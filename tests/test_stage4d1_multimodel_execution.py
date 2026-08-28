"""
Drug-OPT Stage 4D-1: Multi-Model Execution, Storage & Consensus Foundation Tests.

Covers:
1. Multi-model storage and model version isolation
2. Endpoint contract and compatibility enforcement
3. Species isolation gate
4. Classification semantics and vote pattern
5. Failure isolation and transparent weight renormalization
6. Applicability-domain-aware weighting and diversity penalties
7. Regression dispersion and model disagreement metrics
8. Reciprocal Rank Fusion (RRF) for metabolic soft spots
9. Shadow mode operation and backward-compatible API
10. Database schema integrity and artifact validation
"""

import json
from pathlib import Path
import pytest
from sqlalchemy import select, inspect

from backend.database import SessionLocal, engine
from backend.admet import (
    ADMETModelRegistry,
    ADMETPrediction,
    ADMETConsensusPrediction,
    ADMETEndpoint,
    ensure_admet_schema,
)
from backend.endpoint_contracts import (
    ENDPOINT_CONTRACTS,
    EndpointContract,
    OutputType,
    Directionality,
    ExecutionTier,
    ARM64Status,
    get_endpoint_contract,
    check_ensemble_compatibility,
)
from backend.multimodel import (
    ExecutionStatus,
    ModelExecutionPayload,
    AdmeticaChempropAdapter,
    OpenADMETClearanceAdapter,
    ADMETAISafetyAdapter,
    SyGMaMetabolismAdapter,
    list_registered_adapters,
    get_adapters_for_endpoint,
    compute_prediction_cache_key,
)
from backend.consensus import (
    ConsensusMode,
    AggregationType,
    AgreementStatus,
    calculate_static_model_weight,
    compute_endpoint_consensus,
)


def test_schema_migration_and_db_columns():
    """Verify that multi-model columns exist in database tables."""
    ensure_admet_schema(engine)
    inspector = inspect(engine)
    pred_cols = {c["name"] for c in inspector.get_columns("admet_predictions")}
    cons_cols = {c["name"] for c in inspector.get_columns("admet_consensus_predictions")}

    # admet_predictions
    assert "model_version" in pred_cols
    assert "execution_status" in pred_cols
    assert "standardizer_version" in pred_cols
    assert "canonical_smiles" in pred_cols
    assert "runtime_ms" in pred_cols

    # admet_consensus_predictions
    assert "consensus_version" in cons_cols
    assert "consensus_mode" in cons_cols
    assert "model_agreement" in cons_cols
    assert "dispersion_json" in cons_cols
    assert "vote_pattern" in cons_cols


def test_registered_model_adapters():
    """Verify default model adapters are initialized and available."""
    adapters = list_registered_adapters()
    assert len(adapters) >= 19

    # Admetica Solubility
    sol_adapters = get_adapters_for_endpoint("Solubility")
    assert len(sol_adapters) >= 1
    assert sol_adapters[0].model_family == "admetica"
    assert sol_adapters[0].execution_tier == ExecutionTier.TIER_1_LOCAL_FAST
    assert sol_adapters[0].arm64_status == ARM64Status.RUNS_LOCAL_ARM64

    # OpenADMET HLM
    hlm_adapters = get_adapters_for_endpoint("HLM intrinsic clearance")
    assert len(hlm_adapters) >= 1
    assert hlm_adapters[0].model_family == "openadmet_clearance"

    # ADMET-AI Ames
    ames_adapters = get_adapters_for_endpoint("Ames mutagenicity")
    assert len(ames_adapters) >= 1
    assert ames_adapters[0].model_family == "admet_ai_ensemble"

    # SyGMa
    sygma_adapters = get_adapters_for_endpoint("Metabolic soft spots")
    assert len(sygma_adapters) >= 1
    assert sygma_adapters[0].model_family == "rule_based_smarts"


def test_endpoint_contract_and_species_isolation():
    """Verify that incompatible endpoints and cross-species models cannot ensemble."""
    sol_contract = get_endpoint_contract("Solubility")
    hlm_contract = get_endpoint_contract("HLM intrinsic clearance")
    rlm_contract = get_endpoint_contract("RLM intrinsic clearance")

    # Solubility vs HLM (incompatible endpoint IDs / category mismatch)
    compat, reason = check_ensemble_compatibility(sol_contract, hlm_contract)
    assert not compat
    assert "Incompatible" in reason or "Category mismatch" in reason

    # HLM vs RLM (incompatible species / endpoint IDs)
    compat_species, reason_species = check_ensemble_compatibility(hlm_contract, rlm_contract)
    assert not compat_species
    assert "Incompatible" in reason_species or "Species mismatch" in reason_species


def test_fault_tolerant_failure_renormalization():
    """Verify that a failing model is isolated and weights renormalize cleanly."""
    p1 = ModelExecutionPayload(
        model_id="model_1",
        model_name="Model 1",
        model_family="family_a",
        model_version="1.0",
        endpoint_id="ADMET-EXP-001",
        endpoint_name="Solubility",
        canonical_unit="log(mol/L)",
        execution_status=ExecutionStatus.SUCCESS,
        value=-2.50,
        applicability_domain="IN_DOMAIN",
        confidence="HIGH",
        canonical_smiles="CC(=O)Oc1ccccc1C(=O)O",
    )
    p2_failed = ModelExecutionPayload(
        model_id="model_2_fault",
        model_name="Model 2 (Failed)",
        model_family="family_b",
        model_version="1.0",
        endpoint_id="ADMET-EXP-001",
        endpoint_name="Solubility",
        canonical_unit="log(mol/L)",
        execution_status=ExecutionStatus.RUNTIME_ERROR,
        error_message="CUDA/PyTorch tensor allocation error",
        canonical_smiles="CC(=O)Oc1ccccc1C(=O)O",
    )
    p3 = ModelExecutionPayload(
        model_id="model_3",
        model_name="Model 3",
        model_family="family_c",
        model_version="1.0",
        endpoint_id="ADMET-EXP-001",
        endpoint_name="Solubility",
        canonical_unit="log(mol/L)",
        execution_status=ExecutionStatus.SUCCESS,
        value=-2.70,
        applicability_domain="IN_DOMAIN",
        confidence="HIGH",
        canonical_smiles="CC(=O)Oc1ccccc1C(=O)O",
    )

    consensus = compute_endpoint_consensus("Solubility", 101, [p1, p2_failed, p3])

    assert consensus.consensus_mode == ConsensusMode.SHADOW
    assert len(consensus.models_used) == 2
    assert "model_1" in consensus.models_used
    assert "model_3" in consensus.models_used
    assert "model_2_fault" not in consensus.models_used
    assert pytest.approx(sum(consensus.effective_weights.values()), 0.001) == 1.0
    assert pytest.approx(consensus.combined_value, 0.05) == -2.60
    assert any("Failure isolation" in w for w in consensus.warnings)


def test_applicability_domain_and_diversity_penalty():
    """Verify that OOD models are downweighted and shared-architecture models penalized."""
    # Model in-domain vs model OOD
    p_in = ModelExecutionPayload(
        model_id="model_in",
        model_name="Model In",
        model_family="admetica",
        model_version="1.0",
        endpoint_id="ADMET-EXP-001",
        endpoint_name="Solubility",
        canonical_unit="log(mol/L)",
        execution_status=ExecutionStatus.SUCCESS,
        value=-2.0,
        applicability_domain="IN_DOMAIN",
        confidence="HIGH",
    )
    p_ood = ModelExecutionPayload(
        model_id="model_ood",
        model_name="Model OOD",
        model_family="custom_dnn",
        model_version="1.0",
        endpoint_id="ADMET-EXP-001",
        endpoint_name="Solubility",
        canonical_unit="log(mol/L)",
        execution_status=ExecutionStatus.SUCCESS,
        value=-4.0,
        applicability_domain="OUT_OF_DOMAIN",
        confidence="LOW",
    )

    w_in, _ = calculate_static_model_weight(p_in, [p_in, p_ood])
    w_ood, _ = calculate_static_model_weight(p_ood, [p_in, p_ood])
    assert w_in > w_ood * 5.0  # IN_DOMAIN significantly outweighs OUT_OF_DOMAIN

    # Test diversity penalty for shared Admetica + ADMET-AI pair
    p_admetica = ModelExecutionPayload(
        model_id="admetica_cyp",
        model_name="Admetica CYP3A4",
        model_family="admetica",
        model_version="1.0",
        endpoint_id="ADMET-EXP-008",
        endpoint_name="CYP3A4 inhibitor",
        canonical_unit="probability",
        execution_status=ExecutionStatus.SUCCESS,
        probability=0.80,
    )
    p_admet_ai = ModelExecutionPayload(
        model_id="admet_ai_cyp",
        model_name="ADMET-AI CYP3A4",
        model_family="admet_ai",
        model_version="1.0",
        endpoint_id="ADMET-EXP-008",
        endpoint_name="CYP3A4 inhibitor",
        canonical_unit="probability",
        execution_status=ExecutionStatus.SUCCESS,
        probability=0.75,
    )
    w_pen, r_pen = calculate_static_model_weight(p_admetica, [p_admetica, p_admet_ai])
    assert "Diversity(0.55)" in r_pen


def test_regression_dispersion_and_model_disagreement():
    """Verify regression consensus calculates weighted std as model disagreement."""
    p1 = ModelExecutionPayload(
        model_id="m1", model_name="M1", model_family="fam1", model_version="1.0",
        endpoint_id="ADMET-EXP-001", endpoint_name="Solubility", canonical_unit="log(mol/L)",
        execution_status=ExecutionStatus.SUCCESS, value=-2.0,
    )
    p2 = ModelExecutionPayload(
        model_id="m2", model_name="M2", model_family="fam2", model_version="1.0",
        endpoint_id="ADMET-EXP-001", endpoint_name="Solubility", canonical_unit="log(mol/L)",
        execution_status=ExecutionStatus.SUCCESS, value=-2.4,
    )
    consensus = compute_endpoint_consensus("Solubility", 101, [p1, p2])

    assert consensus.aggregation_type == AggregationType.REGRESSION_WEIGHTED
    assert pytest.approx(consensus.combined_value, 0.01) == -2.20
    assert consensus.dispersion["model_disagreement_std"] > 0.0
    assert consensus.dispersion["range"] == 0.40
    assert consensus.model_agreement in {AgreementStatus.HIGH_AGREEMENT, AgreementStatus.MODERATE_AGREEMENT}
    assert "MODEL DISAGREEMENT" in consensus.dispersion["interpretation"]


def test_classification_consensus_and_vote_pattern():
    """Verify classification consensus computes weighted probability and vote pattern."""
    p1 = ModelExecutionPayload(
        model_id="m1", model_name="M1", model_family="fam1", model_version="1.0",
        endpoint_id="ADMET-EXP-008", endpoint_name="CYP3A4 inhibitor", canonical_unit="probability",
        execution_status=ExecutionStatus.SUCCESS, probability=0.85,
    )
    p2 = ModelExecutionPayload(
        model_id="m2", model_name="M2", model_family="fam2", model_version="1.0",
        endpoint_id="ADMET-EXP-008", endpoint_name="CYP3A4 inhibitor", canonical_unit="probability",
        execution_status=ExecutionStatus.SUCCESS, probability=0.75,
    )
    consensus = compute_endpoint_consensus("CYP3A4 inhibitor", 101, [p1, p2])

    assert consensus.aggregation_type == AggregationType.CLASSIFICATION_WEIGHTED
    assert consensus.consensus_classification == "INHIBITOR"
    assert pytest.approx(consensus.combined_probability, 0.01) == 0.80
    assert consensus.model_agreement == AgreementStatus.HIGH_AGREEMENT
    assert "m1:INHIBITOR" in consensus.vote_pattern
    assert "m2:INHIBITOR" in consensus.vote_pattern


def test_metabolic_soft_spots_rank_fusion():
    """Verify that metabolic soft spots use Reciprocal Rank Fusion instead of value averaging."""
    p_sygma = ModelExecutionPayload(
        model_id="sygma_phase1_2",
        model_name="SyGMa Phase I & II",
        model_family="rule_based_smarts",
        model_version="1.1.0",
        endpoint_id="ADMET-MET-001",
        endpoint_name="Metabolic soft spots",
        canonical_unit="atom_index_ranking",
        execution_status=ExecutionStatus.SUCCESS,
        raw_outputs={
            "spots": [
                {"atom_index": 3, "atom_environment": "aliphatic_hydroxylation", "reactions": ["Hydroxylation"]},
                {"atom_index": 5, "atom_environment": "aromatic_hydroxylation", "reactions": ["Aromatic hydroxylation"]},
            ]
        },
    )
    consensus = compute_endpoint_consensus("Metabolic soft spots", 101, [p_sygma])

    assert consensus.aggregation_type == AggregationType.RANK_FUSION
    assert "top_soft_spots" in consensus.dispersion
    top_spots = consensus.dispersion["top_soft_spots"]
    assert len(top_spots) == 2
    assert top_spots[0]["atom_index"] == 3
    assert top_spots[0]["rank"] == 1
    assert top_spots[0]["rrf_score"] == round(1.0 / (60.0 + 1), 4)


def test_prediction_cache_key_generation():
    """Verify cache key uniqueness across compound versions and models."""
    k1 = compute_prediction_cache_key(101, "CC(=O)Oc1ccccc1C(=O)O", "ADMET-EXP-001", "admetica_solubility", "2.1.0")
    k2 = compute_prediction_cache_key(102, "CC(=O)Oc1ccccc1C(=O)O", "ADMET-EXP-001", "admetica_solubility", "2.1.0")
    k3 = compute_prediction_cache_key(101, "CC(=O)Oc1ccccc1C(=O)O", "ADMET-EXP-001", "admetica_solubility", "2.2.0")

    assert k1 != k2  # Different compound version
    assert k1 != k3  # Different model version
    assert len(k1) == 64  # Valid SHA-256 hex string


def test_stage4d1_artifacts_and_documentation():
    """Verify all machine-readable JSON files and documentation markdown files exist."""
    req_json = [
        "validation/stage4d1_adapter_registry.json",
        "validation/stage4d1_runtime_benchmark.json",
        "validation/stage4d1_consensus_policy.json",
        "validation/stage4d1_shadow_validation.json",
    ]
    for r in req_json:
        p = Path(r)
        assert p.exists(), f"Missing required validation artifact: {r}"
        with open(p) as f:
            data = json.load(f)
            assert data is not None

    req_docs = [
        "docs/stage4d1-multimodel-execution.md",
        "docs/stage4d1-consensus-policy.md",
        "docs/stage4d1-shadow-mode.md",
    ]
    for d in req_docs:
        p = Path(d)
        assert p.exists(), f"Missing required documentation: {d}"
        assert p.stat().st_size > 500


def test_multimodel_provenance_api_and_shadow_mode_e2e():
    """Verify the /api/compound-versions/{version_id}/multimodel-provenance endpoint and shadow mode."""
    from fastapi.testclient import TestClient
    from backend.main import app
    client = TestClient(app)

    # 1. Create temporary project & compound
    proj_resp = client.post("/api/projects", json={
        "name": "Stage4D1 Test Project",
        "description": "Multi-model shadow consensus validation",
        "target": "TEST",
        "molecule_type": "Small Molecule"
    })
    assert proj_resp.status_code == 201
    proj_id = proj_resp.json()["id"]

    try:
        comp_resp = client.post(f"/api/projects/{proj_id}/compounds", json={
            "compound_id": "TEST-4D1-001",
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
            "notes": "Stage 4D-1 test compound"
        })
        assert comp_resp.status_code == 201
        comp_data = comp_resp.json()
        compound_row_id = comp_data["row_id"]
        version_id = comp_data["version"]["id"]

        # 2. Run Save & Predict (existing API)
        pred_resp = client.post(f"/api/admet/predict/{version_id}")
        assert pred_resp.status_code == 202
        pred_data = pred_resp.json()
        assert pred_data["status"] in {"COMPLETE", "CACHED", "PARTIAL"}
        assert len(pred_data["predictions"]) >= 18

        # 3. Call new Multi-Model Provenance API
        prov_resp = client.get(f"/api/compound-versions/{version_id}/multimodel-provenance")
        assert prov_resp.status_code == 200
        prov_data = prov_resp.json()

        assert prov_data["compound_version_id"] == version_id
        assert prov_data["consensus_mode"] == "SHADOW"
        assert prov_data["total_endpoints"] >= 18

        # Verify Solubility endpoint structure
        sol_ep = next((e for e in prov_data["endpoints"] if e["endpoint_name"] == "Solubility"), None)
        assert sol_ep is not None
        assert sol_ep["canonical_unit"] == "log10(mol/L)"
        assert len(sol_ep["models"]) >= 1
        assert sol_ep["consensus"] is not None
        assert sol_ep["consensus"]["consensus_mode"] == "SHADOW"

    finally:
        # Clean up temporary test project
        client.delete(f"/api/projects/{proj_id}")
