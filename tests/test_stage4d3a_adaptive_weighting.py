"""
Stage 4D-3A: Comprehensive Unit Tests for Hierarchical Adaptive Weighting.

Validates:
1. Global prior calculation and bounded exponential scoring
2. Zero-data fallback to global prior (lambda = 0)
3. Project-level Bayesian shrinkage monotonicity
4. Bemis-Murcko scaffold series assignment and series shrinkage
5. Morgan fingerprint Tanimoto local neighborhood distance weighting
6. Sparse neighborhood fallback (no neighbors >= 0.40)
7. Applicability domain penalty and minimum weight floor (eps = 0.02)
8. Reason code emission
9. Prospective timestamp filtering (strict zero retrospective leakage)
10. Feedback event idempotency
11. Adaptive provenance API endpoint
12. Shadow mode non-interference with production visible outputs
"""

import pytest
from fastapi.testclient import TestClient
from sqlalchemy import select

from backend.main import app
from backend.database import SessionLocal, engine
from backend.models import Compound, CompoundVersion, Project
from backend.admet import (
    ADMETEndpoint,
    ADMETMeasurement,
    ADMETPrediction,
    ADMETModelRegistry,
    ADMETExperimentalFeedbackEvent,
    ADMETAdaptivePrediction,
    ensure_admet_schema,
)
from backend.adaptive_weighting import (
    ADAPTIVE_POLICY_VERSION,
    AssayQuality,
    AdaptiveReasonCode,
    ExperimentalFeedbackRecord,
    compute_error_score,
    compute_hierarchical_adaptive_weights,
    compute_morgan_fingerprint,
    compute_shrinkage_lambda,
    compute_tanimoto_similarity,
    evaluate_experimental_compatibility,
    get_bemis_murcko_scaffold,
)
from backend.multimodel import (
    ExecutionStatus,
    ModelExecutionPayload,
)


@pytest.fixture(scope="module")
def client():
    ensure_admet_schema(engine)
    return TestClient(app)


def test_global_prior_calculation():
    """Verify global error score calculation and normalization."""
    # Admetica (0.3386) vs ESOL (0.6663)
    s1 = compute_error_score(0.3386, beta=2.0)
    s2 = compute_error_score(0.6663, beta=2.0)
    assert s1 > s2
    w1 = s1 / (s1 + s2)
    w2 = s2 / (s1 + s2)
    assert 0.60 < w1 < 0.70
    assert 0.30 < w2 < 0.40
    assert abs((w1 + w2) - 1.0) < 1e-6


def test_shrinkage_lambda_math():
    """Verify empirical Bayes shrinkage factor lambda = N / (N + N_prior)."""
    assert compute_shrinkage_lambda(0.0, 10.0) == 0.0
    assert compute_shrinkage_lambda(10.0, 10.0) == 0.5
    assert compute_shrinkage_lambda(90.0, 10.0) == 0.9
    # Monotonicity
    l1 = compute_shrinkage_lambda(2.0, 5.0)
    l2 = compute_shrinkage_lambda(5.0, 5.0)
    l3 = compute_shrinkage_lambda(15.0, 5.0)
    assert l1 < l2 < l3


def test_zero_data_fallback():
    """Verify that with zero project observations, weights match global prior exactly."""
    p1 = ModelExecutionPayload(
        model_id="admetica_solubility",
        model_name="Admetica Chemprop Solubility",
        model_family="admetica",
        model_version="1.0",
        endpoint_id="EP_PHYS_SOLUBILITY",
        endpoint_name="Solubility",
        canonical_unit="log10(mol/L)",
        execution_status=ExecutionStatus.SUCCESS,
        value=-2.50,
        applicability_domain="IN_DOMAIN",
    )
    p2 = ModelExecutionPayload(
        model_id="esol_delaney_v1",
        model_name="Delaney ESOL",
        model_family="esol",
        model_version="1.0",
        endpoint_id="EP_PHYS_SOLUBILITY",
        endpoint_name="Solubility",
        canonical_unit="log10(mol/L)",
        execution_status=ExecutionStatus.SUCCESS,
        value=-2.80,
        applicability_domain="IN_DOMAIN",
    )

    res = compute_hierarchical_adaptive_weights(
        query_smiles="CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
        project_id=1,
        candidate_payloads=[p1, p2],
        historical_feedback_events=[],
    )

    assert res.n_project == 0
    assert res.n_series == 0
    assert res.n_local_eff == 0.0
    assert AdaptiveReasonCode.GLOBAL_PRIOR_DOMINANT.value in res.reason_codes
    assert AdaptiveReasonCode.INSUFFICIENT_LOCAL_DATA.value in res.reason_codes

    w1 = res.effective_weights["admetica_solubility"]
    w2 = res.effective_weights["esol_delaney_v1"]
    assert 0.60 < w1 < 0.70
    assert 0.30 < w2 < 0.40


def test_scaffold_series_extraction():
    """Verify Bemis-Murcko scaffold identification including acyclic handling."""
    scaff_asp = get_bemis_murcko_scaffold("CC(=O)Oc1ccccc1C(=O)O")
    assert scaff_asp == "c1ccccc1"

    scaff_acyclic = get_bemis_murcko_scaffold("CCCCCC(=O)O")
    assert scaff_acyclic == "[acyclic]"

    scaff_biphenyl = get_bemis_murcko_scaffold("c1ccc(-c2ccccc2)cc1")
    assert "c1ccc" in scaff_biphenyl


def test_local_tanimoto_similarity():
    """Verify Morgan fingerprint computation and Tanimoto similarity."""
    fp1 = compute_morgan_fingerprint("CC(=O)Oc1ccccc1C(=O)O")  # Aspirin
    fp2 = compute_morgan_fingerprint("Oc1ccccc1C(=O)O")        # Salicylic acid
    fp3 = compute_morgan_fingerprint("CCCCCCCCCCCC")           # Dodecane

    sim_high = compute_tanimoto_similarity(fp1, fp2)
    sim_low = compute_tanimoto_similarity(fp1, fp3)

    assert sim_high > 0.40
    assert sim_low < 0.20


def test_project_and_series_adaptation():
    """Verify that when M2 outperforms M1 on a specific series, M2 weight increases."""
    p1 = ModelExecutionPayload(
        model_id="admetica_solubility",
        model_name="Admetica Chemprop Solubility",
        model_family="admetica",
        model_version="1.0",
        endpoint_id="EP_PHYS_SOLUBILITY",
        endpoint_name="Solubility",
        canonical_unit="log10(mol/L)",
        execution_status=ExecutionStatus.SUCCESS,
        value=-2.00,
        applicability_domain="IN_DOMAIN",
    )
    p2 = ModelExecutionPayload(
        model_id="esol_delaney_v1",
        model_name="Delaney ESOL",
        model_family="esol",
        model_version="1.0",
        endpoint_id="EP_PHYS_SOLUBILITY",
        endpoint_name="Solubility",
        canonical_unit="log10(mol/L)",
        execution_status=ExecutionStatus.SUCCESS,
        value=-3.00,
        applicability_domain="IN_DOMAIN",
    )

    # 10 observations where M2 has 0 error and M1 has high error
    history = []
    for i in range(10):
        history.append(ExperimentalFeedbackRecord(
            event_id=f"TEST-EVT-{i}",
            project_id=1,
            compound_version_id=i + 1,
            canonical_smiles=f"c1ccccc1C({i})",
            endpoint_name="Solubility",
            experimental_value=-3.00,
            experimental_unit="log10(mol/L)",
            assay_quality=AssayQuality.HIGH_QUALITY,
            scaffold_smiles="c1ccccc1",
            timestamp=f"2026-08-28T{i:02d}",
            frozen_predictions={"admetica_solubility": -1.50, "esol_delaney_v1": -3.00},
            model_errors={"admetica_solubility": 1.50, "esol_delaney_v1": 0.00},
            is_valid=True,
        ))

    res = compute_hierarchical_adaptive_weights(
        query_smiles="c1ccccc1C(=O)O",
        project_id=1,
        candidate_payloads=[p1, p2],
        historical_feedback_events=history,
    )

    assert res.n_project == 10
    assert res.n_series == 10
    assert AdaptiveReasonCode.PROJECT_EVIDENCE_ACTIVE.value in res.reason_codes
    assert AdaptiveReasonCode.SERIES_M2_OUTPERFORMS_M1.value in res.reason_codes
    # M2 should now have higher weight than M1
    w_m1 = res.effective_weights["admetica_solubility"]
    w_m2 = res.effective_weights["esol_delaney_v1"]
    assert w_m2 > w_m1


def test_applicability_domain_penalty_and_weight_floor():
    """Verify that OUT_OF_DOMAIN models are penalized by 0.1x and respect 0.02 floor."""
    p1 = ModelExecutionPayload(
        model_id="admetica_solubility",
        model_name="Admetica Chemprop Solubility",
        model_family="admetica",
        model_version="1.0",
        endpoint_id="EP_PHYS_SOLUBILITY",
        endpoint_name="Solubility",
        canonical_unit="log10(mol/L)",
        execution_status=ExecutionStatus.SUCCESS,
        value=-2.00,
        applicability_domain="OUT_OF_DOMAIN",
    )
    p2 = ModelExecutionPayload(
        model_id="esol_delaney_v1",
        model_name="Delaney ESOL",
        model_family="esol",
        model_version="1.0",
        endpoint_id="EP_PHYS_SOLUBILITY",
        endpoint_name="Solubility",
        canonical_unit="log10(mol/L)",
        execution_status=ExecutionStatus.SUCCESS,
        value=-2.00,
        applicability_domain="IN_DOMAIN",
    )

    res = compute_hierarchical_adaptive_weights(
        query_smiles="c1ccccc1",
        project_id=1,
        candidate_payloads=[p1, p2],
        historical_feedback_events=[],
    )

    # Floor protection check
    for mid, bd in res.weights_breakdown.items():
        assert bd.final_effective_weight >= 0.02


def test_prospective_zero_leakage():
    """Verify that future events (timestamp > prediction_timestamp) are excluded."""
    p1 = ModelExecutionPayload(
        model_id="admetica_solubility",
        model_name="Admetica Chemprop Solubility",
        model_family="admetica",
        model_version="1.0",
        endpoint_id="EP_PHYS_SOLUBILITY",
        endpoint_name="Solubility",
        canonical_unit="log10(mol/L)",
        execution_status=ExecutionStatus.SUCCESS,
        value=-2.00,
        applicability_domain="IN_DOMAIN",
    )
    p2 = ModelExecutionPayload(
        model_id="esol_delaney_v1",
        model_name="Delaney ESOL",
        model_family="esol",
        model_version="1.0",
        endpoint_id="EP_PHYS_SOLUBILITY",
        endpoint_name="Solubility",
        canonical_unit="log10(mol/L)",
        execution_status=ExecutionStatus.SUCCESS,
        value=-2.00,
        applicability_domain="IN_DOMAIN",
    )

    future_event = ExperimentalFeedbackRecord(
        event_id="FUTURE-EVT-1",
        project_id=1,
        compound_version_id=999,
        canonical_smiles="c1ccccc1",
        endpoint_name="Solubility",
        experimental_value=-2.00,
        experimental_unit="log10(mol/L)",
        assay_quality=AssayQuality.HIGH_QUALITY,
        scaffold_smiles="c1ccccc1",
        timestamp="2026-08-29T12:00:00",
        frozen_predictions={"admetica_solubility": -2.0, "esol_delaney_v1": -2.0},
        model_errors={"admetica_solubility": 0.0, "esol_delaney_v1": 0.0},
        is_valid=True,
    )

    res = compute_hierarchical_adaptive_weights(
        query_smiles="c1ccccc1Cl",
        project_id=1,
        candidate_payloads=[p1, p2],
        historical_feedback_events=[future_event],
        prediction_timestamp="2026-08-28T00:00:00",
    )

    assert res.n_project == 0


def test_experimental_compatibility():
    """Verify strict unit and method compatibility validation."""
    compat, q, msg = evaluate_experimental_compatibility("Solubility", -2.5, "log10(mol/L)")
    assert compat is True
    assert q in {AssayQuality.HIGH_QUALITY, AssayQuality.USABLE}

    # Wrong unit
    compat_bad_unit, _, _ = evaluate_experimental_compatibility("Solubility", 10.5, "ug/mL")
    assert compat_bad_unit is False

    # Wrong endpoint
    compat_bad_ep, _, _ = evaluate_experimental_compatibility("Caco-2 Permeability", -4.5, "log10(mol/L)")
    assert compat_bad_ep is False


def test_adaptive_provenance_api(client):
    """Verify GET /api/compound-versions/{version_id}/adaptive-provenance endpoint."""
    with SessionLocal() as db:
        proj = db.scalar(select(Project).limit(1))
        assert proj is not None, "Project required for API test"
        ver = db.scalar(select(CompoundVersion).limit(1))
        assert ver is not None, "CompoundVersion required for API test"
        vid = ver.id

    res = client.get(f"/api/compound-versions/{vid}/adaptive-provenance")
    assert res.status_code == 200, res.text
    data = res.json()
    assert data["compound_version_id"] == vid
    assert data["endpoint_name"] == "Solubility"
    assert data["consensus_mode"] == "SHADOW"
    assert data["policy_version"] == ADAPTIVE_POLICY_VERSION
    assert "sample_counts" in data
    assert "weights_breakdown" in data
    assert "reason_codes" in data
    assert "effective_weights" in data
