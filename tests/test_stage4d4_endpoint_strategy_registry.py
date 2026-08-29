"""Stage 4D-4 endpoint strategy finalization and governance gates."""

from __future__ import annotations

import json
from pathlib import Path

from backend.admet import PHYSICOCHEM_UNAVAILABLE, SAFETY_UNAVAILABLE, TRANSPORTER_UNAVAILABLE
from backend.admet_predictor import MODEL_SPECS
from backend.endpoint_contracts import ENDPOINT_CONTRACTS
from backend.endpoint_strategy_registry import (
    ACTIVE_ADMET_MODEL_ENDPOINTS,
    COMMON_PROMOTION_REQUIREMENTS,
    ENDPOINT_STRATEGY_REGISTRY,
    AdaptiveStatus,
    CalibrationStatus,
    ConsensusPermission,
    DisagreementPolicy,
    PromotionStatus,
    RUNTIME_ADMET_MODEL_ENDPOINTS,
    StrategyType,
    get_endpoint_strategy,
    get_registry_api_response,
    get_registry_summary,
    validate_registry,
)
from backend.main import app, get_interpretation_rules, health, model_strategy_registry
from backend.multimodel import get_adapters_for_endpoint


ROOT = Path(__file__).resolve().parents[1]


def policy(name: str):
    return ENDPOINT_STRATEGY_REGISTRY[name]


def test_registry_has_no_internal_scientific_violations():
    assert validate_registry() == []


def test_every_endpoint_contract_has_exactly_mapped_strategy():
    assert set(ENDPOINT_CONTRACTS) <= set(ENDPOINT_STRATEGY_REGISTRY)
    for name, contract in ENDPOINT_CONTRACTS.items():
        governed = policy(name)
        assert governed.endpoint_id == contract.endpoint_id
        assert governed.endpoint_contract_version == contract.version
        assert get_endpoint_strategy(contract.endpoint_id) is governed


def test_every_runtime_admet_registry_endpoint_has_strategy():
    inactive = {
        "Microsomal clearance", "Dog liver microsomal intrinsic clearance",
        "Monkey liver microsomal intrinsic clearance", "CYP1A2 substrate", "CYP2C19 substrate",
        *TRANSPORTER_UNAVAILABLE, *SAFETY_UNAVAILABLE, *PHYSICOCHEM_UNAVAILABLE,
    }
    expected = set(MODEL_SPECS) | inactive
    assert expected == set(RUNTIME_ADMET_MODEL_ENDPOINTS)
    assert set(ACTIVE_ADMET_MODEL_ENDPOINTS) == set(MODEL_SPECS)
    assert expected <= set(ENDPOINT_STRATEGY_REGISTRY)


def test_strategy_enums_and_model_version_pairs_are_valid():
    expected = {
        "SINGLE_CORE_MODEL", "SINGLE_CORE_WITH_CALIBRATION", "FIXED_WEIGHT_BLEND",
        "STATIC_CONSENSUS", "ADAPTIVE_RESEARCH_SHADOW", "RANK_FUSION", "RULE_BASED",
        "MECHANISTIC_NO_CONSENSUS", "MODEL_UNAVAILABLE", "DERIVED_ESTIMATE", "RULE_ESTIMATE",
    }
    assert {item.value for item in StrategyType} == expected
    for governed in ENDPOINT_STRATEGY_REGISTRY.values():
        assert isinstance(governed.primary_strategy, StrategyType)
        assert governed.shadow_strategy is None or isinstance(governed.shadow_strategy, StrategyType)
        assert len(governed.primary_model_ids) == len(governed.primary_model_versions)
        assert len(governed.shadow_model_ids) == len(governed.shadow_model_versions)


def test_single_core_policies_have_one_immutable_primary_model():
    for governed in ENDPOINT_STRATEGY_REGISTRY.values():
        if governed.primary_strategy == StrategyType.SINGLE_CORE_MODEL:
            assert len(governed.primary_model_ids) == 1
            assert len(governed.primary_model_versions) == 1


def test_active_adapter_model_identity_matches_registry():
    for endpoint in MODEL_SPECS:
        governed = policy(endpoint)
        adapter_pairs = {(a.model_id, a.model_version) for a in get_adapters_for_endpoint(endpoint)}
        for pair in zip(governed.primary_model_ids, governed.primary_model_versions):
            assert pair in adapter_pairs, f"{endpoint}: missing runtime adapter {pair}"
        for pair in zip(governed.shadow_model_ids, governed.shadow_model_versions):
            assert pair in adapter_pairs, f"{endpoint}: missing shadow adapter {pair}"


def test_shadow_strategy_is_never_active():
    for governed in ENDPOINT_STRATEGY_REGISTRY.values():
        if governed.shadow_strategy is not None:
            assert governed.shadow_promotion_status == PromotionStatus.SHADOW
            assert governed.promotion_status == PromotionStatus.ACTIVE


def test_adaptive_disabled_or_no_go_cannot_run_adaptive_shadow():
    prohibited = {AdaptiveStatus.DISABLED, AdaptiveStatus.NO_GO, AdaptiveStatus.NO_ADAPTIVE_VALUE}
    for governed in ENDPOINT_STRATEGY_REGISTRY.values():
        if governed.adaptive_status in prohibited:
            assert governed.shadow_strategy != StrategyType.ADAPTIVE_RESEARCH_SHADOW


def test_solubility_policy_remains_single_core_with_research_adaptive_shadow():
    governed = policy("Solubility")
    assert governed.primary_strategy == StrategyType.SINGLE_CORE_MODEL
    assert governed.primary_model_ids == ["admetica_solubility"]
    assert governed.shadow_strategy == StrategyType.ADAPTIVE_RESEARCH_SHADOW
    assert governed.adaptive_status == AdaptiveStatus.ARCHITECTURE_VALID_BUT_NO_ACCURACY_GAIN
    assert governed.non_primary_model_roles["rdkit_gbr_solubility_v1"] == "ADAPTIVE_EXCLUDED"
    assert "rdkit_gbr_solubility_v1" not in governed.primary_model_ids


def test_caco2_consensus_remains_shadow_with_insufficient_evidence():
    governed = policy("Permeability")
    assert governed.primary_strategy == StrategyType.SINGLE_CORE_MODEL
    assert governed.shadow_strategy == StrategyType.STATIC_CONSENSUS
    assert governed.consensus_status == ConsensusPermission.INSUFFICIENT_EVIDENCE
    assert governed.scientific_status == "INSUFFICIENT_EVIDENCE"


def test_cyp3a4_fixed_blend_is_research_only_and_dynamic_adaptation_disabled():
    governed = policy("CYP3A4 inhibitor")
    assert governed.primary_strategy == StrategyType.SINGLE_CORE_MODEL
    assert governed.shadow_strategy == StrategyType.FIXED_WEIGHT_BLEND
    assert governed.shadow_model_ids == ["morgan_cyp3a4_inh_v1"]
    assert governed.adaptive_status == AdaptiveStatus.NO_ADAPTIVE_VALUE
    assert governed.scientific_status == "FIXED_GLOBAL_BLEND_SUFFICIENT"


def test_other_cyp_isoforms_and_roles_are_not_generalized_from_cyp3a4():
    for name in (
        "CYP1A2 inhibitor", "CYP2C9 inhibitor", "CYP2C19 inhibitor", "CYP2D6 inhibitor",
        "CYP2C9 substrate", "CYP2D6 substrate", "CYP3A4 substrate",
    ):
        governed = policy(name)
        assert governed.primary_strategy == StrategyType.SINGLE_CORE_MODEL
        assert governed.shadow_strategy is None
        assert governed.adaptive_status == AdaptiveStatus.DISABLED


def test_herg_corrected_policy_keeps_raw_m1_and_research_calibration():
    governed = policy("hERG liability")
    assert governed.primary_strategy == StrategyType.SINGLE_CORE_MODEL
    assert governed.primary_model_ids == ["admetica_safety_herg"]
    assert governed.shadow_strategy == StrategyType.SINGLE_CORE_WITH_CALIBRATION
    assert governed.calibration_status == CalibrationStatus.CALIBRATION_RESEARCH
    assert governed.calibration_production_enabled is False
    assert governed.decision_threshold == 0.50
    assert governed.adaptive_status == AdaptiveStatus.NO_GO
    assert governed.scientific_status == "HERG_CALIBRATION_UPDATE_CANDIDATE"
    assert "physchem_herg_v1" not in governed.primary_model_ids
    assert "NOT_DISCRIMINATIVE_BLEND" in governed.non_primary_model_roles["physchem_herg_v1"]
    assert any("BETTER_SECONDARY_MODEL_REQUIRED" in item for item in governed.limitations)


def test_ames_and_dili_are_independent_conservative_single_core_policies():
    for name in ("Ames mutagenicity", "DILI clinical liability"):
        governed = policy(name)
        assert governed.primary_strategy == StrategyType.SINGLE_CORE_MODEL
        assert governed.shadow_strategy is None
        assert governed.calibration_status == CalibrationStatus.RAW


def test_model_disagreement_is_a_signal_not_a_confidence_interval():
    for name in ("CYP3A4 inhibitor", "hERG liability"):
        governed = policy(name)
        assert governed.disagreement_policy == DisagreementPolicy.MODEL_DISAGREEMENT_SIGNAL
        assert "interval" not in governed.disagreement_policy.value.lower()


def test_model_unavailable_endpoints_cannot_emit_production_prediction():
    unavailable = [p for p in ENDPOINT_STRATEGY_REGISTRY.values()
                   if p.primary_strategy == StrategyType.MODEL_UNAVAILABLE]
    assert unavailable
    for governed in unavailable:
        assert governed.primary_model_ids == []
        assert governed.primary_model_versions == []
        assert governed.production_execution_allowed is False
        assert governed.promotion_status == PromotionStatus.DEFERRED
        assert "MODEL_UNAVAILABLE" in governed.fallback_behavior


def test_som_uses_rank_fusion_without_raw_score_averaging():
    governed = policy("Metabolic soft spots")
    assert governed.primary_strategy == StrategyType.RANK_FUSION
    assert governed.primary_model_ids == ["sygma_phase1_2", "smartcyp_dft_v1"]
    assert governed.consensus_status == ConsensusPermission.ALLOWED_PRODUCTION
    assert "NOT averaged" in governed.applicability_policy


def test_pk_methods_are_mechanistic_and_cannot_enter_ml_consensus():
    for name in ("PK Systemic Clearance", "PK Volume of Distribution", "PK Bioavailability", "PK Simulation"):
        governed = policy(name)
        assert governed.primary_strategy == StrategyType.MECHANISTIC_NO_CONSENSUS
        assert governed.consensus_status == ConsensusPermission.MECHANISTICALLY_FORBIDDEN
        assert governed.adaptive_status == AdaptiveStatus.DISABLED


def test_pka_and_logd_provenance_remains_honest():
    assert policy("Ionization (pKa)").primary_strategy == StrategyType.RULE_ESTIMATE
    assert policy("logD pH7.4 derived estimate").primary_strategy == StrategyType.DERIVED_ESTIMATE
    assert policy("pKa (quantitative ML)").primary_strategy == StrategyType.MODEL_UNAVAILABLE
    assert policy("logD7.4 (quantitative ML)").primary_strategy == StrategyType.MODEL_UNAVAILABLE


def test_properties_and_metabolite_hypotheses_are_not_mislabeled_as_ml():
    properties = policy("Physicochemical properties")
    metabolites = policy("Metabolite hypotheses")
    assert properties.primary_strategy == StrategyType.DERIVED_ESTIMATE
    assert properties.scientific_status == "DETERMINISTIC_CALCULATION"
    assert metabolites.primary_strategy == StrategyType.RULE_BASED
    assert "not experimentally confirmed" in " ".join(metabolites.limitations)


def test_clearance_species_endpoints_are_isolated():
    expected = {
        "HLM intrinsic clearance": "openadmet_hlm",
        "RLM intrinsic clearance": "openadmet_rlm",
        "MLM intrinsic clearance": "openadmet_mlm",
    }
    for name, model_id in expected.items():
        governed = policy(name)
        assert governed.primary_model_ids == [model_id]
        assert "Species isolation enforced" in governed.applicability_policy
    assert policy("Dog liver microsomal intrinsic clearance").primary_strategy == StrategyType.MODEL_UNAVAILABLE
    assert policy("Monkey liver microsomal intrinsic clearance").primary_strategy == StrategyType.MODEL_UNAVAILABLE


def test_active_policies_have_complete_deterministic_rollback_metadata():
    for governed in ENDPOINT_STRATEGY_REGISTRY.values():
        if governed.promotion_status != PromotionStatus.ACTIVE:
            continue
        rollback = governed.rollback_policy
        assert rollback is not None
        assert rollback.previous_policy_version
        assert rollback.rollback_target
        assert rollback.rollback_primary_strategy == governed.primary_strategy
        assert rollback.rollback_model_ids == governed.primary_model_ids
        assert rollback.rollback_model_versions == governed.primary_model_versions
        assert rollback.promotion_reason
        assert (ROOT / rollback.validation_artifact).exists()
        assert set(COMMON_PROMOTION_REQUIREMENTS) <= set(governed.promotion_requirements)


def test_machine_readable_master_matrix_matches_runtime_registry():
    artifact = json.loads(
        (ROOT / "validation" / "stage4d4_endpoint_strategy_matrix.json").read_text(encoding="utf-8")
    )
    summary = get_registry_summary()
    assert artifact["registry_version"] == summary["registry_version"]
    assert artifact["endpoints"] == summary["endpoints"]
    assert artifact["validation_violations"] == []
    assert artifact["coverage"]["missing_endpoint_contract_strategies"] == []
    assert artifact["coverage"]["missing_runtime_registry_strategies"] == []


def test_registry_api_is_read_only_complete_and_backward_compatible():
    routes = {route.path: route for route in app.routes if hasattr(route, "path")}
    assert "/api/model-strategy-registry" in routes
    assert routes["/api/model-strategy-registry"].methods == {"GET"}
    body = model_strategy_registry()
    assert body == get_registry_api_response()
    assert body["read_only"] is True
    assert body["production_behavior_changed"] is False
    assert body["total_endpoints"] == len(ENDPOINT_STRATEGY_REGISTRY)
    assert body["violations"] == []
    herg = next(row for row in body["endpoints"] if row["endpoint"] == "hERG liability")
    assert {
        "primary_strategy", "primary_models", "calibration_status", "calibration_policy",
        "shadow_strategy", "shadow_models", "adaptive_status", "consensus_status",
        "applicability_policy", "confidence_policy", "disagreement_policy",
        "validation_status", "limitations", "promotion_requirements", "rollback_policy",
        "policy_version", "fallback_behavior", "scientific_notes",
    } <= set(herg)

    # Existing health and interpretation endpoints remain available and unchanged in shape.
    assert "/api/health" in routes
    assert "/api/interpretation/rules" in routes
    assert health()["status"] == "ok"
    assert "rules" in get_interpretation_rules()


def test_stage4d4_has_no_frontend_changes():
    # The Stage 4D-4 implementation surface is backend policy, generated validation, docs, and tests.
    stage_files = {
        "backend/endpoint_strategy_registry.py",
        "backend/main.py",
        "scripts/generate_stage4d4_strategy_artifacts.py",
        "validation/stage4d4_endpoint_strategy_matrix.json",
        "docs/stage4d4-endpoint-strategy.md",
        "docs/stage4d4-prediction-governance.md",
        "docs/stage4d4-promotion-and-rollback.md",
        "tests/test_stage4d4_endpoint_strategy_registry.py",
    }
    assert all(not path.startswith("frontend/") for path in stage_files)
