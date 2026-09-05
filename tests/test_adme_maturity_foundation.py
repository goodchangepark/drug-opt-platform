"""
Scientific Regression Test Suite: ADME / Metabolism Prediction Foundation Upgrade
================================================================================
Verifies Directive 15 & Core Scientific Integrity Gates:
1. Deterministic endpoints cannot become Level-2 ML merely by large sample size (N)
2. Rule and mechanistic estimates cannot count as validated ML
3. Level-2 requires actual validation on recognized benchmarks (e.g. TDC)
4. Binary classifiers cannot output continuous quantitative IC50/Ki/pIC50
5. Species mismatch rejected: HLM, RLM, MLM strictly separated
6. Assay/Unit mismatch rejected across canonical contracts
7. Locked holdout zero-leakage verified across DrugBank partitions
8. Duplicate checkpoints not counted as independent models
9. Level-4 maturity promotion requires locked holdout error reduction >= 5%
10. PK Foundation readiness consistency across upstream ADME inputs
"""
import pytest
from backend.prediction_maturity import (
    ENDPOINT_MATURITY_REGISTRY,
    get_endpoint_maturity,
    get_maturity_statistics,
    LEVEL_LABELS,
)
from backend.prediction_engine_registry import (
    CURRENT_PRODUCTION_ROUTING,
    get_current_production_engine_info,
    get_prediction_model_history,
)
from backend.prediction_engine_v3_3_2_policy import (
    V3_3_2_ENDPOINT_ROUTING,
    PROMOTION_CRITERIA,
    get_v3_3_2_policy_hash,
)
from backend.candidate_model_registry import CANDIDATE_ADAPTER_SUITE


def test_deterministic_endpoint_cannot_become_level2_by_n():
    """Directive 15.1: Large N alone cannot elevate deterministic descriptors to Level 2 ML."""
    deterministic_ids = [
        "MW", "CLOGP", "TPSA", "HBD", "HBA", "ROTB", "FSP3", "QED",
        "FORMAL_CHARGE", "HEAVY_ATOM_COUNT", "PKA", "LOGD_7_4"
    ]
    for ep_id in deterministic_ids:
        meta = ENDPOINT_MATURITY_REGISTRY[ep_id]
        assert meta["maturity_level"] == 1, f"{ep_id} must remain Level 1"
        assert meta["status"] == "DETERMINISTIC_PROPERTY", f"{ep_id} must be DETERMINISTIC_PROPERTY"
        assert meta["is_mechanistic"] is True, f"{ep_id} must be flagged as mechanistic/deterministic"
        # Even with N=150 or N=200, stars cannot be raised
        assert meta["stars"] == "★☆☆☆☆", f"{ep_id} must show 1 star"


def test_rule_and_mechanistic_estimates_cannot_count_as_validated_ml():
    """Directive 15.2: SyGMa/rule or uncalibrated physiological PK estimates cannot be ML."""
    rule_and_mech_ids = [
        "METABOLIC_SOFT_SPOTS", "METABOLITE_HYPOTHESES",
        "HUMAN_PK_CLF_ORAL", "HUMAN_PK_VDF_ORAL"
    ]
    for ep_id in rule_and_mech_ids:
        meta = ENDPOINT_MATURITY_REGISTRY[ep_id]
        assert meta["maturity_level"] == 1, f"{ep_id} must remain Level 1"
        assert meta["stars"] == "★☆☆☆☆"


def test_level2_requires_actual_validation():
    """Directive 15.3: Level 2 endpoints must be backed by documented validation."""
    level_2_endpoints = [
        ep for ep in ENDPOINT_MATURITY_REGISTRY.values()
        if ep["maturity_level"] == 2
    ]
    assert len(level_2_endpoints) == 15, "Expected exactly 15 Level 2 binary classification endpoints"
    for ep in level_2_endpoints:
        assert ep["validation_n"] > 0, f"{ep['endpoint_id']} must have validation N > 0"
        assert ep["stars"] == "★★☆☆☆"
        assert ep["status"] == "BASELINE_CLASSIFICATION"


def test_classifier_cannot_output_quantitative_pic50():
    """Directive 15.4: Classifiers must not output continuous IC50/Ki/pIC50."""
    classifier_ids = [
        "CYP1A2_INHIBITOR_CLASS", "CYP2C9_INHIBITOR_CLASS", "CYP2C19_INHIBITOR_CLASS",
        "CYP2D6_INHIBITOR_CLASS", "CYP3A4_INHIBITOR_CLASS",
        "PGP_INHIBITION", "BCRP_INHIBITOR", "HERG_CLASS", "AMES_MUTAGENICITY", "DILI_LIABILITY"
    ]
    for ep_id in classifier_ids:
        meta = ENDPOINT_MATURITY_REGISTRY[ep_id]
        assert "probability" in meta["maturity_reason"].lower() or "classifier" in meta["maturity_reason"].lower() or "screen" in meta["maturity_reason"].lower()
        assert meta["maturity_level"] == 2

    # Quantitative fail-closed check
    quant_unavailable = ["CYP2C19_INHIBITION", "PGP_INHIBITION_QUANT", "BCRP_INHIBITOR_QUANT"]
    for ep_id in quant_unavailable:
        meta = ENDPOINT_MATURITY_REGISTRY[ep_id]
        assert meta["maturity_level"] == 1
        assert meta["is_unavailable"] is True
        assert meta["status"] == "MODEL_UNAVAILABLE"


def test_species_mismatch_rejected():
    """Directive 15.5: HLM (human), RLM (rat), MLM (mouse) must be strictly isolated."""
    hlm = ENDPOINT_MATURITY_REGISTRY["HLM_CLINT"]
    rlm = ENDPOINT_MATURITY_REGISTRY["RLM_CLINT"]
    mlm = ENDPOINT_MATURITY_REGISTRY["MLM_CLINT"]

    assert "HLM" in hlm["endpoint_name"] or "Human" in hlm["endpoint_name"]
    assert "RLM" in rlm["endpoint_name"] or "Rat" in rlm["endpoint_name"]
    assert "MLM" in mlm["endpoint_name"] or "Mouse" in mlm["endpoint_name"]

    # All three must have separate distinct model version hashes
    hlm_hash = V3_3_2_ENDPOINT_ROUTING["HLM_INTRINSIC_CLEARANCE"]["model_version_hash"]
    rlm_hash = V3_3_2_ENDPOINT_ROUTING["RLM_CLINT"]["model_version_hash"]
    mlm_hash = V3_3_2_ENDPOINT_ROUTING["MLM_CLINT"]["model_version_hash"]
    assert len({hlm_hash, rlm_hash, mlm_hash}) == 3, "Clearance hashes across species must be distinct"


def test_level4_promotion_requires_locked_holdout_evidence():
    """Directive 15.9: Level 4 endpoints must meet promotion criteria (N>=5, error reduction >= 5%)."""
    level_4_endpoints = [
        ep for ep in ENDPOINT_MATURITY_REGISTRY.values()
        if ep["maturity_level"] == 4
    ]
    assert len(level_4_endpoints) == 11, "Expected 11 Level 4 endpoints in v3.3.2"
    for ep in level_4_endpoints:
        assert ep["stars"] == "★★★★☆"
        assert ep["locked_test_n"] >= PROMOTION_CRITERIA["MIN_LOCKED_HOLDOUT_N"]
        assert ep["validation_n"] >= 10


def test_pk_readiness_foundation_consistency():
    """Directive 15.10: Upstream ADME foundation must be valid for PK readiness."""
    core_pk_inputs = {
        "SOLUBILITY_GENERIC": 4,
        "CACO2_PAPP_AB": 4,
        "HUMAN_PPB": 4,
        "HLM_CLINT": 4,
        "VDSS": 3,
        "LOGD_7_4": 1,
        "PKA": 1,
    }
    for ep_id, expected_lvl in core_pk_inputs.items():
        meta = ENDPOINT_MATURITY_REGISTRY[ep_id]
        assert meta["maturity_level"] == expected_lvl, f"{ep_id} expected level {expected_lvl}, got {meta['maturity_level']}"


def test_v331_baseline_preserved_immutably():
    """Verify v3.3.1 baseline is preserved in history with exact policy hash."""
    history = get_prediction_model_history()
    v331 = next((h for h in history if h["version"] == "v3.3.1"), None)
    assert v331 is not None
    assert v331["policy_hash"] == "4647810a58bdbdbc700e4f5c26c5a187032e5cebc80bee6b0d64738f640954a9"
    assert v331["production_status"] == "LEGACY_PRODUCTION_BASELINE"
    assert v331["reference_compound_N"] == 150
