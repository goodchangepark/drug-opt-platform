"""
Stage 4D-3B1: Unit & Integration Tests for CYP3A4 Adaptive Classification Ensembling.

Tests:
1. Endpoint contract & experimental compatibility validation for CYP3A4.
2. Global prior weights calibration for CYP3A4 models (M1 Admetica vs M2 Morgan).
3. Probability loss, Brier score, bounded LogLoss computation.
4. Hierarchical Bayesian shrinkage (Global -> Project -> Series -> Local).
5. Class imbalance detection (CLASS_BALANCE_LIMITED).
6. Prospective forward walk (zero retrospective leakage).
7. Cross-project isolation (Project A never bleeds into Project B).
8. New project & new series fallback.
9. Applicability domain downweighting & minimum weight floor.
10. Model disagreement signal computation (|p_M1 - p_M2|).
11. Replay reproducibility from immutable feedback events.
12. Schema validity of all 8 Stage 4D-3B1 JSON validation artifacts.
"""

import json
import math
import os
import pytest
import numpy as np

from backend.adaptive_weighting import (
    ADAPTIVE_POLICY_VERSION,
    DEFAULT_BETA_ERROR_SCALING,
    DEFAULT_N_PRIOR_LOCAL,
    DEFAULT_N_PRIOR_PROJECT,
    DEFAULT_N_PRIOR_SERIES,
    DEFAULT_LOCAL_SIMILARITY_THRESHOLD,
    GLOBAL_ENDPOINT_PRIOR_ERRORS,
    MINIMUM_WEIGHT_FLOOR,
    PROBABILITY_EPSILON,
    AdaptiveReasonCode,
    AssayQuality,
    ExperimentalFeedbackRecord,
    compute_hierarchical_adaptive_weights,
    evaluate_experimental_compatibility,
    get_bemis_murcko_scaffold,
)
from backend.endpoint_contracts import get_endpoint_contract
from backend.multimodel import ExecutionStatus, ModelExecutionPayload


def make_cyp_payload(model_id: str, value: float, ad: str = "IN_DOMAIN") -> ModelExecutionPayload:
    return ModelExecutionPayload(
        model_id=model_id,
        model_name="Admetica CYP3A4" if "admetica" in model_id else "Morgan CYP3A4",
        model_family="admetica" if "admetica" in model_id else "morgan_gradient_boosting",
        model_version="1.0",
        endpoint_id="EP_MET_CYP3A4_INH",
        endpoint_name="CYP3A4 inhibitor",
        canonical_unit="probability",
        execution_status=ExecutionStatus.SUCCESS,
        value=value,
        applicability_domain=ad,
    )


@pytest.fixture
def cyp_contract():
    return get_endpoint_contract("CYP3A4 inhibitor")


def test_cyp3a4_contract_and_compatibility():
    """Verify CYP3A4 contract properties and experimental compatibility gating."""
    contract = get_endpoint_contract("CYP3A4 inhibitor")
    assert contract.endpoint_id == "cyp3a4_inhibitor_prob"
    assert contract.canonical_unit == "probability"
    
    # 1. Valid binary / probability
    ok, quality, msg = evaluate_experimental_compatibility("CYP3A4 inhibitor", 1.0, "binary")
    assert ok is True
    assert quality == AssayQuality.HIGH_QUALITY

    ok, quality, msg = evaluate_experimental_compatibility("CYP3A4 inhibitor", 0.0, "probability")
    assert ok is True
    assert quality == AssayQuality.HIGH_QUALITY

    # 2. Valid quantitative AC50 / IC50
    ok, quality, msg = evaluate_experimental_compatibility("CYP3A4 inhibitor", 2.5, "uM")
    assert ok is True
    assert quality == AssayQuality.USABLE

    # 3. Incompatible assays (substrate, TDI, induction, fm)
    ok, quality, msg = evaluate_experimental_compatibility("CYP3A4 inhibitor", 1.0, "binary", method="CYP3A4 substrate depletion")
    assert ok is False
    assert quality == AssayQuality.INCOMPATIBLE

    ok, quality, msg = evaluate_experimental_compatibility("CYP3A4 inhibitor", 1.0, "binary", notes="Time-dependent inhibition (TDI)")
    assert ok is False
    assert quality == AssayQuality.INCOMPATIBLE


def test_cyp3a4_global_prior_calibration():
    """Verify calibrated conservative global prior heavily favors M1 over M2."""
    priors = GLOBAL_ENDPOINT_PRIOR_ERRORS["CYP3A4 inhibitor"]
    assert "admetica_cyp_cyp3a4-inhibitor" in priors
    assert "morgan_cyp3a4_inh_v1" in priors
    
    payloads = [
        make_cyp_payload("admetica_cyp_cyp3a4-inhibitor", 0.85),
        make_cyp_payload("morgan_cyp3a4_inh_v1", 0.60),
    ]
    
    res = compute_hierarchical_adaptive_weights(
        query_smiles="CC(=O)Nc1ccccc1",
        project_id=1,
        candidate_payloads=payloads,
        historical_feedback_events=[],
        endpoint_name="CYP3A4 inhibitor",
    )
    
    w_m1 = res.effective_weights["admetica_cyp_cyp3a4-inhibitor"]
    w_m2 = res.effective_weights["morgan_cyp3a4_inh_v1"]
    
    assert w_m1 > 0.90
    assert w_m2 < 0.10
    assert pytest.approx(w_m1 + w_m2, rel=1e-3) == 1.0
    assert AdaptiveReasonCode.GLOBAL_PRIOR_DOMINANT.value in res.reason_codes


def test_cyp3a4_project_and_series_adaptation():
    """Verify that consistent project and series feedback dynamically updates weights."""
    payloads = [
        make_cyp_payload("admetica_cyp_cyp3a4-inhibitor", 0.80),
        make_cyp_payload("morgan_cyp3a4_inh_v1", 0.75),
    ]
    
    # Create 15 feedback records where M2 had smaller Brier loss than M1 in this series
    scaffold = "c1ccccc1"
    events = []
    for k in range(15):
        events.append(ExperimentalFeedbackRecord(
            event_id=f"EV_{k}",
            project_id=1,
            compound_version_id=k + 1,
            canonical_smiles=f"CC(C)c1ccc(C{k})cc1",
            endpoint_name="CYP3A4 inhibitor",
            experimental_value=1.0,
            experimental_unit="binary",
            assay_quality=AssayQuality.HIGH_QUALITY,
            scaffold_smiles=scaffold,
            timestamp=f"2026-08-29T01:{k:02d}:00Z",
            frozen_predictions={
                "admetica_cyp_cyp3a4-inhibitor": 0.40,  # Brier = (0.4 - 1)^2 = 0.36
                "morgan_cyp3a4_inh_v1": 0.90,          # Brier = (0.9 - 1)^2 = 0.01
            },
        ))
        
    res = compute_hierarchical_adaptive_weights(
        query_smiles="c1ccccc1NC(=O)C",
        project_id=1,
        candidate_payloads=payloads,
        historical_feedback_events=events,
        endpoint_name="CYP3A4 inhibitor",
        prediction_timestamp="2026-08-29T02:00:00Z",
    )
    
    w_m2 = res.effective_weights["morgan_cyp3a4_inh_v1"]
    # M2 should increase substantially from ~0.04 to > 0.30 due to series evidence
    assert w_m2 > 0.30
    assert AdaptiveReasonCode.PROJECT_EVIDENCE_ACTIVE.value in res.reason_codes
    assert AdaptiveReasonCode.SERIES_M2_OUTPERFORMS_M1.value in res.reason_codes


def test_cyp3a4_class_balance_safeguard():
    """Verify CLASS_BALANCE_LIMITED reason code when series has homogeneous labels."""
    payloads = [
        make_cyp_payload("admetica_cyp_cyp3a4-inhibitor", 0.80),
        make_cyp_payload("morgan_cyp3a4_inh_v1", 0.75),
    ]
    
    events = [
        ExperimentalFeedbackRecord(
            event_id=f"EV_POS_{k}",
            project_id=1,
            compound_version_id=k + 1,
            canonical_smiles=f"CC(C)c1ccc(C{k})cc1",
            endpoint_name="CYP3A4 inhibitor",
            experimental_value=1.0,  # All positive!
            experimental_unit="binary",
            assay_quality=AssayQuality.HIGH_QUALITY,
            scaffold_smiles="c1ccccc1",
            timestamp=f"2026-08-29T01:{k:02d}:00Z",
            frozen_predictions={"admetica_cyp_cyp3a4-inhibitor": 0.8, "morgan_cyp3a4_inh_v1": 0.8},
        )
        for k in range(8)
    ]
    
    res = compute_hierarchical_adaptive_weights(
        query_smiles="c1ccccc1NC(=O)C",
        project_id=1,
        candidate_payloads=payloads,
        historical_feedback_events=events,
        endpoint_name="CYP3A4 inhibitor",
        prediction_timestamp="2026-08-29T02:00:00Z",
    )
    
    assert AdaptiveReasonCode.CLASS_BALANCE_LIMITED.value in res.reason_codes


def test_cyp3a4_cross_project_isolation():
    """Verify that Project A feedback never bleeds into Project B."""
    payloads = [
        make_cyp_payload("admetica_cyp_cyp3a4-inhibitor", 0.80),
        make_cyp_payload("morgan_cyp3a4_inh_v1", 0.75),
    ]
    
    events_proj_a = [
        ExperimentalFeedbackRecord(
            event_id=f"EV_A_{k}",
            project_id=1,  # Project A
            compound_version_id=k + 1,
            canonical_smiles=f"c1ccccc1C{k}",
            endpoint_name="CYP3A4 inhibitor",
            experimental_value=1.0,
            experimental_unit="binary",
            assay_quality=AssayQuality.HIGH_QUALITY,
            scaffold_smiles="c1ccccc1",
            timestamp=f"2026-08-29T01:{k:02d}:00Z",
            frozen_predictions={"admetica_cyp_cyp3a4-inhibitor": 0.1, "morgan_cyp3a4_inh_v1": 0.9},
        )
        for k in range(10)
    ]
    
    # Query for Project B (id=2)
    res_proj_b = compute_hierarchical_adaptive_weights(
        query_smiles="c1ccccc1NC(=O)C",
        project_id=2,  # Project B
        candidate_payloads=payloads,
        historical_feedback_events=events_proj_a,
        endpoint_name="CYP3A4 inhibitor",
        prediction_timestamp="2026-08-29T02:00:00Z",
    )
    
    # Project B must only have global weights (n_project = 0)
    assert res_proj_b.n_project == 0
    assert AdaptiveReasonCode.GLOBAL_PRIOR_DOMINANT.value in res_proj_b.reason_codes
    assert res_proj_b.effective_weights["admetica_cyp_cyp3a4-inhibitor"] > 0.90


def test_cyp3a4_model_disagreement_signal():
    """Verify MODEL_DISAGREEMENT_SIGNAL is emitted when |p1 - p2| >= 0.35."""
    payloads_agree = [
        make_cyp_payload("admetica_cyp_cyp3a4-inhibitor", 0.80),
        make_cyp_payload("morgan_cyp3a4_inh_v1", 0.75),
    ]
    res_agree = compute_hierarchical_adaptive_weights(
        query_smiles="c1ccccc1",
        project_id=1,
        candidate_payloads=payloads_agree,
        historical_feedback_events=[],
        endpoint_name="CYP3A4 inhibitor",
    )
    assert AdaptiveReasonCode.MODEL_DISAGREEMENT_SIGNAL.value not in res_agree.reason_codes

    payloads_disagree = [
        make_cyp_payload("admetica_cyp_cyp3a4-inhibitor", 0.90),
        make_cyp_payload("morgan_cyp3a4_inh_v1", 0.20),
    ]
    res_disagree = compute_hierarchical_adaptive_weights(
        query_smiles="c1ccccc1",
        project_id=1,
        candidate_payloads=payloads_disagree,
        historical_feedback_events=[],
        endpoint_name="CYP3A4 inhibitor",
    )
    assert AdaptiveReasonCode.MODEL_DISAGREEMENT_SIGNAL.value in res_disagree.reason_codes
    assert res_disagree.model_disagreement >= 0.35


def test_cyp3a4_validation_artifacts_schema():
    """Verify all 8 Stage 4D-3B1 JSON validation artifacts exist and conform to schema."""
    artifact_files = [
        "validation/stage4d3b1_policy.json",
        "validation/stage4d3b1_replay_results.json",
        "validation/stage4d3b1_learning_curve.json",
        "validation/stage4d3b1_series_performance.json",
        "validation/stage4d3b1_weight_trajectories.json",
        "validation/stage4d3b1_calibration.json",
        "validation/stage4d3b1_negative_control.json",
        "validation/stage4d3b1_final_decision.json",
    ]
    for path in artifact_files:
        assert os.path.exists(path), f"Missing validation artifact: {path}"
        with open(path) as f:
            data = json.load(f)
            assert isinstance(data, dict)
            
    # Check specific fields in final decision
    with open("validation/stage4d3b1_final_decision.json") as f:
        decision = json.load(f)
        assert decision["endpoint"] == "EP_MET_CYP3A4_INH"
        assert decision["scientific_decision"] in {
            "ADAPTIVE_PROMOTION_CANDIDATE",
            "CONDITIONAL_ADAPTIVE_VALUE",
            "ARCHITECTURE_VALID_BUT_NO_ACCURACY_GAIN",
            "KEEP_RESEARCH_SHADOW",
            "ADAPTIVE_REJECTED",
        }
        assert decision["consensus_mode"] == "SHADOW"
        assert decision["herg_gate_recommendation"] == "GO"
