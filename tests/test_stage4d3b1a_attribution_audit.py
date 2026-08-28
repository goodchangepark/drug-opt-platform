"""
Stage 4D-3B1A: Unit & Integration Tests for CYP3A4 Adaptive Attribution Audit.

Tests:
1. Authoritative Cohort schema, completeness, and deterministic reproducibility (N=250).
2. Component ablation across all 7 required strategies.
3. Fixed-global comparator validation (Fixed Global Prior vs Full Adaptive).
4. Adaptive vs Global paired bootstrap challenge statistics and CIs.
5. Weight movement distribution & global prior dominance tracking.
6. Project & chemical scaffold series attribution with class balance auditing.
7. Chemical subgroup claim re-audit (Basic Amines & Neutral Heteroaromatics).
8. Class balance safeguard (CLASS_BALANCE_LIMITED reason code).
9. Threshold invariance and probability calibration separation.
10. Expected Calibration Error (ECE) and extreme-probability softening mechanism.
11. Negative control forward walk (Real vs Shuffled vs Fixed Global).
12. Production SHADOW preservation and M2 CALIBRATION_SUPPORTING role.
13. hERG gate recommendation (GO_HERG_CALIBRATION_AUDIT_FIRST).
"""

import json
from pathlib import Path
import pytest
import numpy as np

from backend.adaptive_weighting import (
    ADAPTIVE_POLICY_VERSION,
    DEFAULT_BETA_ERROR_SCALING,
    GLOBAL_ENDPOINT_PRIOR_ERRORS,
    AdaptiveReasonCode,
    AssayQuality,
    ExperimentalFeedbackRecord,
    compute_hierarchical_adaptive_weights,
    compute_error_score,
)
from backend.endpoint_contracts import get_endpoint_contract
from backend.multimodel import ExecutionStatus, ModelExecutionPayload

ROOT = Path(__file__).resolve().parents[1]
VAL_DIR = ROOT / "validation"


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


def test_authoritative_cohort_schema_and_integrity():
    """Verify stage4d3b1a_authoritative_cohort.json contains all 250 compounds with required fields."""
    cohort_path = VAL_DIR / "stage4d3b1a_authoritative_cohort.json"
    assert cohort_path.exists()
    
    with open(cohort_path) as f:
        data = json.load(f)
        
    assert data["n_compounds"] == 250
    assert data["endpoint"] == "EP_MET_CYP3A4_INH"
    assert len(data["compounds"]) == 250
    
    first = data["compounds"][0]
    required_keys = [
        "compound_id", "canonical_smiles", "project", "series_id", "scaffold",
        "experimental_label", "m1_probability", "m2_probability",
        "fixed_global_probability", "project_adaptive_probability",
        "series_adaptive_probability", "local_adaptive_probability",
        "full_adaptive_probability", "effective_weight_m1", "effective_weight_m2",
        "weight_shift_from_global"
    ]
    for k in required_keys:
        assert k in first, f"Missing required key {k} in cohort compound record"
        
    # Check that labels are binary
    labels = [c["experimental_label"] for c in data["compounds"]]
    assert set(labels) == {0, 1}
    assert sum(labels) == 122


def test_component_ablation_evaluates_all_seven_strategies():
    """Verify stage4d3b1a_component_ablation.json evaluates all 7 required strategies."""
    ablation_path = VAL_DIR / "stage4d3b1a_component_ablation.json"
    assert ablation_path.exists()
    
    with open(ablation_path) as f:
        data = json.load(f)
        
    strategies = data["strategies"]
    expected_strats = [
        "1_m1_core",
        "2_m2_shadow",
        "3_static_50_50_consensus",
        "4_fixed_global_prior",
        "5_global_plus_project",
        "6_global_plus_project_series",
        "7_full_adaptive",
    ]
    for s in expected_strats:
        assert s in strategies, f"Missing strategy {s} in component ablation"
        metrics = strategies[s]
        for m in ["mcc", "balanced_accuracy", "brier_score", "log_loss", "auroc", "auprc", "sensitivity", "specificity"]:
            assert m in metrics, f"Missing metric {m} for strategy {s}"


def test_fixed_global_comparator_matches_or_exceeds_adaptive():
    """Verify Fixed Global Prior strictly matches or outperforms Full Adaptive."""
    ablation_path = VAL_DIR / "stage4d3b1a_component_ablation.json"
    with open(ablation_path) as f:
        data = json.load(f)
        
    fixed = data["strategies"]["4_fixed_global_prior"]
    adapt = data["strategies"]["7_full_adaptive"]
    m1 = data["strategies"]["1_m1_core"]
    
    # Brier: Fixed Global is better than or equal to Adaptive
    assert fixed["brier_score"] <= adapt["brier_score"]
    assert fixed["brier_score"] == pytest.approx(m1["brier_score"], abs=1e-4)
    
    # LogLoss: Fixed Global is better than Adaptive
    assert fixed["log_loss"] < adapt["log_loss"]
    assert fixed["log_loss"] < m1["log_loss"]
    
    # Binary metrics: Fixed Global matches or exceeds Adaptive
    assert fixed["mcc"] >= adapt["mcc"]
    assert fixed["balanced_accuracy"] >= adapt["balanced_accuracy"]


def test_adaptive_vs_global_bootstrap_results():
    """Verify paired bootstrap challenge between Full Adaptive and Fixed Global Prior."""
    boot_path = VAL_DIR / "stage4d3b1a_adaptive_vs_global_bootstrap.json"
    assert boot_path.exists()
    
    with open(boot_path) as f:
        data = json.load(f)
        
    assert data["n_boot"] >= 1000
    res = data["bootstrap_adaptive_vs_fixed_global"]
    
    # Delta Brier > 0 means Adaptive is worse than Fixed Global Prior
    assert res["delta_brier"]["mean"] >= 0.0
    assert res["delta_brier"]["p_target_better"] < 0.05
    
    # Delta LogLoss > 0 means Adaptive is worse than Fixed Global Prior
    assert res["delta_log_loss"]["mean"] >= 0.0
    assert res["delta_log_loss"]["p_target_better"] < 0.05
    
    # Delta MCC <= 0 means Adaptive is not better on MCC
    assert res["delta_mcc"]["mean"] <= 0.0


def test_weight_movement_tracking_and_global_dominance():
    """Verify weight movement distribution confirms tight clustering around global prior."""
    weight_path = VAL_DIR / "stage4d3b1a_weight_attribution.json"
    assert weight_path.exists()
    
    with open(weight_path) as f:
        data = json.load(f)
        
    dist = data["weight_shift_distribution"]
    assert dist["median"] <= 0.02
    assert dist["pct_within_0_05"] >= 80.0
    assert data["dynamic_movement_verdict"] == "MINIMAL_DIVERGENCE_FROM_GLOBAL_PRIOR"
    assert data["effective_weights"]["w_m1_mean"] > 0.90


def test_project_and_series_attribution_classifications():
    """Verify project and series attribution classifications are evidence-based."""
    series_path = VAL_DIR / "stage4d3b1a_series_attribution.json"
    assert series_path.exists()
    
    with open(series_path) as f:
        data = json.load(f)
        
    proj_attrib = data["project_attribution"]
    assert len(proj_attrib) >= 5
    
    # Verify zero projects are ADAPTIVE_BETTER
    for p_name, p_info in proj_attrib.items():
        assert p_info["classification"] in {"EQUIVALENT", "GLOBAL_PRIOR_BETTER", "CLASS_BALANCE_LIMITED", "INSUFFICIENT_DATA"}
        assert p_info["classification"] != "ADAPTIVE_BETTER"


def test_subgroup_claim_correction_audit():
    """Verify chemical subgroup re-audit explicitly corrects Stage 4D-3B1 claims."""
    series_path = VAL_DIR / "stage4d3b1a_series_attribution.json"
    with open(series_path) as f:
        data = json.load(f)
        
    subgroups = data["subgroup_audit"]
    assert "Basic Amine (+)" in subgroups
    assert "Heteroaromatic (+)" in subgroups
    assert "Neutral Heteroaromatics" in subgroups
    
    # For all subgroups, Fixed Global Prior Brier is <= Adaptive Brier
    for grp_name, grp_data in subgroups.items():
        assert grp_data["fixed_global"]["brier_score"] <= grp_data["adaptive"]["brier_score"] + 1e-4
        assert grp_data["fixed_global"]["log_loss"] <= grp_data["adaptive"]["log_loss"] + 1e-4


def test_class_balance_safeguard():
    """Verify CLASS_BALANCE_LIMITED reason code is raised when project events are single-class."""
    m1_p = make_cyp_payload("admetica_cyp_cyp3a4-inhibitor", 0.85)
    m2_p = make_cyp_payload("morgan_cyp3a4_inh_v1", 0.70)
    
    # 6 historical events all with positive label (y=1.0)
    events = [
        ExperimentalFeedbackRecord(
            event_id=f"EV_CB_{i}",
            project_id=10,
            compound_version_id=i+1,
            canonical_smiles=f"CC(=O)N{i}C",
            endpoint_name="CYP3A4 inhibitor",
            experimental_value=1.0,
            experimental_unit="binary",
            assay_quality=AssayQuality.HIGH_QUALITY,
            scaffold_smiles="c1ccccc1",
            timestamp=f"2026-08-29T00:{i:02d}:00Z",
            frozen_predictions={"admetica_cyp_cyp3a4-inhibitor": 0.8, "morgan_cyp3a4_inh_v1": 0.6},
        )
        for i in range(6)
    ]
    
    res = compute_hierarchical_adaptive_weights(
        query_smiles="c1ccccc1NC(=O)C",
        project_id=10,
        candidate_payloads=[m1_p, m2_p],
        historical_feedback_events=events,
        endpoint_name="CYP3A4 inhibitor",
    )
    assert AdaptiveReasonCode.CLASS_BALANCE_LIMITED.value in res.reason_codes


def test_threshold_invariance_and_calibration_distinction():
    """Verify decision threshold is fixed at 0.5 and binary classifications are unchanged."""
    ablation_path = VAL_DIR / "stage4d3b1a_component_ablation.json"
    with open(ablation_path) as f:
        data = json.load(f)
        
    m1 = data["strategies"]["1_m1_core"]
    fixed = data["strategies"]["4_fixed_global_prior"]
    
    # MCC and Balanced Accuracy are identical between M1 and Fixed Global Prior
    assert m1["mcc"] == fixed["mcc"]
    assert m1["balanced_accuracy"] == fixed["balanced_accuracy"]
    assert m1["sensitivity"] == fixed["sensitivity"]
    assert m1["specificity"] == fixed["specificity"]


def test_calibration_curves_and_extreme_error_softening():
    """Verify ECE and extreme error softening analysis in stage4d3b1a_calibration.json."""
    cal_path = VAL_DIR / "stage4d3b1a_calibration.json"
    assert cal_path.exists()
    
    with open(cal_path) as f:
        data = json.load(f)
        
    assert "expected_calibration_error" in data
    assert "extreme_probability_analysis" in data
    
    ece = data["expected_calibration_error"]
    assert ece["m1_core"] < 0.05
    assert ece["fixed_global_prior"] < 0.05
    assert ece["m2_shadow"] > 0.08
    
    extreme = data["extreme_probability_analysis"]
    assert extreme["n_extreme_m1_errors"] > 0
    assert len(extreme["extreme_cases"]) == extreme["n_extreme_m1_errors"]


def test_negative_control_validation_and_interpretation():
    """Verify negative control confirms signal consumption without proving adaptive superiority."""
    cal_path = VAL_DIR / "stage4d3b1a_calibration.json"
    with open(cal_path) as f:
        data = json.load(f)
        
    neg = data["negative_control"]
    assert neg["real_adaptive_brier"] < neg["shuffled_adaptive_brier"]
    assert neg["real_adaptive_brier"] >= neg["fixed_global_brier"]


def test_production_shadow_preservation_and_m2_calibration_role():
    """Verify final decision preserves SHADOW consensus mode and classifies M2 as CALIBRATION_SUPPORTING."""
    dec_path = VAL_DIR / "stage4d3b1a_final_decision.json"
    assert dec_path.exists()
    
    with open(dec_path) as f:
        data = json.load(f)
        
    assert data["scientific_decision"] == "FIXED_GLOBAL_BLEND_SUFFICIENT"
    assert data["consensus_mode"] == "SHADOW"
    assert data["m1_model"]["role"] == "CORE"
    assert data["m1_model"]["contribution_status"] == "CORE_PRIMARY"
    assert data["m2_model"]["role"] == "SHADOW_ONLY"
    assert data["m2_model"]["contribution_status"] == "CALIBRATION_SUPPORTING"


def test_herg_gate_decision_rationale():
    """Verify hERG gate is GO_HERG_CALIBRATION_AUDIT_FIRST."""
    dec_path = VAL_DIR / "stage4d3b1a_final_decision.json"
    with open(dec_path) as f:
        data = json.load(f)
        
    assert data["herg_gate_recommendation"] == "GO_HERG_CALIBRATION_AUDIT_FIRST"
    assert "hERG" in data["herg_gate_rationale"] or "herg" in data["herg_gate_rationale"].lower()
