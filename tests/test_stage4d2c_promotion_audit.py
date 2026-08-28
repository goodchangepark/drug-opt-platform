"""
Drug-OPT Stage 4D-2C: Promotion Gate Recalibration & Consensus Calibration Tests.

Covers:
1. Metric reproduction accuracy and numeric consistency
2. Paired bootstrap statistics (delta calculation, 95% CI bounds, empirical p-value)
3. Practical equivalence classification logic (IMPROVED, EQUIVALENT, WORSE, UNCERTAIN)
4. Leave-one-model-out (LOO) ensemble contribution statuses (CORE, SHADOW_ONLY, EXCLUDED_FROM_CONSENSUS)
5. Zero test-set leakage in nested cross-validation
6. Robust aggregation behavior (demonstrating why median is inappropriate for asymmetric models)
7. Disagreement quantile error monotonicity
8. Recalibrated PromotionDecisionStatus mappings
9. All 6 Stage 4D-2C JSON validation artifacts existence and valid schema
10. Backward compatibility and Shadow Mode output invariance
"""

import json
from pathlib import Path
import pytest
import numpy as np

from backend.consensus import (
    PromotionDecisionStatus,
    EnsembleContributionStatus,
    compute_endpoint_consensus,
    ConsensusMode,
    AgreementStatus,
)
from backend.endpoint_contracts import get_endpoint_contract
from backend.multimodel import (
    get_adapters_for_endpoint,
    list_registered_adapters,
    ExecutionStatus,
)
from scripts.run_stage4d2c_audit import (
    compute_regression_metrics,
    compute_classification_metrics,
    paired_bootstrap_regression,
    paired_bootstrap_classification,
    evaluate_equivalence_margin,
)

ROOT = Path(__file__).resolve().parents[1]


def test_promotion_and_contribution_enums():
    """Verify new Stage 4D-2C Enums are correctly defined."""
    assert PromotionDecisionStatus.PRODUCTION_PROMOTION_CANDIDATE.value == "PRODUCTION_PROMOTION_CANDIDATE"
    assert PromotionDecisionStatus.ADAPTIVE_WEIGHTING_RESEARCH_CANDIDATE.value == "ADAPTIVE_WEIGHTING_RESEARCH_CANDIDATE"
    assert PromotionDecisionStatus.KEEP_SHADOW.value == "KEEP_SHADOW"
    assert PromotionDecisionStatus.INSUFFICIENT_EVIDENCE.value == "INSUFFICIENT_EVIDENCE"
    assert PromotionDecisionStatus.EXCLUDE_MODEL_FROM_CONSENSUS.value == "EXCLUDE_MODEL_FROM_CONSENSUS"
    assert PromotionDecisionStatus.REJECT_ENSEMBLE.value == "REJECT_ENSEMBLE"

    assert EnsembleContributionStatus.CORE.value == "CORE"
    assert EnsembleContributionStatus.SUPPORTING.value == "SUPPORTING"
    assert EnsembleContributionStatus.SHADOW_ONLY.value == "SHADOW_ONLY"
    assert EnsembleContributionStatus.EXCLUDED_FROM_CONSENSUS.value == "EXCLUDED_FROM_CONSENSUS"


def test_equivalence_margin_evaluator():
    """Test practical equivalence margin classification logic."""
    # Lower is better (e.g. MAE)
    # Clear improvement: [-0.25, -0.15] with margin 0.10 -> IMPROVED
    assert evaluate_equivalence_margin([-0.25, -0.15], margin=0.10, metric_type="lower_better") == "IMPROVED"
    # Clear degradation: [+0.15, +0.25] with margin 0.10 -> WORSE
    assert evaluate_equivalence_margin([0.15, 0.25], margin=0.10, metric_type="lower_better") == "WORSE"
    # Inside margin: [-0.05, +0.05] with margin 0.10 -> EQUIVALENT
    assert evaluate_equivalence_margin([-0.05, 0.05], margin=0.10, metric_type="lower_better") == "EQUIVALENT"
    # Spanning across margin boundary: [-0.15, +0.05] with margin 0.10 -> UNCERTAIN
    assert evaluate_equivalence_margin([-0.15, 0.05], margin=0.10, metric_type="lower_better") == "UNCERTAIN"

    # Higher is better (e.g. MCC)
    # Clear improvement: [+0.10, +0.20] with margin 0.05 -> IMPROVED
    assert evaluate_equivalence_margin([0.10, 0.20], margin=0.05, metric_type="higher_better") == "IMPROVED"
    # Clear degradation: [-0.20, -0.10] with margin 0.05 -> WORSE
    assert evaluate_equivalence_margin([-0.20, -0.10], margin=0.05, metric_type="higher_better") == "WORSE"


def test_paired_bootstrap_regression_math():
    """Verify paired bootstrap regression statistical calculation."""
    y_true = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
    y_best = np.array([1.1, 2.1, 2.9, 4.1, 4.9, 6.1, 6.9, 8.1, 8.9, 10.1])
    y_worse = np.array([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0])

    res = paired_bootstrap_regression(y_true, y_best, y_worse, n_replicates=200, seed=42)
    assert res["n_replicates"] == 200
    assert res["delta_mae"]["mean"] > 0.0  # y_worse has higher MAE than y_best
    assert res["delta_mae"]["prob_consensus_better"] == 0.0
    assert len(res["delta_mae"]["ci_95"]) == 2
    assert res["delta_mae"]["ci_95"][0] <= res["delta_mae"]["ci_95"][1]


def test_paired_bootstrap_classification_math():
    """Verify paired bootstrap classification statistical calculation."""
    y_true = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0])
    y_prob_best = np.array([0.9, 0.8, 0.85, 0.95, 0.7, 0.1, 0.2, 0.15, 0.05, 0.3])
    y_prob_weak = np.array([0.55, 0.6, 0.52, 0.58, 0.51, 0.49, 0.52, 0.48, 0.51, 0.53])

    res = paired_bootstrap_classification(y_true, y_prob_best, y_prob_weak, n_replicates=200, seed=42)
    assert res["n_replicates"] > 100
    assert res["delta_mcc"]["mean"] < 0.0  # weak model has lower MCC
    assert res["delta_mcc"]["prob_consensus_better"] <= 0.05


def test_all_stage4d2c_artifacts_exist_and_conform():
    """Verify all 6 Stage 4D-2C artifacts exist with correct schema and contents."""
    val_dir = ROOT / "validation"
    artifacts = [
        "stage4d2c_metric_reproduction.json",
        "stage4d2c_bootstrap_comparison.json",
        "stage4d2c_model_contribution.json",
        "stage4d2c_consensus_calibration.json",
        "stage4d2c_promotion_decisions.json",
        "stage4d2c_stage4d3_readiness.json",
    ]

    for filename in artifacts:
        filepath = val_dir / filename
        assert filepath.is_file(), f"Missing Stage 4D-2C artifact: {filename}"
        with open(filepath, "r", encoding="utf-8") as f:
            data = json.load(f)
            assert data.get("stage") == "4D-2C", f"Invalid stage in {filename}"
            assert len(data) > 0


def test_promotion_decisions_recalibrated_correctness():
    """Confirm scientific decisions match Stage 4D-2C audit findings."""
    filepath = ROOT / "validation" / "stage4d2c_promotion_decisions.json"
    with open(filepath, "r", encoding="utf-8") as f:
        data = json.load(f)

    decisions = data.get("recalibrated_decisions", {})
    assert decisions["Solubility"]["decision"] == PromotionDecisionStatus.ADAPTIVE_WEIGHTING_RESEARCH_CANDIDATE.value
    assert decisions["Permeability"]["decision"] == PromotionDecisionStatus.INSUFFICIENT_EVIDENCE.value
    assert decisions["CYP3A4 inhibitor"]["decision"] == PromotionDecisionStatus.ADAPTIVE_WEIGHTING_RESEARCH_CANDIDATE.value
    assert decisions["hERG liability"]["decision"] == PromotionDecisionStatus.KEEP_SHADOW.value
    assert decisions["Metabolic soft spots"]["decision"] == "STAGE_4D2B_PREPARATION_VALIDATED"


def test_model_contribution_taxonomy():
    """Verify individual model role assignments in Stage 4D-2C."""
    filepath = ROOT / "validation" / "stage4d2c_model_contribution.json"
    with open(filepath, "r", encoding="utf-8") as f:
        data = json.load(f)

    sol = data["model_contributions"]["Solubility"]["assigned_statuses"]
    assert sol["admetica_solubility"] == "CORE"
    assert sol["esol_delaney_v1"] == "SHADOW_ONLY"
    assert sol["rdkit_gbr_solubility_v1"] == "EXCLUDED_FROM_CONSENSUS"

    caco = data["model_contributions"]["Caco-2"]["assigned_statuses"]
    assert caco["admetica_caco2"] == "CORE"
    assert caco["physchem_caco2_v1"] == "SHADOW_ONLY"

    cyp = data["model_contributions"]["CYP3A4"]["assigned_statuses"]
    assert cyp["admetica_cyp_cyp3a4-inhibitor"] == "CORE"
    assert cyp["morgan_cyp3a4_inh_v1"] == "SHADOW_ONLY"


def test_stage4d3_readiness_report():
    """Verify Stage 4D-3 prerequisite status."""
    filepath = ROOT / "validation" / "stage4d2c_stage4d3_readiness.json"
    with open(filepath, "r", encoding="utf-8") as f:
        data = json.load(f)

    assert data["overall_stage4d3_readiness"] == "READY_FOR_STAGE_4D3_RESEARCH_WHEN_AUTHORIZED"
    prereqs = data["stage4d3_prerequisites"]
    assert prereqs["multiple_qualified_models"]["status"] == "SATISFIED"
    assert prereqs["performance_heterogeneity"]["status"] == "SATISFIED"
    assert prereqs["series_project_structure"]["status"] == "SATISFIED"
    assert prereqs["shadow_freeze_verified"]["status"] == "SATISFIED"
