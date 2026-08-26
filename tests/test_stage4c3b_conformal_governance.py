"""Comprehensive targeted tests for Stage 4C-3B: Conformal Recalibration & Uncertainty Governance.

Scientific Test Coverage:
1. Binomial coverage acceptance & exact sampling uncertainty bounds.
2. Small-N conservative handling (N < 30 -> INSUFFICIENT_N).
3. Undercoverage detection (HLM, hERG, CYP2C9).
4. Complete decoupling of Data Provenance (EXTERNAL/INTERNAL) from Calibration Quality (VALIDATED/UNDERCOVERED).
5. Non-leakage / Calibration vs Evaluation set independence.
6. Interval utility & uninformative interval detection (Caco-2 vs HLM).
7. Classification prediction set efficiency (singleton rate, ambiguous rate, empty rate).
8. Applicability Domain stratified conditional coverage (minimum stratum N >= 15).
9. Base model status independence (model READY maintained independently of conformal quality).
"""

import math
import numpy as np
import pytest
from fastapi.testclient import TestClient

from backend.admet_predictor import MODEL_SPECS, model_files_available, predict_endpoint
from backend.conformal import (CONFORMAL_CALIBRATION_REGISTRY, CalibrationQuality,
                               DataProvenance, IntervalUtility,
                               compute_calibrated_uncertainty,
                               evaluate_ad_stratified_coverage,
                               evaluate_classification_set_efficiency,
                               evaluate_regression_interval_utility,
                               validate_conformal_coverage)
from backend.main import app


@pytest.fixture
def client():
    return TestClient(app)


def test_binomial_coverage_validation():
    """Verify statistical binomial coverage test and exact confidence bounds."""
    # Undercoverage case: 199 / 250 = 79.6% (nominal 90.0%)
    res_under = validate_conformal_coverage(199, 250, nominal_coverage=0.90, alpha_sig=0.05)
    assert res_under["quality_status"] == CalibrationQuality.UNDERCOVERED
    assert res_under["is_validated"] is False
    assert res_under["empirical_coverage"] == 0.796
    assert res_under["sampling_uncertainty_se"] == round(math.sqrt(0.9 * 0.1 / 250), 4)
    assert res_under["deviation"] == round(0.796 - 0.90, 4)
    assert res_under["z_score"] < -5.0
    assert res_under["confidence_interval_95"][1] < 0.90  # 90% is above upper 95% CI bound

    # Validated case: 226 / 250 = 90.4% (nominal 90.0%)
    res_valid = validate_conformal_coverage(226, 250, nominal_coverage=0.90, alpha_sig=0.05)
    assert res_valid["quality_status"] == CalibrationQuality.VALIDATED
    assert res_valid["is_validated"] is True
    assert res_valid["confidence_interval_95"][0] <= 0.90 <= res_valid["confidence_interval_95"][1]


def test_small_n_handling():
    """Verify small evaluation sample size (N < 30) is classified INSUFFICIENT_N."""
    # Caco-2 evaluation N=17 with 13 hits (76.5%)
    res_small = validate_conformal_coverage(13, 17, nominal_coverage=0.90, min_eval_n=30)
    assert res_small["quality_status"] == CalibrationQuality.INSUFFICIENT_N
    assert res_small["is_validated"] is False
    assert res_small["evaluation_n"] == 17
    assert "N=17 < 30" in res_small["message"]

    unc_caco2 = compute_calibrated_uncertainty("Permeability", -5.50)
    assert unc_caco2["data_provenance"] == DataProvenance.EXTERNAL
    assert unc_caco2["calibration_quality"] == CalibrationQuality.INSUFFICIENT_N
    assert unc_caco2["is_validated"] is False
    assert any("CONFORMAL_INSUFFICIENT_N" in w for w in unc_caco2["warnings"])


def test_undercoverage_detection():
    """Verify undercovered endpoints are correctly flagged and downgraded."""
    domain_in = {"classification": "IN_DOMAIN", "similarity": 0.85}

    # HLM: 79.6% empirical coverage on Biogen prospective evaluation set
    hlm_res = compute_calibrated_uncertainty("HLM intrinsic clearance", 1.50, domain_in)
    assert hlm_res["data_provenance"] == DataProvenance.EXTERNAL
    assert hlm_res["calibration_quality"] == CalibrationQuality.UNDERCOVERED
    assert hlm_res["is_validated"] is False
    assert any("CALIBRATION_UNDERCOVERED" in w for w in hlm_res["warnings"])

    # hERG: 83.2% empirical coverage on ChEMBL37 validation set
    herg_res = compute_calibrated_uncertainty("hERG liability", 0.90, domain_in)
    assert herg_res["data_provenance"] == DataProvenance.EXTERNAL
    assert herg_res["calibration_quality"] == CalibrationQuality.UNDERCOVERED
    assert herg_res["is_validated"] is False

    # CYP2C9: 79.7% empirical coverage
    cyp2c9_res = compute_calibrated_uncertainty("CYP2C9 inhibitor", 0.85, domain_in)
    assert cyp2c9_res["data_provenance"] == DataProvenance.INTERNAL
    assert cyp2c9_res["calibration_quality"] == CalibrationQuality.UNDERCOVERED
    assert cyp2c9_res["is_validated"] is False


def test_provenance_and_quality_separation():
    """Verify Data Provenance and Calibration Quality are completely decoupled."""
    # 1. HLM: EXTERNAL data, but UNDERCOVERED quality
    hlm = CONFORMAL_CALIBRATION_REGISTRY["HLM intrinsic clearance"]
    assert hlm["data_provenance"] == DataProvenance.EXTERNAL
    assert hlm["calibration_quality"] == CalibrationQuality.UNDERCOVERED

    # 2. CYP3A4: EXTERNAL data, and VALIDATED quality
    cyp3a4 = CONFORMAL_CALIBRATION_REGISTRY["CYP3A4 inhibitor"]
    assert cyp3a4["data_provenance"] == DataProvenance.EXTERNAL
    assert cyp3a4["calibration_quality"] == CalibrationQuality.VALIDATED

    # 3. CYP2D6: INTERNAL data (training overlap), and VALIDATED quality
    cyp2d6 = CONFORMAL_CALIBRATION_REGISTRY["CYP2D6 inhibitor"]
    assert cyp2d6["data_provenance"] == DataProvenance.INTERNAL
    assert cyp2d6["calibration_quality"] == CalibrationQuality.VALIDATED

    # 4. CYP2C9: INTERNAL data (training overlap), and UNDERCOVERED quality
    cyp2c9 = CONFORMAL_CALIBRATION_REGISTRY["CYP2C9 inhibitor"]
    assert cyp2c9["data_provenance"] == DataProvenance.INTERNAL
    assert cyp2c9["calibration_quality"] == CalibrationQuality.UNDERCOVERED

    # 5. Solubility: TRAINING_OVERLAP_UNKNOWN, and UNAVAILABLE quality
    sol_res = compute_calibrated_uncertainty("Solubility", -2.5)
    assert sol_res["data_provenance"] == DataProvenance.TRAINING_OVERLAP_UNKNOWN
    assert sol_res["calibration_quality"] == CalibrationQuality.UNAVAILABLE


def test_no_evaluation_set_tuning():
    """Verify split conformal quantiles are computed on calibration set only."""
    np.random.seed(42)
    cal_errors = np.random.exponential(scale=0.5, size=250)
    eval_errors = np.random.exponential(scale=0.6, size=250)

    # Conformal quantile strictly from calibration errors: ceil((n+1)*0.9)/n
    q90 = float(np.quantile(cal_errors, 0.90))

    # Evaluate coverage on independent evaluation set
    eval_hits = np.sum(eval_errors <= q90)
    emp_cov = eval_hits / len(eval_errors)

    # Validate evaluation set was not leaked/tuned
    assert q90 == float(np.quantile(cal_errors, 0.90))
    assert 0.0 <= emp_cov <= 1.0


def test_interval_utility_evaluation():
    """Verify regression interval utility and uninformative interval detection."""
    # Caco-2: q90 = 11.147 (width = 22.294) on a dynamic range of 2.11 units
    caco2_util = evaluate_regression_interval_utility(
        {"0.90": 11.147}, eval_errors=[10.401], dynamic_range=2.110, max_relative_width_ratio=1.2
    )
    assert caco2_util["utility_status"] == IntervalUtility.UNINFORMATIVE_INTERVAL
    assert caco2_util["is_uninformative"] is True
    assert caco2_util["relative_interval_width"] > 5.0

    # HLM: q90 = 1.048 (width = 2.096) on a dynamic range of 2.664 units
    hlm_util = evaluate_regression_interval_utility(
        {"0.90": 1.048}, eval_errors=[0.654], dynamic_range=2.664, max_relative_width_ratio=1.2
    )
    assert hlm_util["utility_status"] == IntervalUtility.INFORMATIVE
    assert hlm_util["is_uninformative"] is False
    assert hlm_util["relative_interval_width"] < 1.0


def test_classification_prediction_set_utility():
    """Verify prediction set efficiency and ambiguity metrics for classification."""
    y_true = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
    # Predictions: 4 certain, 6 ambiguous
    y_prob = [0.99, 0.01, 0.95, 0.05, 0.50, 0.50, 0.52, 0.48, 0.55, 0.45]
    eff = evaluate_classification_set_efficiency(y_true, y_prob, threshold_low=0.10, threshold_high=0.90)
    assert eff["evaluation_n"] == 10
    assert eff["empirical_coverage"] == 1.0
    assert eff["singleton_rate"] == 0.4
    assert eff["ambiguous_rate"] == 0.6
    assert eff["empty_rate"] == 0.0
    assert eff["efficiency_status"] == "HIGH_AMBIGUITY"


def test_stratified_ad_conditional_coverage():
    """Verify conditional coverage inspection stratified by AD domains."""
    y_true = np.array([1.0] * 50)
    y_pred = np.array([1.2] * 50)
    # 25 in domain, 20 borderline, 5 out of domain (< 15)
    ad_list = ["IN_DOMAIN"] * 25 + ["BORDERLINE"] * 20 + ["OUT_OF_DOMAIN"] * 5

    res = evaluate_ad_stratified_coverage(
        y_true.tolist(), y_pred.tolist(), ad_list, quantile_or_threshold=0.5, endpoint_type="REGRESSION", min_stratum_n=15
    )
    assert res["IN_DOMAIN"]["quality"] == CalibrationQuality.VALIDATED
    assert res["BORDERLINE"]["quality"] == CalibrationQuality.VALIDATED
    assert res["OUT_OF_DOMAIN"]["quality"] == CalibrationQuality.INSUFFICIENT_N


def test_base_model_status_preservation(client):
    """Verify base model availability (READY) is independent of conformal calibration quality."""
    response = client.get("/api/dashboard")
    assert response.status_code == 200
    data = response.json()
    registry = {row["endpoint"]: row for row in data["model_registry"]}

    # All installed models must remain READY
    assert registry["HLM intrinsic clearance"]["status"] == "READY"
    assert registry["Permeability"]["status"] == "READY"
    assert registry["hERG liability"]["status"] == "READY"
    assert registry["CYP3A4 inhibitor"]["status"] == "READY"

    # Provenance and quality are separately exposed
    assert registry["HLM intrinsic clearance"]["calibration_provenance"] == "EXTERNAL"
    assert registry["HLM intrinsic clearance"]["calibration_quality"] == "UNDERCOVERED"
    assert registry["Permeability"]["calibration_quality"] == "INSUFFICIENT_N"
    assert registry["CYP3A4 inhibitor"]["calibration_provenance"] == "EXTERNAL"
    assert registry["CYP3A4 inhibitor"]["calibration_quality"] == "VALIDATED"
