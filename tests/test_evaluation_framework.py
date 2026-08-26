"""Unit tests for Scientific Evaluation Framework (Stage 4C-2)."""

import pytest

from backend.evaluation import (
    EVALUATION_REGISTRY,
    aggregate_replicates,
    check_data_leakage,
    compute_classification_metrics,
    compute_regression_metrics,
    evaluate_mmp_directional_accuracy,
    get_rdkit_upgrade_readiness_report,
    parse_censored_observation,
    perform_lightning_security_audit,
)


def test_censored_data_parsing():
    c1 = parse_censored_observation(">100")
    assert c1["numeric_value"] == 100.0
    assert c1["operator"] == ">"
    assert c1["is_censored"] is True

    c2 = parse_censored_observation("<= 0.05")
    assert c2["numeric_value"] == 0.05
    assert c2["operator"] == "<="
    assert c2["is_censored"] is True

    c3 = parse_censored_observation(15.5)
    assert c3["numeric_value"] == 15.5
    assert c3["operator"] == "="
    assert c3["is_censored"] is False


def test_replicate_aggregation():
    # Normal replicates
    r1 = aggregate_replicates([10.0, 12.0, 9.5])
    assert r1["flag"] == "NORMAL"
    assert r1["n_replicates"] == 3
    assert 9.0 <= r1["aggregated_value"] <= 11.0

    # High variability replicates (> 10-fold spread)
    r2 = aggregate_replicates([1.0, 50.0, 2.0])
    assert r2["flag"] == "HIGH_EXPERIMENTAL_VARIABILITY"
    assert r2["fold_spread"] == 50.0


def test_mmp_directional_accuracy():
    pairs = [
        {"exp_A": 10.0, "exp_B": 100.0, "pred_A": 12.0, "pred_B": 80.0},   # Correct increase
        {"exp_A": 50.0, "exp_B": 5.0,   "pred_A": 40.0, "pred_B": 6.0},    # Correct decrease
        {"exp_A": 10.0, "exp_B": 100.0, "pred_A": 50.0, "pred_B": 10.0},   # Incorrect direction
    ]

    res = evaluate_mmp_directional_accuracy(pairs, min_delta_fold=1.5)
    assert res["eligible_pair_count"] == 3
    assert res["correct_direction_count"] == 2
    assert res["incorrect_direction_count"] == 1
    assert round(res["directional_accuracy_pct"], 1) == 66.7


def test_regression_metrics():
    y_true = [1.0, 2.0, 3.0, 4.0, 5.0]
    y_pred = [1.1, 1.9, 3.2, 3.8, 5.1]
    res = compute_regression_metrics(y_true, y_pred, scope="TEST")

    assert res["N"] == 5
    assert res["MAE"] < 0.2
    assert res["RMSE"] < 0.2
    assert res["R2"] > 0.95
    assert res["pct_within_2fold"] == 100.0


def test_classification_metrics():
    y_true = [1, 1, 1, 0, 0, 0]
    y_prob = [0.9, 0.8, 0.7, 0.1, 0.2, 0.3]
    res = compute_classification_metrics(y_true, y_prob, threshold=0.5, scope="TEST")

    assert res["N"] == 6
    assert res["balanced_accuracy"] == 1.0
    assert res["mcc"] == 1.0
    assert res["brier_score"] < 0.05


def test_data_leakage_check():
    train = ["CC(=O)Oc1ccccc1C(=O)O", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"]
    test = ["CC(=O)Oc1ccccc1C(=O)O", "c1ccccc1"]
    leak = check_data_leakage(train, test)

    assert leak["status"] == "EVALUATED"
    assert leak["exact_structure_overlap_count"] == 1
    assert leak["exact_structure_overlap_pct"] == 50.0


def test_lightning_security_audit():
    audit = perform_lightning_security_audit()
    assert audit["status"] == "SECURE"
    assert audit["is_safe"] is True
    assert audit["installed_version"] == "2.6.5"


def test_rdkit_upgrade_readiness():
    report = get_rdkit_upgrade_readiness_report()
    assert report["readiness_status"] == "READY_FOR_CANDIDATE_TESTING"
    assert report["golden_gate_summary"]["gate_passed"] is True
