"""v4.4B public PK benchmark isolation, split, and fail-closed baseline."""
from pathlib import Path
import pytest

from backend.public_pk_benchmark import (
    BENCHMARK_VERSION, DEVELOPMENT, FINAL_EVALUATION, TRACK_B, baseline_rows,
    benchmark_package, coverage, freeze_compound_split, mechanistic_verification,
    records_for_mode, validate_observation,
)


def test_public_benchmark_records_are_context_complete_and_identity_qualified():
    package = benchmark_package()
    assert package["benchmark_version"] == BENCHMARK_VERSION
    assert package["dataset_status"] == "LIMITED_BENCHMARK"
    assert package["observation_count"] >= 20
    assert all(not validate_observation(row) for row in package["observations"])
    assert all(compound["identity_status"] == "EXACT_STRUCTURE_MATCH" for compound in package["compounds"])


def test_compound_level_split_is_deterministic_and_has_no_overlap():
    rows = benchmark_package()["observations"]
    first, second = freeze_compound_split(rows), freeze_compound_split(rows)
    assert first["assignment"] == second["assignment"]
    assert not set(first["development_compounds"]) & set(first["holdout_compounds"])
    assert first["compound_overlap"] == []


def test_development_mode_cannot_read_final_holdout_targets():
    rows = benchmark_package()["observations"]; split = freeze_compound_split(rows)
    with pytest.raises(PermissionError, match="locked"):
        records_for_mode(rows, split, mode=DEVELOPMENT, partition="HOLDOUT")
    assert records_for_mode(rows, split, mode=FINAL_EVALUATION, partition="HOLDOUT")


def test_current_engine_baseline_is_fail_closed_without_target_leakage():
    rows = benchmark_package()["observations"]; split = freeze_compound_split(rows)
    baseline = baseline_rows(rows, split=split, partition=DEVELOPMENT, mode=DEVELOPMENT)
    assert baseline and all(row["prediction_available"] is False for row in baseline)
    assert all(row["status"] == "INSUFFICIENT_INPUT" for row in baseline)
    assert all(row["observed_value"] is not None and row["predicted_value"] is None for row in baseline)


def test_track_b_is_explicit_equation_verification_not_predictive_performance():
    row = mechanistic_verification()[0]
    assert row["track"] == TRACK_B
    assert row["counts_as_predictive_validation"] is False
    assert row["derived_value"] > 0


def test_coverage_is_separate_from_accuracy_and_benchmark_is_not_project_data():
    package = benchmark_package(); profile = coverage(package["observations"])
    assert profile["HUMAN CMAX"]["prediction_available"] == 0
    assert profile["HUMAN CMAX"]["coverage_percent"] == 0.0
    assert package["project_safety"] == {"stored_in_project_database": False, "can_increase_effective_n": False, "can_increase_maturity": False, "can_train_adapter": False}


def test_required_v44b_artifacts_are_present():
    root = Path(__file__).resolve().parents[1] / "validation"
    for name in ("public_pk_benchmark_v1.json", "public_pk_benchmark_sources_v1.json", "public_pk_benchmark_quality_v1.json", "pk_benchmark_review_queue_v4_4b.json", "pk_benchmark_split_v4_4b.json", "pk_benchmark_development_baseline_v4_4b.json", "pk_benchmark_holdout_baseline_v4_4b.json", "pk_benchmark_coverage_v4_4b.json"):
        assert (root / name).is_file()
