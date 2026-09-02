"""v4.2 contracts for PK provenance, numeric error, and Help history."""

from pathlib import Path

from backend.endpoint_comparison import _comparison, _pk_snapshot_values, pk_f_prediction_is_quantitative
from backend.ivive import PKParameterSet
from backend.platform_info import APP_VERSION, version_history
from backend.canonical_endpoints import DIRECT


def test_iv_reference_arm_is_not_an_oral_bioavailability_prediction():
    pset = PKParameterSet(species="Dog", route="IV", f_predicted=100.0)
    assert not pk_f_prediction_is_quantitative(pset)
    assert all(item[0] != "F" for item in _pk_snapshot_values(pset))


def test_oral_f_requires_all_absorption_components_and_preserves_provenance_gate():
    incomplete = PKParameterSet(
        species="Dog", route="PO", f_predicted=100.0,
        provenance_json={"absorption_info": {"fa_value": .8, "fg_value": None, "fh_value": .9}},
    )
    assert not pk_f_prediction_is_quantitative(incomplete)
    assert all(item[0] != "F" for item in _pk_snapshot_values(incomplete))
    complete = PKParameterSet(
        species="Dog", route="PO", f_predicted=72.0,
        provenance_json={"absorption_info": {"fa_value": .8, "fg_value": .9, "fh_value": 1.0}},
    )
    assert pk_f_prediction_is_quantitative(complete)
    assert [item[0] for item in _pk_snapshot_values(complete)] == ["F"]


def test_direct_semantic_match_reports_error_without_accuracy_grade():
    result = _comparison(
        {"available": True, "canonical_endpoint_id": "DOG_PK_F_ORAL", "display_value": 100.0, "unit": "%"},
        [{"comparability": DIRECT, "normalized_value": 51.1, "state": "EXTERNAL_CANDIDATE", "id": 1}],
    )
    assert result["status"] == "DIRECT"
    assert result["error_metric_type"] == "percentage_points"
    assert result["error_value"] == 48.9
    assert result["performance_status"] == "PERFORMANCE_NOT_CALIBRATED"
    assert "accuracy" not in result["performance_status"].lower()


def test_help_history_restores_post_v1_scientific_milestones_without_changing_product_version():
    history = version_history()
    versions = {row["version"] for row in history}
    assert APP_VERSION == "1.0.0"
    assert {"v3.5", "v3.6", "v3.7", "v3.8A", "v3.8B", "v3.9", "v4.0", "v4.1", "v4.2"} <= versions
    assert history[-1]["version"] == "v4.3"
    assert all(row["date"] == row["release_date"] and row["highlights"] == row["improvements"] for row in history)


def test_v42_artifact_contract_files_are_present():
    root = Path(__file__).resolve().parents[1]
    for filename in (
        "pk_prediction_scientific_contract_v4_2.json",
        "sunvozertinib_pk_pair_audit_v4_2.json",
        "prediction_performance_profile_v4_2.json",
        "prediction_run_dedup_audit_v4_2.json",
        "help_version_history_audit_v4_2.json",
    ):
        assert (root / "validation" / filename).is_file()
