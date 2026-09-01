from backend.prediction_experimental_comparison import (
    compare_prediction_experiment,
    generate_pairs,
    independent_compound_count,
    performance_summary,
)
import json
from pathlib import Path


def pred(value=1.0, endpoint="Solubility", created="2025-01-01T00:00:00+00:00", version=1):
    return {"id": 7, "version_id": version, "endpoint": endpoint, "predicted_value": value,
            "unit": "log10(mol/L)", "created_at": created}


def exp(value=-1.0, endpoint="Solubility", status="DIRECTLY_COMPARABLE", imported="2025-02-01T00:00:00+00:00", version=1):
    return {"id": 9, "compound_version_id": version, "endpoint": endpoint, "raw_value": value,
            "raw_unit": "log10(mol/L)", "normalized_value": value, "normalized_unit": "log10(mol/L)",
            "canonical_endpoint_id": "solubility_aqueous_logs", "comparability_status": status,
            "display": {"normalized_value": value, "normalized_unit": "log10(mol/L)", "comparability_status": status},
            "import_eligible": True, "duplicate_status": "DISTINCT_MEASUREMENT", "imported_at": imported,
            "source_quality_class": "A", "reference_status": "REFERENCE_RESOLVED_SOURCE_RECORD"}


def test_log_scale_pair_error_and_prospective_gate():
    pair = compare_prediction_experiment(pred() | {"model_predictions": {"A": 1.2, "B": 0.8}}, exp())
    assert pair.pair_class == "TRUE_PROSPECTIVE"
    assert pair.comparison_metric_type == "LOG_ABSOLUTE_ERROR"
    assert pair.signed_error == 2.0 and pair.absolute_error == 2.0
    assert pair.model_errors == {"A": 2.2, "B": 1.8}
    assert pair.adaptation_eligibility


def test_relation_bounds_are_not_equal_numeric_errors():
    pair = compare_prediction_experiment(pred(-2), exp(-1) | {"raw_relation": ">"})
    assert pair.comparison_metric_type == "BOUND_INCONSISTENT"
    assert pair.absolute_error is None


def test_related_cyp_and_species_or_endpoint_mismatch_do_not_train():
    pair = compare_prediction_experiment(
        pred(.9, "CYP3A4 inhibitor"),
        exp(2.6, "CYP3A4 IC50") | {"canonical_endpoint_id": "", "comparability_status": "RELATED_NOT_SAME_ENDPOINT"},
    )
    assert pair.pair_class == "TRUE_PROSPECTIVE"
    assert not pair.adaptation_eligibility and pair.absolute_error is None


def test_same_compound_multiple_observations_do_not_inflate_compound_n():
    rows = [exp(-1.0, version=1) | {"id": 1}, exp(-1.2, version=1) | {"id": 2}, exp(-1.1, version=2) | {"id": 3}]
    pairs = generate_pairs([pred(version=1), pred(version=2)], rows)
    assert independent_compound_count(pairs) == 2
    assert performance_summary(pairs)["pair_count"] == 3


def test_missing_freeze_is_not_adaptation_eligible_and_historical_is_blocked():
    missing = compare_prediction_experiment(pred(created=None), exp())
    historical = compare_prediction_experiment(pred(created="2025-03-01T00:00:00+00:00"), exp())
    assert missing.pair_class == "HISTORICAL_VISIBLE" and not missing.adaptation_eligibility
    assert historical.pair_class == "HISTORICAL_VISIBLE" and not historical.adaptation_eligibility


def test_five_drug_comparison_artifact_is_present_and_engine_pinned():
    artifact = json.loads((Path(__file__).parents[1] / "validation" / "prediction_experimental_comparison_v3_3.json").read_text())
    assert len(artifact["drugs"]) == 5
    assert artifact["pair_policy"]["same_compound_leakage"] == "BLOCKED"
    assert artifact["engine_hash"] == "12757ab197b5a70d8ea1754678d9a342ab0b6ea0d82f2896bebb767d686bbdeb"
