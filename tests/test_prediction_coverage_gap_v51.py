import json
from pathlib import Path

from backend.prediction_coverage import build_prediction_coverage


def test_no_model_inventory_is_deterministic_and_unknown_free():
    view = {"endpoints": [{
        "endpoint_id": "HUMAN_PK_CMAX_ORAL",
        "section": "PK",
        "prediction": {"available": False, "unavailable_reason": "No compatible persisted prediction"},
        "experimental_external_candidates": [{
            "id": 1, "normalized_value": 412.0, "normalized_unit": "ng/mL",
            "species": "HUMAN", "route": "ORAL", "dose": 200.0, "dose_unit": "mg", "regimen": "SINGLE_DOSE",
        }],
        "comparison": {"status": "NO_MATCHING_PREDICTION_MODEL"},
    }]}
    result = build_prediction_coverage(view)
    assert result["unknown_missing_reason"] == 0
    assert result["existing_capability_not_exposed"] == 0
    assert result["no_model_endpoints"][0]["action"] == "CONTEXT_SPECIFIC_MODEL_REQUIRED"


def test_persisted_calculation_is_not_reported_as_unexposed():
    view = {"endpoints": [{
        "endpoint_id": "HUMAN_PK_CMAX_ORAL", "section": "PK",
        "prediction": {"available": False},
        "experimental_external_candidates": [{"normalized_value": 1.0, "normalized_unit": "ng/mL", "species": "HUMAN"}],
        "comparison": {"status": "NO_MATCHING_PREDICTION_MODEL"},
    }]}
    result = build_prediction_coverage(view, calculated_endpoints={"HUMAN_PK_CMAX_ORAL"})
    assert result["existing_capability_not_exposed"] == 1
    assert result["no_model_endpoints"][0]["action"] == "EXISTING_PK_CALCULATION_NOT_EXPOSED"


def test_sunvozertinib_gap_artifact_qa_contract():
    path = Path(__file__).parents[1] / "validation" / "prediction_coverage_gap_closure_v5_1.json"
    data = json.loads(path.read_text())
    assert data["no_model_inventory"]["reported_no_model_groups"] == 33
    assert data["qa_metrics"] == {
        "EXISTING_CAPABILITY_NOT_EXPOSED": 0,
        "PAIRABLE_EVIDENCE_WITHOUT_PAIR": 0,
        "UNKNOWN_MISSING_REASON": 0,
    }
