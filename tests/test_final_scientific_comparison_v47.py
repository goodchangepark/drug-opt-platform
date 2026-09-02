"""v4.7 automatic evidence display and conservative interpretation tests."""
import json
from pathlib import Path

from backend.endpoint_comparison import _blank, _scientific_rows
from backend.scientific_interpretation import interpret_row, policy_report


def test_qualified_external_is_displayable_without_import_and_prediction_is_explicit_when_absent():
    row = _blank("HUMAN_PPB")
    row["experimental_external_candidates"] = [{"id": 9, "origin": "AUTO_QUALIFIED_EXTERNAL", "state": "AUTO_QUALIFIED_EXTERNAL", "raw_endpoint": "Plasma protein binding", "normalized_value": 91.46, "normalized_unit": "% bound", "comparability": "DIRECT"}]
    result = _scientific_rows([row])[0]
    assert result["primary_experimental_display"]["value"] == 91.46
    assert result["prediction"]["available"] is False
    assert result["prediction"]["unavailable_reason"]


def test_interpretation_keeps_value_and_agreement_independent_and_conservative():
    result = interpret_row(prediction_available=True, direct=True, difference_available=True)
    assert result["value_assessment"] == "CONTEXT_DEPENDENT"
    assert result["agreement"] == "NOT_CALIBRATED"
    assert result["value_assessment_color"] == "neutral"
    assert policy_report()["universal_thresholds"] is False


def test_v47_artifacts_define_stable_rows_and_auto_learning_gates():
    root = Path(__file__).resolve().parents[1] / "validation"
    contract = json.loads((root / "final_scientific_comparison_contract_v4_7.json").read_text())
    learning = json.loads((root / "auto_external_learning_eligibility_v4_7.json").read_text())
    assert contract["row_contract"] == "ScientificResultRow"
    assert contract["frontend_scientific_matching"] is False
    assert learning["automatic_display_without_import"] is True
    assert learning["same_compound_post_prediction_target_used_for_validation"] is False
