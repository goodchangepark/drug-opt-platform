import json
from pathlib import Path

from backend.pk_validation_modes import endpoint_completeness, validate_mode_inputs


def test_hybrid_does_not_relabel_observed_upstream_values():
    result = validate_mode_inputs("HYBRID_PREDICTION", {
        "fu": .1, "HLM_Clint": 2, "VDss": 1, "provenance": {
            "fu": "EXPERIMENTAL", "HLM_Clint": "EXPERIMENTAL", "VDss": "EXPERIMENTAL"
        }
    })
    assert result["complete"] is False
    assert set(result["observed_inputs"]) == {"fu", "HLM_Clint", "VDss"}


def test_endpoint_completeness_is_not_all_or_nothing():
    result = endpoint_completeness(mode="OBSERVED_PARAMETER_ASSISTED", outputs={
        "AUC": {"value": 1, "provenance": "MODEL_PREDICTED"}, "Tmax": None
    })
    assert result["endpoints"]["AUC"]["complete"] is True
    assert result["endpoints"]["Tmax"]["complete"] is False


def test_frozen_baseline_preserves_strict_zero_hybrid_full():
    d = json.loads(Path("validation/end_to_end_pk_baseline_v5_3.json").read_text())
    assert d["modes"]["HYBRID_PREDICTION"]["complete_n"] == 0
    assert d["modes"]["FULL_END_TO_END_PREDICTION"]["complete_n"] == 0
