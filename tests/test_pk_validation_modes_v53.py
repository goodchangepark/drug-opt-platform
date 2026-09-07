from backend.pk_validation_modes import PKValidationMode, validate_mode_inputs, classify_compound_modes


def test_assisted_mode_allows_explicit_observed_upstream_inputs():
    result = validate_mode_inputs(PKValidationMode.OBSERVED_PARAMETER_ASSISTED, {
        "fu": .1, "HLM_Clint": 2, "CL": 1, "VDss": 1, "F": .5, "ka": 1,
        "pKa": 7, "logD7.4": 2, "provenance": {"fu": "EXPERIMENTAL", "HLM_Clint": "EXPERIMENTAL", "CL": "EXPERIMENTAL", "VDss": "EXPERIMENTAL"},
    })
    assert result["complete"] and result["target_leakage"] is False


def test_full_end_to_end_rejects_observed_target_inputs():
    result = validate_mode_inputs(PKValidationMode.FULL_END_TO_END_PREDICTION, {
        "fu": .1, "HLM_Clint": 2, "CL": 1, "VDss": 1, "F": .5, "ka": 1,
        "pKa": 7, "logD7.4": 2,
        "provenance": {"fu": "EXPERIMENTAL", "HLM_Clint": "EXPERIMENTAL", "CL": "EXPERIMENTAL", "VDss": "EXPERIMENTAL"},
    })
    assert result["complete"] is False
    assert set(result["forbidden_observed_inputs"]) == {"fu", "HLM_Clint", "CL", "VDss"}
    assert result["target_leakage"] is True


def test_missing_parameter_propagates_to_all_modes_without_default():
    result = classify_compound_modes(assisted_inputs={"provenance": {}}, hybrid_inputs={"provenance": {}}, full_inputs={"provenance": {}})
    assert all(not value["complete"] for value in result.values())
    assert all(value["status"] == "END_TO_END_INCOMPLETE" for value in result.values())
