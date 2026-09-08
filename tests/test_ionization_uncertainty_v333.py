from backend.ionization import analyze_ionization


def test_rule_pka_and_logd_expose_honest_uncertainty() -> None:
    result = analyze_ionization("CCN(CC)CC")
    assert result["model_provenance"]["version"] == "1.1.0"
    assert result["uncertainty"]["interval_type"] == "RULE_RANGE_NOT_CALIBRATED"
    assert result["uncertainty"]["do_not_interpret_as_confidence"] is True
    interval = result["physiological_state_7_4"]["mechanistic_uncertainty_interval"]
    assert interval["lower"] <= result["physiological_state_7_4"]["estimated_logd74"] <= interval["upper"]
    assert interval["calibrated"] is False
    assert result["micro_pka"]["site_microstates_at_ph_7_4"]


def test_experimental_logd_is_retained_but_not_substituted_into_prediction() -> None:
    baseline = analyze_ionization("CC(=O)O")
    with_observation = analyze_ionization(
        "CC(=O)O",
        experimental_logd_records=[{"value": -1.23, "ph": 7.4, "source": "locked observation"}],
    )
    state = with_observation["physiological_state_7_4"]
    assert state["experimental_logd74"] == {
        "value": -1.23,
        "ph": 7.4,
        "source": "locked observation",
        "evidence_type": "EXPERIMENTAL",
    }
    assert state["estimated_logd74"] == baseline["physiological_state_7_4"]["estimated_logd74"]


def test_polyprotic_uncertainty_is_wider_than_monoprotic_default() -> None:
    result = analyze_ionization("NCC(=O)O")
    interval = result["physiological_state_7_4"]["mechanistic_uncertainty_interval"]
    assert interval["upper"] - interval["lower"] >= 3.0
    assert len(result["micro_pka"]["site_microstates_at_ph_7_4"]) >= 2
