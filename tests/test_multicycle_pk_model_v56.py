import json


def test_three_numerical_cycles_are_recorded_without_production_promotion():
    with open("validation/multicycle_pk_model_development_v5_6.json") as handle:
        artifact = json.load(handle)
    assert set(artifact["cycles"]) == {"cycle_1_total_iv_cl", "cycle_2_vdss", "cycle_3_iv_half_life"}
    assert artifact["cycles"]["cycle_1_total_iv_cl"]["stress_n"] == 187
    assert artifact["cycles"]["cycle_3_iv_half_life"]["metrics"]["n"] == 7
    assert artifact["production"] == "UNCHANGED"


def test_total_iv_and_vdss_are_not_claimed_as_oral_or_hepatic_targets():
    with open("validation/multicycle_pk_model_development_v5_6.json") as handle:
        artifact = json.load(handle)
    assert artifact["cycles"]["cycle_1_total_iv_cl"]["decision"] == "CANDIDATE_NOT_PRODUCTION_QUALIFIED"
    assert artifact["cycles"]["cycle_2_vdss"]["decision"] == "RESEARCH_CANDIDATE_ONLY"
    assert artifact["next_data_requirement"].startswith("Independent source")
