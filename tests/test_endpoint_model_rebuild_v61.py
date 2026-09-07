import json


def _artifact():
    with open("validation/endpoint_model_rebuild_v6_1.json") as handle:
        return json.load(handle)


def test_rebuild_keeps_endpoint_targets_and_source_roles_separate():
    data = _artifact()
    assert data["experiments"]["VDSS"]["target"] == "HUMAN_VDSS"
    assert data["experiments"]["TOTAL_IV_CL"]["target"] == "HUMAN_TOTAL_IV_CL"
    assert data["experiments"]["PPB_FU"]["target"] == "HUMAN_PPB"
    assert data["experiments"]["VDSS"]["source_status"] == "WITHIN_SOURCE_VALIDATION_ONLY"


def test_rebuild_records_locked_denominators_and_no_production_change():
    data = _artifact()
    assert data["experiments"]["VDSS"]["locked_validation_n"] == 30
    assert data["experiments"]["TOTAL_IV_CL"]["locked_validation_n"] == 7
    assert data["production"] == "UNCHANGED"
