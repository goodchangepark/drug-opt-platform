import json


def test_same_resource_holdout_is_not_called_independent_validation():
    with open("validation/pksmart_external_total_iv_cl_v5_6.json") as handle:
        artifact = json.load(handle)
    assert artifact["evaluation_n"] == 187
    assert artifact["production_qualification"] == "NOT_SUFFICIENT_FOR_INDEPENDENT_EXTERNAL_VALIDATION"
    assert artifact["overlap_guard"]["exact_or_morgan_overlap_excluded"] == 115
