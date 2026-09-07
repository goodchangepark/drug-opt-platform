import numpy as np

from scripts.evaluate_human_total_iv_cl_v56 import metric


def test_total_iv_metric_is_fold_based_and_positive():
    result = metric(np.asarray([1.0, 10.0]), np.asarray([2.0, 5.0]))
    assert result["n"] == 2
    assert result["aafe"] == 2.0
    assert result["within_2_fold_pct"] == 100.0


def test_total_iv_target_is_not_renamed_as_hepatic():
    import json

    with open("validation/human_total_iv_cl_v5_6.json") as handle:
        artifact = json.load(handle)
    assert artifact["target"]["id"] == "HUMAN_TOTAL_IV_CL"
    assert "not hepatic CL" in artifact["comparison_guard"]
