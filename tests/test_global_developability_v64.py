import json
from pathlib import Path


def load(name):
    return json.loads((Path("validation") / name).read_text())


def test_v64_pka_and_logd_remain_experimental_only_data_limited():
    pka = load("pka_experimental_dataset_v1.json")
    logd = load("logd74_experimental_dataset_v1.json")
    assert pka["qualified_n"] == 0 and pka["acid_n"] == 0 and pka["base_n"] == 0
    assert logd["qualified_n"] == 0
    assert pka["status"] == logd["status"] == "DATA_LIMITED"


def test_v64_never_substitutes_logp_or_derived_pka_as_experiment():
    logd = load("logd74_experimental_dataset_v1.json")
    assert any(x["reason"] == "LOGP_NOT_LOGD74" for x in logd["rejections"])
    assert any(x["reason"] == "DERIVED_NOT_EXPERIMENTAL" for x in logd["rejections"])
    assert any(x["reason"] == "MODELLED_OR_RULE_DERIVED" for x in load("pka_experimental_dataset_v1.json")["rejections"])


def test_v64_global_baseline_ready_with_explicit_limitations():
    data = load("global_developability_v1_closure_v6_4.json")
    baseline = load("global_developability_v1_baseline.json")
    assert data["global_status"] == "GLOBAL_DEVELOPABILITY_V1_BASELINE_READY_WITH_LIMITATIONS"
    assert baseline["reference_library"]["n"] == 1000
    assert baseline["production_engine"]["status"] == "UNCHANGED"
    assert baseline["production_engine"]["historical_runs_modified"] == 0


def test_v64_rebuild_trigger_and_adapter_contract_are_complete():
    baseline = load("global_developability_v1_baseline.json")
    assert baseline["licensed_data_rebuild_trigger"] == ["NEW_LICENSED_DATA", "canonical qualification", "dataset version increment", "leakage-safe rebuild", "candidate comparison", "independent validation", "production review"]
    assert set(baseline["project_adapter_contract"]) == {"global_prediction", "global_uncertainty", "AD", "model_version", "model_status", "residual-compatible output"}
