import json
from pathlib import Path


def test_ppb_reconciliation_keeps_percent_bound_and_fu_metrics_separate():
    data = json.loads(Path("validation/ppb_canonical_benchmark_v72.json").read_text())
    assert data["reconciliation"]["verdict"] == "C_METRICS_NOT_DIRECTLY_COMPARABLE"
    assert data["legacy_reproduction"]["target"] == "HUMAN_PPB_PERCENT_BOUND"
    assert data["legacy_reproduction"]["dataset_n"] == 232
    assert data["canonical_contract"]["endpoint"] == "HUMAN_FU_PLASMA"
    assert data["canonical_contract"]["dataset_n"] == 476


def test_ppb_canonical_candidate_is_research_only_and_has_required_subgroups():
    data = json.loads(Path("validation/ppb_canonical_benchmark_v72.json").read_text())
    selected = data["selected_research_candidate"]
    assert selected == "extra_trees_logfu"
    assert data["decision"] == "RESEARCH_ONLY_NO_PRODUCTION_CHANGE"
    for candidate in data["candidates"].values():
        assert candidate["pooled_oof"]["n"] == 476
        assert set(candidate["subgroups"]) == {"fu_lt_0_01", "fu_0_01_to_0_1", "fu_ge_0_1"}
