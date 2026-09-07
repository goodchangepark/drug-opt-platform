import json
from pathlib import Path


def test_hepatic_hybrid_contract_prevents_assisted_target_leakage():
    data = json.loads(Path("validation/hepatic_hybrid_benchmark_v73.json").read_text())
    assert data["canonical_contract"]["endpoint"] == "HUMAN_HEPATIC_CL"
    assert data["canonical_contract"]["dataset_n"] == 30
    assert data["decision"] == "RESEARCH_ONLY_NO_PRODUCTION_CHANGE"
    assert all(row["target_fitting_scope"] == "train_fold_only" for row in data["fold_provenance"])
