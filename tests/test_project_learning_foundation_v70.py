import json
from pathlib import Path


def test_real_project_readiness_artifact_is_fail_closed():
    data = json.loads(Path("validation/project_learning_foundation_v7_0.json").read_text())
    assert [x["project_id"] for x in data["projects"]] == [1, 3, 5]
    assert data["real_pilot"]["trained"] is False
    assert data["global_engine"]["status"] == "UNCHANGED"
    assert "TEST_ONLY" in data["synthetic_fixture_policy"]


def test_project_contract_exposes_both_global_and_project_ad():
    data = json.loads(Path("validation/project_learning_foundation_v7_0.json").read_text())
    assert "GLOBAL_AD" in data["adapter_contract"]
    assert "PROJECT_AD" in data["adapter_contract"]
    assert "recommended_source" in data["adapter_contract"]
