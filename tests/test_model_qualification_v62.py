import json
from pathlib import Path


ARTIFACT = Path("validation/model_qualification_data_closure_v6_2.json")


def test_v62_candidate_freezes_do_not_promote_or_change_engine():
    data = json.loads(ARTIFACT.read_text())
    assert data["production"] == "UNCHANGED"
    assert data["reference_library_n"] == 1000
    assert data["frozen_candidate_status"]["PPB_FU"]["status"] == "RESEARCH_CANDIDATE"
    assert data["frozen_candidate_status"]["TOTAL_IV_CL"]["status"] == "RESEARCH_CANDIDATE"


def test_v62_independent_validation_is_not_claimed_without_source_family():
    data = json.loads(ARTIFACT.read_text())
    for endpoint in ("PPB_FU", "HUMAN_TOTAL_IV_CL", "VDSS"):
        assert data["independent_validation"][endpoint]["qualified_n"] == 0


def test_v62_vdss_strict_semantics_and_n30_guard():
    vdss = json.loads(ARTIFACT.read_text())["vdss_strict_audit"]
    assert vdss["original_v6_1_qualified_n"] == 735
    assert vdss["strict_qualified_n"] == 50
    assert vdss["pksmart_marker_n"] == 750
    assert vdss["ambiguous_n"] == 735
    assert vdss["rejected_n"] == 15
    assert vdss["strict_n30_overlap"] == 0
    assert vdss["strict_semantics"]["unit"] == "L/kg"
    assert vdss["strict_semantics"]["context"] == "Clinical IV Pharmacokinetics"


def test_v62_global_baseline_requires_independent_validation():
    data = json.loads(ARTIFACT.read_text())
    assert data["global_status"] == "MORE_FOUNDATION_REQUIRED"
    assert "independent validation" in data["global_status_reason"]
