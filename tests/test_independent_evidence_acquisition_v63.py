import json
from pathlib import Path


ROOT = Path("validation")


def test_v63_all_three_tracks_are_source_family_independent_or_fail_closed():
    data = json.loads((ROOT / "independent_evidence_acquisition_v6_3.json").read_text())
    assert len(data["sources_searched"]) >= 6
    for endpoint in ("HUMAN_VDSS", "HUMAN_TOTAL_IV_CL", "HUMAN_FU_PLASMA"):
        assert data["endpoint_datasets"][endpoint]["qualified_n"] == 0
        assert data["frozen_candidate_evaluation"][endpoint if endpoint != "HUMAN_FU_PLASMA" else "PPB_FU"]["external_n"] == 0


def test_v63_lineage_never_counts_unknown_as_independent():
    lineage = json.loads((ROOT / "independent_evidence_acquisition_v6_3.json").read_text())["study_lineage"]
    assert lineage["independent_studies"] == 0
    assert lineage["unknown_lineage_records_accepted"] == 0


def test_v63_dataset_artifacts_have_exact_units_and_hashes():
    expected = {"human_vdss": "L/kg", "human_total_iv_cl": "mL/min/kg", "human_fu_plasma": "fraction unbound"}
    for stem, unit in expected.items():
        data = json.loads((ROOT / f"{stem}_independent_v1.json").read_text())
        assert data["qualified_n"] == 0
        assert data["canonical_unit"] == unit
        assert len(data["dataset_hash"]) == 64


def test_v63_global_baseline_remains_honest():
    data = json.loads((ROOT / "independent_evidence_acquisition_v6_3.json").read_text())
    assert data["global_status"] == "MORE_FOUNDATION_REQUIRED"
    assert data["production"] == "UNCHANGED"
