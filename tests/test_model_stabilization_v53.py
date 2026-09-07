import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_stabilization_baseline_preserves_engine_and_cohort():
    data = json.loads((ROOT / "validation/model_stabilization_baseline_v5_3.json").read_text())
    assert data["frozen_prediction_engine_v1"]["id"] == "drugopt-prediction-engine-v1@1.0.0"
    assert data["frozen_prediction_engine_v1"]["hash"] == "12757ab197b5a70d8ea1754678d9a342ab0b6ea0d82f2896bebb767d686bbdeb"
    assert data["validation_cohorts"]["clinical_pk"]["n"] == 30
    assert data["drugbank_compounds"] == 250
    assert data["production_unchanged"] is True


def test_inventory_has_required_scientific_endpoints_and_triage():
    data = json.loads((ROOT / "validation/model_performance_inventory_v5_3.json").read_text())
    rows = {row["endpoint"]: row for row in data["rows"]}
    for endpoint in ("SOLUBILITY_GENERIC", "CACO2_PAPP_AB", "HUMAN_PPB", "HLM_CLINT", "RLM_CLINT", "MLM_CLINT", "VDSS", "PKA", "LOGD_7_4", "HUMAN_HEPATIC_CL", "HUMAN_TOTAL_CL", "PK_AUC", "PK_CMAX"):
        assert endpoint in rows
        assert rows[endpoint]["triage"]
    assert rows["HUMAN_TOTAL_CL"]["validation_n"] == 30


def test_information_gap_map_does_not_expand_drugbank_by_default():
    data = json.loads((ROOT / "validation/information_gap_map_v5_3.json").read_text())
    assert data["decision"] == "DO_NOT_EXPAND_YET"
    assert data["gaps"]
