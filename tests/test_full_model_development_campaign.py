"""Scientific governance checks for the restartable full-model campaign."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TERMINAL = {
    "PROMOTED", "RETAIN_EXISTING_MODEL", "CURRENT_DATA_CEILING",
    "MECHANISTIC_ONLY", "MODEL_UNAVAILABLE", "CONTEXT_REQUIRED",
}


def _load(name: str):
    return json.loads((ROOT / "validation" / name).read_text(encoding="utf-8"))


def test_all_scheduled_endpoints_have_evidence_addressed_terminal_decisions():
    state = _load("full_model_development_campaign.json")
    ledger = _load("full_model_optimization_ledger.json")
    tasks = state["tasks"]
    results = ledger["endpoint_results"]
    assert len(tasks) == len(results) == 59
    assert len(state["completed_endpoints"]) == 59
    assert {row["status"] for row in tasks} <= TERMINAL
    assert {row["decision"] for row in results} <= TERMINAL
    assert {row["endpoint"] for row in tasks} == {row["endpoint"] for row in results}
    for result in results:
        assert result["supporting_evidence"]
        for source in result["supporting_evidence"]:
            path = ROOT / source["artifact"]
            assert path.is_file()
            assert hashlib.sha256(path.read_bytes()).hexdigest() == source["sha256"]


def test_no_unqualified_model_was_promoted_or_new_engine_claimed():
    state = _load("full_model_development_campaign.json")
    ledger = _load("full_model_optimization_ledger.json")
    audit = _load("full_model_integrated_admission_audit.json")
    assert ledger["terminal_counts"]["PROMOTED"] == 0
    assert state["accepted_models"] == []
    assert state["engine_baseline"] == "drugopt-prediction-engine-v3@3.3.3"
    assert state["integrated_engine_decision"] == "RETAIN_V3_3_3_NO_QUALIFIED_INTEGRATED_PROMOTION"
    assert audit["summary"]["new_models_promoted"] == 0
    assert audit["summary"]["engine_decision"] == "RETAIN_V3_3_3"


def test_canonical_fu_and_pk_research_limits_are_not_overstated():
    ledger = _load("full_model_optimization_ledger.json")
    rows = {row["endpoint"]: row for row in ledger["endpoint_results"]}
    fu = rows["HUMAN_FU_PLASMA"]
    assert fu["decision"] == "CURRENT_DATA_CEILING"
    assert fu["qualified_n"] == 476
    assert fu["best_metrics"]["aafe"] == 4.760840445074518
    assert fu["best_metrics"]["within_2_fold_pct"] == 50.42016806722689
    assert rows["HUMAN_HEPATIC_CL"]["decision"] == "CURRENT_DATA_CEILING"
    assert rows["HUMAN_HEPATIC_CL"]["best_metrics"]["aafe"] == 2.097108825628647
    assert "ASSISTED" in rows["HUMAN_HEPATIC_CL"]["reason"]
    assert rows["VDSS"]["best_metrics"]["aafe"] == 3.7057366603955706
    assert rows["HUMAN_TOTAL_IV_CL"]["best_metrics"]["aafe"] == 7.341110193181601


def test_only_verified_runtime_routes_are_marked_on_demand_publishable():
    audit = _load("full_model_integrated_admission_audit.json")
    rows = {row["endpoint"]: row for row in audit["rows"]}
    publishable = {
        endpoint for endpoint, row in rows.items() if row["stable_core_publishable_now"]
    }
    assert len(publishable) == 13
    assert {"RLM_CLINT", "MLM_CLINT", "PGP_INHIBITION", "HERG_CLASS", "DILI_LIABILITY"} <= publishable
    assert not rows["SOLUBILITY_GENERIC"]["stable_core_publishable_now"]
    assert rows["SOLUBILITY_GENERIC"]["registry_status"] == "INCOMPLETE_WEIGHTED_ENSEMBLE_BUNDLE"
    assert not rows["HUMAN_PPB"]["stable_core_publishable_now"]
    assert rows["HUMAN_FU_PLASMA"]["scientific_decision"] == "CURRENT_DATA_CEILING"
