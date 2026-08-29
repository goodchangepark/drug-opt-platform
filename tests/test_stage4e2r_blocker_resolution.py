"""Stage 4E-2R may qualify assets, but must never alter production policy."""
from __future__ import annotations

import json
from pathlib import Path

from backend.endpoint_strategy_registry import get_all_strategies

ROOT = Path(__file__).resolve().parents[1]


def load(name: str):
    return json.loads((ROOT / "validation" / name).read_text())


def test_all_original_blockers_have_explicit_outcomes():
    rows = load("stage4e2r_blocker_resolution.json")["records"]
    assert len(rows) == 7
    assert all(row["final_decision"] and row["sources_checked"] for row in rows)


def test_only_licensed_raw_dataset_can_pass_stage4e3_entry():
    gates = load("stage4e2r_stage4e3_entry_gate.json")["assets"]
    passing = [row for row in gates if row["stage4e3_eligible"]]
    assert [row["asset"] for row in passing] == ["DATA_OPENADMET_EXPANSIONRX_CACO2_PAPP_AB"]
    row = passing[0]
    assert all(row[key] for key in ("license_pass", "checkpoint_or_data_pass", "endpoint_pass", "runtime_pass", "overlap_plan_pass", "external_dataset_available"))


def test_expansionrx_intake_preserves_censored_values_and_reproducible_identity():
    row = load("stage4e2r_dataset_candidates.json")["datasets"][0]
    assert row["license"] == "CC-BY-4.0"
    assert row["n"] == 7618
    assert row["acquisition"]["sha256"] == "f674ec74cca1146bc386f832a32d4b8d921d3c312f92cb436cc005901c724a3c"
    assert row["decision"] == "PASS_WITH_EXCLUSIONS"


def test_no_replacement_is_registered_or_eligible_as_runtime_model():
    replacement = load("stage4e2r_replacement_candidates.json")["candidates"]
    policy_models = {model for policy in get_all_strategies().values() for model in policy.primary_model_ids + policy.shadow_model_ids}
    assert not {row["candidate_id"] for row in replacement} & policy_models
    assert load("stage4e2r_stage4e3_plan.json")["eligible_model_ids"] == []


def test_stage4e3_plan_remains_benchmark_only_and_production_is_unchanged():
    plan = load("stage4e2r_stage4e3_plan.json")
    assert plan["ready_state"] == "PARTIAL_READY_FOR_STAGE_4E3"
    assert plan["plans"][0]["purpose"].endswith("no model activation")
    assert load("stage4e2r_stage4e3_entry_gate.json")["production_changed"] is False
