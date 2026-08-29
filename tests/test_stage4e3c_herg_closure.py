from __future__ import annotations

import json
from pathlib import Path

from backend.endpoint_strategy_registry import get_endpoint_strategy

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    return json.loads((ROOT / "validation" / name).read_text())


def test_exact_raw_policy_and_threshold_remain_unchanged():
    p = load("stage4e3c_herg_protocol.json")
    strategy = get_endpoint_strategy("safety_herg_blocker_prob")
    assert p["production"]["model_id"] == "admetica_safety_herg"
    assert p["production"]["model_version"] == "admetica-d4f7056-herg-chemprop-v2.1"
    assert p["production"]["threshold"] == 0.5
    assert strategy.primary_strategy.value == "SINGLE_CORE_MODEL"
    assert strategy.primary_model_ids == ["admetica_safety_herg"]
    assert strategy.decision_threshold == 0.5


def test_platt_is_separated_from_discrimination_and_not_promoted():
    p = load("stage4e3c_herg_platt_audit.json")
    assert p["fit_split"]["heldout_touched_by_fit"] is False
    assert p["heldout_m1_platt"]["Brier"] < p["heldout_raw_m1"]["Brier"]
    assert abs(p["heldout_m1_platt"]["AUROC"] - p["heldout_raw_m1"]["AUROC"]) < 0.01
    assert p["decision"] == "CALIBRATION_RESEARCH_ONLY"


def test_candidate_search_is_small_preregistered_and_fail_closed():
    landscape = load("stage4e3c_herg_candidate_landscape.json")
    prereg = load("stage4e3c_herg_candidate_preregistration.json")
    assert landscape["candidate_count"] <= 3
    assert landscape["no_candidate_benchmarked"] is True
    assert prereg["no_new_candidate_search_after_preregistration"] is True
    assert prereg["preregistered_candidates"] == []


def test_final_closure_and_no_activation():
    final = load("stage4e3c_herg_final_decision.json")
    assert final["primary_decision"] == "HERG_NO_QUALIFIED_REPLACEMENT_RAW_M1_FROZEN"
    assert final["platt_decision"] == "CALIBRATION_RESEARCH_ONLY"
    assert final["secondary_decision"] == "NO_QUALIFIED_SECONDARY_MODEL"
    assert final["production_changed"] is False
    assert final["engine_v1_status"] == "CLOSED"
    assert final["no_further_dedicated_herg_optimization_before_engine_v1_freeze"] is True
