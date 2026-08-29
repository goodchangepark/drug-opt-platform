from __future__ import annotations

import json
from pathlib import Path

from backend.endpoint_strategy_registry import get_endpoint_strategy

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    return json.loads((ROOT / "validation" / name).read_text())


def test_search_frozen_and_small():
    protocol = load("stage4e3b_caco2_candidate_search_protocol.json")
    landscape = load("stage4e3b_caco2_candidate_landscape.json")
    prereg = load("stage4e3b_caco2_candidate_preregistration.json")
    assert protocol["status"] == "FROZEN_BEFORE_CANDIDATE_RESULTS"
    assert landscape["candidate_count"] <= 4
    assert prereg["no_new_candidate_search_after_preregistration"] is True


def test_all_candidates_fail_closed_and_not_registered():
    landscape = load("stage4e3b_caco2_candidate_landscape.json")
    strategy = get_endpoint_strategy("permeability_caco2_logpapp")
    registered = set(strategy.primary_model_ids + strategy.shadow_model_ids)
    assert not registered.intersection({c["candidate_id"] for c in landscape["candidates"]})
    assert all(c["decision"] in {"LEGAL_NOT_QUALIFIED", "NO_GO_LICENSE", "NO_GO_ENDPOINT_MISMATCH"} for c in landscape["candidates"])


def test_no_candidate_benchmark_or_expansionrx_tuning():
    overlap = load("stage4e3b_caco2_overlap_audit.json")
    metrics = load("stage4e3b_caco2_candidate_metrics.json")
    assert overlap["benchmark_started"] is False
    assert overlap["expansionrx_used_for_training_or_tuning"] is False
    assert metrics["no_fitting"] is True


def test_final_caco2_closure_preserves_production():
    final = load("stage4e3b_caco2_final_decision.json")
    strategy = get_endpoint_strategy("permeability_caco2_logpapp")
    assert final["closure_decision"] == "CACO2_NO_QUALIFIED_REPLACEMENT_FOUND_CORE_FROZEN"
    assert final["production_changed"] is False
    assert final["no_further_dedicated_caco2_optimization_before_engine_v1_freeze"] is True
    assert strategy.primary_strategy.value == "SINGLE_CORE_MODEL"
    assert strategy.primary_model_ids == ["admetica_caco2"]
