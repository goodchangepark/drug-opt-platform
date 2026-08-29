"""Stage 4E-3A is a fixed-model, external-data benchmark only."""
from __future__ import annotations

import json
from pathlib import Path

from backend.endpoint_strategy_registry import get_endpoint_strategy


ROOT = Path(__file__).resolve().parents[1]


def load(name: str):
    return json.loads((ROOT / "validation" / name).read_text())


def test_protocol_pins_source_models_units_and_no_consensus():
    protocol = load("stage4e3a_caco2_benchmark_protocol.json")
    assert protocol["revision"] == "6b898ccc43d10d25b230fb09e22a6e30c30022b5"
    assert protocol["sha256"] == "f674ec74cca1146bc386f832a32d4b8d921d3c312f92cb436cc005901c724a3c"
    assert protocol["canonical_unit"] == "log10(cm/s)"
    assert protocol["core"] == "admetica_caco2/admetica-d4f7056-chemprop-v2.1"
    assert protocol["shadow"] == "physchem_caco2_v1/physchem-caco2-v1.0"
    assert protocol["consensus"] == "NO_NUMERIC_CONSENSUS_DEFINED_BY_CURRENT_POLICY"
    assert protocol["no_fitting"] is True


def test_numeric_zero_and_source_censor_provenance_are_distinct():
    flow = load("stage4e3a_caco2_dataset_flow.json")
    assert flow["source_censored_observations"] == 33
    assert flow["non_positive_papp_excluded"] == 2
    assert flow["source_censored_label"] == "SOURCE_CENSORED"
    assert "NON_POSITIVE_PAPP_EXCLUDED" in flow["non_positive_papp_exclusion_label"]
    assert flow["cohort_membership_comparison"]["semantic_correction"] == "METADATA_ONLY_NO_INFERENCE_CHANGE"
    exclusions = load("stage4e3a_caco2_exclusions.json")
    assert exclusions["no_floor_epsilon_imputation_or_replacement"] is True
    assert len(exclusions["source_censored"]) == 33
    assert len(exclusions["non_positive_papp_excluded"]) == 2
    assert all(row["source_row_id"] for row in exclusions["source_censored"])
    assert all(row["category"] == "NON_POSITIVE_PAPP_EXCLUDED" for row in exclusions["non_positive_papp_excluded"])


def test_paired_cache_is_complete_and_exactly_model_identified():
    audit = load("stage4e3a_caco2_prediction_cache_audit.json")
    assert audit["cache_matches_expected_target_set"] is True
    assert audit["expected_unique_target_count"] == 3498
    assert audit["duplicate_structure_hash_count"] == 0
    assert audit["core"]["success"] == 3498
    assert audit["shadow"]["success"] == 3498
    assert audit["core"]["nan_or_inf"] == audit["shadow"]["nan_or_inf"] == 0
    assert audit["production_database_records_created"] is False
    assert audit["core"]["model_id"] == "admetica_caco2"
    assert audit["shadow"]["model_id"] == "physchem_caco2_v1"


def test_bootstrap_is_deterministic_and_core_remains_better():
    bootstrap = load("stage4e3a_caco2_bootstrap.json")
    mae = bootstrap["delta_mae"]
    assert bootstrap["replicates"] == 1000
    assert bootstrap["seed"] == 20260829
    assert mae["point_estimate"] > 0
    assert mae["ci95"][0] > 0
    assert mae["p_shadow_better"] == 0.0
    assert bootstrap["noninferiority_margin"] == "NOT_CONFIGURED"


def test_policy_and_runtime_boundary_stay_unchanged():
    policy = get_endpoint_strategy("permeability_caco2_logpapp")
    assert policy.primary_strategy.value == "SINGLE_CORE_MODEL"
    assert policy.primary_model_ids == ["admetica_caco2"]
    decision = load("stage4e3a_caco2_decision.json")
    assert decision["decision"] == "CURRENT_CORE_CONFIRMED"
    assert decision["production_decision"] == "UNCHANGED"
    assert decision["numeric_consensus"] == "NONE"
    assert decision["promotion"] == "NONE"
