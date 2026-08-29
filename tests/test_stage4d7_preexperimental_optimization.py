"""Stage 4D-7: global pre-experimental optimization review safeguards."""

from __future__ import annotations

import json
from pathlib import Path

from backend.endpoint_strategy_registry import StrategyType, get_all_strategies
from backend.database import SessionLocal
from backend.prediction_orchestrator import _build_execution_plan
from backend.preexperimental_optimization import (
    STAGE4D7_POLICY_VERSION,
    build_candidate_results,
    build_endpoint_accuracy_matrix,
    canonical_hash,
    evaluate_candidate_for_promotion,
)


ROOT = Path(__file__).resolve().parents[1]


def artifact(name: str):
    return json.loads((ROOT / "validation" / name).read_text())


def test_1_every_endpoint_has_one_preexperimental_decision():
    matrix = build_endpoint_accuracy_matrix()
    policies = get_all_strategies()
    assert matrix["endpoint_count"] == len(policies) == 49
    assert {row["endpoint_name"] for row in matrix["endpoints"]} == set(policies)
    assert all(row["production_decision"] for row in matrix["endpoints"])


def test_2_global_review_never_uses_project_or_same_compound_evidence():
    matrix = build_endpoint_accuracy_matrix()
    assert matrix["scope"] == "GLOBAL_PRE_EXPERIMENTAL_ONLY"
    assert matrix["project_or_series_adaptation"] == "NOT_IMPLEMENTED"
    assert matrix["same_compound_experimental_data_used"] is False
    assert matrix["aleniglipron_used_for_optimization"] is False


def test_3_unconfigured_stage4d5_gates_fail_closed():
    review = evaluate_candidate_for_promotion(
        "hERG liability", "CALIBRATED_SINGLE_CORE",
        no_leakage=True,
        independent_validation=True,
        endpoint_requirements_configured=False,
        noninferiority_margin_configured=False,
    )
    assert review.may_promote is False
    assert review.decision == "INSUFFICIENT_EVIDENCE"
    assert "ENDPOINT_SAMPLE_REQUIREMENT_UNCONFIGURED_FAIL_CLOSED" in review.reasons
    assert "NONINFERIORITY_MARGIN_UNCONFIGURED_FAIL_CLOSED" in review.reasons


def test_4_no_candidate_can_activate_from_stage4d7_review():
    candidates = build_candidate_results()
    assert candidates["manual_activation_mode_preserved"] is True
    assert candidates["production_policy_changes"] == []
    assert candidates["reviews"]
    assert not any(review["may_promote"] for review in candidates["reviews"])


def test_5_solubility_bootstrap_preserves_m1_decision():
    bootstrap = artifact("stage4d7_bootstrap_results.json")["solubility"]["delta_mae"]
    assert bootstrap["ci_95"][0] < 0 < bootstrap["ci_95"][1]
    row = next(item for item in build_endpoint_accuracy_matrix()["endpoints"] if item["endpoint_name"] == "Solubility")
    assert row["best_validated_strategy"] == "CURRENT_SINGLE_CORE"
    assert "SHADOW_RETAINED" in row["decision_flags"]


def test_6_caco2_is_conservative_when_data_limited():
    row = next(item for item in build_endpoint_accuracy_matrix()["endpoints"] if item["endpoint_name"] == "Permeability")
    assert row["sample_size"] == 34
    assert "DATA_LIMITED_ACCURACY_CEILING" in row["decision_flags"]
    assert row["best_validated_strategy"] == "CURRENT_SINGLE_CORE"


def test_7_cyp3a4_separates_probability_research_from_production_decision():
    calibration = artifact("stage4d7_calibration_results.json")["cyp3a4"]
    assert calibration["ece"]["m1_core"] < calibration["ece"]["fixed_global_prior"]
    assert "RESEARCH_ONLY" in calibration["decision"]
    policy = get_all_strategies()["CYP3A4 inhibitor"]
    assert policy.primary_strategy == StrategyType.SINGLE_CORE_MODEL
    assert policy.decision_threshold == 0.5


def test_8_herg_calibration_and_discrimination_are_separate():
    calibration = artifact("stage4d7_calibration_results.json")["herg"]
    assert calibration["platt"]["brier_score"] < calibration["raw"]["brier_score"]
    assert abs(calibration["platt"]["auroc"] - calibration["raw"]["auroc"]) < 0.01
    policy = get_all_strategies()["hERG liability"]
    assert policy.primary_strategy == StrategyType.SINGLE_CORE_MODEL
    assert policy.calibration_production_enabled is False
    assert policy.decision_threshold == 0.5


def test_9_model_unavailable_and_mechanistic_endpoints_cannot_enter_ml_promotion():
    matrix = {row["endpoint_name"]: row for row in build_endpoint_accuracy_matrix()["endpoints"]}
    assert matrix["P-gp substrate"]["production_decision"] == "MODEL_UNAVAILABLE"
    assert matrix["PK Systemic Clearance"]["production_decision"] == "MECHANISTIC_NO_CONSENSUS"
    assert matrix["Ionization (pKa)"]["production_decision"] == "RULE_OR_DERIVED_ONLY"


def test_10_domain_analysis_forbids_unvalidated_switching():
    domain = artifact("stage4d7_domain_analysis.json")
    assert domain["domain_aware_selector_fitted"] is False
    assert "automatic model switch" in domain["reason"]


def test_11_stage4d7_baseline_preserves_policy_versions_and_freezes():
    baseline = artifact("stage4d7_preexperimental_baseline.json")
    assert baseline["historical_freezes_mutated"] is False
    assert len(baseline["policies"]) == 49
    actual = get_all_strategies()
    assert {row["endpoint_name"]: row["policy_version"] for row in baseline["policies"]} == {
        name: policy.policy_version for name, policy in actual.items()
    }


def test_12_artifacts_are_reproducibly_bound_to_source_content():
    matrix = artifact("stage4d7_endpoint_accuracy_matrix.json")
    row = next(item for item in matrix["endpoints"] if item["endpoint_name"] == "Solubility")
    source = next(item for item in row["source_artifacts"] if item["path"].endswith("stage4d3a2_m1_bootstrap.json"))
    assert source["sha256"] == canonical_hash(json.loads((ROOT / source["path"]).read_text()))
    assert matrix["review_version"] == STAGE4D7_POLICY_VERSION


def test_13_runtime_execution_plan_still_uses_authoritative_production_core():
    policy = get_all_strategies()["Solubility"]
    with SessionLocal() as session:
        plan = _build_execution_plan("Solubility", session, policy)
    assert plan.core_model_key == policy.primary_model_ids[0]
    assert "esol_delaney_v1" in plan.shadow_adapter_ids
    assert "rdkit_gbr_solubility_v1" not in plan.shadow_adapter_ids
