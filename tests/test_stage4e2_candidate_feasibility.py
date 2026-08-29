"""Stage 4E-2 must remain a fail-closed, non-runtime qualification gate."""
from __future__ import annotations

import json
from pathlib import Path

from backend.endpoint_strategy_registry import get_all_strategies


ROOT = Path(__file__).resolve().parents[1]


def load(name: str):
    return json.loads((ROOT / "validation" / name).read_text())


def test_all_candidates_and_datasets_have_final_decisions():
    decisions = load("stage4e2_candidate_decisions.json")
    assert {row["candidate_id"] for row in decisions["models"]} == {
        "MODEL_CARDIOGENAI_HERG", "MODEL_METABOGNN_CLEARANCE",
        "MODEL_PKASOLVER_LITE", "MODEL_PKALEARN_GNN",
    }
    assert {row["dataset_id"] for row in decisions["datasets"]} == {
        "DATA_BIOGEN_PROSPECTIVE", "DATA_EXPANSIONRX", "DATA_LOGD74_1130",
    }


def test_passes_require_license_contract_asset_and_arm64_evidence():
    decisions = load("stage4e2_candidate_decisions.json")
    licenses = {row["candidate_id"]: row for row in load("stage4e2_license_matrix.json")["entries"]}
    contracts = {row["candidate_id"]: row for row in load("stage4e2_endpoint_contract_matrix.json")["entries"]}
    assets = {row["candidate_id"]: row for row in load("stage4e2_asset_manifest.json")["assets"]}
    arm = {row["candidate_id"]: row for row in load("stage4e2_arm64_runtime_matrix.json")["entries"]}
    for row in decisions["models"]:
        if row["decision"].startswith("PASS_TO_STAGE4E3"):
            assert licenses[row["candidate_id"]]["decision"] == "PASS_INTERNAL_RESEARCH"
            assert contracts[row["candidate_id"]]["compatibility_status"].startswith("ENDPOINT_COMPATIBLE")
            assert assets[row["candidate_id"]]["acquired"] is True
            assert arm[row["candidate_id"]]["CPU_inference"] == "PASS"


def test_no_candidate_can_enter_runtime_or_change_production_policy():
    decisions = load("stage4e2_candidate_decisions.json")
    assert decisions["production_changed"] is False
    assert not any(row["may_enter_runtime"] for row in decisions["models"])
    assert load("stage4e2_endpoint_contract_matrix.json")["production_contracts_changed"] is False
    assert load("stage4e2_asset_manifest.json")["assets_are_not_committed"] is True


def test_no_go_and_review_candidates_are_not_registered_models():
    production_model_ids = {
        model_id
        for policy in get_all_strategies().values()
        for model_id in policy.primary_model_ids + policy.shadow_model_ids
    }
    decisions = load("stage4e2_candidate_decisions.json")
    assert not ({row["candidate_id"] for row in decisions["models"]} & production_model_ids)


def test_dataset_access_failures_do_not_claim_independent_benchmark_status():
    datasets = load("stage4e2_dataset_manifest.json")["datasets"]
    assert all(row["usable_n"] == 0 for row in datasets)
    assert all(row["overlap_status"] == "NOT_AUDITABLE_WITHOUT_RAW_STRUCTURES" for row in datasets)
    assert load("stage4e2_overlap_audit.json")["performed"] is False


def test_stage4e3_plan_is_explicitly_blocked_without_passed_assets():
    plan = load("stage4e2_stage4e3_plan.json")
    assert plan["ready"] is False
    assert plan["pass_model_ids"] == []
    assert plan["pass_dataset_ids"] == []
