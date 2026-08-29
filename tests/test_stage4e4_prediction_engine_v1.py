"""Stage 4E-4 Engine-v1 freeze contract tests (no model search or fitting)."""
from __future__ import annotations

import json
from pathlib import Path

from backend.endpoint_strategy_registry import get_all_strategies
from backend.prediction_engine_v1_policy import (
    ENGINE_V1_POLICY_ID, ENGINE_V1_ROLLBACK_COMMIT, canonical_policy_payload,
    policy_api_response, policy_hash, policy_rows,
)
from backend.prediction_orchestrator import POLICY_VERSION, _build_execution_plan
from backend.admet import ensure_admet_schema
from backend.database import SessionLocal, engine

ROOT = Path(__file__).resolve().parents[1]
V = ROOT / "validation"


def load(name):
    return json.loads((V / name).read_text())


def test_all_registry_policies_are_represented_once():
    rows = policy_rows()
    assert len(rows) == len(get_all_strategies()) == 49
    assert {row["endpoint_name"] for row in rows} == set(get_all_strategies())
    assert len({row["endpoint_id"] for row in rows}) == 49


def test_policy_hash_is_canonical_and_deterministic():
    assert policy_hash() == policy_hash()
    assert policy_hash() == policy_api_response()["policy_hash"]
    assert canonical_policy_payload()["policy_id"] == ENGINE_V1_POLICY_ID


def test_every_available_endpoint_has_a_production_reference():
    for row in policy_rows():
        if row["production_strategy"] != "MODEL_UNAVAILABLE":
            assert row["production_model"], row["endpoint_name"]
            assert len(row["production_model"]) == len(row["model_version"])


def test_closed_caco2_herg_clearance_and_physchem_states_are_preserved():
    rows = {r["endpoint_name"]: r for r in policy_rows()}
    assert rows["Permeability"]["production_strategy"] == "SINGLE_CORE_MODEL"
    assert rows["Permeability"]["production_model"] == ["admetica_caco2"]
    assert rows["hERG liability"]["calibration_state"]["threshold"] == 0.5
    assert rows["hERG liability"]["calibration_state"]["production_enabled"] is False
    assert rows["HLM intrinsic clearance"]["reliability"] == "LOW-MEDIUM"
    assert rows["RLM intrinsic clearance"]["reliability"] == "LOW-MEDIUM"
    assert rows["MLM intrinsic clearance"]["reliability"] == "INSUFFICIENT_EVIDENCE"
    assert rows["Ionization (pKa)"]["evidence_class"] == "RULE_ESTIMATE"
    assert rows["logD pH7.4 derived estimate"]["evidence_class"] == "DERIVED_ESTIMATE"


def test_shadow_roles_cannot_be_production_and_unavailable_stays_unavailable():
    roles = load("stage4e4_engine_v1_model_roles.json")["roles"]
    assert all(not r["production_eligible"] for r in roles if r["role"] in {"SHADOW", "RESEARCH_ONLY", "CALIBRATION_SUPPORTING", "EXCLUDED_FROM_CONSENSUS"})
    unavailable = load("stage4e4_engine_v1_unavailable_endpoints.json")["endpoints"]
    assert len(unavailable) == 22
    assert {r["endpoint_id"] for r in unavailable}.isdisjoint({"transporter_pgp_inhibitor_prob"})


def test_artifact_contracts_and_ai_evidence_are_complete():
    matrix = load("stage4e4_engine_v1_endpoint_policy_matrix.json")
    reliability = load("stage4e4_engine_v1_reliability_matrix.json")
    ai = load("stage4e4_ai_evidence_contract.json")
    assert len(matrix["endpoints"]) == len(reliability["endpoints"]) == len(ai["endpoint_templates"]) == 49
    assert set(ai["required_prediction_evidence_fields"]) >= {"value", "AD", "freeze_id", "policy_version"}
    assert matrix["policy_hash"] == reliability["policy_hash"] == ai["policy_hash"] == policy_hash()


def test_runtime_uses_engine_v1_version_and_registry_plan():
    assert POLICY_VERSION.startswith(ENGINE_V1_POLICY_ID)
    ensure_admet_schema(engine)
    with SessionLocal() as db:
        plan = _build_execution_plan("hERG liability", db, get_all_strategies()["hERG liability"])
    assert plan.production_strategy == "SINGLE_CORE_MODEL"
    assert plan.decision_threshold == 0.5


def test_freeze_rollback_and_acceptance_contract_are_exact():
    freeze = load("stage4e4_prediction_engine_v1_freeze.json")
    acceptance = load("stage4e4_engine_v1_acceptance.json")
    assert freeze["rollback"]["commit"] == ENGINE_V1_ROLLBACK_COMMIT
    assert freeze["policy_hash"] == policy_hash()
    assert all(acceptance["gates"].values())
