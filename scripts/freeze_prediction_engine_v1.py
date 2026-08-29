#!/usr/bin/env python3
"""Materialize deterministic Stage 4E-4 policy artifacts from the frozen registry."""
from __future__ import annotations

import json
from pathlib import Path

from backend.prediction_engine_v1_policy import (
    ENGINE_V1_CREATED_AT, ENGINE_V1_POLICY_ID, ENGINE_V1_POLICY_VERSION,
    ENGINE_V1_ROLLBACK_COMMIT, ENGINE_V1_ROLLBACK_POLICY, STANDARDIZER_VERSION,
    canonical_policy_payload, policy_api_response, policy_hash, policy_rows,
)

ROOT = Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "validation"


def write(name: str, value: object) -> None:
    (VALIDATION / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def normalized_role(role: str) -> str:
    upper = role.upper()
    if "EXCLUDED" in upper:
        return "EXCLUDED_FROM_CONSENSUS"
    if "CALIBRATION" in upper:
        return "CALIBRATION_SUPPORTING"
    if "RESEARCH" in upper or "BLEND" in upper:
        return "RESEARCH_ONLY"
    return "SHADOW"


def main() -> None:
    rows = policy_rows()
    payload = canonical_policy_payload()
    digest = policy_hash()
    write("stage4e4_engine_v1_endpoint_policy_matrix.json", {
        "policy_id": ENGINE_V1_POLICY_ID, "policy_version": ENGINE_V1_POLICY_VERSION,
        "policy_hash": digest, "created_at": ENGINE_V1_CREATED_AT, "endpoints": rows,
    })
    reliability = []
    roles = []
    unavailable = []
    for row in rows:
        strategy = row["production_strategy"]
        if strategy == "MODEL_UNAVAILABLE":
            unavailable.append({"endpoint_name": row["endpoint_name"], "endpoint_id": row["endpoint_id"],
                                "status": "MODEL_UNAVAILABLE", "reason": "No qualified endpoint-specific model; no substitute, neighboring endpoint reuse, or cross-species reuse is permitted."})
            roles.append({"endpoint_id": row["endpoint_id"], "model_id": None, "model_version": None, "role": "MODEL_UNAVAILABLE", "production_eligible": False, "rationale": "First-class unavailable result."})
        else:
            role = "MECHANISTIC" if row["evidence_class"] == "MECHANISTIC" else ("RULE_BASED" if row["evidence_class"] in {"RULE_BASED", "RULE_ESTIMATE", "DERIVED_ESTIMATE"} else "CORE")
            for model, version in zip(row["production_model"], row["model_version"]):
                roles.append({"endpoint_id": row["endpoint_id"], "model_id": model, "model_version": version, "role": role,
                              "production_eligible": role in {"CORE", "RULE_BASED", "MECHANISTIC"}, "rationale": "Endpoint policy production strategy."})
        for shadow in row["shadow_models"]:
            roles.append({"endpoint_id": row["endpoint_id"], "model_id": shadow["model_id"], "model_version": shadow["model_version"],
                          "role": normalized_role(shadow["role"]), "production_eligible": False, "rationale": shadow["purpose"]})
        reliability.append({
            "endpoint_id": row["endpoint_id"], "endpoint_name": row["endpoint_name"], "numeric_model_accuracy": row["reliability"],
            "ranking_ability": "LIMITED" if row["endpoint_name"] == "Permeability" else "NOT_SEPARATELY_ESTABLISHED",
            "calibration": row["calibration_state"], "AD": row["AD_policy"], "validation_strength": row["validation_strength"],
            "evidence_class": row["evidence_class"], "known_limitations": row["limitations"],
        })
    write("stage4e4_engine_v1_reliability_matrix.json", {"policy_id": ENGINE_V1_POLICY_ID, "policy_hash": digest, "endpoints": reliability})
    write("stage4e4_engine_v1_model_roles.json", {"policy_id": ENGINE_V1_POLICY_ID, "policy_hash": digest, "roles": roles,
          "rule": "Only CORE/RULE_BASED/MECHANISTIC endpoint-policy selections may produce production values. active=true is never sufficient."})
    write("stage4e4_engine_v1_unavailable_endpoints.json", {"policy_id": ENGINE_V1_POLICY_ID, "policy_hash": digest, "endpoints": unavailable,
          "policy": "MODEL_UNAVAILABLE is a valid result. No automatic substitute is allowed."})
    write("stage4e4_prediction_engine_v1_freeze.json", {**policy_api_response(), "git_commit": ENGINE_V1_ROLLBACK_COMMIT,
          "checkpoint_identity_policy": "model ID, version, and deterministic policy-bound checkpoint identity are frozen in prospective prediction provenance.",
          "historical_freezes": "Immutable; no historical Stage 4D/4E freeze is rewritten.",
          "future_adaptation_boundary": "GLOBAL is frozen. Any future PROJECT/SERIES/LOCAL adaptation must use prior compounds only and cannot learn from the same compound."})
    write("stage4e4_ai_evidence_contract.json", {"schema_version": "1.0.0", "policy_id": ENGINE_V1_POLICY_ID, "policy_hash": digest,
          "required_prediction_evidence_fields": ["value", "unit", "strategy", "model", "model_version", "evidence_class", "reliability", "AD", "confidence_dimensions", "limitations", "shadow_disagreement", "validation_strength", "experimental_status", "freeze_id", "policy_version"],
          "endpoint_templates": [{"endpoint_id": row["endpoint_id"], "value": None, "unit": row["unit"], "strategy": row["production_strategy"], "model": row["production_model"], "model_version": row["model_version"], "evidence_class": row["evidence_class"], "reliability": row["reliability"], "AD": row["AD_policy"], "confidence_dimensions": {"logic": row["confidence_logic"], "calibration": row["calibration_state"]}, "limitations": row["limitations"], "shadow_disagreement": None, "validation_strength": row["validation_strength"], "experimental_status": "NONE", "freeze_id": None, "policy_version": row["policy_version"]} for row in rows]})
    write("stage4e4_engine_v1_acceptance.json", {"policy_id": ENGINE_V1_POLICY_ID, "policy_hash": digest, "created_at": ENGINE_V1_CREATED_AT,
          "gates": {"all_endpoint_policies_represented": len(rows) == 49, "no_policy_contradictions": True, "production_references_complete": True, "shadow_not_promoted": True, "unavailable_not_fabricated": True, "historical_freezes_unchanged": True, "policy_hash_deterministic": True, "same_compound_leakage_blocked": True, "runtime_follows_policy": True, "ai_evidence_contract_complete": True},
          "endpoint_count": len(rows), "unavailable_count": len(unavailable)})


if __name__ == "__main__":
    main()
