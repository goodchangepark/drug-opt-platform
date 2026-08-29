"""Immutable Stage 4E-4 Prediction Engine v1 policy projection.

The Stage 4D registry remains the executable endpoint-policy source.  This
module projects that registry into the Engine-v1 freeze contract and provides
the one canonical content hash used by runtime provenance and release files.
It deliberately contains no model selection or scientific transformation.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any, Dict

from .endpoint_contracts import ENDPOINT_CONTRACTS
from .endpoint_strategy_registry import StrategyType, get_all_strategies

ENGINE_V1_POLICY_ID = "drugopt-prediction-engine-v1"
ENGINE_V1_POLICY_VERSION = "1.0.0"
ENGINE_V1_CREATED_AT = "2026-08-29T00:00:00+00:00"
ENGINE_V1_ROLLBACK_COMMIT = "91fbf21"
ENGINE_V1_ROLLBACK_POLICY = "stage4d4-endpoint-strategy-v1"
STANDARDIZER_VERSION = "CHEM_STANDARDIZER_V1"

_CLOSED = {
    "Permeability": ("LIMITED / LOW-MEDIUM", "INDEPENDENT_COMPLETE_WITH_LIMITATIONS",
        ["ExpansionRx external MAE 0.5695; RMSE 0.7457; Spearman 0.0410", "Limited ranking ability and AD coverage", "Residual training overlap unknown", "Assay heterogeneity"], "CLOSED"),
    "hERG liability": ("LOW / LOW-MEDIUM", "INDEPENDENT_COMPLETE_WITH_LIMITATIONS",
        ["Raw M1: AUROC 0.6669; MCC 0.1844; BAcc 0.5442; sensitivity 0.9755; specificity 0.113 at threshold 0.50", "Raw M1: Brier 0.2745; LogLoss 1.6901; ECE 0.2651", "Platt is calibration research only; it improved historical-holdout Brier/LogLoss/ECE but did not qualify a production threshold", "Prevalence shift, assay heterogeneity, borderline labels, and limited discrimination"], "CLOSED"),
    "HLM intrinsic clearance": ("LOW-MEDIUM", "EXTERNAL_VALIDATION_INSUFFICIENT",
        ["Independent raw benchmark unavailable", "Residual training overlap unknown", "Strict human species isolation"], "CLOSED"),
    "RLM intrinsic clearance": ("LOW-MEDIUM", "EXTERNAL_VALIDATION_INSUFFICIENT",
        ["Independent raw benchmark unavailable", "Residual training overlap unknown", "Strict rat species isolation"], "CLOSED"),
    "MLM intrinsic clearance": ("INSUFFICIENT_EVIDENCE", "EXTERNAL_VALIDATION_INSUFFICIENT",
        ["Independent raw benchmark unavailable", "Residual training overlap unknown", "Strict mouse species isolation"], "CLOSED"),
    "Ionization (pKa)": ("LOW-MEDIUM", "RULE_ESTIMATE",
        ["Site/rule semantics, not validated quantitative ML", "Approximate ±1–2 pKa-unit uncertainty", "Polyprotic and zwitterion limitations"], "CLOSED"),
    "logD pH7.4 derived estimate": ("LOW-MEDIUM", "DERIVED_ESTIMATE",
        ["Depends on cLogP and pKa quality", "Simplified monoprotic ionization and pH 7.4 assumption", "Polyprotic/zwitterion limitations; not validated quantitative ML"], "CLOSED"),
}


def _evidence(strategy: str, validation: str) -> str:
    if strategy == StrategyType.MODEL_UNAVAILABLE.value:
        return "MODEL_UNAVAILABLE"
    if strategy == StrategyType.RULE_ESTIMATE.value:
        return "RULE_ESTIMATE"
    if strategy == StrategyType.DERIVED_ESTIMATE.value:
        return "DERIVED_ESTIMATE"
    if strategy in {StrategyType.RULE_BASED.value, StrategyType.RANK_FUSION.value}:
        return "RULE_BASED"
    if strategy == StrategyType.MECHANISTIC_NO_CONSENSUS.value or validation == "MECHANISTIC":
        return "MECHANISTIC"
    return "VALIDATED_MODEL_PREDICTION" if validation.startswith("INDEPENDENT") else "MODEL_PREDICTION"


def policy_rows() -> list[Dict[str, Any]]:
    """Return every Engine-v1 endpoint row in deterministic endpoint-id order."""
    rows = []
    for name, policy in get_all_strategies().items():
        contract = ENDPOINT_CONTRACTS.get(name)
        strategy = policy.primary_strategy.value
        default_reliability = "UNAVAILABLE" if strategy == "MODEL_UNAVAILABLE" else (
            "MECHANISTIC" if policy.validation_status.value == "MECHANISTIC" else "LIMITED"
        )
        reliability, validation_strength, extra_limits, status = _CLOSED.get(
            name, (default_reliability, policy.validation_status.value, [], "FROZEN")
        )
        shadows = []
        for model_id, version in zip(policy.shadow_model_ids, policy.shadow_model_versions):
            role = policy.non_primary_model_roles.get(model_id, "SHADOW")
            purpose = "calibration support" if "CALIBRATION" in role else "accuracy/complementarity monitoring"
            shadows.append({"model_id": model_id, "model_version": version, "role": role, "purpose": purpose,
                            "changes_production": False})
        if name == "CYP3A4 inhibitor":
            extra_limits.append("Fixed 0.9578/0.0422 blend and dynamic adaptation remain research/shadow only (NO_GO).")
        calibration = policy.calibration_status.value
        rows.append({
            "endpoint_name": name, "endpoint_id": policy.endpoint_id,
            "species": contract.species if contract else "Not applicable / workflow-level",
            "output_type": contract.output_type.value if contract else "WORKFLOW_OUTPUT",
            "unit": contract.canonical_unit if contract else "context-dependent",
            "endpoint_contract_version": policy.endpoint_contract_version,
            "production_strategy": strategy,
            "production_model": list(policy.primary_model_ids), "model_version": list(policy.primary_model_versions),
            "evidence_class": _evidence(strategy, policy.validation_status.value),
            "validation_state": policy.validation_status.value, "validation_strength": validation_strength,
            "reliability": reliability,
            "AD_policy": policy.applicability_policy,
            "confidence_logic": policy.confidence_policy,
            "shadow_models": shadows, "research_models": [s for s in shadows if "RESEARCH" in s["role"]],
            "calibration_state": {"status": calibration, "production_enabled": policy.calibration_production_enabled,
                                  "threshold": policy.decision_threshold},
            "limitations": list(dict.fromkeys([*policy.limitations, *extra_limits])),
            "rollback_target": {"commit": ENGINE_V1_ROLLBACK_COMMIT, "policy_version": ENGINE_V1_ROLLBACK_POLICY,
                                "endpoint_target": policy.rollback_target},
            "policy_version": ENGINE_V1_POLICY_ID + "@" + ENGINE_V1_POLICY_VERSION,
            "Engine_v1_status": status,
            "freeze_semantics": "Immutable prospective prediction provenance; experimental feedback never rewrites this prediction.",
        })
    return sorted(rows, key=lambda row: row["endpoint_id"])


def canonical_policy_payload() -> Dict[str, Any]:
    return {"policy_id": ENGINE_V1_POLICY_ID, "policy_version": ENGINE_V1_POLICY_VERSION,
            "standardizer_version": STANDARDIZER_VERSION, "registry_version": ENGINE_V1_ROLLBACK_POLICY,
            "rollback": {"commit": ENGINE_V1_ROLLBACK_COMMIT, "policy_version": ENGINE_V1_ROLLBACK_POLICY},
            "evidence_hierarchy": ["EXPERIMENTAL", "VALIDATED_MODEL_PREDICTION", "RULE_ESTIMATE", "DERIVED_ESTIMATE", "MODEL_UNAVAILABLE"],
            "endpoints": policy_rows()}


def policy_hash() -> str:
    payload = json.dumps(canonical_policy_payload(), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def policy_api_response() -> Dict[str, Any]:
    payload = canonical_policy_payload()
    return {**payload, "policy_hash": policy_hash(), "created_at": ENGINE_V1_CREATED_AT,
            "freeze_semantics": "Versioned immutable policy; historical freezes retain their original policy provenance."}
