"""Stage 4D-7 pre-experimental strategy review.

This module deliberately does *not* train, calibrate, or select a model at
prediction time.  It turns the independently produced Stage 4D evidence and
the Stage 4D-5 fail-closed promotion policy into a reproducible endpoint
review.  The resulting decision is deterministic for a new compound and is
independent of project measurements or the compound's own later result.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from .endpoint_strategy_registry import StrategyType, get_all_strategies


STAGE4D7_POLICY_VERSION = "stage4d7-preexperimental-review-v1"
STAGE4D7_DECISION = "CURRENT_PRODUCTION_RETAINED"
ROOT = Path(__file__).resolve().parents[1]


def canonical_hash(value: Any) -> str:
    """Stable content hash used to bind reviews to the exact source artifact."""
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _artifact(path: str) -> dict[str, Any]:
    data = json.loads((ROOT / path).read_text())
    return {"path": path, "sha256": canonical_hash(data), "data": data}


EVIDENCE_BY_ENDPOINT: Mapping[str, Mapping[str, Any]] = {
    "Solubility": {
        "candidates": ["CURRENT_SINGLE_CORE", "FIXED_WEIGHT_BLEND", "DOMAIN_AWARE_MODEL_SELECTOR"],
        "artifact_paths": [
            "validation/stage4d3a2_m1_bootstrap.json",
            "validation/stage4d3a2_final_decision.json",
        ],
        "decision_flags": ["CURRENT_PRODUCTION_RETAINED", "SHADOW_RETAINED"],
        "summary": "M1 MAE 0.4159 on N=250; adaptive M1+M2 delta MAE +0.0074, 95% CI [-0.0061, 0.0207], P(better)=0.159.",
        "dataset": "Stage 4D-3A2 authoritative diverse solubility cohort",
        "n": 250,
        "separation": "Existing replay/bootstrap evidence; no project measurements and no ALENIGLIPRON tuning.",
        "best": "CURRENT_SINGLE_CORE",
    },
    "Permeability": {
        "candidates": ["CURRENT_SINGLE_CORE", "FIXED_WEIGHT_BLEND", "STATIC_CONSENSUS"],
        "artifact_paths": ["validation/stage4d2c_bootstrap_comparison.json"],
        "decision_flags": ["CURRENT_PRODUCTION_RETAINED", "DATA_LIMITED_ACCURACY_CEILING", "INSUFFICIENT_EVIDENCE"],
        "summary": "N≈34; apparent consensus MAE benefit has CI crossing zero and inadequate independent qualification evidence.",
        "dataset": "Stage 4D-2C Caco-2 representative cohort",
        "n": 34,
        "separation": "Existing bootstrap comparison; no selector fitted on the cohort.",
        "best": "CURRENT_SINGLE_CORE",
    },
    "CYP3A4 inhibitor": {
        "candidates": ["CURRENT_SINGLE_CORE", "FIXED_WEIGHT_BLEND", "CALIBRATED_FIXED_BLEND"],
        "artifact_paths": [
            "validation/stage4d3b1a_final_decision.json",
            "validation/stage4d3b1a_calibration.json",
        ],
        "decision_flags": ["CURRENT_PRODUCTION_RETAINED", "SHADOW_RETAINED"],
        "summary": "Fixed 0.957828/0.042172 blend is research evidence for probability behavior; dynamic adaptation is inferior to the fixed prior.",
        "dataset": "Stage 4D-3B1A representative CYP3A4 cohort",
        "n": 250,
        "separation": "Existing attribution/calibration analysis; probability and discrimination are kept separate.",
        "best": "CURRENT_SINGLE_CORE",
    },
    "hERG liability": {
        "candidates": ["CURRENT_SINGLE_CORE", "CALIBRATED_SINGLE_CORE", "FIXED_WEIGHT_BLEND"],
        "artifact_paths": [
            "validation/stage4d3b2a_calibration.json",
            "validation/stage4d3b2a_model_metrics.json",
            "validation/stage4d3b2a_final_decision.json",
        ],
        "decision_flags": ["CURRENT_PRODUCTION_RETAINED", "CALIBRATION_RESEARCH_ONLY", "BETTER_SECONDARY_MODEL_REQUIRED"],
        "summary": "Scaffold-aware 546/182 calibration/holdout split improves probability calibration for Platt M1 but does not improve ranking; M2 discrimination remains inadequate.",
        "dataset": "Stage 4D-3B2A hERG audit cohort",
        "n": 728,
        "separation": "Calibration was fitted on 546 and reported on untouched 182 holdout; threshold remains 0.50.",
        "best": "CURRENT_SINGLE_CORE",
    },
}


@dataclass(frozen=True)
class PromotionReview:
    endpoint_name: str
    current_strategy: str
    candidate_strategy: str
    no_leakage: bool
    independent_validation: bool
    endpoint_requirements_configured: bool
    noninferiority_margin_configured: bool
    decision: str
    reasons: tuple[str, ...]

    @property
    def may_promote(self) -> bool:
        return self.decision == "PRODUCTION_CHANGED_VALIDATED"


def evaluate_candidate_for_promotion(
    endpoint_name: str,
    candidate_strategy: str,
    *,
    no_leakage: bool,
    independent_validation: bool,
    endpoint_requirements_configured: bool,
    noninferiority_margin_configured: bool,
) -> PromotionReview:
    """Apply the Stage 4D-5 conjunctive, fail-closed gate to a reviewed candidate.

    This is intentionally stricter than a point-estimate comparison.  It has
    no input for project observations, so same-compound experimental leakage
    cannot make a candidate eligible.
    """
    policy = get_all_strategies()[endpoint_name]
    reasons: list[str] = []
    if not no_leakage:
        reasons.append("NO_LEAKAGE_GATE_FAILED")
    if not independent_validation:
        reasons.append("INDEPENDENT_VALIDATION_GATE_FAILED")
    if not endpoint_requirements_configured:
        reasons.append("ENDPOINT_SAMPLE_REQUIREMENT_UNCONFIGURED_FAIL_CLOSED")
    if not noninferiority_margin_configured:
        reasons.append("NONINFERIORITY_MARGIN_UNCONFIGURED_FAIL_CLOSED")
    if candidate_strategy == "CURRENT_SINGLE_CORE":
        reasons.append("CURRENT_POLICY_IS_ALREADY_ACTIVE")
    decision = "CURRENT_PRODUCTION_RETAINED" if candidate_strategy == "CURRENT_SINGLE_CORE" else "INSUFFICIENT_EVIDENCE"
    if not reasons and candidate_strategy != "CURRENT_SINGLE_CORE":
        # A future stage may extend this only with the remaining Stage 4D-5
        # gates (calibration, subgroup robustness, rollback, authorization).
        decision = "REVIEW_REQUIRED"
    return PromotionReview(
        endpoint_name=endpoint_name,
        current_strategy=policy.primary_strategy.value,
        candidate_strategy=candidate_strategy,
        no_leakage=no_leakage,
        independent_validation=independent_validation,
        endpoint_requirements_configured=endpoint_requirements_configured,
        noninferiority_margin_configured=noninferiority_margin_configured,
        decision=decision,
        reasons=tuple(reasons),
    )


def build_endpoint_accuracy_matrix() -> dict[str, Any]:
    """Return one honest, deterministic pre-experimental decision per policy."""
    rows: list[dict[str, Any]] = []
    for name, policy in sorted(get_all_strategies().items()):
        evidence = dict(EVIDENCE_BY_ENDPOINT.get(name, {}))
        special_strategy = policy.primary_strategy
        if special_strategy == StrategyType.MODEL_UNAVAILABLE:
            decision_flags = ["MODEL_UNAVAILABLE"]
            candidates = []
            best = "MODEL_UNAVAILABLE"
        elif special_strategy == StrategyType.MECHANISTIC_NO_CONSENSUS:
            decision_flags = ["MECHANISTIC_NO_CONSENSUS"]
            candidates = ["METHOD_QUALIFICATION_ONLY"]
            best = "MECHANISTIC_NO_CONSENSUS"
        elif special_strategy == StrategyType.RANK_FUSION:
            decision_flags = ["CURRENT_PRODUCTION_RETAINED"]
            candidates = ["RANK_FUSION"]
            best = "RANK_FUSION"
        elif special_strategy in {StrategyType.RULE_BASED, StrategyType.RULE_ESTIMATE, StrategyType.DERIVED_ESTIMATE}:
            decision_flags = ["RULE_OR_DERIVED_ONLY"]
            candidates = [special_strategy.value]
            best = special_strategy.value
        else:
            decision_flags = evidence.get("decision_flags", ["CURRENT_PRODUCTION_RETAINED", "INSUFFICIENT_EVIDENCE"])
            candidates = evidence.get("candidates", ["CURRENT_SINGLE_CORE"])
            best = evidence.get("best", "CURRENT_SINGLE_CORE")
        artifacts = [_artifact(path) for path in evidence.get("artifact_paths", [])]
        rows.append({
            "endpoint_name": name,
            "endpoint_id": policy.endpoint_id,
            "endpoint_contract_version": policy.endpoint_contract_version,
            "current_production_strategy": policy.primary_strategy.value,
            "current_production_models": list(policy.primary_model_ids),
            "current_model_versions": list(policy.primary_model_versions),
            "available_qualified_alternatives": list(policy.shadow_model_ids),
            "calibration_candidates": [policy.calibration_status.value] if policy.calibration_status.value not in {"NOT_APPLICABLE", "RAW"} else [],
            "consensus_candidate": policy.shadow_strategy.value if policy.shadow_strategy else None,
            "candidate_strategies_tested": candidates,
            "best_validated_strategy": best,
            "production_decision": decision_flags[0],
            "decision_flags": decision_flags,
            "validation_dataset": evidence.get("dataset"),
            "sample_size": evidence.get("n"),
            "dataset_provenance": [item["path"] for item in artifacts],
            "source_artifacts": [{"path": item["path"], "sha256": item["sha256"]} for item in artifacts],
            "training_overlap_risk": "Not resolved to a promotion-grade independent claim unless explicitly stated in source artifact.",
            "chemical_space_coverage": "See source artifact; no new selector/AD threshold was fitted in Stage 4D-7.",
            "current_validation_status": policy.validation_status.value,
            "promotion_eligibility": "FAIL_CLOSED_UNCONFIGURED_STAGE4D5_REQUIREMENTS",
            "known_limitations": list(policy.limitations),
            "evidence_summary": evidence.get("summary", "No endpoint-specific independent optimization evidence available."),
            "data_separation": evidence.get("separation", "No project or same-compound experimental data is used."),
            "policy_version": policy.policy_version,
            "stage4d7_review_version": STAGE4D7_POLICY_VERSION,
        })
    return {
        "artifact": "STAGE4D7_ENDPOINT_ACCURACY_MATRIX",
        "review_version": STAGE4D7_POLICY_VERSION,
        "scope": "GLOBAL_PRE_EXPERIMENTAL_ONLY",
        "project_or_series_adaptation": "NOT_IMPLEMENTED",
        "same_compound_experimental_data_used": False,
        "aleniglipron_used_for_optimization": False,
        "production_policy_changes": [],
        "endpoint_count": len(rows),
        "endpoints": rows,
    }


def build_candidate_results() -> dict[str, Any]:
    matrix = build_endpoint_accuracy_matrix()
    reviews = []
    for row in matrix["endpoints"]:
        for candidate in row["candidate_strategies_tested"]:
            review = evaluate_candidate_for_promotion(
                row["endpoint_name"], candidate,
                no_leakage=True,
                independent_validation=bool(row["source_artifacts"]),
                endpoint_requirements_configured=False,
                noninferiority_margin_configured=False,
            )
            reviews.append({**review.__dict__, "reasons": list(review.reasons), "may_promote": review.may_promote})
    return {
        "artifact": "STAGE4D7_CANDIDATE_STRATEGY_RESULTS",
        "review_version": STAGE4D7_POLICY_VERSION,
        "promotion_policy": "stage4d5-qualification-policy-v1",
        "manual_activation_mode_preserved": True,
        "production_policy_changes": [],
        "reviews": reviews,
    }
