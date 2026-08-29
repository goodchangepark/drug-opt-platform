#!/usr/bin/env python3
"""Generate deterministic Stage 4D-5 governance artifacts from runtime policy."""

from __future__ import annotations

import json
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from backend.endpoint_strategy_registry import StrategyType
from backend.production_qualification import (
    ACTIVATION_MODE,
    CANDIDATE_TRACKS,
    DEFAULT_DRIFT_POLICY,
    DRIFT_POLICY_VERSION,
    EvidenceTiming,
    EligibilityStatus,
    ExperimentalQualificationResult,
    FrozenModelIdentity,
    ManualPromotionAuthorization,
    MinimumSampleRequirement,
    MonitoringScope,
    NonInferiorityPolicy,
    PerformanceObservation,
    PromotionGateInput,
    ProspectivePredictionFreeze,
    QUALIFICATION_POLICY_REGISTRY,
    QUALIFICATION_POLICY_VERSION,
    QualificationDecision,
    QualificationKind,
    QualificationLifecycle,
    QualificationLifecycleService,
    QualificationPolicy,
    QualificationRecord,
    RollbackReason,
    STATE_MACHINE_VERSION,
    VALIDATION_EVIDENCE_HIERARCHY,
    ValidationType,
    AssayEvidenceQuality,
    CandidateSpecification,
    detect_drift,
    evaluate_experimental_eligibility,
    evaluate_promotion_gate,
    get_production_baseline,
    get_strategy_cards,
    rebuild_performance_snapshot,
    validate_qualification_registry,
)


VALIDATION = ROOT / "validation"


def write_json(name: str, payload: dict | list) -> None:
    (VALIDATION / name).write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=False) + "\n",
        encoding="utf-8",
    )


def build_qualification_policy_artifact() -> dict:
    return {
        "artifact": "stage4d5_qualification_policy",
        "policy_version": QUALIFICATION_POLICY_VERSION,
        "state_machine_version": STATE_MACHINE_VERSION,
        "activation_mode": ACTIVATION_MODE,
        "automatic_shadow_activation": False,
        "automatic_retraining": False,
        "validation_type_taxonomy": [item.value for item in ValidationType],
        "validation_evidence_hierarchy": [
            {"validation_type": item.value, "evidence_rank": VALIDATION_EVIDENCE_HIERARCHY[item]}
            for item in ValidationType
        ],
        "legal_transitions": {
            "RESEARCH_ONLY": ["SHADOW", "RETIRED"],
            "SHADOW": ["VALIDATED", "RETIRED"],
            "VALIDATED": ["SHADOW", "PRODUCTION_CANDIDATE", "RETIRED"],
            "PRODUCTION_CANDIDATE": ["VALIDATED", "ACTIVE", "RETIRED"],
            "ACTIVE": ["RETIRED", "ROLLED_BACK"],
            "RETIRED": [],
            "ROLLED_BACK": ["RETIRED"],
        },
        "direct_shadow_to_active_allowed": False,
        "qualification_policy_count": len(QUALIFICATION_POLICY_REGISTRY),
        "policies": [policy.to_dict() for _, policy in sorted(QUALIFICATION_POLICY_REGISTRY.items())],
        "candidate_tracks": CANDIDATE_TRACKS,
        "registry_violations": validate_qualification_registry(),
    }


def build_promotion_gate_artifact() -> dict:
    return {
        "artifact": "stage4d5_promotion_gates",
        "policy_version": QUALIFICATION_POLICY_VERSION,
        "decision": "NO_CURRENT_CANDIDATE_IS_AUTOMATICALLY_QUALIFIED",
        "activation_mode": ACTIVATION_MODE,
        "required_conjunctive_gates": [
            "endpoint_contract_compatible",
            "model_identity_frozen",
            "no_leakage",
            "independent_validation",
            "endpoint_specific_sample_and_class_balance_sufficient",
            "primary_metric_noninferior",
            "meaningful_or_equivalent_benefit",
            "calibration_acceptable",
            "subgroup_robust",
            "no_unacceptable_safety_tradeoff",
            "rollback_target_exists",
            "artifacts_reproducible",
            "manual_activation_only",
        ],
        "one_metric_promotion_allowed": False,
        "minimum_sample_policy": (
            "Endpoint-specific requirements are intentionally unconfigured until approved with scientific provenance. "
            "An unconfigured requirement fails closed as INSUFFICIENT_EVIDENCE."
        ),
        "noninferiority_policy": (
            "Endpoint-specific margins are versioned and fail closed until approved. Equivalent accuracy may only be "
            "accepted with documented calibration, robustness, applicability, cost, or provenance benefit."
        ),
        "candidate_freeze": [
            "models", "model_versions", "checkpoint_hashes", "weights", "decision_threshold",
            "calibration", "policy_version", "standardizer_version", "endpoint_contract_version",
        ],
        "candidate_change_rule": "ANY_CHANGE_REQUIRES_NEW_CANDIDATE_ID_AND_VERSION",
    }


def build_drift_artifact() -> dict:
    return {
        "artifact": "stage4d5_drift_policy",
        "policy_version": DRIFT_POLICY_VERSION,
        "warning_states": [
            "PERFORMANCE_DRIFT", "CALIBRATION_DRIFT", "DOMAIN_SHIFT", "PRIOR_SHIFT",
            "INSUFFICIENT_DATA", "MODEL_VERSION_CHANGED", "ENDPOINT_MISMATCH",
        ],
        "deterministic_review_thresholds": {
            key: value for key, value in DEFAULT_DRIFT_POLICY.__dict__.items()
        },
        "response": "REVIEW_REQUIRED",
        "automatic_actions": [],
        "explicitly_forbidden": [
            "AUTO_RETRAIN", "THRESHOLD_CHANGE", "ENSEMBLE_WEIGHT_CHANGE",
            "CALIBRATION_ACTIVATION", "PRODUCTION_MODEL_SWITCH",
        ],
        "monitoring_scopes": [item.value for item in MonitoringScope],
        "regression_metrics": ["n", "mae", "rmse", "bias", "spearman", "within_2_fold", "within_3_fold"],
        "classification_metrics": [
            "n", "mcc", "balanced_accuracy", "auroc", "auprc", "brier",
            "log_loss", "sensitivity", "specificity", "ece",
        ],
    }


def build_dry_run() -> dict:
    """TEST-only full lifecycle. No runtime database or production policy is touched."""
    start = datetime(2026, 8, 29, 0, 0, tzinfo=timezone.utc)
    active_model = FrozenModelIdentity("test_herg_active_m1", "test-v1", "sha256:test-active-checkpoint")
    candidate_model = FrozenModelIdentity("test_herg_platt_m1", "test-v1", "sha256:test-candidate-checkpoint")
    active_spec = CandidateSpecification(
        candidate_id="TEST-ACTIVE-HERG-V1", endpoint_id="safety_herg_blocker_prob",
        endpoint_contract_version="1.0.0", candidate_strategy=StrategyType.SINGLE_CORE_MODEL,
        models=(active_model,), policy_version="TEST-POLICY-ACTIVE-V1", candidate_version="TEST-ACTIVE-V1",
        decision_threshold=0.5, calibration="RAW",
    )
    candidate_spec = CandidateSpecification(
        candidate_id="TEST-CANDIDATE-HERG-PLATT-V1", endpoint_id="safety_herg_blocker_prob",
        endpoint_contract_version="1.0.0", candidate_strategy=StrategyType.SINGLE_CORE_WITH_CALIBRATION,
        models=(candidate_model,), policy_version="TEST-POLICY-CANDIDATE-V1", candidate_version="TEST-CANDIDATE-V1",
        decision_threshold=0.5, calibration="PLATT_TEST_ONLY",
    )
    labels = [1, 1, 1, 1, 0, 0, 0, 0]
    active_probabilities = [0.70, 0.60, 0.55, 0.40, 0.60, 0.45, 0.35, 0.20]
    candidate_probabilities = [0.90, 0.85, 0.80, 0.75, 0.25, 0.20, 0.15, 0.10]
    observations: list[PerformanceObservation] = []
    freezes: list[dict] = []
    results: list[dict] = []
    for model_role, spec, probabilities in (
        ("ACTIVE", active_spec, active_probabilities),
        ("SHADOW", candidate_spec, candidate_probabilities),
    ):
        for index, (probability, label) in enumerate(zip(probabilities, labels), start=1):
            frozen_at = start + timedelta(minutes=index)
            frozen = ProspectivePredictionFreeze(
                frozen_prediction_id=f"TEST-FREEZE-{model_role}-{index:02d}",
                compound_version_id=f"TEST-COMPOUND-V{index:02d}",
                endpoint_id=spec.endpoint_id,
                endpoint_contract_version=spec.endpoint_contract_version,
                candidate_id=spec.candidate_id,
                candidate_specification_hash=spec.specification_hash,
                strategy=spec.candidate_strategy,
                models=spec.models,
                prediction_value=probability,
                probability=probability,
                unit="probability",
                frozen_at=frozen_at,
                policy_version=spec.policy_version,
                standardizer_version=spec.standardizer_version,
                project_id="TEST-PROJECT",
                chemical_series_id="TEST-SERIES-A" if index <= 4 else "TEST-SERIES-B",
                applicability_domain="IN_DOMAIN",
                provenance=(("data_scope", "TEST_ONLY_SYNTHETIC"),),
            )
            result = ExperimentalQualificationResult(
                experimental_result_id=f"TEST-RESULT-{model_role}-{index:02d}",
                frozen_prediction_id=frozen.frozen_prediction_id,
                endpoint_id=spec.endpoint_id,
                endpoint_contract_version=spec.endpoint_contract_version,
                experimental_value=float(label),
                unit="binary_label",
                assay_type="Automated patch clamp hERG inhibition label",
                species="Human",
                experiment_date="2026-08-30",
                result_available_at=frozen_at + timedelta(days=1),
                source="TEST_ONLY_SYNTHETIC",
                quality=AssayEvidenceQuality.QUALIFICATION_GRADE,
                protocol_metadata=(("fixture", "stage4d5_dry_run"),),
                recorded_at=frozen_at + timedelta(days=1, minutes=1),
            )
            eligibility = evaluate_experimental_eligibility(frozen, result)
            assert eligibility.status == EligibilityStatus.QUALIFICATION_ELIGIBLE
            observations.append(PerformanceObservation(
                frozen_prediction_id=frozen.frozen_prediction_id,
                endpoint_id=frozen.endpoint_id,
                candidate_id=spec.candidate_id,
                model_version_key=spec.models[0].model_version,
                policy_version=spec.policy_version,
                project_id=frozen.project_id,
                chemical_series_id=frozen.chemical_series_id,
                predicted_value=probability,
                experimental_value=float(label),
                probability=probability,
                eligibility=eligibility.status,
                applicability_domain=frozen.applicability_domain,
            ))
            freezes.append({"role": model_role, **frozen.to_dict()})
            results.append({"role": model_role, **result.to_dict(), "eligibility": eligibility.to_dict()})

    active_snapshot = rebuild_performance_snapshot(
        observations, endpoint_id=active_spec.endpoint_id, candidate_id=active_spec.candidate_id,
        task_type="BINARY_CLASSIFICATION", rebuilt_at=start + timedelta(days=2),
    )
    candidate_snapshot = rebuild_performance_snapshot(
        observations, endpoint_id=candidate_spec.endpoint_id, candidate_id=candidate_spec.candidate_id,
        task_type="BINARY_CLASSIFICATION", rebuilt_at=start + timedelta(days=2),
    )
    test_policy = QualificationPolicy(
        endpoint_name="TEST hERG lifecycle", endpoint_id=candidate_spec.endpoint_id,
        endpoint_contract_version="1.0.0", qualification_kind=QualificationKind.STRATEGY_QUALIFICATION,
        qualification_allowed=True, current_strategy=StrategyType.SINGLE_CORE_MODEL,
        current_policy_version="TEST-POLICY-GATE-V1", current_lifecycle=QualificationLifecycle.ACTIVE,
        output_type="BINARY_CLASSIFICATION",
        minimum_sample_requirement=MinimumSampleRequirement(
            version="TEST-SAMPLE-V1", provenance="TEST_ONLY_SYNTHETIC_DRY_RUN", minimum_n=8,
            minimum_positive_n=4, minimum_negative_n=4, minimum_per_subgroup=4,
        ),
        noninferiority_policy=NonInferiorityPolicy(
            version="TEST-NI-V1", primary_metric="mcc", higher_is_better=True,
            margin=0.05, provenance="TEST_ONLY_SYNTHETIC_DRY_RUN",
        ),
    )
    gate_input = PromotionGateInput(
        endpoint_contract_compatible=True, model_identity_frozen=True, no_leakage=True,
        independent_validation=True, calibration_acceptable=True, subgroup_robust=True,
        no_unacceptable_safety_tradeoff=True, rollback_target_exists=True,
        artifacts_reproducible=True, candidate_metrics=candidate_snapshot.metrics,
        active_metrics=active_snapshot.metrics, n_positive=4, n_negative=4,
        minimum_subgroup_n=4, claimed_improvement=True,
    )
    gate = evaluate_promotion_gate(test_policy, gate_input)
    assert gate.passed

    candidate_metric_values = dict(candidate_snapshot.metrics)
    qualification_record = QualificationRecord(
        qualification_record_id="TEST-QUALIFICATION-HERG-PLATT-V1",
        endpoint_id=candidate_spec.endpoint_id,
        policy_version=candidate_spec.policy_version,
        candidate_strategy=candidate_spec.candidate_strategy,
        candidate_models=tuple(item.model_id for item in candidate_spec.models),
        model_versions=tuple(item.model_version for item in candidate_spec.models),
        checkpoint_hashes=tuple(item.checkpoint_hash for item in candidate_spec.models),
        current_active_strategy=active_spec.candidate_strategy,
        validation_dataset="TEST_ONLY_SYNTHETIC_HERG_DRY_RUN",
        validation_snapshot_hash=candidate_snapshot.statistics_hash,
        validation_type=ValidationType.PROSPECTIVE_INTERNAL,
        prospective_or_retrospective=EvidenceTiming.PROSPECTIVE,
        sample_size=int(candidate_metric_values["n"]),
        primary_metrics=(("mcc", candidate_metric_values["mcc"]),),
        secondary_metrics=tuple(
            (key, candidate_metric_values[key])
            for key in ("balanced_accuracy", "auroc", "auprc", "sensitivity", "specificity")
        ),
        subgroup_metrics=(("TEST-SERIES-A", "PASS_TEST_ONLY"), ("TEST-SERIES-B", "PASS_TEST_ONLY")),
        calibration_metrics=tuple(
            (key, candidate_metric_values[key]) for key in ("brier", "log_loss", "ece")
        ),
        applicability_domain_metrics=(("in_domain_rate", 1.0),),
        known_limitations=("TEST_ONLY_SYNTHETIC", "NOT_SCIENTIFIC_EVIDENCE"),
        qualification_decision=QualificationDecision.QUALIFIED,
        review_timestamp=start + timedelta(days=2, hours=1),
        promotion_status=QualificationLifecycle.VALIDATED,
        rollback_target=active_spec.candidate_id,
        provenance=(("generator", "generate_stage4d5_qualification_artifacts.py"),),
    )

    service = QualificationLifecycleService()
    service.register(active_spec, QualificationLifecycle.ACTIVE)
    service.register(candidate_spec, QualificationLifecycle.SHADOW)
    service.transition(
        candidate_spec.candidate_id, QualificationLifecycle.VALIDATED,
        "TEST gate evidence assembled", qualification_record=qualification_record,
    )
    service.transition(
        candidate_spec.candidate_id, QualificationLifecycle.PRODUCTION_CANDIDATE,
        "TEST promotion gate passed", gate_decision=gate,
    )
    authorization = ManualPromotionAuthorization(
        authorization_id="TEST-MANUAL-AUTHORIZATION-V1", authorized_by="TEST_FIXTURE",
        authorized_at=start + timedelta(days=3), qualification_record_hash=qualification_record.record_hash,
        reason="TEST-only simulated manual activation",
    )
    service.activate(candidate_spec.candidate_id, authorization)
    test_drift_policy = DEFAULT_DRIFT_POLICY.__class__(
        version="TEST-DRIFT-V1", minimum_n=8, performance_absolute_delta=0.10,
        ece_absolute_delta=0.05, domain_out_rate_delta=0.10,
        prevalence_absolute_delta=0.10, chemical_distance_delta=0.10,
        provenance="TEST_ONLY_SYNTHETIC_DRY_RUN",
    )
    drift = detect_drift(
        policy=test_drift_policy,
        reference_metrics=dict(candidate_snapshot.metrics),
        current_metrics={**dict(candidate_snapshot.metrics), "n": 8, "mcc": 0.40, "ece": 0.30},
        primary_metric="mcc", higher_is_better=True,
        reference_model_version="test-v1", current_model_version="test-v1",
        reference_endpoint_contract="1.0.0", current_endpoint_contract="1.0.0",
        reference_out_of_domain_rate=0.0, current_out_of_domain_rate=0.2,
        reference_prevalence=0.5, current_prevalence=0.75, chemical_distance_delta=0.2,
    )
    rollback_event, restored_id = service.rollback(candidate_spec.candidate_id, RollbackReason.PERFORMANCE_REGRESSION)
    assert restored_id == active_spec.candidate_id
    assert service.active_by_endpoint[candidate_spec.endpoint_id] == active_spec.candidate_id

    return {
        "artifact": "stage4d5_dry_run",
        "scope": "TEST_ONLY_SYNTHETIC",
        "real_research_data_used": False,
        "production_state_changed": False,
        "active_and_shadow_stored_separately": True,
        "active_specification": active_spec.to_dict(),
        "candidate_specification": candidate_spec.to_dict(),
        "frozen_predictions": freezes,
        "experimental_results": results,
        "active_metrics": active_snapshot.to_dict(),
        "candidate_metrics": candidate_snapshot.to_dict(),
        "qualification_record": qualification_record.to_dict(),
        "promotion_gate": gate.to_dict(),
        "manual_authorization": {
            "authorization_id": authorization.authorization_id,
            "authorized_by": authorization.authorized_by,
            "authorized_at": authorization.authorized_at.isoformat(),
            "qualification_record_hash": authorization.qualification_record_hash,
            "reason": authorization.reason,
        },
        "lifecycle_events": [
            {
                **{key: (value.value if hasattr(value, "value") else value) for key, value in event.__dict__.items()},
                "occurred_at": event.occurred_at.isoformat(),
            }
            for event in service.events
        ],
        "simulated_drift": drift.to_dict(),
        "rollback": {
            "event_id": rollback_event.event_id,
            "reason": rollback_event.rollback_reason,
            "deterministic_target": restored_id,
            "final_active_candidate_id": service.active_by_endpoint[candidate_spec.endpoint_id],
        },
        "lifecycle_complete": True,
    }


def main() -> None:
    VALIDATION.mkdir(parents=True, exist_ok=True)
    write_json("stage4d5_production_baseline.json", get_production_baseline())
    write_json("stage4d5_qualification_policy.json", build_qualification_policy_artifact())
    write_json("stage4d5_promotion_gates.json", build_promotion_gate_artifact())
    write_json("stage4d5_drift_policy.json", build_drift_artifact())
    write_json("stage4d5_dry_run.json", build_dry_run())
    write_json("stage4d5_strategy_cards.json", {
        "artifact": "stage4d5_strategy_cards",
        "policy_version": QUALIFICATION_POLICY_VERSION,
        "card_count": len(QUALIFICATION_POLICY_REGISTRY),
        "cards": get_strategy_cards(),
    })


if __name__ == "__main__":
    main()
