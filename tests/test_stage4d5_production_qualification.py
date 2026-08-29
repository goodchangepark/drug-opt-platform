"""Stage 4D-5 production qualification and prospective validation gates."""

from __future__ import annotations

import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest
from sqlalchemy import create_engine, select
from sqlalchemy.orm import sessionmaker

from backend.database import Base
from backend.endpoint_strategy_registry import ENDPOINT_STRATEGY_REGISTRY, StrategyType
from backend.main import (
    app,
    health,
    model_strategy_registry,
    qualification_candidates,
    qualification_drift,
    qualification_endpoint,
    qualification_strategies,
)
from backend.production_qualification import (
    ACTIVATION_MODE,
    CANDIDATE_TRACKS,
    DEFAULT_DRIFT_POLICY,
    LEGAL_TRANSITIONS,
    QUALIFICATION_POLICY_REGISTRY,
    STATE_MACHINE_VERSION,
    VALIDATION_EVIDENCE_HIERARCHY,
    AssayEvidenceQuality,
    CandidateSpecification,
    DriftWarning,
    EligibilityStatus,
    EvidenceTiming,
    ExperimentalQualificationResult,
    FrozenModelIdentity,
    ManualPromotionAuthorization,
    MinimumSampleRequirement,
    MonitoringScope,
    NonInferiorityPolicy,
    PerformanceObservation,
    PromotionGateInput,
    ProspectivePredictionFreeze,
    QualificationDecision,
    QualificationEvidenceStore,
    QualificationKind,
    QualificationLifecycle,
    QualificationLifecycleService,
    QualificationPolicy,
    QualificationPredictionFreezeRow,
    QualificationExperimentalResultRow,
    StrategyQualificationRecordRow,
    QualificationLifecycleEventRow,
    QualificationRecord,
    RollbackReason,
    ValidationType,
    canonical_hash,
    classification_metrics,
    detect_drift,
    evaluate_experimental_eligibility,
    evaluate_promotion_gate,
    get_production_baseline,
    get_strategy_cards,
    rebuild_performance_snapshot,
    regression_metrics,
    validate_qualification_registry,
    validate_transition,
)


ROOT = Path(__file__).resolve().parents[1]
UTC = timezone.utc


def identity(name: str = "candidate") -> FrozenModelIdentity:
    return FrozenModelIdentity(name, "model-v1", f"sha256:{name}-checkpoint")


def specification(
    candidate_id: str = "CANDIDATE-V1",
    strategy: StrategyType = StrategyType.SINGLE_CORE_WITH_CALIBRATION,
) -> CandidateSpecification:
    return CandidateSpecification(
        candidate_id=candidate_id,
        endpoint_id="safety_herg_blocker_prob",
        endpoint_contract_version="1.0.0",
        candidate_strategy=strategy,
        models=(identity(candidate_id),),
        policy_version="POLICY-V1",
        candidate_version="CANDIDATE-SPEC-V1",
        decision_threshold=0.5,
        calibration="PLATT_TEST_ONLY",
    )


def freeze(
    frozen_id: str = "FREEZE-1",
    candidate: CandidateSpecification | None = None,
    frozen_at: datetime | None = None,
    probability: float = 0.8,
) -> ProspectivePredictionFreeze:
    candidate = candidate or specification()
    return ProspectivePredictionFreeze(
        frozen_prediction_id=frozen_id,
        compound_version_id="COMPOUND-V1",
        endpoint_id=candidate.endpoint_id,
        endpoint_contract_version=candidate.endpoint_contract_version,
        candidate_id=candidate.candidate_id,
        candidate_specification_hash=candidate.specification_hash,
        strategy=candidate.candidate_strategy,
        models=candidate.models,
        prediction_value=probability,
        probability=probability,
        unit="probability",
        frozen_at=frozen_at or datetime(2026, 8, 1, tzinfo=UTC),
        policy_version=candidate.policy_version,
        standardizer_version=candidate.standardizer_version,
        project_id="PROJECT-1",
        chemical_series_id="SERIES-A",
        applicability_domain="IN_DOMAIN",
        provenance=(("source", "TEST_ONLY"),),
    )


def experimental_result(
    frozen_id: str = "FREEZE-1",
    available_at: datetime | None = None,
    **overrides,
) -> ExperimentalQualificationResult:
    values = {
        "experimental_result_id": f"RESULT-{frozen_id}",
        "frozen_prediction_id": frozen_id,
        "endpoint_id": "safety_herg_blocker_prob",
        "endpoint_contract_version": "1.0.0",
        "experimental_value": 1.0,
        "unit": "binary_label",
        "assay_type": "Automated patch clamp",
        "species": "Human",
        "experiment_date": "2026-08-02",
        "result_available_at": available_at or datetime(2026, 8, 2, tzinfo=UTC),
        "source": "TEST_ONLY",
        "quality": AssayEvidenceQuality.QUALIFICATION_GRADE,
        "protocol_metadata": (("voltage_protocol", "TEST"),),
        "recorded_at": datetime(2026, 8, 2, 1, tzinfo=UTC),
    }
    values.update(overrides)
    return ExperimentalQualificationResult(**values)


def configured_policy() -> QualificationPolicy:
    return QualificationPolicy(
        endpoint_name="TEST hERG",
        endpoint_id="safety_herg_blocker_prob",
        endpoint_contract_version="1.0.0",
        qualification_kind=QualificationKind.STRATEGY_QUALIFICATION,
        qualification_allowed=True,
        current_strategy=StrategyType.SINGLE_CORE_MODEL,
        current_policy_version="TEST-POLICY-V1",
        current_lifecycle=QualificationLifecycle.ACTIVE,
        output_type="BINARY_CLASSIFICATION",
        minimum_sample_requirement=MinimumSampleRequirement(
            "TEST-SAMPLE-V1", "TEST_ONLY", minimum_n=8,
            minimum_positive_n=4, minimum_negative_n=4, minimum_per_subgroup=4,
        ),
        noninferiority_policy=NonInferiorityPolicy(
            "TEST-NI-V1", "mcc", True, 0.05, "TEST_ONLY",
        ),
    )


def passing_gate_input() -> PromotionGateInput:
    return PromotionGateInput(
        endpoint_contract_compatible=True,
        model_identity_frozen=True,
        no_leakage=True,
        independent_validation=True,
        calibration_acceptable=True,
        subgroup_robust=True,
        no_unacceptable_safety_tradeoff=True,
        rollback_target_exists=True,
        artifacts_reproducible=True,
        candidate_metrics=(("n", 8), ("mcc", 0.85), ("ece", 0.04)),
        active_metrics=(("n", 8), ("mcc", 0.80), ("ece", 0.10)),
        n_positive=4,
        n_negative=4,
        minimum_subgroup_n=4,
        claimed_improvement=True,
    )


def qualified_record(candidate: CandidateSpecification) -> QualificationRecord:
    return QualificationRecord(
        qualification_record_id=f"QUAL-{candidate.candidate_id}",
        endpoint_id=candidate.endpoint_id,
        policy_version=candidate.policy_version,
        candidate_strategy=candidate.candidate_strategy,
        candidate_models=tuple(item.model_id for item in candidate.models),
        model_versions=tuple(item.model_version for item in candidate.models),
        checkpoint_hashes=tuple(item.checkpoint_hash for item in candidate.models),
        current_active_strategy=StrategyType.SINGLE_CORE_MODEL,
        validation_dataset="TEST_ONLY",
        validation_snapshot_hash="sha256:test-snapshot",
        validation_type=ValidationType.PROSPECTIVE_INTERNAL,
        prospective_or_retrospective=EvidenceTiming.PROSPECTIVE,
        sample_size=8,
        primary_metrics=(("mcc", 1.0),),
        secondary_metrics=(("auroc", 1.0),),
        subgroup_metrics=(("SERIES-A", "PASS"),),
        calibration_metrics=(("ece", 0.05),),
        applicability_domain_metrics=(("in_domain_rate", 1.0),),
        known_limitations=("TEST_ONLY",),
        qualification_decision=QualificationDecision.QUALIFIED,
        review_timestamp=datetime(2026, 8, 3, tzinfo=UTC),
        promotion_status=QualificationLifecycle.VALIDATED,
        rollback_target="ACTIVE-V1",
        provenance=(("scope", "TEST_ONLY"),),
    )


def observations(candidate_id: str, model_version: str = "model-v1", policy_version: str = "policy-v1"):
    labels = [1, 1, 1, 1, 0, 0, 0, 0]
    probabilities = [0.9, 0.85, 0.8, 0.75, 0.25, 0.2, 0.15, 0.1]
    return [
        PerformanceObservation(
            frozen_prediction_id=f"F-{candidate_id}-{index}",
            endpoint_id="safety_herg_blocker_prob",
            candidate_id=candidate_id,
            model_version_key=model_version,
            policy_version=policy_version,
            project_id="PROJECT-A",
            chemical_series_id="SERIES-A" if index < 4 else "SERIES-B",
            predicted_value=probability,
            experimental_value=float(label),
            probability=probability,
            eligibility=EligibilityStatus.QUALIFICATION_ELIGIBLE,
            applicability_domain="IN_DOMAIN",
        )
        for index, (probability, label) in enumerate(zip(probabilities, labels))
    ]


def test_qualification_registry_covers_all_49_policies_without_violation():
    assert len(QUALIFICATION_POLICY_REGISTRY) == 49
    assert set(QUALIFICATION_POLICY_REGISTRY) == set(ENDPOINT_STRATEGY_REGISTRY)
    assert validate_qualification_registry() == []


def test_lifecycle_enum_and_versioned_legal_transitions_are_complete():
    assert {item.value for item in QualificationLifecycle} == {
        "RESEARCH_ONLY", "SHADOW", "VALIDATED", "PRODUCTION_CANDIDATE",
        "ACTIVE", "RETIRED", "ROLLED_BACK",
    }
    assert STATE_MACHINE_VERSION == "stage4d5-lifecycle-v1"
    assert QualificationLifecycle.ACTIVE not in LEGAL_TRANSITIONS[QualificationLifecycle.SHADOW]
    assert QualificationLifecycle.PRODUCTION_CANDIDATE in LEGAL_TRANSITIONS[QualificationLifecycle.VALIDATED]


def test_shadow_cannot_jump_directly_to_active_even_with_authorization():
    authorization = ManualPromotionAuthorization(
        "AUTH", "reviewer", datetime.now(UTC), "hash", "manual test",
    )
    with pytest.raises(ValueError, match="illegal qualification transition"):
        validate_transition(QualificationLifecycle.SHADOW, QualificationLifecycle.ACTIVE,
                            manual_authorization=authorization)


def test_activation_requires_explicit_manual_authorization():
    with pytest.raises(PermissionError, match="manual promotion authorization"):
        validate_transition(QualificationLifecycle.PRODUCTION_CANDIDATE, QualificationLifecycle.ACTIVE)
    validate_transition(
        QualificationLifecycle.PRODUCTION_CANDIDATE,
        QualificationLifecycle.ACTIVE,
        manual_authorization=ManualPromotionAuthorization(
            "AUTH", "reviewer", datetime.now(UTC), "record-hash", "authorized test",
        ),
    )


def test_candidate_identity_freezes_models_threshold_calibration_and_policy():
    service = QualificationLifecycleService()
    original = specification("FROZEN-ID")
    service.register(original, QualificationLifecycle.SHADOW)
    changed = CandidateSpecification(
        candidate_id="FROZEN-ID", endpoint_id=original.endpoint_id,
        endpoint_contract_version=original.endpoint_contract_version,
        candidate_strategy=original.candidate_strategy,
        models=(FrozenModelIdentity("changed", "model-v2", "sha256:changed"),),
        policy_version=original.policy_version, candidate_version="CHANGED-V2",
        decision_threshold=0.6, calibration="CHANGED",
    )
    with pytest.raises(ValueError, match="new candidate_id"):
        service.register(changed, QualificationLifecycle.SHADOW)


def test_full_manual_lifecycle_and_deterministic_rollback():
    service = QualificationLifecycleService()
    active = specification("ACTIVE-V1", StrategyType.SINGLE_CORE_MODEL)
    candidate = specification("SHADOW-V1")
    service.register(active, QualificationLifecycle.ACTIVE)
    service.register(candidate, QualificationLifecycle.SHADOW)
    service.transition(
        candidate.candidate_id, QualificationLifecycle.VALIDATED, "validated",
        qualification_record=qualified_record(candidate),
    )
    gate = evaluate_promotion_gate(configured_policy(), passing_gate_input())
    service.transition(
        candidate.candidate_id, QualificationLifecycle.PRODUCTION_CANDIDATE, "gate passed",
        gate_decision=gate,
    )
    auth = ManualPromotionAuthorization("AUTH-1", "reviewer", datetime.now(UTC), "record", "manual activation")
    service.activate(candidate.candidate_id, auth)
    assert service.active_by_endpoint[candidate.endpoint_id] == candidate.candidate_id
    assert service.strategies[active.candidate_id].state == QualificationLifecycle.RETIRED
    event, target = service.rollback(candidate.candidate_id, RollbackReason.PERFORMANCE_REGRESSION)
    assert event.rollback_reason == "PERFORMANCE_REGRESSION"
    assert target == active.candidate_id
    assert service.active_by_endpoint[candidate.endpoint_id] == active.candidate_id
    assert service.strategies[candidate.candidate_id].state == QualificationLifecycle.ROLLED_BACK


def test_validated_and_production_candidate_transitions_require_evidence_objects():
    service = QualificationLifecycleService()
    candidate = specification("EVIDENCE-GATED")
    service.register(candidate, QualificationLifecycle.SHADOW)
    with pytest.raises(ValueError, match="structured QUALIFIED"):
        service.transition(candidate.candidate_id, QualificationLifecycle.VALIDATED, "unsupported")
    service.transition(
        candidate.candidate_id, QualificationLifecycle.VALIDATED, "qualified",
        qualification_record=qualified_record(candidate),
    )
    with pytest.raises(ValueError, match="passing conjunctive"):
        service.transition(candidate.candidate_id, QualificationLifecycle.PRODUCTION_CANDIDATE, "unsupported")


def test_prediction_freeze_requires_value_and_timezone_and_has_stable_hash():
    row = freeze()
    assert len(row.record_hash) == 64
    assert row.record_hash == freeze().record_hash
    with pytest.raises(ValueError, match="value or probability"):
        ProspectivePredictionFreeze(**{**row.__dict__, "prediction_value": None, "probability": None})
    with pytest.raises(ValueError, match="timezone-aware"):
        ProspectivePredictionFreeze(**{**row.__dict__, "frozen_at": datetime(2026, 8, 1)})


def test_prospective_experimental_link_is_eligible_only_after_freeze():
    frozen = freeze()
    result = experimental_result()
    eligible = evaluate_experimental_eligibility(frozen, result)
    assert eligible.status == EligibilityStatus.QUALIFICATION_ELIGIBLE
    assert eligible.counts_toward_qualification is True
    post_hoc = evaluate_experimental_eligibility(
        freeze(frozen_at=datetime(2026, 8, 3, tzinfo=UTC)), result,
    )
    assert post_hoc.status == EligibilityStatus.POST_HOC_PREDICTION
    assert post_hoc.counts_toward_qualification is False


def test_missing_freeze_metadata_and_incompatible_assay_fail_closed():
    result = experimental_result()
    assert evaluate_experimental_eligibility(None, result).status == EligibilityStatus.NO_FROZEN_PREDICTION
    missing = experimental_result(source="")
    assert evaluate_experimental_eligibility(freeze(), missing).status == EligibilityStatus.MISSING_METADATA
    wrong_endpoint = experimental_result(endpoint_id="safety_ames_mutagenicity_prob")
    assert evaluate_experimental_eligibility(freeze(), wrong_endpoint).status == EligibilityStatus.INCOMPATIBLE
    limited = experimental_result(quality=AssayEvidenceQuality.LIMITED)
    assert evaluate_experimental_eligibility(freeze(), limited).status == EligibilityStatus.LIMITED
    unrelated = experimental_result(assay_type="Ames bacterial reverse mutation")
    assert evaluate_experimental_eligibility(freeze(), unrelated).status == EligibilityStatus.INCOMPATIBLE


def test_contract_specific_assay_species_and_unit_compatibility():
    sol_spec = CandidateSpecification(
        "SOL-CAND", "solubility_aqueous_logs", "1.0.0", StrategyType.SINGLE_CORE_MODEL,
        (identity("sol"),), "SOL-POLICY", "SOL-V1",
    )
    sol_freeze = ProspectivePredictionFreeze(
        "SOL-FREEZE", "CMP", sol_spec.endpoint_id, "1.0.0", sol_spec.candidate_id,
        sol_spec.specification_hash, sol_spec.candidate_strategy, sol_spec.models, -3.0, None,
        "log10(mol/L)", datetime(2026, 8, 1, tzinfo=UTC), sol_spec.policy_version,
        sol_spec.standardizer_version,
    )
    incompatible = ExperimentalQualificationResult(
        "SOL-RESULT", "SOL-FREEZE", sol_spec.endpoint_id, "1.0.0", -2.8,
        "mg/mL", "Organic solvent solubility", "Chemical / In Vitro", "2026-08-02",
        datetime(2026, 8, 2, tzinfo=UTC), "TEST", AssayEvidenceQuality.QUALIFICATION_GRADE,
    )
    decision = evaluate_experimental_eligibility(sol_freeze, incompatible)
    assert decision.status == EligibilityStatus.INCOMPATIBLE
    assert {"ASSAY_TYPE_INCOMPATIBLE", "UNIT_INCOMPATIBLE"} <= set(decision.reason_codes)


def test_validation_taxonomy_is_explicit_and_non_equivalent():
    assert {item.value for item in ValidationType} == {
        "TRAINING_INTERNAL", "CROSS_VALIDATION", "EXTERNAL_RETROSPECTIVE",
        "PSEUDO_PROSPECTIVE", "PROSPECTIVE_INTERNAL", "PROSPECTIVE_EXTERNAL",
        "CLINICAL_RETROSPECTIVE", "CLINICAL_PROSPECTIVE",
    }
    assert len(set(VALIDATION_EVIDENCE_HIERARCHY.values())) == len(ValidationType)
    assert VALIDATION_EVIDENCE_HIERARCHY[ValidationType.PROSPECTIVE_EXTERNAL] > \
        VALIDATION_EVIDENCE_HIERARCHY[ValidationType.CROSS_VALIDATION]


def test_classification_metrics_include_discrimination_calibration_and_class_balance():
    metrics = classification_metrics([0.9, 0.8, 0.2, 0.1], [1, 1, 0, 0])
    assert metrics["n"] == 4
    assert metrics["mcc"] == pytest.approx(1.0)
    assert metrics["balanced_accuracy"] == pytest.approx(1.0)
    assert metrics["auroc"] == pytest.approx(1.0)
    assert metrics["auprc"] == pytest.approx(1.0)
    assert 0 < metrics["brier"] < 0.1
    assert 0 < metrics["log_loss"] < 0.3
    assert metrics["sensitivity"] == pytest.approx(1.0)
    assert metrics["specificity"] == pytest.approx(1.0)
    assert metrics["ece"] is not None


def test_regression_metrics_include_bias_rank_and_valid_fold_metrics():
    metrics = regression_metrics([-3.0, -4.0, -5.0], [-3.1, -3.9, -5.0])
    assert metrics["n"] == 3
    assert metrics["mae"] == pytest.approx(0.0666666667)
    assert metrics["rmse"] is not None and metrics["bias"] == pytest.approx(0.0)
    assert metrics["spearman"] == pytest.approx(1.0)
    assert metrics["within_2_fold"] == 1.0
    nonlog = regression_metrics([50.0], [55.0], log10_fold_metrics=False)
    assert nonlog["within_2_fold"] is None and nonlog["within_3_fold"] is None


def test_rebuild_is_reproducible_and_ignores_ineligible_observations():
    rows = observations("CAND")
    rows.append(PerformanceObservation(
        "POST-HOC", "safety_herg_blocker_prob", "CAND", "model-v1", "policy-v1",
        "PROJECT-A", "SERIES-A", 0.99, 0.0, 0.99, EligibilityStatus.POST_HOC_PREDICTION,
    ))
    first = rebuild_performance_snapshot(
        rows, endpoint_id="safety_herg_blocker_prob", candidate_id="CAND",
        task_type="BINARY_CLASSIFICATION", rebuilt_at=datetime(2026, 8, 3, tzinfo=UTC),
    )
    second = rebuild_performance_snapshot(
        list(reversed(rows)), endpoint_id="safety_herg_blocker_prob", candidate_id="CAND",
        task_type="BINARY_CLASSIFICATION", rebuilt_at=datetime(2026, 8, 4, tzinfo=UTC),
    )
    assert dict(first.metrics)["n"] == 8
    assert first.metrics == second.metrics
    assert first.source_snapshot_hash == second.source_snapshot_hash
    assert first.statistics_hash == second.statistics_hash


def test_scope_endpoint_model_and_policy_versions_remain_isolated():
    rows = observations("CAND")
    rows.extend(observations("OTHER"))
    project = rebuild_performance_snapshot(
        rows, endpoint_id="safety_herg_blocker_prob", candidate_id="CAND",
        task_type="BINARY_CLASSIFICATION", scope=MonitoringScope.PROJECT, scope_key="PROJECT-A",
    )
    assert dict(project.metrics)["n"] == 8
    mixed_model = rows[:2] + [PerformanceObservation(
        **{**rows[2].__dict__, "model_version_key": "model-v2"}
    )]
    with pytest.raises(ValueError, match="model versions"):
        rebuild_performance_snapshot(
            mixed_model, endpoint_id="safety_herg_blocker_prob", candidate_id="CAND",
            task_type="BINARY_CLASSIFICATION",
        )
    mixed_policy = rows[:2] + [PerformanceObservation(
        **{**rows[2].__dict__, "policy_version": "policy-v2"}
    )]
    with pytest.raises(ValueError, match="policy versions"):
        rebuild_performance_snapshot(
            mixed_policy, endpoint_id="safety_herg_blocker_prob", candidate_id="CAND",
            task_type="BINARY_CLASSIFICATION",
        )


def test_unconfigured_real_endpoint_requirements_fail_closed():
    policy = QUALIFICATION_POLICY_REGISTRY["hERG liability"]
    decision = evaluate_promotion_gate(policy, passing_gate_input())
    assert decision.passed is False
    assert decision.decision == QualificationDecision.INSUFFICIENT_EVIDENCE
    assert "MINIMUM_SAMPLE_REQUIREMENTS_CONFIGURED" in decision.reason_codes
    assert "NONINFERIORITY_MARGIN_CONFIGURED" in decision.reason_codes


def test_promotion_gate_requires_every_gate_not_one_metric():
    passed = evaluate_promotion_gate(configured_policy(), passing_gate_input())
    assert passed.passed is True
    failed_input = PromotionGateInput(**{
        **passing_gate_input().__dict__, "no_leakage": False, "subgroup_robust": False,
    })
    failed = evaluate_promotion_gate(configured_policy(), failed_input)
    assert failed.passed is False
    assert {"NO_LEAKAGE", "SUBGROUP_ROBUST"} <= set(failed.reason_codes)


def test_noninferiority_supports_documented_equivalent_benefit():
    base = passing_gate_input()
    equivalent = PromotionGateInput(**{
        **base.__dict__,
        "candidate_metrics": (("n", 8), ("mcc", 0.78), ("ece", 0.03)),
        "active_metrics": (("n", 8), ("mcc", 0.80), ("ece", 0.10)),
        "claimed_improvement": False,
        "equivalent_benefits": ("BETTER_CALIBRATION",),
    })
    decision = evaluate_promotion_gate(configured_policy(), equivalent)
    assert decision.passed is True


def test_drift_detection_is_deterministic_and_never_acts_automatically():
    assessment = detect_drift(
        policy=DEFAULT_DRIFT_POLICY,
        reference_metrics={"n": 100, "mcc": 0.8, "ece": 0.05},
        current_metrics={"n": 40, "mcc": 0.5, "ece": 0.2},
        primary_metric="mcc", higher_is_better=True,
        reference_model_version="v1", current_model_version="v2",
        reference_endpoint_contract="1.0.0", current_endpoint_contract="2.0.0",
        reference_out_of_domain_rate=0.05, current_out_of_domain_rate=0.3,
        reference_prevalence=0.4, current_prevalence=0.7,
        chemical_distance_delta=0.2,
    )
    assert assessment.status == "REVIEW_REQUIRED"
    assert assessment.review_required is True
    assert assessment.automatic_action == "NONE"
    assert {
        DriftWarning.PERFORMANCE_DRIFT, DriftWarning.CALIBRATION_DRIFT,
        DriftWarning.DOMAIN_SHIFT, DriftWarning.PRIOR_SHIFT,
        DriftWarning.MODEL_VERSION_CHANGED, DriftWarning.ENDPOINT_MISMATCH,
    } <= set(assessment.warnings)


def test_drift_insufficient_data_is_review_not_retraining():
    assessment = detect_drift(
        policy=DEFAULT_DRIFT_POLICY,
        reference_metrics={"n": 100, "mae": 0.4}, current_metrics={"n": 2, "mae": 0.1},
        primary_metric="mae", higher_is_better=False,
        reference_model_version="v1", current_model_version="v1",
        reference_endpoint_contract="1.0.0", current_endpoint_contract="1.0.0",
    )
    assert assessment.warnings == (DriftWarning.INSUFFICIENT_DATA,)
    assert assessment.automatic_action == "NONE"


def test_model_unavailable_endpoints_cannot_enter_promotion_pipeline():
    unavailable = []
    for name, governed in ENDPOINT_STRATEGY_REGISTRY.items():
        if governed.primary_strategy == StrategyType.MODEL_UNAVAILABLE:
            unavailable.append(name)
            policy = QUALIFICATION_POLICY_REGISTRY[name]
            assert policy.qualification_kind == QualificationKind.EXCLUDED
            assert policy.qualification_allowed is False
            assert name not in CANDIDATE_TRACKS
    assert {
        "Dog liver microsomal intrinsic clearance", "Monkey liver microsomal intrinsic clearance",
        "BCRP inhibitor", "BSEP inhibitor", "OATP1B1 inhibitor", "OCT2 inhibitor",
        "P-gp substrate", "pKa (quantitative ML)", "logD7.4 (quantitative ML)",
    } <= set(unavailable)
    with pytest.raises(ValueError, match="MODEL_UNAVAILABLE"):
        specification(strategy=StrategyType.MODEL_UNAVAILABLE)


def test_mechanistic_endpoints_use_method_qualification_not_ml_consensus():
    for name in ("PK Systemic Clearance", "PK Volume of Distribution", "PK Bioavailability", "PK Simulation"):
        policy = QUALIFICATION_POLICY_REGISTRY[name]
        assert policy.qualification_kind == QualificationKind.METHOD_QUALIFICATION
        assert policy.current_strategy == StrategyType.MECHANISTIC_NO_CONSENSUS


def test_priority_tracks_remain_shadow_and_scientifically_bounded():
    assert set(CANDIDATE_TRACKS) == {"Solubility", "Permeability", "CYP3A4 inhibitor", "hERG liability"}
    assert all(track["state"] == "SHADOW" and not track["automatic_activation"]
               for track in CANDIDATE_TRACKS.values())
    assert CANDIDATE_TRACKS["Solubility"]["excluded_models"] == ["rdkit_gbr_solubility_v1"]
    assert CANDIDATE_TRACKS["Permeability"]["qualification_status"] == "INSUFFICIENT_EVIDENCE"
    assert CANDIDATE_TRACKS["CYP3A4 inhibitor"]["weights"] == [0.9578, 0.0422]
    assert CANDIDATE_TRACKS["CYP3A4 inhibitor"]["dynamic_adaptive_promotion"] == "CLOSED_NO_ADAPTIVE_VALUE"
    assert CANDIDATE_TRACKS["hERG liability"]["candidate_models"] == ["admetica_safety_herg"]
    assert CANDIDATE_TRACKS["hERG liability"]["supporting_only_models"] == ["physchem_herg_v1"]


def test_sql_evidence_store_is_append_only_and_links_eligibility():
    engine = create_engine("sqlite+pysqlite:///:memory:")
    Base.metadata.create_all(engine, tables=[
        QualificationPredictionFreezeRow.__table__, QualificationExperimentalResultRow.__table__,
        StrategyQualificationRecordRow.__table__, QualificationLifecycleEventRow.__table__,
    ])
    SessionLocal = sessionmaker(bind=engine)
    session = SessionLocal()
    store = QualificationEvidenceStore(session)
    frozen = freeze()
    freeze_row = store.freeze_prediction(frozen)
    result_row, eligibility = store.link_experimental_result(experimental_result())
    session.commit()
    assert freeze_row.record_hash == frozen.record_hash
    assert result_row.eligibility_status == EligibilityStatus.QUALIFICATION_ELIGIBLE.value
    assert eligibility.counts_toward_qualification is True
    freeze_row.unit = "changed"
    with pytest.raises(ValueError, match="append-only"):
        session.flush()
    session.rollback()
    session.close()


def test_structured_qualification_record_persists_all_required_fields_and_is_immutable():
    engine = create_engine("sqlite+pysqlite:///:memory:")
    Base.metadata.create_all(engine, tables=[StrategyQualificationRecordRow.__table__])
    session = sessionmaker(bind=engine)()
    record = QualificationRecord(
        qualification_record_id="QUAL-1", endpoint_id="safety_herg_blocker_prob",
        policy_version="POLICY-V1", candidate_strategy=StrategyType.SINGLE_CORE_WITH_CALIBRATION,
        candidate_models=("m1",), model_versions=("v1",), checkpoint_hashes=("sha256:x",),
        current_active_strategy=StrategyType.SINGLE_CORE_MODEL,
        validation_dataset="TEST frozen dataset", validation_snapshot_hash="sha256:snapshot",
        validation_type=ValidationType.PROSPECTIVE_INTERNAL,
        prospective_or_retrospective=EvidenceTiming.PROSPECTIVE,
        sample_size=8, primary_metrics=(("mcc", 1.0),), secondary_metrics=(("auroc", 1.0),),
        subgroup_metrics=(("SERIES-A", "ROBUST_TEST_ONLY"),),
        calibration_metrics=(("ece", 0.1),),
        applicability_domain_metrics=(("in_domain_rate", 1.0),),
        known_limitations=("TEST_ONLY",), qualification_decision=QualificationDecision.QUALIFIED,
        review_timestamp=datetime(2026, 8, 3, tzinfo=UTC),
        promotion_status=QualificationLifecycle.VALIDATED,
        rollback_target="ACTIVE-V1", provenance=(("scope", "TEST_ONLY"),),
    )
    row = QualificationEvidenceStore(session).append_qualification_record(record)
    session.commit()
    assert row.validation_snapshot_hash == "sha256:snapshot"
    assert row.checkpoint_hashes_json == ["sha256:x"]
    assert row.record_hash == record.record_hash
    row.sample_size = 9
    with pytest.raises(ValueError, match="append-only"):
        session.flush()
    session.rollback()
    session.close()


def test_strategy_cards_cover_every_endpoint_without_unsupported_confidence_claims():
    cards = get_strategy_cards()
    assert len(cards) == 49
    assert {row["endpoint"] for row in cards} == set(ENDPOINT_STRATEGY_REGISTRY)
    assert all(row["unsupported_confidence_claims"] is False for row in cards)
    unavailable = [row for row in cards if not row["qualification_allowed"]]
    assert unavailable and all(row["qualification_kind"] == "EXCLUDED" for row in unavailable)


def test_production_baseline_preserves_all_current_active_policies():
    baseline = get_production_baseline()
    expected = [row for row in ENDPOINT_STRATEGY_REGISTRY.values() if row.promotion_status.value == "ACTIVE"]
    assert baseline["active_policy_count"] == len(expected) == 27
    assert baseline["production_behavior_changed"] is False
    for row in baseline["active_endpoint_policies"]:
        current = ENDPOINT_STRATEGY_REGISTRY[row["endpoint"]]
        assert row["primary_strategy"] == current.primary_strategy.value
        assert row["decision_threshold"] == current.decision_threshold
        assert row["calibration_status"] == current.calibration_status.value


def test_read_only_qualification_api_and_backward_compatibility():
    routes = {route.path: route for route in app.routes if hasattr(route, "path")}
    for path in (
        "/api/qualification/strategies", "/api/qualification/endpoint/{endpoint_id}",
        "/api/qualification/candidates", "/api/qualification/drift",
    ):
        assert path in routes
        assert routes[path].methods == {"GET"}
    assert not any(
        path.startswith("/api/qualification") and (route.methods or set()) & {"POST", "PUT", "PATCH", "DELETE"}
        for path, route in routes.items()
    )
    assert qualification_strategies()["read_only"] is True
    assert qualification_strategies()["strategy_count"] == 49
    assert qualification_endpoint("safety_herg_blocker_prob")["candidate_track"]["state"] == "SHADOW"
    assert qualification_candidates()["automatic_activation"] is False
    assert qualification_drift()["automatic_action"] == "NONE"
    assert health()["status"] == "ok"
    assert model_strategy_registry()["production_behavior_changed"] is False


def test_required_artifacts_match_runtime_and_dry_run_completes():
    required = {
        "stage4d5_production_baseline.json", "stage4d5_qualification_policy.json",
        "stage4d5_promotion_gates.json", "stage4d5_drift_policy.json",
        "stage4d5_dry_run.json", "stage4d5_strategy_cards.json",
    }
    for name in required:
        assert (ROOT / "validation" / name).exists()
    policy_artifact = json.loads((ROOT / "validation/stage4d5_qualification_policy.json").read_text())
    assert policy_artifact["qualification_policy_count"] == 49
    assert policy_artifact["registry_violations"] == []
    assert policy_artifact["automatic_shadow_activation"] is False
    dry_run = json.loads((ROOT / "validation/stage4d5_dry_run.json").read_text())
    assert dry_run["scope"] == "TEST_ONLY_SYNTHETIC"
    assert dry_run["real_research_data_used"] is False
    assert dry_run["production_state_changed"] is False
    assert dry_run["promotion_gate"]["passed"] is True
    assert dry_run["simulated_drift"]["status"] == "REVIEW_REQUIRED"
    assert dry_run["rollback"]["final_active_candidate_id"] == "TEST-ACTIVE-HERG-V1"
    assert dry_run["lifecycle_complete"] is True


def test_stage4d5_has_no_frontend_or_ai_changes():
    stage_files = {
        "backend/production_qualification.py", "backend/main.py",
        "scripts/generate_stage4d5_qualification_artifacts.py",
        "validation/stage4d5_production_baseline.json",
        "validation/stage4d5_qualification_policy.json",
        "validation/stage4d5_promotion_gates.json",
        "validation/stage4d5_drift_policy.json",
        "validation/stage4d5_dry_run.json",
        "validation/stage4d5_strategy_cards.json",
        "docs/stage4d5-production-qualification.md",
        "docs/stage4d5-prospective-validation.md",
        "docs/stage4d5-drift-monitoring.md",
        "docs/stage4d5-activation-and-rollback.md",
        "tests/test_stage4d5_production_qualification.py",
    }
    assert all(not path.startswith("frontend/") for path in stage_files)
    assert "LLM" not in " ".join(path.lower() for path in stage_files)
