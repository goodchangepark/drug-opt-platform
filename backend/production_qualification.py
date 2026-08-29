"""Stage 4D-5 production qualification and prospective validation governance.

This module is deliberately separate from prediction execution.  It records
immutable evidence, evaluates deterministic scientific gates, and exposes an
internal lifecycle service.  It never selects a production model, retrains a
model, or activates a shadow strategy automatically.
"""

from __future__ import annotations

import enum
import hashlib
import json
import math
import uuid
from dataclasses import asdict, dataclass, field
from datetime import date, datetime, timezone
from typing import Any, Iterable, Mapping, Sequence

from sqlalchemy import JSON, Boolean, DateTime, Float, Integer, String, Text, event, inspect
from sqlalchemy.orm import Mapped, Session, mapped_column

from .database import Base
from .admet_predictor import MODEL_SPECS
from .endpoint_contracts import ENDPOINT_CONTRACTS, OutputType
from .endpoint_strategy_registry import (
    ENDPOINT_STRATEGY_REGISTRY,
    CalibrationStatus,
    PromotionStatus,
    StrategyType,
)
from .models import utcnow


QUALIFICATION_POLICY_VERSION = "stage4d5-qualification-policy-v1"
STATE_MACHINE_VERSION = "stage4d5-lifecycle-v1"
DRIFT_POLICY_VERSION = "stage4d5-drift-review-v1"
ACTIVATION_MODE = "MANUAL_PROMOTION_REQUIRED"


class QualificationLifecycle(str, enum.Enum):
    RESEARCH_ONLY = "RESEARCH_ONLY"
    SHADOW = "SHADOW"
    VALIDATED = "VALIDATED"
    PRODUCTION_CANDIDATE = "PRODUCTION_CANDIDATE"
    ACTIVE = "ACTIVE"
    RETIRED = "RETIRED"
    ROLLED_BACK = "ROLLED_BACK"


class ValidationType(str, enum.Enum):
    TRAINING_INTERNAL = "TRAINING_INTERNAL"
    CROSS_VALIDATION = "CROSS_VALIDATION"
    EXTERNAL_RETROSPECTIVE = "EXTERNAL_RETROSPECTIVE"
    PSEUDO_PROSPECTIVE = "PSEUDO_PROSPECTIVE"
    PROSPECTIVE_INTERNAL = "PROSPECTIVE_INTERNAL"
    PROSPECTIVE_EXTERNAL = "PROSPECTIVE_EXTERNAL"
    CLINICAL_RETROSPECTIVE = "CLINICAL_RETROSPECTIVE"
    CLINICAL_PROSPECTIVE = "CLINICAL_PROSPECTIVE"


class EvidenceTiming(str, enum.Enum):
    RETROSPECTIVE = "RETROSPECTIVE"
    PROSPECTIVE = "PROSPECTIVE"


class QualificationKind(str, enum.Enum):
    STRATEGY_QUALIFICATION = "STRATEGY_QUALIFICATION"
    METHOD_QUALIFICATION = "METHOD_QUALIFICATION"
    EXCLUDED = "EXCLUDED"


class EligibilityStatus(str, enum.Enum):
    QUALIFICATION_ELIGIBLE = "QUALIFICATION_ELIGIBLE"
    LIMITED = "LIMITED"
    INCOMPATIBLE = "INCOMPATIBLE"
    MISSING_METADATA = "MISSING_METADATA"
    NO_FROZEN_PREDICTION = "NO_FROZEN_PREDICTION"
    POST_HOC_PREDICTION = "POST_HOC_PREDICTION"


class QualificationDecision(str, enum.Enum):
    QUALIFIED = "QUALIFIED"
    NOT_QUALIFIED = "NOT_QUALIFIED"
    INSUFFICIENT_EVIDENCE = "INSUFFICIENT_EVIDENCE"
    REVIEW_REQUIRED = "REVIEW_REQUIRED"


class DriftWarning(str, enum.Enum):
    PERFORMANCE_DRIFT = "PERFORMANCE_DRIFT"
    CALIBRATION_DRIFT = "CALIBRATION_DRIFT"
    DOMAIN_SHIFT = "DOMAIN_SHIFT"
    PRIOR_SHIFT = "PRIOR_SHIFT"
    INSUFFICIENT_DATA = "INSUFFICIENT_DATA"
    MODEL_VERSION_CHANGED = "MODEL_VERSION_CHANGED"
    ENDPOINT_MISMATCH = "ENDPOINT_MISMATCH"


class RollbackReason(str, enum.Enum):
    PERFORMANCE_REGRESSION = "PERFORMANCE_REGRESSION"
    CALIBRATION_FAILURE = "CALIBRATION_FAILURE"
    RUNTIME_FAILURE = "RUNTIME_FAILURE"
    DATA_QUALITY_ISSUE = "DATA_QUALITY_ISSUE"
    POLICY_ERROR = "POLICY_ERROR"
    MANUAL_ROLLBACK = "MANUAL_ROLLBACK"


class MonitoringScope(str, enum.Enum):
    GLOBAL = "GLOBAL"
    PROJECT = "PROJECT"
    CHEMICAL_SERIES = "CHEMICAL_SERIES"
    MODEL_VERSION = "MODEL_VERSION"
    POLICY_VERSION = "POLICY_VERSION"


class AssayEvidenceQuality(str, enum.Enum):
    QUALIFICATION_GRADE = "QUALIFICATION_GRADE"
    LIMITED = "LIMITED"
    INCOMPATIBLE = "INCOMPATIBLE"


VALIDATION_EVIDENCE_HIERARCHY: Mapping[ValidationType, int] = {
    ValidationType.TRAINING_INTERNAL: 10,
    ValidationType.CROSS_VALIDATION: 20,
    ValidationType.EXTERNAL_RETROSPECTIVE: 30,
    ValidationType.PSEUDO_PROSPECTIVE: 40,
    ValidationType.PROSPECTIVE_INTERNAL: 50,
    ValidationType.PROSPECTIVE_EXTERNAL: 60,
    ValidationType.CLINICAL_RETROSPECTIVE: 70,
    ValidationType.CLINICAL_PROSPECTIVE: 80,
}


LEGAL_TRANSITIONS: Mapping[QualificationLifecycle, frozenset[QualificationLifecycle]] = {
    QualificationLifecycle.RESEARCH_ONLY: frozenset({
        QualificationLifecycle.SHADOW, QualificationLifecycle.RETIRED,
    }),
    QualificationLifecycle.SHADOW: frozenset({
        QualificationLifecycle.VALIDATED, QualificationLifecycle.RETIRED,
    }),
    QualificationLifecycle.VALIDATED: frozenset({
        QualificationLifecycle.SHADOW,
        QualificationLifecycle.PRODUCTION_CANDIDATE,
        QualificationLifecycle.RETIRED,
    }),
    QualificationLifecycle.PRODUCTION_CANDIDATE: frozenset({
        QualificationLifecycle.VALIDATED,
        QualificationLifecycle.ACTIVE,
        QualificationLifecycle.RETIRED,
    }),
    QualificationLifecycle.ACTIVE: frozenset({
        QualificationLifecycle.RETIRED,
        QualificationLifecycle.ROLLED_BACK,
    }),
    QualificationLifecycle.RETIRED: frozenset(),
    QualificationLifecycle.ROLLED_BACK: frozenset({QualificationLifecycle.RETIRED}),
}


def _jsonable(value: Any) -> Any:
    if isinstance(value, enum.Enum):
        return value.value
    if isinstance(value, (datetime, date)):
        return value.isoformat()
    if hasattr(value, "__dataclass_fields__"):
        return {key: _jsonable(item) for key, item in asdict(value).items()}
    if hasattr(value, "to_dict"):
        return value.to_dict()
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in sorted(value.items(), key=lambda row: str(row[0]))}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_jsonable(item) for item in value]
    return value


def canonical_hash(value: Any) -> str:
    payload = json.dumps(_jsonable(value), sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class FrozenModelIdentity:
    model_id: str
    model_version: str
    checkpoint_hash: str

    def __post_init__(self) -> None:
        if not self.model_id or not self.model_version or not self.checkpoint_hash:
            raise ValueError("model_id, model_version, and checkpoint_hash are required")


@dataclass(frozen=True)
class CandidateSpecification:
    candidate_id: str
    endpoint_id: str
    endpoint_contract_version: str
    candidate_strategy: StrategyType
    models: tuple[FrozenModelIdentity, ...]
    policy_version: str
    candidate_version: str
    weights: tuple[float, ...] = ()
    decision_threshold: float | None = None
    calibration: str = "NOT_APPLICABLE"
    standardizer_version: str = "CHEM_STANDARDIZER_V1"
    qualification_kind: QualificationKind = QualificationKind.STRATEGY_QUALIFICATION

    def __post_init__(self) -> None:
        if not self.candidate_id or not self.endpoint_id or not self.policy_version or not self.candidate_version:
            raise ValueError("candidate, endpoint, policy, and candidate versions are required")
        if self.candidate_strategy == StrategyType.MODEL_UNAVAILABLE:
            raise ValueError("MODEL_UNAVAILABLE cannot enter the qualification pipeline")
        if self.weights and len(self.weights) != len(self.models):
            raise ValueError("weights must match frozen model identities")

    @property
    def specification_hash(self) -> str:
        return canonical_hash(self)

    def to_dict(self) -> dict[str, Any]:
        row = _jsonable(self)
        row["specification_hash"] = self.specification_hash
        return row


@dataclass(frozen=True)
class ProspectivePredictionFreeze:
    frozen_prediction_id: str
    compound_version_id: str
    endpoint_id: str
    endpoint_contract_version: str
    candidate_id: str
    candidate_specification_hash: str
    strategy: StrategyType
    models: tuple[FrozenModelIdentity, ...]
    prediction_value: float | None
    probability: float | None
    unit: str
    frozen_at: datetime
    policy_version: str
    standardizer_version: str
    project_id: str = ""
    chemical_series_id: str = ""
    applicability_domain: str = "UNKNOWN"
    provenance: tuple[tuple[str, str], ...] = ()

    def __post_init__(self) -> None:
        if self.prediction_value is None and self.probability is None:
            raise ValueError("a frozen prediction value or probability is required")
        if self.frozen_at.tzinfo is None:
            raise ValueError("frozen_at must be timezone-aware")

    @property
    def record_hash(self) -> str:
        return canonical_hash(self)

    def to_dict(self) -> dict[str, Any]:
        row = _jsonable(self)
        row["record_hash"] = self.record_hash
        return row


@dataclass(frozen=True)
class ExperimentalQualificationResult:
    experimental_result_id: str
    frozen_prediction_id: str
    endpoint_id: str
    endpoint_contract_version: str
    experimental_value: float
    unit: str
    assay_type: str
    species: str
    experiment_date: str
    result_available_at: datetime
    source: str
    quality: AssayEvidenceQuality
    protocol_metadata: tuple[tuple[str, str], ...] = ()
    recorded_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))

    @property
    def record_hash(self) -> str:
        return canonical_hash(self)

    def to_dict(self) -> dict[str, Any]:
        row = _jsonable(self)
        row["record_hash"] = self.record_hash
        return row


@dataclass(frozen=True)
class EligibilityDecision:
    status: EligibilityStatus
    reason_codes: tuple[str, ...]
    counts_toward_qualification: bool

    def to_dict(self) -> dict[str, Any]:
        return _jsonable(self)


@dataclass(frozen=True)
class MinimumSampleRequirement:
    version: str
    provenance: str
    minimum_n: int | None = None
    minimum_positive_n: int | None = None
    minimum_negative_n: int | None = None
    minimum_per_subgroup: int | None = None

    @property
    def configured(self) -> bool:
        return self.minimum_n is not None

    def to_dict(self) -> dict[str, Any]:
        return {**_jsonable(self), "configured": self.configured}


@dataclass(frozen=True)
class NonInferiorityPolicy:
    version: str
    primary_metric: str
    higher_is_better: bool
    margin: float | None
    provenance: str
    permitted_equivalent_benefits: tuple[str, ...] = (
        "BETTER_CALIBRATION", "BETTER_ROBUSTNESS", "BETTER_APPLICABILITY",
        "LOWER_COMPUTATIONAL_COST", "CLEARER_PROVENANCE",
    )

    @property
    def configured(self) -> bool:
        return self.margin is not None

    def to_dict(self) -> dict[str, Any]:
        return {**_jsonable(self), "configured": self.configured}


@dataclass(frozen=True)
class QualificationPolicy:
    endpoint_name: str
    endpoint_id: str
    endpoint_contract_version: str
    qualification_kind: QualificationKind
    qualification_allowed: bool
    current_strategy: StrategyType
    current_policy_version: str
    current_lifecycle: QualificationLifecycle
    output_type: str
    minimum_sample_requirement: MinimumSampleRequirement
    noninferiority_policy: NonInferiorityPolicy
    activation_mode: str = ACTIVATION_MODE
    required_evidence_timing: EvidenceTiming = EvidenceTiming.PROSPECTIVE
    minimum_validation_type: ValidationType = ValidationType.PROSPECTIVE_INTERNAL
    limitations: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return _jsonable(self)


@dataclass(frozen=True)
class QualificationRecord:
    qualification_record_id: str
    endpoint_id: str
    policy_version: str
    candidate_strategy: StrategyType
    candidate_models: tuple[str, ...]
    model_versions: tuple[str, ...]
    checkpoint_hashes: tuple[str, ...]
    current_active_strategy: StrategyType
    validation_dataset: str
    validation_snapshot_hash: str
    validation_type: ValidationType
    prospective_or_retrospective: EvidenceTiming
    sample_size: int
    primary_metrics: tuple[tuple[str, float | int | None], ...]
    secondary_metrics: tuple[tuple[str, float | int | None], ...]
    subgroup_metrics: tuple[tuple[str, str], ...]
    calibration_metrics: tuple[tuple[str, float | int | None], ...]
    applicability_domain_metrics: tuple[tuple[str, float | int | None], ...]
    known_limitations: tuple[str, ...]
    qualification_decision: QualificationDecision
    review_timestamp: datetime
    promotion_status: QualificationLifecycle
    rollback_target: str
    provenance: tuple[tuple[str, str], ...]

    @property
    def record_hash(self) -> str:
        return canonical_hash(self)

    def to_dict(self) -> dict[str, Any]:
        row = _jsonable(self)
        row["record_hash"] = self.record_hash
        return row


@dataclass(frozen=True)
class PerformanceObservation:
    frozen_prediction_id: str
    endpoint_id: str
    candidate_id: str
    model_version_key: str
    policy_version: str
    project_id: str
    chemical_series_id: str
    predicted_value: float
    experimental_value: float
    probability: float | None
    eligibility: EligibilityStatus
    applicability_domain: str = "UNKNOWN"


@dataclass(frozen=True)
class PerformanceSnapshot:
    endpoint_id: str
    candidate_id: str
    task_type: str
    scope: MonitoringScope
    scope_key: str
    model_version_key: str
    policy_version: str
    metrics: tuple[tuple[str, float | int | None], ...]
    source_record_hashes: tuple[str, ...]
    rebuilt_at: datetime

    @property
    def source_snapshot_hash(self) -> str:
        return canonical_hash(self.source_record_hashes)

    @property
    def statistics_hash(self) -> str:
        return canonical_hash({
            "endpoint_id": self.endpoint_id,
            "candidate_id": self.candidate_id,
            "task_type": self.task_type,
            "scope": self.scope.value,
            "scope_key": self.scope_key,
            "model_version_key": self.model_version_key,
            "policy_version": self.policy_version,
            "metrics": self.metrics,
            "source_record_hashes": self.source_record_hashes,
        })

    def to_dict(self) -> dict[str, Any]:
        row = _jsonable(self)
        row["source_snapshot_hash"] = self.source_snapshot_hash
        row["statistics_hash"] = self.statistics_hash
        return row


@dataclass(frozen=True)
class PromotionGateInput:
    endpoint_contract_compatible: bool
    model_identity_frozen: bool
    no_leakage: bool
    independent_validation: bool
    calibration_acceptable: bool
    subgroup_robust: bool
    no_unacceptable_safety_tradeoff: bool
    rollback_target_exists: bool
    artifacts_reproducible: bool
    candidate_metrics: tuple[tuple[str, float | int | None], ...]
    active_metrics: tuple[tuple[str, float | int | None], ...]
    n_positive: int | None = None
    n_negative: int | None = None
    minimum_subgroup_n: int | None = None
    claimed_improvement: bool = True
    equivalent_benefits: tuple[str, ...] = ()


@dataclass(frozen=True)
class PromotionGateDecision:
    passed: bool
    decision: QualificationDecision
    checks: tuple[tuple[str, bool], ...]
    reason_codes: tuple[str, ...]
    activation_mode: str = ACTIVATION_MODE

    def to_dict(self) -> dict[str, Any]:
        return _jsonable(self)


@dataclass(frozen=True)
class DriftPolicy:
    version: str
    minimum_n: int
    performance_absolute_delta: float
    ece_absolute_delta: float
    domain_out_rate_delta: float
    prevalence_absolute_delta: float
    chemical_distance_delta: float
    provenance: str


@dataclass(frozen=True)
class DriftAssessment:
    warnings: tuple[DriftWarning, ...]
    status: str
    review_required: bool
    details: tuple[tuple[str, str], ...]
    policy_version: str
    automatic_action: str = "NONE"

    def to_dict(self) -> dict[str, Any]:
        return _jsonable(self)


@dataclass(frozen=True)
class ManualPromotionAuthorization:
    authorization_id: str
    authorized_by: str
    authorized_at: datetime
    qualification_record_hash: str
    reason: str


@dataclass
class ManagedStrategy:
    specification: CandidateSpecification
    state: QualificationLifecycle
    previous_active_candidate_id: str | None = None


@dataclass(frozen=True)
class LifecycleEvent:
    event_id: str
    candidate_id: str
    endpoint_id: str
    from_state: QualificationLifecycle
    to_state: QualificationLifecycle
    occurred_at: datetime
    state_machine_version: str
    reason: str
    authorization_id: str = ""
    rollback_reason: str = ""


def validate_transition(
    from_state: QualificationLifecycle,
    to_state: QualificationLifecycle,
    *,
    manual_authorization: ManualPromotionAuthorization | None = None,
) -> None:
    if to_state not in LEGAL_TRANSITIONS[from_state]:
        raise ValueError(f"illegal qualification transition: {from_state.value} -> {to_state.value}")
    if to_state == QualificationLifecycle.ACTIVE and manual_authorization is None:
        raise PermissionError("manual promotion authorization is required for activation")


def evaluate_experimental_eligibility(
    frozen: ProspectivePredictionFreeze | None,
    result: ExperimentalQualificationResult,
) -> EligibilityDecision:
    reasons: list[str] = []
    if frozen is None:
        return EligibilityDecision(EligibilityStatus.NO_FROZEN_PREDICTION, ("NO_FROZEN_PREDICTION",), False)
    if frozen.frozen_at >= result.result_available_at:
        return EligibilityDecision(EligibilityStatus.POST_HOC_PREDICTION, ("PREDICTION_NOT_FROZEN_BEFORE_RESULT",), False)
    if not all((result.endpoint_id, result.endpoint_contract_version, result.unit,
                result.assay_type, result.species, result.experiment_date, result.source)):
        return EligibilityDecision(EligibilityStatus.MISSING_METADATA, ("REQUIRED_EXPERIMENT_METADATA_MISSING",), False)
    if result.quality == AssayEvidenceQuality.INCOMPATIBLE:
        return EligibilityDecision(EligibilityStatus.INCOMPATIBLE, ("ASSAY_MARKED_INCOMPATIBLE",), False)
    if frozen.endpoint_id != result.endpoint_id:
        return EligibilityDecision(EligibilityStatus.INCOMPATIBLE, ("ENDPOINT_ID_MISMATCH",), False)
    if frozen.endpoint_contract_version != result.endpoint_contract_version:
        return EligibilityDecision(EligibilityStatus.INCOMPATIBLE, ("ENDPOINT_CONTRACT_VERSION_MISMATCH",), False)

    contract = next((row for row in ENDPOINT_CONTRACTS.values() if row.endpoint_id == result.endpoint_id), None)
    if contract is None:
        reasons.append("ENDPOINT_CONTRACT_NOT_REGISTERED")
    else:
        compatible_types = (
            contract.experimental_compatibility_rules.get("compatible_experimental_types")
            or contract.experimental_compatibility_rules.get("compatible_types")
            or []
        )
        if compatible_types and result.assay_type.casefold() not in {str(item).casefold() for item in compatible_types}:
            reasons.append("ASSAY_TYPE_INCOMPATIBLE")
        incompatible_types = contract.experimental_compatibility_rules.get("incompatible_types") or []
        if result.assay_type.casefold() in {str(item).casefold() for item in incompatible_types}:
            reasons.append("ASSAY_TYPE_INCOMPATIBLE")
        if result.endpoint_id == "safety_herg_blocker_prob" and not any(
            token in result.assay_type.casefold()
            for token in ("herg", "patch clamp", "dofetilide", "astemizole")
        ):
            reasons.append("ASSAY_TYPE_INCOMPATIBLE")
        if contract.output_type == OutputType.BINARY_CLASSIFICATION:
            accepted_units = {"binary", "binary_label", "class", "0/1", contract.canonical_unit.casefold()}
        else:
            accepted_units = {contract.canonical_unit.casefold(), contract.raw_unit.casefold()}
        if result.unit.casefold() not in accepted_units:
            reasons.append("UNIT_INCOMPATIBLE")
        contract_species = contract.species.casefold()
        result_species = result.species.casefold()
        if "chemical" in contract_species:
            if result_species not in {"chemical", "in vitro", "chemical / in vitro"}:
                reasons.append("SPECIES_INCOMPATIBLE")
        elif result_species != contract_species:
            reasons.append("SPECIES_INCOMPATIBLE")

    if reasons:
        return EligibilityDecision(EligibilityStatus.INCOMPATIBLE, tuple(sorted(set(reasons))), False)
    if result.quality == AssayEvidenceQuality.LIMITED:
        return EligibilityDecision(EligibilityStatus.LIMITED, ("LIMITED_ASSAY_QUALITY",), False)
    return EligibilityDecision(EligibilityStatus.QUALIFICATION_ELIGIBLE, ("ALL_COMPATIBILITY_GATES_PASSED",), True)


def _average_ranks(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    position = 0
    while position < len(order):
        end = position
        while end + 1 < len(order) and values[order[end + 1]] == values[order[position]]:
            end += 1
        average = (position + end + 2) / 2.0
        for index in order[position:end + 1]:
            ranks[index] = average
        position = end + 1
    return ranks


def _pearson(left: Sequence[float], right: Sequence[float]) -> float | None:
    if len(left) < 2 or len(left) != len(right):
        return None
    l_mean = sum(left) / len(left)
    r_mean = sum(right) / len(right)
    numerator = sum((a - l_mean) * (b - r_mean) for a, b in zip(left, right))
    denominator = math.sqrt(sum((a - l_mean) ** 2 for a in left) * sum((b - r_mean) ** 2 for b in right))
    return numerator / denominator if denominator else None


def regression_metrics(
    predicted: Sequence[float],
    observed: Sequence[float],
    *,
    log10_fold_metrics: bool = True,
) -> dict[str, float | int | None]:
    if len(predicted) != len(observed):
        raise ValueError("predicted and observed lengths differ")
    n = len(predicted)
    if not n:
        return {"n": 0, "mae": None, "rmse": None, "bias": None, "spearman": None,
                "within_2_fold": None, "within_3_fold": None}
    errors = [float(p) - float(y) for p, y in zip(predicted, observed)]
    absolute = [abs(item) for item in errors]
    return {
        "n": n,
        "mae": sum(absolute) / n,
        "rmse": math.sqrt(sum(item * item for item in errors) / n),
        "bias": sum(errors) / n,
        "spearman": _pearson(_average_ranks(predicted), _average_ranks(observed)),
        "within_2_fold": (sum(item <= math.log10(2.0) for item in absolute) / n if log10_fold_metrics else None),
        "within_3_fold": (sum(item <= math.log10(3.0) for item in absolute) / n if log10_fold_metrics else None),
    }


def classification_metrics(probabilities: Sequence[float], observed: Sequence[int], threshold: float = 0.5) -> dict[str, float | int | None]:
    if len(probabilities) != len(observed):
        raise ValueError("probability and observed lengths differ")
    n = len(probabilities)
    if not n:
        return {key: (0 if key == "n" else None) for key in (
            "n", "mcc", "balanced_accuracy", "auroc", "auprc", "brier", "log_loss",
            "sensitivity", "specificity", "ece",
        )}
    if any(label not in (0, 1) for label in observed):
        raise ValueError("classification observations must be binary")
    clipped = [min(max(float(value), 1e-15), 1.0 - 1e-15) for value in probabilities]
    predicted = [int(value >= threshold) for value in clipped]
    tp = sum(a == 1 and p == 1 for a, p in zip(observed, predicted))
    tn = sum(a == 0 and p == 0 for a, p in zip(observed, predicted))
    fp = sum(a == 0 and p == 1 for a, p in zip(observed, predicted))
    fn = sum(a == 1 and p == 0 for a, p in zip(observed, predicted))
    sensitivity = tp / (tp + fn) if tp + fn else None
    specificity = tn / (tn + fp) if tn + fp else None
    balanced = (sensitivity + specificity) / 2.0 if sensitivity is not None and specificity is not None else None
    denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / denominator if denominator else None
    positives = [p for p, label in zip(clipped, observed) if label == 1]
    negatives = [p for p, label in zip(clipped, observed) if label == 0]
    auroc = None
    if positives and negatives:
        wins = sum((a > b) + 0.5 * (a == b) for a in positives for b in negatives)
        auroc = wins / (len(positives) * len(negatives))
    ranked = sorted(zip(clipped, observed), reverse=True)
    auprc = None
    if positives:
        hit = 0
        precision_sum = 0.0
        for index, (_, label) in enumerate(ranked, start=1):
            if label:
                hit += 1
                precision_sum += hit / index
        auprc = precision_sum / len(positives)
    brier = sum((p - y) ** 2 for p, y in zip(clipped, observed)) / n
    log_loss = -sum(y * math.log(p) + (1 - y) * math.log(1 - p) for p, y in zip(clipped, observed)) / n
    ece = 0.0
    for bin_index in range(10):
        low, high = bin_index / 10.0, (bin_index + 1) / 10.0
        indexes = [i for i, value in enumerate(clipped) if low <= value < high or (bin_index == 9 and value == 1.0)]
        if indexes:
            confidence = sum(clipped[i] for i in indexes) / len(indexes)
            prevalence = sum(observed[i] for i in indexes) / len(indexes)
            ece += len(indexes) / n * abs(confidence - prevalence)
    return {
        "n": n, "mcc": mcc, "balanced_accuracy": balanced, "auroc": auroc,
        "auprc": auprc, "brier": brier, "log_loss": log_loss,
        "sensitivity": sensitivity, "specificity": specificity, "ece": ece,
    }


def rebuild_performance_snapshot(
    observations: Iterable[PerformanceObservation],
    *,
    endpoint_id: str,
    candidate_id: str,
    task_type: str,
    scope: MonitoringScope = MonitoringScope.GLOBAL,
    scope_key: str = "GLOBAL",
    rebuilt_at: datetime | None = None,
) -> PerformanceSnapshot:
    eligible = [row for row in observations if (
        row.endpoint_id == endpoint_id
        and row.candidate_id == candidate_id
        and row.eligibility == EligibilityStatus.QUALIFICATION_ELIGIBLE
    )]
    if scope == MonitoringScope.PROJECT:
        eligible = [row for row in eligible if row.project_id == scope_key]
    elif scope == MonitoringScope.CHEMICAL_SERIES:
        eligible = [row for row in eligible if row.chemical_series_id == scope_key]
    elif scope == MonitoringScope.MODEL_VERSION:
        eligible = [row for row in eligible if row.model_version_key == scope_key]
    elif scope == MonitoringScope.POLICY_VERSION:
        eligible = [row for row in eligible if row.policy_version == scope_key]
    version_keys = {row.model_version_key for row in eligible}
    policy_versions = {row.policy_version for row in eligible}
    if len(version_keys) > 1:
        raise ValueError("model versions must remain isolated in one performance snapshot")
    if len(policy_versions) > 1:
        raise ValueError("policy versions must remain isolated in one performance snapshot")
    if task_type == OutputType.BINARY_CLASSIFICATION.value:
        metrics = classification_metrics(
            [float(row.probability if row.probability is not None else row.predicted_value) for row in eligible],
            [int(row.experimental_value) for row in eligible],
        )
    elif task_type == OutputType.REGRESSION.value:
        contract = next((row for row in ENDPOINT_CONTRACTS.values() if row.endpoint_id == endpoint_id), None)
        log10_fold_metrics = bool(contract and (
            "log10" in contract.canonical_unit.casefold() or "log10" in contract.transformation.casefold()
        ))
        metrics = regression_metrics(
            [row.predicted_value for row in eligible],
            [row.experimental_value for row in eligible],
            log10_fold_metrics=log10_fold_metrics,
        )
    else:
        raise ValueError(f"metrics are not defined for output type {task_type}")
    source_hashes = tuple(sorted(canonical_hash(row) for row in eligible))
    return PerformanceSnapshot(
        endpoint_id=endpoint_id,
        candidate_id=candidate_id,
        task_type=task_type,
        scope=scope,
        scope_key=scope_key,
        model_version_key=next(iter(version_keys), "NO_DATA"),
        policy_version=next(iter(policy_versions), "NO_DATA"),
        metrics=tuple(sorted(metrics.items())),
        source_record_hashes=source_hashes,
        rebuilt_at=rebuilt_at or datetime.now(timezone.utc),
    )


def evaluate_promotion_gate(
    policy: QualificationPolicy,
    evidence: PromotionGateInput,
) -> PromotionGateDecision:
    candidate = dict(evidence.candidate_metrics)
    active = dict(evidence.active_metrics)
    sample = policy.minimum_sample_requirement
    noninferiority = policy.noninferiority_policy
    n = int(candidate.get("n") or 0)
    sample_sufficient = sample.configured and n >= int(sample.minimum_n or 0)
    if sample.minimum_positive_n is not None:
        sample_sufficient = sample_sufficient and (evidence.n_positive or 0) >= sample.minimum_positive_n
    if sample.minimum_negative_n is not None:
        sample_sufficient = sample_sufficient and (evidence.n_negative or 0) >= sample.minimum_negative_n
    subgroup_sufficient = sample.minimum_per_subgroup is None or (
        evidence.minimum_subgroup_n is not None and evidence.minimum_subgroup_n >= sample.minimum_per_subgroup
    )
    metric_name = noninferiority.primary_metric
    candidate_metric, active_metric = candidate.get(metric_name), active.get(metric_name)
    noninferior = False
    if noninferiority.configured and candidate_metric is not None and active_metric is not None:
        margin = float(noninferiority.margin or 0.0)
        noninferior = (
            float(candidate_metric) >= float(active_metric) - margin
            if noninferiority.higher_is_better
            else float(candidate_metric) <= float(active_metric) + margin
        )
    benefit_supported = bool(evidence.claimed_improvement and candidate_metric is not None and active_metric is not None and (
        float(candidate_metric) > float(active_metric)
        if noninferiority.higher_is_better else float(candidate_metric) < float(active_metric)
    ))
    if not evidence.claimed_improvement:
        benefit_supported = bool(set(evidence.equivalent_benefits) & set(noninferiority.permitted_equivalent_benefits))
    checks = {
        "qualification_pipeline_allowed": policy.qualification_allowed,
        "endpoint_contract_compatible": evidence.endpoint_contract_compatible,
        "model_identity_frozen": evidence.model_identity_frozen,
        "no_leakage": evidence.no_leakage,
        "independent_validation": evidence.independent_validation,
        "minimum_sample_requirements_configured": sample.configured,
        "sample_size_and_class_balance_sufficient": sample_sufficient,
        "subgroup_sample_size_sufficient": subgroup_sufficient,
        "noninferiority_margin_configured": noninferiority.configured,
        "primary_metric_noninferior": noninferior,
        "meaningful_or_equivalent_benefit": benefit_supported,
        "calibration_acceptable": evidence.calibration_acceptable,
        "subgroup_robust": evidence.subgroup_robust,
        "no_unacceptable_safety_tradeoff": evidence.no_unacceptable_safety_tradeoff,
        "rollback_target_exists": evidence.rollback_target_exists,
        "artifacts_reproducible": evidence.artifacts_reproducible,
        "manual_activation_only": policy.activation_mode == ACTIVATION_MODE,
    }
    failed = tuple(key.upper() for key, passed in checks.items() if not passed)
    passed = all(checks.values())
    if passed:
        decision = QualificationDecision.QUALIFIED
    elif not sample.configured or not noninferiority.configured or not sample_sufficient:
        decision = QualificationDecision.INSUFFICIENT_EVIDENCE
    else:
        decision = QualificationDecision.NOT_QUALIFIED
    return PromotionGateDecision(passed, decision, tuple(checks.items()), failed)


def detect_drift(
    *,
    policy: DriftPolicy,
    reference_metrics: Mapping[str, float | int | None],
    current_metrics: Mapping[str, float | int | None],
    primary_metric: str,
    higher_is_better: bool,
    reference_model_version: str,
    current_model_version: str,
    reference_endpoint_contract: str,
    current_endpoint_contract: str,
    reference_out_of_domain_rate: float | None = None,
    current_out_of_domain_rate: float | None = None,
    reference_prevalence: float | None = None,
    current_prevalence: float | None = None,
    chemical_distance_delta: float | None = None,
) -> DriftAssessment:
    warnings: list[DriftWarning] = []
    details: dict[str, str] = {}
    if reference_model_version != current_model_version:
        warnings.append(DriftWarning.MODEL_VERSION_CHANGED)
    if reference_endpoint_contract != current_endpoint_contract:
        warnings.append(DriftWarning.ENDPOINT_MISMATCH)
    current_n = int(current_metrics.get("n") or 0)
    if current_n < policy.minimum_n:
        warnings.append(DriftWarning.INSUFFICIENT_DATA)
        details["sample_size"] = f"{current_n} < {policy.minimum_n}"
    reference_value = reference_metrics.get(primary_metric)
    current_value = current_metrics.get(primary_metric)
    if reference_value is not None and current_value is not None and current_n >= policy.minimum_n:
        degradation = (
            float(reference_value) - float(current_value)
            if higher_is_better else float(current_value) - float(reference_value)
        )
        if degradation > policy.performance_absolute_delta:
            warnings.append(DriftWarning.PERFORMANCE_DRIFT)
            details["performance_delta"] = str(degradation)
    reference_ece, current_ece = reference_metrics.get("ece"), current_metrics.get("ece")
    if reference_ece is not None and current_ece is not None and float(current_ece) - float(reference_ece) > policy.ece_absolute_delta:
        warnings.append(DriftWarning.CALIBRATION_DRIFT)
    if reference_out_of_domain_rate is not None and current_out_of_domain_rate is not None:
        if current_out_of_domain_rate - reference_out_of_domain_rate > policy.domain_out_rate_delta:
            warnings.append(DriftWarning.DOMAIN_SHIFT)
    if chemical_distance_delta is not None and chemical_distance_delta > policy.chemical_distance_delta:
        warnings.append(DriftWarning.DOMAIN_SHIFT)
    if reference_prevalence is not None and current_prevalence is not None:
        if abs(current_prevalence - reference_prevalence) > policy.prevalence_absolute_delta:
            warnings.append(DriftWarning.PRIOR_SHIFT)
    unique = tuple(dict.fromkeys(warnings))
    return DriftAssessment(
        warnings=unique,
        status="REVIEW_REQUIRED" if unique else "STABLE",
        review_required=bool(unique),
        details=tuple(sorted(details.items())),
        policy_version=policy.version,
    )


class QualificationLifecycleService:
    """Internal-only lifecycle coordinator; no public mutation route calls it."""

    def __init__(self) -> None:
        self.strategies: dict[str, ManagedStrategy] = {}
        self.active_by_endpoint: dict[str, str] = {}
        self.events: list[LifecycleEvent] = []

    def register(self, specification: CandidateSpecification, state: QualificationLifecycle) -> ManagedStrategy:
        existing = self.strategies.get(specification.candidate_id)
        if existing and existing.specification.specification_hash != specification.specification_hash:
            raise ValueError("candidate identity is frozen; any specification change requires a new candidate_id")
        if existing:
            return existing
        managed = ManagedStrategy(specification=specification, state=state)
        self.strategies[specification.candidate_id] = managed
        if state == QualificationLifecycle.ACTIVE:
            if specification.endpoint_id in self.active_by_endpoint:
                raise ValueError("an endpoint cannot have contradictory ACTIVE strategies")
            self.active_by_endpoint[specification.endpoint_id] = specification.candidate_id
        return managed

    def transition(
        self,
        candidate_id: str,
        to_state: QualificationLifecycle,
        reason: str,
        authorization: ManualPromotionAuthorization | None = None,
        qualification_record: QualificationRecord | None = None,
        gate_decision: PromotionGateDecision | None = None,
    ) -> LifecycleEvent:
        managed = self.strategies[candidate_id]
        validate_transition(managed.state, to_state, manual_authorization=authorization)
        if managed.state == QualificationLifecycle.SHADOW and to_state == QualificationLifecycle.VALIDATED:
            if qualification_record is None or qualification_record.qualification_decision != QualificationDecision.QUALIFIED:
                raise ValueError("a structured QUALIFIED qualification record is required for VALIDATED")
            expected_models = tuple(item.model_id for item in managed.specification.models)
            expected_versions = tuple(item.model_version for item in managed.specification.models)
            if (
                qualification_record.endpoint_id != managed.specification.endpoint_id
                or qualification_record.candidate_strategy != managed.specification.candidate_strategy
                or qualification_record.candidate_models != expected_models
                or qualification_record.model_versions != expected_versions
            ):
                raise ValueError("qualification record identity does not match the frozen candidate")
        if to_state == QualificationLifecycle.PRODUCTION_CANDIDATE:
            if gate_decision is None or not gate_decision.passed:
                raise ValueError("a passing conjunctive promotion gate is required for PRODUCTION_CANDIDATE")
        if to_state == QualificationLifecycle.ACTIVE:
            previous_id = self.active_by_endpoint.get(managed.specification.endpoint_id)
            if previous_id and previous_id != candidate_id:
                previous = self.strategies[previous_id]
                previous.state = QualificationLifecycle.RETIRED
                managed.previous_active_candidate_id = previous_id
            self.active_by_endpoint[managed.specification.endpoint_id] = candidate_id
        event_row = LifecycleEvent(
            event_id=f"LCE-{uuid.uuid4().hex}", candidate_id=candidate_id,
            endpoint_id=managed.specification.endpoint_id, from_state=managed.state,
            to_state=to_state, occurred_at=datetime.now(timezone.utc),
            state_machine_version=STATE_MACHINE_VERSION, reason=reason,
            authorization_id=authorization.authorization_id if authorization else "",
        )
        managed.state = to_state
        self.events.append(event_row)
        return event_row

    def activate(
        self,
        candidate_id: str,
        authorization: ManualPromotionAuthorization,
    ) -> LifecycleEvent:
        return self.transition(candidate_id, QualificationLifecycle.ACTIVE, authorization.reason, authorization)

    def rollback(self, candidate_id: str, reason: RollbackReason) -> tuple[LifecycleEvent, str]:
        managed = self.strategies[candidate_id]
        if managed.state != QualificationLifecycle.ACTIVE:
            raise ValueError("only an ACTIVE candidate can be rolled back")
        target_id = managed.previous_active_candidate_id
        if not target_id or target_id not in self.strategies:
            raise ValueError("deterministic rollback target is missing")
        target = self.strategies[target_id]
        event_row = LifecycleEvent(
            event_id=f"LCE-{uuid.uuid4().hex}", candidate_id=candidate_id,
            endpoint_id=managed.specification.endpoint_id,
            from_state=QualificationLifecycle.ACTIVE,
            to_state=QualificationLifecycle.ROLLED_BACK,
            occurred_at=datetime.now(timezone.utc),
            state_machine_version=STATE_MACHINE_VERSION,
            reason="Deterministic rollback to frozen previous ACTIVE specification",
            rollback_reason=reason.value,
        )
        managed.state = QualificationLifecycle.ROLLED_BACK
        target.state = QualificationLifecycle.ACTIVE
        self.active_by_endpoint[managed.specification.endpoint_id] = target_id
        self.events.append(event_row)
        return event_row, target_id


class QualificationPredictionFreezeRow(Base):
    __tablename__ = "qualification_prediction_freezes"
    frozen_prediction_id: Mapped[str] = mapped_column(String(80), primary_key=True)
    compound_version_id: Mapped[str] = mapped_column(String(80), index=True)
    project_id: Mapped[str] = mapped_column(String(80), default="", index=True)
    chemical_series_id: Mapped[str] = mapped_column(String(120), default="", index=True)
    endpoint_id: Mapped[str] = mapped_column(String(120), index=True)
    endpoint_contract_version: Mapped[str] = mapped_column(String(80))
    candidate_id: Mapped[str] = mapped_column(String(120), index=True)
    candidate_specification_hash: Mapped[str] = mapped_column(String(64))
    strategy: Mapped[str] = mapped_column(String(80))
    models_json: Mapped[list] = mapped_column(JSON)
    prediction_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    probability: Mapped[float | None] = mapped_column(Float, nullable=True)
    unit: Mapped[str] = mapped_column(String(60))
    frozen_at: Mapped[datetime] = mapped_column(DateTime(timezone=True))
    policy_version: Mapped[str] = mapped_column(String(100), index=True)
    standardizer_version: Mapped[str] = mapped_column(String(100))
    applicability_domain: Mapped[str] = mapped_column(String(60), default="UNKNOWN")
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)
    record_hash: Mapped[str] = mapped_column(String(64), unique=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class QualificationExperimentalResultRow(Base):
    __tablename__ = "qualification_experimental_results"
    experimental_result_id: Mapped[str] = mapped_column(String(80), primary_key=True)
    frozen_prediction_id: Mapped[str] = mapped_column(String(80), index=True)
    endpoint_id: Mapped[str] = mapped_column(String(120), index=True)
    endpoint_contract_version: Mapped[str] = mapped_column(String(80))
    experimental_value: Mapped[float] = mapped_column(Float)
    unit: Mapped[str] = mapped_column(String(60))
    assay_type: Mapped[str] = mapped_column(String(240))
    species: Mapped[str] = mapped_column(String(120))
    experiment_date: Mapped[str] = mapped_column(String(40))
    result_available_at: Mapped[datetime] = mapped_column(DateTime(timezone=True))
    source: Mapped[str] = mapped_column(Text)
    quality: Mapped[str] = mapped_column(String(60))
    protocol_metadata_json: Mapped[dict] = mapped_column(JSON, default=dict)
    eligibility_status: Mapped[str] = mapped_column(String(60), index=True)
    eligibility_reasons_json: Mapped[list] = mapped_column(JSON, default=list)
    counts_toward_qualification: Mapped[bool] = mapped_column(Boolean, default=False)
    recorded_at: Mapped[datetime] = mapped_column(DateTime(timezone=True))
    record_hash: Mapped[str] = mapped_column(String(64), unique=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class StrategyQualificationRecordRow(Base):
    __tablename__ = "strategy_qualification_records"
    qualification_record_id: Mapped[str] = mapped_column(String(80), primary_key=True)
    endpoint_id: Mapped[str] = mapped_column(String(120), index=True)
    policy_version: Mapped[str] = mapped_column(String(100), index=True)
    candidate_strategy: Mapped[str] = mapped_column(String(80))
    candidate_models_json: Mapped[list] = mapped_column(JSON)
    model_versions_json: Mapped[list] = mapped_column(JSON)
    checkpoint_hashes_json: Mapped[list] = mapped_column(JSON)
    current_active_strategy: Mapped[str] = mapped_column(String(80))
    validation_dataset: Mapped[str] = mapped_column(Text)
    validation_snapshot_hash: Mapped[str] = mapped_column(String(64))
    validation_type: Mapped[str] = mapped_column(String(80))
    prospective_or_retrospective: Mapped[str] = mapped_column(String(40))
    sample_size: Mapped[int] = mapped_column(Integer)
    primary_metrics_json: Mapped[dict] = mapped_column(JSON)
    secondary_metrics_json: Mapped[dict] = mapped_column(JSON)
    subgroup_metrics_json: Mapped[dict] = mapped_column(JSON)
    calibration_metrics_json: Mapped[dict] = mapped_column(JSON)
    applicability_domain_metrics_json: Mapped[dict] = mapped_column(JSON)
    known_limitations_json: Mapped[list] = mapped_column(JSON)
    qualification_decision: Mapped[str] = mapped_column(String(80))
    review_timestamp: Mapped[datetime] = mapped_column(DateTime(timezone=True))
    promotion_status: Mapped[str] = mapped_column(String(80))
    rollback_target: Mapped[str] = mapped_column(String(160))
    provenance_json: Mapped[dict] = mapped_column(JSON)
    record_hash: Mapped[str] = mapped_column(String(64), unique=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class QualificationLifecycleEventRow(Base):
    __tablename__ = "qualification_lifecycle_events"
    event_id: Mapped[str] = mapped_column(String(80), primary_key=True)
    candidate_id: Mapped[str] = mapped_column(String(120), index=True)
    endpoint_id: Mapped[str] = mapped_column(String(120), index=True)
    from_state: Mapped[str] = mapped_column(String(80))
    to_state: Mapped[str] = mapped_column(String(80))
    occurred_at: Mapped[datetime] = mapped_column(DateTime(timezone=True))
    state_machine_version: Mapped[str] = mapped_column(String(100))
    reason: Mapped[str] = mapped_column(Text)
    authorization_id: Mapped[str] = mapped_column(String(120), default="")
    rollback_reason: Mapped[str] = mapped_column(String(80), default="")
    record_hash: Mapped[str] = mapped_column(String(64), unique=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


_APPEND_ONLY_MODELS = (
    QualificationPredictionFreezeRow,
    QualificationExperimentalResultRow,
    StrategyQualificationRecordRow,
    QualificationLifecycleEventRow,
)


def _reject_mutation(mapper: Any, connection: Any, target: Any) -> None:
    del mapper, connection, target
    raise ValueError("Stage 4D-5 qualification evidence is append-only and immutable")


for _model in _APPEND_ONLY_MODELS:
    event.listen(_model, "before_update", _reject_mutation)
    event.listen(_model, "before_delete", _reject_mutation)


def ensure_qualification_schema(engine: Any) -> None:
    if "projects" not in inspect(engine).get_table_names():
        return
    Base.metadata.create_all(bind=engine, tables=[row.__table__ for row in _APPEND_ONLY_MODELS])


class QualificationEvidenceStore:
    """Append-only SQL store used by internal qualification workflows."""

    def __init__(self, session: Session):
        self.session = session

    def freeze_prediction(self, frozen: ProspectivePredictionFreeze) -> QualificationPredictionFreezeRow:
        row = QualificationPredictionFreezeRow(
            frozen_prediction_id=frozen.frozen_prediction_id,
            compound_version_id=frozen.compound_version_id,
            project_id=frozen.project_id,
            chemical_series_id=frozen.chemical_series_id,
            endpoint_id=frozen.endpoint_id,
            endpoint_contract_version=frozen.endpoint_contract_version,
            candidate_id=frozen.candidate_id,
            candidate_specification_hash=frozen.candidate_specification_hash,
            strategy=frozen.strategy.value,
            models_json=[_jsonable(item) for item in frozen.models],
            prediction_value=frozen.prediction_value,
            probability=frozen.probability,
            unit=frozen.unit,
            frozen_at=frozen.frozen_at,
            policy_version=frozen.policy_version,
            standardizer_version=frozen.standardizer_version,
            applicability_domain=frozen.applicability_domain,
            provenance_json=dict(frozen.provenance),
            record_hash=frozen.record_hash,
        )
        self.session.add(row)
        self.session.flush()
        return row

    def link_experimental_result(
        self,
        result: ExperimentalQualificationResult,
    ) -> tuple[QualificationExperimentalResultRow, EligibilityDecision]:
        freeze_row = self.session.get(QualificationPredictionFreezeRow, result.frozen_prediction_id)
        frozen = _freeze_from_row(freeze_row) if freeze_row else None
        eligibility = evaluate_experimental_eligibility(frozen, result)
        row = QualificationExperimentalResultRow(
            experimental_result_id=result.experimental_result_id,
            frozen_prediction_id=result.frozen_prediction_id,
            endpoint_id=result.endpoint_id,
            endpoint_contract_version=result.endpoint_contract_version,
            experimental_value=result.experimental_value,
            unit=result.unit,
            assay_type=result.assay_type,
            species=result.species,
            experiment_date=result.experiment_date,
            result_available_at=result.result_available_at,
            source=result.source,
            quality=result.quality.value,
            protocol_metadata_json=dict(result.protocol_metadata),
            eligibility_status=eligibility.status.value,
            eligibility_reasons_json=list(eligibility.reason_codes),
            counts_toward_qualification=eligibility.counts_toward_qualification,
            recorded_at=result.recorded_at,
            record_hash=result.record_hash,
        )
        self.session.add(row)
        self.session.flush()
        return row, eligibility

    def append_qualification_record(self, record: QualificationRecord) -> StrategyQualificationRecordRow:
        row = StrategyQualificationRecordRow(
            qualification_record_id=record.qualification_record_id,
            endpoint_id=record.endpoint_id,
            policy_version=record.policy_version,
            candidate_strategy=record.candidate_strategy.value,
            candidate_models_json=list(record.candidate_models),
            model_versions_json=list(record.model_versions),
            checkpoint_hashes_json=list(record.checkpoint_hashes),
            current_active_strategy=record.current_active_strategy.value,
            validation_dataset=record.validation_dataset,
            validation_snapshot_hash=record.validation_snapshot_hash,
            validation_type=record.validation_type.value,
            prospective_or_retrospective=record.prospective_or_retrospective.value,
            sample_size=record.sample_size,
            primary_metrics_json=dict(record.primary_metrics),
            secondary_metrics_json=dict(record.secondary_metrics),
            subgroup_metrics_json=dict(record.subgroup_metrics),
            calibration_metrics_json=dict(record.calibration_metrics),
            applicability_domain_metrics_json=dict(record.applicability_domain_metrics),
            known_limitations_json=list(record.known_limitations),
            qualification_decision=record.qualification_decision.value,
            review_timestamp=record.review_timestamp,
            promotion_status=record.promotion_status.value,
            rollback_target=record.rollback_target,
            provenance_json=dict(record.provenance),
            record_hash=record.record_hash,
        )
        self.session.add(row)
        self.session.flush()
        return row

    def append_lifecycle_event(self, lifecycle_event: LifecycleEvent) -> QualificationLifecycleEventRow:
        row = QualificationLifecycleEventRow(
            event_id=lifecycle_event.event_id,
            candidate_id=lifecycle_event.candidate_id,
            endpoint_id=lifecycle_event.endpoint_id,
            from_state=lifecycle_event.from_state.value,
            to_state=lifecycle_event.to_state.value,
            occurred_at=lifecycle_event.occurred_at,
            state_machine_version=lifecycle_event.state_machine_version,
            reason=lifecycle_event.reason,
            authorization_id=lifecycle_event.authorization_id,
            rollback_reason=lifecycle_event.rollback_reason,
            record_hash=canonical_hash(lifecycle_event),
        )
        self.session.add(row)
        self.session.flush()
        return row


def _freeze_from_row(row: QualificationPredictionFreezeRow) -> ProspectivePredictionFreeze:
    return ProspectivePredictionFreeze(
        frozen_prediction_id=row.frozen_prediction_id,
        compound_version_id=row.compound_version_id,
        project_id=row.project_id,
        chemical_series_id=row.chemical_series_id,
        endpoint_id=row.endpoint_id,
        endpoint_contract_version=row.endpoint_contract_version,
        candidate_id=row.candidate_id,
        candidate_specification_hash=row.candidate_specification_hash,
        strategy=StrategyType(row.strategy),
        models=tuple(FrozenModelIdentity(**item) for item in row.models_json),
        prediction_value=row.prediction_value,
        probability=row.probability,
        unit=row.unit,
        frozen_at=row.frozen_at if row.frozen_at.tzinfo else row.frozen_at.replace(tzinfo=timezone.utc),
        policy_version=row.policy_version,
        standardizer_version=row.standardizer_version,
        applicability_domain=row.applicability_domain,
        provenance=tuple(sorted((row.provenance_json or {}).items())),
    )


def _default_output_type(endpoint_name: str, strategy: StrategyType) -> str:
    contract = ENDPOINT_CONTRACTS.get(endpoint_name)
    if contract:
        return contract.output_type.value
    if strategy == StrategyType.RANK_FUSION:
        return OutputType.RANKING.value
    if strategy in {StrategyType.MECHANISTIC_NO_CONSENSUS, StrategyType.DERIVED_ESTIMATE, StrategyType.RULE_ESTIMATE}:
        return OutputType.MECHANISTIC_DERIVED.value
    return "NOT_APPLICABLE"


def _metric_policy(output_type: str) -> tuple[str, bool]:
    if output_type == OutputType.BINARY_CLASSIFICATION.value:
        return "mcc", True
    if output_type == OutputType.REGRESSION.value:
        return "mae", False
    if output_type == OutputType.RANKING.value:
        return "top_k_recall", True
    return "method_error", False


def build_qualification_policy_registry() -> dict[str, QualificationPolicy]:
    policies: dict[str, QualificationPolicy] = {}
    for name, governed in ENDPOINT_STRATEGY_REGISTRY.items():
        unavailable = governed.primary_strategy == StrategyType.MODEL_UNAVAILABLE
        mechanistic = governed.primary_strategy in {
            StrategyType.MECHANISTIC_NO_CONSENSUS,
            StrategyType.DERIVED_ESTIMATE,
            StrategyType.RULE_ESTIMATE,
            StrategyType.RULE_BASED,
        }
        kind = QualificationKind.EXCLUDED if unavailable else (
            QualificationKind.METHOD_QUALIFICATION if mechanistic else QualificationKind.STRATEGY_QUALIFICATION
        )
        output_type = _default_output_type(name, governed.primary_strategy)
        metric, higher = _metric_policy(output_type)
        policies[name] = QualificationPolicy(
            endpoint_name=name,
            endpoint_id=governed.endpoint_id,
            endpoint_contract_version=governed.endpoint_contract_version,
            qualification_kind=kind,
            qualification_allowed=not unavailable,
            current_strategy=governed.primary_strategy,
            current_policy_version=governed.policy_version,
            current_lifecycle=(QualificationLifecycle.RESEARCH_ONLY if unavailable else QualificationLifecycle.ACTIVE),
            output_type=output_type,
            minimum_sample_requirement=MinimumSampleRequirement(
                version=f"{QUALIFICATION_POLICY_VERSION}:{governed.endpoint_id}:sample-v1",
                provenance=(
                    "Not configured: an endpoint scientific owner must approve endpoint-, class-, and subgroup-specific minimum counts before qualification."
                    if not unavailable else "Excluded because MODEL_UNAVAILABLE has no candidate."
                ),
            ),
            noninferiority_policy=NonInferiorityPolicy(
                version=f"{QUALIFICATION_POLICY_VERSION}:{governed.endpoint_id}:ni-v1",
                primary_metric=metric,
                higher_is_better=higher,
                margin=None,
                provenance=(
                    "Not configured: an endpoint scientific owner must approve a versioned practical-equivalence margin before qualification; no regulatory claim is implied."
                    if not unavailable else "Excluded because MODEL_UNAVAILABLE has no candidate."
                ),
            ),
            limitations=tuple(governed.limitations),
        )
    return policies


QUALIFICATION_POLICY_REGISTRY = build_qualification_policy_registry()


CANDIDATE_TRACKS: dict[str, dict[str, Any]] = {
    "Solubility": {
        "candidate_id": "solubility-adaptive-shadow-stage4d5-v1",
        "state": QualificationLifecycle.SHADOW.value,
        "candidate_strategy": StrategyType.ADAPTIVE_RESEARCH_SHADOW.value,
        "candidate_models": ["admetica_solubility", "esol_delaney_v1"],
        "excluded_models": ["rdkit_gbr_solubility_v1"],
        "qualification_status": "AWAITING_COMPATIBLE_PROSPECTIVE_EVIDENCE",
        "required_metrics": ["n", "mae", "rmse", "bias", "spearman", "within_2_fold", "within_3_fold"],
        "automatic_activation": False,
    },
    "Permeability": {
        "candidate_id": "caco2-consensus-shadow-stage4d5-v1",
        "state": QualificationLifecycle.SHADOW.value,
        "candidate_strategy": StrategyType.STATIC_CONSENSUS.value,
        "candidate_models": ["admetica_caco2", "physchem_caco2_v1"],
        "qualification_status": "INSUFFICIENT_EVIDENCE",
        "required_metrics": ["n", "mae", "rmse", "bias", "spearman", "within_2_fold", "within_3_fold"],
        "automatic_activation": False,
    },
    "CYP3A4 inhibitor": {
        "candidate_id": "cyp3a4-fixed-9578-0422-shadow-stage4d5-v1",
        "state": QualificationLifecycle.SHADOW.value,
        "candidate_strategy": StrategyType.FIXED_WEIGHT_BLEND.value,
        "candidate_models": ["admetica_cyp_cyp3a4-inhibitor", "morgan_cyp3a4_inh_v1"],
        "weights": [0.9578, 0.0422],
        "dynamic_adaptive_promotion": "CLOSED_NO_ADAPTIVE_VALUE",
        "qualification_status": "AWAITING_COMPATIBLE_PROSPECTIVE_EVIDENCE",
        "required_metrics": ["n", "mcc", "balanced_accuracy", "auroc", "auprc", "brier", "log_loss", "sensitivity", "specificity", "ece"],
        "automatic_activation": False,
    },
    "hERG liability": {
        "candidate_id": "herg-platt-m1-shadow-stage4d5-v1",
        "state": QualificationLifecycle.SHADOW.value,
        "candidate_strategy": StrategyType.SINGLE_CORE_WITH_CALIBRATION.value,
        "candidate_models": ["admetica_safety_herg"],
        "calibration": "PLATT_CALIBRATION_RESEARCH",
        "supporting_only_models": ["physchem_herg_v1"],
        "qualification_status": "AWAITING_COMPATIBLE_PROSPECTIVE_EVIDENCE",
        "required_metrics": ["n", "mcc", "balanced_accuracy", "auroc", "auprc", "brier", "log_loss", "sensitivity", "specificity", "ece"],
        "required_subgroups": ["chemical_series", "applicability_domain", "class_balance"],
        "automatic_activation": False,
    },
}


DEFAULT_DRIFT_POLICY = DriftPolicy(
    version=DRIFT_POLICY_VERSION,
    minimum_n=30,
    performance_absolute_delta=0.10,
    ece_absolute_delta=0.05,
    domain_out_rate_delta=0.10,
    prevalence_absolute_delta=0.10,
    chemical_distance_delta=0.10,
    provenance=(
        "Stage 4D-5 deterministic REVIEW_REQUIRED trigger defaults. These are monitoring thresholds, "
        "not regulatory acceptance criteria, and cannot retrain or change production policy."
    ),
)


def get_production_baseline() -> dict[str, Any]:
    active = []
    for name, governed in sorted(ENDPOINT_STRATEGY_REGISTRY.items()):
        if governed.promotion_status != PromotionStatus.ACTIVE:
            continue
        active.append({
            "endpoint": name,
            "endpoint_id": governed.endpoint_id,
            "endpoint_contract_version": governed.endpoint_contract_version,
            "primary_strategy": governed.primary_strategy.value,
            "models": [
                {"model_id": model_id, "model_version": version}
                for model_id, version in zip(governed.primary_model_ids, governed.primary_model_versions)
            ],
            "decision_threshold": governed.decision_threshold,
            "calibration_status": governed.calibration_status.value,
            "calibration_production_enabled": governed.calibration_production_enabled,
            "policy_version": governed.policy_version,
            "rollback_policy": _jsonable(governed.rollback_policy),
        })
    return {
        "artifact": "stage4d5_production_baseline",
        "baseline_commit": "083163c6362fd2d2d10fb5a473468a1429a6a274",
        "baseline_tag": "stage4d4-endpoint-strategy-finalized",
        "created_for_policy": QUALIFICATION_POLICY_VERSION,
        "active_policy_count": len(active),
        "active_endpoint_policies": active,
        "production_behavior_changed": False,
    }


def get_strategy_cards() -> list[dict[str, Any]]:
    cards = []
    for name, governed in sorted(ENDPOINT_STRATEGY_REGISTRY.items()):
        policy = QUALIFICATION_POLICY_REGISTRY[name]
        contract = ENDPOINT_CONTRACTS.get(name)
        model_spec = MODEL_SPECS.get(name, {})
        evidence_artifact = governed.rollback_policy.validation_artifact if governed.rollback_policy else None
        cards.append({
            "endpoint": name,
            "endpoint_id": governed.endpoint_id,
            "endpoint_contract_version": governed.endpoint_contract_version,
            "scientific_purpose": contract.scientific_definition if contract else governed.scientific_notes,
            "endpoint_definition": contract.to_dict() if contract else None,
            "models": [
                {"model_id": model_id, "model_version": version}
                for model_id, version in zip(governed.primary_model_ids, governed.primary_model_versions)
            ],
            "training_source": {
                "dataset": model_spec.get("training_dataset"),
                "source": model_spec.get("source"),
                "license": model_spec.get("license"),
                "stage_evidence": governed.evidence_stage,
                "stage4d5_training_performed": False,
            },
            "validation_source": {
                "model_registry_validation": model_spec.get("validation"),
                "independent_validation": model_spec.get("independent_validation"),
                "evidence_artifact": evidence_artifact,
            },
            "known_limitations": list(governed.limitations),
            "applicability": governed.applicability_policy,
            "calibration_status": governed.calibration_status.value,
            "production_status": governed.promotion_status.value,
            "production_strategy": governed.primary_strategy.value,
            "shadow_status": governed.shadow_promotion_status.value if governed.shadow_strategy else None,
            "shadow_strategy": governed.shadow_strategy.value if governed.shadow_strategy else None,
            "qualification_kind": policy.qualification_kind.value,
            "qualification_allowed": policy.qualification_allowed,
            "candidate_track": CANDIDATE_TRACKS.get(name),
            "promotion_history": [],
            "rollback_history": [],
            "unsupported_confidence_claims": False,
        })
    return cards


def get_qualification_api_response() -> dict[str, Any]:
    return {
        "policy_version": QUALIFICATION_POLICY_VERSION,
        "state_machine_version": STATE_MACHINE_VERSION,
        "activation_mode": ACTIVATION_MODE,
        "read_only": True,
        "automatic_shadow_activation": False,
        "automatic_retraining": False,
        "strategy_count": len(QUALIFICATION_POLICY_REGISTRY),
        "strategies": [policy.to_dict() for _, policy in sorted(QUALIFICATION_POLICY_REGISTRY.items())],
    }


def get_qualification_endpoint_response(endpoint_id: str) -> dict[str, Any] | None:
    for name, policy in QUALIFICATION_POLICY_REGISTRY.items():
        if endpoint_id in {name, policy.endpoint_id}:
            governed = ENDPOINT_STRATEGY_REGISTRY[name]
            return {
                "read_only": True,
                "policy": policy.to_dict(),
                "active_strategy": governed.to_dict(),
                "candidate_track": CANDIDATE_TRACKS.get(name),
                "activation_mode": ACTIVATION_MODE,
            }
    return None


def get_candidates_api_response() -> dict[str, Any]:
    return {
        "read_only": True,
        "activation_mode": ACTIVATION_MODE,
        "automatic_activation": False,
        "candidate_count": len(CANDIDATE_TRACKS),
        "candidates": [
            {"endpoint": name, "endpoint_id": QUALIFICATION_POLICY_REGISTRY[name].endpoint_id, **track}
            for name, track in sorted(CANDIDATE_TRACKS.items())
        ],
    }


def get_drift_api_response() -> dict[str, Any]:
    return {
        "read_only": True,
        "policy": _jsonable(DEFAULT_DRIFT_POLICY),
        "automatic_action": "NONE",
        "review_action": "REVIEW_REQUIRED",
        "endpoints": [
            {
                "endpoint": name,
                "endpoint_id": policy.endpoint_id,
                "status": DriftWarning.INSUFFICIENT_DATA.value,
                "reason": "No Stage 4D-5 prospective qualification observations have been accumulated.",
            }
            for name, policy in sorted(QUALIFICATION_POLICY_REGISTRY.items())
            if policy.qualification_allowed
        ],
    }


def validate_qualification_registry() -> list[str]:
    violations: list[str] = []
    if set(QUALIFICATION_POLICY_REGISTRY) != set(ENDPOINT_STRATEGY_REGISTRY):
        violations.append("qualification policies do not cover every endpoint strategy")
    for name, policy in QUALIFICATION_POLICY_REGISTRY.items():
        governed = ENDPOINT_STRATEGY_REGISTRY[name]
        if governed.primary_strategy == StrategyType.MODEL_UNAVAILABLE:
            if policy.qualification_allowed or policy.qualification_kind != QualificationKind.EXCLUDED:
                violations.append(f"{name}: MODEL_UNAVAILABLE entered qualification pipeline")
        if governed.primary_strategy == StrategyType.MECHANISTIC_NO_CONSENSUS:
            if policy.qualification_kind != QualificationKind.METHOD_QUALIFICATION:
                violations.append(f"{name}: mechanistic endpoint lacks METHOD_QUALIFICATION")
        if policy.activation_mode != ACTIVATION_MODE:
            violations.append(f"{name}: activation is not manual")
    for name, track in CANDIDATE_TRACKS.items():
        if track["state"] != QualificationLifecycle.SHADOW.value or track["automatic_activation"]:
            violations.append(f"{name}: candidate track is not safely frozen in SHADOW")
    if CANDIDATE_TRACKS["hERG liability"]["candidate_models"] != ["admetica_safety_herg"]:
        violations.append("hERG calibration candidate must be calibrated M1 only")
    if CANDIDATE_TRACKS["CYP3A4 inhibitor"]["dynamic_adaptive_promotion"] != "CLOSED_NO_ADAPTIVE_VALUE":
        violations.append("CYP3A4 dynamic adaptive promotion must remain closed")
    return violations
