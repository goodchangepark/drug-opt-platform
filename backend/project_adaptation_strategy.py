"""Endpoint-aware adaptation overlays above the frozen Engine-v1 output."""
from __future__ import annotations

from dataclasses import dataclass
from statistics import mean

from .project_adaptation_v2 import (
    ENGINE_V1_HASH, ENGINE_V1_POLICY, QualifiedEvidencePair, fit_project_adapter,
)

BASE_ONLY = "BASE_ONLY"
MULTI_MODEL_ENSEMBLE_ADAPTATION = "MULTI_MODEL_ENSEMBLE_ADAPTATION"
SINGLE_MODEL_RESIDUAL_CALIBRATION = "SINGLE_MODEL_RESIDUAL_CALIBRATION"
LOCAL_SERIES_ADAPTATION_RESEARCH = "LOCAL_SERIES_ADAPTATION_RESEARCH"
UNSUPPORTED_FOR_ADAPTATION = "UNSUPPORTED_FOR_ADAPTATION"


@dataclass(frozen=True)
class AdaptationStrategy:
    strategy_type: str
    model_count: int
    reason: str


@dataclass
class StrategyResult:
    endpoint_id: str
    strategy_type: str
    model_count: int
    status: str
    raw_n: int
    effective_n: float
    global_weights: dict
    project_weights: dict
    base_validation_error: float | None
    adapted_validation_error: float | None
    observed_bias: float | None = None
    shrinkage_factor: float = 0.0
    calibration_adjustment: float = 0.0
    calibration_scale: str = ""
    holdouts_improved: float | None = None
    stability: str = "INSUFFICIENT_EVIDENCE"
    activation_decision: str = "BASE_RETAINED"
    reason: str = ""

    def to_dict(self):
        return {
            **vars(self), "base_engine_policy": ENGINE_V1_POLICY, "base_engine_hash": ENGINE_V1_HASH,
            "strategy": self.strategy_type,
        }


def resolve_adaptation_strategy(model_count: int, endpoint_id: str = "", prediction_type: str = "REGRESSION") -> AdaptationStrategy:
    if model_count <= 0:
        return AdaptationStrategy(BASE_ONLY, 0, "No active production model outputs")
    if str(prediction_type).lower() in {"binary_classification", "classification", "probability"}:
        return AdaptationStrategy(UNSUPPORTED_FOR_ADAPTATION, model_count, "Classification calibration is not enabled for small-N overlays")
    if model_count >= 2:
        return AdaptationStrategy(MULTI_MODEL_ENSEMBLE_ADAPTATION, model_count, "At least two active model outputs are available")
    return AdaptationStrategy(SINGLE_MODEL_RESIDUAL_CALIBRATION, 1, "One active production model; residual calibration is evaluated instead of weight fitting")


def _eligible(events):
    return [event for event in events if event.eligible]


def _single_model(events: list[QualifiedEvidencePair], endpoint_id: str, global_weights: dict) -> StrategyResult:
    events = _eligible(events)
    model_ids = list(global_weights)
    if len(model_ids) != 1:
        return StrategyResult(endpoint_id, SINGLE_MODEL_RESIDUAL_CALIBRATION, len(model_ids), "BASE_ONLY", len(events), 0.0, global_weights, global_weights, None, None, reason="Single-model strategy requires exactly one model")
    model = model_ids[0]
    if len(events) < 5:
        return StrategyResult(endpoint_id, SINGLE_MODEL_RESIDUAL_CALIBRATION, 1, "BASE_ONLY", len(events), float(len(events)), global_weights, global_weights, None, None, calibration_scale="endpoint canonical scale", reason="INSUFFICIENT_EVIDENCE")
    residuals = [event.value - event.frozen_predictions[model] for event in events if model in event.frozen_predictions]
    if len(residuals) < 5:
        return StrategyResult(endpoint_id, SINGLE_MODEL_RESIDUAL_CALIBRATION, 1, "BASE_ONLY", len(events), float(len(residuals)), global_weights, global_weights, None, None, reason="Missing frozen model output")
    bias = mean(residuals)
    same_sign = max(sum(value >= 0 for value in residuals), sum(value <= 0 for value in residuals)) / len(residuals)
    # Ridge/shrinkage: the observed bias is estimated from all points, then
    # pulled toward zero. The strength grows smoothly with independent N.
    shrinkage = min(0.55, len(residuals) / (len(residuals) + 8.0))
    correction = shrinkage * bias
    base_errors, adapted_errors = [], []
    for index, held_out in enumerate(events):
        train = events[:index] + events[index + 1:]
        train_residuals = [item.value - item.frozen_predictions[model] for item in train if model in item.frozen_predictions]
        train_bias = mean(train_residuals) if train_residuals else 0.0
        train_shrinkage = min(0.55, len(train_residuals) / (len(train_residuals) + 8.0))
        base = held_out.frozen_predictions[model]
        adapted = base + train_shrinkage * train_bias
        base_errors.append(abs(base - held_out.value)); adapted_errors.append(abs(adapted - held_out.value))
    base_error = mean(base_errors) if base_errors else None
    adapted_error = mean(adapted_errors) if adapted_errors else None
    improved = sum(a < b for a, b in zip(adapted_errors, base_errors)) / len(base_errors) if base_errors else None
    stable = same_sign >= 0.6 and abs(bias) > 1e-12
    activate = stable and adapted_error is not None and base_error is not None and adapted_error < base_error - 1e-12
    if activate:
        return StrategyResult(endpoint_id, SINGLE_MODEL_RESIDUAL_CALIBRATION, 1, "CANDIDATE_VALIDATED_IMPROVEMENT", len(events), float(len(events)), global_weights, global_weights, base_error, adapted_error, bias, shrinkage, correction, "endpoint canonical scale", improved, "STABLE_BIAS", "ACTIVATED", "LOO demonstrated improvement")
    reason = "NO_STABLE_PROJECT_BIAS" if not stable else "NO_IMPROVEMENT"
    return StrategyResult(endpoint_id, SINGLE_MODEL_RESIDUAL_CALIBRATION, 1, "BASE_RETAINED_NO_STABLE_PROJECT_BIAS" if reason == "NO_STABLE_PROJECT_BIAS" else "NO_IMPROVEMENT_BASE_RETAINED", len(events), float(len(events)), global_weights, global_weights, base_error, adapted_error, bias, shrinkage, 0.0, "endpoint canonical scale", improved, "UNSTABLE_BIAS" if not stable else "STABLE_BIAS_NO_GAIN", "BASE_RETAINED", reason)


def fit_project_adaptation_strategy(endpoint_id: str, evidence, global_weights: dict, *, prediction_type: str = "REGRESSION") -> StrategyResult:
    weights = {key: float(value) for key, value in (global_weights or {}).items()}
    strategy = resolve_adaptation_strategy(len(weights), endpoint_id, prediction_type)
    if strategy.strategy_type == SINGLE_MODEL_RESIDUAL_CALIBRATION:
        return _single_model(list(evidence), endpoint_id, weights)
    if strategy.strategy_type == MULTI_MODEL_ENSEMBLE_ADAPTATION:
        result = fit_project_adapter(endpoint_id, evidence, weights)
        return StrategyResult(endpoint_id, strategy.strategy_type, len(weights), result.status, result.raw_n, result.effective_n,
                              result.global_weights, result.project_weights, result.base_validation_error, result.adapted_validation_error,
                              activation_decision=result.activation_decision,
                              reason="Existing regularized non-negative ensemble adapter")
    return StrategyResult(endpoint_id, strategy.strategy_type, len(weights), "BASE_ONLY", 0, 0.0, weights, weights, None, None, reason=strategy.reason)
