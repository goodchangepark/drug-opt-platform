"""Conservative project-specific adapter above frozen Engine v1.

The inputs are pre-existing frozen model predictions plus subsequently
qualified evidence.  No model parameters or historical predictions are ever
changed by this module.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

from .adaptive_weighting import compute_morgan_fingerprint, compute_tanimoto_similarity

ADAPTER_POLICY_VERSION = "drugopt-project-adapter-v2"
ENGINE_V1_POLICY = "drugopt-prediction-engine-v1@1.0.0"
ENGINE_V1_HASH = "12757ab197b5a70d8ea1754678d9a342ab0b6ea0d82f2896bebb767d686bbdeb"


@dataclass(frozen=True)
class QualifiedEvidencePair:
    evidence_id: str
    compound_version_id: int
    smiles: str
    endpoint_id: str
    value: float
    frozen_predictions: dict[str, float]
    origin: str = "EXPERIMENTAL_INTERNAL"
    source_quality: str = "A"
    comparability_status: str = "DIRECTLY_COMPARABLE"
    duplicate_status: str = "DISTINCT_MEASUREMENT"
    # Optional validation metadata; production callers remain backward compatible.
    series_id: str = ""

    @property
    def quality_weight(self) -> float:
        if self.origin == "EXPERIMENTAL_INTERNAL": return 1.0
        return 0.7 if self.source_quality == "A" else 0.5 if self.source_quality == "B" else 0.0

    @property
    def eligible(self) -> bool:
        return self.quality_weight > 0 and self.comparability_status in {"DIRECTLY_COMPARABLE", "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION"} and self.duplicate_status == "DISTINCT_MEASUREMENT"


@dataclass
class ProjectAdapterResult:
    endpoint_id: str
    status: str
    raw_n: int
    effective_n: float
    global_weights: dict[str, float]
    project_weights: dict[str, float]
    beta: float
    base_validation_error: float | None
    adapted_validation_error: float | None
    activation_decision: str
    adapter_version: str = ADAPTER_POLICY_VERSION

    def to_dict(self):
        return vars(self) | {"effective_n": round(self.effective_n, 3), "beta": round(self.beta, 3),
                              "base_validation_error": round(self.base_validation_error, 4) if self.base_validation_error is not None else None,
                              "adapted_validation_error": round(self.adapted_validation_error, 4) if self.adapted_validation_error is not None else None,
                              "base_engine_policy": ENGINE_V1_POLICY, "base_engine_hash": ENGINE_V1_HASH}


def _normalise(weights: dict[str, float]) -> dict[str, float]:
    total = sum(max(0., x) for x in weights.values())
    return {k: (max(0., v) / total if total else 1 / len(weights)) for k, v in weights.items()} if weights else {}

def _tier(effective_n: float) -> str:
    if effective_n < 5: return "BASE_ONLY"
    if effective_n < 10: return "LIGHT_PROJECT_ADAPTATION"
    if effective_n < 20: return "REGULARIZED_PROJECT_ENSEMBLE"
    return "LOCAL_SERIES_ADAPTATION"

def _beta(effective_n: float) -> float:
    # Strong shrinkage at activation; cannot permit an extreme five-point fit.
    if effective_n < 5: return 0.0
    return min(0.65, 0.12 + 0.53 * (effective_n - 5) / (effective_n + 15))

def _fit(events: list[QualifiedEvidencePair], global_weights: dict[str, float], query_smiles: str = "") -> tuple[dict[str, float], float]:
    model_ids = list(global_weights)
    query_fp = compute_morgan_fingerprint(query_smiles) if query_smiles else None
    errors = {m: 0.0 for m in model_ids}; masses = {m: 0.0 for m in model_ids}; eff_n = 0.0
    for event in events:
        sim = compute_tanimoto_similarity(query_fp, compute_morgan_fingerprint(event.smiles)) if query_fp else 1.0
        mass = event.quality_weight * max(0.15, sim)
        eff_n += mass
        for model in model_ids:
            if model in event.frozen_predictions:
                errors[model] += mass * abs(event.frozen_predictions[model] - event.value)
                masses[model] += mass
    scores = {m: 1 / max(0.05, errors[m] / masses[m]) if masses[m] else global_weights[m] for m in model_ids}
    return _normalise(scores), eff_n

def _loo_error(events: list[QualifiedEvidencePair], global_weights: dict[str, float], beta: float) -> tuple[float | None, float | None]:
    if len(events) < 2: return None, None
    base_errors, adapted_errors = [], []
    for index, held_out in enumerate(events):
        train = events[:index] + events[index + 1:]
        fitted, _ = _fit(train, global_weights, held_out.smiles)
        weights = _normalise({m: (1-beta)*global_weights[m] + beta*fitted[m] for m in global_weights})
        base = sum(global_weights[m] * held_out.frozen_predictions[m] for m in global_weights if m in held_out.frozen_predictions)
        adapted = sum(weights[m] * held_out.frozen_predictions[m] for m in weights if m in held_out.frozen_predictions)
        base_errors.append(abs(base - held_out.value)); adapted_errors.append(abs(adapted - held_out.value))
    return sum(base_errors)/len(base_errors), sum(adapted_errors)/len(adapted_errors)

def fit_project_adapter(endpoint_id: str, evidence: Iterable[QualifiedEvidencePair], global_weights: dict[str, float], *, query_smiles: str = "") -> ProjectAdapterResult:
    """Fit only with previously frozen, quality-qualified distinct evidence."""
    selected = [x for x in evidence if x.endpoint_id == endpoint_id and x.eligible and set(global_weights).issubset(x.frozen_predictions) and (not query_smiles or x.smiles != query_smiles)]
    fitted, effective_n = _fit(selected, _normalise(global_weights), query_smiles)
    global_weights = _normalise(global_weights)
    beta = _beta(effective_n)
    weights = _normalise({m: (1-beta)*global_weights[m] + beta*fitted[m] for m in global_weights})
    status = _tier(effective_n)
    base_error, adapted_error = _loo_error(selected, global_weights, beta)
    # Equality is not evidence of learning.  Require a strict, leakage-safe
    # validation improvement before an adapter can be activated.
    activate = status != "BASE_ONLY" and base_error is not None and adapted_error is not None and adapted_error < base_error - 1e-12
    if not activate:
        weights, status = global_weights, "BASE_ONLY" if effective_n < 5 else "NO_IMPROVEMENT_BASE_RETAINED"
    return ProjectAdapterResult(endpoint_id, status, len(selected), effective_n, global_weights, weights, beta if activate else 0.0, base_error, adapted_error, "ACTIVATED" if activate else "BASE_RETAINED")
