"""Leakage-safe project learning-curve validation.

This module is deliberately side-effect free.  It consumes frozen model
predictions and qualified experimental pairs, evaluates an adapter only on
compounds held out from that adapter, and returns auditable JSON-compatible
records.  It does not write a database row, activate an adapter, or alter
Engine v1.
"""
from __future__ import annotations

import math
import random
from collections import defaultdict
from statistics import median
from typing import Iterable

from .adaptive_weighting import compute_morgan_fingerprint, compute_tanimoto_similarity, get_bemis_murcko_scaffold
from .project_adaptation_v2 import ENGINE_V1_HASH, ENGINE_V1_POLICY, QualifiedEvidencePair, fit_project_adapter


LEARNING_CURVE_POLICY_VERSION = "drugopt-project-learning-curve-v1"
MATURITY_POLICY_VERSION = "drugopt-project-maturity-policy-v1"
DEFAULT_ORDERING_SEEDS = tuple(range(20))


def _mean(values):
    return sum(values) / len(values) if values else None


def _rmse(values):
    return math.sqrt(sum(value * value for value in values) / len(values)) if values else None


def _bias(values):
    return _mean(values)


def _spearman(xs, ys):
    """Small dependency-free Spearman calculation; omit it for tiny samples."""
    if len(xs) < 5 or len(xs) != len(ys):
        return None

    def ranks(values):
        order = sorted(range(len(values)), key=lambda index: values[index])
        result = [0.0] * len(values)
        position = 0
        while position < len(order):
            end = position
            while end + 1 < len(order) and values[order[end + 1]] == values[order[position]]:
                end += 1
            rank = (position + end + 2) / 2.0
            for index in order[position:end + 1]:
                result[index] = rank
            position = end + 1
        return result

    rx, ry = ranks(xs), ranks(ys)
    mx, my = _mean(rx), _mean(ry)
    numerator = sum((x - mx) * (y - my) for x, y in zip(rx, ry))
    denominator = math.sqrt(sum((x - mx) ** 2 for x in rx) * sum((y - my) ** 2 for y in ry))
    return numerator / denominator if denominator else None


def _prediction(event: QualifiedEvidencePair, weights: dict[str, float]) -> float | None:
    values = [weights[model] * event.frozen_predictions[model] for model in weights if model in event.frozen_predictions]
    return sum(values) if len(values) == len(weights) else None


def _metrics(events: list[QualifiedEvidencePair], weights: dict[str, float]) -> dict:
    errors, signed, actual, predicted = [], [], [], []
    for event in events:
        prediction = _prediction(event, weights)
        if prediction is None:
            continue
        error = prediction - event.value
        signed.append(error)
        errors.append(abs(error))
        actual.append(event.value)
        predicted.append(prediction)
    return {
        "n": len(errors),
        "mae": _mean(errors),
        "rmse": _rmse([value for value in signed]),
        "bias": _bias(signed),
        "median_absolute_error": median(errors) if errors else None,
        "spearman": _spearman(actual, predicted),
    }


def _independent_events(events: Iterable[QualifiedEvidencePair]) -> list[QualifiedEvidencePair]:
    """Select one distinct, eligible observation per compound, without averaging."""
    grouped = defaultdict(list)
    for event in events:
        if event.eligible:
            grouped[str(event.compound_version_id)].append(event)
    selected = []
    for compound_id, candidates in grouped.items():
        # Internal/high-quality observations are preferred; ties are stable.
        candidates.sort(key=lambda event: (0 if event.origin == "EXPERIMENTAL_INTERNAL" else 1,
                                           0 if event.source_quality == "A" else 1,
                                           str(event.evidence_id)))
        selected.append(candidates[0])
    return sorted(selected, key=lambda event: (str(event.compound_version_id), str(event.evidence_id)))


def _similarity_to_training(event: QualifiedEvidencePair, training: list[QualifiedEvidencePair]) -> float:
    if not training:
        return 0.0
    query = compute_morgan_fingerprint(event.smiles)
    values = [compute_tanimoto_similarity(query, compute_morgan_fingerprint(item.smiles)) for item in training]
    return max(values, default=0.0)


def _similarity_bin(value: float) -> str:
    if value >= 0.70:
        return "high"
    if value >= 0.40:
        return "medium"
    return "low"


def _candidate_for_training(endpoint_id, training, global_weights):
    if not training:
        return None
    return fit_project_adapter(endpoint_id, training, global_weights)


def _evaluate_order(events: list[QualifiedEvidencePair], endpoint_id: str, global_weights: dict[str, float], seed: int) -> list[dict]:
    ordered = list(events)
    random.Random(seed).shuffle(ordered)
    curve = []
    for n in range(0, len(ordered)):
        training = ordered[:n]
        holdout = ordered[n:]
        baseline = _metrics(holdout, global_weights)
        fit = _candidate_for_training(endpoint_id, training, global_weights) if n else None
        candidate_weights = fit.project_weights if fit and fit.activation_decision == "ACTIVATED" else None
        adapted = _metrics(holdout, candidate_weights) if candidate_weights else None
        improved = None
        if adapted and baseline["n"]:
            improved = sum(
                abs(_prediction(event, candidate_weights) - event.value) < abs(_prediction(event, global_weights) - event.value)
                for event in holdout
                if _prediction(event, candidate_weights) is not None and _prediction(event, global_weights) is not None
            ) / baseline["n"]
        decision = "BASE_ONLY"
        maturity = {"level": 1, "label": "Base Prediction", "activated": False}
        if n >= 5 and fit:
            if fit.activation_decision == "ACTIVATED":
                decision = "LIGHT_PROJECT_ADAPTATION_CANDIDATE"
                maturity = {"level": 2, "label": "Early Adaptation", "activated": False}
            else:
                decision = "BASE_RETAINED_NO_IMPROVEMENT"
        curve.append({
            "n": n,
            "seed": seed,
            "training_compounds": [str(event.compound_version_id) for event in training],
            "holdout_compounds": [str(event.compound_version_id) for event in holdout],
            "base": baseline,
            "adapted": adapted,
            "base_mae": baseline["mae"],
            "adapted_mae": adapted["mae"] if adapted else None,
            "delta_mae": (adapted["mae"] - baseline["mae"]) if adapted and adapted["mae"] is not None and baseline["mae"] is not None else None,
            "fraction_holdouts_improved": improved,
            "effective_n": round(fit.effective_n, 4) if fit else 0.0,
            "candidate_weights": dict(fit.project_weights) if fit else None,
            "global_weights": dict(global_weights),
            "validation_decision": decision,
            "maturity_candidate": maturity,
            "holdout_predictions": [
                {
                    "compound_version_id": str(event.compound_version_id),
                    "experimental_value": event.value,
                    "base_prediction": _prediction(event, global_weights),
                    "project_prediction": _prediction(event, candidate_weights) if candidate_weights else None,
                    "similarity_to_training": round(_similarity_to_training(event, training), 4),
                }
                for event in holdout
            ],
        })
    return curve


def _quantile(values, fraction):
    if not values:
        return None
    values = sorted(values)
    index = (len(values) - 1) * fraction
    lower, upper = math.floor(index), math.ceil(index)
    if lower == upper:
        return values[lower]
    return values[lower] + (values[upper] - values[lower]) * (index - lower)


def _aggregate(curves: list[list[dict]]) -> list[dict]:
    by_n = defaultdict(list)
    for curve in curves:
        for point in curve:
            by_n[point["n"]].append(point)
    result = []
    for n in sorted(by_n):
        points = by_n[n]
        base = [point["base_mae"] for point in points if point["base_mae"] is not None]
        adapted = [point["adapted_mae"] for point in points if point["adapted_mae"] is not None]
        improved = [point["fraction_holdouts_improved"] for point in points if point["fraction_holdouts_improved"] is not None]
        result.append({
            "n": n,
            "base_mae": _mean(base),
            "adapted_mae": _mean(adapted),
            "adapted_mae_median": median(adapted) if adapted else None,
            "adapted_mae_iqr": (_quantile(adapted, .25), _quantile(adapted, .75)) if adapted else None,
            "delta_mae": (_mean(adapted) - _mean(base)) if adapted and base else None,
            "fraction_holdouts_improved": _mean(improved),
            "effective_n": _mean([point["effective_n"] for point in points]),
            "maturity_candidate": next((point["maturity_candidate"] for point in points if point["maturity_candidate"]["level"] > 1), {"level": 1, "label": "Base Prediction", "activated": False}),
            "validation_decisions": sorted(set(point["validation_decision"] for point in points)),
        })
    return result


def _split_evaluation(events, endpoint_id, global_weights, mode):
    if len(events) < 2:
        return {"mode": mode, "status": "INSUFFICIENT_DATA"}
    groups = defaultdict(list)
    for event in events:
        key = event.series_id if getattr(event, "series_id", "") and mode == "series" else get_bemis_murcko_scaffold(event.smiles) if mode in {"scaffold", "series"} else str(event.compound_version_id)
        groups[key].append(event)
    if mode == "random":
        shuffled = list(events)
        random.Random(1701).shuffle(shuffled)
        cut = max(1, min(len(shuffled) - 1, round(len(shuffled) * .7)))
        training, holdout = shuffled[:cut], shuffled[cut:]
    else:
        keys = sorted(groups)
        if len(keys) < 2:
            return {"mode": mode, "status": "UNAVAILABLE_ONE_GROUP"}
        holdout = groups[keys[-1]]
        training = [event for key in keys[:-1] for event in groups[key]]
    fit = _candidate_for_training(endpoint_id, training, global_weights) if len(training) >= 5 else None
    baseline = _metrics(holdout, global_weights)
    adapted = _metrics(holdout, fit.project_weights) if fit and fit.activation_decision == "ACTIVATED" else None
    return {
        "mode": mode,
        "status": "COMPLETED",
        "training_compounds": [str(event.compound_version_id) for event in training],
        "holdout_compounds": [str(event.compound_version_id) for event in holdout],
        "base": baseline,
        "adapted": adapted,
        "decision": fit.activation_decision if fit else "INSUFFICIENT_EVIDENCE",
    }


def _similarity_analysis(events, global_weights, curve_point):
    buckets = defaultdict(lambda: {"base": [], "adapted": [], "n": 0})
    training_ids = set(curve_point["training_compounds"])
    training = [event for event in events if str(event.compound_version_id) in training_ids]
    candidate = curve_point.get("candidate_weights")
    for item in curve_point["holdout_predictions"]:
        event = next((event for event in events if str(event.compound_version_id) == item["compound_version_id"]), None)
        if not event:
            continue
        bucket = buckets[_similarity_bin(item["similarity_to_training"])]
        if item["base_prediction"] is not None:
            bucket["base"].append(abs(item["base_prediction"] - item["experimental_value"]))
        if item["project_prediction"] is not None:
            bucket["adapted"].append(abs(item["project_prediction"] - item["experimental_value"]))
        bucket["n"] += 1
    return {key: {"n": value["n"], "base_mae": _mean(value["base"]), "adapted_mae": _mean(value["adapted"])} for key, value in sorted(buckets.items())}


def build_learning_curve(endpoint_id: str, evidence: Iterable[QualifiedEvidencePair], global_weights: dict[str, float], *, ordering_seeds=DEFAULT_ORDERING_SEEDS) -> dict:
    """Build repeated, prospective holdout curves for one endpoint."""
    evidence = list(evidence)
    eligible = _independent_events(event for event in evidence if event.endpoint_id == endpoint_id)
    weights = dict(global_weights)
    curves = [_evaluate_order(eligible, endpoint_id, weights, int(seed)) for seed in ordering_seeds] if eligible else []
    primary = curves[0] if curves else []
    final_point = primary[-1] if primary else None
    return {
        "policy_version": LEARNING_CURVE_POLICY_VERSION,
        "engine_policy": ENGINE_V1_POLICY,
        "engine_hash": ENGINE_V1_HASH,
        "endpoint_id": endpoint_id,
        "raw_observations": sum(1 for event in evidence if event.endpoint_id == endpoint_id),
        "unique_experiments": len(eligible),
        "independent_compounds": len(eligible),
        "ordering_seeds": [int(seed) for seed in ordering_seeds],
        "primary_ordering": primary,
        "aggregate": _aggregate(curves),
        "split_evaluations": [_split_evaluation(eligible, endpoint_id, weights, mode) for mode in ("random", "scaffold", "series")],
        "similarity_analysis": _similarity_analysis(eligible, weights, final_point) if final_point else {},
        "validation_note": "Synthetic/public disposable validation only; no protected real project adapter was activated.",
    }


def build_synthetic_validation_dataset(endpoint_id="Solubility", count=12) -> list[QualifiedEvidencePair]:
    """Controlled acceptance fixture: Model B is intentionally best in-series."""
    smiles = [
        "c1ccccc1C", "c1ccccc1CC", "c1ccncc1C", "c1ccncc1CC", "CCOC(=O)C", "CCOC(=O)CC",
        "CCNCC", "CCNCCC", "CC(C)O", "CC(C)CO", "CC(=O)N", "CC(=O)NC",
        "c1ccccc1Cl", "c1ccccc1F", "c1ccccc1OC", "c1ccccc1OC" + "C",
    ]
    events = []
    for index in range(count):
        value = -4.0 + (index % 4) * 0.11 + index * 0.015
        events.append(QualifiedEvidencePair(
            evidence_id=f"SYN-{endpoint_id}-{index + 1}",
            compound_version_id=index + 1,
            smiles=smiles[index % len(smiles)],
            endpoint_id=endpoint_id,
            value=value,
            frozen_predictions={"model_a": value + 0.40, "model_b": value - 0.08, "model_c": value + 0.55},
            origin="EXPERIMENTAL_INTERNAL",
            source_quality="A",
            series_id="series-a" if index < count // 2 else "series-b",
        ))
    return events


def build_disposable_learning_demo(endpoint_id="Solubility") -> dict:
    events = build_synthetic_validation_dataset(endpoint_id, 6)
    global_weights = {"model_a": .60, "model_b": .25, "model_c": .15}
    curve = build_learning_curve(endpoint_id, events[:5], global_weights, ordering_seeds=(0,))
    fit = fit_project_adapter(endpoint_id, events[:5], global_weights)
    holdout = events[5]
    base = _prediction(holdout, global_weights)
    project = _prediction(holdout, fit.project_weights) if fit.activation_decision == "ACTIVATED" else None
    return {
        "training_compounds": [str(event.compound_version_id) for event in events[:5]],
        "new_compound": str(holdout.compound_version_id),
        "experiment_initially_absent": True,
        "base_prediction": base,
        "project_prediction": project,
        "maturity": {"level": 2, "label": "Early Adaptation", "activated": fit.activation_decision == "ACTIVATED"},
        "effective_n": fit.effective_n,
        "adapter_status": fit.status,
        "adapter_activation_decision": fit.activation_decision,
        "revealed_experiment": holdout.value,
        "base_error": abs(base - holdout.value),
        "project_error": abs(project - holdout.value) if project is not None else None,
        "curve": curve,
        "snapshot_immutable": True,
    }


def maturity_policy() -> dict:
    return {
        "policy_version": MATURITY_POLICY_VERSION,
        "meaning": "Project-specific experimental adaptation maturity, not an absolute guarantee of accuracy.",
        "calibration_basis": {
            "validation_artifact": "project_learning_curve_v3_6.json",
            "validation_type": "controlled synthetic disposable validation",
            "endpoints": ["PPB", "Solubility", "Caco2", "HLM"],
            "compounds_per_endpoint": 12,
            "repeat_orderings": 20,
            "observed_level_2_candidate_gate": "At N=5, candidate weights were evaluated only after leakage-safe LOO and were retained as a candidate only when holdout-compatible validation improved strictly over base.",
            "interpretation": "The controlled fixture demonstrates the gate and learning-curve machinery; it is not a claim about real-drug accuracy.",
        },
        "rules": {
            "level_1": {"label": "Base Prediction", "gate": "Valid prediction with no explicitly activated validated project adapter."},
            "level_2": {"label": "Early Adaptation", "gate": "At least five independent qualified comparable compounds, leakage-safe holdout/LOO validation, non-degradation or improvement, and explicit activation."},
            "level_3": {"label": "Project Adapted", "gate": "Stronger endpoint-specific evidence, repeated validation and stable improvement beyond the Level 2 gate; not N alone."},
            "level_4": {"label": "Series Adapted", "gate": "Validated scaffold/series-aware improvement on held-out compounds in addition to project-global validation."},
            "level_5": {"label": "Mature Project Prediction", "gate": "Substantial independent evidence, repeated successful out-of-sample cycles, low instability, representative coverage, and adapter-history consistency; never N alone."},
        },
        "hard_gates": ["same-compound leakage blocked", "pre-experimental frozen prediction required", "raw/duplicate/related evidence cannot increase maturity", "weights non-negative and sum to one", "explicit user activation required"],
    }
