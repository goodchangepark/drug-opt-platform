"""Frozen prediction ↔ qualified experimental comparison layer.

This module is intentionally side-effect free.  It creates auditable pair
records from already-frozen predictions and qualified observations; it never
recalculates a prediction and never activates a project adapter.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, asdict
from datetime import datetime
from statistics import median
from typing import Any, Iterable

DIRECT = {"DIRECTLY_COMPARABLE", "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION"}


def _dt(value: Any) -> datetime | None:
    if isinstance(value, datetime):
        return value
    if not value:
        return None
    try:
        return datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    except (TypeError, ValueError):
        return None


def _num(value: Any) -> float | None:
    try:
        value = str(value).replace(",", "").strip()
        return float(value) if value else None
    except (TypeError, ValueError):
        return None


def _endpoint_id(prediction: dict, evidence: dict) -> str:
    endpoint = str(prediction.get("endpoint") or prediction.get("endpoint_id") or "")
    canonical = str(evidence.get("canonical_endpoint_id") or evidence.get("routing", {}).get("canonical_endpoint_id") or "")
    mapping = {
        "solubility_aqueous_logs": "Solubility",
        "permeability_caco2_logpapp": "Permeability",
        "ppb_human_percent_bound": "Plasma protein binding",
        "hlm_intrinsic_clearance_scaled_log10": "HLM intrinsic clearance",
        "rlm_intrinsic_clearance_scaled_log10": "RLM intrinsic clearance",
        "mlm_intrinsic_clearance_scaled_log10": "MLM intrinsic clearance",
        "pka": "pKa", "logd_7_4": "logD7.4",
    }
    if canonical in mapping and endpoint.lower() == mapping[canonical].lower():
        return canonical
    if canonical in mapping and not endpoint:
        return canonical
    aliases = {
        "solubility": "solubility_aqueous_logs", "permeability": "permeability_caco2_logpapp",
        "plasma protein binding": "ppb_human_percent_bound", "hlm intrinsic clearance": "hlm_intrinsic_clearance_scaled_log10",
        "rlm intrinsic clearance": "rlm_intrinsic_clearance_scaled_log10", "mlm intrinsic clearance": "mlm_intrinsic_clearance_scaled_log10",
        "pka": "pka", "logd7.4": "logd_7_4",
    }
    return aliases.get(endpoint.lower(), endpoint)


def _relation_status(prediction: float, relation: str, experiment: float) -> str:
    relation = (relation or "=").strip()
    if relation in {"=", "~", ""}:
        return "EQUAL_VALUE"
    if relation in {">", ">="}:
        return "BOUND_CONSISTENT" if prediction >= experiment else "BOUND_INCONSISTENT"
    if relation in {"<", "<="}:
        return "BOUND_CONSISTENT" if prediction <= experiment else "BOUND_INCONSISTENT"
    return "BOUND_INDETERMINATE"


@dataclass(frozen=True)
class PredictionExperimentalPair:
    pair_id: str
    project_id: int | str | None
    compound_id: int | str | None
    compound_version_id: int | str | None
    endpoint_id: str
    prediction_record_id: int | str | None
    experimental_evidence_id: int | str | None
    prediction_value: float | None
    prediction_unit: str
    experimental_raw_value: Any
    experimental_raw_unit: str
    experimental_normalized_value: float | None
    experimental_normalized_unit: str
    relation: str
    comparability_status: str
    comparison_metric_type: str
    absolute_error: float | None
    signed_error: float | None
    relative_error: float | None
    fold_error: float | None
    classification_match: bool | None
    adaptation_eligibility: bool
    evidence_quality: str
    independent_experiment_group_id: str
    pair_class: str
    eligibility_reason: str
    prediction_created_at: str | None
    experimental_created_at: str | None

    def to_dict(self) -> dict:
        return asdict(self)


def compare_prediction_experiment(prediction: dict, evidence: dict, *, project_id=None,
                                  compound_id=None, compound_version_id=None) -> PredictionExperimentalPair:
    """Build one pair, including a hard pre-experimental freeze gate."""
    endpoint_id = _endpoint_id(prediction, evidence)
    display = evidence.get("display") or {}
    status = str(evidence.get("comparability_status") or display.get("comparability_status") or "UNSUPPORTED")
    pred = _num(prediction.get("predicted_value", prediction.get("prediction_value")))
    exp = _num(display.get("normalized_value", evidence.get("normalized_value")))
    relation = str(evidence.get("raw_relation", evidence.get("relation", "=")) or "=")
    pred_at = _dt(prediction.get("created_at") or prediction.get("prediction_created_at"))
    exp_at = _dt(evidence.get("imported_at") or evidence.get("experimental_created_at") or evidence.get("created_at"))
    if pred_at and exp_at and pred_at < exp_at:
        pair_class, time_reason = "TRUE_PROSPECTIVE", "Frozen prediction predates experiment"
    elif pred_at and exp_at:
        pair_class, time_reason = "HISTORICAL_VISIBLE", "Experiment was visible before or at prediction time"
    else:
        pair_class, time_reason = "HISTORICAL_VISIBLE", "No complete historical freeze timestamps"

    metric = "NONE"
    absolute = signed = relative = fold = None
    class_match = None
    if status in DIRECT and pred is not None and exp is not None and endpoint_id:
        if "probability" in str(prediction.get("unit", "")).lower() or evidence.get("measurement_type") == "classification":
            metric = "CLASSIFICATION"
        elif "bound" in str(prediction.get("unit", "")).lower():
            metric = "PERCENTAGE_POINTS"
        elif "log" in (str(prediction.get("unit", "")) + str(display.get("normalized_unit", ""))).lower():
            metric = "LOG_ABSOLUTE_ERROR"
        else:
            metric = "LINEAR_ERROR"
        relation_result = _relation_status(pred, relation, exp)
        if relation_result == "EQUAL_VALUE":
            signed, absolute = pred - exp, abs(pred - exp)
            if metric == "LINEAR_ERROR" and exp != 0:
                relative, fold = abs(pred - exp) / abs(exp), max(pred / exp, exp / pred) if pred > 0 and exp > 0 else None
            elif metric == "PERCENTAGE_POINTS":
                relative, fold = None, None
        else:
            metric = relation_result
    eligible = bool(
        status in DIRECT and evidence.get("duplicate_status", "DISTINCT_MEASUREMENT") == "DISTINCT_MEASUREMENT"
        and evidence.get("import_eligible", evidence.get("adaptation_eligibility", False))
        and pair_class == "TRUE_PROSPECTIVE" and pred is not None and exp is not None
    )
    reason = "Eligible prospective directly comparable pair" if eligible else time_reason
    if status not in DIRECT:
        reason = evidence.get("gap_reason") or evidence.get("routing", {}).get("routing_reason") or "Evidence is not directly comparable"
    elif evidence.get("duplicate_status") != "DISTINCT_MEASUREMENT":
        reason = "Duplicate measurement excluded"
    elif pred is None or exp is None:
        reason = "Prediction or normalized experimental value is unavailable"
    return PredictionExperimentalPair(
        pair_id=f"PAIR-{prediction.get('id', 'p')}-{evidence.get('id', evidence.get('source_record_id', 'e'))}",
        project_id=project_id, compound_id=compound_id, compound_version_id=compound_version_id or prediction.get("version_id") or evidence.get("compound_version_id"),
        endpoint_id=endpoint_id, prediction_record_id=prediction.get("id"), experimental_evidence_id=evidence.get("id"),
        prediction_value=pred, prediction_unit=str(prediction.get("unit", "")),
        experimental_raw_value=evidence.get("raw_value", evidence.get("value")), experimental_raw_unit=str(evidence.get("raw_unit", evidence.get("unit", ""))),
        experimental_normalized_value=exp, experimental_normalized_unit=str(display.get("normalized_unit", evidence.get("normalized_unit", ""))),
        relation=relation, comparability_status=status, comparison_metric_type=metric,
        absolute_error=absolute, signed_error=signed, relative_error=relative, fold_error=fold,
        classification_match=class_match, adaptation_eligibility=eligible, evidence_quality=str(evidence.get("source_quality_class", "D")),
        independent_experiment_group_id=str(evidence.get("independent_experiment_group_id") or evidence.get("source_document_id") or evidence.get("source_record_id") or "unknown"),
        pair_class=pair_class, eligibility_reason=reason,
        prediction_created_at=pred_at.isoformat() if pred_at else None, experimental_created_at=exp_at.isoformat() if exp_at else None,
    )


def generate_pairs(predictions: Iterable[dict], evidence: Iterable[dict], **scope) -> list[PredictionExperimentalPair]:
    result = []
    for prediction in predictions:
        for item in evidence:
            if prediction.get("version_id") and item.get("compound_version_id") and prediction["version_id"] != item["compound_version_id"]:
                continue
            pair = compare_prediction_experiment(prediction, item, **scope)
            if pair.endpoint_id and pair.comparability_status not in {"UNSUPPORTED", "NOT_COMPARABLE"}:
                result.append(pair)
    return result


def independent_compound_count(pairs: Iterable[PredictionExperimentalPair], *, eligible_only=True) -> int:
    selected = [p for p in pairs if (p.adaptation_eligibility if eligible_only else True)]
    return len({p.compound_version_id for p in selected})


def performance_summary(pairs: Iterable[PredictionExperimentalPair]) -> dict:
    selected = [p for p in pairs if p.absolute_error is not None]
    errors = [p.absolute_error for p in selected]
    signed = [p.signed_error for p in selected if p.signed_error is not None]
    compounds = independent_compound_count(selected)
    return {
        "pair_count": len(selected), "independent_compounds": compounds,
        "effective_n": float(compounds),
        "mae": sum(errors) / len(errors) if errors else None,
        "rmse": math.sqrt(sum(e * e for e in errors) / len(errors)) if errors else None,
        "bias": sum(signed) / len(signed) if signed else None,
        "median_absolute_error": median(errors) if errors else None,
        "status": "Collecting" if compounds < 5 else "Eligible — validation required",
    }
