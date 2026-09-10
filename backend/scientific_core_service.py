"""Single Stable Core scientific row and comparison authority."""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any

from sqlalchemy import select

from .canonical_endpoints import REGISTRY
from .prediction_engine_registry import CURRENT_ENGINE_ID, get_current_production_routing
from .prediction_maturity import ENDPOINT_MATURITY_REGISTRY, get_endpoint_maturity
from .stable_core import (
    CurrentPredictionSnapshot,
    ExperimentalObservation,
    UNKNOWN_PROVENANCE,
    admit_current_prediction,
    canonical_prediction_context,
)
from .unit_normalization import (
    clean_unit_str,
    convert_auc,
    convert_clearance,
    convert_concentration,
    convert_ppb,
    convert_time,
    convert_volume,
)

CONTEXTUAL_EXPOSURE_TOKENS = ("CMAX", "AUC", "TMAX", "_F_", "_KA_")
MISSING_CONTEXT_VALUES = {"", "UNKNOWN", "UNSPECIFIED", "NOT_REPORTED", "NONE", "N/A"}


def _relevant_context(endpoint: str, context: dict | None) -> dict:
    return canonical_prediction_context(endpoint, context)


def _compatible_context(endpoint: str, experimental: dict, prediction: dict) -> tuple[bool, str]:
    left = _relevant_context(endpoint, experimental)
    right = _relevant_context(endpoint, prediction)
    if any(token in endpoint.upper() for token in CONTEXTUAL_EXPOSURE_TOKENS):
        required = ("route", "dose", "dose_unit")
        missing = [field for field in required if str(left.get(field, "")).strip().upper() in MISSING_CONTEXT_VALUES or str(right.get(field, "")).strip().upper() in MISSING_CONTEXT_VALUES]
        if missing:
            return False, "CONTEXT_REQUIRED:" + ",".join(missing)
    for field in sorted(set(left) & set(right)):
        if str(left[field]).strip().upper() != str(right[field]).strip().upper():
            return False, f"CONTEXT_MISMATCH:{field}"
    one_sided = sorted(field for field in set(left) ^ set(right) if field in {"route", "dose", "dose_unit", "formulation", "regimen", "fed_state"})
    if one_sided:
        return False, "CONTEXT_MISSING_ON_ONE_SIDE:" + ",".join(one_sided)
    return True, "COMPATIBLE_CONTEXT"


def _compatible_numeric_values(endpoint: str, exp: float, exp_unit: str, pred: float, pred_unit: str):
    """Return deterministic same-semantics values, never assumption-based conversions."""
    left, right = clean_unit_str(exp_unit), clean_unit_str(pred_unit)
    if left == right:
        return exp, pred, exp_unit or pred_unit
    ep = endpoint.upper()
    try:
        if "T_HALF" in ep or "HALF_LIFE" in ep or "TMAX" in ep:
            return convert_time(exp, exp_unit, "hours").normalized_value, convert_time(pred, pred_unit, "hours").normalized_value, "hours"
        if "CMAX" in ep:
            return convert_concentration(exp, exp_unit, "ng/mL").normalized_value, convert_concentration(pred, pred_unit, "ng/mL").normalized_value, "ng/mL"
        if "AUC" in ep:
            return convert_auc(exp, exp_unit, "ng*h/mL").normalized_value, convert_auc(pred, pred_unit, "ng*h/mL").normalized_value, "ng*h/mL"
        if "PPB" in ep:
            return convert_ppb(exp, exp_unit, "% bound").normalized_value, convert_ppb(pred, pred_unit, "% bound").normalized_value, "% bound"
        if "_F_" in ep:
            def percent(value, unit):
                cleaned = clean_unit_str(unit)
                if "%" in cleaned or "percent" in cleaned:
                    return value
                if cleaned in {"fraction", "ratio", ""} and 0 <= value <= 1:
                    return value * 100
                raise ValueError("ambiguous bioavailability unit")
            return percent(exp, exp_unit), percent(pred, pred_unit), "%"
        if any(token in ep for token in ("CLINT", "_CL_", "CLEARANCE")) and "/kg" in left and "/kg" in right:
            return convert_clearance(exp, exp_unit, "mL/min/kg").normalized_value, convert_clearance(pred, pred_unit, "mL/min/kg").normalized_value, "mL/min/kg"
        if any(token in ep for token in ("VDSS", "_VD_", "_VDF_", "_VSSF_")) and "/kg" in left and "/kg" in right:
            return convert_volume(exp, exp_unit, "L/kg").normalized_value, convert_volume(pred, pred_unit, "L/kg").normalized_value, "L/kg"
    except (TypeError, ValueError):
        return None
    return None


def compare_scientific_pair(experimental: ExperimentalObservation, prediction: CurrentPredictionSnapshot) -> dict:
    display = experimental.canonical_endpoint == prediction.canonical_endpoint
    if not display:
        return {"display_comparable": False, "numeric_pairable": False, "learning_eligible": False, "fold_error": None, "signed_log_error": None, "reason": "ENDPOINT_SEMANTICS_MISMATCH"}
    if experimental.species != prediction.species:
        return {"display_comparable": True, "numeric_pairable": False, "learning_eligible": False, "fold_error": None, "signed_log_error": None, "reason": "SPECIES_MISMATCH"}
    context_ok, context_reason = _compatible_context(experimental.canonical_endpoint, experimental.context_json, prediction.context_json)
    if not context_ok:
        return {"display_comparable": True, "numeric_pairable": False, "learning_eligible": False, "fold_error": None, "signed_log_error": None, "reason": context_reason}
    endpoint = REGISTRY.get(experimental.canonical_endpoint)
    if endpoint and endpoint.canonical_scale in {"PROBABILITY", "CLASSIFICATION", "RANKING", "MIXED"}:
        return {"display_comparable": True, "numeric_pairable": False, "learning_eligible": False, "fold_error": None, "signed_log_error": None, "reason": "CLASSIFICATION_OR_NON_CONTINUOUS_ENDPOINT"}
    if experimental.numeric_value is None or prediction.value is None:
        return {"display_comparable": True, "numeric_pairable": False, "learning_eligible": False, "fold_error": None, "signed_log_error": None, "reason": "NON_NUMERIC_VALUE"}
    normalized = _compatible_numeric_values(
        experimental.canonical_endpoint,
        experimental.numeric_value,
        experimental.unit,
        prediction.value,
        prediction.unit,
    )
    if normalized is None:
        return {"display_comparable": True, "numeric_pairable": False, "learning_eligible": False, "fold_error": None, "signed_log_error": None, "reason": "UNIT_MISMATCH"}
    if experimental.qualifier not in {"=", "~"}:
        return {"display_comparable": True, "numeric_pairable": False, "learning_eligible": False, "fold_error": None, "signed_log_error": None, "reason": "CENSORED_VALUE_REQUIRES_INTERVAL_COMPARISON"}
    exp, pred, comparison_unit = normalized
    if exp <= 0 or pred <= 0:
        return {"display_comparable": True, "numeric_pairable": False, "learning_eligible": False, "fold_error": None, "signed_log_error": None, "reason": "FOLD_ERROR_REQUIRES_POSITIVE_VALUES"}
    ratio = pred / exp
    return {
        "display_comparable": True,
        "numeric_pairable": True,
        "learning_eligible": bool(experimental.learning_eligible),
        "prediction_experimental_ratio": ratio,
        "fold_error": max(ratio, 1.0 / ratio),
        "signed_log_error": math.log10(ratio),
        "comparison_unit": comparison_unit,
        "reason": "NUMERIC_PAIRABLE",
    }


def model_authority(db, endpoint: str, prediction: CurrentPredictionSnapshot | None = None) -> dict:
    route = next((row for row in get_current_production_routing() if row.get("endpoint_id") == endpoint), None)
    unknown_maturity = {
        "level": 0,
        "label": "Unknown maturity",
        "stars": "",
        "reason": "No admitted current model artifact is selected for this scientific row.",
        "model_route": "MODEL_UNAVAILABLE",
        "engine_version": UNKNOWN_PROVENANCE,
        "validation_n": 0,
        "locked_test_n": 0,
        "real_world_n": 0,
        "is_mechanistic": False,
        "is_unavailable": True,
        "ad_status": "UNKNOWN",
    }
    maturity = unknown_maturity
    if prediction is not None:
        admission = admit_current_prediction(db, prediction)
        if admission.eligible and endpoint in ENDPOINT_MATURITY_REGISTRY:
            maturity = get_endpoint_maturity(endpoint)
            maturity = {
                **maturity,
                "artifact_hash": prediction.model_artifact_hash,
                "model_id": prediction.model_id,
                "model_version": prediction.model_version,
                "validation_artifact": admission.model_registration.validation_artifact,
            }
    return {
        "route": route or {"endpoint_id": endpoint, "status": "MODEL_UNAVAILABLE"},
        "maturity": maturity,
    }


def build_scientific_endpoint_rows(db, version_id: int, *, category: str | None = None) -> dict:
    observations = list(db.scalars(select(ExperimentalObservation).where(
        ExperimentalObservation.compound_version_id == version_id,
        ExperimentalObservation.curation_status == "ACCEPTED",
    ).order_by(ExperimentalObservation.canonical_endpoint, ExperimentalObservation.id)))
    predictions = list(db.scalars(select(CurrentPredictionSnapshot).where(
        CurrentPredictionSnapshot.compound_version_id == version_id,
        CurrentPredictionSnapshot.is_current.is_(True),
        CurrentPredictionSnapshot.engine_release == CURRENT_ENGINE_ID,
    ).order_by(CurrentPredictionSnapshot.created_at.desc(), CurrentPredictionSnapshot.id.desc())))
    grouped_observations: dict[tuple, list] = defaultdict(list)
    grouped_predictions: dict[tuple, list] = defaultdict(list)
    for row in observations:
        relevant = _relevant_context(row.canonical_endpoint, row.context_json)
        grouped_observations[(row.canonical_endpoint, row.species, tuple(sorted(relevant.items())))].append(row)
    for row in predictions:
        # Revalidate persisted flags and provenance at read time.  A legacy
        # publisher or migration cannot make malformed data current merely by
        # leaving is_current=true.
        if not admit_current_prediction(db, row).eligible:
            continue
        relevant = _relevant_context(row.canonical_endpoint, row.context_json)
        grouped_predictions[(row.canonical_endpoint, row.species, tuple(sorted(relevant.items())))].append(row)
    keys = sorted(
        set(grouped_observations) | set(grouped_predictions),
        key=lambda item: (item[0], item[1], repr(item[2])),
    )
    rows = []
    for endpoint_id, species, context_key in keys:
        definition = REGISTRY.get(endpoint_id)
        section = definition.section if definition else "UNCLASSIFIED"
        if category and section.upper() != category.upper():
            continue
        exps = grouped_observations.get((endpoint_id, species, context_key), [])
        preds = grouped_predictions.get((endpoint_id, species, context_key), [])
        # Current means the explicitly routed current engine.  A legacy
        # artifact with missing/older provenance remains inspectable but must
        # never be silently promoted into the Current Prediction column.
        prediction = preds[0] if preds else None
        experimental = exps[0] if exps else None
        comparison = (
            compare_scientific_pair(experimental, prediction)
            if experimental and prediction else {
                "display_comparable": bool(experimental or prediction), "numeric_pairable": False,
                "learning_eligible": False, "fold_error": None, "signed_log_error": None,
                "reason": "NO_CURRENT_PREDICTION" if experimental else "NO_ACCEPTED_EXPERIMENT",
            }
        )
        authority = model_authority(db, endpoint_id, prediction)
        rows.append({
            "canonical_endpoint": endpoint_id,
            "display_name": definition.display_name if definition else endpoint_id,
            "category": section,
            "species": species,
            "context_identity": dict(context_key),
            "context": (experimental.context_json if experimental else prediction.context_json if prediction else {}),
            "experimental": None if not experimental else {
                "observation_id": experimental.id, "value": experimental.numeric_value,
                "display_value": experimental.value_text, "unit": experimental.unit,
                "qualifier": experimental.qualifier, "source": experimental.source,
                "source_url": experimental.source_url, "status": experimental.curation_status,
            },
            "prediction": None if not prediction else {
                "snapshot_id": prediction.id, "value": prediction.value, "unit": prediction.unit,
                "classification": prediction.classification, "model_id": prediction.model_id,
                "model_version": prediction.model_version,
                "model_artifact_hash": prediction.model_artifact_hash,
                "engine_version": prediction.engine_release or UNKNOWN_PROVENANCE,
                "applicability_domain": prediction.applicability_domain_json,
                "uncertainty": prediction.uncertainty_json, "mode": prediction.prediction_mode,
            },
            "comparison": comparison,
            "maturity": authority["maturity"],
            "model_status": authority["route"].get("status", "MODEL_UNAVAILABLE"),
        })
    return {
        "contract": "ScientificEndpointRow/stable-core-v1",
        "compound_version_id": version_id,
        "category": category,
        "current_engine": CURRENT_ENGINE_ID,
        "rows": rows,
    }
