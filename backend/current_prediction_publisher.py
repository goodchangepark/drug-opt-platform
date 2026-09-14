"""Publish verified runtime predictions into the Stable Core authority.

Legacy calculation tables remain immutable artifacts.  This module is the
single producer adapter that admits only outputs whose exact installed model,
scientific identity, unit, and structure revision can be proven.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Iterable

from sqlalchemy import select

from .model_artifact_authority import artifact_bundle_sha256, model_artifact_registration
from .prediction_engine_registry import CURRENT_ENGINE_ID
from .request_context import current_request_trace
from .stable_core import (
    CurrentPredictionSnapshot,
    ScientificMutationAudit,
    admit_current_prediction,
    expected_structure_revision,
    prediction_context_identity,
    publish_current_prediction,
)


# Explicitly distinguish binary screening outputs from quantitative inhibition.
# The general alias normalizer intentionally cannot make this decision because
# the same legacy label has historically been used for both tasks.
ADMET_CURRENT_ENDPOINTS: dict[str, tuple[str, str, dict[str, Any]]] = {
    "RLM intrinsic clearance": ("RLM_CLINT", "RAT", {"matrix": "microsomes"}),
    "MLM intrinsic clearance": ("MLM_CLINT", "MOUSE", {"matrix": "microsomes"}),
    "CYP1A2 inhibitor": ("CYP1A2_INHIBITOR_CLASS", "HUMAN", {"assay_type": "binary inhibitor classification"}),
    "CYP2C9 inhibitor": ("CYP2C9_INHIBITOR_CLASS", "HUMAN", {"assay_type": "binary inhibitor classification"}),
    "CYP2C19 inhibitor": ("CYP2C19_INHIBITOR_CLASS", "HUMAN", {"assay_type": "binary inhibitor classification"}),
    "CYP2D6 inhibitor": ("CYP2D6_INHIBITOR_CLASS", "HUMAN", {"assay_type": "binary inhibitor classification"}),
    "CYP3A4 inhibitor": ("CYP3A4_INHIBITOR_CLASS", "HUMAN", {"assay_type": "binary inhibitor classification"}),
    "CYP2C9 substrate": ("CYP2C9_SUBSTRATE", "HUMAN", {"assay_type": "binary substrate classification"}),
    "CYP2D6 substrate": ("CYP2D6_SUBSTRATE", "HUMAN", {"assay_type": "binary substrate classification"}),
    "CYP3A4 substrate": ("CYP3A4_SUBSTRATE", "HUMAN", {"assay_type": "binary substrate classification"}),
    "P-gp inhibitor": ("PGP_INHIBITION", "HUMAN", {"assay_type": "binary inhibitor classification"}),
    "hERG liability": ("HERG_CLASS", "HUMAN", {"assay_type": "binary blocker classification"}),
    "DILI clinical liability": ("DILI_LIABILITY", "HUMAN", {"assay_type": "binary clinical liability classification"}),
}

ADMET_UNQUALIFIED_ENDPOINTS: dict[str, tuple[str, str]] = {
    # The legacy ADMET panel still calculates these base-model outputs.  They
    # are reported as ineligible here; the full Predict workflow separately
    # executes and publishes the exact current-production routes.
    "Solubility": ("SOLUBILITY_GENERIC", "ROUTED_IMPLEMENTATION_MISMATCH"),
    "Permeability": ("CACO2_PAPP_AB", "ROUTED_IMPLEMENTATION_MISMATCH"),
    "Plasma protein binding": ("HUMAN_PPB", "ROUTED_IMPLEMENTATION_MISMATCH"),
    "HLM intrinsic clearance": ("HLM_CLINT", "ROUTED_IMPLEMENTATION_MISMATCH"),
}


def _audit(db, *, snapshot: CurrentPredictionSnapshot, before: dict, after: dict, operation: str) -> None:
    from .models import Compound, CompoundVersion

    trace = current_request_trace()
    project_id = db.scalar(
        select(Compound.project_id)
        .join(CompoundVersion, CompoundVersion.compound_row_id == Compound.id)
        .where(CompoundVersion.id == snapshot.compound_version_id)
    )
    db.add(ScientificMutationAudit(
        actor=trace.caller_name,
        process=trace.route_action,
        operation=operation,
        object_type="CurrentPredictionSnapshot",
        object_id=str(snapshot.id or ""),
        project_id_snapshot=project_id,
        transaction_id=trace.transaction_id,
        caller_type=trace.caller_type,
        caller_name=trace.caller_name,
        request_id=trace.request_id,
        workflow_id=trace.workflow_id,
        execution_id=trace.execution_id,
        conversation_id=trace.conversation_id,
        route_action=trace.route_action,
        before_identity_json=before,
        after_identity_json=after,
        reason="Verified prediction workflow publication into Stable Core authority",
        outcome="COMMITTED",
    ))


def _identity(snapshot: CurrentPredictionSnapshot) -> dict:
    return {
        "id": snapshot.id,
        "compound_version_id": snapshot.compound_version_id,
        "canonical_endpoint": snapshot.canonical_endpoint,
        "species": snapshot.species,
        "context_identity": snapshot.context_identity,
        "engine_release": snapshot.engine_release,
        "value": snapshot.value,
        "unit": snapshot.unit,
        "classification": snapshot.classification,
        "model_id": snapshot.model_id,
        "model_version": snapshot.model_version,
        "model_artifact_hash": snapshot.model_artifact_hash,
        "prediction_mode": snapshot.prediction_mode,
        "source_artifact_type": snapshot.source_artifact_type,
        "source_artifact_id": snapshot.source_artifact_id,
        "is_current": snapshot.is_current,
    }


def _upsert_verified_snapshot(db, candidate: CurrentPredictionSnapshot) -> CurrentPredictionSnapshot:
    existing = db.scalar(select(CurrentPredictionSnapshot).where(
        CurrentPredictionSnapshot.compound_version_id == candidate.compound_version_id,
        CurrentPredictionSnapshot.canonical_endpoint == candidate.canonical_endpoint,
        CurrentPredictionSnapshot.species == candidate.species,
        CurrentPredictionSnapshot.context_identity == candidate.context_identity,
        CurrentPredictionSnapshot.engine_release == candidate.engine_release,
    ))
    if existing is None:
        publish_current_prediction(db, candidate)
        _audit(db, snapshot=candidate, before={}, after=_identity(candidate), operation="PUBLISH")
        return candidate

    fields = (
        "context_json", "value", "unit", "classification", "model_id", "model_version",
        "model_artifact_hash", "applicability_domain_json", "uncertainty_json",
        "prediction_mode", "source_artifact_type", "source_artifact_id", "structure_revision",
        "created_at", "replaced_at", "is_current",
    )
    before_values = {field: getattr(existing, field) for field in fields}
    before = _identity(existing)
    for field in fields[:-3]:
        setattr(existing, field, getattr(candidate, field))
    existing.is_current = True
    existing.created_at = datetime.now(timezone.utc)
    existing.replaced_at = None
    admission = admit_current_prediction(db, existing)
    if not admission.eligible:
        for field, value in before_values.items():
            setattr(existing, field, value)
        db.flush()
        raise ValueError(f"CURRENT_PREDICTION_REJECTED:{admission.reason}")
    db.flush()
    _audit(db, snapshot=existing, before=before, after=_identity(existing), operation="REPLACE")
    return existing


def publish_admet_current_predictions(db, version, predictions: Iterable[Any]) -> list[dict[str, Any]]:
    """Publish the newest exact CORE result per supported scientific endpoint."""
    predictions = list(predictions)
    candidates: dict[str, list[Any]] = {}
    for prediction in predictions:
        model = getattr(prediction, "model", None)
        endpoint_name = str(getattr(model, "endpoint_name", "") or "")
        if endpoint_name in ADMET_CURRENT_ENDPOINTS:
            candidates.setdefault(endpoint_name, []).append(prediction)

    results: list[dict[str, Any]] = []
    seen_names = {str(getattr(getattr(row, "model", None), "endpoint_name", "") or "") for row in predictions}
    for endpoint_name, (endpoint_id, reason) in sorted(ADMET_UNQUALIFIED_ENDPOINTS.items()):
        if endpoint_name in seen_names:
            results.append({"endpoint": endpoint_id, "status": "CALCULATED_BUT_NOT_ELIGIBLE", "reason": reason})
    for endpoint_name, endpoint_predictions in sorted(candidates.items()):
        endpoint_id, species, context = ADMET_CURRENT_ENDPOINTS[endpoint_name]
        registration = model_artifact_registration(endpoint_id)
        if registration is None:
            results.append({"endpoint": endpoint_id, "status": "CALCULATED_BUT_NOT_ELIGIBLE", "reason": "MODEL_NOT_REGISTERED"})
            continue
        matching = [
            row for row in endpoint_predictions
            if row.model.model_name == registration.model_id and row.model.model_version == registration.model_version
        ]
        if not matching:
            results.append({"endpoint": endpoint_id, "status": "CALCULATED_BUT_NOT_ELIGIBLE", "reason": "REGISTRY_MISMATCH"})
            continue
        prediction = max(matching, key=lambda row: row.id)
        model = prediction.model
        if model.model_name != registration.model_id or model.model_version != registration.model_version:
            results.append({"endpoint": endpoint_id, "status": "CALCULATED_BUT_NOT_ELIGIBLE", "reason": "REGISTRY_MISMATCH"})
            continue
        artifact_hash = artifact_bundle_sha256(registration)
        if not artifact_hash:
            results.append({"endpoint": endpoint_id, "status": "CALCULATED_BUT_NOT_ELIGIBLE", "reason": "MISSING_MODEL_ARTIFACT"})
            continue
        outputs = dict(prediction.outputs_json or {})
        classification = str(outputs.get("classification") or "")
        candidate = CurrentPredictionSnapshot(
            compound_version_id=version.id,
            canonical_endpoint=endpoint_id,
            species=species,
            context_json=context,
            context_identity=prediction_context_identity(endpoint_id, context),
            engine_release=CURRENT_ENGINE_ID,
            value=prediction.predicted_value,
            unit=prediction.unit,
            classification=classification,
            model_id=registration.model_id,
            model_version=registration.model_version,
            model_artifact_hash=artifact_hash,
            applicability_domain_json={
                "classification": prediction.applicability_domain,
                "confidence": prediction.confidence,
            },
            uncertainty_json={"value": prediction.uncertainty, "source": "runtime model output"},
            prediction_mode="FULL_PREDICTION",
            source_artifact_type="ADMETPrediction",
            source_artifact_id=prediction.id,
            structure_revision=expected_structure_revision(version),
            is_current=True,
        )
        try:
            published = _upsert_verified_snapshot(db, candidate)
        except ValueError as exc:
            reason = str(exc).split(":", 1)[-1]
            results.append({"endpoint": endpoint_id, "status": "CALCULATED_BUT_NOT_ELIGIBLE", "reason": reason})
            continue
        results.append({
            "endpoint": endpoint_id,
            "status": "CALCULATED_AND_PUBLISHED",
            "snapshot_id": published.id,
            "source_artifact_id": prediction.id,
        })
    return results


def publish_cached_admet_current_predictions(db, version) -> list[dict[str, Any]]:
    from .admet import ADMETPrediction

    predictions = list(db.scalars(select(ADMETPrediction).where(
        ADMETPrediction.version_id == version.id,
        ADMETPrediction.execution_status == "SUCCESS",
    )))
    return publish_admet_current_predictions(db, version, predictions)


PROPERTY_CURRENT_ENDPOINTS: dict[str, tuple[str, str]] = {
    "MW": ("exact_molecular_weight", "g/mol"),
    "CLOGP": ("clogp", "logP"),
    "TPSA": ("tpsa", "Å²"),
    "HBD": ("hbd", "count"),
    "HBA": ("hba", "count"),
    "ROTB": ("rotatable_bonds", "count"),
    "FSP3": ("fraction_csp3", "fraction"),
    "QED": ("qed", "score (0-1)"),
    "FORMAL_CHARGE": ("formal_charge", "charge"),
    "HEAVY_ATOM_COUNT": ("heavy_atom_count", "count"),
}


def _registered_candidate(
    version,
    endpoint_id: str,
    *,
    value: float | None,
    unit: str,
    species: str = "HUMAN",
    context: dict[str, Any] | None = None,
    prediction_mode: str = "FULL_PREDICTION",
    ad: dict[str, Any] | None = None,
    uncertainty: dict[str, Any] | None = None,
    source_artifact_type: str,
    source_artifact_id: int | None = None,
) -> CurrentPredictionSnapshot | None:
    registration = model_artifact_registration(endpoint_id)
    if registration is None or value is None:
        return None
    artifact_hash = artifact_bundle_sha256(registration)
    if not artifact_hash:
        return None
    context = context or {}
    return CurrentPredictionSnapshot(
        compound_version_id=version.id,
        canonical_endpoint=endpoint_id,
        species=species,
        context_json=context,
        context_identity=prediction_context_identity(endpoint_id, context),
        engine_release=CURRENT_ENGINE_ID,
        value=float(value),
        unit=unit,
        classification="",
        model_id=registration.model_id,
        model_version=registration.model_version,
        model_artifact_hash=artifact_hash,
        applicability_domain_json=ad or {"classification": "IN_DOMAIN"},
        uncertainty_json=uncertainty or {},
        prediction_mode=prediction_mode,
        source_artifact_type=source_artifact_type,
        source_artifact_id=source_artifact_id,
        structure_revision=expected_structure_revision(version),
        is_current=True,
    )


def publish_property_current_predictions(db, version) -> list[dict[str, Any]]:
    """Publish deterministic Stage-1 outputs through unchanged admission."""
    properties = dict(version.properties_json or {})
    values: list[tuple[str, Any, str, dict[str, Any]]] = [
        (endpoint, properties.get(key), unit, {})
        for endpoint, (key, unit) in PROPERTY_CURRENT_ENDPOINTS.items()
    ]
    results = []
    for endpoint_id, value, unit, context in values:
        candidate = _registered_candidate(
            version, endpoint_id, value=value, unit=unit, context=context,
            source_artifact_type="PropertyCalculation",
        )
        if candidate is None:
            results.append({"endpoint": endpoint_id, "status": "NOT_EMITTED", "reason": "VALUE_OR_REGISTERED_ARTIFACT_UNAVAILABLE"})
            continue
        try:
            published = _upsert_verified_snapshot(db, candidate)
            results.append({"endpoint": endpoint_id, "status": "CALCULATED_AND_PUBLISHED", "snapshot_id": published.id})
        except ValueError as exc:
            results.append({"endpoint": endpoint_id, "status": "CALCULATED_BUT_NOT_ELIGIBLE", "reason": str(exc).split(":", 1)[-1]})
    return results


GLOBAL_CURRENT_CONTEXT: dict[str, tuple[str, dict[str, Any], str]] = {
    "SOLUBILITY_GENERIC": ("HUMAN", {}, "log10(mol/L)"),
    "CACO2_PAPP_AB": ("HUMAN", {"matrix": "Caco-2", "direction": "A→B"}, "log10(cm/s)"),
    "HUMAN_PPB": ("HUMAN", {"matrix": "plasma"}, "% bound"),
    "HLM_CLINT": ("HUMAN", {"matrix": "microsomes"}, "log10(mL/min/kg)"),
    "RLM_CLINT": ("RAT", {"matrix": "microsomes"}, "log10(mL/min/kg)"),
    "MLM_CLINT": ("MOUSE", {"matrix": "microsomes"}, "log10(mL/min/kg)"),
    "CYP1A2_INHIBITION": ("HUMAN", {"assay_type": "quantitative inhibition"}, "pIC50"),
    "CYP2C9_INHIBITION": ("HUMAN", {"assay_type": "quantitative inhibition"}, "pIC50"),
    "CYP2D6_INHIBITION": ("HUMAN", {"assay_type": "quantitative inhibition"}, "pIC50"),
    "CYP3A4_INHIBITION": ("HUMAN", {"assay_type": "quantitative inhibition"}, "pIC50"),
    "HERG_LIABILITY": ("HUMAN", {"assay_type": "quantitative inhibition"}, "pIC50"),
}

GLOBAL_CURRENT_ENDPOINT_ALIASES = {
    "CACO2_PERMEABILITY": "CACO2_PAPP_AB",
    "HLM_INTRINSIC_CLEARANCE": "HLM_CLINT",
}


def publish_global_current_predictions(db, version, predictions: dict[str, dict[str, Any]], source_artifact_id: int | None = None) -> list[dict[str, Any]]:
    """Publish exact routed v3.3.3 values, never legacy-model substitutes."""
    from .current_production_executor import EXECUTION_CONTRACT

    results = []
    for raw_endpoint_id, output in predictions.items():
        endpoint_id = GLOBAL_CURRENT_ENDPOINT_ALIASES.get(raw_endpoint_id, raw_endpoint_id)
        if endpoint_id not in GLOBAL_CURRENT_CONTEXT:
            continue
        registration = model_artifact_registration(endpoint_id)
        if (
            output.get("execution_contract") != EXECUTION_CONTRACT
            or output.get("execution_status") != "SUCCESS"
            or registration is None
            or output.get("model_id") != registration.model_id
            or output.get("model_version") != registration.model_version
            or output.get("route") != registration.model_route
        ):
            results.append({
                "endpoint": endpoint_id,
                "status": "CALCULATED_BUT_NOT_ELIGIBLE",
                "reason": "ROUTED_IMPLEMENTATION_MISMATCH",
            })
            continue
        species, context, unit = GLOBAL_CURRENT_CONTEXT[endpoint_id]
        value = output.get("production_prediction")
        candidate = _registered_candidate(
            version, endpoint_id, value=value, unit=unit, species=species, context=context,
            ad={
                "classification": output.get("applicability_domain", "UNKNOWN"),
                "nearest_neighbor_similarity": output.get("nearest_neighbor_similarity"),
                "guard_applied": output.get("ad_extrapolation_guard_applied", False),
            },
            uncertainty={"value": output.get("prediction_uncertainty"), "source": "current engine routed output"},
            source_artifact_type="ADMETPredictionRun",
            source_artifact_id=source_artifact_id,
        )
        if candidate is None:
            results.append({"endpoint": endpoint_id, "status": "NOT_EMITTED", "reason": "VALUE_OR_REGISTERED_ARTIFACT_UNAVAILABLE"})
            continue
        try:
            published = _upsert_verified_snapshot(db, candidate)
            results.append({"endpoint": endpoint_id, "status": "CALCULATED_AND_PUBLISHED", "snapshot_id": published.id})
        except ValueError as exc:
            results.append({"endpoint": endpoint_id, "status": "CALCULATED_BUT_NOT_ELIGIBLE", "reason": str(exc).split(":", 1)[-1]})
    return results
