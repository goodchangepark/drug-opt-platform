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
    "Solubility": ("SOLUBILITY_GENERIC", "ROUTED_IMPLEMENTATION_MISMATCH"),
    "Permeability": ("CACO2_PAPP_AB", "ROUTED_IMPLEMENTATION_MISMATCH"),
    "Plasma protein binding": ("HUMAN_PPB", "ROUTED_IMPLEMENTATION_MISMATCH"),
    "HLM intrinsic clearance": ("HLM_CLINT", "ROUTED_IMPLEMENTATION_MISMATCH"),
    "Ames mutagenicity": ("AMES_MUTAGENICITY", "INVALID_SPECIES"),
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
