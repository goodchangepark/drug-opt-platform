"""Canonical internal experimental-evidence ingestion.

Manual capture deliberately writes the same durable scientific-evidence row
used for searched/imported observations.  It never updates a prediction.
"""
from __future__ import annotations

import hashlib
import json
import uuid
from datetime import datetime, timezone

from fastapi import HTTPException
from sqlalchemy import select
from sqlalchemy.orm import Session

from .canonical_endpoints import (CANONICAL_ENDPOINT_VERSION,
    COMPARISON_UNIT_VERSION, REGISTRY, normalize_experimental_observation)
from .models import Compound, CompoundVersion, ExternalExperimentalEvidence
from .activity_models import AssayDefinition


def _current_version(compound: Compound) -> CompoundVersion:
    version = next((item for item in compound.versions if item.version_number == compound.current_version), None)
    if not version:
        raise HTTPException(status_code=400, detail="Compound has no current version")
    return version


def entry_options() -> dict:
    """Shared form metadata; frontend has no independent endpoint taxonomy."""
    endpoints = []
    for endpoint in REGISTRY.values():
        required = []
        if endpoint.species_requirement: required.append("species")
        if endpoint.matrix_requirement: required.append("matrix")
        if endpoint.direction_requirement: required.append("direction")
        if endpoint.route_requirement: required.append("route")
        endpoints.append({"canonical_endpoint_id": endpoint.canonical_endpoint_id,
                          "section": endpoint.section, "display_name": endpoint.display_name,
                          "canonical_unit": endpoint.canonical_unit,
                          "value_type": endpoint.value_type,
                          "required_fields": required})
    endpoints.extend([
        {"canonical_endpoint_id": "", "raw_endpoint": kind, "section": "ACTIVITY", "display_name": kind,
         "canonical_unit": "nM", "value_type": "numeric", "required_fields": ["assay"]}
        for kind in ("IC50", "EC50", "Ki", "Kd")
    ])
    # PK endpoint identity is species/route-aware and therefore expanded from
    # the same canonical naming contract rather than duplicated in JS.
    pk_parameters = ("CL", "CL/F", "Vd", "Vss", "Vd/F", "F", "Cmax", "Tmax", "AUC0-t", "AUC0-inf", "AUCtau", "t1/2")
    return {"canonical_endpoint_version": CANONICAL_ENDPOINT_VERSION,
            "comparison_unit_version": COMPARISON_UNIT_VERSION,
            "endpoints": endpoints,
            "pk_parameters": pk_parameters,
            "units": {"PPB": ["% bound", "fraction bound", "fu", "% unbound"],
                      "CACO2": ["×10^-6 cm/s", "cm/s", "log10(cm/s)"],
                      "PK_CMAX": ["ng/mL", "µg/L", "µg/mL", "mg/L"],
                      "PK_AUC": ["ng*h/mL", "µg*h/L", "mg*h/L"],
                      "PK_TIME": ["min", "h", "day"],
                      "PK_CL": ["mL/min/kg", "L/h/kg", "L/h"],
                      "PK_VOLUME": ["L/kg", "mL/kg"], "PK_F": ["%", "fraction"]}}


def _raw_endpoint(payload: dict) -> str:
    endpoint = str(payload.get("canonical_endpoint_id") or "").upper()
    # Do not smuggle a required context field into raw semantics merely because
    # the selected canonical label contains it.  Caco-2 A→B remains reviewable
    # until the user explicitly records direction.
    if endpoint in {"CACO2_PAPP_AB", "CACO2_PAPP_BA"}:
        return "Caco-2 Papp"
    if endpoint in REGISTRY:
        return REGISTRY[endpoint].display_name
    parameter = str(payload.get("parameter") or "").strip()
    if parameter:
        return parameter
    return str(payload.get("raw_endpoint") or "").strip()


def _context(payload: dict) -> dict:
    keys = ("species", "matrix", "assay", "direction", "route", "dose", "dose_unit", "regimen",
            "analyte", "measurement_subtype", "temperature", "ph", "concentration",
            "concentration_unit", "study_id", "batch_id", "assay_id", "experiment_date", "fed_fasted", "steady_state")
    context = {key: payload[key] for key in keys if payload.get(key) not in (None, "")}
    context["internal_notes"] = str(payload.get("notes") or "")
    return context


def _fingerprint(version: CompoundVersion, endpoint: str, value: object, unit: str, context: dict) -> str:
    material = {"version": version.inchikey, "endpoint": endpoint, "value": str(value), "unit": unit,
                "species": context.get("species"), "route": context.get("route"), "dose": context.get("dose"),
                "dose_unit": context.get("dose_unit"), "regimen": context.get("regimen"),
                "analyte": context.get("analyte"), "measurement_subtype": context.get("measurement_subtype"),
                "study_id": context.get("study_id"), "batch_id": context.get("batch_id")}
    return hashlib.sha256(json.dumps(material, sort_keys=True, default=str).encode()).hexdigest()


def _qualification(mapped: dict, context: dict) -> dict:
    status = mapped.get("comparability_status")
    direct = status in {"DIRECTLY_COMPARABLE", "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION"}
    context_ok = bool(mapped.get("normalized_value") is not None and status not in {"UNSUPPORTED", "NOT_COMPARABLE", "CONDITIONALLY_COMPARABLE"})
    reason = mapped.get("reason") or ("Required endpoint context is incomplete" if not context_ok else "")
    return {"qualification_version": "drugopt-experimental-qualification-v4",
            "stages": {"IDENTITY_QUALIFIED": True, "REFERENCE_QUALIFIED": True,
                       "NUMERIC_QUALIFIED": mapped.get("normalized_value") is not None,
                       "ENDPOINT_QUALIFIED": mapped.get("canonical_endpoint_id") not in {"", "UNRESOLVED"},
                       "CONTEXT_QUALIFIED": context_ok, "DIRECTLY_COMPARABLE": direct,
                       "IMPORTABLE": False, "ADAPTATION_ELIGIBLE": False},
            "context_status": "CONTEXT_QUALIFIED" if context_ok else "CONTEXT_NOT_QUALIFIED",
            "primary_gap_reason": reason}


def save_internal_evidence(db: Session, project_id: int, compound_id: int, payload: dict, *, supersedes: ExternalExperimentalEvidence | None = None) -> ExternalExperimentalEvidence:
    compound = db.get(Compound, compound_id)
    if not compound or compound.project_id != project_id:
        raise HTTPException(status_code=404, detail="Compound not found in project")
    version = _current_version(compound)
    raw_endpoint = _raw_endpoint(payload)
    raw_value = payload.get("raw_value", payload.get("value"))
    raw_unit = str(payload.get("raw_unit", payload.get("unit", ""))).strip()
    if not raw_endpoint or raw_value in (None, ""):
        raise HTTPException(status_code=400, detail="Canonical endpoint and raw value are required")
    context = _context(payload)
    # Preserve explicit endpoint-specific semantics in the raw source label.
    subtype = str(payload.get("measurement_type") or context.get("measurement_subtype") or "").strip()
    if subtype and subtype.lower() not in raw_endpoint.lower():
        raw_endpoint = f"{raw_endpoint} {subtype}"
    mapped = normalize_experimental_observation(raw_endpoint, raw_value, raw_unit,
        species=context.get("species", ""), context=context, assay_type=context.get("assay", ""),
        canonical_hint=payload.get("canonical_endpoint_id", ""))
    # Activity values are only scientifically meaningful against a project
    # assay with the same measurement semantics.  Keeping this association in
    # the shared evidence row lets the canonical comparison view place the
    # manual value beside the correct frozen activity prediction.
    if mapped.get("section") == "ACTIVITY":
        assay_id = context.get("assay_id")
        assay = db.get(AssayDefinition, int(assay_id)) if str(assay_id or "").isdigit() else None
        if not assay or assay.project_id != project_id:
            raise HTTPException(status_code=400, detail="A project assay is required for manual activity evidence")
        expected = str(assay.measurement_type or "").upper()
        observed = str(mapped.get("measurement_subtype") or "").upper()
        if expected and observed and expected != observed:
            raise HTTPException(status_code=400, detail=f"Selected assay measures {expected}; cannot save {observed} as the same activity endpoint")
        context["assay"] = assay.name
        context["target"] = assay.target
    qualification = _qualification(mapped, context)
    fingerprint = _fingerprint(version, mapped["canonical_endpoint_id"], raw_value, raw_unit, context)
    if supersedes is None:
        existing = db.scalar(select(ExternalExperimentalEvidence).where(
            ExternalExperimentalEvidence.provenance_key == fingerprint,
            ExternalExperimentalEvidence.lifecycle_status == "ACTIVE"))
        if existing:
            return existing
    now = datetime.now(timezone.utc)
    row = ExternalExperimentalEvidence(
        compound_version_id=version.id, provenance_key=hashlib.sha256((fingerprint + uuid.uuid4().hex).encode()).hexdigest() if supersedes else fingerprint,
        cas_number=compound.cas_number or "", raw_endpoint_name=raw_endpoint, raw_value=str(raw_value),
        raw_relation=str(payload.get("relation") or "="), raw_unit=raw_unit,
        assay_type=str(payload.get("measurement_type") or context.get("assay") or ""),
        assay_conditions_json=context, species=str(context.get("species") or ""),
        source_database="Internal Experimental", source_record_id=str(context.get("study_id") or context.get("batch_id") or f"internal-{uuid.uuid4().hex}"),
        source_assay_id=str(context.get("assay_id") or context.get("batch_id") or ""), source_document_id=str(context.get("study_id") or ""),
        reference_text=str(context.get("study_id") or "Internal experimental study"), source_url="",
        identity_match_status="EXACT_STRUCTURE_MATCH", endpoint_match_status="CANONICAL_MANUAL_ENTRY",
        mapping_status="CANONICAL_SCIENTIFIC_EVIDENCE", evidence_origin="INTERNAL_EXPERIMENTAL",
        canonical_endpoint_id=mapped["canonical_endpoint_id"], normalized_value="" if mapped.get("normalized_value") is None else str(mapped["normalized_value"]),
        normalized_unit=str(mapped.get("normalized_unit") or ""), normalization_rule=str(mapped.get("normalization_rule") or ""),
        normalization_version="drugopt-experimental-normalization-v1", comparability_status=str(mapped.get("comparability_status") or "UNSUPPORTED"),
        source_quality_class="A", duplicate_status="DISTINCT_MEASUREMENT", provenance_fingerprint=fingerprint,
        evidence_state="INTERNAL_EXPERIMENTAL", first_seen_at=now, last_seen_at=now, accepted_at=now,
        qualification_version="drugopt-experimental-qualification-v4", routing_version="drugopt-canonical-endpoint-v1",
        canonical_endpoint_version=CANONICAL_ENDPOINT_VERSION, unit_normalization_version=COMPARISON_UNIT_VERSION,
        display_evidence_group_id=f"internal-{fingerprint[:20]}", independent_experiment_group_id=str(context.get("study_id") or context.get("batch_id") or f"internal-{fingerprint[:20]}"),
        qualification_json=qualification, qualification_status=qualification["context_status"], routing_section=mapped.get("section", "UNCLASSIFIED"), routing_reason=qualification["primary_gap_reason"],
        retrieved_at=now, imported_at=now, lifecycle_status="ACTIVE", revision_number=(supersedes.revision_number + 1 if supersedes else 1),
        supersedes_evidence_id=supersedes.id if supersedes else None, updated_at=now)
    if supersedes:
        supersedes.lifecycle_status = "SUPERSEDED"
    db.add(row); db.flush()
    return row
