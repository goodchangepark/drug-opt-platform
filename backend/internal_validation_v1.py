"""Internal Prospective Validation Framework — Engine v1.

Stage 6 — Internal Prospective Validation of Prediction Engine v1.

Scientific principle:
    STRUCTURE → ENGINE V1 PREDICTION → IMMUTABLE FREEZE → EXPERIMENT
    → EXPERIMENTAL RESULT → COMPARISON

This module provides:
  - InternalValidationCampaign: campaign lifecycle entity
  - InternalValidationObservation: immutable prediction↔experiment linkage
  - InternalValidationExperimentalRecord: raw experimental data (immutable)
  - InternalValidationPredictionFreeze: v1-campaign-scoped prediction freeze
  - Campaign/observation registration helpers
  - Prediction-before-experiment ordering enforcement
  - Evidence classification (TRUE_PROSPECTIVE / BLINDED_RETROSPECTIVE /
    HISTORICAL_VISIBLE)
  - Experimental import with endpoint compatibility verification
  - Blinding enforcement: experimental values are NOT readable during
    prediction freeze creation

Forbidden in this module:
  - Model fitting, retraining, recalibration
  - Threshold modification
  - AD threshold modification
  - Shadow promotion
  - Bias correction fitting
  - Policy hash change
"""

from __future__ import annotations

import hashlib
import json
import uuid
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

from sqlalchemy import Boolean, DateTime, Float, Integer, String, Text
from sqlalchemy.orm import Mapped, mapped_column
from sqlalchemy.types import JSON

from .database import Base

ENGINE_V1_POLICY_ID = "drugopt-prediction-engine-v1"
ENGINE_V1_POLICY_VERSION = "drugopt-prediction-engine-v1@1.0.0"
ENGINE_V1_POLICY_HASH = (
    "12757ab197b5a70d8ea1754678d9a342ab0b6ea0d82f2896bebb767d686bbdeb"
)
STANDARDIZER_VERSION = "CHEM_STANDARDIZER_V1"

# Internal Validation Campaign v1
CAMPAIGN_ID = "IVC-engine-v1-2026-08-29"
CAMPAIGN_PROTOCOL_ID = "internal-validation-v1-protocol-2026-08-29"

# Evidence classification
EVIDENCE_TRUE_PROSPECTIVE = "TRUE_PROSPECTIVE"
EVIDENCE_BLINDED_RETROSPECTIVE = "BLINDED_RETROSPECTIVE"
EVIDENCE_HISTORICAL_VISIBLE = "HISTORICAL_VISIBLE"

# Endpoint compatibility
COMPAT_DIRECT_MATCH = "DIRECT_MATCH"
COMPAT_UNIT_TRANSFORM = "DETERMINISTIC_UNIT_TRANSFORM"
COMPAT_ASSAY_LIMITED = "ASSAY_CONTEXT_LIMITED"
COMPAT_MISMATCH = "ENDPOINT_MISMATCH"
COMPAT_NOT_COMPARABLE = "NOT_COMPARABLE"

# Campaign statuses
CAMPAIGN_STATUS_ACTIVE = "ACTIVE"
CAMPAIGN_STATUS_COLLECTING = "COLLECTING"
CAMPAIGN_STATUS_PARTIAL = "PARTIAL"
CAMPAIGN_STATUS_SUFFICIENT = "SUFFICIENT_FOR_BASELINE_ASSESSMENT"


def utcnow() -> datetime:
    return datetime.now(timezone.utc)


def _sha256(payload: str) -> str:
    return hashlib.sha256(payload.encode()).hexdigest()


def _new_id(prefix: str) -> str:
    return f"{prefix}-{uuid.uuid4().hex[:16]}"


# ---------------------------------------------------------------------------
# ORM Models
# ---------------------------------------------------------------------------


class InternalValidationCampaignRow(Base):
    """Campaign lifecycle entity for Engine v1 internal validation.

    Distinct from ordinary research projects.  Status is driven by data
    collection progress, not by code deployment.
    """

    __tablename__ = "internal_validation_campaigns"

    campaign_id: Mapped[str] = mapped_column(String(120), primary_key=True)
    name: Mapped[str] = mapped_column(String(240))
    protocol_id: Mapped[str] = mapped_column(String(120))
    engine_policy_id: Mapped[str] = mapped_column(String(120))
    engine_policy_version: Mapped[str] = mapped_column(String(120))
    engine_policy_hash: Mapped[str] = mapped_column(String(64))
    standardizer_version: Mapped[str] = mapped_column(String(80))
    # Framework/scientific status
    framework_status: Mapped[str] = mapped_column(String(60), default="READY")
    scientific_status: Mapped[str] = mapped_column(
        String(80), default="COLLECTING"
    )
    status: Mapped[str] = mapped_column(String(60), default=CAMPAIGN_STATUS_ACTIVE)
    compound_count: Mapped[int] = mapped_column(Integer, default=0)
    endpoint_count: Mapped[int] = mapped_column(Integer, default=0)
    prediction_freeze_complete: Mapped[bool] = mapped_column(Boolean, default=False)
    experiment_import_complete: Mapped[bool] = mapped_column(Boolean, default=False)
    analysis_complete: Mapped[bool] = mapped_column(Boolean, default=False)
    bootstrap_seed: Mapped[int] = mapped_column(Integer, default=42)
    notes: Mapped[str] = mapped_column(Text, default="")
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=utcnow
    )
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=utcnow, onupdate=utcnow
    )


class InternalValidationCohortEntryRow(Base):
    """One compound enrolled in a validation campaign.

    Records compound identity at enrollment time (structure hash, version).
    """

    __tablename__ = "internal_validation_cohort_entries"

    entry_id: Mapped[str] = mapped_column(String(80), primary_key=True)
    campaign_id: Mapped[str] = mapped_column(String(120), index=True)
    compound_version_id: Mapped[str] = mapped_column(String(80), index=True)
    # De-identified label (no confidential names in public artifacts)
    compound_label: Mapped[str] = mapped_column(String(120))
    compound_identifier: Mapped[str] = mapped_column(String(120))
    # Structure fingerprint (no raw SMILES stored in the model)
    inchikey: Mapped[str] = mapped_column(String(32))
    structure_hash: Mapped[str] = mapped_column(String(64))
    standardizer_version: Mapped[str] = mapped_column(String(80))
    project_label: Mapped[str] = mapped_column(String(120), default="")
    chemical_series_label: Mapped[str] = mapped_column(String(120), default="")
    murcko_scaffold_hash: Mapped[str] = mapped_column(String(64), default="")
    enrolled_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=utcnow
    )
    eligibility_status: Mapped[str] = mapped_column(String(60), default="ELIGIBLE")
    eligibility_notes: Mapped[str] = mapped_column(Text, default="")
    record_hash: Mapped[str] = mapped_column(String(64), unique=True)


class InternalValidationPredictionFreezeRow(Base):
    """Engine v1 prediction freeze created for a validation campaign compound.

    This is a campaign-scoped view of the qualification_prediction_freezes
    record.  It must be created BEFORE experimental values are imported.
    The freeze_timestamp ordering is enforced by the framework.
    """

    __tablename__ = "internal_validation_prediction_freezes"

    vfreeze_id: Mapped[str] = mapped_column(String(120), primary_key=True)
    campaign_id: Mapped[str] = mapped_column(String(120), index=True)
    entry_id: Mapped[str] = mapped_column(String(80), index=True)
    compound_version_id: Mapped[str] = mapped_column(String(80), index=True)
    # Link to upstream qualification_prediction_freezes record
    upstream_frozen_prediction_id: Mapped[str] = mapped_column(String(80), index=True)
    inchikey: Mapped[str] = mapped_column(String(32))
    structure_hash: Mapped[str] = mapped_column(String(64))
    standardizer_version: Mapped[str] = mapped_column(String(80))
    engine_policy_version: Mapped[str] = mapped_column(String(120))
    engine_policy_hash: Mapped[str] = mapped_column(String(64))
    endpoint_id: Mapped[str] = mapped_column(String(120), index=True)
    strategy: Mapped[str] = mapped_column(String(80))
    core_model_id: Mapped[str] = mapped_column(String(240), default="")
    core_model_version: Mapped[str] = mapped_column(String(120), default="")
    evidence_class: Mapped[str] = mapped_column(String(80))
    prediction_value: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    probability: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    unit: Mapped[str] = mapped_column(String(60))
    applicability_domain: Mapped[str] = mapped_column(String(60))
    reliability: Mapped[str] = mapped_column(String(60), default="")
    limitations_json: Mapped[list] = mapped_column(JSON, default=list)
    shadow_outputs_json: Mapped[dict] = mapped_column(JSON, default=dict)
    # Immutable timestamp — must exist before any experiment import
    freeze_timestamp: Mapped[datetime] = mapped_column(DateTime(timezone=True))
    record_hash: Mapped[str] = mapped_column(String(64), unique=True)
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=utcnow
    )


class InternalValidationExperimentalRecordRow(Base):
    """Raw experimental measurement.

    Raw values are immutable once stored.  Normalized values are computed
    and stored separately as pairing fields.  Censored/qualifier values
    are preserved exactly.

    result_available_at is the timestamp at which the experimental result
    first became available to the team (not just assay date).
    """

    __tablename__ = "internal_validation_experimental_records"

    exp_record_id: Mapped[str] = mapped_column(String(80), primary_key=True)
    campaign_id: Mapped[str] = mapped_column(String(120), index=True)
    compound_version_id: Mapped[str] = mapped_column(String(80), index=True)
    inchikey: Mapped[str] = mapped_column(String(32))
    structure_hash: Mapped[str] = mapped_column(String(64))
    endpoint_id: Mapped[str] = mapped_column(String(120), index=True)
    # Raw observation (never overwritten)
    raw_value: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    raw_unit: Mapped[str] = mapped_column(String(60))
    qualifier: Mapped[str] = mapped_column(
        String(20), default=""
    )  # "<", ">", "=", "~", "BLQ", "ULOQ"
    species: Mapped[str] = mapped_column(String(120), default="")
    assay_type: Mapped[str] = mapped_column(String(240), default="")
    assay_direction: Mapped[str] = mapped_column(String(60), default="")
    assay_ph: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    assay_protocol: Mapped[str] = mapped_column(String(240), default="")
    replicate_id: Mapped[str] = mapped_column(String(80), default="")
    assay_date: Mapped[str] = mapped_column(String(40), default="")
    result_available_at: Mapped[Optional[datetime]] = mapped_column(
        DateTime(timezone=True), nullable=True
    )
    source: Mapped[str] = mapped_column(Text, default="")
    censor_flag: Mapped[bool] = mapped_column(Boolean, default=False)
    # Endpoint compatibility classification
    endpoint_compatibility: Mapped[str] = mapped_column(
        String(60), default=COMPAT_DIRECT_MATCH
    )
    compatibility_notes: Mapped[str] = mapped_column(Text, default="")
    # Import audit
    imported_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=utcnow
    )
    record_hash: Mapped[str] = mapped_column(String(64), unique=True)


class InternalValidationObservationRow(Base):
    """Immutable linkage between one prediction freeze and one experimental result.

    Evidence classification:
      TRUE_PROSPECTIVE        — freeze_timestamp < result_available_at
      BLINDED_RETROSPECTIVE   — result existed but was hidden during freeze
      HISTORICAL_VISIBLE      — result was visible during Engine v1 development

    Only DIRECT_MATCH and DETERMINISTIC_UNIT_TRANSFORM enter primary metrics.
    """

    __tablename__ = "internal_validation_observations"

    observation_id: Mapped[str] = mapped_column(String(80), primary_key=True)
    campaign_id: Mapped[str] = mapped_column(String(120), index=True)
    entry_id: Mapped[str] = mapped_column(String(80), index=True)
    compound_version_id: Mapped[str] = mapped_column(String(80), index=True)
    inchikey: Mapped[str] = mapped_column(String(32))
    structure_hash: Mapped[str] = mapped_column(String(64))
    endpoint_id: Mapped[str] = mapped_column(String(120), index=True)
    # Prediction side
    vfreeze_id: Mapped[str] = mapped_column(String(120), index=True)
    upstream_frozen_prediction_id: Mapped[str] = mapped_column(String(80))
    engine_policy_version: Mapped[str] = mapped_column(String(120))
    engine_policy_hash: Mapped[str] = mapped_column(String(64))
    strategy: Mapped[str] = mapped_column(String(80))
    evidence_class: Mapped[str] = mapped_column(String(80))
    prediction_value: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    prediction_unit: Mapped[str] = mapped_column(String(60))
    core_model_id: Mapped[str] = mapped_column(String(240), default="")
    core_model_version: Mapped[str] = mapped_column(String(120), default="")
    applicability_domain: Mapped[str] = mapped_column(String(60))
    reliability: Mapped[str] = mapped_column(String(60), default="")
    freeze_timestamp: Mapped[datetime] = mapped_column(DateTime(timezone=True))
    # Experimental side
    exp_record_id: Mapped[str] = mapped_column(String(80), index=True)
    raw_value: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    raw_unit: Mapped[str] = mapped_column(String(60))
    qualifier: Mapped[str] = mapped_column(String(20), default="")
    censor_flag: Mapped[bool] = mapped_column(Boolean, default=False)
    assay_type: Mapped[str] = mapped_column(String(240), default="")
    species: Mapped[str] = mapped_column(String(120), default="")
    result_available_at: Mapped[Optional[datetime]] = mapped_column(
        DateTime(timezone=True), nullable=True
    )
    # Pairing classification
    endpoint_compatibility: Mapped[str] = mapped_column(String(60))
    prospective_evidence_class: Mapped[str] = mapped_column(
        String(60)
    )  # TRUE_PROSPECTIVE | BLINDED_RETROSPECTIVE | HISTORICAL_VISIBLE
    enters_primary_metrics: Mapped[bool] = mapped_column(Boolean, default=False)
    # Comparison (computed at analysis time — not at import time)
    comparison_status: Mapped[str] = mapped_column(
        String(60), default="PENDING"
    )
    absolute_error: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    signed_error: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    classification_correct: Mapped[Optional[bool]] = mapped_column(
        Boolean, nullable=True
    )
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=utcnow
    )
    record_hash: Mapped[str] = mapped_column(String(64), unique=True)


class InternalValidationAnalysisRow(Base):
    """Snapshot of analysis results for a campaign + endpoint.

    Overwritten each time analysis is re-run (not append-only).
    The raw observations are the immutable source of truth.
    """

    __tablename__ = "internal_validation_analyses"

    analysis_id: Mapped[str] = mapped_column(String(80), primary_key=True)
    campaign_id: Mapped[str] = mapped_column(String(120), index=True)
    endpoint_id: Mapped[str] = mapped_column(String(120), index=True)
    analysis_type: Mapped[str] = mapped_column(String(60))  # REGRESSION | CLASSIFICATION | DESCRIPTIVE
    n_total: Mapped[int] = mapped_column(Integer, default=0)
    n_prospective: Mapped[int] = mapped_column(Integer, default=0)
    n_blinded_retro: Mapped[int] = mapped_column(Integer, default=0)
    n_historical: Mapped[int] = mapped_column(Integer, default=0)
    n_primary_metrics: Mapped[int] = mapped_column(Integer, default=0)
    metrics_json: Mapped[dict] = mapped_column(JSON, default=dict)
    bootstrap_json: Mapped[dict] = mapped_column(JSON, default=dict)
    ad_analysis_json: Mapped[dict] = mapped_column(JSON, default=dict)
    reliability_analysis_json: Mapped[dict] = mapped_column(JSON, default=dict)
    shadow_disagreement_json: Mapped[dict] = mapped_column(JSON, default=dict)
    scaffold_analysis_json: Mapped[dict] = mapped_column(JSON, default=dict)
    data_insufficiency_flag: Mapped[bool] = mapped_column(Boolean, default=False)
    data_insufficiency_reason: Mapped[str] = mapped_column(Text, default="")
    computed_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=utcnow
    )


# ---------------------------------------------------------------------------
# Schema helpers
# ---------------------------------------------------------------------------

_VALIDATION_TABLES = (
    InternalValidationCampaignRow,
    InternalValidationCohortEntryRow,
    InternalValidationPredictionFreezeRow,
    InternalValidationExperimentalRecordRow,
    InternalValidationObservationRow,
    InternalValidationAnalysisRow,
)


def ensure_validation_schema(engine_obj) -> None:
    """Create validation tables if they do not exist.

    Safe to call with a running service — uses CREATE TABLE IF NOT EXISTS.
    Does not drop or alter existing tables.
    """
    from sqlalchemy import inspect, text

    inspector = inspect(engine_obj)
    existing = set(inspector.get_table_names())
    for model in _VALIDATION_TABLES:
        tname = model.__tablename__
        if tname not in existing:
            model.__table__.create(engine_obj, checkfirst=True)


# ---------------------------------------------------------------------------
# Blinding enforcement
# ---------------------------------------------------------------------------


def assert_no_experimental_access(
    session,
    compound_version_id: str,
    endpoint_id: str,
) -> None:
    """Raise RuntimeError if an experimental record already exists for this
    compound+endpoint in a way that would contaminate prediction generation.

    This enforces: prediction MUST be stored before experimental import.
    """
    from sqlalchemy import select

    existing = session.scalars(
        select(InternalValidationExperimentalRecordRow).where(
            InternalValidationExperimentalRecordRow.compound_version_id
            == compound_version_id,
            InternalValidationExperimentalRecordRow.endpoint_id == endpoint_id,
        )
    ).first()
    if existing is not None:
        raise RuntimeError(
            f"BLINDING VIOLATION: experimental record {existing.exp_record_id} "
            f"for compound_version_id={compound_version_id}, endpoint={endpoint_id} "
            f"already exists before prediction freeze. "
            f"Prediction must precede experimental import."
        )


# ---------------------------------------------------------------------------
# Evidence classification
# ---------------------------------------------------------------------------


def classify_evidence(
    freeze_timestamp: datetime,
    result_available_at: Optional[datetime],
    blinded_retrospective: bool = False,
) -> str:
    """Classify a paired observation by evidence type.

    TRUE_PROSPECTIVE: freeze happened before result was available.
    BLINDED_RETROSPECTIVE: result existed but was hidden during freeze,
        and that blinding is documented.
    HISTORICAL_VISIBLE: result was available and potentially visible during
        model development.
    """
    if result_available_at is None:
        # No timestamp for result — cannot classify as prospective
        return EVIDENCE_HISTORICAL_VISIBLE

    # Normalize to UTC-aware
    ft = freeze_timestamp.replace(tzinfo=timezone.utc) if freeze_timestamp.tzinfo is None else freeze_timestamp
    rat = result_available_at.replace(tzinfo=timezone.utc) if result_available_at.tzinfo is None else result_available_at

    if ft < rat:
        return EVIDENCE_TRUE_PROSPECTIVE
    elif blinded_retrospective:
        return EVIDENCE_BLINDED_RETROSPECTIVE
    else:
        return EVIDENCE_HISTORICAL_VISIBLE


# ---------------------------------------------------------------------------
# Endpoint compatibility verification
# ---------------------------------------------------------------------------

# Known endpoint compatibility contracts
_ENDPOINT_COMPAT_RULES: Dict[str, Dict[str, Any]] = {
    "solubility_aqueous_logs": {
        "required_unit_raw": ["log10(mol/L)", "log S", "logS", "mol/L"],
        "required_assay_type_keywords": [],  # kinetic vs thermodynamic noted
        "notes": "Kinetic vs thermodynamic distinction. Assay pH should match.",
    },
    "permeability_caco2_logpapp": {
        "required_unit_raw": ["log10(cm/s)", "cm/s", "log(cm/s)"],
        "required_assay_direction": ["A→B", "A-B", "AB", "A to B"],
        "notes": "Papp A→B only. Unit must be log10(cm/s) or convertible.",
    },
    "ppb_human_percent_bound": {
        "required_species": ["human", "Human"],
        "required_unit_raw": ["% bound", "%bound", "fraction bound"],
        "notes": "Human PPB only. fu vs percent bound conversion is deterministic.",
    },
    "hlm_intrinsic_clearance_scaled_log10": {
        "required_species": ["human", "Human"],
        "required_assay_type_keywords": ["microsom", "HLM"],
        "notes": "Human liver microsomal. Clint log10(mL/min/kg). Hepatocyte ≠ microsomal.",
    },
    "rlm_intrinsic_clearance_scaled_log10": {
        "required_species": ["rat", "Rat"],
        "required_assay_type_keywords": ["microsom", "RLM"],
        "notes": "Rat liver microsomal. Strict species isolation.",
    },
    "mlm_intrinsic_clearance_scaled_log10": {
        "required_species": ["mouse", "Mouse"],
        "required_assay_type_keywords": ["microsom", "MLM"],
        "notes": "Mouse liver microsomal. Strict species isolation.",
    },
    "safety_herg_blocker_prob": {
        "notes": "hERG IC50 or patch-clamp. Binary threshold at 1 µM or 10 µM must be specified.",
    },
    "cyp3a4_inhibitor_prob": {
        "notes": "CYP3A4 inhibitor (not substrate). IC50-based or direct inhibition assay.",
    },
    "ionization_pka_estimated": {
        "notes": "pKa: compare measured vs rule estimate only. ±1–2 pKa unit uncertainty expected.",
    },
    "physchem_logd_ph74_derived_estimate": {
        "notes": "logD pH 7.4 comparison: measured vs derived estimate. Not quantitative ML.",
    },
}


def check_endpoint_compatibility(
    endpoint_id: str,
    raw_unit: str,
    species: str = "",
    assay_type: str = "",
    assay_direction: str = "",
) -> Tuple[str, str]:
    """Return (compatibility_status, notes) for an experimental observation.

    Used during experimental import to classify each record before it enters
    the pairing pipeline.
    """
    rule = _ENDPOINT_COMPAT_RULES.get(endpoint_id)
    if rule is None:
        return COMPAT_ASSAY_LIMITED, f"No specific compatibility contract for {endpoint_id}"

    notes_parts = [rule.get("notes", "")]

    # Species check
    req_species = rule.get("required_species", [])
    if req_species and species and not any(s.lower() in species.lower() for s in req_species):
        return COMPAT_MISMATCH, f"Species mismatch: got '{species}', expected one of {req_species}"

    # Unit check
    req_units = rule.get("required_unit_raw", [])
    if req_units and raw_unit:
        unit_ok = any(u.lower() in raw_unit.lower() or raw_unit.lower() in u.lower() for u in req_units)
        if not unit_ok:
            return COMPAT_UNIT_TRANSFORM, f"Unit '{raw_unit}' may need conversion. Expected: {req_units}"

    # Assay direction check
    req_dir = rule.get("required_assay_direction", [])
    if req_dir and assay_direction:
        dir_ok = any(d.lower() in assay_direction.lower() for d in req_dir)
        if not dir_ok:
            return COMPAT_MISMATCH, f"Direction mismatch: got '{assay_direction}', expected one of {req_dir}"

    return COMPAT_DIRECT_MATCH, "; ".join(p for p in notes_parts if p)


# ---------------------------------------------------------------------------
# Non-positive log handling
# ---------------------------------------------------------------------------


def safe_log_error(
    prediction_log: Optional[float],
    experimental_raw: Optional[float],
    endpoint_id: str,
) -> Tuple[Optional[float], Optional[float], str]:
    """Return (absolute_error, signed_error, status) for a log-scale endpoint.

    If experimental value is non-positive (cannot be log-transformed),
    returns (None, None, 'NON_POSITIVE_EXCLUDED_FROM_LOG_METRIC').
    """
    if prediction_log is None or experimental_raw is None:
        return None, None, "MISSING_VALUE"

    # Endpoints that store already-log-transformed predictions
    log_endpoints = {
        "solubility_aqueous_logs",
        "permeability_caco2_logpapp",
        "hlm_intrinsic_clearance_scaled_log10",
        "rlm_intrinsic_clearance_scaled_log10",
        "mlm_intrinsic_clearance_scaled_log10",
    }

    if endpoint_id in log_endpoints:
        import math

        # Experimental value: check if it's already log-transformed or raw
        # Convention: if abs(experimental_raw) < 20 it is likely already log
        # But we rely on unit metadata. Here we assume the caller has verified.
        if experimental_raw <= 0:
            return None, None, "NON_POSITIVE_EXCLUDED_FROM_LOG_METRIC"
        try:
            exp_log = math.log10(experimental_raw)
        except ValueError:
            return None, None, "NON_POSITIVE_EXCLUDED_FROM_LOG_METRIC"

        signed = prediction_log - exp_log
        return abs(signed), signed, "OK"

    # For probability endpoints — direct difference (not log)
    signed = prediction_log - experimental_raw
    return abs(signed), signed, "OK"


# ---------------------------------------------------------------------------
# Campaign registration helpers
# ---------------------------------------------------------------------------


def get_or_create_campaign(session) -> InternalValidationCampaignRow:
    """Return the canonical Engine v1 internal validation campaign row,
    creating it if it does not yet exist.
    """
    from sqlalchemy import select

    row = session.scalars(
        select(InternalValidationCampaignRow).where(
            InternalValidationCampaignRow.campaign_id == CAMPAIGN_ID
        )
    ).first()

    if row is None:
        row = InternalValidationCampaignRow(
            campaign_id=CAMPAIGN_ID,
            name="Engine v1 Internal Prospective Validation",
            protocol_id=CAMPAIGN_PROTOCOL_ID,
            engine_policy_id=ENGINE_V1_POLICY_ID,
            engine_policy_version=ENGINE_V1_POLICY_VERSION,
            engine_policy_hash=ENGINE_V1_POLICY_HASH,
            standardizer_version=STANDARDIZER_VERSION,
            framework_status="READY",
            scientific_status="COLLECTING",
            status=CAMPAIGN_STATUS_ACTIVE,
            bootstrap_seed=42,
        )
        session.add(row)
        session.commit()
        session.refresh(row)

    return row


def register_cohort_entry(
    session,
    campaign_id: str,
    compound_version_id: str,
    compound_label: str,
    compound_identifier: str,
    inchikey: str,
    structure_hash: str,
    project_label: str = "",
    chemical_series_label: str = "",
    murcko_scaffold_hash: str = "",
) -> InternalValidationCohortEntryRow:
    """Enroll a compound in the validation campaign.

    Idempotent: returns existing row if compound_version_id already enrolled.
    """
    from sqlalchemy import select

    existing = session.scalars(
        select(InternalValidationCohortEntryRow).where(
            InternalValidationCohortEntryRow.campaign_id == campaign_id,
            InternalValidationCohortEntryRow.compound_version_id
            == compound_version_id,
        )
    ).first()

    if existing is not None:
        return existing

    entry_id = _new_id("CE")
    rec_hash = _sha256(
        f"{campaign_id}|{compound_version_id}|{inchikey}|{structure_hash}"
    )

    row = InternalValidationCohortEntryRow(
        entry_id=entry_id,
        campaign_id=campaign_id,
        compound_version_id=compound_version_id,
        compound_label=compound_label,
        compound_identifier=compound_identifier,
        inchikey=inchikey,
        structure_hash=structure_hash,
        standardizer_version=STANDARDIZER_VERSION,
        project_label=project_label,
        chemical_series_label=chemical_series_label,
        murcko_scaffold_hash=murcko_scaffold_hash,
        eligibility_status="ELIGIBLE",
        record_hash=rec_hash,
    )
    session.add(row)
    session.commit()
    session.refresh(row)
    return row


def register_prediction_freeze(
    session,
    campaign_id: str,
    entry_id: str,
    upstream_frozen_prediction_id: str,
    compound_version_id: str,
    inchikey: str,
    structure_hash: str,
    endpoint_id: str,
    strategy: str,
    evidence_class: str,
    prediction_value: Optional[float],
    probability: Optional[float],
    unit: str,
    applicability_domain: str,
    reliability: str,
    freeze_timestamp: datetime,
    core_model_id: str = "",
    core_model_version: str = "",
    limitations: Optional[List[str]] = None,
    shadow_outputs: Optional[Dict] = None,
) -> InternalValidationPredictionFreezeRow:
    """Register a prediction freeze for a validation compound+endpoint.

    Enforces blinding: raises if experimental record already exists.
    Idempotent by upstream_frozen_prediction_id.
    """
    from sqlalchemy import select

    # Blinding check
    assert_no_experimental_access(session, compound_version_id, endpoint_id)

    # Idempotency
    existing = session.scalars(
        select(InternalValidationPredictionFreezeRow).where(
            InternalValidationPredictionFreezeRow.upstream_frozen_prediction_id
            == upstream_frozen_prediction_id,
            InternalValidationPredictionFreezeRow.campaign_id == campaign_id,
        )
    ).first()
    if existing is not None:
        return existing

    vfreeze_id = _new_id("VF")
    rec_hash = _sha256(
        f"{campaign_id}|{upstream_frozen_prediction_id}|{endpoint_id}"
        f"|{ENGINE_V1_POLICY_HASH}"
    )

    row = InternalValidationPredictionFreezeRow(
        vfreeze_id=vfreeze_id,
        campaign_id=campaign_id,
        entry_id=entry_id,
        compound_version_id=compound_version_id,
        upstream_frozen_prediction_id=upstream_frozen_prediction_id,
        inchikey=inchikey,
        structure_hash=structure_hash,
        standardizer_version=STANDARDIZER_VERSION,
        engine_policy_version=ENGINE_V1_POLICY_VERSION,
        engine_policy_hash=ENGINE_V1_POLICY_HASH,
        endpoint_id=endpoint_id,
        strategy=strategy,
        evidence_class=evidence_class,
        prediction_value=prediction_value,
        probability=probability,
        unit=unit,
        applicability_domain=applicability_domain,
        reliability=reliability,
        freeze_timestamp=freeze_timestamp,
        core_model_id=core_model_id,
        core_model_version=core_model_version,
        limitations_json=limitations or [],
        shadow_outputs_json=shadow_outputs or {},
        record_hash=rec_hash,
    )
    session.add(row)
    session.commit()
    session.refresh(row)
    return row


def import_experimental_record(
    session,
    campaign_id: str,
    compound_version_id: str,
    inchikey: str,
    structure_hash: str,
    endpoint_id: str,
    raw_value: Optional[float],
    raw_unit: str,
    qualifier: str = "=",
    species: str = "",
    assay_type: str = "",
    assay_direction: str = "",
    assay_ph: Optional[float] = None,
    assay_protocol: str = "",
    replicate_id: str = "",
    assay_date: str = "",
    result_available_at: Optional[datetime] = None,
    source: str = "",
    censor_flag: bool = False,
) -> InternalValidationExperimentalRecordRow:
    """Import one raw experimental observation.

    Endpoint compatibility is verified before storage.
    Duplicate replicates are allowed (different replicate_id).
    """
    compat, compat_notes = check_endpoint_compatibility(
        endpoint_id, raw_unit, species, assay_type, assay_direction
    )

    exp_record_id = _new_id("EXP")
    rec_hash = _sha256(
        f"{campaign_id}|{compound_version_id}|{endpoint_id}"
        f"|{raw_value}|{raw_unit}|{replicate_id}|{assay_date}"
    )

    row = InternalValidationExperimentalRecordRow(
        exp_record_id=exp_record_id,
        campaign_id=campaign_id,
        compound_version_id=compound_version_id,
        inchikey=inchikey,
        structure_hash=structure_hash,
        endpoint_id=endpoint_id,
        raw_value=raw_value,
        raw_unit=raw_unit,
        qualifier=qualifier,
        species=species,
        assay_type=assay_type,
        assay_direction=assay_direction,
        assay_ph=assay_ph,
        assay_protocol=assay_protocol,
        replicate_id=replicate_id,
        assay_date=assay_date,
        result_available_at=result_available_at,
        source=source,
        censor_flag=censor_flag,
        endpoint_compatibility=compat,
        compatibility_notes=compat_notes,
        record_hash=rec_hash,
    )
    session.add(row)
    session.commit()
    session.refresh(row)
    return row


def pair_observation(
    session,
    campaign_id: str,
    vfreeze: InternalValidationPredictionFreezeRow,
    exp_record: InternalValidationExperimentalRecordRow,
    blinded_retrospective_documented: bool = False,
) -> InternalValidationObservationRow:
    """Create an immutable prediction↔experiment pairing.

    Enforces freeze_timestamp < result_available_at for TRUE_PROSPECTIVE.
    Evidence class and primary metrics eligibility are computed here.
    """
    from sqlalchemy import select

    # Idempotency
    existing = session.scalars(
        select(InternalValidationObservationRow).where(
            InternalValidationObservationRow.vfreeze_id == vfreeze.vfreeze_id,
            InternalValidationObservationRow.exp_record_id
            == exp_record.exp_record_id,
        )
    ).first()
    if existing is not None:
        return existing

    ev_class = classify_evidence(
        vfreeze.freeze_timestamp,
        exp_record.result_available_at,
        blinded_retrospective=blinded_retrospective_documented,
    )

    # Primary metrics eligibility
    compat = exp_record.endpoint_compatibility
    enters_primary = compat in (COMPAT_DIRECT_MATCH, COMPAT_UNIT_TRANSFORM) and not exp_record.censor_flag

    obs_id = _new_id("OBS")
    rec_hash = _sha256(
        f"{campaign_id}|{vfreeze.vfreeze_id}|{exp_record.exp_record_id}"
    )

    row = InternalValidationObservationRow(
        observation_id=obs_id,
        campaign_id=campaign_id,
        entry_id=vfreeze.entry_id,
        compound_version_id=vfreeze.compound_version_id,
        inchikey=vfreeze.inchikey,
        structure_hash=vfreeze.structure_hash,
        endpoint_id=vfreeze.endpoint_id,
        vfreeze_id=vfreeze.vfreeze_id,
        upstream_frozen_prediction_id=vfreeze.upstream_frozen_prediction_id,
        engine_policy_version=ENGINE_V1_POLICY_VERSION,
        engine_policy_hash=ENGINE_V1_POLICY_HASH,
        strategy=vfreeze.strategy,
        evidence_class=vfreeze.evidence_class,
        prediction_value=vfreeze.prediction_value,
        prediction_unit=vfreeze.unit,
        core_model_id=vfreeze.core_model_id,
        core_model_version=vfreeze.core_model_version,
        applicability_domain=vfreeze.applicability_domain,
        reliability=vfreeze.reliability,
        freeze_timestamp=vfreeze.freeze_timestamp,
        exp_record_id=exp_record.exp_record_id,
        raw_value=exp_record.raw_value,
        raw_unit=exp_record.raw_unit,
        qualifier=exp_record.qualifier,
        censor_flag=exp_record.censor_flag,
        assay_type=exp_record.assay_type,
        species=exp_record.species,
        result_available_at=exp_record.result_available_at,
        endpoint_compatibility=compat,
        prospective_evidence_class=ev_class,
        enters_primary_metrics=enters_primary,
        comparison_status="PENDING",
        record_hash=rec_hash,
    )
    session.add(row)
    session.commit()
    session.refresh(row)
    return row


# ---------------------------------------------------------------------------
# Campaign summary helper
# ---------------------------------------------------------------------------


def campaign_summary(session, campaign_id: str) -> Dict[str, Any]:
    """Return a summary dict for the campaign (no experimental values exposed)."""
    from sqlalchemy import func, select

    campaign = session.scalars(
        select(InternalValidationCampaignRow).where(
            InternalValidationCampaignRow.campaign_id == campaign_id
        )
    ).first()
    if campaign is None:
        return {"error": f"Campaign {campaign_id} not found"}

    n_compounds = session.scalar(
        select(func.count(InternalValidationCohortEntryRow.entry_id)).where(
            InternalValidationCohortEntryRow.campaign_id == campaign_id
        )
    )
    n_freezes = session.scalar(
        select(
            func.count(InternalValidationPredictionFreezeRow.vfreeze_id)
        ).where(
            InternalValidationPredictionFreezeRow.campaign_id == campaign_id
        )
    )
    n_experiments = session.scalar(
        select(
            func.count(
                InternalValidationExperimentalRecordRow.exp_record_id
            )
        ).where(
            InternalValidationExperimentalRecordRow.campaign_id == campaign_id
        )
    )
    n_observations = session.scalar(
        select(
            func.count(InternalValidationObservationRow.observation_id)
        ).where(
            InternalValidationObservationRow.campaign_id == campaign_id
        )
    )

    return {
        "campaign_id": campaign.campaign_id,
        "name": campaign.name,
        "protocol_id": campaign.protocol_id,
        "engine_policy_version": campaign.engine_policy_version,
        "engine_policy_hash": campaign.engine_policy_hash,
        "framework_status": campaign.framework_status,
        "scientific_status": campaign.scientific_status,
        "status": campaign.status,
        "n_compounds_enrolled": n_compounds,
        "n_prediction_freezes": n_freezes,
        "n_experimental_records": n_experiments,
        "n_paired_observations": n_observations,
        "prediction_freeze_complete": campaign.prediction_freeze_complete,
        "experiment_import_complete": campaign.experiment_import_complete,
        "analysis_complete": campaign.analysis_complete,
        "created_at": campaign.created_at.isoformat() if campaign.created_at else None,
    }
