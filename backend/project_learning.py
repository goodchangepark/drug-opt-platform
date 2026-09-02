"""Continuous project-learning lifecycle helpers.

Policy Version: drugopt-learning-observation-policy-v1

The functions in this module record comparisons, validation classification, and
learning eligibility. They never rewrite a frozen prediction and never activate
an adapter without explicit user action.

Validation types:
- PROSPECTIVE_VALIDATION: Prediction frozen strictly before experimental result existed.
- RETROSPECTIVE_OUT_OF_FOLD_VALIDATION: Historical public/retrospective evidence evaluated
  via leakage-safe leave-one-compound-out (LOCO) or out-of-fold (OOF) evaluation.
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Set, Tuple

from sqlalchemy import select
from sqlalchemy.orm import Session

from .admet import ADMETMeasurement, ADMETModelRegistry, ADMETPrediction, PredictionExperimentalPairRecord
from .models import Compound, CompoundVersion, ExternalExperimentalEvidence

LEARNING_OBSERVATION_POLICY_VERSION = "drugopt-learning-observation-policy-v1"
DIRECT_COMPARABILITY = {"DIRECTLY_COMPARABLE", "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION"}


def _aware(value):
    if value is None:
        return None
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value


def _number(value):
    try:
        return float(value) if value is not None and str(value).strip() else None
    except (TypeError, ValueError):
        return None


def _model_key(prediction):
    name = str(prediction.model.model_name or "").lower()
    if "admetica" in name:
        return "admetica_solubility"
    if "esol" in name:
        return "esol_delaney_v1"
    if prediction.model.endpoint_name == "Solubility":
        return "rdkit_gbr_solubility_v1"
    return f"{prediction.model.endpoint_name}:{prediction.model_id}"


def _prediction_bundle(db: Session, version_id: int, endpoint_name: str, before=None):
    """Return the latest pre-experiment model rows (or frozen base rows) and bundle."""
    query = select(ADMETPrediction).join(ADMETModelRegistry).where(
        ADMETPrediction.version_id == version_id,
        ADMETModelRegistry.endpoint_name == endpoint_name,
        ADMETPrediction.execution_status == "SUCCESS",
    )
    if before is not None:
        # Check if pre-experiment predictions exist
        pre_query = query.where(ADMETPrediction.created_at < before).order_by(ADMETPrediction.created_at.desc())
        rows = list(db.scalars(pre_query).all())
        if not rows:
            # Fall back to latest base prediction for retrospective evaluation
            rows = list(db.scalars(query.order_by(ADMETPrediction.created_at.desc())).all())
    else:
        rows = list(db.scalars(query.order_by(ADMETPrediction.created_at.desc())).all())

    latest = {}
    for row in rows:
        latest.setdefault(_model_key(row), row)
    values = {key: float(row.predicted_value) for key, row in latest.items() if row.predicted_value is not None}
    if not values:
        return None, {}, {}
    snapshot = dict(next(iter(latest.values())).outputs_json or {}).get("prediction_snapshot") or {}
    base = _number(snapshot.get("base_prediction"))
    if base is None:
        base = sum(values.values()) / len(values)
    project = _number(snapshot.get("project_prediction"))
    first = next(iter(latest.values()))
    return first, values, {"base_prediction": base, "project_prediction": project,
                           "snapshot": snapshot, "unit": first.unit}


def _pair_key(project_id, version_id, experiment_kind, experiment_id, endpoint, prediction_id):
    raw = f"{project_id}|{version_id}|{experiment_kind}|{experiment_id}|{endpoint}|{prediction_id}"
    return hashlib.sha256(raw.encode()).hexdigest()


def _write_pair(db: Session, *, project_id, version, endpoint_name, prediction, values,
                bundle, experiment_id, external_id, evidence_origin, experiment_at,
                experimental_value, experimental_unit, raw_value, raw_unit, relation,
                comparability_status, independent_group, evidence_quality="A", reason=""):
    pair_key = _pair_key(project_id, version.id, "external" if external_id else "internal",
                         external_id or experiment_id, endpoint_name, prediction.id if prediction else 0)
    if db.scalar(select(PredictionExperimentalPairRecord.id).where(PredictionExperimentalPairRecord.pair_key == pair_key)):
        return None
    prediction_at = _aware(prediction.created_at) if prediction else None
    experiment_at = _aware(experiment_at) or datetime.now(timezone.utc)
    prospective = bool(prediction_at and prediction_at < experiment_at)
    if prospective:
        pair_class = "TRUE_PROSPECTIVE"
    elif prediction:
        pair_class = "RETROSPECTIVE_OOF"
    else:
        pair_class = "HISTORICAL_VISIBLE"
    exp = _number(experimental_value)
    base = _number(bundle.get("base_prediction")) if bundle else None
    project = _number(bundle.get("project_prediction")) if bundle else None
    comparable = comparability_status in DIRECT_COMPARABILITY
    eligible = bool(comparable and exp is not None and base is not None and evidence_origin in {"INTERNAL_EXPERIMENTAL", "EXTERNAL_IMPORTED", "AUTO_QUALIFIED_EXTERNAL", "EXPERIMENTAL_EXTERNAL_AUTO", "EXPERIMENTAL_EXTERNAL"})
    
    if not reason:
        if not prediction:
            reason = "NO_PREEXPERIMENTAL_FREEZE"
        elif not comparable:
            reason = "CONTEXT_MISMATCH_OR_RELATED_ENDPOINT"
        elif exp is None:
            reason = "NON_NUMERIC_EXPERIMENT"
        elif eligible:
            reason = "Eligible prospective pair" if prospective else "Eligible retrospective OOF pair"
            
    return PredictionExperimentalPairRecord(
        pair_key=pair_key, project_id=project_id, compound_version_id=version.id,
        endpoint_name=endpoint_name, prediction_record_id=prediction.id if prediction else None,
        experiment_id=experiment_id, external_evidence_id=external_id,
        evidence_origin=evidence_origin, prediction_created_at=prediction_at,
        experiment_created_at=experiment_at, pair_class=pair_class,
        comparability_status=comparability_status or "NOT_COMPARABLE",
        adaptation_eligibility=eligible, independent_experiment_group_id=independent_group or pair_key,
        base_prediction=base, project_prediction=project, experimental_value=exp,
        experimental_unit=experimental_unit or "", raw_value=str(raw_value or ""), raw_unit=raw_unit or "",
        relation=relation or "=", absolute_error=abs(base - exp) if eligible else None,
        signed_error=base - exp if eligible else None,
        project_absolute_error=abs(project - exp) if eligible and project is not None else None,
        adapter_version=str((bundle or {}).get("snapshot", {}).get("adapter_version") or ""),
        included_in_future_adapter=eligible, exclusion_reason="" if eligible else reason,
        snapshot_json=(bundle or {}).get("snapshot", {}),
    )


def record_internal_measurement_pair(db: Session, project_id: int, version: CompoundVersion,
                                     measurement: ADMETMeasurement):
    endpoint_name = measurement.endpoint.name if measurement.endpoint else ""
    if not endpoint_name:
        return None
    prediction, values, bundle = _prediction_bundle(db, version.id, endpoint_name, measurement.created_at)
    status = "DIRECTLY_COMPARABLE" if prediction and _number(measurement.value) is not None and prediction.unit == measurement.unit else "CONDITIONALLY_COMPARABLE"
    row = _write_pair(
        db, project_id=project_id, version=version, endpoint_name=endpoint_name,
        prediction=prediction, values=values, bundle=bundle, experiment_id=measurement.id,
        external_id=None, evidence_origin="INTERNAL_EXPERIMENTAL", experiment_at=measurement.created_at,
        experimental_value=measurement.value if measurement.value is not None else measurement.mean_value,
        experimental_unit=measurement.unit, raw_value=measurement.value if measurement.value is not None else measurement.qualitative_value,
        raw_unit=measurement.unit, relation=measurement.qualifier, comparability_status=status,
        independent_group=(measurement.provenance_json or {}).get("independent_experiment_group_id") or f"MEAS-{measurement.id}",
        evidence_quality="A",
    )
    if row:
        db.add(row)
    return row


def record_external_evidence_pair(db: Session, project_id: int, evidence: ExternalExperimentalEvidence):
    version = db.get(CompoundVersion, evidence.compound_version_id)
    if not version:
        return None
    mapping = {
        "solubility_aqueous_logs": "Solubility", "SOLUBILITY_GENERIC": "Solubility",
        "SOLUBILITY_THERMODYNAMIC": "Solubility", "SOLUBILITY_KINETIC": "Solubility", "SOLUBILITY_INTRINSIC": "Solubility",
        "permeability_caco2_logpapp": "Permeability", "CACO2_PAPP_AB": "Permeability",
        "ppb_human_percent_bound": "Plasma protein binding", "HUMAN_PPB": "Plasma protein binding",
        "hlm_intrinsic_clearance_scaled_log10": "HLM intrinsic clearance", "HLM_CLINT": "HLM intrinsic clearance",
        "rlm_intrinsic_clearance_scaled_log10": "RLM intrinsic clearance", "RLM_CLINT": "RLM intrinsic clearance",
        "mlm_intrinsic_clearance_scaled_log10": "MLM intrinsic clearance", "MLM_CLINT": "MLM intrinsic clearance",
    }
    endpoint_name = mapping.get(evidence.canonical_endpoint_id, evidence.raw_endpoint_name)
    experiment_at = evidence.imported_at or evidence.retrieved_at or datetime.now(timezone.utc)
    prediction, values, bundle = _prediction_bundle(db, version.id, endpoint_name, experiment_at)
    row = _write_pair(
        db, project_id=project_id, version=version, endpoint_name=endpoint_name,
        prediction=prediction, values=values, bundle=bundle, experiment_id=None,
        external_id=evidence.id, evidence_origin="EXTERNAL_IMPORTED", experiment_at=experiment_at,
        experimental_value=evidence.normalized_value, experimental_unit=evidence.normalized_unit,
        raw_value=evidence.raw_value, raw_unit=evidence.raw_unit, relation=evidence.raw_relation,
        comparability_status=evidence.comparability_status,
        independent_group=evidence.source_document_id or evidence.source_record_id,
        evidence_quality=evidence.source_quality_class,
    )
    if row:
        db.add(row)
    return row


def record_canonical_evidence_pair(db: Session, project_id: int, evidence: ExternalExperimentalEvidence):
    """Pair shared internal/imported canonical evidence with a frozen model run."""
    version = db.get(CompoundVersion, evidence.compound_version_id)
    if not version:
        return None
    endpoint_name = {
        "SOLUBILITY_GENERIC": "Solubility", "SOLUBILITY_KINETIC": "Solubility",
        "SOLUBILITY_THERMODYNAMIC": "Solubility", "SOLUBILITY_INTRINSIC": "Solubility",
        "CACO2_PAPP_AB": "Permeability", "HUMAN_PPB": "Plasma protein binding",
        "HLM_CLINT": "HLM intrinsic clearance", "RLM_CLINT": "RLM intrinsic clearance",
        "MLM_CLINT": "MLM intrinsic clearance",
    }.get(evidence.canonical_endpoint_id, "")
    if not endpoint_name:
        return None
    experiment_at = evidence.imported_at or evidence.retrieved_at or datetime.now(timezone.utc)
    prediction, values, bundle = _prediction_bundle(db, version.id, endpoint_name, experiment_at)
    origin = "INTERNAL_EXPERIMENTAL" if evidence.evidence_state == "INTERNAL_EXPERIMENTAL" else ("AUTO_QUALIFIED_EXTERNAL" if evidence.evidence_state == "AUTO_QUALIFIED_EXTERNAL" else "EXTERNAL_IMPORTED")
    row = _write_pair(
        db, project_id=project_id, version=version, endpoint_name=endpoint_name,
        prediction=prediction, values=values, bundle=bundle, experiment_id=None,
        external_id=evidence.id, evidence_origin=origin, experiment_at=experiment_at,
        experimental_value=evidence.normalized_value, experimental_unit=evidence.normalized_unit,
        raw_value=evidence.raw_value, raw_unit=evidence.raw_unit, relation=evidence.raw_relation,
        comparability_status=evidence.comparability_status,
        independent_group=evidence.independent_experiment_group_id or evidence.source_record_id,
        evidence_quality=evidence.source_quality_class,
    )
    if row:
        db.add(row)
    return row


def ledger_out(row: PredictionExperimentalPairRecord):
    return {
        "id": row.id, "pair_key": row.pair_key, "project_id": row.project_id,
        "compound_version_id": row.compound_version_id, "endpoint": row.endpoint_name,
        "prediction_record_id": row.prediction_record_id, "experiment_id": row.experiment_id,
        "external_evidence_id": row.external_evidence_id, "evidence_origin": row.evidence_origin,
        "prediction_created_at": row.prediction_created_at.isoformat() if row.prediction_created_at else None,
        "experiment_created_at": row.experiment_created_at.isoformat() if row.experiment_created_at else None,
        "pair_class": row.pair_class, "comparability": row.comparability_status,
        "adaptation_eligibility": row.adaptation_eligibility,
        "independent_experiment_group_id": row.independent_experiment_group_id,
        "base_prediction": row.base_prediction, "project_prediction": row.project_prediction,
        "experimental_value": row.experimental_value, "experimental_unit": row.experimental_unit,
        "absolute_error": row.absolute_error, "signed_error": row.signed_error,
        "project_absolute_error": row.project_absolute_error,
        "adapter_version": row.adapter_version, "included_in_future_adapter": row.included_in_future_adapter,
        "exclusion_reason": row.exclusion_reason,
    }


def project_learning_summary(db: Session, project_id: int):
    rows = list(db.scalars(select(PredictionExperimentalPairRecord).where(
        PredictionExperimentalPairRecord.project_id == project_id
    ).order_by(PredictionExperimentalPairRecord.created_at, PredictionExperimentalPairRecord.id)).all())
    result = {}
    for row in rows:
        entry = result.setdefault(row.endpoint_name, {
            "endpoint": row.endpoint_name,
            "pairs": 0,
            "eligible_pairs": 0,
            "prospective_pairs": 0,
            "retrospective_oof_pairs": 0,
            "independent_compounds": set(),
            "effective_n": 0.0,
            "base_errors": [],
            "project_errors": [],
        })
        entry["pairs"] += 1
        if row.pair_class == "TRUE_PROSPECTIVE":
            entry["prospective_pairs"] += 1
        else:
            entry["retrospective_oof_pairs"] += 1

        if row.adaptation_eligibility:
            entry["eligible_pairs"] += 1
            entry["independent_compounds"].add(row.compound_version_id)
            if row.absolute_error is not None:
                entry["base_errors"].append(row.absolute_error)
            if row.project_absolute_error is not None:
                entry["project_errors"].append(row.project_absolute_error)

    for entry in result.values():
        entry["independent_compounds"] = len(entry.pop("independent_compounds"))
        entry["effective_n"] = float(entry["independent_compounds"])
        base_errors = entry.pop("base_errors")
        project_errors = entry.pop("project_errors")
        entry["base_mae"] = sum(base_errors) / len(base_errors) if base_errors else None
        entry["project_mae"] = sum(project_errors) / len(project_errors) if project_errors else None
        if entry["effective_n"] >= 5:
            entry["status"] = "ELIGIBLE_FOR_LIGHT_PROJECT_ADAPTATION"
            entry["reason"] = f"Validated on N={int(entry['effective_n'])} independent compounds"
        else:
            entry["status"] = "COLLECTING"
            entry["reason"] = f"INSUFFICIENT_INDEPENDENT_COMPOUNDS (N={int(entry['effective_n'])} < 5 required)"

    return list(result.values()), [ledger_out(row) for row in rows]
