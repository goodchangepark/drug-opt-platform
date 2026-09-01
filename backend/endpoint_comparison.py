"""Canonical, persisted endpoint aggregation for the compound detail UI."""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime

from sqlalchemy import select

from .admet import ADMETEndpoint, ADMETMeasurement, ADMETPrediction, ADMETModelRegistry, PredictionEndpointSnapshot
from .activity_models import ActivityMeasurement, ActivityPrediction, AssayDefinition
from .experimental_evidence_router import route_evidence
from .models import Compound, CompoundVersion, ExternalExperimentalEvidence


CANONICAL_ENDPOINTS = {
    "Solubility": "solubility_aqueous_logs",
    "Permeability": "permeability_caco2_logpapp",
    "Plasma protein binding": "ppb_human_percent_bound",
    "HLM intrinsic clearance": "hlm_intrinsic_clearance_scaled_log10",
    "RLM intrinsic clearance": "rlm_intrinsic_clearance_scaled_log10",
    "MLM intrinsic clearance": "mlm_intrinsic_clearance_scaled_log10",
}


def _iso(value):
    return value.isoformat() if value else None


def _number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _section(endpoint: str, canonical: str = ""):
    text = f"{endpoint} {canonical}".lower()
    if any(token in text for token in ("ic50", "ec50", " ki", " kd", "activity", "potency")):
        return "ACTIVITY"
    if any(token in text for token in ("cyp", "metabol", "hepatocyte", "microsom", "clint", "mass balance")):
        return "METABOLISM"
    if any(token in text for token in ("cmax", "tmax", "auc", "half-life", "clearance", "cl/f", "v d", "vd/f", "bioavailability", "excretion")):
        return "PK"
    if any(token in text for token in ("herg", "ames", "dili", "tox")):
        return "TOXICITY"
    if canonical or any(token in text for token in ("solubility", "permeability", "caco", "papp", "ppb", "protein", "pka", "logd", "transporter")):
        return "ADMET"
    return "UNCLASSIFIED"


def _reference(row):
    return {
        "source": row.source_database,
        "source_record_id": row.source_record_id,
        "document_id": row.source_document_id,
        "reference": row.reference_text,
        "url": row.source_url,
        "page": (row.assay_conditions_json or {}).get("page"),
        "section": (row.assay_conditions_json or {}).get("section"),
        "table": (row.assay_conditions_json or {}).get("table"),
    }


def _external(row):
    state = row.evidence_state or ("EXTERNAL_IMPORTED" if row.accepted_at else "EXTERNAL_CANDIDATE")
    comparable = row.comparability_status in {"DIRECTLY_COMPARABLE", "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION"}
    endpoint_id = row.canonical_endpoint_id or ""
    section = _section(row.raw_endpoint_name, endpoint_id)
    routing_reason = f"Matched {section} endpoint semantics" if section != "UNCLASSIFIED" else row.routing_reason
    return {
        "id": row.id, "origin": state, "state": state,
        "raw_value": row.raw_value, "raw_unit": row.raw_unit, "relation": row.raw_relation,
        "normalized_value": _number(row.normalized_value), "normalized_unit": row.normalized_unit,
        "species": row.species, "context": row.assay_conditions_json or {},
        "reference": _reference(row), "qualification": row.qualification_status,
        "comparability": row.comparability_status, "importable": comparable and state != "EXTERNAL_IMPORTED",
        "identity_match_status": row.identity_match_status, "reference_status": "REFERENCE_RESOLVED_IMPORTED" if state == "EXTERNAL_IMPORTED" else "REFERENCE_RESOLVED_CANDIDATE",
        "endpoint": row.raw_endpoint_name, "assay_type": row.assay_type, "assay_id": row.source_assay_id,
        "adaptation_eligibility": bool(row.accepted_at and comparable),
        "display_evidence_group_id": row.provenance_fingerprint or f"evidence-{row.id}",
        "independent_experiment_group_id": row.source_document_id or row.source_record_id or f"evidence-{row.id}",
        "canonical_endpoint_id": row.canonical_endpoint_id,
        "display_source": row.source_database,
        "routing_reason": routing_reason,
    }


def _comparison(prediction, experiments):
    if not prediction or not prediction.get("available"):
        return None
    for experiment in experiments:
        pv = prediction.get("display_value")
        ev = experiment.get("normalized_value")
        if pv is None or ev is None:
            continue
        status = experiment.get("comparability")
        if status not in {"DIRECTLY_COMPARABLE", "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION"}:
            continue
        diff = float(pv) - float(ev)
        return {
            "status": "DIRECT" if status == "DIRECTLY_COMPARABLE" else "CONVERTED",
            "comparability": status, "prediction_value": pv, "experimental_value": ev,
            "difference": diff, "signed_error": diff, "absolute_error": abs(diff),
            "preview": experiment.get("state") == "EXTERNAL_CANDIDATE",
        }
    return None


def build_endpoint_comparison(db, version_id: int) -> dict:
    version = db.get(CompoundVersion, version_id)
    if not version:
        raise ValueError("CompoundVersion not found")
    compound = db.get(Compound, version.compound_row_id)
    endpoint_names = {row.id: row.name for row in db.scalars(select(ADMETEndpoint).where(ADMETEndpoint.project_id == compound.project_id)).all()}
    endpoint_rows = {}

    predictions = db.scalars(
        select(ADMETPrediction).join(ADMETModelRegistry).where(
            ADMETPrediction.version_id == version_id,
            ADMETPrediction.execution_status == "SUCCESS",
        ).order_by(ADMETPrediction.created_at.desc())
    ).all()
    # PredictionEndpointSnapshot is the durable run index.  Prefer its latest
    # row for the default display so a second Predict creates a new visible
    # run without rewriting the older model output rows.
    latest_snapshots = {}
    for snapshot_row in db.scalars(
        select(PredictionEndpointSnapshot)
        .where(PredictionEndpointSnapshot.compound_version_id == version_id)
        .order_by(PredictionEndpointSnapshot.created_at.desc())
    ).all():
        latest_snapshots.setdefault(snapshot_row.endpoint_name, snapshot_row)
    # Keep the newest immutable prediction per model/endpoint, but join all
    # endpoint models into one row rather than rendering parallel lists.
    by_model = {}
    for row in predictions:
        key = (row.model.endpoint_name, row.model_id)
        by_model.setdefault(key, row)
    grouped = defaultdict(list)
    for row in by_model.values():
        grouped[row.model.endpoint_name].append(row)
    for endpoint, rows in grouped.items():
        first = rows[0]
        values = [float(row.predicted_value) for row in rows if row.predicted_value is not None]
        if not values:
            continue
        snapshot = dict(first.outputs_json or {}).get("prediction_snapshot") or {}
        persisted_snapshot = latest_snapshots.get(endpoint)
        if persisted_snapshot is not None:
            snapshot = dict(persisted_snapshot.snapshot_json or snapshot)
        base_value = (persisted_snapshot.base_value if persisted_snapshot is not None and persisted_snapshot.base_value is not None
                      else snapshot.get("base_prediction", sum(values) / len(values)))
        project_value = (persisted_snapshot.project_value if persisted_snapshot is not None
                         else snapshot.get("project_prediction"))
        display_value = project_value if project_value is not None else base_value
        endpoint_id = CANONICAL_ENDPOINTS.get(endpoint, endpoint.lower().replace(" ", "_"))
        endpoint_rows[endpoint_id] = {
            "endpoint_id": endpoint_id, "section": _section(endpoint, endpoint_id), "display_name": endpoint,
            "prediction": {
                "available": True, "prediction_snapshot_id": persisted_snapshot.id if persisted_snapshot is not None else None, "prediction_run_id": persisted_snapshot.prediction_run_id if persisted_snapshot is not None else first.run_id, "base_value": base_value,
                "project_value": project_value, "display_value": display_value,
                "unit": persisted_snapshot.base_unit if persisted_snapshot is not None else first.unit, "prediction_type": (snapshot.get("prediction_type") or "BASE_PREDICTION"),
                "adapter": snapshot.get("adapter_version") or "", "maturity": snapshot.get("maturity") or {
                    "level": 1, "label": "Base Prediction", "stars": "★☆☆☆☆"
                }, "ood": first.applicability_domain, "timestamp": _iso(first.created_at),
                "model_count": len(rows), "model_predictions": {str(row.model_id): row.predicted_value for row in rows},
            },
            "experimental_internal": [], "experimental_external_imported": [],
            "experimental_external_candidates": [], "related_evidence": [], "needs_review": [],
            "references": [], "project_learning": {},
        }

    activity_assays = {row.id: row for row in db.scalars(select(AssayDefinition).where(AssayDefinition.project_id == compound.project_id)).all()}
    activity_predictions = db.scalars(select(ActivityPrediction).where(ActivityPrediction.version_id == version_id).order_by(ActivityPrediction.created_at.desc())).all()
    activity_measurements = db.scalars(select(ActivityMeasurement).where(ActivityMeasurement.version_id == version_id).order_by(ActivityMeasurement.created_at.desc())).all()
    for assay_id, assay in activity_assays.items():
        pred = next((item for item in activity_predictions if item.assay_id == assay_id), None)
        measured = [item for item in activity_measurements if item.assay_id == assay_id]
        if pred is None and not measured:
            continue
        endpoint_id = f"activity:{assay.measurement_type.lower()}:{assay_id}"
        row = endpoint_rows.setdefault(endpoint_id, {
            "endpoint_id": endpoint_id, "section": "ACTIVITY", "display_name": f"{assay.name} ({assay.measurement_type})",
            "prediction": {"available": False}, "experimental_internal": [], "experimental_external_imported": [],
            "experimental_external_candidates": [], "related_evidence": [], "needs_review": [], "references": [], "project_learning": {},
        })
        if pred is not None:
            row["prediction"] = {"available": True, "prediction_run_id": None, "base_value": pred.predicted_value_nm,
                "project_value": None, "display_value": pred.predicted_value_nm, "unit": pred.assay.unit,
                "prediction_type": pred.prediction_type, "adapter": "", "maturity": {"level": 1, "label": "Base Prediction", "stars": "★☆☆☆☆"},
                "ood": pred.applicability_domain, "timestamp": _iso(pred.created_at), "model_count": 1}
        for measurement in measured:
            row["experimental_internal"].append({"id": measurement.id, "origin": "INTERNAL_EXPERIMENTAL", "state": "INTERNAL_EXPERIMENTAL",
                "raw_value": measurement.raw_value, "normalized_value": measurement.normalized_value_nm, "raw_unit": measurement.original_unit,
                "normalized_unit": "nM", "relation": measurement.qualifier, "species": assay.species,
                "context": {"target": assay.target, "cell_line": assay.cell_line, "assay": assay.name},
                "reference": {"source": measurement.source, "reference": measurement.notes}, "qualification": "QUALIFIED_DIRECT",
                "comparability": "DIRECTLY_COMPARABLE", "importable": False, "adaptation_eligibility": True})

    measurements = db.scalars(select(ADMETMeasurement).where(ADMETMeasurement.version_id == version_id)).all()
    for measurement in measurements:
        endpoint = endpoint_names.get(measurement.endpoint_id, "Unknown endpoint")
        endpoint_id = CANONICAL_ENDPOINTS.get(endpoint, endpoint.lower().replace(" ", "_"))
        row = endpoint_rows.setdefault(endpoint_id, {
            "endpoint_id": endpoint_id, "section": _section(endpoint, endpoint_id), "display_name": endpoint,
            "prediction": {"available": False}, "experimental_internal": [], "experimental_external_imported": [],
            "experimental_external_candidates": [], "related_evidence": [], "needs_review": [], "references": [], "project_learning": {},
        })
        row["experimental_internal"].append({
            "id": measurement.id, "origin": "INTERNAL_EXPERIMENTAL", "state": "INTERNAL_EXPERIMENTAL",
            "raw_value": measurement.value if measurement.value is not None else measurement.qualitative_value,
            "normalized_value": measurement.value, "raw_unit": measurement.unit, "normalized_unit": measurement.unit,
            "species": measurement.species, "context": {"matrix": measurement.matrix, "method": measurement.method},
            "reference": {"source": measurement.source, "reference": measurement.notes},
            "qualification": "QUALIFIED_DIRECT", "comparability": "DIRECTLY_COMPARABLE", "importable": False,
            "adaptation_eligibility": True,
        })

    evidence_rows = db.scalars(select(ExternalExperimentalEvidence).where(ExternalExperimentalEvidence.compound_version_id == version_id)).all()
    for evidence in evidence_rows:
        endpoint = evidence.canonical_endpoint_id or evidence.raw_endpoint_name
        endpoint_id = endpoint if endpoint in CANONICAL_ENDPOINTS.values() else CANONICAL_ENDPOINTS.get(evidence.raw_endpoint_name, endpoint.lower().replace(" ", "_"))
        # Re-derive the section from the persisted canonical endpoint.  This
        # keeps records harvested by an older router display-safe after a
        # routing-version upgrade (without deleting or rewriting raw values).
        section = _section(evidence.raw_endpoint_name, endpoint_id)
        row = endpoint_rows.setdefault(endpoint_id, {
            "endpoint_id": endpoint_id, "section": section, "display_name": evidence.raw_endpoint_name,
            "prediction": {"available": False}, "experimental_internal": [], "experimental_external_imported": [],
            "experimental_external_candidates": [], "related_evidence": [], "needs_review": [], "references": [], "project_learning": {},
        })
        item = _external(evidence)
        if evidence.evidence_state == "EXTERNAL_IMPORTED":
            row["experimental_external_imported"].append(item)
        elif evidence.comparability_status == "RELATED_NOT_SAME_ENDPOINT":
            row["related_evidence"].append(item)
        elif evidence.qualification_status in {"NEEDS_REVIEW", "UNSUPPORTED"}:
            row["needs_review"].append(item)
        else:
            row["experimental_external_candidates"].append(item)
        row["references"].append(item["reference"])

    for row in endpoint_rows.values():
        experiments = row["experimental_internal"] + row["experimental_external_imported"] + row["experimental_external_candidates"]
        row["comparison"] = _comparison(row["prediction"], experiments)
        row["summary"] = {
            "both": int(bool(row["prediction"].get("available") and experiments)),
            "prediction_only": int(bool(row["prediction"].get("available") and not experiments)),
            "experimental_only": int(bool(experiments and not row["prediction"].get("available"))),
            "related": len(row["related_evidence"]), "needs_review": len(row["needs_review"]),
            "ready_to_import": sum(bool(item.get("importable")) for item in row["experimental_external_candidates"]),
        }
    endpoints = sorted(endpoint_rows.values(), key=lambda item: (item["section"], item["display_name"]))
    summary = {key: sum(row["summary"][key] for row in endpoints) for key in ("both", "prediction_only", "experimental_only", "related", "needs_review", "ready_to_import")}
    summary["imported_pairs"] = 0
    return {"version_id": version_id, "project_id": compound.project_id, "compound_id": compound.id, "endpoints": endpoints, "summary": summary}
