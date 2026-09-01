"""Canonical endpoint aggregation for the persisted comparison API.

Raw source labels are audit fields, never comparison keys. Both prediction
snapshots and experimental evidence pass through the versioned registry.
"""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
import hashlib
from sqlalchemy import select

from .admet import ADMETEndpoint, ADMETMeasurement, ADMETPrediction, ADMETModelRegistry, ADMETPredictionRun, PredictionEndpointSnapshot
from .activity_models import ActivityMeasurement, ActivityPrediction, AssayDefinition
from .canonical_endpoints import (
    CANONICAL_ENDPOINT_VERSION, COMPARISON_UNIT_VERSION, CONVERTED, DIRECT,
    RELATED, REGISTRY, UNSUPPORTED, canonicalize_prediction_endpoint,
    endpoint_contract, normalize_experimental_observation, normalize_species,
    prediction_source_label, prediction_source_type, PREDICTION_DERIVED,
    PREDICTION_MECHANISTIC, PREDICTION_MODEL, PREDICTION_RULE, PREDICTION_UNAVAILABLE,
)
from .ivive import PKParameterSet
from .pk import PKNCAResult, PKStudy
from .simulation import PKSimulationRun
from .metabolism import MetabolicPredictionRun
from .qualification_contract import aggregate_qualification, qualify_record
from .models import Compound, CompoundVersion, ExternalExperimentalEvidence


CANONICAL_ENDPOINTS = {
    "Solubility": "SOLUBILITY_GENERIC", "Permeability": "CACO2_PAPP_AB",
    "Plasma protein binding": "HUMAN_PPB", "HLM intrinsic clearance": "HLM_CLINT",
    "RLM intrinsic clearance": "RLM_CLINT", "MLM intrinsic clearance": "MLM_CLINT",
}


def _iso(value):
    return value.isoformat() if value else None


def _number(value):
    try: return float(value)
    except (TypeError, ValueError): return None


def _context(row):
    value = getattr(row, "assay_conditions_json", None)
    return value if isinstance(value, dict) else {"conditions": value or ""}


def _dose_mg_kg(value, unit=""):
    number = _number(value)
    if number is None: return None
    text = str(unit or "mg/kg").lower().replace("µ", "u")
    if "ug/kg" in text or "mcg/kg" in text: return number / 1000.0
    if "g/kg" in text and "mg/kg" not in text: return number * 1000.0
    return number


def _section(endpoint_id):
    contract = endpoint_contract(endpoint_id)
    if contract: return contract.section
    text = str(endpoint_id or "").upper()
    if text.startswith("ACTIVITY_"): return "ACTIVITY"
    if text.startswith("PK_"): return "PK"
    if text.startswith(("CYP", "METABOLITE", "HEPATOCYTE", "EXCRETION")): return "METABOLISM"
    if text.startswith(("HERG", "AMES", "DILI")): return "TOXICITY"
    if text.startswith(("SOLUBILITY", "CACO2", "HUMAN_PPB", "RAT_PPB", "MOUSE_PPB", "HLM", "RLM", "MLM", "PKA", "LOG")): return "ADMET"
    return "UNCLASSIFIED"


def _display_name(endpoint_id, fallback=""):
    contract = endpoint_contract(endpoint_id)
    if contract: return contract.display_name
    if str(endpoint_id).startswith("ACTIVITY_"): return str(endpoint_id).removeprefix("ACTIVITY_").split(":", 1)[0]
    return fallback or str(endpoint_id).replace("_", " ").title()


def _reference(row):
    context = _context(row)
    return {"source": row.source_database, "source_record_id": row.source_record_id, "document_id": row.source_document_id, "reference": row.reference_text, "url": row.source_url, "page": context.get("page"), "section": context.get("section"), "table": context.get("table")}


def _mapped_external(row):
    context = _context(row)
    mapped = normalize_experimental_observation(row.raw_endpoint_name, row.raw_value, row.raw_unit, species=row.species, context=context, assay_type=row.assay_type, target=context.get("target", ""), canonical_hint=row.canonical_endpoint_id)
    state = row.evidence_state or ("EXTERNAL_IMPORTED" if row.accepted_at else "EXTERNAL_CANDIDATE")
    endpoint_id = mapped["canonical_endpoint_id"]
    comparable = mapped["comparability_status"] in {DIRECT, CONVERTED}
    qualification = {DIRECT: "QUALIFIED_DIRECT", CONVERTED: "QUALIFIED_DETERMINISTIC_CONVERSION", RELATED: "QUALIFIED_RELATED", "CONDITIONALLY_COMPARABLE": "QUALIFIED_CONDITIONAL"}.get(mapped["comparability_status"], "NEEDS_REVIEW")
    item = {"id": row.id, "origin": state, "state": state, "raw_endpoint": row.raw_endpoint_name, "endpoint": _display_name(endpoint_id, row.raw_endpoint_name), "raw_value": row.raw_value, "raw_unit": row.raw_unit, "relation": row.raw_relation, "normalized_value": mapped.get("normalized_value"), "normalized_unit": mapped.get("normalized_unit", ""), "species": mapped.get("species", normalize_species(row.species, context)), "route": mapped.get("route", "UNSPECIFIED"), "context": context, "dose": context.get("dose"), "dose_unit": context.get("dose_unit") or context.get("dose_units"), "reference": _reference(row), "qualification": qualification, "comparability": mapped["comparability_status"], "importable": comparable and state != "EXTERNAL_IMPORTED", "identity_match_status": row.identity_match_status, "reference_status": "REFERENCE_RESOLVED_IMPORTED" if state == "EXTERNAL_IMPORTED" else "REFERENCE_RESOLVED_CANDIDATE", "assay_type": row.assay_type, "assay_id": row.source_assay_id, "adaptation_eligibility": bool(row.accepted_at and comparable), "display_evidence_group_id": row.display_evidence_group_id or row.provenance_fingerprint or f"evidence-{row.id}", "independent_experiment_group_id": row.independent_experiment_group_id or row.source_document_id or row.source_record_id or f"evidence-{row.id}", "canonical_endpoint_id": endpoint_id, "canonical_comparison_key": mapped["comparison_key"], "display_source": row.source_database, "routing_reason": mapped.get("reason", "") or row.routing_reason, "normalization_rule": mapped.get("normalization_rule", ""), "raw_persisted_canonical_endpoint_id": row.canonical_endpoint_id, "qualification_details": row.qualification_json or {}}
    return item, mapped


def _prediction_object(snapshot, rows, endpoint_id, raw_endpoint, species="", route=""):
    first = rows[0]
    snap = dict(snapshot.snapshot_json or {}) if snapshot else dict(first.outputs_json or {}).get("prediction_snapshot") or {}
    base = snapshot.base_value if snapshot and snapshot.base_value is not None else snap.get("base_prediction", first.predicted_value)
    project = snapshot.project_value if snapshot else snap.get("project_prediction")
    display = project if project is not None else base
    unit = snapshot.base_unit if snapshot and snapshot.base_unit else first.unit
    source_type = snap.get("source_type") or prediction_source_type(
        source=first.model.model_name, prediction_type=snap.get("prediction_type"), endpoint=raw_endpoint
    )
    return {"available": True, "prediction_snapshot_id": snapshot.id if snapshot else None, "prediction_run_id": snapshot.prediction_run_id if snapshot else first.run_id, "raw_endpoint": raw_endpoint, "canonical_endpoint_id": endpoint_id, "canonical_comparison_key": f"{endpoint_id}|{species}|{route}", "base_value": base, "project_value": project, "display_value": display, "unit": unit, "prediction_type": snap.get("prediction_type", "MODEL"), "source_type": source_type, "source_label": prediction_source_label(source_type), "adapter": snap.get("adapter_version", ""), "maturity": snap.get("maturity") or {"level": 1, "label": "Base Prediction", "stars": "★☆☆☆☆"}, "ood": first.applicability_domain, "timestamp": _iso(snapshot.created_at if snapshot else first.created_at), "model_count": len(rows), "model_predictions": {str(row.model_id): row.predicted_value for row in rows}, "species": species, "route": route}


def _snapshot_prediction(snapshot, endpoint_id, raw_endpoint, *, species="", route="", dose=None, dose_unit=""):
    """Build the same public prediction contract for non-ADMET snapshots."""
    data = dict(snapshot.snapshot_json or {})
    source_type = data.get("source_type") or prediction_source_type(
        source=data.get("source", ""), prediction_type=snapshot.prediction_type, endpoint=raw_endpoint
    )
    base = snapshot.base_value
    project = snapshot.project_value
    return {
        "available": base is not None or project is not None,
        "prediction_snapshot_id": snapshot.id,
        "prediction_run_id": snapshot.prediction_run_id,
        "raw_endpoint": raw_endpoint,
        "canonical_endpoint_id": endpoint_id,
        "canonical_comparison_key": f"{endpoint_id}|{species}|{route}|PARENT",
        "base_value": base, "project_value": project,
        "display_value": project if project is not None else base,
        "unit": snapshot.project_unit or snapshot.base_unit,
        "prediction_type": snapshot.prediction_type,
        "source_type": source_type,
        "source_label": prediction_source_label(source_type),
        "adapter": snapshot.adapter_version or "",
        "maturity": data.get("maturity") or {"level": snapshot.maturity_level or 1, "label": snapshot.maturity_label or "Base Prediction", "stars": "★☆☆☆☆"},
        "ood": data.get("ood_applicability", "UNKNOWN"),
        "timestamp": _iso(snapshot.created_at),
        "model_count": data.get("model_count", 0),
        "model_predictions": data.get("model_predictions", {}),
        "species": species, "route": route,
        "dose": dose if dose is not None else data.get("dose"),
        "dose_unit": dose_unit or data.get("dose_unit", ""),
        "provenance": data.get("provenance", {}),
    }


def _pk_prediction_source(source: str, *, simulation: bool = False) -> str:
    if simulation:
        return PREDICTION_MECHANISTIC
    return prediction_source_type(source=source, endpoint="PK foundation", default=PREDICTION_DERIVED)


def _pk_snapshot_values(pset):
    """Yield only actual non-experimental PK foundation outputs."""
    route = "ORAL" if str(pset.route).upper() == "PO" else str(pset.route).upper()
    v_type = str(pset.v_type or "").upper()
    values = [
        ("CLF_ORAL" if route == "ORAL" else "CL", pset.cl_value, pset.cl_unit, pset.cl_source_type),
        ("VDF_ORAL" if route == "ORAL" else ("VSS" if "VSS" in v_type else "VD"), pset.v_value, pset.v_unit, pset.v_source_type),
        ("F", pset.f_predicted, "%" if pset.f_predicted is not None and abs(float(pset.f_predicted)) > 1 else "fraction", (pset.provenance_json or {}).get("f_source_type", "MECHANISTIC_ASSEMBLY")),
    ]
    for parameter, value, unit, source in values:
        if value is None or prediction_source_type(source=source, endpoint=parameter) == PREDICTION_UNAVAILABLE:
            continue
        yield parameter, value, unit, _pk_prediction_source(source)


def _simulation_values(sim):
    metrics = sim.output_metrics or {}
    return [
        ("CMAX", metrics.get("cmax_ng_ml"), "ng/mL"),
        ("TMAX", metrics.get("tmax_hours"), "hours"),
        ("AUC0_T", metrics.get("auc_last_ng_h_ml"), "ng*h/mL"),
        ("AUC0_INF", metrics.get("auc_inf_analytical_ng_h_ml"), "ng*h/mL"),
        ("T_HALF", metrics.get("half_life_hours"), "hours"),
    ]


def persist_pk_prediction_snapshots(db, version_id: int, prediction_run_id: int | None = None, *, reuse_existing: bool = False) -> dict:
    """Index persisted IVIVE/Stage-5 outputs as immutable endpoint snapshots.

    This does not execute a model or alter an existing PK result.  It only
    makes already persisted calculations available to the canonical endpoint
    API and gives them a durable prediction-run identity.
    """
    version = db.get(CompoundVersion, version_id)
    if not version:
        raise ValueError("CompoundVersion not found")
    psets = db.scalars(select(PKParameterSet).where(PKParameterSet.version_id == version_id)).all()
    sims = db.scalars(select(PKSimulationRun).where(PKSimulationRun.version_id == version_id)).all()
    metabolism_runs = db.scalars(select(MetabolicPredictionRun).where(MetabolicPredictionRun.version_id == version_id, MetabolicPredictionRun.status == "COMPLETE")).all()
    if not psets and not sims and not metabolism_runs:
        return {"prediction_run_id": prediction_run_id, "created": 0, "existing": 0}
    run = db.get(ADMETPredictionRun, prediction_run_id) if prediction_run_id else None
    if run is None and reuse_existing:
        run = db.scalar(select(ADMETPredictionRun).where(
            ADMETPredictionRun.version_id == version_id,
            ADMETPredictionRun.requested_by == "system:canonical-prediction-index",
        ).order_by(ADMETPredictionRun.started_at.desc()))
    if run is None:
        digest = hashlib.sha256(f"canonical-pk|{version_id}|{datetime.now(timezone.utc).isoformat()}".encode()).hexdigest()
        run = ADMETPredictionRun(version_id=version_id, requested_by="system:canonical-prediction-index", inputs_hash=digest, status="COMPLETE", message="Persisted IVIVE and Stage-5 PK calculation index.", started_at=datetime.now(timezone.utc), completed_at=datetime.now(timezone.utc))
        db.add(run); db.flush()
    existing = {(row.prediction_run_id, row.endpoint_id) for row in db.scalars(select(PredictionEndpointSnapshot).where(PredictionEndpointSnapshot.prediction_run_id == run.id)).all()}
    created = 0
    for pset in psets:
        species = normalize_species(pset.species); route = "ORAL" if str(pset.route).upper() == "PO" else str(pset.route).upper()
        for parameter, value, unit, source_type in _pk_snapshot_values(pset):
            # ``f_predicted`` in the IVIVE foundation is oral bioavailability
            # relative to the IV reference arm.  It is not an IV
            # bioavailability endpoint, even when the foundation row is
            # stored under its IV route context.
            eid = f"{species}_PK_F_ORAL" if parameter == "F" else f"{species}_PK_{parameter}_{route}"
            if (run.id, eid) in existing: continue
            if parameter == "F":
                normalized_value, normalized_unit = value, "%"
                snapshot_route = "ORAL"
            else:
                mapped = normalize_experimental_observation(parameter, value, unit, species=species, context={"route": route, "dose": pset.dose_value, "dose_unit": pset.dose_unit})
                normalized_value, normalized_unit, snapshot_route = mapped.get("normalized_value", value), mapped.get("normalized_unit", unit), route
            snapshot = {"source": "IVIVE PK foundation", "source_type": source_type, "prediction_type": source_type, "species": species, "route": snapshot_route, "reference_route": "IV" if parameter == "F" else "", "dose": pset.dose_value, "dose_unit": pset.dose_unit, "v_type": pset.v_type, "provenance": pset.provenance_json or {}, "pk_parameter_set_id": pset.id, "canonical_endpoint_version": CANONICAL_ENDPOINT_VERSION, "comparison_unit_version": COMPARISON_UNIT_VERSION, "maturity": {"level": 1, "label": "Base Prediction", "stars": "★☆☆☆☆"}}
            db.add(PredictionEndpointSnapshot(prediction_run_id=run.id, project_id=version.compound.project_id, compound_version_id=version_id, endpoint_id=eid, endpoint_name=eid, base_value=normalized_value, base_unit=normalized_unit, prediction_type=source_type, maturity_level=1, maturity_label="Base Prediction", snapshot_json=snapshot, created_at=pset.created_at))
            existing.add((run.id, eid)); created += 1
    for sim in sims:
        species = normalize_species(sim.species); route = "ORAL" if str(sim.route).upper() == "PO" else str(sim.route).upper()
        for parameter, value, unit in _simulation_values(sim):
            if value is None: continue
            eid = f"{species}_PK_{parameter}_{route}"
            if (run.id, eid) in existing: continue
            snapshot = {"source": "Stage-5 PK simulation", "source_type": PREDICTION_MECHANISTIC, "prediction_type": PREDICTION_MECHANISTIC, "species": species, "route": route, "dose": sim.dose, "dose_unit": sim.dose_unit, "simulation_run_id": sim.id, "provenance": sim.provenance or {}, "canonical_endpoint_version": CANONICAL_ENDPOINT_VERSION, "comparison_unit_version": COMPARISON_UNIT_VERSION, "maturity": {"level": 1, "label": "Base Prediction", "stars": "★☆☆☆☆"}}
            db.add(PredictionEndpointSnapshot(prediction_run_id=run.id, project_id=version.compound.project_id, compound_version_id=version_id, endpoint_id=eid, endpoint_name=eid, base_value=value, base_unit=unit, prediction_type=PREDICTION_MECHANISTIC, maturity_level=1, maturity_label="Base Prediction", snapshot_json=snapshot, created_at=sim.created_at))
            existing.add((run.id, eid)); created += 1
    # SyGMa/SMARTCyp output is a ranked rule/derived hypothesis set, not a
    # scalar ML endpoint.  Persist its counts so the canonical UI cannot lose
    # the output between the metabolism page and a comparison reload.
    if metabolism_runs:
        met_run = sorted(metabolism_runs, key=lambda item: item.started_at or datetime.min.replace(tzinfo=timezone.utc), reverse=True)[0]
        for eid, value, unit, label in (
            ("METABOLIC_SOFT_SPOTS", len(met_run.spots), "ranked sites", "Metabolic soft spots"),
            ("METABOLITE_HYPOTHESES", len(met_run.metabolites), "hypotheses", "Metabolite hypotheses"),
        ):
            if (run.id, eid) in existing or value == 0: continue
            snapshot = {"source": "SyGMa/SMARTCyp metabolism calculation", "source_type": PREDICTION_RULE, "prediction_type": PREDICTION_RULE, "metabolic_prediction_run_id": met_run.id, "provenance": {"engine_name": met_run.engine_name, "engine_version": met_run.engine_version}, "canonical_endpoint_version": CANONICAL_ENDPOINT_VERSION, "comparison_unit_version": COMPARISON_UNIT_VERSION, "maturity": {"level": 1, "label": "Base Prediction", "stars": "★☆☆☆☆"}}
            db.add(PredictionEndpointSnapshot(prediction_run_id=run.id, project_id=version.compound.project_id, compound_version_id=version_id, endpoint_id=eid, endpoint_name=eid, base_value=float(value), base_unit=unit, prediction_type=PREDICTION_RULE, maturity_level=1, maturity_label="Base Prediction", snapshot_json=snapshot, created_at=met_run.completed_at or met_run.started_at))
            existing.add((run.id, eid)); created += 1
    db.flush()
    return {"prediction_run_id": run.id, "created": created, "existing": len(existing) - created}


def ensure_pk_prediction_snapshot_index(db) -> dict:
    """Backfill the durable index for old, already persisted calculations."""
    versions = db.scalars(select(CompoundVersion.id)).all()
    total = 0
    runs = 0
    for version_id in versions:
        result = persist_pk_prediction_snapshots(db, version_id, reuse_existing=True)
        total += result["created"]
        runs += int(result["prediction_run_id"] is not None and result["created"] > 0)
    return {"versions_examined": len(versions), "snapshots_created": total, "runs_created": runs}


def ensure_admet_prediction_snapshot_index(db, version_id: int | None = None) -> dict:
    """Index every successful user-facing ADMET prediction exactly once.

    Older runs froze only the production-core rows while the workspace could
    still expose additional successful model rows.  This backfill creates the
    missing immutable per-run/per-endpoint index without recomputing a model
    or changing a prediction value.
    """
    query = select(ADMETPrediction).where(ADMETPrediction.execution_status == "SUCCESS")
    if version_id is not None:
        query = query.where(ADMETPrediction.version_id == version_id)
    predictions = db.scalars(query.order_by(ADMETPrediction.created_at.asc())).all()
    # SQLite stores the snapshot endpoint key as text while the legacy ADMET
    # endpoint foreign key is integer-valued.  Normalize both sides before
    # checking the immutable uniqueness constraint.
    existing = {(row.prediction_run_id, str(row.endpoint_id)) for row in db.scalars(select(PredictionEndpointSnapshot)).all()}
    # Prediction workflows add the primary snapshot to the SQLAlchemy unit of
    # work immediately before this compatibility index runs.  Include pending
    # objects as well as committed rows; otherwise the index attempts to
    # insert the same (run, endpoint) pair a second time during one request.
    existing.update(
        (row.prediction_run_id, str(row.endpoint_id))
        for row in db.new
        if isinstance(row, PredictionEndpointSnapshot)
    )
    created = 0
    grouped = defaultdict(list)
    for prediction in predictions:
        grouped[(prediction.run_id, prediction.endpoint_id)].append(prediction)
    for (run_id, endpoint_id), rows in grouped.items():
        if (run_id, str(endpoint_id)) in existing:
            continue
        first = rows[0]
        raw_endpoint = first.model.endpoint_name
        values = {str(row.model_id): float(row.predicted_value) for row in rows if row.predicted_value is not None}
        if not values:
            continue
        base = sum(values.values()) / len(values)
        source_type = prediction_source_type(source=first.model.model_name, endpoint=raw_endpoint, default=PREDICTION_MODEL)
        snapshot = dict(first.outputs_json or {}).get("prediction_snapshot") or {
            "compound_version_id": first.version_id,
            "project_id": first.version.compound.project_id if first.version and first.version.compound else None,
            "endpoint": raw_endpoint, "base_prediction": base,
            "project_prediction": None, "project_adjustment": 0.0,
            "model_predictions": values, "source_type": source_type,
            "prediction_type": source_type, "canonical_endpoint_version": CANONICAL_ENDPOINT_VERSION,
            "comparison_unit_version": COMPARISON_UNIT_VERSION,
            "maturity": {"level": 1, "label": "Base Prediction", "stars": "★☆☆☆☆"},
        }
        project_id = first.version.compound.project_id if first.version and first.version.compound else 0
        db.add(PredictionEndpointSnapshot(
            prediction_run_id=run_id, project_id=project_id, compound_version_id=first.version_id,
            endpoint_id=str(endpoint_id), endpoint_name=raw_endpoint, base_value=base,
            base_unit=first.unit or "", prediction_type=source_type,
            maturity_level=(snapshot.get("maturity") or {}).get("level", 1),
            maturity_label=(snapshot.get("maturity") or {}).get("label", "Base Prediction"),
            snapshot_json=snapshot, created_at=first.created_at,
        ))
        existing.add((run_id, str(endpoint_id))); created += 1
    if created:
        db.flush()
    return {"predictions_examined": len(predictions), "snapshots_created": created}


def requalify_persisted_evidence(db, version_id: int | None = None) -> dict:
    """Recompute v4 qualification from stored raw evidence, without search."""
    versions = [version_id] if version_id is not None else db.scalars(select(CompoundVersion.id)).all()
    changed = 0
    for current_version_id in versions:
        prediction_ids = {row.endpoint_id for row in db.scalars(select(PredictionEndpointSnapshot).where(PredictionEndpointSnapshot.compound_version_id == current_version_id)).all()}
        for evidence in db.scalars(select(ExternalExperimentalEvidence).where(ExternalExperimentalEvidence.compound_version_id == current_version_id)).all():
            item, mapped = _mapped_external(evidence)
            q = qualify_record(item, prediction_endpoints=prediction_ids, imported=evidence.evidence_state == "EXTERNAL_IMPORTED")
            evidence.canonical_endpoint_id = mapped.get("canonical_endpoint_id", evidence.canonical_endpoint_id)
            evidence.normalized_value = "" if mapped.get("normalized_value") is None else str(mapped.get("normalized_value"))
            evidence.normalized_unit = str(mapped.get("normalized_unit") or evidence.normalized_unit or "")
            evidence.normalization_rule = str(mapped.get("normalization_rule") or evidence.normalization_rule or "")
            evidence.comparability_status = str(mapped.get("comparability_status") or evidence.comparability_status or UNSUPPORTED)
            evidence.qualification_version = q["qualification_version"]
            evidence.qualification_status = q["endpoint_status"]
            evidence.routing_section = mapped.get("section", evidence.routing_section or "")
            evidence.routing_reason = q.get("primary_gap_reason", "")
            evidence.qualification_json = q
            evidence.display_evidence_group_id = evidence.display_evidence_group_id or evidence.provenance_fingerprint or f"evidence-{evidence.id}"
            evidence.independent_experiment_group_id = evidence.independent_experiment_group_id or evidence.source_document_id or evidence.source_record_id or f"evidence-{evidence.id}"
            evidence.canonical_endpoint_version = CANONICAL_ENDPOINT_VERSION
            evidence.unit_normalization_version = COMPARISON_UNIT_VERSION
            changed += 1
    db.flush()
    return {"versions_examined": len(versions), "evidence_requalified": changed}


def _blank(endpoint_id, display_name=""):
    return {"endpoint_id": endpoint_id, "canonical_comparison_key": endpoint_id, "section": _section(endpoint_id), "display_name": _display_name(endpoint_id, display_name), "species": "UNSPECIFIED", "route": "UNSPECIFIED", "prediction": {"available": False}, "experimental_internal": [], "experimental_external_imported": [], "experimental_external_candidates": [], "related_evidence": [], "needs_review": [], "references": [], "project_learning": {}}


def _comparison(prediction, experiments):
    if not prediction or not prediction.get("available"):
        if any(item.get("comparability") == RELATED for item in experiments): return {"status": "RELATED_SAME_SCIENTIFIC_GROUP", "reason": "Related measurement semantics; no numeric error calculated", "difference": None}
        return None
    direct, related = [], []
    for experiment in experiments:
        if experiment.get("comparability") in {DIRECT, CONVERTED} and experiment.get("normalized_value") is not None:
            # PK rows include analyte in their comparison key. A parent
            # prediction must never be scored against a metabolite result.
            if str(prediction.get("canonical_endpoint_id", "")).startswith("PK_") or "_PK_" in str(prediction.get("canonical_endpoint_id", "")):
                pred_key = str(prediction.get("canonical_comparison_key", ""))
                exp_key = str(experiment.get("canonical_comparison_key", ""))
                if pred_key and exp_key and pred_key.rsplit("|", 1)[-1] != exp_key.rsplit("|", 1)[-1]:
                    continue
                pred_route, exp_route = prediction.get("route", "UNSPECIFIED"), experiment.get("route", "UNSPECIFIED")
                if pred_route != "UNSPECIFIED" and exp_route != "UNSPECIFIED" and pred_route != exp_route:
                    continue
                pred_dose = _dose_mg_kg(prediction.get("dose"), prediction.get("dose_unit"))
                exp_dose = _dose_mg_kg(experiment.get("dose"), experiment.get("dose_unit"))
                if pred_dose is not None and exp_dose is not None and abs(pred_dose - exp_dose) > max(1e-9, 1e-6 * max(abs(pred_dose), abs(exp_dose), 1.0)):
                    continue
            pv, ev = _number(prediction.get("display_value")), _number(experiment.get("normalized_value"))
            if pv is None or ev is None: continue
            diff = pv - ev
            direct.append({"status": "DIRECT" if experiment["comparability"] == DIRECT else "CONVERTED", "comparability": experiment["comparability"], "prediction_value": pv, "experimental_value": ev, "difference": diff, "signed_error": diff, "absolute_error": abs(diff), "preview": experiment.get("state") == "EXTERNAL_CANDIDATE", "experimental_id": experiment.get("id"), "unit": prediction.get("unit") or experiment.get("normalized_unit")})
        elif experiment.get("comparability") == RELATED: related.append(experiment)
    if direct: return {"status": direct[0]["status"], "comparability": direct[0]["comparability"], "matches": direct, **direct[0]}
    if related: return {"status": "RELATED_SAME_SCIENTIFIC_GROUP", "reason": "Related measurement semantics; no numeric error calculated", "difference": None, "related_observation_count": len(related)}
    return None


def _pk_internal_item(study, nca, raw_endpoint, value, unit, mapped):
    """Convert a persisted NCA result into the same evidence contract as search data."""
    context = {
        "study": study.study_name,
        "route": study.route,
        "dose": study.dose,
        "dose_unit": study.dose_unit,
        "dosing_frequency": study.dosing_frequency,
        "fed_fasted": study.fed_fasted,
        "matrix": study.matrix,
        "strain": study.strain,
    }
    endpoint = mapped["canonical_endpoint_id"]
    return {
        "id": nca.id,
        "origin": "INTERNAL_EXPERIMENTAL",
        "state": "INTERNAL_EXPERIMENTAL",
        "raw_endpoint": raw_endpoint,
        "endpoint": _display_name(endpoint, raw_endpoint),
        "raw_value": value,
        "normalized_value": mapped.get("normalized_value"),
        "raw_unit": unit,
        "normalized_unit": mapped.get("normalized_unit", unit),
        "relation": "=",
        "species": mapped.get("species"),
        "route": mapped.get("route"),
        "context": context,
        "reference": {"source": study.source or "Internal PK study", "reference": study.notes or study.study_name, "study_id": study.id},
        "qualification": "QUALIFIED_DIRECT" if mapped.get("comparability_status") in {DIRECT, CONVERTED} else "NEEDS_REVIEW",
        "comparability": mapped.get("comparability_status", UNSUPPORTED),
        "importable": False,
        "adaptation_eligibility": mapped.get("comparability_status") in {DIRECT, CONVERTED},
        "canonical_endpoint_id": endpoint,
        "canonical_comparison_key": mapped["comparison_key"],
    }


def _add_experiment(row, item):
    if item["comparability"] == RELATED: row["related_evidence"].append(item)
    elif item["qualification"] in {"NEEDS_REVIEW", "UNSUPPORTED"} or item["comparability"] == UNSUPPORTED: row["needs_review"].append(item)
    elif item["state"] == "EXTERNAL_IMPORTED": row["experimental_external_imported"].append(item)
    else: row["experimental_external_candidates"].append(item)
    row["references"].append(item["reference"])


def build_endpoint_comparison(db, version_id: int) -> dict:
    version = db.get(CompoundVersion, version_id)
    if not version: raise ValueError("CompoundVersion not found")
    compound = db.get(Compound, version.compound_row_id)
    endpoint_rows = {}
    predictions = db.scalars(select(ADMETPrediction).join(ADMETModelRegistry).where(ADMETPrediction.version_id == version_id, ADMETPrediction.execution_status == "SUCCESS").order_by(ADMETPrediction.created_at.desc())).all()
    latest_snapshots = {}
    for snap in db.scalars(select(PredictionEndpointSnapshot).where(PredictionEndpointSnapshot.compound_version_id == version_id).order_by(PredictionEndpointSnapshot.created_at.desc())).all(): latest_snapshots.setdefault(snap.endpoint_name, snap)
    latest_canonical_snapshots = {}
    for snap in db.scalars(select(PredictionEndpointSnapshot).where(PredictionEndpointSnapshot.compound_version_id == version_id).order_by(PredictionEndpointSnapshot.created_at.desc())).all(): latest_canonical_snapshots.setdefault(snap.endpoint_id, snap)
    by_raw = defaultdict(list)
    for pred in predictions: by_raw[pred.model.endpoint_name].append(pred)
    for raw_endpoint, rows in by_raw.items():
        mapping = canonicalize_prediction_endpoint(raw_endpoint, species=(rows[0].model.species or ""))
        eid = mapping["canonical_endpoint_id"]; row = endpoint_rows.setdefault(eid, _blank(eid, raw_endpoint)); row["section"] = _section(eid); row["display_name"] = _display_name(eid, raw_endpoint); row["species"] = mapping["species"]; row["route"] = mapping["route"]; row["canonical_comparison_key"] = mapping["comparison_key"]; row["prediction"] = _prediction_object(latest_snapshots.get(raw_endpoint), rows, eid, raw_endpoint, species=mapping["species"], route=mapping["route"])

    endpoint_names = {row.id: row.name for row in db.scalars(select(ADMETEndpoint).where(ADMETEndpoint.project_id == compound.project_id)).all()}
    for measurement in db.scalars(select(ADMETMeasurement).where(ADMETMeasurement.version_id == version_id)).all():
        raw_endpoint = endpoint_names.get(measurement.endpoint_id, "Unknown endpoint")
        mapped = normalize_experimental_observation(raw_endpoint, measurement.value if measurement.value is not None else measurement.qualitative_value, measurement.unit, species=measurement.species, context={"matrix": measurement.matrix, "method": measurement.method})
        eid = mapped["canonical_endpoint_id"]; row = endpoint_rows.setdefault(eid, _blank(eid, raw_endpoint)); row["experimental_internal"].append({"id": measurement.id, "origin": "INTERNAL_EXPERIMENTAL", "state": "INTERNAL_EXPERIMENTAL", "raw_endpoint": raw_endpoint, "endpoint": row["display_name"], "raw_value": measurement.value if measurement.value is not None else measurement.qualitative_value, "normalized_value": mapped.get("normalized_value", measurement.value), "raw_unit": measurement.unit, "normalized_unit": mapped.get("normalized_unit", measurement.unit), "relation": measurement.qualifier, "species": mapped.get("species"), "route": mapped.get("route"), "context": {"matrix": measurement.matrix, "method": measurement.method}, "reference": {"source": measurement.source, "reference": measurement.notes}, "qualification": "QUALIFIED_DIRECT", "comparability": mapped.get("comparability_status", DIRECT), "importable": False, "adaptation_eligibility": True, "canonical_endpoint_id": eid, "canonical_comparison_key": mapped["comparison_key"]})

    assays = {row.id: row for row in db.scalars(select(AssayDefinition).where(AssayDefinition.project_id == compound.project_id)).all()}
    activity_preds = db.scalars(select(ActivityPrediction).where(ActivityPrediction.version_id == version_id).order_by(ActivityPrediction.created_at.desc())).all(); activity_meas = db.scalars(select(ActivityMeasurement).where(ActivityMeasurement.version_id == version_id).order_by(ActivityMeasurement.created_at.desc())).all()
    for assay_id, assay in assays.items():
        pred = next((p for p in activity_preds if p.assay_id == assay_id), None); measured = [m for m in activity_meas if m.assay_id == assay_id]
        if pred is None and not measured: continue
        eid = f"ACTIVITY_{str(assay.measurement_type).upper()}:{assay_id}"; row = endpoint_rows.setdefault(eid, _blank(eid, f"{assay.name} ({assay.measurement_type})")); row["section"] = "ACTIVITY"; row["display_name"] = f"{assay.name} ({assay.measurement_type})"; row["species"] = normalize_species(assay.species)
        if pred is not None: row["prediction"] = {"available": True, "raw_endpoint": assay.measurement_type, "canonical_endpoint_id": eid, "base_value": pred.predicted_value_nm, "project_value": None, "display_value": pred.predicted_value_nm, "unit": assay.unit, "prediction_type": pred.prediction_type, "maturity": {"level": 1, "label": "Base Prediction", "stars": "★☆☆☆☆"}, "ood": pred.applicability_domain, "timestamp": _iso(pred.created_at), "model_count": 1}
        for m in measured: row["experimental_internal"].append({"id": m.id, "origin": "INTERNAL_EXPERIMENTAL", "state": "INTERNAL_EXPERIMENTAL", "raw_endpoint": assay.measurement_type, "raw_value": m.raw_value, "normalized_value": m.normalized_value_nm, "raw_unit": m.original_unit, "normalized_unit": "nM", "relation": m.qualifier, "species": normalize_species(assay.species), "context": {"target": assay.target, "cell_line": assay.cell_line, "assay": assay.name}, "reference": {"source": m.source, "reference": m.notes}, "qualification": "QUALIFIED_DIRECT", "comparability": DIRECT, "importable": False, "adaptation_eligibility": True})

    for evidence in db.scalars(select(ExternalExperimentalEvidence).where(ExternalExperimentalEvidence.compound_version_id == version_id)).all():
        item, mapped = _mapped_external(evidence); eid = mapped["canonical_endpoint_id"]; row = endpoint_rows.setdefault(eid, _blank(eid, evidence.raw_endpoint_name)); row["section"] = mapped["section"]; row["display_name"] = _display_name(eid, evidence.raw_endpoint_name); row["species"] = mapped.get("species", "UNSPECIFIED"); row["route"] = mapped.get("route", "UNSPECIFIED"); row["canonical_comparison_key"] = mapped["comparison_key"]; _add_experiment(row, item)

    # Non-scalar metabolism predictions are still canonical prediction output.
    # Keep them in the scientific metabolism section, but do not pretend that
    # a count of hypotheses is a quantitative validation endpoint.
    metabolic_runs = db.scalars(select(MetabolicPredictionRun).where(MetabolicPredictionRun.version_id == version_id, MetabolicPredictionRun.status == "COMPLETE").order_by(MetabolicPredictionRun.started_at.desc())).all()
    if metabolic_runs:
        for eid, label in (("METABOLIC_SOFT_SPOTS", "Metabolic soft spots"), ("METABOLITE_HYPOTHESES", "Metabolite hypotheses")):
            snapshot = latest_canonical_snapshots.get(eid)
            if snapshot is None: continue
            row = endpoint_rows.setdefault(eid, _blank(eid, label)); row["section"] = "METABOLISM"; row["display_name"] = label
            row["prediction"] = _snapshot_prediction(snapshot, eid, label)

    # PK foundations are persisted outside ADMET. Expose only actual values.
    for pset in db.scalars(select(PKParameterSet).where(PKParameterSet.version_id == version_id).order_by(PKParameterSet.created_at.desc())).all():
        species = normalize_species(pset.species); route = "ORAL" if str(pset.route).upper() == "PO" else str(pset.route).upper()
        for parameter, value, unit, source_type in _pk_snapshot_values(pset):
            eid = f"{species}_PK_F_ORAL" if parameter == "F" else f"{species}_PK_{parameter}_{route}"; row = endpoint_rows.setdefault(eid, _blank(eid, eid)); row["section"] = "PK"; row["display_name"] = f"{species.title()} {'Oral Bioavailability F' if parameter == 'F' else parameter.replace('_', ' ')}"; row["species"] = species; row["route"] = "ORAL" if parameter == "F" else route
            snapshot = latest_canonical_snapshots.get(eid)
            if snapshot is not None:
                row["prediction"] = _snapshot_prediction(snapshot, eid, parameter, species=species, route="ORAL" if parameter == "F" else route, dose=pset.dose_value, dose_unit=pset.dose_unit)
            elif not row["prediction"].get("available"):
                mapped_pk = normalize_experimental_observation(parameter, value, unit, species=species, context={"route": "ORAL" if parameter == "F" else route, "dose": pset.dose_value, "dose_unit": pset.dose_unit})
                row["prediction"] = {"available": True, "raw_endpoint": parameter, "canonical_endpoint_id": eid, "canonical_comparison_key": f"{eid}|{species}|{'ORAL' if parameter == 'F' else route}|PARENT", "base_value": mapped_pk.get("normalized_value", value), "project_value": None, "display_value": mapped_pk.get("normalized_value", value), "unit": mapped_pk.get("normalized_unit", unit), "prediction_type": source_type, "source_type": source_type, "source_label": prediction_source_label(source_type), "maturity": {"level": 1, "label": "Base Prediction", "stars": "★☆☆☆☆"}, "timestamp": _iso(pset.created_at), "model_count": 1, "species": species, "route": "ORAL" if parameter == "F" else route, "dose": pset.dose_value, "dose_unit": pset.dose_unit}

    # Concentration-time simulations provide the Stage-5 Cmax/Tmax/AUC/t1/2
    # predictions. They are joined by the same species/route key as external
    # PK observations, rather than being hidden in the simulation tab.
    for sim in db.scalars(select(PKSimulationRun).where(PKSimulationRun.version_id == version_id).order_by(PKSimulationRun.created_at.desc())).all():
        species = normalize_species(sim.species); route = "ORAL" if str(sim.route).upper() == "PO" else str(sim.route).upper()
        for parameter, value, unit in _simulation_values(sim):
            if value is None: continue
            eid = f"{species}_PK_{parameter}_{route}"; row = endpoint_rows.setdefault(eid, _blank(eid, eid)); row["section"] = "PK"; row["display_name"] = f"{species.title()} {parameter.replace('_', ' ')}"; row["species"] = species; row["route"] = route
            snapshot = latest_canonical_snapshots.get(eid)
            if snapshot is not None:
                row["prediction"] = _snapshot_prediction(snapshot, eid, parameter, species=species, route=route, dose=sim.dose, dose_unit=sim.dose_unit)
            elif not row["prediction"].get("available"):
                row["prediction"] = {"available": True, "raw_endpoint": parameter, "canonical_endpoint_id": eid, "canonical_comparison_key": f"{eid}|{species}|{route}|PARENT", "base_value": value, "project_value": None, "display_value": value, "unit": unit, "prediction_type": "MECHANISTIC_ESTIMATE", "source_type": PREDICTION_MECHANISTIC, "source_label": prediction_source_label(PREDICTION_MECHANISTIC), "maturity": {"level": 1, "label": "Base Prediction", "stars": "★☆☆☆☆"}, "timestamp": _iso(sim.created_at), "model_count": 1, "species": species, "route": route, "dose": sim.dose, "dose_unit": sim.dose_unit, "prediction_run_id": sim.id}

    # NCA studies/results are the persisted internal experimental PK stream.
    # They must share the same species/route/parameter ontology as foundations
    # and simulations; raw study provenance remains attached to each row.
    studies = {study.id: study for study in db.scalars(select(PKStudy).where(PKStudy.version_id == version_id)).all()}
    for nca in db.scalars(select(PKNCAResult).where(PKNCAResult.version_id == version_id, PKNCAResult.is_latest.is_(True))).all():
        study = studies.get(nca.pk_study_id)
        if not study: continue
        values = [("Cmax", nca.cmax, nca.cmax_unit), ("Tmax", nca.tmax, nca.tmax_unit), ("AUC0-t", nca.auclast, nca.auclast_unit), ("AUC0-inf", nca.aucinf, nca.aucinf_unit), ("t1/2", nca.terminal_half_life, "hours"), ("CL/F" if str(study.route).upper() != "IV" and nca.cl_f is not None else "CL", nca.cl_f if str(study.route).upper() != "IV" and nca.cl_f is not None else nca.cl, nca.cl_f_unit if str(study.route).upper() != "IV" and nca.cl_f is not None else nca.cl_unit), ("Vd/F" if str(study.route).upper() != "IV" and nca.vz_f is not None else "Vd", nca.vz_f if str(study.route).upper() != "IV" and nca.vz_f is not None else nca.vz, nca.vz_f_unit if str(study.route).upper() != "IV" and nca.vz_f is not None else nca.vz_unit)]
        for raw_endpoint, value, unit in values:
            if value is None: continue
            mapped = normalize_experimental_observation(raw_endpoint, value, unit, species=study.species, context={"route": study.route, "dose": study.dose, "dose_unit": study.dose_unit, "dosing_frequency": study.dosing_frequency, "fed_fasted": study.fed_fasted, "matrix": study.matrix})
            eid = mapped["canonical_endpoint_id"]
            row = endpoint_rows.setdefault(eid, _blank(eid, raw_endpoint)); row["section"] = "PK"; row["display_name"] = mapped.get("display_name", _display_name(eid, raw_endpoint)); row["species"] = mapped.get("species", "UNSPECIFIED"); row["route"] = mapped.get("route", "UNSPECIFIED"); row["canonical_comparison_key"] = mapped["comparison_key"]
            row["experimental_internal"].append(_pk_internal_item(study, nca, raw_endpoint, value, unit, mapped))

    prediction_endpoint_ids = {eid for eid, row in endpoint_rows.items() if row["prediction"].get("available")}
    evidence_lists = ("experimental_internal", "experimental_external_imported", "experimental_external_candidates", "related_evidence", "needs_review")
    for row in endpoint_rows.values():
        for list_name in evidence_lists:
            for item in row[list_name]:
                item["qualification_details"] = qualify_record(
                    item, prediction_endpoints=prediction_endpoint_ids,
                    imported=item.get("state") == "EXTERNAL_IMPORTED",
                )
        experiments = row["experimental_internal"] + row["experimental_external_imported"] + row["experimental_external_candidates"] + row["related_evidence"]
        row["comparison"] = _comparison(row["prediction"], experiments)
        row["summary"] = {"both": int(bool(row["prediction"].get("available") and experiments)), "prediction_only": int(bool(row["prediction"].get("available") and not experiments)), "experimental_only": int(bool(experiments and not row["prediction"].get("available"))), "related": len(row["related_evidence"]), "needs_review": len(row["needs_review"]), "ready_to_import": sum(bool((item.get("qualification_details") or {}).get("stages", {}).get("IMPORTABLE", item.get("importable"))) for item in row["experimental_external_candidates"])}
    endpoints = sorted(endpoint_rows.values(), key=lambda item: (item["section"], item["display_name"], item.get("canonical_comparison_key", ""))); summary = {key: sum(row["summary"][key] for row in endpoints) for key in ("both", "prediction_only", "experimental_only", "related", "needs_review", "ready_to_import")}; summary["imported_pairs"] = 0
    qualification_items = [item for row in endpoints for list_name in evidence_lists for item in row[list_name]]
    qualification = aggregate_qualification(qualification_items, prediction_endpoints=prediction_endpoint_ids)
    summary["qualification"] = qualification["global"]
    summary["source_qualification"] = qualification["sources"]
    return {"version_id": version_id, "project_id": compound.project_id, "compound_id": compound.id, "canonical_endpoint_version": CANONICAL_ENDPOINT_VERSION, "comparison_unit_version": COMPARISON_UNIT_VERSION, "qualification_version": qualification["qualification_version"], "endpoints": endpoints, "summary": summary}
