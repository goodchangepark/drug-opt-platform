"""Canonical endpoint aggregation for the persisted comparison API.

Raw source labels are audit fields, never comparison keys. Both prediction
snapshots and experimental evidence pass through the versioned registry.
"""
from __future__ import annotations

from collections import defaultdict
from sqlalchemy import select

from .admet import ADMETEndpoint, ADMETMeasurement, ADMETPrediction, ADMETModelRegistry, PredictionEndpointSnapshot
from .activity_models import ActivityMeasurement, ActivityPrediction, AssayDefinition
from .canonical_endpoints import (
    CANONICAL_ENDPOINT_VERSION, COMPARISON_UNIT_VERSION, CONVERTED, DIRECT,
    RELATED, REGISTRY, UNSUPPORTED, canonicalize_prediction_endpoint,
    endpoint_contract, normalize_experimental_observation, normalize_species,
)
from .ivive import PKParameterSet
from .pk import PKNCAResult, PKStudy
from .simulation import PKSimulationRun
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
    item = {"id": row.id, "origin": state, "state": state, "raw_endpoint": row.raw_endpoint_name, "endpoint": _display_name(endpoint_id, row.raw_endpoint_name), "raw_value": row.raw_value, "raw_unit": row.raw_unit, "relation": row.raw_relation, "normalized_value": mapped.get("normalized_value"), "normalized_unit": mapped.get("normalized_unit", ""), "species": mapped.get("species", normalize_species(row.species, context)), "route": mapped.get("route", "UNSPECIFIED"), "context": context, "reference": _reference(row), "qualification": qualification, "comparability": mapped["comparability_status"], "importable": comparable and state != "EXTERNAL_IMPORTED", "identity_match_status": row.identity_match_status, "reference_status": "REFERENCE_RESOLVED_IMPORTED" if state == "EXTERNAL_IMPORTED" else "REFERENCE_RESOLVED_CANDIDATE", "assay_type": row.assay_type, "assay_id": row.source_assay_id, "adaptation_eligibility": bool(row.accepted_at and comparable), "display_evidence_group_id": row.provenance_fingerprint or f"evidence-{row.id}", "independent_experiment_group_id": row.source_document_id or row.source_record_id or f"evidence-{row.id}", "canonical_endpoint_id": endpoint_id, "canonical_comparison_key": mapped["comparison_key"], "display_source": row.source_database, "routing_reason": mapped.get("reason", ""), "normalization_rule": mapped.get("normalization_rule", ""), "raw_persisted_canonical_endpoint_id": row.canonical_endpoint_id}
    return item, mapped


def _prediction_object(snapshot, rows, endpoint_id, raw_endpoint, species="", route=""):
    first = rows[0]
    snap = dict(snapshot.snapshot_json or {}) if snapshot else dict(first.outputs_json or {}).get("prediction_snapshot") or {}
    base = snapshot.base_value if snapshot and snapshot.base_value is not None else snap.get("base_prediction", first.predicted_value)
    project = snapshot.project_value if snapshot else snap.get("project_prediction")
    display = project if project is not None else base
    unit = snapshot.base_unit if snapshot and snapshot.base_unit else first.unit
    return {"available": True, "prediction_snapshot_id": snapshot.id if snapshot else None, "prediction_run_id": snapshot.prediction_run_id if snapshot else first.run_id, "raw_endpoint": raw_endpoint, "canonical_endpoint_id": endpoint_id, "canonical_comparison_key": f"{endpoint_id}|{species}|{route}", "base_value": base, "project_value": project, "display_value": display, "unit": unit, "prediction_type": snap.get("prediction_type", "BASE_PREDICTION"), "adapter": snap.get("adapter_version", ""), "maturity": snap.get("maturity") or {"level": 1, "label": "Base Prediction", "stars": "★☆☆☆☆"}, "ood": first.applicability_domain, "timestamp": _iso(snapshot.created_at if snapshot else first.created_at), "model_count": len(rows), "model_predictions": {str(row.model_id): row.predicted_value for row in rows}, "species": species, "route": route}


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

    # PK foundations are persisted outside ADMET. Expose only actual values.
    for pset in db.scalars(select(PKParameterSet).where(PKParameterSet.version_id == version_id).order_by(PKParameterSet.created_at.desc())).all():
        species = normalize_species(pset.species); route = "ORAL" if str(pset.route).upper() == "PO" else str(pset.route).upper(); f_unit = "%" if pset.f_predicted is not None and abs(float(pset.f_predicted)) > 1 else "fraction"; values = [("CL/F" if route == "ORAL" else "CL", pset.cl_value, pset.cl_unit), ("Vd/F" if route == "ORAL" else "Vd", pset.v_value, pset.v_unit), ("F", pset.f_predicted, f_unit)]
        for parameter, value, unit in values:
            if value is None: continue
            eid = f"{species}_PK_{parameter}_{route}"; row = endpoint_rows.setdefault(eid, _blank(eid, eid)); row["section"] = "PK"; row["display_name"] = f"{species.title()} {parameter.replace('_', ' ')}"; row["species"] = species; row["route"] = route
            mapped_pk = normalize_experimental_observation(parameter, value, unit, species=species, context={"route": route, "dose": pset.dose_value, "dose_unit": pset.dose_unit})
            if not row["prediction"].get("available"): row["prediction"] = {"available": True, "raw_endpoint": parameter, "canonical_endpoint_id": eid, "canonical_comparison_key": f"{eid}|{species}|{route}|PARENT", "base_value": mapped_pk.get("normalized_value", value), "project_value": None, "display_value": mapped_pk.get("normalized_value", value), "unit": mapped_pk.get("normalized_unit", unit), "prediction_type": "PK_FOUNDATION", "maturity": {"level": 1, "label": "Base Prediction", "stars": "★☆☆☆☆"}, "timestamp": _iso(pset.created_at), "model_count": 1, "species": species, "route": route, "dose": pset.dose_value, "dose_unit": pset.dose_unit}

    # Concentration-time simulations provide the Stage-5 Cmax/Tmax/AUC/t1/2
    # predictions. They are joined by the same species/route key as external
    # PK observations, rather than being hidden in the simulation tab.
    for sim in db.scalars(select(PKSimulationRun).where(PKSimulationRun.version_id == version_id).order_by(PKSimulationRun.created_at.desc())).all():
        species = normalize_species(sim.species); route = "ORAL" if str(sim.route).upper() == "PO" else str(sim.route).upper(); metrics = sim.output_metrics or {}
        values = [("CMAX", metrics.get("cmax_ng_ml"), "ng/mL"), ("TMAX", metrics.get("tmax_hours"), "hours"), ("AUC0_T", metrics.get("auc_last_ng_h_ml"), "ng*h/mL"), ("AUC0_INF", metrics.get("auc_inf_analytical_ng_h_ml"), "ng*h/mL"), ("T_HALF", metrics.get("half_life_hours"), "hours")]
        for parameter, value, unit in values:
            if value is None: continue
            eid = f"{species}_PK_{parameter}_{route}"; row = endpoint_rows.setdefault(eid, _blank(eid, eid)); row["section"] = "PK"; row["display_name"] = f"{species.title()} {parameter.replace('_', ' ')}"; row["species"] = species; row["route"] = route
            if not row["prediction"].get("available"): row["prediction"] = {"available": True, "raw_endpoint": parameter, "canonical_endpoint_id": eid, "canonical_comparison_key": f"{eid}|{species}|{route}|PARENT", "base_value": value, "project_value": None, "display_value": value, "unit": unit, "prediction_type": "PK_SIMULATION", "maturity": {"level": 1, "label": "Base Prediction", "stars": "★☆☆☆☆"}, "timestamp": _iso(sim.created_at), "model_count": 1, "species": species, "route": route, "dose": sim.dose, "dose_unit": sim.dose_unit, "prediction_run_id": sim.id}

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

    for row in endpoint_rows.values():
        experiments = row["experimental_internal"] + row["experimental_external_imported"] + row["experimental_external_candidates"] + row["related_evidence"]
        row["comparison"] = _comparison(row["prediction"], experiments)
        row["summary"] = {"both": int(bool(row["prediction"].get("available") and experiments)), "prediction_only": int(bool(row["prediction"].get("available") and not experiments)), "experimental_only": int(bool(experiments and not row["prediction"].get("available"))), "related": len(row["related_evidence"]), "needs_review": len(row["needs_review"]), "ready_to_import": sum(bool(item.get("importable")) for item in row["experimental_external_candidates"])}
    endpoints = sorted(endpoint_rows.values(), key=lambda item: (item["section"], item["display_name"], item.get("canonical_comparison_key", ""))); summary = {key: sum(row["summary"][key] for row in endpoints) for key in ("both", "prediction_only", "experimental_only", "related", "needs_review", "ready_to_import")}; summary["imported_pairs"] = 0
    return {"version_id": version_id, "project_id": compound.project_id, "compound_id": compound.id, "canonical_endpoint_version": CANONICAL_ENDPOINT_VERSION, "comparison_unit_version": COMPARISON_UNIT_VERSION, "endpoints": endpoints, "summary": summary}
