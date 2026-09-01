"""Generate the v3.9 prediction coverage audit from persisted application data.

This is an audit/report generator.  It never runs a prediction, imports
evidence, or changes project/adaptation state.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

from sqlalchemy import select

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.admet import ADMETModelRegistry, ADMETPrediction
from backend.activity_models import ActivityPrediction
from backend.canonical_endpoints import (
    CANONICAL_ENDPOINT_VERSION, COMPARISON_UNIT_VERSION, PREDICTION_UNAVAILABLE,
    canonicalize_prediction_endpoint, prediction_source_type,
)
from backend.database import SessionLocal
from backend.endpoint_comparison import build_endpoint_comparison, _pk_snapshot_values, _simulation_values
from backend.ivive import PKParameterSet
from backend.metabolism import MetabolicPredictionRun
from backend.models import Compound, CompoundVersion, PropertyCalculation
from backend.simulation import PKSimulationRun

OUT = ROOT / "validation"
DRUGS = ("Sunvozertinib", "Osimertinib", "Midazolam", "Warfarin", "Metformin")


def _json(value):
    return json.loads(json.dumps(value, default=str))


def _coverage(db):
    rows = {}
    for compound in db.scalars(select(Compound)).all():
        version = next((item for item in compound.versions if item.version_number == compound.current_version), None)
        if not version:
            continue
        rows[version.id] = (compound, version, build_endpoint_comparison(db, version.id))
    return rows


def _row_match_counts(coverage):
    counts = defaultdict(int)
    for _, _, payload in coverage.values():
        for row in payload["endpoints"]:
            if row.get("prediction", {}).get("available"):
                counts[row["endpoint_id"]] += sum(
                    len(row.get(key, [])) for key in
                    ("experimental_internal", "experimental_external_imported", "experimental_external_candidates", "related_evidence", "needs_review")
                )
    return counts


def inventory(db, coverage):
    match_counts = _row_match_counts(coverage)
    records = []
    rendered_endpoints = {
        row["endpoint_id"]
        for _, _, payload in coverage.values()
        for row in payload["endpoints"]
        if row.get("prediction", {}).get("available")
    }

    def add(**record):
        record.setdefault("persisted", True)
        record.setdefault("canonicalized", bool(record.get("canonical_endpoint")))
        record.setdefault("rendered", record.get("canonical_endpoint") in rendered_endpoints)
        record.setdefault("experimental_match_count", match_counts.get(record.get("canonical_endpoint"), 0))
        records.append(_json(record))

    for model in db.scalars(select(ADMETModelRegistry).order_by(ADMETModelRegistry.endpoint_name, ADMETModelRegistry.id)).all():
        endpoint = canonicalize_prediction_endpoint(model.endpoint_name, species=model.species)
        preds = db.scalars(select(ADMETPrediction).where(ADMETPrediction.model_id == model.id, ADMETPrediction.execution_status == "SUCCESS")).all()
        source_type = "MODEL" if preds else PREDICTION_UNAVAILABLE
        add(source_module="ADMETPrediction/Prediction Engine v1", raw_endpoint=model.endpoint_name,
            canonical_endpoint=endpoint["canonical_endpoint_id"], prediction_type=source_type,
            species=endpoint["species"], route=endpoint["route"], unit=model.output_unit,
            value_type=model.output_type or "numeric", output_count=len(preds),
            model_id=model.id, active=bool(model.is_active),
            unavailable_reason=(model.provenance_json or {}).get("reason") if not preds else "")

    for row in db.scalars(select(PKParameterSet)).all():
        for parameter, value, unit, source_type in _pk_snapshot_values(row):
            species = str(row.species or "").upper()
            route = "ORAL" if str(row.route).upper() == "PO" else str(row.route).upper()
            endpoint = f"{species}_PK_{parameter}_{route}"
            add(source_module="IVIVE/PKParameterSet", raw_endpoint=parameter,
                canonical_endpoint=endpoint, prediction_type=source_type,
                species=species, route=route, unit=unit, value_type="numeric",
                output_count=1, dose=row.dose_value, dose_unit=row.dose_unit,
                pk_parameter_set_id=row.id)

    for row in db.scalars(select(PKSimulationRun)).all():
        species = str(row.species or "").upper()
        route = "ORAL" if str(row.route).upper() == "PO" else str(row.route).upper()
        for parameter, value, unit in _simulation_values(row):
            if value is None:
                continue
            add(source_module="Stage-5 PK simulation", raw_endpoint=parameter,
                canonical_endpoint=f"{species}_PK_{parameter}_{route}",
                prediction_type="MECHANISTIC_ESTIMATE", species=species,
                route=route, unit=unit, value_type="numeric", output_count=1,
                dose=row.dose, dose_unit=row.dose_unit, simulation_run_id=row.id)

    for row in db.scalars(select(MetabolicPredictionRun).where(MetabolicPredictionRun.status == "COMPLETE")).all():
        add(source_module="Metabolism/SyGMa+SMARTCyp", raw_endpoint="Metabolic soft spots",
            canonical_endpoint="METABOLIC_SOFT_SPOTS", prediction_type="RULE_ESTIMATE",
            species="UNSPECIFIED", route="UNSPECIFIED", unit="ranked sites",
            value_type="ranking", output_count=len(row.spots), metabolic_run_id=row.id)
        add(source_module="Metabolism/SyGMa", raw_endpoint="Predicted metabolites",
            canonical_endpoint="METABOLITE_HYPOTHESES", prediction_type="RULE_ESTIMATE",
            species="UNSPECIFIED", route="UNSPECIFIED", unit="hypotheses",
            value_type="ranking", output_count=len(row.metabolites), metabolic_run_id=row.id)

    for row in db.scalars(select(ActivityPrediction)).all():
        endpoint = f"ACTIVITY_{str(row.assay.measurement_type if row.assay else 'UNKNOWN').upper()}"
        add(source_module="Activity/QSAR", raw_endpoint=row.assay.measurement_type if row.assay else "Activity",
            canonical_endpoint=endpoint, prediction_type="MODEL", species=str(row.assay.species if row.assay else "UNSPECIFIED"),
            route="UNSPECIFIED", unit="nM", value_type="numeric", output_count=1)

    for row in db.scalars(select(PropertyCalculation)).all():
        add(source_module="Properties/RDKit", raw_endpoint=row.endpoint,
            canonical_endpoint="PROPERTY_" + str(row.endpoint).upper(), prediction_type="DERIVED_ESTIMATE",
            species="CHEMICAL", route="UNSPECIFIED", unit="", value_type="numeric", output_count=1,
            method=row.method, engine=row.engine, rendered=False)

    return records


def endpoint_rows_for(compound, payload):
    return [{
        "endpoint": row["endpoint_id"], "display_name": row["display_name"],
        "section": row["section"], "species": row.get("species"), "route": row.get("route"),
        "experimental_available": bool(row.get("experimental_internal") or row.get("experimental_external_imported") or row.get("experimental_external_candidates") or row.get("related_evidence") or row.get("needs_review")),
        "prediction_available": bool(row.get("prediction", {}).get("available")),
        "prediction_type": row.get("prediction", {}).get("source_type") or row.get("prediction", {}).get("prediction_type"),
        "prediction_value": row.get("prediction", {}).get("display_value"),
        "prediction_unit": row.get("prediction", {}).get("unit"),
        "comparison_status": (row.get("comparison") or {}).get("status"),
        "direct_comparison": (row.get("comparison") or {}).get("status") in {"DIRECT", "CONVERTED"},
        "related": len(row.get("related_evidence", [])),
        "experimental_observations": sum(len(row.get(key, [])) for key in ("experimental_internal", "experimental_external_imported", "experimental_external_candidates")),
        "needs_review": len(row.get("needs_review", [])),
        "gap_reasons": sorted({item.get("routing_reason", "") for item in row.get("needs_review", []) if item.get("routing_reason")}),
    } for row in payload["endpoints"]]


def main():
    with SessionLocal() as db:
        coverage = _coverage(db)
        records = inventory(db, coverage)

        selected = {}
        for name in DRUGS:
            matches = [(compound, version, payload) for compound, version, payload in coverage.values() if compound.name.lower() == name.lower()]
            if matches:
                compound, version, payload = matches[0]
                selected[name] = {"compound_id": compound.id, "compound_version_id": version.id, "rows": endpoint_rows_for(compound, payload), "summary": payload["summary"]}
            else:
                selected[name] = {"status": "NOT_PRESENT_IN_CURRENT_DATABASE", "rows": [], "summary": {}}

        sun = selected["Sunvozertinib"]
        sun_rows = sun.get("rows", [])
        # v3.8B is retained as a supplied comparison baseline; all v3.9
        # values below are read from the persisted database at generation time.
        sun_audit = {
            "drug": "Sunvozertinib", "compound_id": sun.get("compound_id"),
            "compound_version_id": sun.get("compound_version_id"),
            "canonical_endpoint_version": CANONICAL_ENDPOINT_VERSION,
            "comparison_unit_version": COMPARISON_UNIT_VERSION,
            "before_v3_8b": {"both": 3, "direct": 1, "related": 3, "experimental_only": 17, "prediction_only": 44, "needs_review": 46},
            "after_v3_9": {"summary": sun.get("summary", {}), "endpoints": sun_rows},
            "species_pk_pair_status": "NO_CURRENT_VALID_SPECIES_PK_PAIR",
            "species_pk_pair_reason": "Persisted Sunvozertinib PK foundation has no Cmax/AUC/Tmax/half-life simulation rows; available external PK observations lack a complete species+route+dose+parameter match to foundation predictions.",
        }

        experimental_only = defaultdict(lambda: {"observations": 0, "compounds": set(), "sections": set(), "reasons": set()})
        for name, value in selected.items():
            if not value.get("rows"): continue
            for row in value["rows"]:
                if row["experimental_available"] and not row["prediction_available"]:
                    item = experimental_only[row["endpoint"]]
                    item["observations"] += row["experimental_observations"]
                    item["compounds"].add(value.get("compound_id"))
                    item["sections"].add(row["section"])
                    item["reasons"].update(row.get("gap_reasons", []))
        gaps = []
        for endpoint, item in sorted(experimental_only.items()):
            if item["observations"] == 0:
                continue
            gaps.append({"endpoint": endpoint, "experimental_observations": item["observations"], "compounds": len(item["compounds"]), "sections": sorted(item["sections"]), "source_quality": "persisted qualified/candidate evidence audit", "current_prediction_availability": "NO_MATCHING_PREDICTION", "potential_future_model_or_data_path": "Add or qualify an endpoint-specific validated method only after semantic/context contract is defined.", "priority": "HIGH" if endpoint.startswith(("PK_", "HUMAN_PK", "RAT_PK", "DOG_PK")) else "MEDIUM", "remaining_reasons": sorted(item["reasons"])})

        OUT.mkdir(exist_ok=True)
        (OUT / "prediction_endpoint_inventory_v3_9.json").write_text(json.dumps({"version": "v3.9", "canonical_endpoint_version": CANONICAL_ENDPOINT_VERSION, "comparison_unit_version": COMPARISON_UNIT_VERSION, "generated_from": "persisted database", "records": records}, indent=2, sort_keys=True) + "\n")
        (OUT / "prediction_model_gap_v3_9.json").write_text(json.dumps({"version": "v3.9", "policy": "no fabricated prediction values", "gaps": gaps}, indent=2, sort_keys=True) + "\n")
        (OUT / "sunvozertinib_prediction_coverage_v3_9.json").write_text(json.dumps({"version": "v3.9", "audit": sun_audit, "validation_set": selected}, indent=2, sort_keys=True) + "\n")
        print(json.dumps({"inventory": len(records), "model_gaps": len(gaps), "sunvozertinib": sun.get("summary", {}), "validation_drugs": {key: value.get("status", "PRESENT") for key, value in selected.items()}}, indent=2))


if __name__ == "__main__":
    main()
