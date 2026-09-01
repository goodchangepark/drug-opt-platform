"""Generate read-only v4.2 scientific validation artifacts from persisted data."""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from sqlalchemy import select

from backend.database import SessionLocal
from backend.endpoint_comparison import build_endpoint_comparison, pk_f_prediction_is_quantitative
from backend.ivive import PKParameterSet
from backend.models import Compound, CompoundVersion
from backend.admet import ADMETPredictionRun, PredictionEndpointSnapshot
from backend.platform_info import version_history


OUT = ROOT / "validation"
ENGINE = "drugopt-prediction-engine-v1@1.0.0"
ENGINE_HASH = "12757ab197b5a70d8ea1754678d9a342ab0b6ea0d82f2896bebb767d686bbdeb"


def dump(name, payload):
    (OUT / name).write_text(json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True, default=str) + "\n", encoding="utf-8")


PK_CONTRACT = {
    "contract_version": "drugopt-pk-scientific-contract-v4.2",
    "engine_policy": ENGINE,
    "engine_hash": ENGINE_HASH,
    "endpoints": [
        {"parameter": "CL", "prediction_type": "MECHANISTIC_ESTIMATE", "source_module": "backend.ivive", "scientific_definition": "Systemic clearance for the assembled route; not CL/F.", "required_inputs": ["species physiology", "scaled Clint", "fu,b", "hepatic blood flow"], "missing_input_behavior": "MODEL_UNAVAILABLE", "species": "endpoint-specific", "route": "systemic route context", "dose_requirement": "not dose-normalized", "unit": "mL/min/kg", "fallback_possible": False, "bounded_or_clipped": False, "comparison_allowed": "only matching species, parameter, route and context", "known_limitations": ["renal/non-hepatic clearance is not modeled"]},
        {"parameter": "CL/F", "prediction_type": "MECHANISTIC_ESTIMATE", "source_module": "backend.ivive", "scientific_definition": "Apparent oral clearance only when the oral route is explicit.", "required_inputs": ["systemic CL", "oral route"], "missing_input_behavior": "MODEL_UNAVAILABLE", "species": "endpoint-specific", "route": "PO", "dose_requirement": "context retained; no implicit dose scaling", "unit": "mL/min/kg", "fallback_possible": False, "bounded_or_clipped": False, "comparison_allowed": "matching oral apparent-clearance context only", "known_limitations": ["not interchangeable with systemic CL"]},
        {"parameter": "Vd/Vss", "prediction_type": "DERIVED_ESTIMATE", "source_module": "backend.ivive", "scientific_definition": "Systemic volume/distribution estimate.", "required_inputs": ["lipophilicity/ionization", "binding", "species physiology"], "missing_input_behavior": "MODEL_UNAVAILABLE", "species": "endpoint-specific", "route": "systemic", "dose_requirement": "not dose-normalized", "unit": "L/kg", "fallback_possible": False, "bounded_or_clipped": False, "comparison_allowed": "matching parameter and species only", "known_limitations": ["do not merge with Vd/F"]},
        {"parameter": "Vd/F", "prediction_type": "DERIVED_ESTIMATE", "source_module": "backend.ivive", "scientific_definition": "Apparent extravascular volume only when route semantics support it.", "required_inputs": ["systemic volume", "oral route"], "missing_input_behavior": "MODEL_UNAVAILABLE", "species": "endpoint-specific", "route": "PO", "dose_requirement": "context retained", "unit": "L/kg", "fallback_possible": False, "bounded_or_clipped": False, "comparison_allowed": "matching oral apparent-volume context only", "known_limitations": ["not interchangeable with systemic Vd/Vss"]},
        {"parameter": "F", "prediction_type": "MECHANISTIC_ESTIMATE", "source_module": "backend.ivive", "scientific_definition": "Absolute oral bioavailability F = Fa × Fg × Fh relative to an IV reference.", "required_inputs": ["quantitative Fa", "quantitative Fg", "quantitative Fh"], "missing_input_behavior": "MODEL_UNAVAILABLE / INSUFFICIENT_INPUT", "species": "endpoint-specific", "route": "PO; IV is reference arm only", "dose_requirement": "route/dose context retained", "unit": "%", "fallback_possible": False, "bounded_or_clipped": True, "comparison_allowed": "matching species, oral route, analyte and context", "known_limitations": ["missing Fg must not be replaced by 1.0", "IV F=100% is not an oral prediction"]},
        {"parameter": "Cmax", "prediction_type": "MECHANISTIC_ESTIMATE", "source_module": "backend.simulation", "scientific_definition": "Maximum concentration from persisted concentration-time simulation.", "required_inputs": ["PK simulation run", "dose", "route", "concentration-time model"], "missing_input_behavior": "MODEL_UNAVAILABLE", "species": "simulation-specific", "route": "simulation-specific", "dose_requirement": "required", "unit": "ng/mL", "fallback_possible": False, "bounded_or_clipped": False, "comparison_allowed": "matching dose/regimen/context", "known_limitations": ["not fabricated from static foundation parameters"]},
        {"parameter": "AUC", "prediction_type": "MECHANISTIC_ESTIMATE", "source_module": "backend.simulation", "scientific_definition": "Exposure integral from persisted concentration-time simulation.", "required_inputs": ["PK simulation run", "dose", "route", "concentration-time model"], "missing_input_behavior": "MODEL_UNAVAILABLE", "species": "simulation-specific", "route": "simulation-specific", "dose_requirement": "required", "unit": "ng*h/mL", "fallback_possible": False, "bounded_or_clipped": False, "comparison_allowed": "matching AUC subtype and regimen", "known_limitations": ["AUC0-t, AUC0-inf and AUCtau remain distinct"]},
        {"parameter": "Tmax", "prediction_type": "MECHANISTIC_ESTIMATE", "source_module": "backend.simulation", "scientific_definition": "Time of maximum concentration from simulation.", "required_inputs": ["PK simulation run", "dose", "route", "absorption model"], "missing_input_behavior": "MODEL_UNAVAILABLE", "species": "simulation-specific", "route": "simulation-specific", "dose_requirement": "required", "unit": "hours", "fallback_possible": False, "bounded_or_clipped": False, "comparison_allowed": "matching dose/regimen/context", "known_limitations": ["not derived from an unpersisted curve"]},
        {"parameter": "t1/2", "prediction_type": "DERIVED_ESTIMATE", "source_module": "backend.simulation", "scientific_definition": "Terminal or model half-life with explicit subtype.", "required_inputs": ["persisted simulation or valid CL/V context"], "missing_input_behavior": "MODEL_UNAVAILABLE", "species": "context-specific", "route": "context-specific", "dose_requirement": "model-specific", "unit": "hours", "fallback_possible": False, "bounded_or_clipped": False, "comparison_allowed": "matching half-life subtype", "known_limitations": ["terminal and other half-life semantics must not be collapsed"]},
    ],
}


def main():
    with SessionLocal() as db:
        compound = db.scalar(select(Compound).where(Compound.name.ilike("%sunvozertinib%")))
        if not compound:
            raise SystemExit("Sunvozertinib was not found in the persisted database")
        version = next(v for v in compound.versions if v.version_number == compound.current_version)
        comparison = build_endpoint_comparison(db, version.id)
        psets = db.scalars(select(PKParameterSet).where(PKParameterSet.version_id == version.id).order_by(PKParameterSet.created_at)).all()
        snapshots = db.scalars(select(PredictionEndpointSnapshot).where(PredictionEndpointSnapshot.compound_version_id == version.id).order_by(PredictionEndpointSnapshot.created_at)).all()
        runs = db.scalars(select(ADMETPredictionRun).where(ADMETPredictionRun.version_id == version.id).order_by(ADMETPredictionRun.started_at)).all()

        pk_rows = []
        for row in comparison["endpoints"]:
            if row.get("section") != "PK":
                continue
            pk_rows.append({
                "endpoint": row["endpoint_id"], "display_name": row["display_name"],
                "species": row.get("species"), "route": row.get("route"),
                "experimental": [item for name in ("experimental_internal", "experimental_external_imported", "experimental_external_candidates", "related_evidence", "needs_review") for item in row.get(name, [])],
                "prediction": row.get("prediction"), "comparison": row.get("comparison"),
            })

        f_audit = []
        pset_by_id = {row.id: row for row in psets}
        for snap in snapshots:
            if not str(snap.endpoint_id).endswith("_PK_F_ORAL"):
                continue
            data = snap.snapshot_json or {}
            pset = pset_by_id.get(data.get("pk_parameter_set_id"))
            valid = bool(pset and pk_f_prediction_is_quantitative(pset))
            f_audit.append({
                "snapshot_id": snap.id, "endpoint": snap.endpoint_id, "value": snap.base_value,
                "unit": snap.base_unit, "prediction_type": snap.prediction_type,
                "route_in_snapshot": data.get("route"), "parameter_set_id": data.get("pk_parameter_set_id"),
                "parameter_set_route": pset.route if pset else None,
                "valid_current_oral_f": valid,
                "reason": "VALID_QUANTITATIVE_ORAL_F" if valid else "IV_REFERENCE_OR_FG_MISSING_OR_LEGACY_FALLBACK",
                "provenance": data.get("provenance", {}),
            })

        direct_pairs = []
        for row in comparison["endpoints"]:
            if (row.get("comparison") or {}).get("status") in {"DIRECT", "CONVERTED"}:
                direct_pairs.append({"endpoint": row["endpoint_id"], "comparison": row["comparison"], "prediction": row.get("prediction")})

        dump("pk_prediction_scientific_contract_v4_2.json", PK_CONTRACT)
        dump("sunvozertinib_pk_pair_audit_v4_2.json", {
            "generated_at": datetime.now(timezone.utc).isoformat(), "compound": compound.name, "version_id": version.id,
            "engine": ENGINE, "engine_hash": ENGINE_HASH, "rows": pk_rows, "bioavailability_snapshot_audit": f_audit,
            "concentration_time_gap": {"cmax": False, "auc": False, "tmax": False, "t_half": False, "reason": "No persisted concentration-time simulation output for this compound/context."},
        })
        dump("prediction_performance_profile_v4_2.json", {
            "generated_at": datetime.now(timezone.utc).isoformat(), "compound": compound.name,
            "performance_policy": "PERFORMANCE_NOT_CALIBRATED", "semantic_match_is_not_accuracy": True,
            "direct_pairs": direct_pairs, "direct_pair_count": len(direct_pairs),
            "note": "Numeric error is reported without arbitrary good/moderate/bad labels; independent external observations remain individually preserved.",
        })
        hashes = [run.inputs_hash for run in runs]
        dump("prediction_run_dedup_audit_v4_2.json", {
            "generated_at": datetime.now(timezone.utc).isoformat(), "compound": compound.name, "version_id": version.id,
            "prediction_runs": len(runs), "unique_fingerprints": len(set(hashes)), "identical_fingerprint_duplicates": len(hashes) - len(set(hashes)),
            "runs": [{"run_id": run.id, "requested_by": run.requested_by, "status": run.status, "fingerprint": run.inputs_hash, "started_at": run.started_at, "completed_at": run.completed_at} for run in runs],
            "future_policy": "Identical workflow requests reuse a completed workflow fingerprint by default; force_rerun=true creates an explicit immutable rerun.",
        })
        dump("help_version_history_audit_v4_2.json", {
            "current_product_version": "Drug-OPT v1.0",
            "engine_version": ENGINE, "first_post_release_milestone": "v3.5",
            "entries": version_history(), "latest_entry": version_history()[-1]["version"],
            "basis": "Curated user-visible milestones supported by git history and validation artifacts; individual maintenance commits are grouped.",
        })


if __name__ == "__main__":
    main()
