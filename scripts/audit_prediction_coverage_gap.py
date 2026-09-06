"""Generate the v5.1 prediction-coverage gap audit without searching or predicting."""
from __future__ import annotations

import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from sqlalchemy import select

from backend.database import SessionLocal
from backend.endpoint_comparison import build_endpoint_comparison, _pk_snapshot_values, _simulation_values
from backend.ivive import PKParameterSet
from backend.models import Compound, CompoundVersion
from backend.prediction_coverage import build_prediction_coverage
from backend.simulation import PKSimulationRun

OUT = Path(__file__).resolve().parents[1] / "validation" / "prediction_coverage_gap_closure_v5_1.json"


def _calculated_ids(db, version_id: int) -> set[str]:
    result = set()
    for p in db.scalars(select(PKParameterSet).where(PKParameterSet.version_id == version_id)).all():
        species = str(p.species).upper()
        route = "ORAL" if str(p.route).upper() == "PO" else str(p.route).upper()
        for parameter, value, _unit, _source in _pk_snapshot_values(p):
            if value is not None:
                result.add(f"{species}_PK_{parameter}_{'ORAL' if parameter == 'F' else route}")
    for sim in db.scalars(select(PKSimulationRun).where(PKSimulationRun.version_id == version_id)).all():
        species = str(sim.species).upper()
        route = "ORAL" if str(sim.route).upper() == "PO" else str(sim.route).upper()
        for parameter, value, _unit in _simulation_values(sim):
            if value is not None:
                result.add(f"{species}_PK_{parameter}_{route}")
    return result


def main() -> None:
    db = SessionLocal()
    try:
        version = db.scalar(select(CompoundVersion).join(Compound).where(Compound.project_id == 3, Compound.name.ilike("Sunvozertinib")).order_by(CompoundVersion.id.desc()))
        if not version:
            raise SystemExit("Sunvozertinib not found in Project 3")
        golden = build_endpoint_comparison(db, version.id)
        calculated = _calculated_ids(db, version.id)
        gap = build_prediction_coverage(golden, calculated_endpoints=calculated)
        compounds = []
        for name in ("Osimertinib", "Midazolam", "Warfarin", "Metformin", "Mobocertinib"):
            row = db.scalar(select(Compound).where(Compound.name.ilike(name)).order_by(Compound.id.desc()))
            if not row:
                continue
            version_row = db.scalar(select(CompoundVersion).where(CompoundVersion.compound_row_id == row.id).order_by(CompoundVersion.id.desc()))
            if not version_row:
                continue
            view = build_endpoint_comparison(db, version_row.id)
            coverage = view.get("summary", {}).get("comparison_coverage", {})
            compounds.append({"compound": row.name, "compound_id": row.id, "version_id": version_row.id,
                              "evidence_n": coverage.get("total_evidence", 0),
                              "model_capable_evidence_n": coverage.get("model_capable_evidence", 0),
                              "direct_pairs": coverage.get("direct_pairs", 0), "converted_pairs": coverage.get("converted_pairs", 0),
                              "context_only": coverage.get("context_only", 0), "no_model": coverage.get("no_model", 0),
                              "pairable_evidence_without_pair": coverage.get("pairable_evidence_without_pair", 0)})
    finally:
        db.close()
    no_model_count = sum(int(item["evidence_n"]) for item in gap["no_model_endpoints"])
    action_counts = Counter(item["action"] for item in gap["no_model_endpoints"])
    out = {
        "artifact": "prediction_coverage_gap_closure_v5_1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "golden_compound": {"project_id": 3, "compound": "Sunvozertinib", "version_id": version.id},
        "no_model_inventory": {"reported_no_model_groups": golden["summary"]["comparison_coverage"].get("no_model", 0), "evidence_rows_in_no_model_groups": no_model_count, "actions": dict(action_counts), "rows": gap["no_model_endpoints"]},
        "section_coverage": gap["section_coverage"],
        "qa_metrics": {"EXISTING_CAPABILITY_NOT_EXPOSED": gap["existing_capability_not_exposed"], "PAIRABLE_EVIDENCE_WITHOUT_PAIR": golden["summary"]["comparison_coverage"].get("pairable_evidence_without_pair", 0), "UNKNOWN_MISSING_REASON": gap["unknown_missing_reason"]},
        "high_evidence_compounds": compounds,
        "priority_order": [
            {"endpoint": "PK scenario Cmax/AUC/Tmax/t1/2", "score_basis": "dose/context-qualified evidence plus existing simulation layer", "action": "CONTEXT_SPECIFIC_MODEL_REQUIRED"},
            {"endpoint": "CL/CLF and renal clearance", "score_basis": "PK importance; strict systemic versus apparent semantics", "action": "CONTEXT_SPECIFIC_MODEL_REQUIRED"},
            {"endpoint": "VDss", "score_basis": "distribution importance; species-specific model availability", "action": "NEW_MODEL_REQUIRED"},
            {"endpoint": "pKa/logD7.4", "score_basis": "ADME importance; endpoint-specific validation required", "action": "NEW_MODEL_REQUIRED"},
            {"endpoint": "quantitative CYP/transporters", "score_basis": "measurement semantics and assay context", "action": "SCIENTIFICALLY_NOT_PREDICTABLE"},
        ],
        "scientific_policy": {"observed_target_used_as_prediction_input": False, "new_engine_version": False, "raw_evidence_modified": False, "historical_prediction_runs_modified": False},
    }
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False, default=str) + "\n")
    print(json.dumps({"artifact": str(OUT), "reported_no_model_groups": out["no_model_inventory"]["reported_no_model_groups"], "evidence_rows": no_model_count, "actions": dict(action_counts), "qa": out["qa_metrics"], "compounds": len(compounds)}, indent=2))


if __name__ == "__main__":
    main()
