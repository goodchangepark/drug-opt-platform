"""Generate qualification-contract v4 audit artifacts from persisted data."""
from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

from sqlalchemy import select

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.canonical_endpoints import CANONICAL_ENDPOINT_VERSION, COMPARISON_UNIT_VERSION
from backend.database import SessionLocal
from backend.endpoint_comparison import build_endpoint_comparison
from backend.models import Compound, CompoundVersion, ExternalExperimentalEvidence, ExperimentalSearchRun
from backend.qualification_contract import QUALIFICATION_VERSION, qualification_contract_report
from backend.ivive import PKParameterSet


def _json(value):
    return json.loads(json.dumps(value, default=str))


def _items(view):
    keys = ("experimental_internal", "experimental_external_imported", "experimental_external_candidates", "related_evidence", "needs_review")
    return [item for row in view.get("endpoints", []) for key in keys for item in row.get(key, [])]


def _persisted_summary(view, latest_run=None):
    q = dict(view.get("summary", {}).get("qualification", {}))
    if latest_run:
        q["raw_source_records"] = latest_run.raw_count
        q["search_run_unique_records"] = latest_run.unique_count
        q["search_run_endpoint_qualified"] = latest_run.qualified_count
    q["canonical_endpoint_rows"] = len(view.get("endpoints", []))
    return q


def main():
    with SessionLocal() as db:
        views = {}
        for compound in db.scalars(select(Compound)).all():
            version = next((v for v in compound.versions if v.version_number == compound.current_version), None)
            if version:
                views[compound.id] = (compound, version, build_endpoint_comparison(db, version.id))

        sun = next((value for value in views.values() if value[0].name.lower() == "sunvozertinib"), None)
        if sun:
            compound, version, view = sun
            latest = db.scalar(select(ExperimentalSearchRun).where(ExperimentalSearchRun.compound_id == compound.id).order_by(ExperimentalSearchRun.started_at.desc()))
            source = view.get("summary", {}).get("source_qualification", {})
            rows = []
            for row in view.get("endpoints", []):
                q_items = _items({"endpoints": [row]})
                rows.append({
                    "endpoint": row["endpoint_id"], "display_name": row["display_name"], "section": row["section"],
                    "prediction_available": bool(row.get("prediction", {}).get("available")),
                    "prediction_type": (row.get("prediction") or {}).get("source_type"),
                    "experimental_observations": len(q_items),
                    "qualification": [item.get("qualification_details", {}) for item in q_items],
                    "comparison": row.get("comparison"),
                })
            sun_artifact = {
                "version": "v4.0", "qualification_version": QUALIFICATION_VERSION,
                "compound": {"id": compound.id, "version_id": version.id, "name": compound.name},
                "search_run": _json({"id": latest.id, "search_run_id": latest.search_run_id, "raw_count": latest.raw_count, "unique_count": latest.unique_count, "legacy_summary": latest.summary_json}) if latest else None,
                "persisted_requalification": {"global": _persisted_summary(view, latest), "source": source, "endpoint_rows": rows},
                "policy": qualification_contract_report(),
            }
        else:
            sun_artifact = {"version": "v4.0", "qualification_version": QUALIFICATION_VERSION, "status": "SUNVOZERTINIB_NOT_PRESENT"}

        reconciliation = []
        for compound, version, view in views.values():
            latest = db.scalar(select(ExperimentalSearchRun).where(ExperimentalSearchRun.compound_id == compound.id).order_by(ExperimentalSearchRun.started_at.desc()))
            q = _persisted_summary(view, latest)
            source = view.get("summary", {}).get("source_qualification", {})
            reconciliation.append({"compound": compound.name, "version_id": version.id, "global": q, "sources": source,
                                   "invariants": {"unique_le_raw": q.get("unique_scientific_observations", 0) <= q.get("raw_source_records", 0),
                                                   "direct_le_pairable": q.get("direct", 0) <= q.get("prediction_pairable", 0),
                                                   "importable_le_endpoint": q.get("ready_to_import", 0) <= q.get("endpoint_qualified", 0)}})

        pk_rows = []
        for row in db.scalars(select(PKParameterSet)).all():
            route = "ORAL" if str(row.route).upper() == "PO" else str(row.route).upper()
            values = []
            if row.cl_value is not None: values.append({"parameter": "CL/F" if route == "ORAL" else "CL", "value": row.cl_value, "unit": row.cl_unit, "semantics": "apparent oral clearance" if route == "ORAL" else "systemic/route-context clearance"})
            if row.v_value is not None: values.append({"parameter": "Vd/F" if route == "ORAL" else ("Vss" if "VSS" in str(row.v_type or "").upper() else "Vd"), "value": row.v_value, "unit": row.v_unit, "semantics": "apparent oral volume" if route == "ORAL" else ("systemic steady-state volume" if "VSS" in str(row.v_type or "").upper() else "systemic volume")})
            if row.f_predicted is not None: values.append({"parameter": "F", "value": row.f_predicted, "unit": "%", "semantics": "oral bioavailability relative to IV reference", "reference_route": "IV", "canonical_display": "Oral Bioavailability F"})
            if values: pk_rows.append({"pk_parameter_set_id": row.id, "species": row.species, "stored_route": row.route, "dose": row.dose_value, "dose_unit": row.dose_unit, "values": values})

        out = ROOT / "validation"
        out.mkdir(exist_ok=True)
        (out / "qualification_contract_v4.json").write_text(json.dumps({**qualification_contract_report(), "canonical_endpoint_version": CANONICAL_ENDPOINT_VERSION, "comparison_unit_version": COMPARISON_UNIT_VERSION, "count_policy": "global counts are unique scientific observations; source found counts retain raw search records"}, indent=2, sort_keys=True) + "\n")
        (out / "sunvozertinib_qualification_v4.json").write_text(json.dumps(_json(sun_artifact), indent=2, sort_keys=True) + "\n")
        (out / "source_qualification_reconciliation_v4.json").write_text(json.dumps(_json({"version": "v4.0", "compounds": reconciliation}), indent=2, sort_keys=True) + "\n")
        (out / "pk_prediction_semantic_audit_v4.json").write_text(json.dumps(_json({"version": "v4.0", "policy": {"F": "Oral Bioavailability F relative to IV reference; never F/IV", "CL": "systemic clearance and CL/F remain distinct", "Vd": "systemic Vd/Vss and apparent Vd/F remain distinct"}, "rows": pk_rows}), indent=2, sort_keys=True) + "\n")
        print(json.dumps({"compound_views": len(views), "sunvozertinib": sun_artifact.get("persisted_requalification", {}).get("global", {}), "pk_rows": len(pk_rows)}, indent=2, default=str))


if __name__ == "__main__":
    main()
