"""Generate the v3.8B registry and persisted Sunvozertinib match audit."""
from __future__ import annotations

import json
import sys
from pathlib import Path

from sqlalchemy import select

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.canonical_endpoints import registry_report
from backend.database import SessionLocal
from backend.endpoint_comparison import build_endpoint_comparison
from backend.models import Compound


VALIDATION = ROOT / "validation"


def _experimental(row):
    return (
        list(row.get("experimental_internal", []))
        + list(row.get("experimental_external_imported", []))
        + list(row.get("experimental_external_candidates", []))
        + list(row.get("related_evidence", []))
        + list(row.get("needs_review", []))
    )


def _reason(row, has_exp, has_pred):
    comparison = row.get("comparison") or {}
    if comparison.get("status") in {"DIRECT", "CONVERTED"}:
        return ""
    if comparison.get("status") == "RELATED_SAME_SCIENTIFIC_GROUP":
        return "Same scientific group, but measurement semantics do not support numeric error"
    if has_exp and has_pred:
        return "; ".join(item.get("routing_reason", "") for item in _experimental(row) if item.get("routing_reason")) or "Context, unit, or qualification mismatch"
    if has_exp:
        return "No persisted prediction endpoint with the same canonical semantics"
    if has_pred:
        return "No qualified persisted experimental observation for this endpoint"
    return ""


def build_audit(comparison):
    rows = []
    for row in comparison["endpoints"]:
        experiments = _experimental(row)
        prediction = row.get("prediction") or {}
        has_exp = bool(experiments)
        has_pred = bool(prediction.get("available"))
        units = sorted({item.get("normalized_unit", "") for item in experiments if item.get("normalized_unit")})
        prediction_unit = prediction.get("unit") or ""
        comparison_status = (row.get("comparison") or {}).get("status")
        rows.append({
            "canonical_endpoint_id": row["endpoint_id"],
            "canonical_comparison_key": row.get("canonical_comparison_key"),
            "section": row["section"],
            "display_name": row["display_name"],
            "species": row.get("species"),
            "route": row.get("route"),
            "experimental_available": has_exp,
            "prediction_available": has_pred,
            "canonical_match": bool(has_exp and has_pred),
            "unit_match": bool(has_exp and has_pred and (not units or prediction_unit in units)),
            "direct_comparison_possible": comparison_status in {"DIRECT", "CONVERTED"},
            "match_status": comparison_status or ("EXPERIMENTAL_ONLY" if has_exp else "PREDICTION_ONLY"),
            "reason": _reason(row, has_exp, has_pred),
            "prediction": {
                "raw_endpoint": prediction.get("raw_endpoint"),
                "value": prediction.get("display_value"),
                "unit": prediction_unit,
                "prediction_run_id": prediction.get("prediction_run_id"),
            } if has_pred else None,
            "experimental": [{
                "id": item.get("id"),
                "raw_endpoint": item.get("raw_endpoint"),
                "raw_value": item.get("raw_value"),
                "raw_unit": item.get("raw_unit"),
                "normalized_value": item.get("normalized_value"),
                "normalized_unit": item.get("normalized_unit"),
                "origin": item.get("origin"),
                "state": item.get("state"),
                "comparability": item.get("comparability"),
                "importable": item.get("importable"),
                "reference": item.get("reference"),
            } for item in experiments],
            "comparison": row.get("comparison"),
        })
    summary = {}
    for section in ("ACTIVITY", "ADMET", "METABOLISM", "PK", "TOXICITY"):
        subset = [row for row in rows if row["section"] == section]
        summary[section] = {
            "experimental_endpoints": sum(row["experimental_available"] for row in subset),
            "prediction_endpoints": sum(row["prediction_available"] for row in subset),
            "both_scientific_group": sum(row["canonical_match"] for row in subset),
            "directly_comparable": sum(row["direct_comparison_possible"] for row in subset),
            "related_semantic_match": sum(row["match_status"] == "RELATED_SAME_SCIENTIFIC_GROUP" for row in subset),
            "experimental_only": sum(row["match_status"] == "EXPERIMENTAL_ONLY" for row in subset),
            "prediction_only": sum(row["match_status"] == "PREDICTION_ONLY" for row in subset),
        }
    return {
        "artifact": "sunvozertinib_endpoint_match_audit_v3_8b",
        "canonical_endpoint_version": comparison["canonical_endpoint_version"],
        "comparison_unit_version": comparison["comparison_unit_version"],
        "source": "persisted compound endpoint comparison; no new search or prediction",
        "summary": summary,
        "rows": rows,
    }


def main():
    VALIDATION.mkdir(exist_ok=True)
    units = {
        "artifact": "comparison_unit_mapping_v1",
        "comparison_unit_version": "drugopt-comparison-unit-v1",
        "rules": [
            {"endpoint": "HUMAN_PPB / RAT_PPB / MOUSE_PPB", "canonical_unit": "% bound", "conversions": ["fu -> (1-fu)*100", "fraction bound -> fraction*100"], "forbidden": ["ambiguous binding fraction"]},
            {"endpoint": "CACO2_PAPP_AB / CACO2_PAPP_BA", "canonical_unit": "log10(cm/s)", "conversions": ["cm/s -> log10(cm/s)", "10^-6 cm/s -> log10(cm/s)"], "forbidden": ["A->B with B->A", "efflux ratio with Papp"]},
            {"endpoint": "PK Cmax", "canonical_unit": "ng/mL", "conversions": ["µg/L -> ng/mL", "mg/L -> ng/mL"], "forbidden": ["different route/species/analyte without contract"]},
            {"endpoint": "PK AUC", "canonical_unit": "ng*h/mL", "conversions": ["µg*h/L -> ng*h/mL"], "forbidden": ["AUC0-inf with AUCtau"]},
            {"endpoint": "PK CL / CLF", "canonical_unit": "mL/min/kg", "conversions": ["L/h/kg -> mL/min/kg"], "forbidden": ["CL with CL/F"]},
            {"endpoint": "PK Vd / VdF", "canonical_unit": "L/kg", "conversions": ["mL/kg -> L/kg"], "forbidden": ["Vd with Vd/F", "L without body weight -> L/kg"]},
            {"endpoint": "PK Tmax / half-life", "canonical_unit": "hours", "conversions": ["minutes -> hours", "days -> hours"]},
            {"endpoint": "PK F", "canonical_unit": "%", "conversions": ["fraction -> percent"]},
        ],
    }
    (VALIDATION / "comparison_unit_mapping_v1.json").write_text(json.dumps(units, indent=2, ensure_ascii=False) + "\n")
    (VALIDATION / "canonical_endpoint_mapping_v1.json").write_text(json.dumps(registry_report(), indent=2, ensure_ascii=False) + "\n")
    with SessionLocal() as db:
        compound = db.scalar(select(Compound).where(Compound.name.ilike("%sunvozertinib%")))
        if compound is None:
            raise SystemExit("Sunvozertinib compound not found")
        version = next((item for item in compound.versions if item.version_number == compound.current_version), None)
        if version is None:
            raise SystemExit("Sunvozertinib current version not found")
        audit = build_audit(build_endpoint_comparison(db, version.id))
    (VALIDATION / "sunvozertinib_endpoint_match_audit_v3_8b.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False, default=str) + "\n")
    print(json.dumps({"version_id": version.id, "summary": audit["summary"], "rows": len(audit["rows"])}, ensure_ascii=False))


if __name__ == "__main__":
    main()
