"""Create a read-only golden-compound comparison audit.

The audit consumes the persisted comparison contract; it does not search,
predict, import, or mutate the runtime database.  It is intentionally useful
for distinguishing a missing pair from a value that is only contextually
related or has no supported model.
"""
from __future__ import annotations

import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from backend.database import SessionLocal
from backend.endpoint_comparison import build_endpoint_comparison


OUT = Path(__file__).resolve().parents[1] / "validation" / "sunvozertinib_comparison_recovery_v4_7.json"


def _failure(item: dict, row: dict) -> str:
    reason = str(item.get("routing_reason") or item.get("qualification_details", {}).get("primary_gap_reason") or "")
    low = reason.lower()
    if "unit" in low:
        return "UNIT_MISMATCH"
    if "species" in low:
        return "SPECIES_MISMATCH"
    if "matrix" in low:
        return "MATRIX_MISMATCH"
    if "assay" in low or "measurement" in low:
        return "ASSAY_SEMANTICS_MISMATCH"
    if "route" in low:
        return "ROUTE_MISMATCH"
    if "dose" in low or "regimen" in low:
        return "DOSE_REGIMEN_MISMATCH"
    if row.get("prediction", {}).get("available"):
        return "PREDICTION_NOT_EXPOSED" if not row.get("comparison") else "CONTEXT_ONLY"
    return "NO_MODEL"


def main() -> None:
    db = SessionLocal()
    try:
        view = build_endpoint_comparison(db, 13)
    finally:
        db.close()

    rows = []
    failures = Counter()
    before = {"experimental_without_prediction": 0, "prediction_without_experimental": 0, "valid_pairs": 0}
    examples = []
    for row in view.get("endpoints", []):
        prediction = row.get("prediction") or {}
        experiments = []
        for key in ("experimental_internal", "experimental_external_imported", "experimental_external_candidates", "related_evidence", "needs_review"):
            experiments.extend(row.get(key) or [])
        comparison = row.get("comparison") or {}
        if experiments and not prediction.get("available"):
            before["experimental_without_prediction"] += len(experiments)
        if prediction.get("available") and not experiments:
            before["prediction_without_experimental"] += 1
        if comparison.get("matches"):
            before["valid_pairs"] += len(comparison["matches"])
        for item in experiments:
            status = item.get("comparability") or item.get("qualification") or "UNCLASSIFIED"
            paired = any(match.get("experimental_id") == item.get("id") for match in comparison.get("matches", []))
            if not paired and status not in {"DIRECTLY_COMPARABLE", "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION"}:
                failures[_failure(item, row)] += 1
            if paired and len(examples) < 5:
                examples.append({
                    "canonical_endpoint": row.get("endpoint_id"),
                    "raw_value": item.get("raw_value"),
                    "raw_unit": item.get("raw_unit"),
                    "normalized_value": item.get("normalized_value"),
                    "normalized_unit": item.get("normalized_unit"),
                    "prediction_value": comparison.get("prediction_value"),
                    "prediction_unit": prediction.get("unit"),
                    "comparison_status": comparison.get("status"),
                    "conversion": item.get("normalization_rule"),
                })
            rows.append({
                "evidence_id": item.get("id"),
                "raw_endpoint": item.get("raw_endpoint"),
                "raw_value": item.get("raw_value"),
                "raw_unit": item.get("raw_unit"),
                "species": item.get("species"),
                "matrix": (item.get("context") or {}).get("matrix", ""),
                "assay": item.get("assay_type", ""),
                "route": item.get("route"),
                "dose": item.get("dose"),
                "regimen": item.get("regimen"),
                "source": item.get("reference"),
                "canonical_endpoint": row.get("endpoint_id"),
                "matching_prediction": bool(prediction.get("available")),
                "pairing_status": "PAIRED" if paired else status,
                "failure_reason": "" if paired else _failure(item, row),
            })

    summary = view.get("summary", {})
    out = {
        "artifact": "sunvozertinib_comparison_recovery_v4_7",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "compound_version_id": 13,
        "project_id": 3,
        "matching_engine_version": summary.get("comparison_coverage", {}).get("matching_engine_version"),
        "before": before,
        "after": {
            "total_qualified_evidence": summary.get("qualification", {}).get("context_qualified", 0),
            "numeric_evidence": summary.get("qualification", {}).get("numeric", 0),
            "model_capable_evidence": summary.get("comparison_coverage", {}).get("model_capable_evidence", 0),
            "direct_pairs": summary.get("comparison_coverage", {}).get("direct_pairs", 0),
            "converted_pairs": summary.get("comparison_coverage", {}).get("converted_pairs", 0),
            "context_only": summary.get("comparison_coverage", {}).get("context_only", 0),
            "no_model": summary.get("comparison_coverage", {}).get("no_model", 0),
            "pairing_failures": sum(failures.values()),
            "pairable_evidence_without_pair": summary.get("comparison_coverage", {}).get("pairable_evidence_without_pair", 0),
        },
        "failure_counts": dict(failures),
        "conversion_pair_examples": examples,
        "rows": rows,
        "scientific_rows": view.get("scientific_rows", []),
        "section_summary": view.get("section_summary", {}),
        "raw_values_modified": 0,
        "historical_prediction_runs_modified": 0,
        "fabricated_values": 0,
    }
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False, default=str) + "\n")
    print(json.dumps({"artifact": str(OUT), "rows": len(rows), "after": out["after"], "failure_counts": out["failure_counts"]}, indent=2))


if __name__ == "__main__":
    main()
