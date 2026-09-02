"""Create read-only v4.3 scientific-result audit artifacts from persisted data."""
from __future__ import annotations

import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from sqlalchemy import select

from backend.database import SessionLocal
from backend.endpoint_comparison import build_endpoint_comparison
from backend.models import Compound, ExternalExperimentalEvidence


OUT = ROOT / "validation"
ENGINE = "drugopt-prediction-engine-v1@1.0.0"
ENGINE_HASH = "12757ab197b5a70d8ea1754678d9a342ab0b6ea0d82f2896bebb767d686bbdeb"


def dump(name: str, value: dict) -> None:
    (OUT / name).write_text(json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True, default=str) + "\n", encoding="utf-8")


def observation_bucket(qualification: dict) -> str:
    stages = qualification.get("stages", {})
    if stages.get("DIRECTLY_COMPARABLE"):
        return "DIRECT"
    if stages.get("CONDITIONALLY_COMPARABLE"):
        return "CONVERTED_OR_CONDITIONAL"
    if stages.get("RELATED_SAME_GROUP"):
        return "RELATED"
    if stages.get("CONTEXT_QUALIFIED"):
        return "EXPERIMENTAL_ONLY"
    if stages.get("ENDPOINT_QUALIFIED"):
        return "NEEDS_REVIEW"
    return "NON_COMPARISON_SCIENTIFIC_EVIDENCE"


def main() -> None:
    with SessionLocal() as db:
        compound = db.scalar(select(Compound).where(Compound.name.ilike("%sunvozertinib%")))
        if compound is None:
            raise SystemExit("Sunvozertinib was not found")
        version = next(item for item in compound.versions if item.version_number == compound.current_version)
        comparison = build_endpoint_comparison(db, version.id)
        evidence = db.scalars(select(ExternalExperimentalEvidence).where(ExternalExperimentalEvidence.compound_version_id == version.id).order_by(ExternalExperimentalEvidence.id)).all()
        rows = comparison["scientific_rows"]
        context_sources = Counter()
        remaining_reasons = Counter()
        buckets = Counter()
        for item in evidence:
            context = item.assay_conditions_json or {}
            for key in ("species_source", "route_source", "dose_source", "regimen_source", "analyte_source"):
                context_sources[f"{key}:{context.get(key, 'UNRESOLVED')}"] += 1
            qualification = item.qualification_json or {}
            buckets[observation_bucket(qualification)] += 1
            if qualification.get("primary_gap_reason"):
                remaining_reasons[qualification["primary_gap_reason"]] += 1

        dump("scientific_results_ux_contract_v4_3.json", {
            "contract": "ScientificComparisonRow-v4.3",
            "engine": ENGINE, "engine_hash": ENGINE_HASH,
            "fields": ["section", "group", "canonical_endpoint", "display_name", "species", "route", "dose", "dose_unit", "regimen", "matrix", "assay", "direction", "analyte", "experimental_observations", "primary_experimental_display", "prediction", "display_unit", "difference", "semantic_status", "qualification_status", "prediction_type", "maturity", "references", "unmatched_reason"],
            "rules": [
                "Frontend renders persisted scientific_rows and does not join evidence/prediction/PK streams.",
                "Caco-2 and solubility retain raw log scale while exposing exact scientist-facing display conversions.",
                "CYP/P-gp heterogeneous measurement types are listed separately and never pooled into one numeric range.",
                "PK primary rows require species/context identity; identical systemic Vd/Vss snapshots are collapsed only when values are identical.",
            ],
        })
        dump("sunvozertinib_final_scientific_results_v4_3.json", {
            "generated_at": datetime.now(timezone.utc).isoformat(), "compound": compound.name,
            "compound_id": compound.id, "version_id": version.id, "summary": comparison["summary"],
            "scientific_rows": rows,
        })
        human = [row for row in rows if row["section"] == "PK" and row["species"] == "HUMAN"]
        dump("sunvozertinib_human_clinical_pk_v4_3.json", {
            "generated_at": datetime.now(timezone.utc).isoformat(), "compound": compound.name,
            "section": "HUMAN CLINICAL PK", "rows": human,
            "note": "No prediction is fabricated when no persisted prediction snapshot shares the clinical species/parameter/context.",
        })
        dump("pk_context_qualification_v4_3.json", {
            "generated_at": datetime.now(timezone.utc).isoformat(), "compound": compound.name,
            "context_source_counts": dict(sorted(context_sources.items())),
            "pk_observations": [
                {"id": item.id, "raw_endpoint": item.raw_endpoint_name, "raw_value": item.raw_value, "raw_unit": item.raw_unit,
                 "canonical_endpoint": item.canonical_endpoint_id, "species": (item.assay_conditions_json or {}).get("species", item.species),
                 "route": (item.assay_conditions_json or {}).get("route", "UNSPECIFIED"), "dose": (item.assay_conditions_json or {}).get("dose"),
                 "dose_unit": (item.assay_conditions_json or {}).get("dose_unit", ""), "regimen": (item.assay_conditions_json or {}).get("regimen", "UNSPECIFIED"),
                 "sources": {key: (item.assay_conditions_json or {}).get(key, "UNRESOLVED") for key in ("species_source", "route_source", "dose_source", "regimen_source", "analyte_source")},
                 "qualification": item.qualification_json or {}, "reference": item.reference_text}
                for item in evidence if "_PK_" in str(item.canonical_endpoint_id or "")
            ],
        })
        dump("qualification_completion_v4_3.json", {
            "generated_at": datetime.now(timezone.utc).isoformat(), "compound": compound.name,
            "persistent_observations": len(evidence), "qualification": comparison["summary"]["qualification"],
            "observation_classification_counts": dict(sorted(buckets.items())),
            "needs_review_reasons": dict(sorted(remaining_reasons.items())),
            "complete_definition": "All deterministic table/document context has been resolved; remaining records retain an explicit, persisted reason instead of being forced into a comparison.",
        })


if __name__ == "__main__":
    main()
