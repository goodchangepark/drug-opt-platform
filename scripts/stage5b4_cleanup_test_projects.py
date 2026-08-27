#!/usr/bin/env python3
"""Audit and delete only positively identified development projects."""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

from sqlalchemy import func, select

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.admet import ADMETPrediction
from backend.database import SessionLocal, engine
from backend.human_pk import ensure_human_pk_schema
from backend.main import _confirmed_project_delete
from backend.models import Compound, CompoundVersion, Project
from backend.optimization import OptimizationRun
from backend.pk import PKStudy
from backend.stabilization import classify_project


OUTPUT = ROOT / "validation" / "test_projects_cleanup.json"


def project_rows(db) -> list[dict]:
    rows = []
    for project in db.scalars(select(Project).order_by(Project.id)):
        compounds = int(db.scalar(select(func.count(Compound.id)).where(Compound.project_id == project.id)) or 0)
        pk_studies = int(db.scalar(select(func.count(PKStudy.id)).where(PKStudy.project_id == project.id)) or 0)
        predictions = int(db.scalar(select(func.count(ADMETPrediction.id)).join(
            CompoundVersion, ADMETPrediction.version_id == CompoundVersion.id).join(
            Compound, CompoundVersion.compound_row_id == Compound.id).where(Compound.project_id == project.id)) or 0)
        optimization = int(db.scalar(select(func.count(OptimizationRun.id)).where(OptimizationRun.project_id == project.id)) or 0)
        row = {"project_id": project.id, "project_name": project.name, "target": project.target,
               "description": project.description, "compound_count": compounds, "pk_study_count": pk_studies,
               "prediction_count": predictions, "optimization_run_count": optimization,
               "created_at": project.created_at.isoformat() if project.created_at else None}
        classification, reason = classify_project(row)
        row.update({"classification": classification, "reason_classified_as_test": reason,
                    "deletion_status": "PENDING" if classification == "CONFIRMED_TEST" else "PRESERVED"})
        rows.append(row)
    return rows


def write(payload: dict) -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    # Match production lifespan initialization when this maintenance command is
    # run while the service is stopped.
    ensure_human_pk_schema(engine)
    with SessionLocal() as db:
        rows = project_rows(db)
    confirmed = [row for row in rows if row["classification"] == "CONFIRMED_TEST"]
    payload = {"generated_at": datetime.now(timezone.utc).isoformat(), "projects_before": len(rows),
               "summary": {"CONFIRMED_TEST": len(confirmed),
                           "AMBIGUOUS": sum(row["classification"] == "AMBIGUOUS" for row in rows),
                           "KEEP": sum(row["classification"] == "KEEP" for row in rows)},
               "projects": rows}
    write(payload)  # required machine-readable pre-deletion decision record
    if confirmed:
        confirmations = [{"id": row["project_id"], "confirmation_name": row["project_name"]} for row in confirmed]
        with SessionLocal() as db:
            result = _confirmed_project_delete(db, confirmations)
        deleted = set(result["deleted_project_ids"])
        for row in rows:
            if row["project_id"] in deleted:
                row["deletion_status"] = "DELETED"
    with SessionLocal() as db:
        remaining = set(db.scalars(select(Project.id)))
    for row in rows:
        if row["classification"] != "CONFIRMED_TEST" and row["project_id"] not in remaining:
            raise RuntimeError(f"Preserved project unexpectedly deleted: {row['project_id']}")
    payload["projects_after"] = len(remaining)
    payload["deleted_count"] = sum(row["deletion_status"] == "DELETED" for row in rows)
    payload["preserved_count"] = sum(row["deletion_status"] == "PRESERVED" for row in rows)
    payload["completed_at"] = datetime.now(timezone.utc).isoformat()
    write(payload)
    print(json.dumps({"before": len(rows), "deleted": payload["deleted_count"],
                      "preserved": payload["preserved_count"], "after": len(remaining)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
