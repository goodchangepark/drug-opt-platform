#!/usr/bin/env python3
"""Audit and clean confirmed test projects, retaining only the single most recent test project."""

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


OUTPUT = ROOT / "validation" / "ui_refinement_test_project_cleanup.json"


def project_rows(db) -> list[dict]:
    rows = []
    for project in db.scalars(select(Project).order_by(Project.created_at.desc(), Project.id.desc())):
        compounds = int(db.scalar(select(func.count(Compound.id)).where(Compound.project_id == project.id)) or 0)
        pk_studies = int(db.scalar(select(func.count(PKStudy.id)).where(PKStudy.project_id == project.id)) or 0)
        predictions = int(db.scalar(select(func.count(ADMETPrediction.id)).join(
            CompoundVersion, ADMETPrediction.version_id == CompoundVersion.id).join(
            Compound, CompoundVersion.compound_row_id == Compound.id).where(Compound.project_id == project.id)) or 0)
        optimization = int(db.scalar(select(func.count(OptimizationRun.id)).where(OptimizationRun.project_id == project.id)) or 0)
        row = {
            "project_id": project.id,
            "project_name": project.name,
            "target": project.target,
            "description": project.description,
            "compound_count": compounds,
            "pk_study_count": pk_studies,
            "prediction_count": predictions,
            "optimization_run_count": optimization,
            "created_at": project.created_at.isoformat() if project.created_at else None,
        }
        classification, reason = classify_project(row)
        row["classification"] = classification
        row["reason_classified_as_test"] = reason
        rows.append(row)
    return rows


def write(payload: dict) -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    ensure_human_pk_schema(engine)
    with SessionLocal() as db:
        rows = project_rows(db)

    # Separate confirmed test vs other projects
    confirmed_tests = [r for r in rows if r["classification"] == "CONFIRMED_TEST"]
    
    # Sort confirmed tests by created_at DESC, id DESC
    # The first one is the MOST RECENT test project
    most_recent_test = confirmed_tests[0] if confirmed_tests else None
    
    for r in rows:
        if r["classification"] == "CONFIRMED_TEST":
            if most_recent_test and r["project_id"] == most_recent_test["project_id"]:
                r["action"] = "KEEP_MOST_RECENT_TEST"
                r["deletion_status"] = "PRESERVED"
            else:
                r["action"] = "DELETE_OLD_TEST"
                r["deletion_status"] = "PENDING"
        elif r["classification"] == "KEEP":
            r["action"] = "PRESERVE_REAL"
            r["deletion_status"] = "PRESERVED"
        else:  # AMBIGUOUS
            r["action"] = "PRESERVE_AMBIGUOUS"
            r["deletion_status"] = "PRESERVED"

    to_delete = [r for r in rows if r["action"] == "DELETE_OLD_TEST"]

    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "projects_before": len(rows),
        "most_recent_test_project": {
            "project_id": most_recent_test["project_id"],
            "project_name": most_recent_test["project_name"],
            "created_at": most_recent_test["created_at"],
        } if most_recent_test else None,
        "summary": {
            "CONFIRMED_TEST_TOTAL": len(confirmed_tests),
            "CONFIRMED_TEST_RETAINED": 1 if most_recent_test else 0,
            "CONFIRMED_TEST_DELETED": len(to_delete),
            "PRESERVE_REAL": sum(r["action"] == "PRESERVE_REAL" for r in rows),
            "PRESERVE_AMBIGUOUS": sum(r["action"] == "PRESERVE_AMBIGUOUS" for r in rows),
        },
        "projects": rows,
    }
    write(payload)

    if to_delete:
        confirmations = [{"id": r["project_id"], "confirmation_name": r["project_name"]} for r in to_delete]
        with SessionLocal() as db:
            result = _confirmed_project_delete(db, confirmations)
        deleted_ids = set(result["deleted_project_ids"])
        for r in rows:
            if r["project_id"] in deleted_ids:
                r["deletion_status"] = "DELETED"

    with SessionLocal() as db:
        remaining = set(db.scalars(select(Project.id)))

    # Validation assertions
    for r in rows:
        if r["action"] in ("KEEP_MOST_RECENT_TEST", "PRESERVE_REAL", "PRESERVE_AMBIGUOUS"):
            if r["project_id"] not in remaining:
                raise RuntimeError(f"Preserved project unexpectedly deleted: {r['project_id']} ({r['project_name']})")
        elif r["action"] == "DELETE_OLD_TEST":
            if r["project_id"] in remaining:
                raise RuntimeError(f"Test project failed to delete: {r['project_id']} ({r['project_name']})")

    payload["projects_after"] = len(remaining)
    payload["deleted_count"] = sum(r["deletion_status"] == "DELETED" for r in rows)
    payload["preserved_count"] = sum(r["deletion_status"] == "PRESERVED" for r in rows)
    payload["completed_at"] = datetime.now(timezone.utc).isoformat()
    write(payload)

    print(json.dumps({
        "before": len(rows),
        "deleted": payload["deleted_count"],
        "preserved": payload["preserved_count"],
        "retained_test_project": payload["most_recent_test_project"],
        "after": len(remaining)
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
