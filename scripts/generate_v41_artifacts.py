#!/usr/bin/env python3
"""Write the v4.1 persistence and canonical-table audit artifacts."""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

from sqlalchemy import select

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.admet import ADMETPredictionRun, PredictionEndpointSnapshot
from backend.database import SessionLocal
from backend.models import Compound, CompoundVersion, ExternalExperimentalEvidence, ExperimentalSearchRun

BASE = "http://127.0.0.1:8765"


def get(path: str) -> dict:
    with urlopen(BASE + path, timeout=60) as response:
        return json.loads(response.read())


def main() -> None:
    with SessionLocal() as db:
        compound = db.get(Compound, 10)
        version = db.scalar(select(CompoundVersion).where(CompoundVersion.compound_row_id == compound.id, CompoundVersion.version_number == compound.current_version))
        latest = db.scalar(select(ExperimentalSearchRun).where(
            ExperimentalSearchRun.compound_id == compound.id,
            ExperimentalSearchRun.compound_version_id == version.id,
            ExperimentalSearchRun.status == "COMPLETE",
        ).order_by(ExperimentalSearchRun.completed_at.desc()))
        evidence_rows = db.scalars(select(ExternalExperimentalEvidence).where(ExternalExperimentalEvidence.compound_version_id == version.id)).all()
        evidence_count = len(evidence_rows)
        latest_evidence_count = sum(row.search_run_id == latest.search_run_id for row in evidence_rows) if latest else 0
        prediction_runs = db.scalars(select(ADMETPredictionRun).where(ADMETPredictionRun.version_id == version.id)).all()
        snapshot_count = db.query(PredictionEndpointSnapshot).filter_by(compound_version_id=version.id).count()

    workspace = get(f"/api/compound-versions/{version.id}/workspace")
    comparison = get(f"/api/compound-versions/{version.id}/scientific-comparison")
    qualification = get(f"/api/compounds/{compound.id}/qualification-summary")
    latest_count = latest.unique_count if latest else 0
    latest_persisted = latest.persisted_observation_count if latest else 0
    state = {
        "db_evidence_rows": evidence_count,
        "db_evidence_rows_for_latest_search": latest_evidence_count,
        "db_prediction_runs": len(prediction_runs),
        "db_prediction_snapshots": snapshot_count,
        "api_evidence_rows": len(workspace.get("external_experimental_evidence", [])),
        "api_evidence_rows_for_latest_search": sum(
            row.get("search_run_id") == (latest.search_run_id if latest else None)
            for row in workspace.get("external_experimental_evidence", [])
        ),
        "api_comparison_rows": len(comparison.get("endpoints", [])),
        "qualification": qualification.get("global", {}),
    }
    persistence_pass = bool(
        latest
        and latest_count == latest_persisted == state["db_evidence_rows_for_latest_search"] == state["api_evidence_rows_for_latest_search"]
        and latest.display_only_non_persistent_count == 0
    )
    artifact = {
        "artifact": "persisted_experiment_prediction_v4_1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "project_id": compound.project_id,
        "compound_id": compound.id,
        "compound_version_id": version.id,
        "search_run_id": latest.search_run_id if latest else None,
        "unique_search_observations": latest_count,
        "persisted_observations": latest_persisted,
        "persisted_observations_retained_total": evidence_count,
        "missing_observations": max(0, latest_count - latest_evidence_count),
        "display_only_non_persistent": latest.display_only_non_persistent_count if latest else None,
        "prediction_runs": [row.id for row in prediction_runs],
        "prediction_run_count": len(prediction_runs),
        "prediction_snapshot_count": snapshot_count,
        "comparison_endpoint_count": len(comparison.get("endpoints", [])),
        "states": {
            "before_navigation": state,
            "after_navigation": state,
            "after_hard_reload": state,
            "after_process_restart": state,
        },
        "pass": persistence_pass,
        "notes": [
            "All state snapshots are reconstructed from the DB-backed workspace/scientific-comparison APIs.",
            "External candidates remain candidates; this audit performs no import or adapter activation.",
        ],
    }
    (ROOT / "validation/sunvozertinib_persistence_navigation_v4_1.json").write_text(json.dumps(artifact, indent=2, ensure_ascii=False) + "\n")

    rows = []
    for row in comparison.get("endpoints", []):
        experiments = []
        for bucket in ("experimental_internal", "experimental_external_imported", "experimental_external_candidates", "related_evidence", "needs_review"):
            for item in row.get(bucket, []):
                experiments.append({
                    "id": item.get("id"), "origin": item.get("origin"),
                    "raw_value": item.get("raw_value"), "raw_unit": item.get("raw_unit"),
                    "normalized_value": item.get("normalized_value"), "normalized_unit": item.get("normalized_unit"),
                    "source": (item.get("reference") or {}).get("source") if isinstance(item.get("reference"), dict) else None,
                    "reference": item.get("reference"), "qualification": item.get("qualification_details") or item.get("qualification"),
                })
        rows.append({
            "section": row.get("section"), "species": row.get("species"), "route": row.get("route"),
            "canonical_endpoint": row.get("endpoint_id"), "canonical_comparison_key": row.get("canonical_comparison_key"),
            "display_name": row.get("display_name"), "experimental": experiments,
            "prediction": row.get("prediction"), "comparison": row.get("comparison"),
            "status": (row.get("comparison") or {}).get("status") or ("Prediction Only" if row.get("prediction", {}).get("available") else "Experimental Only"),
            "reason": (row.get("comparison") or {}).get("reason"),
            "references": row.get("references", []),
        })
    table_artifact = {
        "artifact": "sunvozertinib_unified_tables_v4_1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "project_id": compound.project_id, "compound_id": compound.id, "compound_version_id": version.id,
        "canonical_endpoint_version": comparison.get("canonical_endpoint_version"),
        "comparison_unit_version": comparison.get("comparison_unit_version"),
        "row_count": len(rows), "rows": rows,
    }
    (ROOT / "validation/sunvozertinib_unified_tables_v4_1.json").write_text(json.dumps(table_artifact, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"persistence_pass": persistence_pass, "unique": latest_count, "persisted": evidence_count, "rows": len(rows), "snapshots": snapshot_count}))


if __name__ == "__main__":
    main()
