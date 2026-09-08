#!/usr/bin/env python3
"""Materialize current v3.3.3 snapshot indexes from persisted predictions.

This is deliberately an indexing pass: it never overwrites an ADMET/PK value
or creates a historical PredictionRun.  Expensive prediction workflows remain
explicit and are not triggered by opening a compound.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from sqlalchemy import func, select

from backend.database import SessionLocal
from backend.endpoint_comparison import ensure_admet_prediction_snapshot_index, persist_pk_prediction_snapshots
from backend.models import Compound, CompoundVersion
from backend.admet import PredictionEndpointSnapshot
from backend.prediction_engine_registry import CURRENT_ENGINE_ID, CURRENT_POLICY_HASH, CURRENT_ENGINE_VERSION

PROJECT_IDS = (1, 3, 5, 300)
OUT = ROOT / "validation/v333_snapshot_generation.json"


def main() -> None:
    with SessionLocal() as db:
        versions = db.scalars(select(CompoundVersion).join(Compound).where(Compound.project_id.in_(PROJECT_IDS))).all()
        before = db.scalar(select(func.count(PredictionEndpointSnapshot.id))) or 0
        # ADMET indexing performs one grouped scan over persisted predictions;
        # invoking it once per version would be O(versions × predictions).
        global_admet = ensure_admet_prediction_snapshot_index(db)
        created_admet = int(global_admet.get("snapshots_created", 0))
        created_pk = 0
        by_project = {project_id: {"versions": 0, "snapshots_created": 0} for project_id in PROJECT_IDS}
        for version in versions:
            pk = persist_pk_prediction_snapshots(db, version.id, reuse_existing=True)
            created = int(pk.get("created", 0))
            created_pk += int(pk.get("created", 0))
            project_id = version.compound.project_id
            by_project[project_id]["versions"] += 1
            by_project[project_id]["snapshots_created"] += created
        db.commit()
        after = db.scalar(select(func.count(PredictionEndpointSnapshot.id))) or 0
    payload = {
        "artifact": "V333_CURRENT_SNAPSHOT_GENERATION",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "engine_id": CURRENT_ENGINE_ID,
        "engine_version": CURRENT_ENGINE_VERSION,
        "policy_hash": CURRENT_POLICY_HASH,
        "projects": list(PROJECT_IDS),
        "versions_examined": len(versions),
        "snapshots_before": before,
        "admet_snapshots_created": created_admet,
        "pk_snapshots_created": created_pk,
        "snapshots_after": after,
        "historical_runs_mutated": False,
        "prediction_values_recomputed": False,
        "on_demand_external_search": False,
        "by_project": by_project,
        "contract": "Current indexes point to persisted values; historical PredictionRuns and values remain immutable.",
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
