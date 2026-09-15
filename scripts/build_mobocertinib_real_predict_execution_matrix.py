#!/usr/bin/env python3
"""Build the production Mobocertinib endpoint execution matrix."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from sqlalchemy import select

from backend.database import SessionLocal
from backend.developability_profile import build_developability_profile
from backend.model_artifact_authority import (
    artifact_bundle_sha256,
    model_artifact_registration,
    validation_registration_error,
)
from backend.models import Compound, CompoundVersion
from backend.predict_all_core_contract import PREDICT_ALL_EXECUTION_ENDPOINTS
from backend.prediction_engine_registry import get_current_production_routing
from backend.stable_core import CurrentPredictionSnapshot, admit_current_prediction


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--response", required=True)
    args = parser.parse_args()
    response = json.loads(Path(args.response).read_text())
    response_rows = {row["endpoint"]: row for row in response.get("endpoint_execution", [])}
    routes = {row["endpoint_id"]: row for row in get_current_production_routing()}
    db = SessionLocal()
    try:
        compound = db.scalar(select(Compound).where(Compound.project_id == 3, Compound.name.ilike("Mobocertinib")))
        version = db.scalar(select(CompoundVersion).where(
            CompoundVersion.compound_row_id == compound.id,
            CompoundVersion.version_number == compound.current_version,
        ))
        profile = build_developability_profile(db, version.id)
        snapshots = {
            row.canonical_endpoint: row
            for row in db.scalars(select(CurrentPredictionSnapshot).where(
                CurrentPredictionSnapshot.compound_version_id == version.id,
                CurrentPredictionSnapshot.is_current.is_(True),
            ))
        }
        rows = []
        for profile_row in profile["availability_catalog"]:
            endpoint = profile_row["query_endpoint"]
            snapshot = snapshots.get(endpoint)
            registration = model_artifact_registration(endpoint)
            route = routes.get(endpoint)
            admission = admit_current_prediction(db, snapshot) if snapshot else None
            prediction = profile_row.get("prediction") or {}
            response_row = response_rows.get(endpoint, {})
            status = profile_row.get("status") or profile_row.get("availability")
            if snapshot and admission and admission.eligible:
                terminal = "PUBLISHED"
            elif status == "MECHANISTIC_ONLY":
                terminal = "MECHANISTIC_ONLY"
            elif status == "MODEL_NOT_REGISTERED":
                terminal = "MODEL_NOT_REGISTERED"
            elif status == "MODEL_UNAVAILABLE":
                terminal = "MODEL_UNAVAILABLE"
            elif status == "CURRENT_DATA_CEILING":
                terminal = "CURRENT_DATA_CEILING"
            elif status == "CONTEXT_REQUIRED":
                terminal = "CONTEXT_REQUIRED"
            elif response_row.get("terminal_state") == "ADMISSION_FAILED":
                terminal = "ADMISSION_FAILED"
            elif status == "ON_DEMAND":
                terminal = "NOT_IN_WORKFLOW"
            else:
                terminal = "NOT_IN_WORKFLOW"
            artifact_paths = list(registration.artifact_paths) if registration else []
            artifact_missing = bool(registration and validation_registration_error(registration))
            reason = (
                admission.reason if admission and not admission.eligible else
                profile_row.get("reason") or response_row.get("reason") or status
            )
            rows.append({
                "endpoint": endpoint,
                "display_name": profile_row.get("display_name"),
                "category": profile_row.get("category"),
                "expected_core": endpoint in PREDICT_ALL_EXECUTION_ENDPOINTS,
                "species": profile_row.get("species"),
                "unit": profile_row.get("unit"),
                "model_route_exists": route is not None,
                "model_registry_entry": registration is not None,
                "registered_model_id": snapshot.model_id if snapshot else (registration.model_id if registration else None),
                "model_version": snapshot.model_version if snapshot else (registration.model_version if registration else None),
                "artifact": artifact_paths,
                "artifact_exists": bool(registration) and not artifact_missing,
                "artifact_hash": snapshot.model_artifact_hash if snapshot else (artifact_bundle_sha256(registration) if registration else None),
                "execution_attempted": bool(response_row.get("execution_attempted")) or endpoint in PREDICT_ALL_EXECUTION_ENDPOINTS,
                "execution_succeeded": bool(snapshot and admission and admission.eligible) or response_row.get("execution_succeeded", False),
                "raw_result": prediction.get("value", prediction.get("classification")) if prediction else None,
                "publisher_invoked": bool(response_row.get("publisher_invoked")) or bool(snapshot),
                "stable_core_admission": {
                    "eligible": admission.eligible if admission else False,
                    "reason": admission.reason if admission else "NOT_SUBMITTED",
                },
                "current_snapshot_id": snapshot.id if snapshot else None,
                "canonical_profile_value": prediction.get("value") if prediction else None,
                "canonical_profile_unit": prediction.get("unit") if prediction else profile_row.get("unit"),
                "ui_value": prediction.get("value") if prediction else None,
                "ui_status": status,
                "terminal_state": terminal,
                "failure_reason": reason,
            })
        payload = {
            "matrix": "MobocertinibRealPredictExecutionMatrix/1",
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "production_data_modified_by_explicit_predict": True,
            "compound": {
                "project_id": compound.project_id,
                "compound_row_id": compound.id,
                "compound_version_id": version.id,
                "version_number": version.version_number,
                "inchikey": version.inchikey,
            },
            "predict_response_summary": response.get("summary"),
            "endpoint_count": len(rows),
            "rows": rows,
        }
        Path(args.output).write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str) + "\n")
    finally:
        db.close()


if __name__ == "__main__":
    main()
