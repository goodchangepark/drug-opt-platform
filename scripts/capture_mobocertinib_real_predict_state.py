#!/usr/bin/env python3
"""Read-only production Mobocertinib state capture for the P0 repair trace."""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

from sqlalchemy import select

from backend.admet import ADMETPrediction, ADMETPredictionRun, PredictionEndpointSnapshot
from backend.database import SessionLocal
from backend.models import Compound, CompoundVersion, PredictionRun
from backend.stable_core import CurrentPredictionSnapshot, admit_current_prediction


def iso(value):
    return value.isoformat() if isinstance(value, datetime) else value


def capture(output: str, phase: str, response_path: str | None = None) -> None:
    db = SessionLocal()
    try:
        compound = db.scalar(select(Compound).where(Compound.project_id == 3, Compound.name.ilike("Mobocertinib")))
        if compound is None:
            raise RuntimeError("EGFR Mobocertinib was not found")
        version = db.scalar(select(CompoundVersion).where(
            CompoundVersion.compound_row_id == compound.id,
            CompoundVersion.version_number == compound.current_version,
        ))
        if version is None:
            raise RuntimeError("Current Mobocertinib version was not found")

        current = list(db.scalars(select(CurrentPredictionSnapshot).where(
            CurrentPredictionSnapshot.compound_version_id == version.id,
        ).order_by(CurrentPredictionSnapshot.id)))
        endpoint_snapshots = list(db.scalars(select(PredictionEndpointSnapshot).where(
            PredictionEndpointSnapshot.compound_version_id == version.id,
        ).order_by(PredictionEndpointSnapshot.id)))
        prediction_runs = list(db.scalars(select(PredictionRun).where(
            PredictionRun.version_id == version.id,
        ).order_by(PredictionRun.id)))
        admet_runs = list(db.scalars(select(ADMETPredictionRun).where(
            ADMETPredictionRun.version_id == version.id,
        ).order_by(ADMETPredictionRun.id)))
        admet_predictions = list(db.scalars(select(ADMETPrediction).where(
            ADMETPrediction.version_id == version.id,
        ).order_by(ADMETPrediction.id)))

        snapshot_rows = []
        for row in current:
            admission = admit_current_prediction(db, row)
            snapshot_rows.append({
                "id": row.id,
                "endpoint": row.canonical_endpoint,
                "species": row.species,
                "value": row.value,
                "unit": row.unit,
                "classification": row.classification,
                "model_id": row.model_id,
                "model_version": row.model_version,
                "artifact_hash": row.model_artifact_hash,
                "prediction_mode": row.prediction_mode,
                "engine_release": row.engine_release,
                "context": row.context_json,
                "structure_revision": row.structure_revision,
                "is_current": bool(row.is_current),
                "admission": {"eligible": admission.eligible, "reason": admission.reason},
                "created_at": iso(row.created_at),
            })
        endpoint_rows = [{
            "id": row.id, "prediction_run_id": row.prediction_run_id,
            "endpoint": row.endpoint_id, "endpoint_name": row.endpoint_name,
            "base_value": row.base_value, "base_unit": row.base_unit,
            "project_value": row.project_value, "project_unit": row.project_unit,
            "prediction_type": row.prediction_type, "adapter_version": row.adapter_version,
            "maturity_level": row.maturity_level, "maturity_label": row.maturity_label,
            "created_at": iso(row.created_at),
        } for row in endpoint_snapshots]
        run_rows = [{
            "id": row.id, "stage": row.stage, "model_name": row.model_name,
            "model_version": row.model_version, "inputs_hash": row.inputs_hash,
            "outputs": row.outputs_json, "provenance": row.provenance_json,
            "confidence": row.confidence, "created_at": iso(row.created_at),
        } for row in prediction_runs]
        admet_run_rows = [{
            "id": row.id, "status": row.status, "requested_by": row.requested_by,
            "inputs_hash": row.inputs_hash, "message": row.message,
            "started_at": iso(row.started_at), "completed_at": iso(row.completed_at),
        } for row in admet_runs]
        admet_rows = [{
            "id": row.id, "run_id": row.run_id, "endpoint_id": row.endpoint_id,
            "endpoint": row.endpoint.name if row.endpoint else None,
            "version_id": row.version_id, "model_id": row.model_id,
            "model_version": row.model_version, "execution_status": row.execution_status,
            "predicted_value": row.predicted_value, "unit": row.unit,
            "confidence": row.confidence, "applicability_domain": row.applicability_domain,
            "uncertainty": row.uncertainty, "outputs": row.outputs_json,
            "created_at": iso(row.created_at),
        } for row in admet_predictions]
        payload = {
            "capture": "MobocertinibRealPredictState/1",
            "phase": phase,
            "captured_at": datetime.utcnow().isoformat() + "Z",
            "production_database": "drug_opt.db",
            "compound": {
                "project_id": compound.project_id,
                "project_name": compound.project.name if compound.project else "EGFR",
                "compound_row_id": compound.id,
                "compound_id": compound.compound_id,
                "name": compound.name,
                "compound_version_id": version.id,
                "version_number": version.version_number,
                "current_structure_revision": f"v{version.version_number}",
                "canonical_smiles": version.canonical_smiles,
                "inchi": version.inchi,
                "inchikey": version.inchikey,
            },
            "predict_route": {
                "method": "POST",
                "route": f"/api/compounds/{compound.id}/predict-all",
                "payload": {},
                "frontend_source": "frontend/static/app.js:4538",
                "workflow_route_alias": f"/api/compounds/{compound.id}/predict-workflow",
            },
            "current_prediction_snapshots": snapshot_rows,
            "prediction_endpoint_snapshots": endpoint_rows,
            "admet_predictions": admet_rows,
            "admet_prediction_runs": admet_run_rows,
            "prediction_runs": run_rows,
            "counts": {
                "current_prediction_snapshots": len(snapshot_rows),
                "current_prediction_snapshots_current": sum(row["is_current"] for row in snapshot_rows),
                "prediction_endpoint_snapshots": len(endpoint_rows),
                "admet_predictions": len(admet_rows),
                "admet_prediction_runs": len(admet_run_rows),
                "prediction_runs": len(run_rows),
            },
        }
        if response_path:
            payload["predict_response"] = json.loads(Path(response_path).read_text())
        Path(output).write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str) + "\n")
    finally:
        db.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--phase", choices=("before", "after"), required=True)
    parser.add_argument("--response")
    args = parser.parse_args()
    capture(args.output, args.phase, args.response)
