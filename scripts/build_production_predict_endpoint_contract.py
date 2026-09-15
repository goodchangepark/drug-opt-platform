#!/usr/bin/env python3
"""Persist the authoritative structure-only Predict contract."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from sqlalchemy import select

from backend.database import SessionLocal
from backend.developability_profile import build_developability_profile
from backend.model_artifact_authority import artifact_bundle_sha256, model_artifact_registration
from backend.models import Compound, CompoundVersion
from backend.predict_all_core_contract import (
    CONTEXT_REQUIRED_ENDPOINTS,
    CURRENT_DATA_CEILING_ENDPOINTS,
    MECHANISTIC_EXECUTION_ENDPOINTS,
    MODEL_NOT_REGISTERED_ENDPOINTS,
    MODEL_UNAVAILABLE_ENDPOINTS,
    PREDICT_ALL_EXECUTION_ENDPOINTS,
)
from backend.prediction_engine_registry import get_current_production_routing


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    db = SessionLocal()
    try:
        compound = db.scalar(select(Compound).where(Compound.project_id == 3, Compound.name.ilike("Mobocertinib")))
        version = db.scalar(select(CompoundVersion).where(
            CompoundVersion.compound_row_id == compound.id,
            CompoundVersion.version_number == compound.current_version,
        ))
        profile = build_developability_profile(db, version.id)
        catalog = {row["query_endpoint"]: row for row in profile["availability_catalog"]}
        routes = {row["endpoint_id"]: row for row in get_current_production_routing()}
        endpoint_ids = sorted(
            PREDICT_ALL_EXECUTION_ENDPOINTS
            | MODEL_NOT_REGISTERED_ENDPOINTS
            | MODEL_UNAVAILABLE_ENDPOINTS
            | CURRENT_DATA_CEILING_ENDPOINTS
            | CONTEXT_REQUIRED_ENDPOINTS
        )
        rows = []
        for endpoint in endpoint_ids:
            row = catalog.get(endpoint, {})
            registration = model_artifact_registration(endpoint)
            route = routes.get(endpoint)
            if endpoint in MECHANISTIC_EXECUTION_ENDPOINTS:
                support = "MECHANISTIC"
            elif endpoint in CONTEXT_REQUIRED_ENDPOINTS:
                support = "CONTEXT_REQUIRED"
            elif endpoint in MODEL_UNAVAILABLE_ENDPOINTS or endpoint in MODEL_NOT_REGISTERED_ENDPOINTS:
                support = "UNAVAILABLE"
            elif endpoint in CURRENT_DATA_CEILING_ENDPOINTS:
                support = "RESEARCH_ONLY"
            else:
                support = "PRODUCTION_EXECUTABLE"
            rows.append({
                "canonical_endpoint": endpoint,
                "display_name": row.get("display_name"),
                "category": row.get("category"),
                "species": row.get("species"),
                "unit": row.get("unit"),
                "support": support,
                "expected_after_predict": endpoint in PREDICT_ALL_EXECUTION_ENDPOINTS,
                "model_route_exists": route is not None,
                "model_id": registration.model_id if registration else None,
                "model_version": registration.model_version if registration else None,
                "artifact": list(registration.artifact_paths) if registration else [],
                "artifact_hash": artifact_bundle_sha256(registration) if registration else None,
                "availability_reason": row.get("reason") or row.get("status"),
            })
        payload = {
            "contract": "DrugOPTProductionPredictEndpointContract/3.3.3-stable-core-v1.2",
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "engine": "drugopt-prediction-engine-v3@3.3.3",
            "stable_core": "v1.2",
            "compound_reference": {"project_id": compound.project_id, "compound_row_id": compound.id, "compound_version_id": version.id},
            "endpoint_count": len(rows),
            "endpoints": rows,
        }
        Path(args.output).write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str) + "\n")
    finally:
        db.close()


if __name__ == "__main__":
    main()
