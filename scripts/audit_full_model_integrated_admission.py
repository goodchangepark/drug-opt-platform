#!/usr/bin/env python3
"""Audit endpoint routing against executable Stable Core admission provenance.

This is a filesystem/registry audit only.  It does not open a database, run a
model, change routing, or generate a prediction.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.current_prediction_publisher import ADMET_CURRENT_ENDPOINTS  # noqa: E402
from backend.model_artifact_authority import (  # noqa: E402
    artifact_bundle_sha256,
    model_artifact_registration,
    validation_registration_error,
)
from backend.prediction_engine_registry import get_current_production_routing  # noqa: E402

LEDGER = ROOT / "validation/full_model_optimization_ledger.json"
OUTPUT = ROOT / "validation/full_model_integrated_admission_audit.json"


def main() -> int:
    ledger = json.loads(LEDGER.read_text(encoding="utf-8"))
    routes = {row["endpoint_id"]: row for row in get_current_production_routing()}
    producer_endpoints = {value[0] for value in ADMET_CURRENT_ENDPOINTS.values()}
    rows = []
    for result in ledger["endpoint_results"]:
        endpoint = result["endpoint"]
        registration = model_artifact_registration(endpoint)
        route = routes.get(endpoint)
        artifact_hash = artifact_bundle_sha256(registration) if registration else None
        validation_error = validation_registration_error(registration) if registration else "MODEL_NOT_REGISTERED"
        weighted_components = {
            key: value for key, value in (route or {}).get("weights", {}).items()
            if isinstance(value, (int, float)) and value > 0
        }
        # A weighted route needs an immutable composition manifest that binds
        # every named component.  The existing per-endpoint registration only
        # binds its listed files; never infer that a single base checkpoint is
        # the advertised multi-model ensemble.
        ensemble_bundle_complete = len(weighted_components) <= 1 or (
            registration is not None
            and len(registration.artifact_paths) >= len(weighted_components)
        )
        producer_status = (
            "CANONICAL_PUBLISHER_CONNECTED"
            if endpoint in producer_endpoints else
            "NOT_CONNECTED_TO_ON_DEMAND_PUBLISHER"
        )
        registry_status = "VERIFIED" if registration and artifact_hash and not validation_error else "UNQUALIFIED"
        if len(weighted_components) > 1 and not ensemble_bundle_complete:
            registry_status = "INCOMPLETE_WEIGHTED_ENSEMBLE_BUNDLE"
        rows.append({
            "endpoint": endpoint,
            "scientific_decision": result["decision"],
            "production_route": (route or {}).get("route"),
            "registered_model_id": registration.model_id if registration else None,
            "registered_model_version": registration.model_version if registration else None,
            "registered_artifacts": list(registration.artifact_paths) if registration else [],
            "artifact_hash": artifact_hash,
            "validation_error": validation_error,
            "weighted_components": weighted_components,
            "ensemble_bundle_complete": ensemble_bundle_complete,
            "registry_status": registry_status,
            "on_demand_producer_status": producer_status,
            "stable_core_publishable_now": bool(
                endpoint in producer_endpoints
                and registry_status == "VERIFIED"
                and result["decision"] == "RETAIN_EXISTING_MODEL"
            ),
        })
    valid_publishers = [row["endpoint"] for row in rows if row["stable_core_publishable_now"]]
    incomplete_ensembles = [
        row["endpoint"] for row in rows
        if row["registry_status"] == "INCOMPLETE_WEIGHTED_ENSEMBLE_BUNDLE"
    ]
    payload = {
        "artifact": "FULL_MODEL_INTEGRATED_ADMISSION_AUDIT_V1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "engine": "drugopt-prediction-engine-v3@3.3.3",
        "scope": "filesystem, endpoint registry, immutable artifacts, on-demand publisher; no runtime DB access",
        "rows": rows,
        "summary": {
            "endpoint_count": len(rows),
            "canonical_on_demand_publishable": valid_publishers,
            "canonical_on_demand_publishable_count": len(valid_publishers),
            "incomplete_weighted_ensemble_bundles": incomplete_ensembles,
            "incomplete_weighted_ensemble_bundle_count": len(incomplete_ensembles),
            "new_models_promoted": 0,
            "engine_decision": "RETAIN_V3_3_3",
        },
        "scientific_interpretation": [
            "A RETAIN_EXISTING_MODEL decision preserves the prior validated route; it does not claim that the on-demand producer currently emits that route.",
            "Only endpoints with an exact runtime producer plus verified artifact registration can create Stable Core CurrentPredictionSnapshots.",
            "Missing canonical publication remains blank/MODEL_UNAVAILABLE and is safer than substituting a base-model or legacy value for an advertised ensemble.",
            "No v3.4.0 candidate is justified because this campaign promoted no new endpoint model.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload["summary"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
