#!/usr/bin/env python3
"""Build the deterministic v3.3.3 release manifest from current registries."""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.prediction_engine_registry import (
    CURRENT_ENGINE_ID,
    CURRENT_ENGINE_VERSION,
    CURRENT_POLICY_HASH,
    CURRENT_PRODUCTION_ROUTING,
    CURRENT_RELEASE_DATE,
    PREVIOUS_ENGINE_ID,
    PREVIOUS_POLICY_HASH,
)

INVENTORY = ROOT / "validation" / "v333_endpoint_inventory.json"
CAMPAIGN = ROOT / "validation" / "v333_completion_campaign.json"
OUTPUT = ROOT / "validation" / "prediction_engine_v3_3_3_manifest.json"


def build() -> dict:
    inventory = json.loads(INVENTORY.read_text(encoding="utf-8"))
    campaign = json.loads(CAMPAIGN.read_text(encoding="utf-8"))
    endpoint_routes = [{
        "endpoint": row["endpoint_id"],
        "route": row["route"],
        "endpoint_model": row["model_or_ensemble"],
        "endpoint_model_version": row["model_version_hash"],
        "routing_status": row["status"],
        "unit": row["unit"],
    } for row in CURRENT_PRODUCTION_ROUTING]
    unavailable = [row["endpoint"] for row in inventory["endpoints"] if row["classification"] == "MODEL_UNAVAILABLE"]
    return {
        "artifact": "PREDICTION_ENGINE_V3_3_3_RELEASE_MANIFEST",
        "engine_id": CURRENT_ENGINE_ID,
        "release_version": CURRENT_ENGINE_VERSION,
        "policy_hash": CURRENT_POLICY_HASH,
        "release_date": CURRENT_RELEASE_DATE,
        "release_status": "PRODUCTION_DEFAULT",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "scientific_scope": "Engine/workflow release; endpoint artifact replacement remains endpoint-gated.",
        "endpoint_inventory_artifact": "validation/v333_endpoint_inventory.json",
        "endpoint_inventory_count": inventory["total_endpoints"],
        "endpoint_routes": endpoint_routes,
        "changes_vs_v3_3_2": [
            "Durable qualified evidence and accepted/candidate separation",
            "Canonical endpoint-semantic comparison and pairability",
            "Bounded reference-library and compound-core navigation",
            "Persisted current prediction snapshot provenance",
            "Explicit species registry and context-aware PK architecture",
        ],
        "improved_endpoint_routes": [],
        "unchanged_endpoint_routes": [row["endpoint_id"] for row in CURRENT_PRODUCTION_ROUTING],
        "research_only_candidates": [
            {
                "endpoint": "HUMAN_FU_PLASMA",
                "artifact": "validation/ppb_canonical_benchmark_v72.json",
                "status": "RESEARCH_BASELINE_RETAINED",
                "rejected_in_v74": ["descriptor_morgan_hybrid", "soft_binding_regime_experts"],
            }
        ],
        "model_unavailable_endpoints": unavailable,
        "rollback": {
            "engine_id": PREVIOUS_ENGINE_ID,
            "policy_hash": PREVIOUS_POLICY_HASH,
            "artifacts_overwritten": False,
            "historical_prediction_runs_mutated": False,
            "procedure": "Select the preserved v3.3.2 registry/route manifest; do not rewrite PredictionRuns.",
        },
        "completion_campaign": {
            "artifact": "validation/v333_completion_campaign.json",
            "status": campaign["status"],
            "exact_next_action": campaign["exact_next_action"],
        },
        "validation": {
            "platform_baseline": "953 passed, 0 failed before v3.3.3 changes",
            "v3_3_3_final_full_pytest": "PENDING",
            "browser_e2e": "PENDING",
            "release_gate": "IN_PROGRESS",
        },
    }


if __name__ == "__main__":
    payload = build()
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
