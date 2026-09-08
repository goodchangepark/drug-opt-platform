#!/usr/bin/env python3
"""Build the v3.3.3 endpoint inventory from authoritative registries.

This is an audit projection only.  It does not train, promote, or route a
model, and it deliberately preserves endpoint artifact versions independently
from the current engine release identifier.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.endpoint_inventory_audit import ENDPOINT_INVENTORY_47
from backend.prediction_engine_registry import (
    CURRENT_ENGINE_ID,
    CURRENT_POLICY_HASH,
    CURRENT_PRODUCTION_ROUTING,
)
from backend.prediction_maturity import get_endpoint_maturity

OUTPUT = ROOT / "validation" / "v333_endpoint_inventory.json"
BOARD = ROOT / "validation" / "endpoint_completion_board.json"

BOARD_ALIASES = {
    "HUMAN_PPB": "HUMAN_FU_PLASMA",
    "VDSS": "HUMAN_VDSS",
}


def _species_support(endpoint_id: str) -> list[str]:
    if endpoint_id == "RLM_CLINT":
        return ["RAT"]
    if endpoint_id == "MLM_CLINT":
        return ["MOUSE"]
    if endpoint_id == "HLM_CLINT" or endpoint_id.startswith("HUMAN_"):
        return ["HUMAN"]
    if endpoint_id in {"MW", "CLOGP", "TPSA", "HBD", "HBA", "ROTB", "FSP3", "QED", "FORMAL_CHARGE", "HEAVY_ATOM_COUNT", "PKA", "LOGD_7_4"}:
        return ["SPECIES_INDEPENDENT"]
    return ["HUMAN"]


def _scientific_class(record) -> str:
    if record.category == "DETERMINISTIC":
        return "DETERMINISTIC"
    if record.category == "MODEL_UNAVAILABLE":
        return "MODEL_UNAVAILABLE"
    if record.output_type in {"BINARY_CLASSIFICATION", "MULTI_CLASSIFICATION"}:
        return "CLASSIFICATION_ML"
    if "MECHANISTIC" in record.output_type or record.primary_model_family in {"ivive_well_stirred", "tissue_composition"}:
        return "MECHANISTIC"
    return "NUMERICAL_ML"


def build() -> dict:
    board_rows = {}
    if BOARD.exists():
        board = json.loads(BOARD.read_text(encoding="utf-8"))
        board_rows = {row["endpoint"]: row for row in board.get("rows", [])}
    routing = {row["endpoint_id"]: row for row in CURRENT_PRODUCTION_ROUTING}
    rows = []
    for record in ENDPOINT_INVENTORY_47:
        route = routing[record.endpoint_id]
        maturity = get_endpoint_maturity(record.endpoint_id)
        board_row = board_rows.get(BOARD_ALIASES.get(record.endpoint_id, record.endpoint_id), {})
        rows.append({
            "endpoint": record.endpoint_id,
            "display_name": record.display_name,
            "domain": record.domain,
            "classification": _scientific_class(record),
            "output_type": record.output_type,
            "species_support": _species_support(record.endpoint_id),
            "unit": route.get("unit") or record.unit,
            "model": route.get("model_or_ensemble") or record.primary_model_family,
            "model_version": route.get("model_version_hash"),
            "engine_id": CURRENT_ENGINE_ID,
            "current_benchmark": board_row.get("last_experiment"),
            "benchmark_n": board_row.get("development_N", maturity.get("validation_n", 0)),
            "primary_metric": board_row.get("metric"),
            "maturity": maturity,
            "ad_capability": route.get("ad_status") or maturity.get("ad_status"),
            "training_evidence_count": None,
            "validation_evidence_count": maturity.get("validation_n", 0),
            "locked_validation_count": maturity.get("locked_test_n", 0),
            "production_status": route.get("status"),
            "route": route.get("route"),
            "optimization_target": board_row.get("target_value"),
            "next_action": board_row.get("next_action") or route.get("known_limitation"),
            "conversion_guard": record.conversion_guard,
        })
    classes = {}
    for row in rows:
        classes[row["classification"]] = classes.get(row["classification"], 0) + 1
    return {
        "artifact": "V333_ENDPOINT_INVENTORY",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "engine_id": CURRENT_ENGINE_ID,
        "policy_hash": CURRENT_POLICY_HASH,
        "total_endpoints": len(rows),
        "classification_counts": classes,
        "invariants": [
            "Endpoint artifact versions are independent of engine release version.",
            "Classifier probabilities are never converted to quantitative potency.",
            "MODEL_UNAVAILABLE and context-dependent endpoints fail closed.",
            "Maturity derives from the endpoint maturity registry.",
        ],
        "endpoints": rows,
    }


if __name__ == "__main__":
    payload = build()
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUTPUT.relative_to(ROOT)} ({payload['total_endpoints']} endpoints)")
