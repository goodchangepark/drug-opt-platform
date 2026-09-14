"""Execute the exact quantitative routes frozen in the production registry.

This module is deliberately an execution bridge, not a model-development
surface.  It resolves the already-approved component identifiers, runs those
installed implementations, and applies only the frozen non-negative route
weights.  No legacy value store or experimental observation participates.
"""

from __future__ import annotations

from typing import Any

from .candidate_model_registry import CANDIDATE_ADAPTER_SUITE
from .endpoint_contracts import get_endpoint_contract
from .multimodel import (
    ExecutionStatus,
    _ADAPTER_REGISTRY,
    _V2_ADAPTER_REGISTRY,
)
from .prediction_engine_registry import (
    ROUTE_BEST_SINGLE,
    ROUTE_WEIGHTED_ENSEMBLE,
    get_current_production_routing,
)


EXECUTION_CONTRACT = "DrugOPTCurrentProductionRouting/3.3.3"

# Only continuous, structure-based routes are handled here.  Deterministic
# properties and classifier endpoints have their own existing executors.
QUANTITATIVE_ENDPOINT_NAMES: dict[str, str] = {
    "SOLUBILITY_GENERIC": "Solubility",
    "CACO2_PAPP_AB": "Permeability",
    "HUMAN_PPB": "Plasma protein binding",
    "HLM_CLINT": "HLM intrinsic clearance",
    "CYP1A2_INHIBITION": "CYP1A2 quantitative inhibition",
    "CYP2C9_INHIBITION": "CYP2C9 quantitative inhibition",
    "CYP2D6_INHIBITION": "CYP2D6 quantitative inhibition",
    "CYP3A4_INHIBITION": "CYP3A4 quantitative inhibition",
    "HERG_LIABILITY": "hERG liability",
}


def _adapter_index() -> dict[str, Any]:
    adapters = [
        *_ADAPTER_REGISTRY.values(),
        *_V2_ADAPTER_REGISTRY.values(),
        *CANDIDATE_ADAPTER_SUITE,
    ]
    return {adapter.model_id: adapter for adapter in adapters}


def execute_current_production_quantitative(
    canonical_smiles: str,
) -> dict[str, dict[str, Any]]:
    """Run every qualified quantitative route using its frozen components."""
    adapters = _adapter_index()
    outputs: dict[str, dict[str, Any]] = {}
    routes = {
        row["endpoint_id"]: row
        for row in get_current_production_routing()
        if row.get("endpoint_id") in QUANTITATIVE_ENDPOINT_NAMES
        and row.get("route") in {ROUTE_WEIGHTED_ENSEMBLE, ROUTE_BEST_SINGLE}
    }
    for endpoint_id, endpoint_name in QUANTITATIVE_ENDPOINT_NAMES.items():
        route = routes.get(endpoint_id)
        if route is None:
            continue
        contract = get_endpoint_contract(endpoint_name)
        component_rows: list[dict[str, Any]] = []
        weighted_values: list[tuple[float, float]] = []
        for model_id, weight in dict(route.get("weights") or {}).items():
            adapter = adapters.get(model_id)
            if adapter is None:
                component_rows.append({
                    "model_id": model_id,
                    "weight": float(weight),
                    "status": "MODEL_NOT_REGISTERED",
                })
                continue
            available, unavailable_reason = adapter.is_available()
            if not available:
                component_rows.append({
                    "model_id": model_id,
                    "model_version": adapter.model_version,
                    "weight": float(weight),
                    "status": "MODEL_UNAVAILABLE",
                    "reason": unavailable_reason,
                })
                continue
            payload = adapter.execute(canonical_smiles, contract)
            status = (
                payload.execution_status.value
                if hasattr(payload.execution_status, "value")
                else str(payload.execution_status)
            )
            component = {
                "model_id": payload.model_id,
                "model_version": payload.model_version,
                "weight": float(weight),
                "status": status,
                "value": payload.value,
                "unit": payload.canonical_unit,
                "applicability_domain": payload.applicability_domain,
                "confidence": payload.confidence,
                "runtime_ms": payload.runtime_ms,
            }
            if payload.error_message:
                component["reason"] = payload.error_message
            component_rows.append(component)
            if payload.execution_status == ExecutionStatus.SUCCESS and payload.value is not None:
                weighted_values.append((float(weight), float(payload.value)))

        expected = set(dict(route.get("weights") or {}))
        succeeded = {row["model_id"] for row in component_rows if row["status"] == ExecutionStatus.SUCCESS.value}
        if succeeded != expected or not weighted_values:
            outputs[endpoint_id] = {
                "execution_contract": EXECUTION_CONTRACT,
                "execution_status": "EXECUTION_FAILED",
                "model_id": route.get("model_or_ensemble"),
                "model_version": route.get("model_version_hash"),
                "route": route.get("route"),
                "components": component_rows,
                "reason": "Every frozen production-route component must succeed; partial ensembles are not published.",
            }
            continue

        total_weight = sum(weight for weight, _ in weighted_values)
        value = sum(weight * component_value for weight, component_value in weighted_values) / total_weight
        domains = [str(row.get("applicability_domain") or "UNKNOWN") for row in component_rows]
        domain = (
            "OUT_OF_DOMAIN" if "OUT_OF_DOMAIN" in domains
            else "BORDERLINE" if "BORDERLINE" in domains
            else "IN_DOMAIN"
        )
        outputs[endpoint_id] = {
            "execution_contract": EXECUTION_CONTRACT,
            "execution_status": "SUCCESS",
            "production_prediction": round(value, 6),
            "unit": route.get("unit"),
            "model_id": route.get("model_or_ensemble"),
            "model_version": route.get("model_version_hash"),
            "route": route.get("route"),
            "applicability_domain": domain,
            "components": component_rows,
        }
    return outputs
