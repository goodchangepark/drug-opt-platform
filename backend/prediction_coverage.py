"""Auditable prediction-coverage classification for the comparison layer.

This module does not create predictions or alter matching.  It explains why a
persisted experimental group has no current prediction and reports whether a
persisted PK calculation is represented in the canonical result rows.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from typing import Any, Iterable


PK_PARAMETER_LABELS = {
    "CL", "CLF", "VDF", "VSSF", "VD", "VSS", "F", "CMAX", "TMAX",
    "AUC", "AUC0_T", "AUC0_INF", "AUCTAU", "T_HALF",
}


def _evidence(row: dict) -> list[dict]:
    result = []
    for key in ("experimental_internal", "experimental_external_imported",
                "experimental_external_candidates", "related_evidence", "needs_review"):
        result.extend(row.get(key) or [])
    return result


def _pk_parameter(endpoint_id: str) -> str:
    text = str(endpoint_id or "").upper()
    if "_PK_" not in text:
        return ""
    value = text.split("_PK_", 1)[1]
    value = value.split("_ORAL", 1)[0].split("_IV", 1)[0].split("_SC", 1)[0].split("_IP", 1)[0]
    value = value.rsplit("_DAY", 1)[0]
    return value


def _action(row: dict, evidence: list[dict], calculated_endpoints: set[str]) -> str:
    endpoint = str(row.get("endpoint_id") or "")
    prediction = row.get("prediction") or {}
    if endpoint in calculated_endpoints:
        return "EXISTING_PK_CALCULATION_NOT_EXPOSED"
    if prediction.get("available"):
        return "EXISTING_MODEL_NOT_EXPOSED"
    if not evidence:
        return "OUT_OF_SCOPE"
    if endpoint.startswith("PK_") or "_PK_" in endpoint:
        # PK outputs are scenario-dependent.  A missing compatible scenario
        # is not permission to borrow another dose, route, species, or analyte.
        return "CONTEXT_SPECIFIC_MODEL_REQUIRED"
    if endpoint in {"METABOLITE_OBSERVATION", "EXCRETION_FECAL", "EXCRETION_URINARY",
                    "CYP3A4_METABOLIC_CONTRIBUTION", "CYP2D6_METABOLIC_CONTRIBUTION"}:
        return "SCIENTIFICALLY_NOT_PREDICTABLE"
    return "NEW_MODEL_REQUIRED"


def build_prediction_coverage(view: dict, *, calculated_endpoints: Iterable[str] = ()) -> dict[str, Any]:
    calculated = {str(x) for x in calculated_endpoints}
    no_model = []
    section = defaultdict(lambda: Counter(qualified_numeric=0, model_capable=0, direct=0,
                                           converted=0, context_only=0, no_model=0, total=0))
    unknown = []
    for row in view.get("endpoints", []):
        evidence = _evidence(row)
        prediction = row.get("prediction") or {}
        status = (row.get("comparison") or {}).get("status") or ""
        sec = row.get("section") or "UNSPECIFIED"
        c = section[sec]
        c["total"] += len(evidence)
        for item in evidence:
            if item.get("normalized_value") is not None:
                c["qualified_numeric"] += 1
        if prediction.get("available"):
            c["model_capable"] += len(evidence) or 1
        matches = (row.get("comparison") or {}).get("matches") or []
        direct_matches = sum(1 for match in matches if match.get("status") in {"DIRECT", "DIRECTLY_COMPARABLE"})
        converted_matches = sum(1 for match in matches if match.get("status") in {"CONVERTED", "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION"})
        c["direct"] += direct_matches
        c["converted"] += converted_matches
        if evidence and not prediction.get("available"):
            c["no_model"] += len(evidence)
        elif evidence:
            c["context_only"] += max(0, len(evidence) - direct_matches - converted_matches)
        if evidence and not prediction.get("available"):
            action = _action(row, evidence, calculated)
            if action not in {
                "EXISTING_MODEL_NOT_EXPOSED", "EXISTING_PK_CALCULATION_NOT_EXPOSED",
                "NEW_MODEL_REQUIRED", "CONTEXT_SPECIFIC_MODEL_REQUIRED",
                "SCIENTIFICALLY_NOT_PREDICTABLE", "OUT_OF_SCOPE",
            }:
                unknown.append(row.get("endpoint_id"))
                action = "UNKNOWN_MISSING_REASON"
            no_model.append({
                "canonical_endpoint": row.get("endpoint_id"),
                "evidence_n": len(evidence),
                "species": sorted({x.get("species", "UNSPECIFIED") for x in evidence}),
                "contexts": [{"route": x.get("route"), "dose": x.get("dose"),
                              "dose_unit": x.get("dose_unit"), "regimen": x.get("regimen"),
                              "analyte": x.get("analyte")} for x in evidence[:10]],
                "unit": sorted({x.get("normalized_unit") or x.get("raw_unit") or "UNRESOLVED" for x in evidence}),
                "current_prediction_route": prediction.get("source_type") or "MODEL_UNAVAILABLE",
                "why_no_model": prediction.get("unavailable_reason") or "No compatible persisted prediction",
                "existing_backend_capability": endpoint_capability(row.get("endpoint_id"), calculated),
                "action": action,
            })
    section_out = {}
    for name, values in section.items():
        d = dict(values)
        d["comparison_coverage_percent"] = round(100.0 * (d["direct"] + d["converted"]) / d["total"], 2) if d["total"] else 0.0
        section_out[name] = d
    return {"no_model_endpoints": no_model, "section_coverage": section_out,
            "unknown_missing_reason": len(unknown),
            "existing_capability_not_exposed": sum(1 for x in no_model if x["action"] in {"EXISTING_MODEL_NOT_EXPOSED", "EXISTING_PK_CALCULATION_NOT_EXPOSED"})}


def endpoint_capability(endpoint_id: str, calculated_endpoints: set[str] | None = None) -> str:
    endpoint = str(endpoint_id or "")
    if endpoint in (calculated_endpoints or set()):
        return "PERSISTED_CALCULATION"
    if "_PK_" in endpoint:
        return "PK_SCENARIO_LAYER_AVAILABLE_BUT_CONTEXT_MUST_MATCH"
    return "NO_ACTIVE_COMPATIBLE_PREDICTION_ROUTE"
