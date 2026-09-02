"""Deterministic experimental representative selection for scientific rows."""
from __future__ import annotations

REPRESENTATIVE_EXPERIMENTAL_VERSION = "drugopt-representative-experimental-v1"

_ORIGIN_PRIORITY = {
    "INTERNAL_EXPERIMENTAL": 0,
    "EXTERNAL_IMPORTED": 1,
    "EXTERNAL_CANDIDATE": 2,
}
_SEMANTIC_PRIORITY = {
    "DIRECTLY_COMPARABLE": 0,
    "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION": 1,
    "QUALIFIED_DIRECT": 0,
    "QUALIFIED_DETERMINISTIC_CONVERSION": 1,
    "RELATED_SAME_SCIENTIFIC_GROUP": 3,
}


def representative_rank(item: dict) -> tuple:
    """Rank without accepting or consulting any predicted numeric value."""
    context = item.get("context") or {}
    qualification = item.get("qualification_details") or {}
    stages = qualification.get("stages") or {}
    origin = item.get("origin") or item.get("state") or "EXTERNAL_CANDIDATE"
    semantic = item.get("comparability") or item.get("qualification") or ""
    complete_context = bool(stages.get("CONTEXT_QUALIFIED")) or all(
        context.get(key) not in (None, "", "UNSPECIFIED", "UNRESOLVED")
        for key in ("species",) if key in context
    )
    source_quality = str(context.get("source_quality") or item.get("source_quality") or "D").upper()
    reference = item.get("reference") or {}
    reference_present = bool(reference.get("reference") or reference.get("url") or reference.get("source_record_id"))
    stable_id = str(item.get("display_evidence_group_id") or item.get("independent_experiment_group_id") or item.get("id") or "")
    return (
        _ORIGIN_PRIORITY.get(origin, 9),
        _SEMANTIC_PRIORITY.get(semantic, 5),
        0 if complete_context else 1,
        source_quality,
        0 if reference_present else 1,
        stable_id,
    )


def select_representative(items: list[dict]) -> tuple[dict | None, str]:
    eligible = [item for item in items if item.get("display", {}).get("value") is not None]
    if not eligible:
        return None, "NO_DISPLAYABLE_NUMERIC_OBSERVATION"
    selected = min(eligible, key=representative_rank)
    origin = selected.get("origin") or selected.get("state") or "EXTERNAL_CANDIDATE"
    return selected, f"{REPRESENTATIVE_EXPERIMENTAL_VERSION}: origin, semantic compatibility, context completeness, source/reference quality, stable identity"
