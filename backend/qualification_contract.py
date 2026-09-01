"""Versioned, stage-specific qualification contract for scientific evidence.

The contract intentionally distinguishes endpoint qualification from the
availability of a Drug-OPT prediction.  A valid experimental observation may
therefore be endpoint-qualified while remaining prediction-unavailable.
"""
from __future__ import annotations

import re
from collections import defaultdict
from typing import Any, Iterable

from .canonical_endpoints import CONVERTED, DIRECT, RELATED, UNSUPPORTED

QUALIFICATION_VERSION = "drugopt-experimental-qualification-v4"

IDENTITY_QUALIFIED = "IDENTITY_QUALIFIED"
REFERENCE_QUALIFIED = "REFERENCE_QUALIFIED"
NUMERIC_QUALIFIED = "NUMERIC_QUALIFIED"
ENDPOINT_QUALIFIED = "ENDPOINT_QUALIFIED"
CONTEXT_QUALIFIED = "CONTEXT_QUALIFIED"
PREDICTION_PAIRABLE = "PREDICTION_PAIRABLE"
DIRECTLY_COMPARABLE = "DIRECTLY_COMPARABLE"
CONDITIONALLY_COMPARABLE = "CONDITIONALLY_COMPARABLE"
RELATED_SAME_GROUP = "RELATED_SAME_GROUP"
IMPORTABLE = "IMPORTABLE"
ADAPTATION_ELIGIBLE = "ADAPTATION_ELIGIBLE"

_NUMBER = re.compile(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


def _bool(value: Any) -> bool:
    return bool(value is True or str(value).strip().lower() in {"true", "yes", "1", "qualified", "resolved"})


def _numeric(record: dict) -> bool:
    if record.get("numeric_observation") is not None:
        return _bool(record.get("numeric_observation"))
    display = record.get("display") or {}
    if display.get("normalized_value") is not None or record.get("normalized_value") is not None:
        return True
    value = record.get("value", record.get("raw_value", ""))
    return bool(_NUMBER.search(str(value or "")))


def _canonical(record: dict) -> str:
    display = record.get("display") or {}
    routing = record.get("routing") or {}
    return str(
        display.get("canonical_endpoint_id")
        or record.get("canonical_endpoint_id")
        or routing.get("canonical_endpoint_id")
        or record.get("canonical_endpoint_candidate")
        or "UNRESOLVED"
    )


def _comparison(record: dict) -> str:
    display = record.get("display") or {}
    routing = record.get("routing") or {}
    return str(display.get("comparability_status") or record.get("comparability_status") or record.get("comparability") or routing.get("comparability_status") or UNSUPPORTED)


def _reference(record: dict) -> bool:
    if str(record.get("origin") or record.get("state") or "") == "INTERNAL_EXPERIMENTAL":
        return True
    status = str(record.get("reference_status") or "")
    return status.startswith("REFERENCE_RESOLVED") or _bool(record.get("reference_qualified"))


def _identity(record: dict) -> bool:
    if str(record.get("origin") or record.get("state") or "") == "INTERNAL_EXPERIMENTAL":
        return True
    status = str(record.get("identity_match_status") or "")
    return status in {"EXACT_STRUCTURE_MATCH", "PUBLIC_IDENTIFIER_RESOLVED", "IDENTITY_QUALIFIED"} or _bool(record.get("identity_qualified"))


def _endpoint(record: dict, numeric: bool) -> bool:
    endpoint = _canonical(record)
    return bool(endpoint and endpoint not in {"UNRESOLVED", "", "NONE"} and (record.get("endpoint_qualified") is not False or numeric))


def _context(record: dict, endpoint: bool, comparison: str) -> tuple[bool, str]:
    if not endpoint:
        return False, "ENDPOINT_UNRESOLVED"
    reason = str((record.get("display") or {}).get("reason") or record.get("routing_reason") or "")
    reason_lower = reason.lower()
    for token, gap in (("direction", "DIRECTION_MISSING"), ("species", "SPECIES_MISSING"), ("route", "ROUTE_MISSING"), ("dose", "DOSE_MISSING"), ("assay", "ASSAY_SYSTEM_MISMATCH"), ("matrix", "MATRIX_MISMATCH"), ("analyte", "ANALYTE_MISMATCH")):
        if token in reason_lower:
            reason = gap
            break
    # PK observations have endpoint-specific context requirements.  Source
    # parsers often leave an explicit UNSPECIFIED sentinel rather than an
    # empty value, so derive an auditable reason from the persisted fields.
    # This does not infer a route or dose; it only explains why comparison is
    # blocked.
    endpoint_id = _canonical(record)
    if endpoint and "_PK_" in endpoint_id:
        route = str(record.get("route") or "").strip().upper()
        if not route or route in {"UNSPECIFIED", "UNKNOWN", "N/A", "NA"}:
            return False, "ROUTE_MISSING"
        if record.get("dose") is None and any(token in endpoint_id for token in ("CMAX", "AUC", "TMAX", "CLF", "VDF", "F_ORAL")):
            return False, "DOSE_MISSING"
        raw_unit = str(record.get("raw_unit") or record.get("unit") or "").strip()
        if not raw_unit:
            return False, "UNIT_UNSUPPORTED"
        normalized = (record.get("display") or {}).get("normalized_value")
        if normalized is None and record.get("normalized_value") is None:
            return False, "UNIT_UNSUPPORTED"
    if comparison in {DIRECT, CONVERTED, RELATED}:
        return True, ""
    if comparison == "CONDITIONALLY_COMPARABLE":
        return False, reason or "CONTEXT_INCOMPLETE"
    for token, gap in (
        ("direction", "DIRECTION_MISSING"),
        ("species", "SPECIES_MISSING"),
        ("route", "ROUTE_MISSING"),
        ("dose", "DOSE_MISSING"),
        ("assay", "ASSAY_SYSTEM_MISMATCH"),
        ("matrix", "MATRIX_MISMATCH"),
        ("analyte", "ANALYTE_MISMATCH"),
    ):
        if token in reason.lower():
            return False, gap
    return False, reason or "CONTEXT_INCOMPLETE"


def qualify_record(record: dict, *, prediction_endpoints: Iterable[str] = (), imported: bool | None = None) -> dict:
    """Return the complete stage state for one raw/display evidence record."""
    prediction_set = {str(x) for x in prediction_endpoints}
    endpoint = _canonical(record)
    numeric = _numeric(record)
    identity = _identity(record)
    reference = _reference(record)
    endpoint_ok = _endpoint(record, numeric)
    comparison = _comparison(record)
    context, context_reason = _context(record, endpoint_ok, comparison)
    pairable = endpoint in prediction_set
    direct = pairable and comparison in {DIRECT, CONVERTED}
    conditional = pairable and comparison == "CONDITIONALLY_COMPARABLE"
    related = pairable and comparison == RELATED
    duplicate = str(record.get("duplicate_status") or "") == "SAME_MEASUREMENT"
    state = str(record.get("evidence_state") or record.get("state") or "")
    is_imported = state == "EXTERNAL_IMPORTED" if imported is None else imported
    policy_importable = record.get("importable") if "importable" in record else None
    importable = bool(policy_importable) if policy_importable is not None else (identity and reference and numeric and endpoint_ok and context and not duplicate and not is_imported)
    importable = importable and identity and reference and numeric and endpoint_ok and context and not duplicate and not is_imported
    adaptation = bool(record.get("adaptation_eligibility")) or (is_imported and direct and bool(record.get("pre_experimental_freeze")) and not duplicate)
    gaps = []
    if not identity: gaps.append("IDENTITY_UNRESOLVED")
    if not reference: gaps.append("REFERENCE_INSUFFICIENT")
    if not numeric: gaps.append("NUMERIC_PARSE")
    if not endpoint_ok: gaps.append("ENDPOINT_UNRESOLVED")
    if not context and context_reason: gaps.append(context_reason)
    if endpoint_ok and not pairable: gaps.append("NO_CURRENT_PREDICTION_ENDPOINT")
    if duplicate: gaps.append("DUPLICATE_DISPLAY_OBSERVATION")
    if pairable and not direct and not conditional and not related: gaps.append("PREDICTION_CONTEXT_INVALID")
    if direct: primary = ""
    elif related: primary = "RELATED_MEASUREMENT_SEMANTICS"
    elif conditional: primary = context_reason or "CONTEXT_INCOMPLETE"
    elif gaps: primary = gaps[0]
    else: primary = ""
    return {
        "identity_status": IDENTITY_QUALIFIED if identity else "IDENTITY_NOT_QUALIFIED",
        "reference_status": REFERENCE_QUALIFIED if reference else "REFERENCE_NOT_QUALIFIED",
        "numeric_status": NUMERIC_QUALIFIED if numeric else "NUMERIC_NOT_QUALIFIED",
        "endpoint_status": ENDPOINT_QUALIFIED if endpoint_ok else "ENDPOINT_NOT_QUALIFIED",
        "context_status": CONTEXT_QUALIFIED if context else "CONTEXT_NOT_QUALIFIED",
        "canonical_endpoint_id": endpoint,
        "prediction_pairability_status": PREDICTION_PAIRABLE if pairable else "PREDICTION_NOT_PAIRABLE",
        "comparability_status": DIRECTLY_COMPARABLE if direct else (CONDITIONALLY_COMPARABLE if conditional else (RELATED_SAME_GROUP if related else comparison)),
        "importability_status": IMPORTABLE if importable else "NOT_IMPORTABLE",
        "adaptation_eligibility": bool(adaptation),
        "primary_gap_reason": primary,
        "secondary_gap_reasons": sorted(set(gaps[1:] if primary in gaps else gaps)),
        "qualification_version": QUALIFICATION_VERSION,
        "stages": {
            IDENTITY_QUALIFIED: identity, REFERENCE_QUALIFIED: reference,
            NUMERIC_QUALIFIED: numeric, ENDPOINT_QUALIFIED: endpoint_ok,
            CONTEXT_QUALIFIED: context, PREDICTION_PAIRABLE: pairable,
            DIRECTLY_COMPARABLE: direct, CONDITIONALLY_COMPARABLE: conditional,
            RELATED_SAME_GROUP: related, IMPORTABLE: importable,
            ADAPTATION_ELIGIBLE: bool(adaptation),
        },
    }


def _source_names(record: dict) -> list[str]:
    provenance = record.get("display_provenance") or []
    names = {str(item.get("source") or "External") for item in provenance if isinstance(item, dict)}
    names.add(str(record.get("source") or record.get("display_source") or "External"))
    return sorted(names)


def aggregate_qualification(records: list[dict], *, prediction_endpoints: Iterable[str] = (), raw_source_counts: dict[str, int] | None = None) -> dict:
    """Aggregate only the v4 stages; source and global counts share this path."""
    qualified = []
    source = defaultdict(lambda: {"found": 0, "unique": 0, "numeric": 0, "identity_qualified": 0, "reference_qualified": 0, "endpoint_qualified": 0, "context_qualified": 0, "prediction_pairable": 0, "direct": 0, "conditional": 0, "related": 0, "ready_to_import": 0, "adaptation_eligible": 0})
    for record in records:
        q = qualify_record(record, prediction_endpoints=prediction_endpoints)
        record["qualification"] = q
        qualified.append(q)
        for name in _source_names(record):
            item = source[name]
            item["unique"] += 1
            item["numeric"] += int(q["stages"][NUMERIC_QUALIFIED])
            item["identity_qualified"] += int(q["stages"][IDENTITY_QUALIFIED])
            item["reference_qualified"] += int(q["stages"][REFERENCE_QUALIFIED])
            item["endpoint_qualified"] += int(q["stages"][ENDPOINT_QUALIFIED])
            item["context_qualified"] += int(q["stages"][CONTEXT_QUALIFIED])
            item["prediction_pairable"] += int(q["stages"][PREDICTION_PAIRABLE])
            item["direct"] += int(q["stages"][DIRECTLY_COMPARABLE])
            item["conditional"] += int(q["stages"][CONDITIONALLY_COMPARABLE])
            item["related"] += int(q["stages"][RELATED_SAME_GROUP])
            item["ready_to_import"] += int(q["stages"][IMPORTABLE])
            item["adaptation_eligible"] += int(q["stages"][ADAPTATION_ELIGIBLE])
    if raw_source_counts:
        for name, count in raw_source_counts.items(): source[name]["found"] = int(count)
    else:
        for name, item in source.items(): item["found"] = item["unique"]
    totals = {
        "numeric": sum(q["stages"][NUMERIC_QUALIFIED] for q in qualified),
        "identity_qualified": sum(q["stages"][IDENTITY_QUALIFIED] for q in qualified),
        "reference_qualified": sum(q["stages"][REFERENCE_QUALIFIED] for q in qualified),
        "endpoint_qualified": sum(q["stages"][ENDPOINT_QUALIFIED] for q in qualified),
        "context_qualified": sum(q["stages"][CONTEXT_QUALIFIED] for q in qualified),
        "prediction_pairable": sum(q["stages"][PREDICTION_PAIRABLE] for q in qualified),
        "direct": sum(q["stages"][DIRECTLY_COMPARABLE] for q in qualified),
        "conditional": sum(q["stages"][CONDITIONALLY_COMPARABLE] for q in qualified),
        "related": sum(q["stages"][RELATED_SAME_GROUP] for q in qualified),
        "ready_to_import": sum(q["stages"][IMPORTABLE] for q in qualified),
        "adaptation_eligible": sum(q["stages"][ADAPTATION_ELIGIBLE] for q in qualified),
    }
    totals["raw_source_records"] = sum(raw_source_counts.values()) if raw_source_counts else len(records)
    totals["unique_scientific_observations"] = len(records)
    totals["manual_review"] = sum(not q["stages"][ENDPOINT_QUALIFIED] or not q["stages"][CONTEXT_QUALIFIED] for q in qualified)
    totals["related_endpoint_groups"] = len({q["canonical_endpoint_id"] for q in qualified if q["stages"][RELATED_SAME_GROUP]})
    totals["qualification_version"] = QUALIFICATION_VERSION
    return {"qualification_version": QUALIFICATION_VERSION, "global": totals, "sources": dict(source), "records": records}


def qualification_contract_report() -> dict:
    return {
        "version": QUALIFICATION_VERSION,
        "stages": [IDENTITY_QUALIFIED, REFERENCE_QUALIFIED, NUMERIC_QUALIFIED, ENDPOINT_QUALIFIED, CONTEXT_QUALIFIED, PREDICTION_PAIRABLE, DIRECTLY_COMPARABLE, CONDITIONALLY_COMPARABLE, RELATED_SAME_GROUP, IMPORTABLE, ADAPTATION_ELIGIBLE],
        "principle": "Endpoint-qualified evidence may be prediction-unavailable; direct comparison requires a persisted prediction and compatible semantics.",
    }
