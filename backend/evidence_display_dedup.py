"""Display-only deduplication for cross-source evidence representations."""
from __future__ import annotations

import hashlib
import json
from collections import OrderedDict


def _text(value):
    if isinstance(value, dict):
        return json.dumps(value, sort_keys=True, ensure_ascii=False)
    return str(value or "").strip().lower()


def display_group_key(row: dict) -> str:
    """Fingerprint the scientific measurement, deliberately excluding source IDs."""
    display = row.get("display") or {}
    context = row.get("conditions") or row.get("raw_context") or ""
    provenance = row.get("provenance") or {}
    parts = (
        row.get("identity_inchikey") or row.get("identity_key") or row.get("compound_id") or "",
        row.get("canonical_endpoint_id") or row.get("routing", {}).get("canonical_endpoint_id") or row.get("endpoint") or "",
        row.get("measurement_type") or row.get("assay_type") or "",
        row.get("target") or "", row.get("assay_id") or row.get("source_assay_id") or "",
        row.get("species") or "", row.get("direction") or "",
        display.get("normalized_value", row.get("normalized_value", row.get("value", ""))),
        display.get("normalized_unit", row.get("normalized_unit", row.get("unit", ""))),
        row.get("relation", row.get("raw_relation", "=")),
        row.get("concentration") or row.get("dose") or context,
        row.get("doi") or provenance.get("doi") or row.get("pmid") or provenance.get("pmid") or row.get("document_id") or "",
        row.get("table") or row.get("section") or provenance.get("table") or provenance.get("section") or "",
    )
    return "DISPLAY-" + hashlib.sha256("|".join(_text(part) for part in parts).encode()).hexdigest()[:24]


def deduplicate_for_display(records: list[dict]) -> tuple[list[dict], int]:
    """Collapse only provenance-equivalent rows and aggregate source provenance."""
    groups = OrderedDict()
    for row in records:
        key = row.get("display_evidence_group_id") or display_group_key(row)
        row["display_evidence_group_id"] = key
        row.setdefault("independent_experiment_group_id", row.get("source_document_id") or row.get("source_record_id") or key)
        if key not in groups:
            item = dict(row)
            item["display_provenance"] = [{
                "source": row.get("source"), "source_record_id": row.get("source_record_id"),
                "source_url": row.get("source_url"), "reference": row.get("reference"),
            }]
            item["display_source_count"] = 1
            groups[key] = item
        else:
            item = groups[key]
            item["display_provenance"].append({
                "source": row.get("source"), "source_record_id": row.get("source_record_id"),
                "source_url": row.get("source_url"), "reference": row.get("reference"),
            })
            item["display_source_count"] += 1
            item["display_sources"] = sorted({str(p.get("source") or "External") for p in item["display_provenance"]})
            if row.get("import_eligible") and not item.get("import_eligible"):
                for field in ("import_eligible", "qualification_state", "reference_status", "identity_match_status", "endpoint_match_status"):
                    if field in row:
                        item[field] = row[field]
    result = list(groups.values())
    for row in result:
        row.setdefault("display_sources", sorted({str(p.get("source") or "External") for p in row["display_provenance"]}))
    return result, max(0, len(records) - len(result))
