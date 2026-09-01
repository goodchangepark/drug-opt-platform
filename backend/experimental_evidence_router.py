"""Single scientific routing decision for external experimental evidence.

The router is deliberately conservative: routing says where an observation
belongs; it never turns related evidence into a prediction target.
"""
from __future__ import annotations

import re

ROUTER_VERSION = "experimental-evidence-routing-v3.2"
SECTIONS = ("ACTIVITY", "ADMET", "METABOLISM", "PK", "TOXICITY", "UNCLASSIFIED")

_PK = re.compile(r"\b(cmax|tmax|auc(?:0[- ]?(?:t|inf)|tau)?|half[- ]life|clearance|cl/f|volume of distribution|vd/f|bioavailability|accumulation|excretion|feces|urine)\b", re.I)
_ACTIVITY = re.compile(r"\b(ic50|ec50|ki|kd|activity|potency|binding)\b", re.I)
_TOX = re.compile(r"\b(herg|ames|dili|cytotox|toxicity|mutagen)\b", re.I)
_METABOLISM = re.compile(r"\b(cyp\w*|cytochrome|metabol|hepatocyte|microsom|clint|mass balance|enzyme involvement)\b", re.I)
_ADMET = re.compile(r"\b(solubility|caco[- ]?2|papp|permeability|ppb|protein binding|fraction unbound|\bfu\b|pka|logd|logp|transporter|p[- ]?gp|bcrp|bsep|oatp|oct|mate)\b", re.I)


def _status(display: dict, row: dict) -> tuple[str, str]:
    value = display.get("comparability_status", "UNSUPPORTED") if display else "UNSUPPORTED"
    return {
        "DIRECTLY_COMPARABLE": ("QUALIFIED_DIRECT", "Directly Comparable"),
        "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION": ("QUALIFIED_DETERMINISTIC_CONVERSION", "Comparable after Conversion"),
        "CONDITIONALLY_COMPARABLE": ("QUALIFIED_CONDITIONAL", "Condition-dependent"),
        "RELATED_NOT_SAME_ENDPOINT": ("QUALIFIED_RELATED", "Related Evidence"),
        "NOT_COMPARABLE": ("NEEDS_REVIEW", display.get("reason") or "Endpoint semantics are not comparable"),
        "UNSUPPORTED": ("UNSUPPORTED", display.get("reason") or "No supported canonical endpoint"),
    }.get(value, ("NEEDS_REVIEW", "Qualification requires review"))


def route_evidence(row: dict, display: dict | None = None) -> dict:
    """Return exactly one primary scientific section for an observation."""
    display = display or row.get("display") or {}
    primary_text = " ".join(str(row.get(k, "")) for k in ("endpoint", "measurement_type", "assay_type")).lower()
    context_text = " ".join(str(row.get(k, "")) for k in ("target", "conditions", "raw_context", "document_type")).lower()
    text = primary_text + " " + context_text
    # A persisted canonical endpoint is authoritative.  In particular,
    # "protein binding" can occur in a narrative containing a target/binding
    # term, but it is ADMET PPB rather than Activity.
    canonical = str(display.get("canonical_endpoint_id") or row.get("canonical_endpoint_id") or "").lower()
    explicit_sections = {
        "ppb_human_percent_bound": "ADMET",
        "solubility_aqueous_logs": "ADMET",
        "permeability_caco2_logpapp": "ADMET",
        "pka": "ADMET",
        "logd_7_4": "ADMET",
        "hlm_intrinsic_clearance_scaled_log10": "METABOLISM",
        "rlm_intrinsic_clearance_scaled_log10": "METABOLISM",
        "mlm_intrinsic_clearance_scaled_log10": "METABOLISM",
    }
    if canonical in explicit_sections:
        section = explicit_sections[canonical]
    elif canonical.startswith(("cyp", "pgp", "bcrp", "bsep", "oatp", "oct", "mate", "metabolite")):
        section = "METABOLISM"
    elif canonical.startswith(("herg", "ames", "dili")):
        section = "TOXICITY"
    elif canonical.startswith("activity:"):
        section = "ACTIVITY"
    elif canonical.startswith(("cmax", "tmax", "auc", "half", "clearance", "cl/", "vd", "bioavailability", "excretion")):
        section = "PK"
    # Prefer the explicitly classified endpoint.  Context may mention a
    # reference value such as Cmax next to a distinct PPB observation.
    elif _PK.search(primary_text):
        section = "PK"
    elif _TOX.search(primary_text):
        section = "TOXICITY"
    elif _ACTIVITY.search(primary_text) and not _METABOLISM.search(primary_text):
        section = "ACTIVITY"
    elif _METABOLISM.search(primary_text):
        section = "METABOLISM"
    elif _ADMET.search(primary_text):
        section = "ADMET"
    elif _PK.search(context_text):
        section = "PK"
    elif _TOX.search(context_text):
        section = "TOXICITY"
    elif _METABOLISM.search(context_text):
        section = "METABOLISM"
    elif _ADMET.search(context_text):
        section = "ADMET"
    else:
        section = "UNCLASSIFIED"

    qualification, label = _status(display, row)
    reference = str(row.get("reference_status", ""))
    if not reference.startswith("REFERENCE_RESOLVED"):
        qualification, label = "NEEDS_REVIEW", "Reference required"
    ready = bool(row.get("import_eligible"))
    return {
        "router_version": ROUTER_VERSION,
        "section": section,
        "canonical_endpoint_id": display.get("canonical_endpoint_id") or row.get("canonical_endpoint_candidate", ""),
        "display_group": str(row.get("endpoint") or row.get("canonical_endpoint_candidate") or "External evidence"),
        "qualification_status": qualification,
        "qualification_label": label,
        "comparability_status": display.get("comparability_status", "UNSUPPORTED"),
        "comparability_label": display.get("comparability_label") or label,
        "routing_reason": f"Matched {section} endpoint semantics" if section != "UNCLASSIFIED" else "No canonical section matched",
        "adaptation_eligibility": ready and qualification in {"QUALIFIED_DIRECT", "QUALIFIED_DETERMINISTIC_CONVERSION"},
        "importability": ready,
    }


def route_records(records: list[dict]) -> list[dict]:
    for row in records:
        row["routing"] = route_evidence(row, row.get("display"))
    return records
