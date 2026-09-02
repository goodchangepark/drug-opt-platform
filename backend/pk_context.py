"""Conservative PK study-context recovery for persisted public evidence.

PK values in regulatory documents often inherit their study semantics from a
table title, header, footnote, or the immediately surrounding paragraph.  The
harvester keeps that text verbatim; this module derives only deterministic
context and records where every resolved field came from.
"""
from __future__ import annotations

import re
from typing import Any

from .canonical_endpoints import normalize_route, normalize_species


PK_CONTEXT_QUALIFICATION_VERSION = "drugopt-pk-context-v4.3"


def _text(value: Any) -> str:
    if isinstance(value, dict):
        return " ".join(f"{key} {item}" for key, item in value.items() if item is not None)
    return str(value or "")


def _first_match(patterns: list[tuple[str, str]], texts: list[tuple[str, str]]) -> tuple[str, str]:
    for source, text in texts:
        for pattern, value in patterns:
            if re.search(pattern, text, re.I):
                return value, source
    return "", ""


def _dose(texts: list[tuple[str, str]]) -> tuple[float | None, str, str]:
    patterns = (
        r"(?:single |multiple |repeat(?:ed)? |daily |recommended(?: clinical)? |therapeutic )?(?:oral )?(?:dose|dosing|administration|given)\s*(?:of|at)?\s*(\d+(?:\.\d+)?)\s*(mg/kg|mg/day|mg)\b",
        r"(?:single |multiple |repeat(?:ed)? )?(?:oral )?(\d+(?:\.\d+)?)\s*(mg/kg|mg/day|mg)\s*dose\b",
        r"\b(\d+(?:\.\d+)?)\s*(mg/kg|mg/day|mg)\s*(?:qd|once daily|daily)\b",
    )
    for source, text in texts:
        for pattern in patterns:
            match = re.search(pattern, text, re.I)
            if match:
                return float(match.group(1)), match.group(2), source
    return None, "", ""


def _measurement_issue(raw_endpoint: str, raw_unit: str, text: str) -> str:
    endpoint = str(raw_endpoint or "").lower()
    unit = str(raw_unit or "").lower()
    # Percentage changes, CV/RSE values and exposure margins are not absolute
    # parent-drug PK parameters, even when a source sentence happens to name
    # Cmax, AUC or CL/F.
    if any(term in endpoint for term in ("cmax", "auc", "cl/f", "clearance", "half-life", "tmax")):
        if "%" in unit and re.search(r"(?:%\s*(?:cv|rse)|increased|decreased|lower|higher|ratio|slope|margin|confidence interval|ci|variability)", text, re.I):
            return "MEASUREMENT_SEMANTICS_DIFFER"
    if "auc" in endpoint and unit.strip() in {"h", "hour", "hours"}:
        return "UNIT_UNRESOLVED"
    return ""


def resolve_pk_study_context(
    *, raw_endpoint: Any, raw_value: Any, raw_unit: Any, species: Any = "",
    context: Any = None, source_database: Any = "", source_record_id: Any = "",
) -> dict:
    """Return inherited PK context plus auditable field sources.

    Precedence is observation text, row/column header, table title, footnote,
    immediate paragraph, then a document-level regulatory administration
    context.  The last level is limited to explicit FDA/NDA Sunvozertinib
    clinical parent-drug records; it never supplies an analyte, dose, or a
    route for interaction-substrate observations.
    """
    base = dict(context or {}) if isinstance(context, dict) else {"conditions": _text(context)}
    observation = _text(base.get("conditions") or base.get("observation_text") or "")
    texts = [
        ("OBSERVATION_TEXT", observation),
        ("ROW_HEADER", _text(base.get("row_header"))),
        ("COLUMN_HEADER", _text(base.get("column_header"))),
        ("TABLE_TITLE", _text(base.get("table_title"))),
        ("TABLE_FOOTNOTE", _text(base.get("table_footnote"))),
        ("IMMEDIATE_PARAGRAPH", _text(base.get("paragraph") or base.get("section_context"))),
        ("SECTION_HEADER", _text(base.get("section") or base.get("study_description"))),
    ]
    all_text = " ".join(text for _, text in texts if text)
    explicit_species = normalize_species(species, "")
    inferred_species, species_source = _first_match([
        (r"\b(?:human|patients?|healthy volunteers?|nsclc)\b", "HUMAN"),
        (r"\b(?:rats?|wistar|sprague[- ]dawley)\b", "RAT"),
        (r"\b(?:mice|mouse)\b", "MOUSE"),
        (r"\b(?:dogs?|beagle)\b", "DOG"),
        (r"\b(?:monkeys?|cynomolgus)\b", "MONKEY"),
    ], texts)
    resolved_species = explicit_species if explicit_species != "UNSPECIFIED" else (inferred_species or "UNSPECIFIED")
    if explicit_species != "UNSPECIFIED":
        species_source = "OBSERVATION_FIELD"

    explicit_route = normalize_route(base)
    inferred_route, route_source = _first_match([
        (r"\b(?:po|oral|orally|per os)\b", "ORAL"),
        (r"\b(?:iv|intravenous)\b", "IV"),
        (r"\b(?:sc|subcutaneous)\b", "SC"),
    ], texts)
    endpoint_lower = str(raw_endpoint or "").lower()
    if not inferred_route and ("cl/f" in endpoint_lower or "vd/f" in endpoint_lower or "apparent (oral)" in all_text.lower()):
        inferred_route, route_source = "ORAL", "PARAMETER_SEMANTICS"

    # Explicitly distinguish co-administered probe/substrate exposure from
    # parent-drug PK.  Those records are still retained as scientific DDI
    # evidence but cannot become a Sunvozertinib PK row.
    analyte, analyte_source = "PARENT", "DEFAULT_PARENT"
    if re.search(r"\b(?:digoxin|rosuvastatin|substrate)\b", all_text, re.I):
        analyte, analyte_source = "OTHER_ANALYTE", "OBSERVATION_TEXT"

    # FDA's Sunvozertinib label describes an oral product.  This document
    # context is used only for parent human clinical exposure records where
    # the observation itself establishes patients/clinical dosing but omits
    # the route.  It is deliberately not used for animal or DDI observations.
    regulatory_oral = (
        str(source_database or "").startswith("FDA")
        and "NDA219839" in str(source_record_id or "")
        and resolved_species == "HUMAN"
        and analyte == "PARENT"
        and re.search(r"\b(?:patients?|clinical|recommended|therapeutic|nsclc|healthy)\b", all_text, re.I)
    )
    route = explicit_route if explicit_route != "UNSPECIFIED" else inferred_route
    if not route and regulatory_oral:
        route, route_source = "ORAL", "REGULATORY_DOCUMENT_CONTEXT"
    route = route or "UNSPECIFIED"

    dose, dose_unit, dose_source = _dose(texts)
    regimen, regimen_source = _first_match([
        (r"\b(?:steady[ -]?state|css|multiple|repeat(?:ed)?|qd|once daily)\b", "MULTIPLE_DOSE"),
        (r"\bsingle(?: dose| oral administration)?\b", "SINGLE_DOSE"),
    ], texts)
    resolved_unit, unit_source = "", ""
    if "auc" in endpoint_lower and str(raw_unit or "").strip().lower() in {"h", "hr", "hour", "hours"} and re.search(r"(?:h|hr)\s*\*?\s*ng\s*/\s*ml", all_text, re.I):
        resolved_unit, unit_source = "ng*h/mL", "OBSERVATION_TEXT"
    issue = _measurement_issue(str(raw_endpoint), "" if resolved_unit else str(raw_unit), all_text)
    base.update({
        "species": resolved_species,
        "route": route,
        "dose": dose if dose is not None else base.get("dose"),
        "dose_unit": dose_unit or base.get("dose_unit") or base.get("dose_units") or "",
        "regimen": regimen or base.get("regimen") or "UNSPECIFIED",
        "analyte": analyte,
        "pk_context_version": PK_CONTEXT_QUALIFICATION_VERSION,
        "species_source": species_source or "UNRESOLVED",
        "route_source": route_source or "UNRESOLVED",
        "dose_source": dose_source or "UNRESOLVED",
        "regimen_source": regimen_source or "UNRESOLVED",
        "analyte_source": analyte_source,
        "resolved_unit": resolved_unit,
        "unit_source": unit_source or "OBSERVATION_FIELD",
        "measurement_semantics_issue": issue,
    })
    return base
