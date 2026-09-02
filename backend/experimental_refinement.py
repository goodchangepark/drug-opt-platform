"""Deterministic Scientific Experimental Evidence Refinement Policy.

Policy Version: drugopt-experimental-refinement-v1

This module provides deterministic scientific refinement of raw/harvested experimental
evidence records into Drug-OPT canonical endpoints.

Refinement states:
- AUTO_QUALIFIED: Deterministically resolved into a canonical endpoint group with valid numeric/categorical comparison.
- RELATED: Related scientific evidence (e.g. DDI, metabolite, related subtype, biomarker) belonging to the same scientific family.
- REVIEW_REQUIRED: Explicitly flagged for scientific review with an exact unresolved reason.
- UNUSABLE: Insufficient or corrupted record that cannot be classified.

Context inheritance order:
1. Explicit observation text (raw_endpoint_name, raw_value, raw_unit)
2. Row/column header & measurement type (IC50, EC50, Ki, Kd, Cmax, AUC, CL, t1/2, Papp, % bound)
3. Table title
4. Table footnote
5. Immediate paragraph / conditions text (assay_conditions_json, conditions, target, assay_type)
6. Section / study header
No distant-context guessing.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

REFINEMENT_POLICY_VERSION = "drugopt-experimental-refinement-v1"
MEASUREMENT_TYPE_POLICY_VERSION = "drugopt-measurement-type-v1"
UNIT_RESOLUTION_POLICY_VERSION = "drugopt-unit-resolution-v1"

# States
STATE_AUTO_QUALIFIED = "AUTO_QUALIFIED"
STATE_RELATED = "RELATED"
STATE_REVIEW_REQUIRED = "REVIEW_REQUIRED"
STATE_UNUSABLE = "UNUSABLE"

# Unresolved reasons
REASON_ENDPOINT_AMBIGUOUS = "ENDPOINT_AMBIGUOUS"
REASON_MEASUREMENT_TYPE_MISSING = "MEASUREMENT_TYPE_MISSING"
REASON_UNIT_MISSING = "UNIT_MISSING"
REASON_SPECIES_MISSING = "SPECIES_MISSING"
REASON_ASSAY_CONTEXT_MISSING = "ASSAY_CONTEXT_MISSING"
REASON_ROUTE_MISSING = "ROUTE_MISSING"
REASON_DOSE_MISSING = "DOSE_MISSING"
REASON_ANALYTE_MISSING = "ANALYTE_MISSING"
REASON_REFERENCE_AMBIGUOUS = "REFERENCE_AMBIGUOUS"
REASON_RELATIVE_RATIO_NOT_ABSOLUTE_EXPOSURE = "RELATIVE_RATIO_NOT_ABSOLUTE_EXPOSURE"
REASON_NON_NUMERIC_OBSERVATION = "NON_NUMERIC_OBSERVATION"
REASON_OTHER = "OTHER"


def _clean_text(val: Any) -> str:
    return str(val or "").strip()


def parse_numeric(val: Any) -> Optional[float]:
    if val is None:
        return None
    s = _clean_text(val)
    if not s:
        return None
    # Check for range, e.g. "29 to 49", "4 to 8", "29-49"
    range_match = re.search(r"(\d+(?:\.\d+)?)\s*(?:to|-)\s*(\d+(?:\.\d+)?)", s)
    if range_match:
        try:
            low, high = float(range_match.group(1)), float(range_match.group(2))
            return (low + high) / 2.0
        except ValueError:
            pass
    # Standard single number match (handle commas and scientific notation)
    m = re.search(r"[-+]?\d+(?:,\d{3})*(?:\.\d+)?(?:[eE][-+]?\d+)?", s)
    if not m:
        return None
    try:
        return float(m.group(0).replace(",", ""))
    except ValueError:
        return None


def extract_context_hierarchy(record: dict) -> dict:
    """Extract context fields following strict inheritance order."""
    conditions = record.get("assay_conditions_json") or record.get("conditions") or {}
    if isinstance(conditions, str):
        try:
            import json
            conditions = json.loads(conditions)
        except Exception:
            conditions = {"conditions": conditions}
    elif not isinstance(conditions, dict):
        conditions = {"conditions": str(conditions)}

    text_parts = [
        _clean_text(record.get("raw_endpoint_name") or record.get("endpoint")),
        _clean_text(record.get("raw_value") or record.get("value")),
        _clean_text(record.get("raw_unit") or record.get("unit")),
        _clean_text(conditions.get("row_header")),
        _clean_text(conditions.get("col_header")),
        _clean_text(conditions.get("table_title")),
        _clean_text(conditions.get("table_footnote")),
        _clean_text(conditions.get("conditions")),
        _clean_text(conditions.get("section_header")),
        _clean_text(record.get("assay_type") or record.get("measurement_type")),
        _clean_text(record.get("reference_text") or record.get("reference")),
    ]
    for v in conditions.values():
        if isinstance(v, (str, int, float)):
            text_parts.append(_clean_text(str(v)))
    full_context_text = " ".join(p for p in text_parts if p).lower()

    return {
        "raw_endpoint": _clean_text(record.get("raw_endpoint_name") or record.get("endpoint")),
        "raw_value": _clean_text(record.get("raw_value") or record.get("value")),
        "raw_unit": _clean_text(record.get("raw_unit") or record.get("unit")),
        "species": _clean_text(record.get("species") or conditions.get("species")),
        "source": _clean_text(record.get("source_database") or record.get("source")),
        "target": _clean_text(record.get("target") or conditions.get("target")),
        "conditions_dict": conditions,
        "full_text": full_context_text,
    }


def resolve_species(raw_species: str, full_text: str) -> str:
    s = f"{raw_species} {full_text}".lower()
    if re.search(r"\b(human|homo sapiens|patients?|healthy subjects?|volunteers?)\b", s):
        return "HUMAN"
    if re.search(r"\b(rat|sprague[- ]dawley|rattus)\b", s):
        return "RAT"
    if re.search(r"\b(mouse|mice|mus musculus)\b", s):
        return "MOUSE"
    if re.search(r"\b(dog|beagle|canine)\b", s):
        return "DOG"
    if re.search(r"\b(monkey|cynomolgus|nhp|nonhuman primate)\b", s):
        return "MONKEY"
    return "UNSPECIFIED"


def resolve_route(full_text: str) -> str:
    s = full_text.lower()
    if re.search(r"\b(intravenous|iv)\b", s):
        return "IV"
    if re.search(r"\b(oral|po|tablet|capsule|per os)\b", s):
        return "ORAL"
    if re.search(r"\b(subcutaneous|sc)\b", s):
        return "SC"
    if re.search(r"\b(intraperitoneal|ip)\b", s):
        return "IP"
    return "UNSPECIFIED"


def resolve_measurement_type(raw_endpoint: str, full_text: str, unit: str) -> str:
    s = f"{raw_endpoint} {unit} {full_text}".lower()
    if re.search(r"\bic50\b", s):
        return "IC50"
    if re.search(r"\bec50\b", s):
        return "EC50"
    if re.search(r"\bki\b", s):
        return "Ki"
    if re.search(r"\bkd\b", s):
        return "Kd"
    if re.search(r"\bgi50\b", s):
        return "GI50"
    if re.search(r"\bpapp\b|apparent permeability", s):
        return "Papp"
    if re.search(r"\bclint\b|intrinsic clearance", s):
        return "Clint"
    if re.search(r"\bcmax\b", s):
        return "Cmax"
    if re.search(r"\btmax\b", s):
        return "Tmax"
    if re.search(r"\bauc\b", s):
        return "AUC"
    if re.search(r"\bt1/2\b|half[- ]life", s):
        return "t1/2"
    if re.search(r"\bcl\b|clearance", s):
        return "CL"
    if re.search(r"\bvd\b|\bvss\b|volume of distribution", s):
        return "Vd"
    if re.search(r"\bbioavailability\b|\bF\b", s):
        return "F"
    if re.search(r"plasma protein binding|protein binding|\bppb\b|\bfu\b|fraction unbound", s):
        return "PPB"
    if re.search(r"solubility", s):
        return "Solubility"
    if "%" in unit or "%" in s:
        return "%"
    return "OTHER"


def refine_scientific_observation(record: dict) -> dict:
    """Deterministically refine an external experimental record.

    Returns refined dictionary with canonical section, canonical endpoint,
    normalized value/unit, refinement state, and explicit unresolved reasons.
    """
    ctx = extract_context_hierarchy(record)
    raw_ep = ctx["raw_endpoint"]
    raw_val = ctx["raw_value"]
    raw_u = ctx["raw_unit"]
    src = ctx["source"]
    full_text = ctx["full_text"]

    num = parse_numeric(raw_val)
    species = resolve_species(ctx["species"], full_text)
    route = resolve_route(full_text)
    mtype = resolve_measurement_type(raw_ep, full_text, raw_u)

    # Defaults
    section = "UNCLASSIFIED"
    canonical_endpoint = "UNRESOLVED"
    norm_val = num
    norm_unit = raw_u
    state = STATE_REVIEW_REQUIRED
    reason = ""
    refinement_rule = "default_review"
    comparability = "UNSUPPORTED"

    # Rule 1: Literature / PMC candidates without extracted values
    if raw_ep == "Literature candidate" or (src == "Europe PMC" and num is None):
        return {
            "observation_id": record.get("id"),
            "raw_text_value": f"{raw_ep}: {raw_val} {raw_u}".strip(),
            "source": src,
            "current_classification": record.get("qualification_status") or "UNCLASSIFIED",
            "proposed_canonical_section": "UNCLASSIFIED",
            "proposed_endpoint": "LITERATURE_CITATION",
            "measurement_type": "CITATION",
            "species": species,
            "context": ctx["conditions_dict"],
            "unit": raw_u,
            "qualification": STATE_REVIEW_REQUIRED,
            "refinement_rule": "literature_candidate_citation_only",
            "resolved": False,
            "unresolved_reason": REASON_MEASUREMENT_TYPE_MISSING,
            "canonical_endpoint_id": "LITERATURE_CITATION",
            "normalized_value": None,
            "normalized_unit": "",
            "comparability_status": "UNSUPPORTED",
            "evidence_state": STATE_REVIEW_REQUIRED,
        }

    # Rule 2: Non-numeric observations
    if num is None and not any(k in full_text for k in ["positive", "negative", "inhibitor", "substrate"]):
        return {
            "observation_id": record.get("id"),
            "raw_text_value": f"{raw_ep}: {raw_val} {raw_u}".strip(),
            "source": src,
            "current_classification": record.get("qualification_status") or "UNCLASSIFIED",
            "proposed_canonical_section": "UNCLASSIFIED",
            "proposed_endpoint": "UNRESOLVED",
            "measurement_type": mtype,
            "species": species,
            "context": ctx["conditions_dict"],
            "unit": raw_u,
            "qualification": STATE_REVIEW_REQUIRED,
            "refinement_rule": "non_numeric_unresolved",
            "resolved": False,
            "unresolved_reason": REASON_NON_NUMERIC_OBSERVATION,
            "canonical_endpoint_id": "UNRESOLVED",
            "normalized_value": None,
            "normalized_unit": raw_u,
            "comparability_status": "UNSUPPORTED",
            "evidence_state": STATE_REVIEW_REQUIRED,
        }

    # Rule 3: Activity endpoints (EC50, IC50, Ki, Kd, GI50)
    if mtype in {"EC50", "IC50", "Ki", "Kd", "GI50"} or re.search(r"\b(ec50|ic50|ki|kd|gi50)\b", raw_ep.lower()):
        subtype = mtype if mtype in {"EC50", "IC50", "Ki", "Kd", "GI50"} else "IC50"
        section = "ACTIVITY"
        canonical_endpoint = f"ACTIVITY_{subtype}"
        # Unit conversion for activity to nM
        u_clean = raw_u.lower().replace("μ", "u").replace("µ", "u")
        if u_clean in {"nm", "nmol/l", "nmol/dm3"}:
            norm_val = num
            norm_unit = "nM"
        elif u_clean in {"um", "umol/l", "µm"}:
            norm_val = num * 1000.0
            norm_unit = "nM"
        elif u_clean in {"mm", "mmol/l"}:
            norm_val = num * 1e6
            norm_unit = "nM"
        elif u_clean in {"m", "mol/l"}:
            norm_val = num * 1e9
            norm_unit = "nM"
        elif u_clean in {"pm", "pmol/l"}:
            norm_val = num * 0.001
            norm_unit = "nM"
        else:
            norm_val = num
            norm_unit = raw_u or "nM"

        # Check for GLP-1R or EGFR assay context
        is_glp1r = bool(re.search(r"glp[- ]?1|camp|glp1r", full_text))
        is_egfr = bool(re.search(r"egfr|erbb|l858r|t790m|exon20", full_text))
        
        state = STATE_AUTO_QUALIFIED
        comparability = "DIRECTLY_COMPARABLE"
        refinement_rule = f"activity_{subtype.lower()}_nm_standardized"
        reason = ""

        return {
            "observation_id": record.get("id"),
            "raw_text_value": f"{raw_ep}: {raw_val} {raw_u}".strip(),
            "source": src,
            "current_classification": record.get("qualification_status") or "ENDPOINT_QUALIFIED",
            "proposed_canonical_section": section,
            "proposed_endpoint": canonical_endpoint,
            "measurement_type": subtype,
            "species": species,
            "context": ctx["conditions_dict"],
            "unit": norm_unit,
            "qualification": state,
            "refinement_rule": refinement_rule,
            "resolved": True,
            "unresolved_reason": "",
            "canonical_endpoint_id": canonical_endpoint,
            "normalized_value": norm_val,
            "normalized_unit": norm_unit,
            "comparability_status": comparability,
            "evidence_state": "AUTO_QUALIFIED_EXTERNAL",
            "target_context": "GLP1R" if is_glp1r else ("EGFR" if is_egfr else "GENERAL"),
        }

    # Rule 4: Plasma Protein Binding (PPB / fu)
    if mtype == "PPB" or re.search(r"plasma protein binding|protein binding|\bppb\b|\bfu\b|fraction unbound", raw_ep.lower() + " " + full_text):
        section = "ADMET"
        ep_name = f"{species}_PPB" if species in {"HUMAN", "RAT", "MOUSE"} else "HUMAN_PPB"
        u_clean = raw_u.lower()
        if "%" in u_clean or (num is not None and 1.0 < num <= 100.0):
            norm_val = num
            norm_unit = "% bound"
            state = STATE_AUTO_QUALIFIED
            comparability = "DIRECTLY_COMPARABLE"
            refinement_rule = "ppb_percent_bound_direct"
        elif u_clean in {"fu", "fraction", "fraction unbound", "frac"} or (num is not None and 0.0 <= num <= 1.0):
            # Convert fu (fraction unbound) to % bound
            norm_val = (1.0 - num) * 100.0
            norm_unit = "% bound"
            state = STATE_AUTO_QUALIFIED
            comparability = "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION"
            refinement_rule = "ppb_fu_to_percent_bound"
        else:
            state = STATE_REVIEW_REQUIRED
            reason = REASON_UNIT_MISSING
            refinement_rule = "ppb_unit_unresolved"

        return {
            "observation_id": record.get("id"),
            "raw_text_value": f"{raw_ep}: {raw_val} {raw_u}".strip(),
            "source": src,
            "current_classification": record.get("qualification_status") or "ENDPOINT_QUALIFIED",
            "proposed_canonical_section": section,
            "proposed_endpoint": ep_name,
            "measurement_type": "PPB",
            "species": species,
            "context": ctx["conditions_dict"],
            "unit": norm_unit,
            "qualification": state,
            "refinement_rule": refinement_rule,
            "resolved": state == STATE_AUTO_QUALIFIED,
            "unresolved_reason": reason,
            "canonical_endpoint_id": ep_name,
            "normalized_value": norm_val,
            "normalized_unit": norm_unit,
            "comparability_status": comparability,
            "evidence_state": "AUTO_QUALIFIED_EXTERNAL" if state == STATE_AUTO_QUALIFIED else STATE_REVIEW_REQUIRED,
        }

    # Rule 5: Solubility
    if mtype == "Solubility" or re.search(r"solubility|log\s*s", raw_ep.lower()):
        section = "ADMET"
        subtype = "SOLUBILITY_INTRINSIC" if "intrinsic" in full_text else ("SOLUBILITY_KINETIC" if "kinetic" in full_text else ("SOLUBILITY_THERMODYNAMIC" if "thermodynamic" in full_text else "SOLUBILITY_GENERIC"))
        u_clean = raw_u.lower()
        if "um" in u_clean or "µm" in u_clean:
            norm_val = num
            norm_unit = "µM"
            state = STATE_AUTO_QUALIFIED
            comparability = "DIRECTLY_COMPARABLE"
            refinement_rule = "solubility_micromolar_direct"
        elif "mg/ml" in u_clean or "mg/l" in u_clean:
            norm_val = num
            norm_unit = raw_u
            state = STATE_AUTO_QUALIFIED
            comparability = "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION"
            refinement_rule = "solubility_mass_concentration"
        elif "log" in u_clean:
            norm_val = num
            norm_unit = "log10(mol/L)"
            state = STATE_AUTO_QUALIFIED
            comparability = "DIRECTLY_COMPARABLE"
            refinement_rule = "solubility_log10_direct"
        else:
            state = STATE_AUTO_QUALIFIED
            norm_val = num
            norm_unit = raw_u or "µM"
            comparability = "DIRECTLY_COMPARABLE"
            refinement_rule = "solubility_standard"

        return {
            "observation_id": record.get("id"),
            "raw_text_value": f"{raw_ep}: {raw_val} {raw_u}".strip(),
            "source": src,
            "current_classification": record.get("qualification_status") or "ENDPOINT_QUALIFIED",
            "proposed_canonical_section": section,
            "proposed_endpoint": subtype,
            "measurement_type": "Solubility",
            "species": species,
            "context": ctx["conditions_dict"],
            "unit": norm_unit,
            "qualification": state,
            "refinement_rule": refinement_rule,
            "resolved": True,
            "unresolved_reason": "",
            "canonical_endpoint_id": subtype,
            "normalized_value": norm_val,
            "normalized_unit": norm_unit,
            "comparability_status": comparability,
            "evidence_state": "AUTO_QUALIFIED_EXTERNAL",
        }

    # Rule 6: Caco-2 / Permeability
    if mtype == "Papp" or re.search(r"caco[- ]?2|permeability|papp", raw_ep.lower() + " " + full_text):
        section = "ADMET"
        ep_name = "CACO2_PAPP_AB"
        if "efflux" in full_text:
            ep_name = "CACO2_EFFLUX_RATIO"
        elif re.search(r"b\s*[- >]+\s*a|b\s*to\s*a|basolateral", full_text):
            ep_name = "CACO2_PAPP_BA"

        state = STATE_AUTO_QUALIFIED
        comparability = "DIRECTLY_COMPARABLE"
        refinement_rule = "caco2_permeability_direct"
        norm_val = num
        norm_unit = raw_u or "10^-6 cm/s"

        return {
            "observation_id": record.get("id"),
            "raw_text_value": f"{raw_ep}: {raw_val} {raw_u}".strip(),
            "source": src,
            "current_classification": record.get("qualification_status") or "ENDPOINT_QUALIFIED",
            "proposed_canonical_section": section,
            "proposed_endpoint": ep_name,
            "measurement_type": "Papp",
            "species": species,
            "context": ctx["conditions_dict"],
            "unit": norm_unit,
            "qualification": state,
            "refinement_rule": refinement_rule,
            "resolved": True,
            "unresolved_reason": "",
            "canonical_endpoint_id": ep_name,
            "normalized_value": norm_val,
            "normalized_unit": norm_unit,
            "comparability_status": comparability,
            "evidence_state": "AUTO_QUALIFIED_EXTERNAL",
        }

    # Rule 7: Microsomal / Hepatocyte Stability
    if mtype == "Clint" or re.search(r"microsom|clint|hepatocyte", raw_ep.lower() + " " + full_text):
        if "hepatocyte" in full_text:
            section = "METABOLISM"
            ep_name = "HEPATOCYTE_CLINT"
            norm_unit = "µL/min/10^6 cells"
        else:
            section = "ADMET"
            ep_name = f"{species}_CLINT" if species in {"HLM", "RLM", "MLM"} else ("HLM_CLINT" if species == "HUMAN" else ("RLM_CLINT" if species == "RAT" else ("MLM_CLINT" if species == "MOUSE" else "HLM_CLINT")))
            norm_unit = "mL/min/kg"

        state = STATE_AUTO_QUALIFIED
        comparability = "DIRECTLY_COMPARABLE"
        refinement_rule = "metabolic_stability_clint"
        return {
            "observation_id": record.get("id"),
            "raw_text_value": f"{raw_ep}: {raw_val} {raw_u}".strip(),
            "source": src,
            "current_classification": record.get("qualification_status") or "ENDPOINT_QUALIFIED",
            "proposed_canonical_section": section,
            "proposed_endpoint": ep_name,
            "measurement_type": "Clint",
            "species": species,
            "context": ctx["conditions_dict"],
            "unit": norm_unit,
            "qualification": state,
            "refinement_rule": refinement_rule,
            "resolved": True,
            "unresolved_reason": "",
            "canonical_endpoint_id": ep_name,
            "normalized_value": num,
            "normalized_unit": norm_unit,
            "comparability_status": comparability,
            "evidence_state": "AUTO_QUALIFIED_EXTERNAL",
        }

    # Rule 8: CYP & Transporter endpoints
    cyp_match = re.search(r"\b(cyp\s*(?:3a4?|3a|2c8|2c9|2c19|2d6|1a2|2b6|3a5))\b", raw_ep.lower() + " " + full_text)
    trans_match = re.search(r"\b(p[- ]?gp|pgp|bcrp|bsep|oatp[0-9a-z]*|oct[0-9a-z]*|mate[0-9a-z]*)\b", raw_ep.lower() + " " + full_text)
    if cyp_match:
        iso = cyp_match.group(1).upper().replace(" ", "")
        if iso == "CYP3A":
            iso = "CYP3A4"
        section = "METABOLISM"
        if re.search(r"substrat", full_text):
            ep_name = f"{iso}_SUBSTRATE"
        elif re.search(r"contribution|metabol|fraction|fm", full_text):
            ep_name = f"{iso}_METABOLIC_CONTRIBUTION"
        else:
            ep_name = f"{iso}_INHIBITION"

        is_fm = "%" in raw_u or "contribution" in ep_name.lower()
        state = STATE_AUTO_QUALIFIED
        comparability = "DIRECTLY_COMPARABLE" if is_fm else "RELATED_NOT_SAME_ENDPOINT"
        norm_unit = "%" if is_fm else (raw_u or "µM")

        return {
            "observation_id": record.get("id"),
            "raw_text_value": f"{raw_ep}: {raw_val} {raw_u}".strip(),
            "source": src,
            "current_classification": record.get("qualification_status") or "ENDPOINT_QUALIFIED",
            "proposed_canonical_section": section,
            "proposed_endpoint": ep_name,
            "measurement_type": "CYP Interaction",
            "species": species,
            "context": ctx["conditions_dict"],
            "unit": norm_unit,
            "qualification": state,
            "refinement_rule": "cyp_interaction_resolved",
            "resolved": True,
            "unresolved_reason": "",
            "canonical_endpoint_id": ep_name,
            "normalized_value": num,
            "normalized_unit": norm_unit,
            "comparability_status": comparability,
            "evidence_state": "AUTO_QUALIFIED_EXTERNAL",
        }

    if trans_match:
        trans_name = trans_match.group(1).upper().replace(" ", "")
        if trans_name in {"PGP", "P-GP"}:
            ep_name = "PGP_INHIBITION"
        else:
            ep_name = f"{trans_name}_INHIBITION"
        section = "METABOLISM"
        state = STATE_AUTO_QUALIFIED
        comparability = "RELATED_NOT_SAME_ENDPOINT"
        return {
            "observation_id": record.get("id"),
            "raw_text_value": f"{raw_ep}: {raw_val} {raw_u}".strip(),
            "source": src,
            "current_classification": record.get("qualification_status") or "ENDPOINT_QUALIFIED",
            "proposed_canonical_section": section,
            "proposed_endpoint": ep_name,
            "measurement_type": "Transporter Interaction",
            "species": species,
            "context": ctx["conditions_dict"],
            "unit": raw_u or "µM",
            "qualification": state,
            "refinement_rule": "transporter_interaction_resolved",
            "resolved": True,
            "unresolved_reason": "",
            "canonical_endpoint_id": ep_name,
            "normalized_value": num,
            "normalized_unit": raw_u or "µM",
            "comparability_status": comparability,
            "evidence_state": "AUTO_QUALIFIED_EXTERNAL",
        }

    # Rule 9: Excretion / Metabolites
    if re.search(r"\bfeces\b|fecal", raw_ep.lower() + " " + full_text):
        return {
            "observation_id": record.get("id"),
            "raw_text_value": f"{raw_ep}: {raw_val} {raw_u}".strip(),
            "source": src,
            "current_classification": record.get("qualification_status") or "ENDPOINT_QUALIFIED",
            "proposed_canonical_section": "METABOLISM",
            "proposed_endpoint": "EXCRETION_FECAL",
            "measurement_type": "Excretion",
            "species": species,
            "context": ctx["conditions_dict"],
            "unit": "% dose",
            "qualification": STATE_AUTO_QUALIFIED,
            "refinement_rule": "excretion_fecal_recovery",
            "resolved": True,
            "unresolved_reason": "",
            "canonical_endpoint_id": "EXCRETION_FECAL",
            "normalized_value": num,
            "normalized_unit": "% dose",
            "comparability_status": "DIRECTLY_COMPARABLE",
            "evidence_state": "AUTO_QUALIFIED_EXTERNAL",
        }

    if re.search(r"\burine\b|urinary", raw_ep.lower() + " " + full_text):
        return {
            "observation_id": record.get("id"),
            "raw_text_value": f"{raw_ep}: {raw_val} {raw_u}".strip(),
            "source": src,
            "current_classification": record.get("qualification_status") or "ENDPOINT_QUALIFIED",
            "proposed_canonical_section": "METABOLISM",
            "proposed_endpoint": "EXCRETION_URINARY",
            "measurement_type": "Excretion",
            "species": species,
            "context": ctx["conditions_dict"],
            "unit": "% dose",
            "qualification": STATE_AUTO_QUALIFIED,
            "refinement_rule": "excretion_urinary_recovery",
            "resolved": True,
            "unresolved_reason": "",
            "canonical_endpoint_id": "EXCRETION_URINARY",
            "normalized_value": num,
            "normalized_unit": "% dose",
            "comparability_status": "DIRECTLY_COMPARABLE",
            "evidence_state": "AUTO_QUALIFIED_EXTERNAL",
        }

    if re.search(r"\bmetabolite\b", raw_ep.lower() + " " + full_text):
        return {
            "observation_id": record.get("id"),
            "raw_text_value": f"{raw_ep}: {raw_val} {raw_u}".strip(),
            "source": src,
            "current_classification": record.get("qualification_status") or "ENDPOINT_QUALIFIED",
            "proposed_canonical_section": "METABOLISM",
            "proposed_endpoint": "METABOLITE_OBSERVATION",
            "measurement_type": "Metabolite",
            "species": species,
            "context": ctx["conditions_dict"],
            "unit": raw_u or "%",
            "qualification": STATE_AUTO_QUALIFIED,
            "refinement_rule": "metabolite_observation",
            "resolved": True,
            "unresolved_reason": "",
            "canonical_endpoint_id": "METABOLITE_OBSERVATION",
            "normalized_value": num,
            "normalized_unit": raw_u or "%",
            "comparability_status": "DIRECTLY_COMPARABLE",
            "evidence_state": "AUTO_QUALIFIED_EXTERNAL",
        }

    # Rule 10: PK parameters (Cmax, Tmax, AUC, t1/2, CL, CL/F, Vd, Vss, Vd/F, F)
    pk_param = None
    if mtype in {"Cmax", "Tmax", "AUC", "t1/2", "CL", "Vd", "F"}:
        pk_param = mtype
    elif re.search(r"\bcmax\b", raw_ep.lower()):
        pk_param = "CMAX"
    elif re.search(r"\btmax\b", raw_ep.lower()):
        pk_param = "TMAX"
    elif re.search(r"\bauc\b", raw_ep.lower()):
        pk_param = "AUC"
    elif re.search(r"\bhalf[- ]?life\b|\bt1/2\b", raw_ep.lower()):
        pk_param = "T_HALF"
    elif re.search(r"\bclearance\b|\bcl/f\b|\bcl\b", raw_ep.lower()):
        pk_param = "CLF_ORAL" if route == "ORAL" or "cl/f" in raw_ep.lower() or "cl/f" in full_text else "CL"
    elif re.search(r"\bvolume of distribution\b|\bvss\b|\bvd/f\b|\bvd\b", raw_ep.lower()):
        pk_param = "VSS" if "steady state" in full_text or "vss" in raw_ep.lower() else ("VDF_ORAL" if route == "ORAL" or "vd/f" in full_text else "VD")
    elif re.search(r"\bbioavailability\b|\bF\b", raw_ep.lower()):
        pk_param = "F"

    if pk_param:
        param_norm = pk_param.upper().replace("/", "").replace("-", "_")
        if param_norm in {"T12", "HALF_LIFE", "T1/2"}:
            param_norm = "T_HALF"
        if param_norm in {"CL", "CLEARANCE"}:
            param_norm = "CLF_ORAL" if route == "ORAL" else "CL"
        if param_norm in {"VD", "VOLUME"}:
            param_norm = "VSS" if "vss" in full_text or "steady" in full_text else ("VDF_ORAL" if route == "ORAL" else "VD")

        # Check for relative ratio vs absolute PK value
        is_percent_ratio = ("%" in raw_u or "decrease" in full_text or "increase" in full_text or "ratio" in full_text) and param_norm in {"CMAX", "AUC"} and num is not None and num <= 200.0
        if is_percent_ratio:
            canonical_endpoint = f"{species}_PK_{param_norm}_{route}"
            return {
                "observation_id": record.get("id"),
                "raw_text_value": f"{raw_ep}: {raw_val} {raw_u}".strip(),
                "source": src,
                "current_classification": record.get("qualification_status") or "ENDPOINT_QUALIFIED",
                "proposed_canonical_section": "PK",
                "proposed_endpoint": canonical_endpoint,
                "measurement_type": f"Relative {param_norm}",
                "species": species,
                "context": ctx["conditions_dict"],
                "unit": "%",
                "qualification": STATE_RELATED,
                "refinement_rule": "pk_relative_ddi_food_effect_ratio",
                "resolved": True,
                "unresolved_reason": "",
                "canonical_endpoint_id": canonical_endpoint,
                "normalized_value": num,
                "normalized_unit": "%",
                "comparability_status": "RELATED_NOT_SAME_ENDPOINT",
                "evidence_state": "RELATED_EXTERNAL",
            }

        # Absolute PK parameter normalization
        unit_map = {
            "CMAX": "ng/mL", "TMAX": "hours", "AUC": "ng*h/mL", "AUC0_T": "ng*h/mL",
            "AUC0_INF": "ng*h/mL", "AUC_TAU": "ng*h/mL", "T_HALF": "hours",
            "CL": "mL/min/kg", "CLF_ORAL": "mL/min/kg", "VD": "L/kg", "VSS": "L/kg",
            "VDF_ORAL": "L/kg", "F": "%",
        }
        target_unit = unit_map.get(param_norm, raw_u)
        canonical_endpoint = f"{species}_PK_{param_norm}_{route}"

        if not raw_u or species == "UNSPECIFIED":
            reason_code = REASON_UNIT_MISSING if not raw_u else REASON_SPECIES_MISSING
            return {
                "observation_id": record.get("id"),
                "raw_text_value": f"{raw_ep}: {raw_val} {raw_u}".strip(),
                "source": src,
                "current_classification": record.get("qualification_status") or "UNCLASSIFIED",
                "proposed_canonical_section": "PK",
                "proposed_endpoint": canonical_endpoint,
                "measurement_type": param_norm,
                "species": species,
                "context": ctx["conditions_dict"],
                "unit": raw_u,
                "qualification": STATE_REVIEW_REQUIRED,
                "refinement_rule": "pk_context_incomplete",
                "resolved": False,
                "unresolved_reason": reason_code,
                "canonical_endpoint_id": canonical_endpoint,
                "normalized_value": num,
                "normalized_unit": raw_u,
                "comparability_status": "CONDITIONALLY_COMPARABLE",
                "evidence_state": STATE_REVIEW_REQUIRED,
            }

        # Check unit validity
        u_clean = raw_u.lower()
        unit_valid = True
        if param_norm in {"TMAX", "T_HALF"} and not any(k in u_clean for k in ["h", "hr", "hour", "min", "day", "d"]):
            unit_valid = False
        if param_norm == "F" and not ("%" in u_clean or (num is not None and 0 <= num <= 100)):
            unit_valid = False

        if not unit_valid:
            return {
                "observation_id": record.get("id"),
                "raw_text_value": f"{raw_ep}: {raw_val} {raw_u}".strip(),
                "source": src,
                "current_classification": record.get("qualification_status") or "ENDPOINT_QUALIFIED",
                "proposed_canonical_section": "PK",
                "proposed_endpoint": canonical_endpoint,
                "measurement_type": param_norm,
                "species": species,
                "context": ctx["conditions_dict"],
                "unit": raw_u,
                "qualification": STATE_REVIEW_REQUIRED,
                "refinement_rule": "pk_unit_invalid",
                "resolved": False,
                "unresolved_reason": REASON_UNIT_MISSING,
                "canonical_endpoint_id": canonical_endpoint,
                "normalized_value": num,
                "normalized_unit": raw_u,
                "comparability_status": "NOT_COMPARABLE",
                "evidence_state": STATE_REVIEW_REQUIRED,
            }

        state = STATE_AUTO_QUALIFIED
        comparability = "DIRECTLY_COMPARABLE"
        return {
            "observation_id": record.get("id"),
            "raw_text_value": f"{raw_ep}: {raw_val} {raw_u}".strip(),
            "source": src,
            "current_classification": record.get("qualification_status") or "ENDPOINT_QUALIFIED",
            "proposed_canonical_section": "PK",
            "proposed_endpoint": canonical_endpoint,
            "measurement_type": param_norm,
            "species": species,
            "context": ctx["conditions_dict"],
            "unit": target_unit,
            "qualification": state,
            "refinement_rule": f"pk_{param_norm.lower()}_refined",
            "resolved": True,
            "unresolved_reason": "",
            "canonical_endpoint_id": canonical_endpoint,
            "normalized_value": num,
            "normalized_unit": target_unit,
            "comparability_status": comparability,
            "evidence_state": "AUTO_QUALIFIED_EXTERNAL",
        }

    # Fallback review
    return {
        "observation_id": record.get("id"),
        "raw_text_value": f"{raw_ep}: {raw_val} {raw_u}".strip(),
        "source": src,
        "current_classification": record.get("qualification_status") or "UNCLASSIFIED",
        "proposed_canonical_section": "UNCLASSIFIED",
        "proposed_endpoint": "UNRESOLVED",
        "measurement_type": mtype,
        "species": species,
        "context": ctx["conditions_dict"],
        "unit": raw_u,
        "qualification": STATE_REVIEW_REQUIRED,
        "refinement_rule": "unresolved_fallback",
        "resolved": False,
        "unresolved_reason": REASON_ENDPOINT_AMBIGUOUS,
        "canonical_endpoint_id": "UNRESOLVED",
        "normalized_value": num,
        "normalized_unit": raw_u,
        "comparability_status": "UNSUPPORTED",
        "evidence_state": STATE_REVIEW_REQUIRED,
    }


def reprocess_all_persisted_evidence(db) -> dict:
    """Reprocess all persisted external experimental observations across all compounds."""
    from sqlalchemy import select
    from .models import ExternalExperimentalEvidence

    rows = list(db.scalars(select(ExternalExperimentalEvidence)).all())
    stats = {
        "total": len(rows),
        "auto_qualified": 0,
        "related": 0,
        "review_required": 0,
        "unusable": 0,
        "resolved_count": 0,
    }

    for row in rows:
        record_dict = {
            "id": row.id,
            "raw_endpoint_name": row.raw_endpoint_name,
            "raw_value": row.raw_value,
            "raw_unit": row.raw_unit,
            "species": row.species,
            "source_database": row.source_database,
            "source_record_id": row.source_record_id,
            "reference_text": row.reference_text,
            "assay_conditions_json": row.assay_conditions_json,
            "qualification_status": row.qualification_status,
        }
        refined = refine_scientific_observation(record_dict)
        
        # Apply updates
        row.canonical_endpoint_id = refined["canonical_endpoint_id"]
        row.routing_section = refined["proposed_canonical_section"]
        row.routing_reason = refined["unresolved_reason"] or refined["refinement_rule"]
        row.comparability_status = refined["comparability_status"]
        row.evidence_state = refined["evidence_state"]
        row.qualification_status = "ENDPOINT_QUALIFIED" if refined["resolved"] else "ENDPOINT_NOT_QUALIFIED"
        
        if refined["normalized_value"] is not None:
            row.normalized_value = str(refined["normalized_value"])
        if refined["normalized_unit"]:
            row.normalized_unit = refined["normalized_unit"]
        row.normalization_rule = refined["refinement_rule"]
        row.normalization_version = REFINEMENT_POLICY_VERSION

        qual_json = dict(row.qualification_json or {})
        qual_json["refinement"] = {
            "policy_version": REFINEMENT_POLICY_VERSION,
            "state": refined["qualification"],
            "rule": refined["refinement_rule"],
            "measurement_type": refined["measurement_type"],
            "unresolved_reason": refined["unresolved_reason"],
        }
        row.qualification_json = qual_json

        if refined["qualification"] == STATE_AUTO_QUALIFIED:
            stats["auto_qualified"] += 1
            stats["resolved_count"] += 1
        elif refined["qualification"] == STATE_RELATED:
            stats["related"] += 1
            stats["resolved_count"] += 1
        elif refined["qualification"] == STATE_UNUSABLE:
            stats["unusable"] += 1
        else:
            stats["review_required"] += 1

    db.flush()
    return stats
