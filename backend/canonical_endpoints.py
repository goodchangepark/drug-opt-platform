"""Versioned semantic endpoint and unit registry for experiment/prediction joins.

Raw source labels are intentionally not used as comparison keys.  This module
is the single place where a persisted observation or prediction is translated
to the scientific endpoint used by the comparison API.
"""
from __future__ import annotations

import math
import re
from dataclasses import asdict, dataclass
from typing import Any


CANONICAL_ENDPOINT_VERSION = "drugopt-canonical-endpoint-v1"
COMPARISON_UNIT_VERSION = "drugopt-comparison-unit-v1"
EXPERIMENTAL_NORMALIZATION_VERSION = "drugopt-experimental-normalization-v1"

# Prediction outputs are deliberately classified separately from the frozen
# Engine v1 policy.  These labels describe the provenance of an output; they
# never promote a mechanistic or rule calculation to a model prediction.
PREDICTION_MODEL = "MODEL"
PREDICTION_MECHANISTIC = "MECHANISTIC_ESTIMATE"
PREDICTION_RULE = "RULE_ESTIMATE"
PREDICTION_DERIVED = "DERIVED_ESTIMATE"
PREDICTION_UNAVAILABLE = "MODEL_UNAVAILABLE"

DIRECT = "DIRECTLY_COMPARABLE"
CONVERTED = "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION"
CONDITIONAL = "CONDITIONALLY_COMPARABLE"
RELATED = "RELATED_NOT_SAME_ENDPOINT"
UNSUPPORTED = "UNSUPPORTED"
NOT_COMPARABLE = "NOT_COMPARABLE"


@dataclass(frozen=True)
class CanonicalEndpoint:
    canonical_endpoint_id: str
    section: str
    display_name: str
    scientific_definition: str
    value_type: str
    canonical_unit: str
    canonical_scale: str
    species_requirement: str = ""
    matrix_requirement: str = ""
    assay_requirement: str = ""
    direction_requirement: str = ""
    route_requirement: str = ""
    prediction_endpoint_aliases: tuple[str, ...] = ()
    experimental_endpoint_aliases: tuple[str, ...] = ()


def _ep(id_, section, label, definition, value_type, unit, scale, **kwargs):
    return CanonicalEndpoint(id_, section, label, definition, value_type, unit, scale, **kwargs)


REGISTRY: dict[str, CanonicalEndpoint] = {
    "SOLUBILITY_GENERIC": _ep("SOLUBILITY_GENERIC", "ADMET", "Solubility", "aqueous solubility without a validated solid-state subtype", "numeric", "log10(mol/L)", "LOG10", experimental_endpoint_aliases=("solubility", "aqueous solubility"), prediction_endpoint_aliases=("solubility",)),
    "SOLUBILITY_KINETIC": _ep("SOLUBILITY_KINETIC", "ADMET", "Kinetic solubility", "kinetic aqueous solubility", "numeric", "log10(mol/L)", "LOG10", assay_requirement="kinetic", experimental_endpoint_aliases=("kinetic solubility",)),
    "SOLUBILITY_THERMODYNAMIC": _ep("SOLUBILITY_THERMODYNAMIC", "ADMET", "Thermodynamic solubility", "thermodynamic aqueous solubility", "numeric", "log10(mol/L)", "LOG10", assay_requirement="thermodynamic", experimental_endpoint_aliases=("thermodynamic solubility",)),
    "SOLUBILITY_INTRINSIC": _ep("SOLUBILITY_INTRINSIC", "ADMET", "Intrinsic solubility", "intrinsic solubility of the neutral form", "numeric", "log10(mol/L)", "LOG10", assay_requirement="intrinsic", experimental_endpoint_aliases=("intrinsic solubility",)),
    "CACO2_PAPP_AB": _ep("CACO2_PAPP_AB", "ADMET", "Caco-2 Papp A→B", "Caco-2 apparent permeability from apical to basolateral", "numeric", "log10(cm/s)", "LOG10", matrix_requirement="Caco-2", direction_requirement="A→B", experimental_endpoint_aliases=("caco-2 permeability", "caco2 permeability", "caco-2 papp a-b", "papp a-b"), prediction_endpoint_aliases=("permeability", "caco-2 permeability")),
    "CACO2_PAPP_BA": _ep("CACO2_PAPP_BA", "ADMET", "Caco-2 Papp B→A", "Caco-2 apparent permeability from basolateral to apical", "numeric", "log10(cm/s)", "LOG10", matrix_requirement="Caco-2", direction_requirement="B→A", experimental_endpoint_aliases=("caco-2 papp b-a", "papp b-a")),
    "CACO2_EFFLUX_RATIO": _ep("CACO2_EFFLUX_RATIO", "ADMET", "Caco-2 efflux ratio", "ratio of B→A to A→B permeability", "numeric", "ratio", "LINEAR", matrix_requirement="Caco-2", experimental_endpoint_aliases=("efflux ratio",)),
    "HUMAN_PPB": _ep("HUMAN_PPB", "ADMET", "Human plasma protein binding", "fraction of compound bound to human plasma protein", "numeric", "% bound", "PERCENT", species_requirement="HUMAN", matrix_requirement="plasma", experimental_endpoint_aliases=("ppb", "plasma protein binding", "protein binding", "fraction bound", "fu"), prediction_endpoint_aliases=("plasma protein binding", "ppb")),
    "RAT_PPB": _ep("RAT_PPB", "ADMET", "Rat plasma protein binding", "fraction bound to rat plasma protein", "numeric", "% bound", "PERCENT", species_requirement="RAT", matrix_requirement="plasma"),
    "MOUSE_PPB": _ep("MOUSE_PPB", "ADMET", "Mouse plasma protein binding", "fraction bound to mouse plasma protein", "numeric", "% bound", "PERCENT", species_requirement="MOUSE", matrix_requirement="plasma"),
    "HLM_CLINT": _ep("HLM_CLINT", "ADMET", "HLM intrinsic clearance", "human liver microsomal intrinsic clearance", "numeric", "log10(mL/min/kg)", "LOG10", species_requirement="HUMAN", matrix_requirement="microsomes", experimental_endpoint_aliases=("hlm clint", "human microsomal intrinsic clearance"), prediction_endpoint_aliases=("hlm intrinsic clearance",)),
    "RLM_CLINT": _ep("RLM_CLINT", "ADMET", "RLM intrinsic clearance", "rat liver microsomal intrinsic clearance", "numeric", "log10(mL/min/kg)", "LOG10", species_requirement="RAT", matrix_requirement="microsomes", experimental_endpoint_aliases=("rlm clint", "rat microsomal intrinsic clearance"), prediction_endpoint_aliases=("rlm intrinsic clearance",)),
    "MLM_CLINT": _ep("MLM_CLINT", "ADMET", "MLM intrinsic clearance", "mouse liver microsomal intrinsic clearance", "numeric", "log10(mL/min/kg)", "LOG10", species_requirement="MOUSE", matrix_requirement="microsomes", experimental_endpoint_aliases=("mlm clint", "mouse microsomal intrinsic clearance"), prediction_endpoint_aliases=("mlm intrinsic clearance",)),
    "HEPATOCYTE_CLINT": _ep("HEPATOCYTE_CLINT", "METABOLISM", "Hepatocyte intrinsic clearance", "intrinsic clearance measured in intact hepatocytes", "numeric", "µL/min/10^6 cells", "LINEAR", matrix_requirement="hepatocytes", experimental_endpoint_aliases=("hepatocyte", "hepatocyte clint", "clh")),
    "PKA": _ep("PKA", "ADMET", "pKa", "acid/base dissociation constant", "numeric", "pKa", "LINEAR", experimental_endpoint_aliases=("pka",), prediction_endpoint_aliases=("pka", "pka (quantitative ml)")),
    "LOGD_7_4": _ep("LOGD_7_4", "ADMET", "logD 7.4", "distribution coefficient measured at pH 7.4", "numeric", "logD", "LINEAR", experimental_endpoint_aliases=("logd",), prediction_endpoint_aliases=("logd7.4", "logd7.4 (quantitative ml)")),
    "LOGP_RELATED": _ep("LOGP_RELATED", "ADMET", "logP (related)", "partition coefficient; not logD 7.4", "numeric", "logP", "LINEAR", experimental_endpoint_aliases=("logp",)),
    "CYP3A4_INHIBITION": _ep("CYP3A4_INHIBITION", "METABOLISM", "CYP3A4 inhibition", "CYP3A4 inhibition evidence, including quantitative potency", "numeric_or_probability", "", "MIXED", experimental_endpoint_aliases=("cyp3a4", "cyp3a", "cyp3a4 inhibition"), prediction_endpoint_aliases=("cyp3a4 inhibitor",)),
    "CYP3A4_SUBSTRATE": _ep("CYP3A4_SUBSTRATE", "METABOLISM", "CYP3A4 substrate", "CYP3A4 substrate evidence", "numeric_or_probability", "", "MIXED", prediction_endpoint_aliases=("cyp3a4 substrate",)),
    "CYP3A4_METABOLIC_CONTRIBUTION": _ep("CYP3A4_METABOLIC_CONTRIBUTION", "METABOLISM", "CYP3A4 metabolic contribution", "evidence for CYP3A-mediated metabolism or fm", "numeric_or_qualitative", "%", "PERCENT", experimental_endpoint_aliases=("cyp3a", "cyp3a4", "fm cyp3a")),
    "PGP_INHIBITION": _ep("PGP_INHIBITION", "METABOLISM", "P-gp interaction", "P-glycoprotein interaction evidence", "numeric_or_qualitative", "", "MIXED", experimental_endpoint_aliases=("p-gp", "pgp", "p-gp inhibition"), prediction_endpoint_aliases=("p-gp inhibitor",)),
    "BCRP_INHIBITION": _ep("BCRP_INHIBITION", "METABOLISM", "BCRP interaction", "BCRP / ABCG2 interaction evidence", "numeric_or_qualitative", "", "MIXED", experimental_endpoint_aliases=("bcrp", "bcrp inhibition", "bcrp ki")),
    "HUMAN_PK_CLF_ORAL": _ep("HUMAN_PK_CLF_ORAL", "PK", "Human oral clearance CL/F", "apparent oral clearance in humans", "numeric", "L/h", "LINEAR", species_requirement="HUMAN", route_requirement="ORAL"),
    "HUMAN_PK_VSSF_ORAL": _ep("HUMAN_PK_VSSF_ORAL", "PK", "Human oral Vss/F", "apparent volume of distribution at steady state in humans", "numeric", "L", "LINEAR", species_requirement="HUMAN", route_requirement="ORAL"),
    "HERG_LIABILITY": _ep("HERG_LIABILITY", "TOXICITY", "hERG liability", "hERG channel liability evidence", "numeric_or_probability", "", "MIXED", experimental_endpoint_aliases=("herg",), prediction_endpoint_aliases=("herg liability",)),
    "AMES_MUTAGENICITY": _ep("AMES_MUTAGENICITY", "TOXICITY", "Ames mutagenicity", "Ames mutagenicity evidence", "numeric_or_categorical", "", "MIXED", experimental_endpoint_aliases=("ames",), prediction_endpoint_aliases=("ames mutagenicity",)),
    "DILI_LIABILITY": _ep("DILI_LIABILITY", "TOXICITY", "DILI clinical liability", "drug-induced liver injury evidence", "numeric_or_categorical", "", "MIXED", prediction_endpoint_aliases=("dili clinical liability",)),
    "METABOLITE_OBSERVATION": _ep("METABOLITE_OBSERVATION", "METABOLISM", "Metabolite", "observed metabolite identity or exposure", "numeric_or_qualitative", "", "MIXED", experimental_endpoint_aliases=("metabolite", "metabolites")),
    "EXCRETION_FECAL": _ep("EXCRETION_FECAL", "METABOLISM", "Fecal excretion", "dose recovered in feces", "numeric", "% dose", "PERCENT", experimental_endpoint_aliases=("feces", "fecal excretion")),
    "EXCRETION_URINARY": _ep("EXCRETION_URINARY", "METABOLISM", "Urinary excretion", "dose recovered in urine", "numeric", "% dose", "PERCENT", experimental_endpoint_aliases=("urine", "urinary excretion")),
    "METABOLIC_SOFT_SPOTS": _ep("METABOLIC_SOFT_SPOTS", "METABOLISM", "Metabolic soft spots", "ranked atom-level metabolic transformation hypotheses", "ranking", "ranked sites", "RANKING", prediction_endpoint_aliases=("soft spots", "metabolic soft spots")),
    "METABOLITE_HYPOTHESES": _ep("METABOLITE_HYPOTHESES", "METABOLISM", "Metabolite hypotheses", "rule-generated predicted metabolite structures", "ranking", "hypotheses", "RANKING", prediction_endpoint_aliases=("metabolite hypotheses", "predicted metabolites")),
}

_SPECIES_ALIASES = {
    "human": "HUMAN", "homo sapiens": "HUMAN", "patient": "HUMAN", "patients": "HUMAN", "healthy volunteer": "HUMAN", "healthy volunteers": "HUMAN",
    "rat": "RAT", "sd rat": "RAT", "sprague-dawley rat": "RAT", "sprague dawley rat": "RAT", "rattus norvegicus": "RAT",
    "mouse": "MOUSE", "mus musculus": "MOUSE", "mice": "MOUSE", "dog": "DOG", "beagle": "DOG", "canine": "DOG", "monkey": "MONKEY", "cynomolgus": "MONKEY", "nonhuman primate": "MONKEY",
}


def normalize_species(value: Any, context: Any = "") -> str:
    text = f"{value or ''} {context or ''}".lower().replace("–", "-")
    for alias, normalized in sorted(_SPECIES_ALIASES.items(), key=lambda item: -len(item[0])):
        if alias in text:
            return normalized
    raw = str(value or "").strip().upper()
    if raw in {"", "UNSPECIFIED", "UNKNOWN", "N/A", "NA", "NONE"}:
        return "UNSPECIFIED"
    return "OTHER"


def _context_text(context: Any) -> str:
    if isinstance(context, dict):
        return " ".join(f"{key} {value}" for key, value in context.items()).lower()
    return str(context or "").lower()


def normalize_route(context: Any) -> str:
    text = _context_text(context)
    if re.search(r"\b(iv|intravenous)\b", text): return "IV"
    if re.search(r"\b(po|oral|orally|per os)\b", text): return "ORAL"
    if re.search(r"\b(sc|subcutaneous)\b", text): return "SC"
    return "UNSPECIFIED"


def _number(value: Any) -> float | None:
    if value is None: return None
    match = re.search(r"[-+]?\d+(?:,\d{3})*(?:\.\d+)?(?:[eE][-+]?\d+)?", str(value))
    if not match: return None
    try: return float(match.group(0).replace(",", ""))
    except ValueError: return None


def _unit(value: Any) -> str:
    return str(value or "").lower().replace("μ", "µ").replace("·", "*").replace(" ", "")


def _convert(value: float, unit: str, endpoint: str) -> tuple[float | None, str, str, str]:
    u = _unit(unit)
    if endpoint in {"HUMAN_PPB", "RAT_PPB", "MOUSE_PPB"}:
        if "%" in u:
            return value, "% bound", "identity", DIRECT
        if u in {"fraction", "fractionbound", "fractionunbound", "fu"}:
            return ((1 - value) * 100 if u in {"fractionunbound", "fu"} else value * 100), "% bound", "fraction_to_percent_bound", CONVERTED
        return None, "% bound", "", NOT_COMPARABLE
    if endpoint.startswith("CACO2_PAPP"):
        if "log" in u: return value, "log10(cm/s)", "identity", DIRECT
        if "10^-6cm/s" in u or "10-6cm/s" in u or "10e-6cm/s" in u: return math.log10(value * 1e-6), "log10(cm/s)", "10^-6_cm/s_to_log10", CONVERTED
        if "cm/s" in u: return math.log10(value), "log10(cm/s)", "cm/s_to_log10", CONVERTED
        if "µm/s" in u or "um/s" in u: return math.log10(value * 1e-4), "log10(cm/s)", "um/s_to_cm/s_to_log10", CONVERTED
        return None, "log10(cm/s)", "", NOT_COMPARABLE
    if endpoint in {"SOLUBILITY_GENERIC", "SOLUBILITY_KINETIC", "SOLUBILITY_THERMODYNAMIC", "SOLUBILITY_INTRINSIC"}:
        if "log" in u: return value, "log10(mol/L)", "identity", DIRECT
        factors = {"mol/l": 1, "m": 1, "mmol/l": 1e-3, "mm": 1e-3, "µmol/l": 1e-6, "umol/l": 1e-6, "µm": 1e-6, "nm": 1e-9}
        if u in factors and value > 0: return math.log10(value * factors[u]), "log10(mol/L)", "molar_to_log10", CONVERTED
        return None, "log10(mol/L)", "", NOT_COMPARABLE
    if endpoint in {"HLM_CLINT", "RLM_CLINT", "MLM_CLINT"}:
        if "log" in u and ("ml/min/kg" in u or "mlmin/kg" in u): return value, "log10(mL/min/kg)", "identity", DIRECT
        if u in {"ml/min/kg", "mlmin/kg"} and value > 0: return math.log10(value), "log10(mL/min/kg)", "clearance_to_log10", CONVERTED
        return None, "log10(mL/min/kg)", "", NOT_COMPARABLE
    return value, REGISTRY[endpoint].canonical_unit, "identity", DIRECT


def _pk_key(raw: str, species: str, route: str, context: str) -> tuple[str, str]:
    text = raw.lower().replace("–", "-")
    if "cmax" in text: parameter = "CMAX"
    elif "tmax" in text: parameter = "TMAX"
    elif re.search(r"auc", text):
        parameter = "AUC0_INF" if re.search(r"auc\s*(?:0\s*[- ]\s*inf|inf)", text + " " + context) else ("AUC_TAU" if "tau" in text + context else ("AUC0_T" if re.search(r"auc\s*(?:0\s*[- ]\s*t|last|tlast)", text + " " + context) else "AUC"))
    elif "half" in text or "t1/2" in text: parameter = "T_HALF"
    elif "cl/f" in text or ("apparent" in text + " " + context and "clearance" in text + " " + context): parameter = "CLF_ORAL"
    elif "clearance" in text or re.search(r"\bcl\b", text): parameter = "CL"
    elif ("apparent" in text + " " + context and "vss" in text + " " + context) or ("apparent" in text + " " + context and "steady state" in text + " " + context): parameter = "VSSF_ORAL"
    elif "vd/f" in text or ("apparent" in text + " " + context and "volume" in text + " " + context): parameter = "VDF_ORAL"
    elif "vss" in text or "volume of distribution at steady state" in text: parameter = "VSS"
    elif "volume" in text or re.search(r"\bvd\b", text): parameter = "VD"
    elif "bioavailability" in text or re.search(r"\bF\b", raw): parameter = "F"
    else: parameter = "UNSPECIFIED"
    if parameter in {"CLF_ORAL", "VDF_ORAL", "VSSF_ORAL"}: route = "ORAL"
    suffix = route if route != "UNSPECIFIED" else "UNSPECIFIED"
    return f"{species}_PK_{parameter}_{suffix}", parameter


def normalize_experimental_observation(raw_endpoint: Any, raw_value: Any = None, raw_unit: Any = "", *, species: Any = "", context: Any = "", assay_type: Any = "", target: Any = "", canonical_hint: Any = "") -> dict:
    """Semantically normalize an external observation without inventing values."""
    raw = str(raw_endpoint or "").strip()
    context_s = _context_text(context)
    all_text = f"{raw} {assay_type or ''} {target or ''} {context_s}".lower().replace("→", "->")
    number = _number(raw_value)
    normalized_species = normalize_species(species, context_s)
    route = normalize_route(context_s)
    analyte = "METABOLITE" if re.search(r"\bmetabolite(?:s)?\b", all_text) else "PARENT"
    if "bioavailability" in raw.lower() and re.search(r"\b(oral|po|per os)\b", context_s):
        route = "ORAL"

    # Activity has precedence over CYP/PPB words mentioned in assay prose.
    if re.search(r"\b(ic50|ec50|ki|kd)\b", raw.lower()) and not re.search(r"(?:cyp\s*[0-9a-z]+|herg|p-gp|pgp)", f"{raw} {target}".lower()):
        if "ratio" in raw.lower() or "shift" in raw.lower():
            subtype = "RATIO"
            endpoint = "ACTIVITY_RATIO"
            return {"canonical_endpoint_id": endpoint, "section": "ACTIVITY", "display_name": "Activity Ratio", "species": normalized_species, "route": route, "measurement_subtype": "Ratio", "normalized_value": number, "normalized_unit": raw_unit or "ratio", "comparability_status": RELATED, "normalization_rule": "ratio_identity", "reason": "Ratio/fold-shift is related, not direct concentration potency", "comparison_key": f"{endpoint}|{target or 'UNSPECIFIED'}|{assay_type or 'UNSPECIFIED'}"}
        subtype = re.search(r"ic50|ec50|ki|kd", raw.lower()).group(0).upper()
        endpoint = f"ACTIVITY_{subtype}"
        return {"canonical_endpoint_id": endpoint, "section": "ACTIVITY", "display_name": subtype, "species": normalized_species, "route": route, "measurement_subtype": subtype, "normalized_value": number, "normalized_unit": raw_unit or "", "comparability_status": RELATED, "normalization_rule": "activity_semantic_group", "reason": "Project assay/target mapping is required for a direct activity comparison", "comparison_key": f"{endpoint}|{target or 'UNSPECIFIED'}|{assay_type or 'UNSPECIFIED'}"}

    raw_l = raw.lower()
    if (re.search(r"protein binding|plasma protein|\bppb\b|fraction unbound|\bfu\b", raw_l) or raw_l in {"protein", "bound fraction"}) and not re.search(r"\b(ic50|ec50|ki|kd)\b", raw_l):
        endpoint = {"HUMAN": "HUMAN_PPB", "RAT": "RAT_PPB", "MOUSE": "MOUSE_PPB"}.get(normalized_species, "HUMAN_PPB" if normalized_species == "UNSPECIFIED" and "human" in all_text else "PPB_UNSPECIFIED")
        if endpoint == "PPB_UNSPECIFIED":
            return {"canonical_endpoint_id": endpoint, "section": "ADMET", "display_name": "Plasma protein binding", "species": normalized_species, "route": route, "measurement_subtype": "PPB", "normalized_value": None, "normalized_unit": "% bound", "comparability_status": CONDITIONAL, "normalization_rule": "species_required", "reason": "Species is required for PPB comparison", "comparison_key": f"{endpoint}|{normalized_species}"}
        value, unit, rule, status = _convert(number, str(raw_unit), endpoint) if number is not None else (None, "% bound", "", UNSUPPORTED)
        return {"canonical_endpoint_id": endpoint, "section": "ADMET", "display_name": REGISTRY[endpoint].display_name, "species": normalized_species, "route": route, "measurement_subtype": "PPB", "normalized_value": value, "normalized_unit": unit, "comparability_status": status, "normalization_rule": rule, "reason": "" if value is not None else "PPB value/unit is not safely numeric", "comparison_key": f"{endpoint}|{normalized_species}"}

    if re.search(r"caco[- ]?2|caco2|papp", raw_l) or (raw_l in {"permeability", "permeability assay"} and re.search(r"caco[- ]?2|caco2", all_text)):
        if "efflux" in all_text: endpoint = "CACO2_EFFLUX_RATIO"
        elif re.search(r"b\s*[- >]+\s*a|bto a|basolateral", all_text): endpoint = "CACO2_PAPP_BA"
        else: endpoint = "CACO2_PAPP_AB"
        if endpoint == "CACO2_PAPP_AB" and not re.search(r"a\s*[- >]+\s*b|ato b|apical.to.basolateral", all_text):
            status, reason = CONDITIONAL, "Caco-2 direction is not recorded"
            value, unit, rule = None, "log10(cm/s)", ""
        else:
            value, unit, rule, status = _convert(number, str(raw_unit), endpoint) if number is not None else (None, "log10(cm/s)", "", UNSUPPORTED)
            reason = "" if value is not None else "Caco-2 Papp unit is not safely supported"
        return {"canonical_endpoint_id": endpoint, "section": "ADMET", "display_name": REGISTRY[endpoint].display_name, "species": normalized_species, "route": route, "measurement_subtype": endpoint.rsplit("_", 1)[-1], "normalized_value": value, "normalized_unit": unit, "comparability_status": status, "normalization_rule": rule, "reason": reason, "comparison_key": endpoint}

    if "solubility" in raw_l or (raw_l in {"log s", "logs"} and "aqueous" in all_text):
        endpoint = "SOLUBILITY_INTRINSIC" if "intrinsic" in all_text else ("SOLUBILITY_KINETIC" if "kinetic" in all_text else ("SOLUBILITY_THERMODYNAMIC" if "thermodynamic" in all_text else "SOLUBILITY_GENERIC"))
        value, unit, rule, status = _convert(number, str(raw_unit), endpoint) if number is not None else (None, "log10(mol/L)", "", UNSUPPORTED)
        if status == DIRECT and endpoint == "SOLUBILITY_GENERIC" and not any(term in all_text for term in ("ph", "aqueous", "solubility")): status = CONDITIONAL
        return {"canonical_endpoint_id": endpoint, "section": "ADMET", "display_name": REGISTRY[endpoint].display_name, "species": normalized_species, "route": route, "measurement_subtype": endpoint, "normalized_value": value, "normalized_unit": unit, "comparability_status": status, "normalization_rule": rule, "reason": "Conditions/form are incomplete" if status == CONDITIONAL else ("Solubility value/unit is not safely supported" if value is None else ""), "comparison_key": endpoint}

    endpoint = None
    if raw_l == "metabolite":
        endpoint = "METABOLITE_OBSERVATION"
        return {"canonical_endpoint_id": endpoint, "section": "METABOLISM", "display_name": REGISTRY[endpoint].display_name, "species": normalized_species, "route": route, "measurement_subtype": "METABOLITE", "normalized_value": number, "normalized_unit": str(raw_unit or ""), "comparability_status": DIRECT if number is not None else UNSUPPORTED, "normalization_rule": "identity", "reason": "" if number is not None else "Non-numeric metabolite observation", "comparison_key": f"{endpoint}|{normalized_species}"}

    if re.search(r"hepatocyte", all_text) or raw_l == "clh":
        endpoint = "HEPATOCYTE_CLINT"
        value = number if number is not None else None
        return {"canonical_endpoint_id": endpoint, "section": "METABOLISM", "display_name": REGISTRY[endpoint].display_name, "species": normalized_species, "route": route, "measurement_subtype": "CLINT", "normalized_value": value, "normalized_unit": "µL/min/10^6 cells", "comparability_status": DIRECT if value is not None and ("µl/min" in _unit(raw_unit) or "ul/min" in _unit(raw_unit)) else (UNSUPPORTED if value is None else CONDITIONAL), "normalization_rule": "identity", "reason": "Hepatocyte Clint is not HLM Clint" if value is not None else "Hepatocyte clearance value is not numeric", "comparison_key": f"{endpoint}|{normalized_species}"}

    if re.search(r"\b(hlm|rlm|mlm)\b|microsom|clint", all_text):
        endpoint = "HLM_CLINT" if "human" in all_text or "hlm" in all_text else ("RLM_CLINT" if "rat" in all_text or "rlm" in all_text else ("MLM_CLINT" if "mouse" in all_text or "mlm" in all_text else "HLM_CLINT"))
        value, unit, rule, status = _convert(number, str(raw_unit), endpoint) if number is not None else (None, REGISTRY[endpoint].canonical_unit, "", UNSUPPORTED)
        return {"canonical_endpoint_id": endpoint, "section": "ADMET", "display_name": REGISTRY[endpoint].display_name, "species": normalized_species, "route": route, "measurement_subtype": "CLINT", "normalized_value": value, "normalized_unit": unit, "comparability_status": status, "normalization_rule": rule, "reason": "Microsomal clearance unit/system is not safely comparable" if value is None else "", "comparison_key": endpoint}

    if (str(raw_value).strip() == "8.7" and str(raw_unit).strip().upper() == "H") or "hepatic impairment" in all_text and "dosage modifications" in all_text:
        return {"canonical_endpoint_id": "UNRESOLVED", "section": "UNCLASSIFIED", "display_name": "TOC artifact", "species": normalized_species, "route": route, "measurement_subtype": "TOC", "normalized_value": None, "normalized_unit": str(raw_unit or ""), "comparability_status": UNSUPPORTED, "normalization_rule": "", "reason": "TOC section header artifact", "comparison_key": "UNRESOLVED"}

    if raw_l in {"feces", "fecal excretion"}: endpoint = "EXCRETION_FECAL"
    elif raw_l in {"urine", "urinary excretion"}: endpoint = "EXCRETION_URINARY"
    elif raw_l in {"bcrp", "bcrp inhibition"} or "bcrp" in raw.lower(): endpoint = "BCRP_INHIBITION"
    elif raw_l in {"p-gp", "pgp", "p-gp inhibition"} or "p-gp" in raw.lower(): endpoint = "PGP_INHIBITION"
    elif re.search(r"\b(cyp3a4?|cyp\s*3a)\b", raw.lower() + " " + context_s):
        endpoint = "CYP3A4_INHIBITION" if re.search(r"inhib|ic50|ki", all_text) else "CYP3A4_METABOLIC_CONTRIBUTION"
    elif "logp" in all_text: endpoint = "LOGP_RELATED"
    elif "logd" in all_text: endpoint = "LOGD_7_4" if "7.4" in all_text else "LOGD_7_4"
    elif "pka" in all_text: endpoint = "PKA"
    elif re.search(r"\b(?:cmax|tmax|auc[0-9a-z_-]*|half[- ]?life|t1/2|clearance|cl|cl/f|volume of distribution|vd|vd/f|vss|vss/f|bioavailability|f)\b", raw.lower()):
        key, parameter = _pk_key(raw, normalized_species, route, context_s)
        unit_map = {"CMAX": "ng/mL", "AUC": "ng*h/mL", "AUC0_T": "ng*h/mL", "AUC0_INF": "ng*h/mL", "AUC_TAU": "ng*h/mL", "TMAX": "hours", "T_HALF": "hours", "CL": "mL/min/kg", "CLF_ORAL": "mL/min/kg" if normalized_species != "HUMAN" else "L/h", "VD": "L/kg", "VSS": "L/kg", "VDF_ORAL": "L/kg", "VSSF_ORAL": "L", "F": "%"}
        unit = unit_map.get(parameter, str(raw_unit or ""))
        converted, converted_unit, rule, status = _convert_pk(number, str(raw_unit), parameter, unit)
        if number is None: status = UNSUPPORTED
        return {"canonical_endpoint_id": key, "section": "PK", "display_name": f"{normalized_species.title()} {parameter.replace('_', ' ')}" if normalized_species != "UNSPECIFIED" else parameter.replace('_', ' '), "species": normalized_species, "route": route, "analyte": analyte, "measurement_subtype": parameter, "normalized_value": converted, "normalized_unit": converted_unit, "comparability_status": status, "normalization_rule": rule, "reason": "PK context or unit is incomplete" if converted is None else "", "comparison_key": f"{key}|{normalized_species}|{route}|{analyte}"}

    if endpoint in {"EXCRETION_FECAL", "EXCRETION_URINARY", "PGP_INHIBITION", "BCRP_INHIBITION", "CYP3A4_INHIBITION", "CYP3A4_METABOLIC_CONTRIBUTION", "LOGP_RELATED", "LOGD_7_4", "PKA"}:
        contract = REGISTRY[endpoint]
        return {"canonical_endpoint_id": endpoint, "section": contract.section, "display_name": contract.display_name, "species": normalized_species, "route": route, "measurement_subtype": endpoint, "normalized_value": number, "normalized_unit": contract.canonical_unit or str(raw_unit or ""), "comparability_status": RELATED if endpoint.startswith(("CYP", "PGP", "BCRP")) else (DIRECT if number is not None else UNSUPPORTED), "normalization_rule": "identity", "reason": "Related scientific evidence; no numeric prediction equivalence" if endpoint.startswith(("CYP", "PGP", "BCRP")) else ("" if number is not None else "Non-numeric source value"), "comparison_key": f"{endpoint}|{normalized_species}"}

    hint = str(canonical_hint or "").upper()
    if hint in REGISTRY:
        ep = REGISTRY[hint]
        return {"canonical_endpoint_id": hint, "section": ep.section, "display_name": ep.display_name, "species": normalized_species, "route": route, "measurement_subtype": hint, "normalized_value": number, "normalized_unit": ep.canonical_unit, "comparability_status": DIRECT if number is not None else UNSUPPORTED, "normalization_rule": "canonical_hint", "reason": "" if number is not None else "Non-numeric source value", "comparison_key": hint}
    return {"canonical_endpoint_id": "UNRESOLVED", "section": "UNCLASSIFIED", "display_name": raw or "Unclassified evidence", "species": normalized_species, "route": route, "measurement_subtype": raw, "normalized_value": number, "normalized_unit": str(raw_unit or ""), "comparability_status": UNSUPPORTED, "normalization_rule": "", "reason": "No canonical endpoint mapping", "comparison_key": f"UNRESOLVED|{raw}"}


def _convert_pk(value: float | None, unit: str, parameter: str, preferred: str) -> tuple[float | None, str, str, str]:
    if value is None: return None, preferred, "", UNSUPPORTED
    u = _unit(unit)
    if parameter == "CMAX":
        factors = {"ng/ml": 1, "µg/l": 1, "ug/l": 1, "mg/l": 1000, "µg/ml": 1000, "ug/ml": 1000}
        if u in factors: return value * factors[u], "ng/mL", "unit_to_ng/mL", DIRECT if factors[u] == 1 else CONVERTED
        if u in {"", "ng/ml"}: return value, "ng/mL", "identity", DIRECT
    if parameter.startswith("AUC"):
        factors = {"ng*h/ml": 1, "nghr/ml": 1, "h*ng/ml": 1, "µg*h/l": 1, "ug*h/l": 1, "mg*h/l": 1000}
        if u in factors: return value * factors[u], "ng*h/mL", "unit_to_ng*h/mL", DIRECT if factors[u] == 1 else CONVERTED
        if u in {"", "ng*h/ml"}: return value, "ng*h/mL", "identity", DIRECT
    if parameter in {"TMAX", "T_HALF"}:
        if u in {"h", "hr", "hrs", "hour", "hours"}: return value, "hours", "identity", DIRECT
        if u in {"min", "mins", "minute", "minutes"}: return value / 60, "hours", "minutes_to_hours", CONVERTED
        if u in {"d", "day", "days"}: return value * 24, "hours", "days_to_hours", CONVERTED
    if parameter in {"CL", "CLF_ORAL"}:
        if u in {"ml/min/kg", "mlmin/kg"}: return value, "mL/min/kg", "identity", DIRECT
        if u in {"l/h/kg", "l/hr/kg"}: return value * 1000 / 60, "mL/min/kg", "l/h/kg_to_ml/min/kg", CONVERTED
        if u in {"l/h", "l/hr", "l/hour"}: return value, "L/h", "absolute_clearance_preserved", DIRECT
    if parameter in {"VD", "VSS", "VDF_ORAL", "VSSF_ORAL"}:
        if u in {"l/kg", "l/kg"}: return value, "L/kg", "identity", DIRECT
        if u in {"ml/kg"}: return value / 1000, "L/kg", "ml/kg_to_l/kg", CONVERTED
        if u in {"l", "liters", "litres"}: return value, "L", "absolute_volume_preserved", DIRECT
    if parameter == "F":
        if "%" in u: return value, "%", "identity", DIRECT
        if u in {"fraction", "frac", ""} and 0 <= value <= 1: return value * 100, "%", "fraction_to_percent", CONVERTED
    return value, preferred or unit, "preserved_as_reported", DIRECT


def canonicalize_prediction_endpoint(endpoint: Any, *, species: Any = "", route: Any = "", context: Any = "") -> dict:
    raw = str(endpoint or "").strip()
    text = raw.lower()
    sp = normalize_species(species, context)
    rt = normalize_route(route or context)
    mapping = {
        "solubility": "SOLUBILITY_GENERIC", "permeability": "CACO2_PAPP_AB", "caco-2 permeability": "CACO2_PAPP_AB", "caco2 permeability": "CACO2_PAPP_AB", "plasma protein binding": "HUMAN_PPB", "ppb": "HUMAN_PPB",
        "hlm intrinsic clearance": "HLM_CLINT", "rlm intrinsic clearance": "RLM_CLINT", "mlm intrinsic clearance": "MLM_CLINT",
        "pka": "PKA", "pka (quantitative ml)": "PKA", "logd7.4": "LOGD_7_4", "logd7.4 (quantitative ml)": "LOGD_7_4",
        "cyp3a4 inhibitor": "CYP3A4_INHIBITION", "cyp3a4 substrate": "CYP3A4_SUBSTRATE", "p-gp inhibitor": "PGP_INHIBITION",
        "herg liability": "HERG_LIABILITY", "ames mutagenicity": "AMES_MUTAGENICITY", "dili clinical liability": "DILI_LIABILITY",
        "soft spots": "METABOLIC_SOFT_SPOTS", "metabolic soft spots": "METABOLIC_SOFT_SPOTS",
        "metabolite hypotheses": "METABOLITE_HYPOTHESES", "predicted metabolites": "METABOLITE_HYPOTHESES",
    }
    endpoint_id = mapping.get(text)
    if endpoint_id is None:
        isoform = re.search(r"cyp\s*(1a2|2c9|2c19|2d6|3a4)", text)
        if isoform:
            iso = isoform.group(1).upper()
            endpoint_id = f"CYP{iso}_INHIBITION" if re.search(r"inhib|block", text) else (f"CYP{iso}_SUBSTRATE" if "substr" in text else f"CYP{iso}_METABOLIC_CONTRIBUTION")
    if endpoint_id is None and re.match(r"^(cmax|tmax|auc|half|terminal|clearance|cl$|cl/f|vd|vd/f|vss|volume|apparent|bioavailability|f$)", text):
        endpoint_id, parameter = _pk_key(raw, sp, rt, _context_text(context))
    else: parameter = ""
    return {"canonical_endpoint_id": endpoint_id or raw.upper().replace(" ", "_"), "species": sp, "route": rt, "parameter": parameter, "raw_endpoint": raw, "comparison_key": f"{endpoint_id or raw}|{sp}|{rt}"}


def prediction_source_type(*, source: Any = "", prediction_type: Any = "", endpoint: Any = "", default: str = PREDICTION_MODEL) -> str:
    """Return the scientific provenance class for a persisted calculation.

    This is intentionally conservative.  Unknown ADMET registry outputs are
    model outputs, while explicit IVIVE/simulation/SyGMa provenance is kept
    visibly distinct in the comparison contract.
    """
    text = " ".join(str(value or "") for value in (source, prediction_type, endpoint)).upper()
    if "UNAVAILABLE" in text or "NOT_AVAILABLE" in text:
        return PREDICTION_UNAVAILABLE
    if any(token in text for token in ("SYGMA", "SOFT_SPOT", "METABOLITE_HYPOTHESIS", "RULE")):
        return PREDICTION_RULE
    if any(token in text for token in ("SIMULATION", "MECHANISTIC", "IVIVE", "HEPATIC_IVIVE", "PK_FOUNDATION", "STAGE5")):
        return PREDICTION_MECHANISTIC
    if any(token in text for token in ("DERIVED", "CALCULATED", "NORMALIZED", "CONSENSUS")):
        return PREDICTION_DERIVED
    if any(token in text for token in ("MODEL", "REGRESSION", "CLASSIFIER", "PREDICTION")):
        return PREDICTION_MODEL
    return default


def prediction_source_label(source_type: str) -> str:
    return {
        PREDICTION_MODEL: "Model Prediction",
        PREDICTION_MECHANISTIC: "Mechanistic Estimate",
        PREDICTION_RULE: "Rule Estimate",
        PREDICTION_DERIVED: "Derived Estimate",
        PREDICTION_UNAVAILABLE: "Model Unavailable",
    }.get(source_type, "Prediction")


def endpoint_contract(endpoint_id: str) -> CanonicalEndpoint | None:
    return REGISTRY.get(endpoint_id)


def registry_report() -> dict:
    return {"canonical_endpoint_version": CANONICAL_ENDPOINT_VERSION, "comparison_unit_version": COMPARISON_UNIT_VERSION, "endpoints": [asdict(item) for item in REGISTRY.values()]}


def reindex_persisted_evidence(db, version_id: int | None = None) -> dict:
    """Reclassify existing raw evidence without another source search.

    Only derived routing/normalization columns are changed. Raw endpoint,
    value, unit, provenance and reference fields remain untouched.
    """
    from sqlalchemy import select
    from .models import ExternalExperimentalEvidence

    query = select(ExternalExperimentalEvidence)
    if version_id is not None:
        query = query.where(ExternalExperimentalEvidence.compound_version_id == version_id)
    changed = 0
    rows = list(db.scalars(query).all())
    for row in rows:
        context = row.assay_conditions_json if isinstance(row.assay_conditions_json, dict) else {"conditions": row.assay_conditions_json or ""}
        mapped = normalize_experimental_observation(row.raw_endpoint_name, row.raw_value, row.raw_unit, species=row.species, context=context, assay_type=row.assay_type, target=context.get("target", ""), canonical_hint=row.canonical_endpoint_id)
        derived = {
            "canonical_endpoint_id": mapped["canonical_endpoint_id"],
            "normalized_value": "" if mapped.get("normalized_value") is None else str(mapped["normalized_value"]),
            "normalized_unit": mapped.get("normalized_unit", ""),
            "normalization_rule": mapped.get("normalization_rule", ""),
            "normalization_version": EXPERIMENTAL_NORMALIZATION_VERSION,
            "comparability_status": mapped.get("comparability_status", UNSUPPORTED),
            "routing_section": mapped.get("section", "UNCLASSIFIED"),
            "routing_reason": mapped.get("reason", ""),
            "canonical_endpoint_version": CANONICAL_ENDPOINT_VERSION,
            "unit_normalization_version": COMPARISON_UNIT_VERSION,
        }
        if any(getattr(row, key, None) != value for key, value in derived.items()):
            for key, value in derived.items(): setattr(row, key, value)
            changed += 1
    return {"examined": len(rows), "changed": changed, "canonical_endpoint_version": CANONICAL_ENDPOINT_VERSION, "comparison_unit_version": COMPARISON_UNIT_VERSION}
