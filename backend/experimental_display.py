"""Frozen-prediction display contract for experimental evidence.

This module never changes model outputs.  It only derives a reproducible,
scientifically bounded display representation from immutable source evidence.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Any

NORMALIZATION_VERSION = "drugopt-experimental-normalization-v1"

DIRECTLY_COMPARABLE = "DIRECTLY_COMPARABLE"
COMPARABLE_AFTER_DETERMINISTIC_CONVERSION = "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION"
CONDITIONALLY_COMPARABLE = "CONDITIONALLY_COMPARABLE"
RELATED_NOT_SAME_ENDPOINT = "RELATED_NOT_SAME_ENDPOINT"
NOT_COMPARABLE = "NOT_COMPARABLE"
UNSUPPORTED = "UNSUPPORTED"

COMPARABILITY_LABELS = {
    DIRECTLY_COMPARABLE: "Directly Comparable",
    COMPARABLE_AFTER_DETERMINISTIC_CONVERSION: "Comparable",
    CONDITIONALLY_COMPARABLE: "Condition-dependent",
    RELATED_NOT_SAME_ENDPOINT: "Related Evidence",
    NOT_COMPARABLE: "Not Directly Comparable",
    UNSUPPORTED: "Unsupported",
}


@dataclass(frozen=True)
class EndpointDisplayContract:
    endpoint_id: str
    display_label: str
    prediction_key: str
    prediction_kind: str
    canonical_display_unit: str
    canonical_scale: str
    display_precision: int
    species_requirement: str = ""
    assay_requirement: str = ""
    direction_requirement: str = ""
    accepted_experimental_measurement_types: tuple[str, ...] = ()


CONTRACTS = {
    "solubility_aqueous_logs": EndpointDisplayContract("solubility_aqueous_logs", "Solubility", "Solubility", "QUANTITATIVE", "log10(mol/L)", "LOG10", 2, assay_requirement="aqueous solubility", accepted_experimental_measurement_types=("solubility",)),
    "permeability_caco2_logpapp": EndpointDisplayContract("permeability_caco2_logpapp", "Caco-2 Papp A→B", "Permeability", "QUANTITATIVE", "log10(cm/s)", "LOG10", 2, species_requirement="Caco-2", direction_requirement="A→B", accepted_experimental_measurement_types=("caco-2", "papp")),
    "ppb_human_percent_bound": EndpointDisplayContract("ppb_human_percent_bound", "Human PPB", "Plasma protein binding", "QUANTITATIVE", "% bound", "PERCENT", 1, species_requirement="Human", assay_requirement="plasma", accepted_experimental_measurement_types=("ppb", "plasma protein binding", "fraction unbound", "fu")),
    "hlm_intrinsic_clearance_scaled_log10": EndpointDisplayContract("hlm_intrinsic_clearance_scaled_log10", "HLM intrinsic clearance", "HLM intrinsic clearance", "QUANTITATIVE", "log10(mL/min/kg)", "LOG10", 2, species_requirement="Human", assay_requirement="microsomes"),
    "rlm_intrinsic_clearance_scaled_log10": EndpointDisplayContract("rlm_intrinsic_clearance_scaled_log10", "RLM intrinsic clearance", "RLM intrinsic clearance", "QUANTITATIVE", "log10(mL/min/kg)", "LOG10", 2, species_requirement="Rat", assay_requirement="microsomes"),
    "mlm_intrinsic_clearance_scaled_log10": EndpointDisplayContract("mlm_intrinsic_clearance_scaled_log10", "MLM intrinsic clearance", "MLM intrinsic clearance", "QUANTITATIVE", "log10(mL/min/kg)", "LOG10", 2, species_requirement="Mouse", assay_requirement="microsomes"),
    "pka": EndpointDisplayContract("pka", "pKa", "pKa", "RULE_ESTIMATE", "pKa", "LINEAR", 2),
    "logd_7_4": EndpointDisplayContract("logd_7_4", "logD 7.4", "logD7.4", "DERIVED_ESTIMATE", "logD", "LINEAR", 2),
}
for isoform in ("CYP1A2", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A4"):
    CONTRACTS[f"{isoform.lower()}_inhibitor_prob"] = EndpointDisplayContract(f"{isoform.lower()}_inhibitor_prob", f"{isoform} Inhibition", f"{isoform} inhibitor", "CLASSIFICATION", "probability", "PROBABILITY", 2, species_requirement="Human")
for endpoint, label, key in (("pgp_inhibitor_prob", "P-gp Inhibition", "P-gp inhibitor"), ("herg_liability_prob", "hERG", "hERG liability"), ("ames_mutagenicity_prob", "Ames", "Ames mutagenicity"), ("dili_clinical_liability_prob", "DILI", "DILI clinical liability")):
    CONTRACTS[endpoint] = EndpointDisplayContract(endpoint, label, key, "CLASSIFICATION", "probability", "PROBABILITY", 2)


def evidence_label(origin: str | None, prediction_kind: str | None = None) -> str:
    if origin == "EXPERIMENTAL_INTERNAL": return "Experimental — Internal"
    if origin in {"EXTERNAL_CANDIDATE", "EXTERNAL_EVIDENCE_CANDIDATE"}: return "External Evidence — Candidate"
    if origin in {"EXTERNAL_IMPORTED", "IMPORTED_EXPERIMENTAL"}: return "Experimental — External — Imported"
    if origin == "RELATED_EVIDENCE": return "External Evidence — Related"
    if origin == "NEEDS_REVIEW": return "External Evidence — Needs Review"
    if origin == "EXPERIMENTAL_EXTERNAL": return "Experimental — External"
    if prediction_kind == "RULE_ESTIMATE": return "Rule Estimate"
    if prediction_kind == "DERIVED_ESTIMATE": return "Derived Estimate"
    return "Prediction" if origin else "Unavailable"


def _number(value: Any) -> float | None:
    try: return float(str(value).replace(",", "").strip())
    except (TypeError, ValueError): return None

def _contains(text: str, *terms: str) -> bool:
    low = text.lower().replace("→", "->")
    return all(term.lower().replace("→", "->") in low for term in terms)

def _result(contract: EndpointDisplayContract | None, raw_value: Any, raw_unit: str, status: str, value: float | None = None, rule: str = "", reason: str = "") -> dict:
    return {"canonical_endpoint_id": contract.endpoint_id if contract else "", "normalized_value": value,
            "normalized_unit": contract.canonical_display_unit if contract and value is not None else "",
            "normalization_rule": rule, "normalization_version": NORMALIZATION_VERSION,
            "comparability_status": status, "comparability_label": COMPARABILITY_LABELS[status],
            "reason": reason, "raw_value": str(raw_value), "raw_unit": raw_unit}

def normalize_experimental(endpoint: str, value: Any, unit: str = "", *, species: str = "", conditions: str = "", measurement_type: str = "", mw: float | None = None, target: str = "") -> dict:
    """Return display metadata; unsupported data is retained without a value transform."""
    text = " ".join((endpoint, measurement_type, conditions, target)).lower()
    unit_l = unit.lower().replace("μ", "µ")
    categorical = str(value).strip().lower()
    if "ames" in text and categorical in {"positive", "negative", "mutagenic", "non-mutagenic"}:
        return _result(CONTRACTS["ames_mutagenicity_prob"], value, unit, DIRECTLY_COMPARABLE, rule="explicit_categorical_result")
    if "herg" in text and categorical in {"blocker", "non-blocker", "positive", "negative"}:
        return _result(CONTRACTS["herg_liability_prob"], value, unit, DIRECTLY_COMPARABLE, rule="explicit_categorical_result")
    number = _number(value)
    if number is None: return _result(None, value, unit, UNSUPPORTED, reason="Non-numeric source value")
    if "pampa" in text or "mdck" in text:
        return _result(CONTRACTS["permeability_caco2_logpapp"], value, unit, NOT_COMPARABLE, reason="PAMPA and MDCK are not Caco-2 A→B Papp")
    if _contains(text, "caco-2") or _contains(text, "caco2"):
        c = CONTRACTS["permeability_caco2_logpapp"]
        if "b->a" in text or "b>a" in text or "efflux" in text or "pampa" in text or "mdck" in text:
            return _result(c, value, unit, NOT_COMPARABLE, reason="Only Caco-2 A→B Papp is eligible")
        if "a->b" not in text and "a>b" not in text:
            return _result(c, value, unit, CONDITIONALLY_COMPARABLE, reason="Caco-2 direction is not recorded")
        if "log" in unit_l: return _result(c, value, unit, DIRECTLY_COMPARABLE, number, "identity")
        if "cm/s" in unit_l and number > 0: return _result(c, value, unit, COMPARABLE_AFTER_DETERMINISTIC_CONVERSION, math.log10(number), "log10(Papp_cm_per_s)")
        return _result(c, value, unit, NOT_COMPARABLE, reason="Papp unit is not supported")
    if "solubility" in text:
        c = CONTRACTS["solubility_aqueous_logs"]
        if "log" in unit_l and ("mol" in unit_l or not unit_l):
            status = CONDITIONALLY_COMPARABLE if not any(x in text for x in ("ph", "aqueous", "thermodynamic", "kinetic")) else DIRECTLY_COMPARABLE
            return _result(c, value, unit, status, number, "identity", "Conditions incomplete" if status != DIRECTLY_COMPARABLE else "")
        factors = {"mol/l": 1., "m": 1., "mmol/l": 1e-3, "mm": 1e-3, "µmol/l": 1e-6, "umol/l": 1e-6, "µm": 1e-6, "nm": 1e-9}
        factor = next((f for u, f in factors.items() if u == unit_l.replace(" ", "")), None)
        if factor and number > 0: return _result(c, value, unit, CONDITIONALLY_COMPARABLE, math.log10(number * factor), "log10(molar_concentration)", "Solubility conditions must be reviewed")
        if unit_l in {"mg/ml", "g/l", "mg/l", "µg/ml", "ug/ml"} and mw and mw > 0 and number > 0:
            g_l = number if unit_l == "g/l" else (number if unit_l == "mg/ml" else number / 1000)
            return _result(c, value, unit, CONDITIONALLY_COMPARABLE, math.log10(g_l / mw), "mass_to_molar_using_exact_MW_then_log10", "Material form and conditions must be reviewed")
        return _result(c, value, unit, NOT_COMPARABLE, reason="Mass solubility requires exact compatible molecular weight")
    if "plasma protein" in text or re.search(r"\bppb\b|\bfu\b", text):
        c = CONTRACTS["ppb_human_percent_bound"]
        if species and species.lower() != "human": return _result(c, value, unit, NOT_COMPARABLE, reason="Human PPB cannot be compared with another species")
        if "fu" in text or "unbound" in text: return _result(c, value, unit, COMPARABLE_AFTER_DETERMINISTIC_CONVERSION, (1-number)*100 if "%" not in unit_l else 100-number, "percent_bound=100*(1-fu)")
        if "%" in unit_l: return _result(c, value, unit, DIRECTLY_COMPARABLE, number, "identity")
        if unit_l in {"fraction", "fraction bound"}: return _result(c, value, unit, COMPARABLE_AFTER_DETERMINISTIC_CONVERSION, number*100, "fraction_bound*100")
        return _result(c, value, unit, NOT_COMPARABLE, reason="PPB unit is not supported")
    if "ic50" in text or "ec50" in text or re.search(r"\bki\b|\bkd\b", text):
        return _result(None, value, unit, RELATED_NOT_SAME_ENDPOINT, reason="Potency measurement is not a model probability; IC50, EC50, Ki, and Kd remain separate")
    if "herg" in text or "cyp" in text or "p-gp" in text or "pgp" in text:
        return _result(None, value, unit, RELATED_NOT_SAME_ENDPOINT, reason="Quantitative assay evidence is related to, not numerically equivalent to, classifier probability")
    if "logd" in text:
        c = CONTRACTS["logd_7_4"]
        if "7.4" in text: return _result(c, value, unit, DIRECTLY_COMPARABLE, number, "identity")
        return _result(c, value, unit, CONDITIONALLY_COMPARABLE, reason="logD pH must be 7.4")
    if "logp" in text: return _result(CONTRACTS["logd_7_4"], value, unit, RELATED_NOT_SAME_ENDPOINT, reason="logP is not logD 7.4")
    if "pka" in text:
        c = CONTRACTS["pka"]
        return _result(c, value, unit, CONDITIONALLY_COMPARABLE if "micro" in text or "site" not in text else DIRECTLY_COMPARABLE, number, "identity", "pKa site/type requires review" if "micro" in text else "")
    return _result(None, value, unit, UNSUPPORTED, reason="No canonical endpoint mapping")

def contract_report() -> dict:
    return {"normalization_version": NORMALIZATION_VERSION, "contracts": [vars(c) | {"accepted_experimental_measurement_types": list(c.accepted_experimental_measurement_types)} for c in CONTRACTS.values()], "comparability_states": COMPARABILITY_LABELS}
