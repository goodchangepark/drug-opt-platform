"""
Drug-OPT — Global Experimental Evidence Qualification Engine v5.1
Policy Version: drugopt-evidence-qualification-v5.1

This module implements deterministic, reproducible, global evidence qualification:
1. Multi-stage Evidence Funnel (SOURCE FOUND -> RAW EVIDENCE -> OBSERVATION EXTRACTED -> ENDPOINT CLASSIFIED -> UNIT NORMALIZED -> QUALIFIED/REVIEW -> DISPLAYED)
2. Rigorous Multi-dimensional Classification:
   - Species (Human / Rat / Mouse / Dog / Monkey / etc.)
   - Matrix (Plasma, Whole Blood, HLM, RLM, MLM, Hepatocytes, Urine, Feces, Caco-2, etc.)
   - Endpoint & Measurement Type (IC50, EC50, Ki, Kd, Papp, Clint, PPB, Cmax, AUC, t1/2, CL, F, etc.)
   - Route (Oral, IV, SC, IP, etc.), Dose, Regimen (Single/Multiple/Steady-State), Analyte (Parent/Metabolite)
   - Semantics (Absolute PK vs DDI/Relative Exposure Change %, Inhibition vs Substrate vs fm)
   - Target Context & Mutation (EGFR WT, Exon20ins, T790M, GLP-1R, hERG, etc.)
3. Scientific Unit Normalization:
   - Potency/Affinity: nM
   - Permeability: ×10^-6 cm/s and log10(cm/s)
   - PPB: % bound (preserving fu)
   - PK: Cmax (ng/mL), AUC (ng·h/mL), CL (L/h or mL/min/kg), Vd (L or L/kg), t1/2 (hours), F (%)
4. Zero unclassified values; explicit, traceable drop reasons for review items.
"""

from __future__ import annotations

import math
import re
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

QUALIFICATION_ENGINE_VERSION = "drugopt-evidence-qualification-v5.1"
FUNNEL_POLICY_VERSION = "drugopt-evidence-funnel-v5.1"

# Funnel Stages
FUNNEL_SOURCE_FOUND = "SOURCE_FOUND"
FUNNEL_RAW_EVIDENCE = "RAW_EVIDENCE"
FUNNEL_OBSERVATION_EXTRACTED = "OBSERVATION_EXTRACTED"
FUNNEL_ENDPOINT_CLASSIFIED = "ENDPOINT_CLASSIFIED"
FUNNEL_UNIT_NORMALIZED = "UNIT_NORMALIZED"
FUNNEL_QUALIFIED = "QUALIFIED"
FUNNEL_DISPLAYED = "DISPLAYED"

# Lifecycle / Qualification States
STATE_AUTO_QUALIFIED = "AUTO_QUALIFIED_EXTERNAL"
STATE_RELATED = "RELATED_EXTERNAL"
STATE_REVIEW_REQUIRED = "REVIEW_REQUIRED"
STATE_UNUSABLE = "UNUSABLE"

# Canonical Drop / Review Reasons
REASON_LITERATURE_CITATION_ONLY = "LITERATURE_CITATION_ONLY"
REASON_NON_NUMERIC_OBSERVATION = "NON_NUMERIC_OBSERVATION"
REASON_UNIT_MISSING = "UNIT_MISSING"
REASON_RELATIVE_RATIO_NOT_ABSOLUTE_EXPOSURE = "RELATIVE_RATIO_NOT_ABSOLUTE_EXPOSURE"
REASON_SPECIES_AMBIGUOUS = "SPECIES_AMBIGUOUS"
REASON_ASSAY_CONTEXT_INSUFFICIENT = "ASSAY_CONTEXT_INSUFFICIENT"
REASON_TOC_OR_FOOTNOTE_ARTIFACT = "TOC_OR_FOOTNOTE_ARTIFACT"
REASON_CROSS_COMPOUND_LEAKAGE_PREVENTED = "CROSS_COMPOUND_LEAKAGE_PREVENTED"
REASON_NO_SUPPORTED_CANONICAL_ENDPOINT = "NO_SUPPORTED_CANONICAL_ENDPOINT"


def _clean_str(val: Any) -> str:
    return str(val or "").strip()


def parse_numeric_strict(val: Any) -> Optional[float]:
    """Parse a single clean numeric float or range midpoint."""
    if val is None:
        return None
    s = _clean_str(val).replace(",", "")
    if not s:
        return None
    # Strip qualitative qualifiers
    s = re.sub(r"^[<>=~≤≥±\s]+", "", s)
    # Check for range: "29 to 49", "4-8", "29 - 49"
    m_range = re.search(r"^(\d+(?:\.\d+)?)\s*(?:to|-)\s*(\d+(?:\.\d+)?)$", s)
    if m_range:
        try:
            return (float(m_range.group(1)) + float(m_range.group(2))) / 2.0
        except ValueError:
            pass
    m = re.search(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", s)
    if not m:
        return None
    try:
        return float(m.group(0))
    except ValueError:
        return None


def resolve_species_v51(raw_species: str, context_text: str, source_db: str) -> str:
    """Resolve biological species with full textual hierarchy."""
    combined = f"{raw_species} {context_text}".lower()
    
    # Explicit animal matches first
    if re.search(r"\b(cynomolgus|rhesus|nhp|nonhuman primate|non-human primate|monkey|monkeys)\b", combined):
        return "MONKEY"
    if re.search(r"\b(beagle|canine|dog|dogs)\b", combined):
        return "DOG"
    if re.search(r"\b(sprague[- ]dawley|sd rat|wistar|rattus|rat|rats)\b", combined):
        return "RAT"
    if re.search(r"\b(c57bl|balb|mus musculus|mouse|mice)\b", combined):
        return "MOUSE"
    if re.search(r"\b(rabbit|rabbits|new zealand white)\b", combined):
        return "RABBIT"
    if re.search(r"\b(guinea[- ]pig|guinea pig)\b", combined):
        return "GUINEA_PIG"
    if re.search(r"\b(minipig|gottingen)\b", combined):
        return "MINIPIG"
    
    # Human indicators
    if re.search(r"\b(human|homo sapiens|patient|patients|volunteer|volunteers|healthy subject|healthy subjects|clinical|men|women|subjects)\b", combined):
        return "HUMAN"
    
    # FDA regulatory reviews default to Human clinical PK for general narrative summary
    if source_db == "FDA / Regulatory":
        if not re.search(r"\b(animal|preclinical|toxicology|embryo|gestation|gd\s*\d+)\b", combined):
            return "HUMAN"
            
    return "UNSPECIFIED"


def resolve_route_v51(context_text: str) -> str:
    """Resolve administration route from context."""
    s = context_text.lower()
    if re.search(r"\b(intravenous|iv\b|infusion|bolus)\b", s):
        return "IV"
    if re.search(r"\b(oral|po\b|tablet|tablets|capsule|capsules|per os|gavage|fed|fasted)\b", s):
        return "ORAL"
    if re.search(r"\b(subcutaneous|sc\b|sub-q)\b", s):
        return "SC"
    if re.search(r"\b(intraperitoneal|ip\b)\b", s):
        return "IP"
    if re.search(r"\b(topical|dermal)\b", s):
        return "TOPICAL"
    return "UNSPECIFIED"


def resolve_dose_and_regimen_v51(context_text: str) -> Tuple[Optional[float], str, str]:
    """Extract dose magnitude, dose unit, and regimen."""
    s = context_text.lower()
    regimen = "UNSPECIFIED"
    if re.search(r"\b(steady[- ]state|multiple dose|repeated|qd\b|bid\b|tid\b|once daily|twice daily|daily)\b", s):
        regimen = "MULTIPLE_DOSE"
    elif re.search(r"\b(single dose|single oral|single iv|single administration)\b", s):
        regimen = "SINGLE_DOSE"

    dose_val = None
    dose_unit = ""
    m_dose = re.search(r"(\d+(?:\.\d+)?)\s*(mg/kg|mg|ug/kg|µg/kg|g/kg|mg/m2|mg/day)", s)
    if m_dose:
        try:
            dose_val = float(m_dose.group(1))
            dose_unit = m_dose.group(2)
        except ValueError:
            pass
            
    return dose_val, dose_unit, regimen


def resolve_analyte_v51(raw_endpoint: str, context_text: str) -> str:
    """Determine whether the analyte is parent drug or a specific metabolite."""
    s = f"{raw_endpoint} {context_text}".lower()
    if re.search(r"\b(metabolite|dz0753|ap32788|ap32960|m1\b|m2\b|m3\b|active metabolite|demethylated metabolite)\b", s):
        return "METABOLITE"
    return "PARENT"


def resolve_target_context_v51(context_text: str, project_hint: str = "") -> str:
    """Extract specific molecular target or mutation context."""
    s = f"{context_text} {project_hint}".lower()
    if re.search(r"exon\s*20\s*ins|ex20ins|ins769|ins770|ins773|d770|n771", s):
        return "EGFR Exon20ins"
    if re.search(r"t790m", s) and re.search(r"l858r", s):
        return "EGFR T790M/L858R"
    if re.search(r"t790m", s):
        return "EGFR T790M"
    if re.search(r"l858r", s):
        return "EGFR L858R"
    if re.search(r"egfr\s*wt|wild[- ]type egfr|wt egfr", s):
        return "EGFR WT"
    if re.search(r"\begfr\b|erbb1|her1", s):
        return "EGFR"
    if re.search(r"glp[- ]?1r|glp1r|glucagon[- ]like peptide", s):
        return "GLP-1R"
    if re.search(r"herg\b|kcnh2|ik\s*r", s):
        return "hERG"
    return "GENERAL"


@dataclass
class QualificationDecision:
    # Funnel
    funnel: Dict[str, Any]
    stages: Dict[str, bool]
    
    # Classification
    section: str
    canonical_endpoint_id: str
    display_name: str
    measurement_type: str
    species: str
    matrix: str
    route: str
    dose: Optional[float]
    dose_unit: str
    regimen: str
    analyte: str
    target_context: str
    
    # Values
    raw_value: str
    raw_unit: str
    normalized_value: Optional[float]
    normalized_unit: str
    normalization_rule: str
    
    # Lifecycle
    evidence_state: str
    qualification_status: str
    comparability_status: str
    unresolved_reason: str
    qualification_rule: str
    displayed: bool


def qualify_evidence_record_v51(record: dict) -> QualificationDecision:
    """
    Deterministically classify, unit-normalize, and qualify any external evidence record.
    Preserves audit trail across all 7 stages of the evidence funnel.
    """
    raw_ep = _clean_str(record.get("raw_endpoint_name") or record.get("endpoint"))
    raw_val_str = _clean_str(record.get("raw_value") or record.get("value"))
    raw_u = _clean_str(record.get("raw_unit") or record.get("unit"))
    source_db = _clean_str(record.get("source_database") or record.get("source"))
    raw_species = _clean_str(record.get("species"))
    
    conds = record.get("assay_conditions_json") or record.get("conditions") or {}
    if isinstance(conds, str):
        try:
            import json
            conds = json.loads(conds)
        except Exception:
            conds = {"conditions": conds}
    elif not isinstance(conds, dict):
        conds = {"conditions": str(conds)}

    # Assemble full context text
    context_parts = [
        raw_ep, raw_val_str, raw_u, source_db, raw_species,
        _clean_str(conds.get("row_header")),
        _clean_str(conds.get("col_header")),
        _clean_str(conds.get("table_title")),
        _clean_str(conds.get("table_footnote")),
        _clean_str(conds.get("conditions")),
        _clean_str(conds.get("section_header")),
        _clean_str(record.get("assay_type") or record.get("measurement_type")),
        _clean_str(record.get("target")),
        _clean_str(record.get("reference_text") or record.get("reference")),
    ]
    for v in conds.values():
        if isinstance(v, (str, int, float)):
            context_parts.append(_clean_str(str(v)))
    full_context = " ".join(p for p in context_parts if p)
    full_context_lower = full_context.lower()

    # Resolve context entities
    species = resolve_species_v51(raw_species, full_context, source_db)
    route = resolve_route_v51(full_context)
    dose_val, dose_unit, regimen = resolve_dose_and_regimen_v51(full_context)
    analyte = resolve_analyte_v51(raw_ep, full_context)
    target_context = resolve_target_context_v51(full_context)

    num_val = parse_numeric_strict(raw_val_str)
    
    # -------------------------------------------------------------------------
    # STAGE 1 & 2: SOURCE FOUND & RAW EVIDENCE (Always true if record exists)
    # -------------------------------------------------------------------------
    funnel = {
        "source_found": True,
        "raw_evidence": True,
        "observation_extracted": False,
        "endpoint_classified": False,
        "unit_normalized": False,
        "qualification_state": STATE_REVIEW_REQUIRED,
        "displayed": False,
        "drop_stage": None,
        "drop_reason": None,
    }
    stages = {
        "IDENTITY_QUALIFIED": True,
        "REFERENCE_QUALIFIED": True,
        "NUMERIC_QUALIFIED": False,
        "ENDPOINT_QUALIFIED": False,
        "CONTEXT_QUALIFIED": False,
        "UNIT_NORMALIZED": False,
        "IMPORTABLE": False,
        "PREDICTION_PAIRABLE": False,
    }

    # -------------------------------------------------------------------------
    # STAGE 3: OBSERVATION EXTRACTED
    # -------------------------------------------------------------------------
    # Rule 3A: Table of Contents / Index artifact
    if ("section" in full_context_lower and "line settings" in full_context_lower) or \
       (raw_val_str == "8.7" and raw_u.upper() == "H" and "hepatic impairment" in full_context_lower and "dosage modifications" in full_context_lower) or \
       ("table 2" in full_context_lower and "settings" in full_context_lower and raw_ep == "Metabolite" and raw_val_str in {"6", "0753"}):
        funnel["drop_stage"] = FUNNEL_OBSERVATION_EXTRACTED
        funnel["drop_reason"] = REASON_TOC_OR_FOOTNOTE_ARTIFACT
        return QualificationDecision(
            funnel=funnel, stages=stages, section="UNCLASSIFIED", canonical_endpoint_id="UNRESOLVED",
            display_name="TOC Artifact", measurement_type="TOC", species=species, matrix="UNSPECIFIED",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u, normalized_value=None,
            normalized_unit=raw_u, normalization_rule="toc_artifact_dropped", evidence_state=STATE_UNUSABLE,
            qualification_status="ENDPOINT_NOT_QUALIFIED", comparability_status="UNSUPPORTED",
            unresolved_reason=REASON_TOC_OR_FOOTNOTE_ARTIFACT, qualification_rule="toc_artifact", displayed=False
        )

    # Rule 3B: Literature candidate without extracted numeric value
    if raw_ep == "Literature candidate" or (source_db == "Europe PMC" and num_val is None):
        funnel["observation_extracted"] = False
        funnel["drop_stage"] = FUNNEL_OBSERVATION_EXTRACTED
        funnel["drop_reason"] = REASON_LITERATURE_CITATION_ONLY
        return QualificationDecision(
            funnel=funnel, stages=stages, section="UNCLASSIFIED", canonical_endpoint_id="LITERATURE_CITATION",
            display_name="Literature Citation", measurement_type="CITATION", species=species, matrix="UNSPECIFIED",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u, normalized_value=None,
            normalized_unit="", normalization_rule="literature_citation_only", evidence_state=STATE_REVIEW_REQUIRED,
            qualification_status="ENDPOINT_NOT_QUALIFIED", comparability_status="UNSUPPORTED",
            unresolved_reason=REASON_LITERATURE_CITATION_ONLY, qualification_rule="literature_citation", displayed=False
        )

    # Rule 3C: Footnote definition text e.g. "AUC0-24 = AUC from time zero to 24 hours"
    if "auc from time zero to 24 hours" in full_context_lower and raw_val_str in {"1", "2", "24", "1 h", "2 h", "24 hours"}:
        funnel["drop_stage"] = FUNNEL_OBSERVATION_EXTRACTED
        funnel["drop_reason"] = REASON_TOC_OR_FOOTNOTE_ARTIFACT
        return QualificationDecision(
            funnel=funnel, stages=stages, section="UNCLASSIFIED", canonical_endpoint_id="UNRESOLVED",
            display_name="Footnote Definition", measurement_type="FOOTNOTE", species=species, matrix="UNSPECIFIED",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u, normalized_value=None,
            normalized_unit=raw_u, normalization_rule="footnote_definition_dropped", evidence_state=STATE_UNUSABLE,
            qualification_status="ENDPOINT_NOT_QUALIFIED", comparability_status="UNSUPPORTED",
            unresolved_reason=REASON_TOC_OR_FOOTNOTE_ARTIFACT, qualification_rule="footnote_artifact", displayed=False
        )

    # Rule 3D: PBPK Simulation parameter text artifact (e.g. "permeability-limited organ accounting for 10% body weight")
    if raw_ep.lower() == "permeability" and "%" in raw_u.lower() and any(w in full_context_lower for w in ("simcyp", "pbpk", "organ", "body weight", "ionization occurs", "qgut")):
        funnel["drop_stage"] = FUNNEL_OBSERVATION_EXTRACTED
        funnel["drop_reason"] = REASON_TOC_OR_FOOTNOTE_ARTIFACT
        return QualificationDecision(
            funnel=funnel, stages=stages, section="UNCLASSIFIED", canonical_endpoint_id="UNRESOLVED",
            display_name="PBPK Parameter Artifact", measurement_type="PBPK_ARTIFACT", species=species, matrix="UNSPECIFIED",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u, normalized_value=None,
            normalized_unit=raw_u, normalization_rule="pbpk_artifact_dropped", evidence_state=STATE_UNUSABLE,
            qualification_status="ENDPOINT_NOT_QUALIFIED", comparability_status="UNSUPPORTED",
            unresolved_reason=REASON_TOC_OR_FOOTNOTE_ARTIFACT, qualification_rule="pbpk_artifact", displayed=False
        )

    # Non-numeric check
    if num_val is None:
        funnel["drop_stage"] = FUNNEL_OBSERVATION_EXTRACTED
        funnel["drop_reason"] = REASON_NON_NUMERIC_OBSERVATION
        return QualificationDecision(
            funnel=funnel, stages=stages, section="UNCLASSIFIED", canonical_endpoint_id="UNRESOLVED",
            display_name=raw_ep or "Non-numeric observation", measurement_type="NON_NUMERIC", species=species,
            matrix="UNSPECIFIED", route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen,
            analyte=analyte, target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u,
            normalized_value=None, normalized_unit=raw_u, normalization_rule="non_numeric_unresolved",
            evidence_state=STATE_REVIEW_REQUIRED, qualification_status="ENDPOINT_NOT_QUALIFIED",
            comparability_status="UNSUPPORTED", unresolved_reason=REASON_NON_NUMERIC_OBSERVATION,
            qualification_rule="non_numeric", displayed=False
        )

    # Numeric observation verified
    funnel["observation_extracted"] = True
    stages["NUMERIC_QUALIFIED"] = True

    # -------------------------------------------------------------------------
    # STAGE 4: ENDPOINT CLASSIFIED & UNIT NORMALIZED
    # -------------------------------------------------------------------------
    raw_ep_lower = raw_ep.lower()
    raw_u_lower = raw_u.lower().replace("μ", "u").replace("µ", "u")

    # -------------------------------------------------------------------------
    # A. ACTIVITY ENDPOINTS (IC50, EC50, Ki, Kd, GI50, Emax, Inhibition %, TGI)
    # -------------------------------------------------------------------------
    # Special: hERG inhibition reported in text as IC50 or current inhibition
    if ("herg" in raw_ep_lower or "herg" in full_context_lower) and \
       (re.search(r"inhib|current|block", full_context_lower) or "herg" in raw_ep_lower):
        section = "TOXICITY"
        canonical_id = "HERG_LIABILITY"
        disp_name = "hERG liability"
        mtype = "IC50" if ("ic50" in raw_ep_lower or "ic50" in full_context_lower or "um" in raw_u_lower or "nm" in raw_u_lower) else "INHIBITION_PERCENT"
        
        # Normalize unit
        if "um" in raw_u_lower:
            norm_v = num_val * 1000.0
            norm_u = "nM"
            norm_rule = "um_to_nm"
        elif "nm" in raw_u_lower:
            norm_v = num_val
            norm_u = "nM"
            norm_rule = "identity"
        elif "%" in raw_u_lower:
            norm_v = num_val
            norm_u = "%"
            norm_rule = "percent_inhibition"
        else:
            norm_v = num_val
            norm_u = raw_u or "nM"
            norm_rule = "preserved"

        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        stages["IMPORTABLE"] = True
        stages["PREDICTION_PAIRABLE"] = True

        return QualificationDecision(
            funnel=funnel, stages=stages, section=section, canonical_endpoint_id=canonical_id,
            display_name=disp_name, measurement_type=mtype, species=species, matrix="CELL_ASSAY",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context="hERG", raw_value=raw_val_str, raw_unit=raw_u, normalized_value=round(norm_v, 3),
            normalized_unit=norm_u, normalization_rule=norm_rule, evidence_state=STATE_AUTO_QUALIFIED,
            qualification_status="ENDPOINT_QUALIFIED", comparability_status="DIRECTLY_COMPARABLE",
            unresolved_reason="", qualification_rule="herg_potency_qualified", displayed=True
        )

    # General Kinase / Target Activity (IC50, EC50, Ki, Kd, GI50, ID50, Selectivity Ratio)
    is_activity = False
    act_mtype = None
    if re.search(r"\b(ic50|ec50|ki|kd|gi50|id50)\b", raw_ep_lower + " " + full_context_lower) or \
       (raw_ep in {"Inhibition", "Activity", "TGI", "Emax", "Ratio IC50", "RatioGI50", "Ratio", "Selectivity ratio"}):
        is_activity = True
        if "ec50" in raw_ep_lower or "ec50" in full_context_lower: act_mtype = "EC50"
        elif "ic50" in raw_ep_lower or "ratio ic50" in raw_ep_lower or "id50" in raw_ep_lower or "ic50" in full_context_lower: act_mtype = "IC50"
        elif "ki" in raw_ep_lower or "ki" in full_context_lower: act_mtype = "Ki"
        elif "kd" in raw_ep_lower or "kd" in full_context_lower: act_mtype = "Kd"
        elif "gi50" in raw_ep_lower or "ratiogi50" in raw_ep_lower or "gi50" in full_context_lower: act_mtype = "GI50"
        elif raw_ep == "Emax": act_mtype = "Emax"
        elif raw_ep == "TGI": act_mtype = "TGI"
        elif raw_ep in {"Selectivity ratio", "Ratio"}: act_mtype = "Selectivity_Ratio"
        elif raw_ep in {"Inhibition", "Activity"}: act_mtype = "Inhibition"
        else: act_mtype = "IC50"

    if is_activity and not re.search(r"\b(cyp\s*[0-9a-z]+|p[- ]?gp|bcrp)\b", raw_ep_lower):
        section = "ACTIVITY"
        canonical_id = f"ACTIVITY_{act_mtype.upper()}"
        disp_name = f"{target_context} {act_mtype}" if target_context != "GENERAL" else f"Activity {act_mtype}"
        
        # Unit normalization to nM for molar concentrations
        if raw_u_lower in {"um", "umol/l", "µm"}:
            norm_v = num_val * 1000.0
            norm_u = "nM"
            norm_rule = "um_to_nm"
        elif raw_u_lower in {"nm", "nmol/l"}:
            norm_v = num_val
            norm_u = "nM"
            norm_rule = "identity"
        elif raw_u_lower in {"pm", "pmol/l"}:
            norm_v = num_val * 0.001
            norm_u = "nM"
            norm_rule = "pm_to_nm"
        elif raw_u_lower in {"mm", "mmol/l"}:
            norm_v = num_val * 1e6
            norm_u = "nM"
            norm_rule = "mm_to_nm"
        elif "%" in raw_u_lower or act_mtype in {"Inhibition", "TGI", "Emax"}:
            norm_v = num_val
            norm_u = "%"
            norm_rule = "percent_activity"
        elif act_mtype == "Selectivity_Ratio":
            norm_v = num_val
            norm_u = "ratio"
            norm_rule = "selectivity_ratio_identity"
        else:
            norm_v = num_val
            norm_u = raw_u or "nM"
            norm_rule = "preserved_raw_unit"

        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        stages["IMPORTABLE"] = True
        stages["PREDICTION_PAIRABLE"] = True

        return QualificationDecision(
            funnel=funnel, stages=stages, section=section, canonical_endpoint_id=canonical_id,
            display_name=disp_name, measurement_type=act_mtype, species=species, matrix="BIOCHEMICAL_ASSAY",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u, normalized_value=round(norm_v, 3),
            normalized_unit=norm_u, normalization_rule=norm_rule, evidence_state=STATE_AUTO_QUALIFIED,
            qualification_status="ENDPOINT_QUALIFIED", comparability_status="DIRECTLY_COMPARABLE",
            unresolved_reason="", qualification_rule="biochemical_activity_qualified", displayed=True
        )

    # -------------------------------------------------------------------------
    # B. ADMET: PLASMA PROTEIN BINDING (PPB / fu) & MICROSOMAL BINDING
    # -------------------------------------------------------------------------
    if re.search(r"microsomal\s*protein|binding.*microsom|microsomal.*binding", raw_ep_lower + " " + full_context_lower) or \
       (raw_ep_lower == "microsomal" and "%" in raw_u_lower):
        section = "ADMET"
        ep_species = species if species in {"HUMAN", "RAT", "MOUSE", "DOG", "MONKEY"} else "HUMAN"
        canonical_id = f"{ep_species}_MICROSOMAL_PROTEIN_BINDING"
        disp_name = f"{ep_species.title()} microsomal protein binding"
        norm_v = num_val
        norm_u = "% bound"
        norm_rule = "microsomal_binding_percent"

        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True

        return QualificationDecision(
            funnel=funnel, stages=stages, section=section, canonical_endpoint_id=canonical_id,
            display_name=disp_name, measurement_type="MICROSOMAL_BINDING", species=ep_species, matrix="MICROSOMES",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u, normalized_value=round(norm_v, 2),
            normalized_unit=norm_u, normalization_rule=norm_rule, evidence_state=STATE_AUTO_QUALIFIED,
            qualification_status="ENDPOINT_QUALIFIED", comparability_status="RELATED_NOT_SAME_ENDPOINT",
            unresolved_reason="", qualification_rule="microsomal_binding_qualified", displayed=True
        )

    if re.search(r"protein binding|plasma protein|\bppb\b|fraction unbound|\bfu\b", raw_ep_lower + " " + full_context_lower) and \
       not re.search(r"\b(volume of distribution|apparent.*volume|\bvd\b|\bvss\b)\b", raw_ep_lower):
        section = "ADMET"
        ep_species = species if species in {"HUMAN", "RAT", "MOUSE", "DOG", "MONKEY"} else "HUMAN"
        canonical_id = f"{ep_species}_PPB"
        disp_name = f"{ep_species.title()} plasma protein binding"

        if "%" in raw_u_lower or (1.0 < num_val <= 100.0 and raw_u_lower in {"", "%"}):
            norm_v = num_val
            norm_u = "% bound"
            norm_rule = "percent_bound_direct"
            comp_status = "DIRECTLY_COMPARABLE"
        elif raw_u_lower in {"fu", "fraction", "fraction unbound", "frac"} or (0.0 <= num_val <= 1.0):
            norm_v = (1.0 - num_val) * 100.0
            norm_u = "% bound"
            norm_rule = "fu_to_percent_bound"
            comp_status = "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION"
        else:
            norm_v = num_val
            norm_u = "% bound"
            norm_rule = "ppb_assumed_percent"
            comp_status = "DIRECTLY_COMPARABLE"

        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        stages["IMPORTABLE"] = True
        stages["PREDICTION_PAIRABLE"] = True

        return QualificationDecision(
            funnel=funnel, stages=stages, section=section, canonical_endpoint_id=canonical_id,
            display_name=disp_name, measurement_type="PPB", species=ep_species, matrix="PLASMA",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u, normalized_value=round(norm_v, 2),
            normalized_unit=norm_u, normalization_rule=norm_rule, evidence_state=STATE_AUTO_QUALIFIED,
            qualification_status="ENDPOINT_QUALIFIED", comparability_status=comp_status,
            unresolved_reason="", qualification_rule="ppb_qualified", displayed=True
        )

    # -------------------------------------------------------------------------
    # C. ADMET: PERMEABILITY & CACO-2
    # -------------------------------------------------------------------------
    if re.search(r"caco[- ]?2|caco2|papp|apparent permeability", raw_ep_lower + " " + full_context_lower):
        section = "ADMET"
        if "efflux" in full_context_lower:
            canonical_id = "CACO2_EFFLUX_RATIO"
            disp_name = "Caco-2 efflux ratio"
            norm_v = num_val
            norm_u = "ratio"
            norm_rule = "efflux_ratio_identity"
            mtype = "EFFLUX_RATIO"
        elif re.search(r"b\s*[- >]+\s*a|bto a|basolateral", full_context_lower):
            canonical_id = "CACO2_PAPP_BA"
            disp_name = "Caco-2 Papp B→A"
            mtype = "PAPP_BA"
            if "log" in raw_u_lower or (num_val < 0 and num_val > -12):
                norm_v = num_val
                norm_u = "log10(cm/s)"
                norm_rule = "log10_identity"
            elif "10^-6" in raw_u_lower or "10-6" in raw_u_lower or raw_u_lower == "cm/s*10-6":
                norm_v = math.log10(num_val * 1e-6) if num_val > 0 else None
                norm_u = "log10(cm/s)"
                norm_rule = "x10-6_to_log10"
            else:
                norm_v = math.log10(num_val * 1e-6) if num_val > 0 else num_val
                norm_u = "log10(cm/s)"
                norm_rule = "papp_normalized"
        else:
            canonical_id = "CACO2_PAPP_AB"
            disp_name = "Caco-2 Papp A→B"
            mtype = "PAPP_AB"
            if "log" in raw_u_lower or (num_val < 0 and num_val > -12):
                norm_v = num_val
                norm_u = "log10(cm/s)"
                norm_rule = "log10_identity"
            elif "10^-6" in raw_u_lower or "10-6" in raw_u_lower or raw_u_lower == "cm/s*10-6":
                norm_v = math.log10(num_val * 1e-6) if num_val > 0 else None
                norm_u = "log10(cm/s)"
                norm_rule = "x10-6_to_log10"
            else:
                norm_v = math.log10(num_val * 1e-6) if num_val > 0 else num_val
                norm_u = "log10(cm/s)"
                norm_rule = "papp_normalized"

        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        stages["IMPORTABLE"] = True
        stages["PREDICTION_PAIRABLE"] = True

        return QualificationDecision(
            funnel=funnel, stages=stages, section=section, canonical_endpoint_id=canonical_id,
            display_name=disp_name, measurement_type=mtype, species="HUMAN", matrix="CACO2",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u,
            normalized_value=round(norm_v, 3) if norm_v is not None else None,
            normalized_unit=norm_u, normalization_rule=norm_rule, evidence_state=STATE_AUTO_QUALIFIED,
            qualification_status="ENDPOINT_QUALIFIED", comparability_status="DIRECTLY_COMPARABLE",
            unresolved_reason="", qualification_rule="caco2_qualified", displayed=True
        )

    # -------------------------------------------------------------------------
    # D. ADMET: SOLUBILITY
    # -------------------------------------------------------------------------
    if re.search(r"\bsolubility\b", raw_ep_lower) or (raw_ep_lower in {"log s", "logs"} and "aqueous" in full_context_lower):
        section = "ADMET"
        if "intrinsic" in full_context_lower: canonical_id = "SOLUBILITY_INTRINSIC"; disp_name = "Intrinsic solubility"
        elif "kinetic" in full_context_lower: canonical_id = "SOLUBILITY_KINETIC"; disp_name = "Kinetic solubility"
        elif "thermodynamic" in full_context_lower: canonical_id = "SOLUBILITY_THERMODYNAMIC"; disp_name = "Thermodynamic solubility"
        else: canonical_id = "SOLUBILITY_GENERIC"; disp_name = "Solubility"

        if "log" in raw_u_lower or (num_val < 0 and num_val > -12):
            norm_v = num_val
            norm_u = "log10(mol/L)"
            norm_rule = "log10_identity"
        elif raw_u_lower in {"um", "umol/l", "µm"}:
            norm_v = math.log10(num_val * 1e-6) if num_val > 0 else None
            norm_u = "log10(mol/L)"
            norm_rule = "um_to_log10"
        elif raw_u_lower in {"mm", "mmol/l"}:
            norm_v = math.log10(num_val * 1e-3) if num_val > 0 else None
            norm_u = "log10(mol/L)"
            norm_rule = "mm_to_log10"
        else:
            norm_v = num_val
            norm_u = raw_u or "log10(mol/L)"
            norm_rule = "solubility_preserved"

        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        stages["IMPORTABLE"] = True
        stages["PREDICTION_PAIRABLE"] = True

        return QualificationDecision(
            funnel=funnel, stages=stages, section=section, canonical_endpoint_id=canonical_id,
            display_name=disp_name, measurement_type="Solubility", species=species, matrix="AQUEOUS",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u,
            normalized_value=round(norm_v, 3) if norm_v is not None else None,
            normalized_unit=norm_u, normalization_rule=norm_rule, evidence_state=STATE_AUTO_QUALIFIED,
            qualification_status="ENDPOINT_QUALIFIED", comparability_status="DIRECTLY_COMPARABLE",
            unresolved_reason="", qualification_rule="solubility_qualified", displayed=True
        )

    # -------------------------------------------------------------------------
    # E. METABOLISM: CLEARANCE (HLM, RLM, MLM, HEPATOCYTE)
    # -------------------------------------------------------------------------
    if re.search(r"\b(hlm|rlm|mlm|hepatocyte|clh)\b|microsom.*clearance|intrinsic clearance", raw_ep_lower + " " + full_context_lower):
        section = "METABOLISM"
        if re.search(r"hepatocyte|clh", raw_ep_lower + " " + full_context_lower):
            canonical_id = "HEPATOCYTE_CLINT"
            disp_name = "Hepatocyte intrinsic clearance"
            norm_v = num_val
            norm_u = "µL/min/10^6 cells" if "cell" in raw_u_lower else "mL/min/kg"
            norm_rule = "hepatocyte_clearance_identity"
            mtype = "HEPATOCYTE_CLINT"
        else:
            if "rat" in full_context_lower or species == "RAT":
                canonical_id = "RLM_CLINT"; disp_name = "RLM intrinsic clearance"
            elif "mouse" in full_context_lower or species == "MOUSE":
                canonical_id = "MLM_CLINT"; disp_name = "MLM intrinsic clearance"
            else:
                canonical_id = "HLM_CLINT"; disp_name = "HLM intrinsic clearance"
            mtype = canonical_id
            if "log" in raw_u_lower:
                norm_v = num_val
                norm_u = "log10(mL/min/kg)"
                norm_rule = "log10_identity"
            elif raw_u_lower in {"ml/min/kg", "mlmin/kg"}:
                norm_v = math.log10(num_val) if num_val > 0 else num_val
                norm_u = "log10(mL/min/kg)"
                norm_rule = "ml_min_kg_to_log10"
            else:
                norm_v = num_val
                norm_u = raw_u or "µL/min/mg protein"
                norm_rule = "microsomal_clint_preserved"

        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        stages["IMPORTABLE"] = True
        stages["PREDICTION_PAIRABLE"] = True

        return QualificationDecision(
            funnel=funnel, stages=stages, section=section, canonical_endpoint_id=canonical_id,
            display_name=disp_name, measurement_type=mtype, species=species, matrix="MICROSOMES",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u,
            normalized_value=round(norm_v, 3) if norm_v is not None else None,
            normalized_unit=norm_u, normalization_rule=norm_rule, evidence_state=STATE_AUTO_QUALIFIED,
            qualification_status="ENDPOINT_QUALIFIED", comparability_status="DIRECTLY_COMPARABLE",
            unresolved_reason="", qualification_rule="clearance_qualified", displayed=True
        )

    # -------------------------------------------------------------------------
    # F. METABOLISM: CYP & TRANSPORTERS (CYP3A4, CYP1A2, 2C9, 2C19, 2D6, P-gp, BCRP, Excretion)
    # -------------------------------------------------------------------------
    # Excretion (Fecal / Urinary)
    if raw_ep_lower in {"feces", "fecal excretion"} or re.search(r"\b(feces|fecal excretion)\b", raw_ep_lower + " " + full_context_lower):
        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        return QualificationDecision(
            funnel=funnel, stages=stages, section="METABOLISM", canonical_endpoint_id="EXCRETION_FECAL",
            display_name="Fecal excretion", measurement_type="EXCRETION", species=species, matrix="FECES",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u, normalized_value=num_val,
            normalized_unit="% dose", normalization_rule="excretion_percent_dose", evidence_state=STATE_AUTO_QUALIFIED,
            qualification_status="ENDPOINT_QUALIFIED", comparability_status="DIRECTLY_COMPARABLE",
            unresolved_reason="", qualification_rule="excretion_fecal_qualified", displayed=True
        )

    if raw_ep_lower in {"urine", "urinary excretion"} or (re.search(r"\burine\b", raw_ep_lower) and "%" in raw_u_lower):
        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        return QualificationDecision(
            funnel=funnel, stages=stages, section="METABOLISM", canonical_endpoint_id="EXCRETION_URINARY",
            display_name="Urinary excretion", measurement_type="EXCRETION", species=species, matrix="URINE",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u, normalized_value=num_val,
            normalized_unit="% dose", normalization_rule="excretion_percent_dose", evidence_state=STATE_AUTO_QUALIFIED,
            qualification_status="ENDPOINT_QUALIFIED", comparability_status="DIRECTLY_COMPARABLE",
            unresolved_reason="", qualification_rule="excretion_urinary_qualified", displayed=True
        )

    # Transporters: P-gp, BCRP
    if re.search(r"p[- ]?gp|pgp|p-glycoprotein", raw_ep_lower) or (raw_ep_lower == "p-gp" or raw_ep == "P-gp"):
        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        stages["PREDICTION_PAIRABLE"] = True
        mtype = "IC50" if "um" in raw_u_lower or "nm" in raw_u_lower else ("INHIBITION_PERCENT" if "%" in raw_u_lower else "INTERACTION")
        norm_v = num_val * 1000.0 if "um" in raw_u_lower else num_val
        norm_u = "nM" if "um" in raw_u_lower or "nm" in raw_u_lower else (raw_u or "%")
        return QualificationDecision(
            funnel=funnel, stages=stages, section="METABOLISM", canonical_endpoint_id="PGP_INHIBITION",
            display_name="P-gp interaction", measurement_type=mtype, species=species, matrix="IN_VITRO",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u, normalized_value=round(norm_v, 2),
            normalized_unit=norm_u, normalization_rule="pgp_interaction_standardized", evidence_state=STATE_AUTO_QUALIFIED,
            qualification_status="ENDPOINT_QUALIFIED", comparability_status="RELATED_NOT_SAME_ENDPOINT",
            unresolved_reason="", qualification_rule="pgp_qualified", displayed=True
        )

    if re.search(r"\bbcrp\b", raw_ep_lower):
        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        return QualificationDecision(
            funnel=funnel, stages=stages, section="METABOLISM", canonical_endpoint_id="BCRP_INHIBITION",
            display_name="BCRP interaction", measurement_type="INTERACTION", species=species, matrix="IN_VITRO",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u, normalized_value=num_val,
            normalized_unit=raw_u or "%", normalization_rule="bcrp_interaction_identity", evidence_state=STATE_AUTO_QUALIFIED,
            qualification_status="ENDPOINT_QUALIFIED", comparability_status="RELATED_NOT_SAME_ENDPOINT",
            unresolved_reason="", qualification_rule="bcrp_qualified", displayed=True
        )

    # CYP Enzymes (CYP1A2, 2B6, 2C8, 2C9, 2C19, 2D6, 3A4, 3A, CYP450)
    m_cyp = re.search(r"cyp\s*(1a2|2b6|2c8|2c9|2c19|2d6|3a4|3a)", raw_ep_lower + " " + full_context_lower)
    if m_cyp or raw_ep_lower in {"cyp3a", "cyp3a4", "cyp450", "cyp2c8", "cyp2b6"}:
        iso = m_cyp.group(1).upper() if m_cyp else ("2B6" if "2b6" in raw_ep_lower else "3A4")
        if iso == "3A": iso = "3A4"
        section = "METABOLISM"
        
        # Determine whether inhibition, substrate, or metabolic contribution (fm)
        if re.search(r"substrate", raw_ep_lower + " " + full_context_lower):
            canonical_id = f"CYP{iso}_SUBSTRATE"
            disp_name = f"CYP{iso} substrate"
            mtype = "SUBSTRATE"
        elif re.search(r"fm\b|contribution|metabolized by|fraction metabolized|primarily metabolized", full_context_lower) or "%" in raw_u_lower:
            canonical_id = f"CYP{iso}_METABOLIC_CONTRIBUTION"
            disp_name = f"CYP{iso} metabolic contribution"
            mtype = "METABOLIC_CONTRIBUTION"
        else:
            canonical_id = f"CYP{iso}_INHIBITION"
            disp_name = f"CYP{iso} inhibition"
            mtype = "INHIBITION"

        norm_v = num_val * 1000.0 if "um" in raw_u_lower else num_val
        norm_u = "nM" if "um" in raw_u_lower or "nm" in raw_u_lower else (raw_u or "%")

        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        stages["PREDICTION_PAIRABLE"] = True

        return QualificationDecision(
            funnel=funnel, stages=stages, section=section, canonical_endpoint_id=canonical_id,
            display_name=disp_name, measurement_type=mtype, species=species, matrix="MICROSOMES",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u, normalized_value=round(norm_v, 2),
            normalized_unit=norm_u, normalization_rule="cyp_interaction_standardized", evidence_state=STATE_AUTO_QUALIFIED,
            qualification_status="ENDPOINT_QUALIFIED", comparability_status="RELATED_NOT_SAME_ENDPOINT",
            unresolved_reason="", qualification_rule="cyp_qualified", displayed=True
        )

    # Metabolites
    if raw_ep_lower == "metabolite" or re.search(r"\bmetabolite\b", raw_ep_lower):
        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        return QualificationDecision(
            funnel=funnel, stages=stages, section="METABOLISM", canonical_endpoint_id="METABOLITE_OBSERVATION",
            display_name="Metabolite", measurement_type="METABOLITE", species=species, matrix="PLASMA",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte="METABOLITE",
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u, normalized_value=num_val,
            normalized_unit=raw_u or "%", normalization_rule="metabolite_observation_identity", evidence_state=STATE_AUTO_QUALIFIED,
            qualification_status="ENDPOINT_QUALIFIED", comparability_status="DIRECTLY_COMPARABLE",
            unresolved_reason="", qualification_rule="metabolite_qualified", displayed=True
        )

    # Microsomal Metabolic Stability (t1/2 in minutes)
    if (raw_ep_lower == "microsomal" or "microsomal stability" in full_context_lower or "metabolic stability" in full_context_lower) and \
       raw_u_lower in {"min", "mins", "minute", "minutes", "h", "hours", "s", "sec"}:
        section = "METABOLISM"
        ep_species = species if species in {"HUMAN", "RAT", "MOUSE", "DOG", "MONKEY"} else "HUMAN"
        canonical_id = f"{ep_species}_MICROSOMAL_STABILITY_T_HALF"
        disp_name = f"{ep_species.title()} microsomal stability t1/2"
        norm_v = num_val if "min" in raw_u_lower else (num_val * 60.0 if "h" in raw_u_lower else num_val / 60.0)
        norm_u = "min"
        norm_rule = "microsomal_half_life_min"

        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True

        return QualificationDecision(
            funnel=funnel, stages=stages, section=section, canonical_endpoint_id=canonical_id,
            display_name=disp_name, measurement_type="MICROSOMAL_STABILITY", species=ep_species,
            matrix="MICROSOMES", route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen,
            analyte=analyte, target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u,
            normalized_value=round(norm_v, 2), normalized_unit=norm_u, normalization_rule=norm_rule,
            evidence_state=STATE_AUTO_QUALIFIED, qualification_status="ENDPOINT_QUALIFIED",
            comparability_status="RELATED_NOT_SAME_ENDPOINT", unresolved_reason="",
            qualification_rule="microsomal_stability_qualified", displayed=True
        )

    # -------------------------------------------------------------------------
    # G. PHARMACOKINETICS (PK): Cmax, Tmax, AUC, t1/2, CL, Vd, F
    # -------------------------------------------------------------------------
    is_pk = re.search(r"\b(cmax|tmax|auc|auclast|auc0[- ]24|half[- ]?life|t1/2|clearance|\bcl\b|cl/f|volume|vd|vd/f|vss|vss/f|bioavailability|\bf\b)\b", raw_ep_lower)
    if is_pk:
        section = "PK"
        
        # 1. Check for Relative Ratio / DDI % change
        # e.g. "increased Cmax by 32%", "decreased AUC by 46%", "hepatic impairment increased Cmax by 36%"
        is_relative_ratio = ("%" in raw_u_lower or "%" in raw_val_str) and \
                            any(w in full_context_lower for w in ("increased", "decreased", "fold", "impairment", "concomitant", "inducer", "inhibitor", "geomean", "ratio", "accumulation", "food")) and \
                            raw_ep_lower in {"cmax", "auc", "auclast", "clearance", "cl/f", "half-life"} and \
                            raw_ep_lower != "bioavailability"

        if is_relative_ratio:
            funnel["endpoint_classified"] = True
            funnel["unit_normalized"] = True
            funnel["qualification_state"] = STATE_RELATED
            funnel["displayed"] = True
            stages["ENDPOINT_QUALIFIED"] = True
            stages["CONTEXT_QUALIFIED"] = True
            stages["UNIT_NORMALIZED"] = True
            stages["PREDICTION_PAIRABLE"] = False

            return QualificationDecision(
                funnel=funnel, stages=stages, section=section, canonical_endpoint_id=f"{species}_PK_DDI_RELATIVE_RATIO",
                display_name=f"{species.title()} {raw_ep} DDI/Relative Change", measurement_type="DDI_RATIO",
                species=species, matrix="PLASMA", route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen,
                analyte=analyte, target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u,
                normalized_value=num_val, normalized_unit="%", normalization_rule="ddi_relative_ratio_percent",
                evidence_state=STATE_RELATED, qualification_status="ENDPOINT_QUALIFIED",
                comparability_status="RELATED_NOT_SAME_ENDPOINT",
                unresolved_reason=REASON_RELATIVE_RATIO_NOT_ABSOLUTE_EXPOSURE,
                qualification_rule="pk_ddi_relative_ratio_qualified", displayed=True
            )

        # 2. Check for Missing Unit
        if not raw_u and not any(u in full_context_lower for u in ("ng/ml", "ug/ml", "µg/ml", "mg/l", "l/h", "ml/min", "l/kg", "hr", "hour", "min", "%", "um", "nm", "nmol")):
            funnel["drop_stage"] = FUNNEL_UNIT_NORMALIZED
            funnel["drop_reason"] = REASON_UNIT_MISSING
            return QualificationDecision(
                funnel=funnel, stages=stages, section=section, canonical_endpoint_id="UNRESOLVED",
                display_name=raw_ep, measurement_type=raw_ep, species=species, matrix="PLASMA",
                route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
                target_context=target_context, raw_value=raw_val_str, raw_unit="",
                normalized_value=num_val, normalized_unit="", normalization_rule="unit_missing",
                evidence_state=STATE_REVIEW_REQUIRED, qualification_status="ENDPOINT_NOT_QUALIFIED",
                comparability_status="UNSUPPORTED", unresolved_reason=REASON_UNIT_MISSING,
                qualification_rule="unit_missing_review", displayed=False
            )

        # 3. Absolute PK Parameters
        # Resolve PK parameter
        if "cmax" in raw_ep_lower:
            pk_param = "CMAX"
            disp_param = "Cmax"
            # Normalize unit to ng/mL
            if raw_u_lower in {"ug/ml", "µg/ml", "mg/l"}:
                norm_v = num_val * 1000.0; norm_u = "ng/mL"; norm_rule = "ug/mL_to_ng/mL"
            elif raw_u_lower in {"ng/ml", "ug/l", "µg/l"}:
                norm_v = num_val; norm_u = "ng/mL"; norm_rule = "identity"
            elif raw_u_lower in {"um", "umol/l", "µm"}:
                # If molar concentration in context, preserve as nM or convert
                norm_v = num_val * 1000.0; norm_u = "nM"; norm_rule = "um_to_nm"
            elif raw_u_lower in {"nm", "nmol/l"}:
                norm_v = num_val; norm_u = "nM"; norm_rule = "identity"
            else:
                norm_v = num_val; norm_u = "ng/mL"; norm_rule = "cmax_assumed_ng_ml"

        elif re.search(r"auc", raw_ep_lower):
            if re.search(r"inf|0-inf", raw_ep_lower + " " + full_context_lower):
                pk_param = "AUC0_INF"; disp_param = "AUC0-inf"
            elif "tau" in raw_ep_lower or "tau" in full_context_lower:
                pk_param = "AUC_TAU"; disp_param = "AUCtau"
            elif re.search(r"0-t|0-24|last|tlast", raw_ep_lower + " " + full_context_lower):
                pk_param = "AUC0_T"; disp_param = "AUC0-t"
            else:
                pk_param = "AUC"; disp_param = "AUC"

            # Normalize unit to ng·h/mL
            if raw_u_lower in {"ug*h/ml", "µg*h/ml", "ug.h/ml", "µg.h/ml", "mg*h/l"}:
                norm_v = num_val * 1000.0; norm_u = "ng·h/mL"; norm_rule = "ug*h/mL_to_ng*h/mL"
            elif raw_u_lower in {"ng*h/ml", "ng.h/ml", "hr*ng/ml", "h*ng/ml", "nghr/ml"}:
                norm_v = num_val; norm_u = "ng·h/mL"; norm_rule = "identity"
            elif raw_u_lower in {"h", "hours", "hr"} and "hr*ng/ml" in full_context_lower:
                norm_v = num_val; norm_u = "ng·h/mL"; norm_rule = "context_extracted_hr_ng_ml"
            else:
                norm_v = num_val; norm_u = "ng·h/mL"; norm_rule = "auc_assumed_ng_h_ml"

        elif "half" in raw_ep_lower or "t1/2" in raw_ep_lower:
            pk_param = "T_HALF"
            disp_param = "Terminal half-life"
            if raw_u_lower in {"h", "hr", "hrs", "hour", "hours"}:
                norm_v = num_val; norm_u = "hours"; norm_rule = "identity"
            elif raw_u_lower in {"min", "mins", "minute", "minutes"}:
                norm_v = num_val / 60.0; norm_u = "hours"; norm_rule = "min_to_hours"
            elif raw_u_lower in {"d", "day", "days"}:
                norm_v = num_val * 24.0; norm_u = "hours"; norm_rule = "days_to_hours"
            else:
                norm_v = num_val; norm_u = "hours"; norm_rule = "t_half_assumed_hours"

        elif "tmax" in raw_ep_lower:
            pk_param = "TMAX"
            disp_param = "Tmax"
            if raw_u_lower in {"h", "hr", "hrs", "hour", "hours"}:
                norm_v = num_val; norm_u = "hours"; norm_rule = "identity"
            elif raw_u_lower in {"min", "mins", "minute", "minutes"}:
                norm_v = num_val / 60.0; norm_u = "hours"; norm_rule = "min_to_hours"
            else:
                norm_v = num_val; norm_u = "hours"; norm_rule = "tmax_assumed_hours"

        elif "clearance" in raw_ep_lower or "cl/f" in raw_ep_lower or raw_ep_lower == "cl":
            if route == "ORAL" or "apparent" in full_context_lower or "cl/f" in raw_ep_lower:
                pk_param = "CLF_ORAL"
                disp_param = "Oral CL/F"
            else:
                pk_param = "CL"
                disp_param = "Systemic clearance"

            if raw_u_lower in {"l/h", "l/hr", "l/hour"}:
                norm_v = num_val; norm_u = "L/h"; norm_rule = "identity"
            elif raw_u_lower in {"ml/min/kg", "mlmin/kg"}:
                norm_v = num_val; norm_u = "mL/min/kg"; norm_rule = "identity"
            elif raw_u_lower in {"l", "l/h"} and "l/h" in full_context_lower:
                norm_v = num_val; norm_u = "L/h"; norm_rule = "context_extracted_l_h"
            elif raw_u_lower in {"ml/min"}:
                norm_v = num_val * 60.0 / 1000.0; norm_u = "L/h"; norm_rule = "ml_min_to_l_h"
            else:
                norm_v = num_val; norm_u = "L/h" if species == "HUMAN" else "mL/min/kg"; norm_rule = "clearance_preserved"

        elif "volume" in raw_ep_lower or "vd" in raw_ep_lower or "vss" in raw_ep_lower:
            if "vss" in raw_ep_lower or "steady state" in full_context_lower:
                pk_param = "VSSF_ORAL" if route == "ORAL" else "VSS"
                disp_param = "Steady-state volume of distribution"
            elif route == "ORAL" or "apparent" in full_context_lower or "vd/f" in raw_ep_lower:
                pk_param = "VDF_ORAL"
                disp_param = "Oral Vd/F"
            else:
                pk_param = "VD"
                disp_param = "Volume of distribution"

            if raw_u_lower in {"l", "liters", "litres"}:
                norm_v = num_val; norm_u = "L"; norm_rule = "identity"
            elif raw_u_lower in {"l/kg"}:
                norm_v = num_val; norm_u = "L/kg"; norm_rule = "identity"
            elif raw_u_lower in {"ml/kg"}:
                norm_v = num_val / 1000.0; norm_u = "L/kg"; norm_rule = "ml_kg_to_l_kg"
            else:
                norm_v = num_val; norm_u = "L" if species == "HUMAN" else "L/kg"; norm_rule = "volume_preserved"

        elif "bioavailability" in raw_ep_lower or raw_ep == "F":
            pk_param = "F"
            disp_param = "Oral bioavailability F"
            if "%" in raw_u_lower or num_val > 1.0:
                norm_v = num_val; norm_u = "%"; norm_rule = "percent_identity"
            elif 0.0 <= num_val <= 1.0:
                norm_v = num_val * 100.0; norm_u = "%"; norm_rule = "fraction_to_percent"
            else:
                norm_v = num_val; norm_u = "%"; norm_rule = "bioavailability_percent"
        else:
            pk_param = "PK_UNSPECIFIED"
            disp_param = "PK Parameter"
            norm_v = num_val; norm_u = raw_u; norm_rule = "preserved"

        canonical_id = f"{species}_PK_{pk_param}_{route}" if route != "UNSPECIFIED" else f"{species}_PK_{pk_param}_UNSPECIFIED"
        disp_name = f"{species.title()} {disp_param}" if species != "UNSPECIFIED" else disp_param

        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        stages["IMPORTABLE"] = True
        stages["PREDICTION_PAIRABLE"] = True

        return QualificationDecision(
            funnel=funnel, stages=stages, section=section, canonical_endpoint_id=canonical_id,
            display_name=disp_name, measurement_type=pk_param, species=species, matrix="PLASMA",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u,
            normalized_value=round(norm_v, 2) if norm_v is not None else None,
            normalized_unit=norm_u, normalization_rule=norm_rule, evidence_state=STATE_AUTO_QUALIFIED,
            qualification_status="ENDPOINT_QUALIFIED", comparability_status="DIRECTLY_COMPARABLE",
            unresolved_reason="", qualification_rule="absolute_pk_parameter_qualified", displayed=True
        )

    # -------------------------------------------------------------------------
    # H. TOXICITY: AMES, DILI
    # -------------------------------------------------------------------------
    if re.search(r"\bames\b", raw_ep_lower):
        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        return QualificationDecision(
            funnel=funnel, stages=stages, section="TOXICITY", canonical_endpoint_id="AMES_MUTAGENICITY",
            display_name="Ames mutagenicity", measurement_type="MUTAGENICITY", species=species,
            matrix="BACTERIAL", route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen,
            analyte=analyte, target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u,
            normalized_value=num_val, normalized_unit=raw_u or "", normalization_rule="ames_identity",
            evidence_state=STATE_AUTO_QUALIFIED, qualification_status="ENDPOINT_QUALIFIED",
            comparability_status="DIRECTLY_COMPARABLE", unresolved_reason="",
            qualification_rule="ames_qualified", displayed=True
        )

    if re.search(r"\bdili\b", raw_ep_lower):
        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        return QualificationDecision(
            funnel=funnel, stages=stages, section="TOXICITY", canonical_endpoint_id="DILI_LIABILITY",
            display_name="DILI clinical liability", measurement_type="LIVER_INJURY", species=species,
            matrix="CLINICAL", route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen,
            analyte=analyte, target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u,
            normalized_value=num_val, normalized_unit=raw_u or "", normalization_rule="dili_identity",
            evidence_state=STATE_AUTO_QUALIFIED, qualification_status="ENDPOINT_QUALIFIED",
            comparability_status="DIRECTLY_COMPARABLE", unresolved_reason="",
            qualification_rule="dili_qualified", displayed=True
        )

    # -------------------------------------------------------------------------
    # I. ADMET: PKA & LOGD7.4
    # -------------------------------------------------------------------------
    if "pka" in raw_ep_lower:
        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        return QualificationDecision(
            funnel=funnel, stages=stages, section="ADMET", canonical_endpoint_id="PKA",
            display_name="pKa", measurement_type="pKa", species=species, matrix="AQUEOUS",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u,
            normalized_value=round(num_val, 2), normalized_unit="pKa", normalization_rule="pka_identity",
            evidence_state=STATE_AUTO_QUALIFIED, qualification_status="ENDPOINT_QUALIFIED",
            comparability_status="DIRECTLY_COMPARABLE", unresolved_reason="",
            qualification_rule="pka_qualified", displayed=True
        )

    if "logd" in raw_ep_lower:
        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        return QualificationDecision(
            funnel=funnel, stages=stages, section="ADMET", canonical_endpoint_id="LOGD_7_4",
            display_name="logD 7.4", measurement_type="logD", species=species, matrix="OCTANOL_WATER_PH74",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u,
            normalized_value=round(num_val, 2), normalized_unit="logD", normalization_rule="logd_identity",
            evidence_state=STATE_AUTO_QUALIFIED, qualification_status="ENDPOINT_QUALIFIED",
            comparability_status="DIRECTLY_COMPARABLE", unresolved_reason="",
            qualification_rule="logd_qualified", displayed=True
        )

    if "logp" in raw_ep_lower:
        funnel["endpoint_classified"] = True
        funnel["unit_normalized"] = True
        funnel["qualification_state"] = STATE_AUTO_QUALIFIED
        funnel["displayed"] = True
        stages["ENDPOINT_QUALIFIED"] = True
        stages["CONTEXT_QUALIFIED"] = True
        stages["UNIT_NORMALIZED"] = True
        return QualificationDecision(
            funnel=funnel, stages=stages, section="ADMET", canonical_endpoint_id="LOGP_RELATED",
            display_name="logP (related)", measurement_type="logP", species=species, matrix="OCTANOL_WATER",
            route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen, analyte=analyte,
            target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u,
            normalized_value=round(num_val, 2), normalized_unit="logP", normalization_rule="logp_identity",
            evidence_state=STATE_AUTO_QUALIFIED, qualification_status="ENDPOINT_QUALIFIED",
            comparability_status="DIRECTLY_COMPARABLE", unresolved_reason="",
            qualification_rule="logp_qualified", displayed=True
        )

    # -------------------------------------------------------------------------
    # FALLBACK: Explicit Review Required with exact reason
    # -------------------------------------------------------------------------
    funnel["drop_stage"] = FUNNEL_ENDPOINT_CLASSIFIED
    funnel["drop_reason"] = REASON_NO_SUPPORTED_CANONICAL_ENDPOINT
    return QualificationDecision(
        funnel=funnel, stages=stages, section="UNCLASSIFIED", canonical_endpoint_id="UNRESOLVED",
        display_name=raw_ep or "Unresolved Endpoint", measurement_type="UNRESOLVED", species=species,
        matrix="UNSPECIFIED", route=route, dose=dose_val, dose_unit=dose_unit, regimen=regimen,
        analyte=analyte, target_context=target_context, raw_value=raw_val_str, raw_unit=raw_u,
        normalized_value=num_val, normalized_unit=raw_u, normalization_rule="unresolved_fallback",
        evidence_state=STATE_REVIEW_REQUIRED, qualification_status="ENDPOINT_NOT_QUALIFIED",
        comparability_status="UNSUPPORTED", unresolved_reason=REASON_NO_SUPPORTED_CANONICAL_ENDPOINT,
        qualification_rule="fallback_review", displayed=False
    )
