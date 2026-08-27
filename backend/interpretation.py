"""Centralized Scientific Interpretation Registry.

Provides deterministic, versioned, and reference-backed interpretation rules for:
- Physicochemical properties & Lipinski / Veber drug-likeness criteria
- ADME prediction endpoints (Solubility, Permeability, PPB/fu, Microsomal stability)
- CYP & Transporter liability classifiers
- Safety toxicology endpoints (hERG, DILI, Ames)
- In vivo Pharmacokinetic & Translational parameters (CL, Vd, t1/2, F)
"""

from __future__ import annotations

from typing import Any


INTERPRETATION_REGISTRY_VERSION = "1.0.0"


# Centralized interpretation rules with directionality, reference thresholds, and scientific citations
INTERPRETATION_RULES: dict[str, dict[str, Any]] = {
    # Physicochemical & Drug-Likeness
    "mw": {
        "endpoint": "Molecular Weight",
        "unit": "g/mol",
        "rule_system": "Lipinski Rule of 5",
        "reference_threshold": "≤ 500",
        "thresholds": {"favorable_max": 500.0, "borderline_max": 600.0},
        "evaluate": lambda v: (
            ("FAVORABLE", "favorable", "Optimal small molecule size (≤500)")
            if v <= 500.0
            else (("INTERMEDIATE", "intermediate", "Borderline size (500–600)") if v <= 600.0 else ("UNFAVORABLE", "liability", "High molecular weight (>600)"))
        ),
        "source": "Lipinski et al., Adv Drug Deliv Rev 2001",
    },
    "clogp": {
        "endpoint": "cLogP (Crippen)",
        "unit": "log10",
        "rule_system": "Lipinski Rule of 5",
        "reference_threshold": "≤ 5.0",
        "thresholds": {"optimal_min": 0.0, "optimal_max": 5.0, "borderline_max": 6.0},
        "evaluate": lambda v: (
            ("FAVORABLE", "favorable", "Optimal lipophilicity (0–5.0)")
            if 0.0 <= v <= 5.0
            else (("INTERMEDIATE", "intermediate", "Borderline lipophilicity") if -1.0 <= v <= 6.0 else ("UNFAVORABLE", "liability", "Extreme lipophilicity / hydrophobicity"))
        ),
        "source": "Lipinski et al., Adv Drug Deliv Rev 2001",
    },
    "tpsa": {
        "endpoint": "Topological Polar Surface Area",
        "unit": "Å²",
        "rule_system": "Veber Rule",
        "reference_threshold": "≤ 140 Å²",
        "thresholds": {"favorable_max": 140.0, "cns_optimal_max": 90.0},
        "evaluate": lambda v: (
            ("FAVORABLE", "favorable", "Optimal polar surface area for oral permeability (≤140 Å²)")
            if v <= 140.0
            else ("UNFAVORABLE", "liability", "High polar surface area (>140 Å²); restricted oral absorption")
        ),
        "source": "Veber et al., J Med Chem 2002",
    },
    "hbd": {
        "endpoint": "Hydrogen Bond Donors",
        "unit": "count",
        "rule_system": "Lipinski Rule of 5",
        "reference_threshold": "≤ 5",
        "evaluate": lambda v: (
            ("FAVORABLE", "favorable", "Optimal H-bond donors (≤5)")
            if v <= 5
            else ("UNFAVORABLE", "liability", "Excessive H-bond donors (>5); desolvation penalty")
        ),
        "source": "Lipinski et al., Adv Drug Deliv Rev 2001",
    },
    "hba": {
        "endpoint": "Hydrogen Bond Acceptors",
        "unit": "count",
        "rule_system": "Lipinski Rule of 5",
        "reference_threshold": "≤ 10",
        "evaluate": lambda v: (
            ("FAVORABLE", "favorable", "Optimal H-bond acceptors (≤10)")
            if v <= 10
            else ("UNFAVORABLE", "liability", "Excessive H-bond acceptors (>10)")
        ),
        "source": "Lipinski et al., Adv Drug Deliv Rev 2001",
    },
    "rotb": {
        "endpoint": "Rotatable Bonds",
        "unit": "count",
        "rule_system": "Veber Rule",
        "reference_threshold": "≤ 10",
        "evaluate": lambda v: (
            ("FAVORABLE", "favorable", "Optimal molecular flexibility (≤10)")
            if v <= 10
            else ("UNFAVORABLE", "liability", "High conformational flexibility (>10); entropic penalty")
        ),
        "source": "Veber et al., J Med Chem 2002",
    },
    "fsp3": {
        "endpoint": "Fraction Csp3",
        "unit": "ratio",
        "rule_system": "Lovering Complexity",
        "reference_threshold": "≥ 0.42",
        "evaluate": lambda v: (
            ("FAVORABLE", "favorable", "Good 3D saturation and complexity (≥0.42)")
            if v >= 0.42
            else ("INTERMEDIATE", "intermediate", "Flat/aromatic dominant structure (<0.42)")
        ),
        "source": "Lovering et al., J Med Chem 2009",
    },
    "qed": {
        "endpoint": "Quantitative Estimate of Drug-likeness",
        "unit": "score (0–1)",
        "rule_system": "Bickerton QED",
        "reference_threshold": "≥ 0.67 (attractive)",
        "evaluate": lambda v: (
            ("FAVORABLE", "favorable", "High drug-likeness (≥0.67)")
            if v >= 0.67
            else (("INTERMEDIATE", "intermediate", "Moderate drug-likeness (0.49–0.67)") if v >= 0.49 else ("UNFAVORABLE", "liability", "Low drug-likeness (<0.49)"))
        ),
        "source": "Bickerton et al., Nature Chem 2012",
    },

    # ADMET Prediction Endpoints
    "solubility": {
        "endpoint": "Aqueous Solubility",
        "unit": "µM",
        "thresholds": {"high": 60.0, "moderate": 10.0},
        "evaluate": lambda v: (
            ("HIGH SOLUBILITY", "favorable", "High aqueous solubility (>60 µM); low dissolution risk")
            if v >= 60.0
            else (("MODERATE SOLUBILITY", "intermediate", "Moderate solubility (10–60 µM)") if v >= 10.0 else ("LOW SOLUBILITY", "liability", "Low solubility (<10 µM); absorption/formulation risk"))
        ),
        "source": "FDA Biopharmaceutics Classification System (BCS)",
    },
    "caco2": {
        "endpoint": "Caco-2 Permeability",
        "unit": "log10 cm/s",
        "thresholds": {"high_log": -5.0, "moderate_log": -6.0},
        "evaluate": lambda v: (
            ("HIGH PERMEABILITY", "favorable", "High apparent permeability (> -5.0 log cm/s)")
            if v >= -5.0
            else (("MODERATE PERMEABILITY", "intermediate", "Moderate permeability (-6.0 to -5.0)") if v >= -6.0 else ("LOW PERMEABILITY", "liability", "Poor membrane permeability (< -6.0)"))
        ),
        "source": "Artursson & Karlsson, Biochem Biophys Res Commun 1991",
    },
    "ppb": {
        "endpoint": "Plasma Protein Binding (fu)",
        "unit": "fraction unbound (fu)",
        "thresholds": {"low_binding_fu": 0.20, "moderate_binding_fu": 0.05},
        "evaluate": lambda fu: (
            ("HIGH UNBOUND FRACTION", "favorable", f"fu = {fu:.3f} (PPB {(1-fu)*100:.1f}%); ample free drug")
            if fu >= 0.20
            else (("MODERATE BINDING", "intermediate", f"fu = {fu:.3f} (PPB {(1-fu)*100:.1f}%)") if fu >= 0.05 else ("HIGHLY BOUND", "intermediate", f"fu = {fu:.3f} (PPB {(1-fu)*100:.1f}%); extensive plasma binding"))
        ),
        "source": "Smith et al., J Med Chem 2010",
    },
    "hlm_clint": {
        "endpoint": "Human Liver Microsomal Clearance",
        "unit": "mL/min/kg",
        "thresholds": {"low": 15.0, "high": 45.0},
        "evaluate": lambda v: (
            ("LOW CLEARANCE", "favorable", "Low metabolic turnover (<15 mL/min/kg); high metabolic stability")
            if v <= 15.0
            else (("MODERATE CLEARANCE", "intermediate", "Moderate metabolic turnover (15–45 mL/min/kg)") if v <= 45.0 else ("HIGH CLEARANCE", "liability", "Rapid hepatic microsomal turnover (>45 mL/min/kg); metabolic liability"))
        ),
        "source": "Obach et al., J Pharmacol Exp Ther 1997",
    },
    "rlm_clint": {
        "endpoint": "Rat Liver Microsomal Clearance",
        "unit": "mL/min/kg",
        "thresholds": {"low": 20.0, "high": 60.0},
        "evaluate": lambda v: (
            ("LOW CLEARANCE", "favorable", "Low rat microsomal turnover (<20 mL/min/kg)")
            if v <= 20.0
            else (("MODERATE CLEARANCE", "intermediate", "Moderate rat clearance (20–60 mL/min/kg)") if v <= 60.0 else ("HIGH CLEARANCE", "liability", "High rat microsomal clearance (>60 mL/min/kg)"))
        ),
        "source": "Obach et al., J Pharmacol Exp Ther 1997",
    },
    "mlm_clint": {
        "endpoint": "Mouse Liver Microsomal Clearance",
        "unit": "mL/min/kg",
        "thresholds": {"low": 30.0, "high": 90.0},
        "evaluate": lambda v: (
            ("LOW CLEARANCE", "favorable", "Low mouse microsomal turnover (<30 mL/min/kg)")
            if v <= 30.0
            else (("MODERATE CLEARANCE", "intermediate", "Moderate mouse clearance (30–90 mL/min/kg)") if v <= 90.0 else ("HIGH CLEARANCE", "liability", "High mouse microsomal clearance (>90 mL/min/kg)"))
        ),
        "source": "Obach et al., J Pharmacol Exp Ther 1997",
    },

    # Safety & Liability Classifiers (Positive = Liability / RED, Negative = Favorable / BLUE)
    "herg": {
        "endpoint": "hERG Liability",
        "unit": "probability (0–1)",
        "evaluate": lambda prob: (
            ("HIGH LIABILITY", "liability", f"Positive inhibitor liability ({prob:.2f} prob); cardiotoxicity risk")
            if prob >= 0.50
            else ("LOW LIABILITY", "favorable", f"Negative liability ({1.0-prob:.2f} safety prob); low cardiotoxicity risk")
        ),
        "source": "Redfern et al., Cardiovasc Res 2003",
    },
    "dili": {
        "endpoint": "DILI Clinical Liability",
        "unit": "probability (0–1)",
        "evaluate": lambda prob: (
            ("HIGH LIABILITY", "liability", f"Positive clinical DILI concern ({prob:.2f} prob); hepatotoxicity risk")
            if prob >= 0.50
            else ("LOW LIABILITY", "favorable", f"Negative liability ({1.0-prob:.2f} safety prob); low hepatotoxicity concern")
        ),
        "source": "Chen et al., Drug Discov Today 2011",
    },
    "ames": {
        "endpoint": "Ames Mutagenicity",
        "unit": "probability (0–1)",
        "evaluate": lambda prob: (
            ("HIGH LIABILITY", "liability", f"Positive Ames mutagenicity ({prob:.2f} prob); genotoxicity risk")
            if prob >= 0.50
            else ("LOW LIABILITY", "favorable", f"Negative Ames mutagenicity ({1.0-prob:.2f} safety prob); non-mutagenic")
        ),
        "source": "Hansen et al., J Chem Inf Model 2009",
    },
    "pgp_inhibitor": {
        "endpoint": "P-gp Inhibitor",
        "unit": "probability (0–1)",
        "evaluate": lambda prob: (
            ("POTENTIAL LIABILITY", "liability", f"Predicted P-gp inhibitor ({prob:.2f} prob); transporter DDI concern")
            if prob >= 0.50
            else ("LOW LIABILITY", "favorable", f"Non-inhibitor ({1.0-prob:.2f} prob); low transporter interaction")
        ),
        "source": "Giacomini et al., Nat Rev Drug Discov 2010",
    },
    "cyp_inhibitor": {
        "endpoint": "CYP Inhibitor",
        "unit": "probability (0–1)",
        "evaluate": lambda prob: (
            ("POTENTIAL LIABILITY", "liability", f"Predicted CYP inhibitor ({prob:.2f} prob); potential DDI liability")
            if prob >= 0.50
            else ("LOW LIABILITY", "favorable", f"Low inhibition concern ({1.0-prob:.2f} prob)")
        ),
        "source": "FDA Guidance: In Vitro Drug Interaction Studies 2020",
    },
    "cyp_substrate": {
        "endpoint": "CYP Substrate",
        "unit": "probability (0–1)",
        "evaluate": lambda prob: (
            ("POTENTIAL SUBSTRATE", "intermediate", f"Predicted substrate ({prob:.2f} prob); active metabolic pathway")
            if prob >= 0.50
            else ("NON-SUBSTRATE", "favorable", f"Low substrate likelihood ({1.0-prob:.2f} prob)")
        ),
        "source": "FDA Guidance: In Vitro Drug Interaction Studies 2020",
    },

    # In Vivo Pharmacokinetics
    "vd": {
        "endpoint": "Volume of Distribution (Vss)",
        "unit": "L/kg",
        "thresholds": {"low": 0.7, "high": 2.0},
        "evaluate": lambda v: (
            ("LOW DISTRIBUTION", "intermediate", "Confined primarily to plasma/extracellular water (<0.7 L/kg)")
            if v < 0.7
            else (("MODERATE DISTRIBUTION", "favorable", "Well-balanced tissue distribution (0.7–2.0 L/kg)") if v <= 2.0 else ("HIGH DISTRIBUTION", "favorable", "Extensive tissue binding / partitioning (>2.0 L/kg)"))
        ),
        "source": "Smith et al., J Med Chem 2015",
    },
    "bioavailability": {
        "endpoint": "Oral Bioavailability (F)",
        "unit": "%",
        "thresholds": {"low": 30.0, "high": 60.0},
        "evaluate": lambda f: (
            ("HIGH BIOAVAILABILITY", "favorable", "High oral bioavailability (>60%)")
            if f >= 60.0
            else (("MODERATE BIOAVAILABILITY", "intermediate", "Moderate oral absorption (30–60%)") if f >= 30.0 else ("LOW BIOAVAILABILITY", "liability", "Low oral bioavailability (<30%); formulation/first-pass limitation"))
        ),
        "source": "Veber et al., J Med Chem 2002",
    },
}


def interpret_property(key: str, value: Any) -> dict[str, Any]:
    """Retrieve deterministic scientific evaluation and color mapping for a property or endpoint."""
    rule = INTERPRETATION_RULES.get(key.lower())
    if not rule or value is None:
        return {
            "status": "UNAVAILABLE",
            "assessment": "UNAVAILABLE",
            "color_class": "unavailable",
            "interpretation": "No reference interpretation available for this property.",
            "rule_system": rule.get("rule_system") if rule else "None",
            "reference_threshold": rule.get("reference_threshold") if rule else "—",
        }

    try:
        val_float = float(value)
        assessment, color_class, interpretation = rule["evaluate"](val_float)
        return {
            "status": "EVALUATED",
            "assessment": assessment,
            "color_class": color_class,
            "interpretation": interpretation,
            "rule_system": rule.get("rule_system", "Scientific Threshold"),
            "reference_threshold": rule.get("reference_threshold", "—"),
            "source": rule.get("source", "Literature Guideline"),
        }
    except (ValueError, TypeError):
        return {
            "status": "UNPARSEABLE",
            "assessment": "UNPARSEABLE",
            "color_class": "unavailable",
            "interpretation": f"Cannot numerically evaluate value: {value}",
            "rule_system": rule.get("rule_system", "None"),
            "reference_threshold": rule.get("reference_threshold", "—"),
        }


def get_interpretation_registry_summary() -> dict[str, Any]:
    """Return public metadata of the interpretation registry."""
    return {
        "version": INTERPRETATION_REGISTRY_VERSION,
        "rules_count": len(INTERPRETATION_RULES),
        "rules": {
            k: {
                "endpoint": v["endpoint"],
                "unit": v.get("unit", ""),
                "rule_system": v.get("rule_system", "Threshold Rule"),
                "reference_threshold": v.get("reference_threshold", "—"),
                "source": v.get("source", ""),
            }
            for k, v in INTERPRETATION_RULES.items()
        },
    }
