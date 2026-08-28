"""
Drug-OPT Stage 4D-0: Authoritative Scientific Endpoint Contracts & Ensemble Governance.

This module defines deterministic, scientifically bounded endpoint contracts,
compatibility verification gates, unit conversion rules, and model adapter
interfaces for future multi-model ensemble architectures.

DESIGN PRINCIPLES:
1. Scientific Semantics First: Two models are NEVER ensemble-compatible merely
   because they share a keyword name (e.g. 'solubility', 'clearance', 'hERG').
2. Strict Unit & Transformation Consistency: Units must match exactly or be
   transformable via mathematically exact transformations.
3. Species & Role Isolation: Interspecies parameter reuse and role confusion
   (e.g. inhibitor vs substrate, in vivo clearance vs in vitro rate) are hard failures.
4. Failure Isolation: A failure in one ensemble member must never crash the
   endpoint; valid remaining predictions continue cleanly.
5. Deterministic Prediction Freeze: Prospective predictions must be frozen
   with exact model, version, input representation, and timestamp before
   experimental results arrive.
"""

from __future__ import annotations

import enum
import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Callable, Dict, List, Optional, Set, Tuple


class EndpointCategory(str, enum.Enum):
    PHYSICOCHEMICAL = "PHYSICOCHEMICAL"
    ABSORPTION = "ABSORPTION"
    DISTRIBUTION = "DISTRIBUTION"
    METABOLISM = "METABOLISM"
    EXCRETION = "EXCRETION"
    CYP_PANEL = "CYP_PANEL"
    TRANSPORTER = "TRANSPORTER"
    SAFETY = "SAFETY"
    PHARMACOKINETICS = "PHARMACOKINETICS"


class OutputType(str, enum.Enum):
    REGRESSION = "REGRESSION"
    BINARY_CLASSIFICATION = "BINARY_CLASSIFICATION"
    MULTI_CLASSIFICATION = "MULTI_CLASSIFICATION"
    RANKING = "RANKING"
    MECHANISTIC_DERIVED = "MECHANISTIC_DERIVED"


class Directionality(str, enum.Enum):
    HIGHER_BETTER = "HIGHER_BETTER"
    LOWER_BETTER = "LOWER_BETTER"
    NEUTRAL = "NEUTRAL"
    CLASSIFICATION_BINARY = "CLASSIFICATION_BINARY"


class QualificationStatus(str, enum.Enum):
    QUALIFIED = "QUALIFIED"
    QUALIFIED_WITH_LIMITATIONS = "QUALIFIED_WITH_LIMITATIONS"
    RESEARCH_ONLY = "RESEARCH_ONLY"
    REJECTED = "REJECTED"
    UNAVAILABLE = "UNAVAILABLE"


class ARM64Status(str, enum.Enum):
    RUNS_LOCAL_ARM64 = "RUNS_LOCAL_ARM64"
    RUNS_LOCAL_CPU_SLOW = "RUNS_LOCAL_CPU_SLOW"
    REQUIRES_ISOLATED_ENV = "REQUIRES_ISOLATED_ENV"
    GPU_OPTIONAL = "GPU_OPTIONAL"
    GPU_REQUIRED = "GPU_REQUIRED"
    EXTERNAL_SERVICE_REQUIRED = "EXTERNAL_SERVICE_REQUIRED"
    NOT_COMPATIBLE = "NOT_COMPATIBLE"


class ExecutionTier(str, enum.Enum):
    TIER_1_LOCAL_FAST = "TIER_1_LOCAL_FAST"      # Small, sub-second CPU models (Descriptor GBDT, lightweight FFN)
    TIER_2_LOCAL_HEAVY = "TIER_2_LOCAL_HEAVY"    # Slower CPU/GPU models (Deep MPNN, Chemprop Ensembles)
    TIER_3_OPTIONAL_EXTERNAL = "TIER_3_OPTIONAL_EXTERNAL"  # Heavy services / distributed workers


class LicenseClassification(str, enum.Enum):
    PERMISSIVE_OPEN_SOURCE = "PERMISSIVE_OPEN_SOURCE"        # MIT, Apache-2.0, BSD
    COPYLEFT_OPEN_SOURCE = "COPYLEFT_OPEN_SOURCE"            # LGPL, GPL-3.0
    ACADEMIC_NON_COMMERCIAL = "ACADEMIC_NON_COMMERCIAL"      # Non-commercial only
    PROPRIETARY = "PROPRIETARY"                              # Paid commercial license required
    LICENSE_REVIEW_REQUIRED = "LICENSE_REVIEW_REQUIRED"      # Ambiguous or upstream dataset restrictions


@dataclass(frozen=True)
class EndpointContract:
    """
    Authoritative scientific contract defining the physical and biological
    meaning of a Drug-OPT prediction endpoint.
    """
    endpoint_id: str
    display_name: str
    category: EndpointCategory
    scientific_definition: str
    species: str
    assay_type: str
    raw_unit: str
    canonical_unit: str
    transformation: str
    output_type: OutputType
    directionality: Directionality
    classification_semantics: Optional[Dict[str, Any]] = None
    cutoff_definition: Optional[Dict[str, Any]] = None
    required_metadata: Tuple[str, ...] = ()
    experimental_compatibility_rules: Dict[str, Any] = field(default_factory=dict)
    comparison_rules: Dict[str, Any] = field(default_factory=dict)
    ensemble_compatibility_rules: Dict[str, Any] = field(default_factory=dict)
    version: str = "1.0.0"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "endpoint_id": self.endpoint_id,
            "display_name": self.display_name,
            "category": self.category.value,
            "scientific_definition": self.scientific_definition,
            "species": self.species,
            "assay_type": self.assay_type,
            "raw_unit": self.raw_unit,
            "canonical_unit": self.canonical_unit,
            "transformation": self.transformation,
            "output_type": self.output_type.value,
            "directionality": self.directionality.value,
            "classification_semantics": self.classification_semantics,
            "cutoff_definition": self.cutoff_definition,
            "required_metadata": list(self.required_metadata),
            "experimental_compatibility_rules": self.experimental_compatibility_rules,
            "comparison_rules": self.comparison_rules,
            "ensemble_compatibility_rules": self.ensemble_compatibility_rules,
            "version": self.version,
        }


# ==============================================================================
# AUTHORITATIVE ENDPOINT CONTRACT REGISTRY
# ==============================================================================

ENDPOINT_CONTRACTS: Dict[str, EndpointContract] = {
    # 1. SOLUBILITY
    "Solubility": EndpointContract(
        endpoint_id="solubility_aqueous_logs",
        display_name="Aqueous Solubility (LogS)",
        category=EndpointCategory.PHYSICOCHEMICAL,
        scientific_definition="Aqueous solubility expressed as LogS = log10(S [mol/L]). Represents room-temperature aggregate aqueous solubility across neutral/buffered media, not a pH-dependent micro-state or intrinsic solubility S0.",
        species="Chemical / In Vitro",
        assay_type="Aggregate Aqueous Shake-Flask / Nephelometric / Potentiometric",
        raw_unit="log10(mol/L)",
        canonical_unit="log10(mol/L)",
        transformation="Identity (already in log10 mol/L)",
        output_type=OutputType.REGRESSION,
        directionality=Directionality.HIGHER_BETTER,
        cutoff_definition={
            "high_solubility_min": -2.0,       # > 10 mM (high)
            "moderate_solubility_min": -4.0,   # 100 uM - 10 mM (moderate)
            "low_solubility_max": -4.0,        # < 100 uM (low)
            "poor_solubility_max": -6.0,       # < 1 uM (very low / insoluble)
        },
        required_metadata=("temperature_celsius", "solvent_system"),
        experimental_compatibility_rules={
            "compatible_experimental_types": ["Thermodynamic aqueous solubility", "Kinetic aqueous solubility (low DMSO <1%)"],
            "incompatible_types": ["Organic solvent solubility", "Mass solubility without MW", "High DMSO cosolvent (>5%)"],
            "conversion_from_mass": "log10(S_mg_per_mL / (MW * 1e-3))",
        },
        comparison_rules={"max_compatible_ph_difference": 1.0},
        ensemble_compatibility_rules={
            "required_species": "Chemical / In Vitro",
            "required_canonical_unit": "log10(mol/L)",
            "forbidden_mix_types": ["intrinsic_s0_without_pka", "kinetic_dmso_turbidity_raw"],
        },
    ),

    # 2. CACO-2 PERMEABILITY
    "Permeability": EndpointContract(
        endpoint_id="permeability_caco2_logpapp",
        display_name="Caco-2 Apparent Permeability (LogPapp)",
        category=EndpointCategory.ABSORPTION,
        scientific_definition="Apparent membrane permeability measured in human colon carcinoma (Caco-2) monolayers, expressed as LogPapp = log10(Papp [cm/s]). Represents overall apical-to-basolateral (A->B) transcellular and paracellular flux.",
        species="Human Caco-2 Cell Line",
        assay_type="Transwell Monolayer Permeability (pH 7.4/7.4 or pH 6.5/7.4)",
        raw_unit="log10(cm/s)",
        canonical_unit="log10(cm/s)",
        transformation="Identity (log10 cm/s)",
        output_type=OutputType.REGRESSION,
        directionality=Directionality.HIGHER_BETTER,
        cutoff_definition={
            "high_permeability_min": -5.0,    # > 10 x 10^-6 cm/s (High)
            "moderate_permeability_min": -6.0, # 1 to 10 x 10^-6 cm/s (Medium)
            "low_permeability_max": -6.0,      # < 1 x 10^-6 cm/s (Low)
        },
        required_metadata=("direction", "ph_apical", "ph_basolateral"),
        experimental_compatibility_rules={
            "compatible_experimental_types": ["Caco-2 A->B apparent permeability", "Caco-2 Papp"],
            "incompatible_types": ["PAMPA", "MDCK-WT", "MDCK-MDR1", "RBE4 BBB permeability", "Efflux Ratio (B->A / A->B)"],
            "unit_conversions": {
                "10^-6 cm/s": "log10(val * 1e-6)",
                "nm/s": "log10(val * 1e-7)",
                "um/s": "log10(val * 1e-4)",
            },
        },
        comparison_rules={"direction_must_match": True, "strict_assay_isolation": True},
        ensemble_compatibility_rules={
            "required_species": "Human Caco-2 Cell Line",
            "required_canonical_unit": "log10(cm/s)",
            "forbidden_mix_types": ["pampa_artificial_membrane", "mdck_canine_cells", "efflux_ratio"],
        },
    ),

    # 3. PLASMA PROTEIN BINDING (PPB)
    "Plasma protein binding": EndpointContract(
        endpoint_id="ppb_human_percent_bound",
        display_name="Human Plasma Protein Binding (% Bound)",
        category=EndpointCategory.DISTRIBUTION,
        scientific_definition="Percent of compound bound to plasma proteins in pooled human plasma at equilibrium. Fraction unbound fu is derived strictly as fu = (100 - % bound) / 100.",
        species="Human",
        assay_type="Equilibrium Dialysis / Rapid Equilibrium Dialysis (RED) / Ultrafiltration",
        raw_unit="% bound",
        canonical_unit="% bound",
        transformation="Identity (% bound, range [0.0, 100.0])",
        output_type=OutputType.REGRESSION,
        directionality=Directionality.NEUTRAL,
        cutoff_definition={
            "low_binding_max": 80.0,          # fu > 0.20 (Low binding / High fu)
            "moderate_binding_min": 80.0,
            "moderate_binding_max": 95.0,      # fu = 0.05 - 0.20 (Moderate)
            "high_binding_min": 95.0,
            "high_binding_max": 99.0,          # fu = 0.01 - 0.05 (High)
            "very_high_binding_min": 99.0,     # fu < 0.01 (Very high binding / Extreme fu restriction)
        },
        required_metadata=("matrix", "species", "protein_concentration_g_per_L"),
        experimental_compatibility_rules={
            "compatible_experimental_types": ["Human plasma equilibrium dialysis % bound", "Human plasma RED % bound", "Human plasma ultrafiltration % bound"],
            "incompatible_types": ["Bovine Serum Albumin (BSA) binding", "Alpha-1-acid glycoprotein (AAG) isolated binding", "Rat/Mouse plasma binding"],
            "unit_conversions": {
                "fraction_bound": "val * 100.0",
                "fraction_unbound_fu": "(1.0 - val) * 100.0",
                "%_unbound": "100.0 - val",
            },
        },
        comparison_rules={"species_must_match": "Human", "matrix_must_match": "Plasma"},
        ensemble_compatibility_rules={
            "required_species": "Human",
            "required_canonical_unit": "% bound",
            "derived_fu_formula": "fu = max(0.0001, (100.0 - percent_bound) / 100.0)",
            "forbidden_mix_types": ["animal_plasma_ppb", "isolated_albumin_binding"],
        },
    ),

    # 4. HLM INTRINSIC CLEARANCE
    "HLM intrinsic clearance": EndpointContract(
        endpoint_id="hlm_intrinsic_clearance_scaled_log10",
        display_name="Human Liver Microsomal Clearance (HLM Clint)",
        category=EndpointCategory.METABOLISM,
        scientific_definition="Human liver microsomal intrinsic clearance scaled to in vivo body weight clearance and expressed in log10(mL/min/kg). Represents Phase I oxidative/hydrolytic metabolic rate under pooled human hepatic microsomes.",
        species="Human",
        assay_type="Pooled Human Liver Microsomes (HLM) Substrate Depletion + NADPH",
        raw_unit="log10(mL/min/kg)",
        canonical_unit="log10(mL/min/kg)",
        transformation="Identity (log10 mL/min/kg)",
        output_type=OutputType.REGRESSION,
        directionality=Directionality.LOWER_BETTER,
        cutoff_definition={
            "stable_max": 0.903090,      # <= 8.0 mL/min/kg (Stable / Low clearance)
            "moderate_max": 1.741998,    # 8.0 to 55.2 mL/min/kg (Moderate clearance)
            "unstable_min": 1.741998,    # >= 55.2 mL/min/kg (Unstable / High clearance)
        },
        required_metadata=("species", "matrix", "mppgl_mg_per_g", "liver_weight_g_per_kg"),
        experimental_compatibility_rules={
            "compatible_experimental_types": ["HLM in vitro substrate depletion Clint (scaled)"],
            "incompatible_types": ["Cryopreserved Human Hepatocyte Clint (without hepatocyte scaling)", "Rat/Mouse liver microsomes", "Recombinant CYP isoform Clint"],
            "scaling_from_raw_in_vitro": "Clint_scaled_mL_min_kg = Clint_uL_min_mg * MPPGL(45.0) * LiverWeight(25.7) / 1000.0",
        },
        comparison_rules={"species_must_match": "Human", "matrix_must_match": "Microsomes"},
        ensemble_compatibility_rules={
            "required_species": "Human",
            "required_canonical_unit": "log10(mL/min/kg)",
            "forbidden_mix_types": ["raw_uL_min_mg_unscaled", "microsomal_half_life_min", "hepatocyte_clearance"],
        },
    ),

    # 5. RLM INTRINSIC CLEARANCE
    "RLM intrinsic clearance": EndpointContract(
        endpoint_id="rlm_intrinsic_clearance_scaled_log10",
        display_name="Rat Liver Microsomal Clearance (RLM Clint)",
        category=EndpointCategory.METABOLISM,
        scientific_definition="Rat liver microsomal intrinsic clearance scaled to in vivo body weight clearance in log10(mL/min/kg). Represents Phase I metabolic rate in pooled Sprague-Dawley rat liver microsomes.",
        species="Rat",
        assay_type="Pooled Rat Liver Microsomes (RLM) Substrate Depletion + NADPH",
        raw_unit="log10(mL/min/kg)",
        canonical_unit="log10(mL/min/kg)",
        transformation="Identity (log10 mL/min/kg)",
        output_type=OutputType.REGRESSION,
        directionality=Directionality.LOWER_BETTER,
        cutoff_definition={
            "stable_max": 1.301030,      # <= 20.0 mL/min/kg
            "moderate_max": 2.171360,    # 20.0 to 148.4 mL/min/kg
            "unstable_min": 2.171360,    # >= 148.4 mL/min/kg
        },
        required_metadata=("species", "matrix", "mppgl_mg_per_g", "liver_weight_g_per_kg"),
        experimental_compatibility_rules={
            "compatible_experimental_types": ["RLM in vitro substrate depletion Clint (scaled)"],
            "incompatible_types": ["Human/Mouse microsomes", "In vivo rat PK clearance without IVIVE"],
            "scaling_from_raw_in_vitro": "Clint_scaled_mL_min_kg = Clint_uL_min_mg * MPPGL(45.0) * LiverWeight(40.0) / 1000.0",
        },
        comparison_rules={"species_must_match": "Rat", "matrix_must_match": "Microsomes"},
        ensemble_compatibility_rules={
            "required_species": "Rat",
            "required_canonical_unit": "log10(mL/min/kg)",
            "forbidden_mix_types": ["human_hlm", "mouse_mlm", "unscaled_half_life"],
        },
    ),

    # 6. MLM INTRINSIC CLEARANCE
    "MLM intrinsic clearance": EndpointContract(
        endpoint_id="mlm_intrinsic_clearance_scaled_log10",
        display_name="Mouse Liver Microsomal Clearance (MLM Clint)",
        category=EndpointCategory.METABOLISM,
        scientific_definition="Mouse liver microsomal intrinsic clearance scaled to in vivo body weight clearance in log10(mL/min/kg). Represents Phase I metabolic rate in pooled CD-1 mouse liver microsomes.",
        species="Mouse",
        assay_type="Pooled Mouse Liver Microsomes (MLM) Substrate Depletion + NADPH",
        raw_unit="log10(mL/min/kg)",
        canonical_unit="log10(mL/min/kg)",
        transformation="Identity (log10 mL/min/kg)",
        output_type=OutputType.REGRESSION,
        directionality=Directionality.LOWER_BETTER,
        cutoff_definition={
            "stable_max": 1.930057,      # <= 85.1 mL/min/kg
            "moderate_max": 2.818209,    # 85.1 to 658.0 mL/min/kg
            "unstable_min": 2.818209,    # >= 658.0 mL/min/kg
        },
        required_metadata=("species", "matrix", "mppgl_mg_per_g", "liver_weight_g_per_kg"),
        experimental_compatibility_rules={
            "compatible_experimental_types": ["MLM in vitro substrate depletion Clint (scaled)"],
            "incompatible_types": ["Human/Rat microsomes"],
            "scaling_from_raw_in_vitro": "Clint_scaled_mL_min_kg = Clint_uL_min_mg * MPPGL(45.0) * LiverWeight(87.5) / 1000.0",
        },
        comparison_rules={"species_must_match": "Mouse", "matrix_must_match": "Microsomes"},
        ensemble_compatibility_rules={
            "required_species": "Mouse",
            "required_canonical_unit": "log10(mL/min/kg)",
            "forbidden_mix_types": ["human_hlm", "rat_rlm"],
        },
    ),

    # 7-11. CYP INHIBITION PANEL (1A2, 2C9, 2C19, 2D6, 3A4)
    "CYP1A2 inhibitor": EndpointContract(
        endpoint_id="cyp1a2_inhibitor_prob",
        display_name="CYP1A2 Inhibition Probability",
        category=EndpointCategory.CYP_PANEL,
        scientific_definition="Calibrated binary probability of inhibiting human CYP1A2 at AC50 <= 10 uM in recombinant fluorogenic/luminescent dealkylation assay (PubChem AID 1851 protocol).",
        species="Human",
        assay_type="Recombinant Human CYP1A2 Functional Dealkylation Assay",
        raw_unit="probability",
        canonical_unit="probability",
        transformation="Identity (probability in [0.0, 1.0])",
        output_type=OutputType.BINARY_CLASSIFICATION,
        directionality=Directionality.LOWER_BETTER,
        classification_semantics={
            "positive_class": "INHIBITOR",
            "negative_class": "NON_INHIBITOR",
            "decision_threshold": 0.50,
            "potency_threshold_um": 10.0,
        },
        required_metadata=("isoform", "role"),
        comparison_rules={"isoform_must_match": "CYP1A2", "role_must_match": "INHIBITOR"},
        ensemble_compatibility_rules={"required_isoform": "CYP1A2", "required_role": "INHIBITOR"},
    ),
    "CYP2C9 inhibitor": EndpointContract(
        endpoint_id="cyp2c9_inhibitor_prob",
        display_name="CYP2C9 Inhibition Probability",
        category=EndpointCategory.CYP_PANEL,
        scientific_definition="Calibrated binary probability of inhibiting human CYP2C9 at AC50 <= 10 uM in recombinant fluorogenic/luminescent assay.",
        species="Human",
        assay_type="Recombinant Human CYP2C9 Functional Assay",
        raw_unit="probability",
        canonical_unit="probability",
        transformation="Identity (probability in [0.0, 1.0])",
        output_type=OutputType.BINARY_CLASSIFICATION,
        directionality=Directionality.LOWER_BETTER,
        classification_semantics={
            "positive_class": "INHIBITOR",
            "negative_class": "NON_INHIBITOR",
            "decision_threshold": 0.50,
            "potency_threshold_um": 10.0,
        },
        required_metadata=("isoform", "role"),
        comparison_rules={"isoform_must_match": "CYP2C9", "role_must_match": "INHIBITOR"},
        ensemble_compatibility_rules={"required_isoform": "CYP2C9", "required_role": "INHIBITOR"},
    ),
    "CYP2C19 inhibitor": EndpointContract(
        endpoint_id="cyp2c19_inhibitor_prob",
        display_name="CYP2C19 Inhibition Probability",
        category=EndpointCategory.CYP_PANEL,
        scientific_definition="Calibrated binary probability of inhibiting human CYP2C19 at AC50 <= 10 uM in recombinant assay.",
        species="Human",
        assay_type="Recombinant Human CYP2C19 Functional Assay",
        raw_unit="probability",
        canonical_unit="probability",
        transformation="Identity (probability in [0.0, 1.0])",
        output_type=OutputType.BINARY_CLASSIFICATION,
        directionality=Directionality.LOWER_BETTER,
        classification_semantics={
            "positive_class": "INHIBITOR",
            "negative_class": "NON_INHIBITOR",
            "decision_threshold": 0.50,
            "potency_threshold_um": 10.0,
        },
        required_metadata=("isoform", "role"),
        comparison_rules={"isoform_must_match": "CYP2C19", "role_must_match": "INHIBITOR"},
        ensemble_compatibility_rules={"required_isoform": "CYP2C19", "required_role": "INHIBITOR"},
    ),
    "CYP2D6 inhibitor": EndpointContract(
        endpoint_id="cyp2d6_inhibitor_prob",
        display_name="CYP2D6 Inhibition Probability",
        category=EndpointCategory.CYP_PANEL,
        scientific_definition="Calibrated binary probability of inhibiting human CYP2D6 at AC50 <= 10 uM in recombinant assay.",
        species="Human",
        assay_type="Recombinant Human CYP2D6 Functional Assay",
        raw_unit="probability",
        canonical_unit="probability",
        transformation="Identity (probability in [0.0, 1.0])",
        output_type=OutputType.BINARY_CLASSIFICATION,
        directionality=Directionality.LOWER_BETTER,
        classification_semantics={
            "positive_class": "INHIBITOR",
            "negative_class": "NON_INHIBITOR",
            "decision_threshold": 0.50,
            "potency_threshold_um": 10.0,
        },
        required_metadata=("isoform", "role"),
        comparison_rules={"isoform_must_match": "CYP2D6", "role_must_match": "INHIBITOR"},
        ensemble_compatibility_rules={"required_isoform": "CYP2D6", "required_role": "INHIBITOR"},
    ),
    "CYP3A4 inhibitor": EndpointContract(
        endpoint_id="cyp3a4_inhibitor_prob",
        display_name="CYP3A4 Inhibition Probability",
        category=EndpointCategory.CYP_PANEL,
        scientific_definition="Calibrated binary probability of inhibiting human CYP3A4 at AC50 <= 10 uM in recombinant assay.",
        species="Human",
        assay_type="Recombinant Human CYP3A4 Functional Assay",
        raw_unit="probability",
        canonical_unit="probability",
        transformation="Identity (probability in [0.0, 1.0])",
        output_type=OutputType.BINARY_CLASSIFICATION,
        directionality=Directionality.LOWER_BETTER,
        classification_semantics={
            "positive_class": "INHIBITOR",
            "negative_class": "NON_INHIBITOR",
            "decision_threshold": 0.50,
            "potency_threshold_um": 10.0,
        },
        required_metadata=("isoform", "role"),
        comparison_rules={"isoform_must_match": "CYP3A4", "role_must_match": "INHIBITOR"},
        ensemble_compatibility_rules={"required_isoform": "CYP3A4", "required_role": "INHIBITOR"},
    ),

    # 12-14. CYP SUBSTRATE PANEL (2C9, 2D6, 3A4)
    "CYP2C9 substrate": EndpointContract(
        endpoint_id="cyp2c9_substrate_prob",
        display_name="CYP2C9 Substrate Probability",
        category=EndpointCategory.CYP_PANEL,
        scientific_definition="Binary classification probability of a compound being a metabolic turnover substrate for human CYP2C9.",
        species="Human",
        assay_type="Literature In Vitro CYP2C9 Turnover Annotation",
        raw_unit="probability",
        canonical_unit="probability",
        transformation="Identity (probability in [0.0, 1.0])",
        output_type=OutputType.BINARY_CLASSIFICATION,
        directionality=Directionality.NEUTRAL,
        classification_semantics={
            "positive_class": "SUBSTRATE",
            "negative_class": "NON_SUBSTRATE",
            "decision_threshold": 0.50,
        },
        required_metadata=("isoform", "role"),
        comparison_rules={"isoform_must_match": "CYP2C9", "role_must_match": "SUBSTRATE"},
        ensemble_compatibility_rules={"required_isoform": "CYP2C9", "required_role": "SUBSTRATE"},
    ),
    "CYP2D6 substrate": EndpointContract(
        endpoint_id="cyp2d6_substrate_prob",
        display_name="CYP2D6 Substrate Probability",
        category=EndpointCategory.CYP_PANEL,
        scientific_definition="Binary classification probability of a compound being a metabolic turnover substrate for human CYP2D6.",
        species="Human",
        assay_type="Literature In Vitro CYP2D6 Turnover Annotation",
        raw_unit="probability",
        canonical_unit="probability",
        transformation="Identity (probability in [0.0, 1.0])",
        output_type=OutputType.BINARY_CLASSIFICATION,
        directionality=Directionality.NEUTRAL,
        classification_semantics={
            "positive_class": "SUBSTRATE",
            "negative_class": "NON_SUBSTRATE",
            "decision_threshold": 0.50,
        },
        required_metadata=("isoform", "role"),
        comparison_rules={"isoform_must_match": "CYP2D6", "role_must_match": "SUBSTRATE"},
        ensemble_compatibility_rules={"required_isoform": "CYP2D6", "required_role": "SUBSTRATE"},
    ),
    "CYP3A4 substrate": EndpointContract(
        endpoint_id="cyp3a4_substrate_prob",
        display_name="CYP3A4 Substrate Probability",
        category=EndpointCategory.CYP_PANEL,
        scientific_definition="Binary classification probability of a compound being a metabolic turnover substrate for human CYP3A4.",
        species="Human",
        assay_type="Literature In Vitro CYP3A4 Turnover Annotation",
        raw_unit="probability",
        canonical_unit="probability",
        transformation="Identity (probability in [0.0, 1.0])",
        output_type=OutputType.BINARY_CLASSIFICATION,
        directionality=Directionality.NEUTRAL,
        classification_semantics={
            "positive_class": "SUBSTRATE",
            "negative_class": "NON_SUBSTRATE",
            "decision_threshold": 0.50,
        },
        required_metadata=("isoform", "role"),
        comparison_rules={"isoform_must_match": "CYP3A4", "role_must_match": "SUBSTRATE"},
        ensemble_compatibility_rules={"required_isoform": "CYP3A4", "required_role": "SUBSTRATE"},
    ),

    # 15. TRANSPORTERS: P-GP INHIBITOR
    "P-gp inhibitor": EndpointContract(
        endpoint_id="transporter_pgp_inhibitor_prob",
        display_name="P-glycoprotein (P-gp/ABCB1) Inhibitor Probability",
        category=EndpointCategory.TRANSPORTER,
        scientific_definition="Binary probability of inhibiting human P-glycoprotein (P-gp/ABCB1) efflux activity (IC50 <= 15 uM or >25% transport inhibition in calcein-AM / rhodamine-123 dye efflux assays).",
        species="Human",
        assay_type="Heterogeneous Human P-gp Transwell / Fluorescent Dye Efflux Assays",
        raw_unit="probability",
        canonical_unit="probability",
        transformation="Identity (probability in [0.0, 1.0])",
        output_type=OutputType.BINARY_CLASSIFICATION,
        directionality=Directionality.LOWER_BETTER,
        classification_semantics={
            "positive_class": "INHIBITOR",
            "negative_class": "NON_INHIBITOR",
            "decision_threshold": 0.50,
        },
        required_metadata=("transporter", "role"),
        comparison_rules={"transporter_must_match": "P-gp", "role_must_match": "INHIBITOR"},
        ensemble_compatibility_rules={"required_transporter": "P-gp / ABCB1", "required_role": "INHIBITOR"},
    ),

    # 16. HERG CARDIAC LIABILITY
    "hERG liability": EndpointContract(
        endpoint_id="safety_herg_blocker_prob",
        display_name="hERG Cardiac Blocker Liability Probability",
        category=EndpointCategory.SAFETY,
        scientific_definition="Binary probability of exhibiting human ether-a-go-go-related gene (hERG / Kv11.1 / KCNH2) potassium channel blocker liability (screening threshold IC50 <= 10 uM or pIC50 >= 5.0).",
        species="Human",
        assay_type="Heterogeneous Patch-Clamp & Radioligand Binding ([3H]-dofetilide / [3H]-astemizole)",
        raw_unit="probability",
        canonical_unit="probability",
        transformation="Identity (probability in [0.0, 1.0])",
        output_type=OutputType.BINARY_CLASSIFICATION,
        directionality=Directionality.LOWER_BETTER,
        classification_semantics={
            "positive_class": "BLOCKER",
            "negative_class": "NON_BLOCKER",
            "decision_threshold": 0.50,
            "potency_threshold_um": 10.0,
        },
        required_metadata=("target", "species"),
        comparison_rules={"target_must_match": "hERG", "species_must_match": "Human"},
        ensemble_compatibility_rules={
            "required_target": "hERG / KCNH2",
            "required_canonical_unit": "probability",
            "forbidden_mix_types": ["patch_clamp_ic50_regression_without_thresholding", "qt_prolongation_clinical"],
        },
    ),

    # 17. AMES MUTAGENICITY
    "Ames mutagenicity": EndpointContract(
        endpoint_id="safety_ames_mutagenicity_prob",
        display_name="Ames Bacterial Reverse Mutation Mutagenicity",
        category=EndpointCategory.SAFETY,
        scientific_definition="Binary probability of inducing bacterial reverse mutation in Salmonella typhimurium strains (TA98, TA100, TA1535, TA1537, TA102) with or without S9 metabolic activation.",
        species="Salmonella typhimurium",
        assay_type="Bacterial Reverse Mutation Test (Ames Assay / OECD 471)",
        raw_unit="probability",
        canonical_unit="probability",
        transformation="Identity (probability in [0.0, 1.0])",
        output_type=OutputType.BINARY_CLASSIFICATION,
        directionality=Directionality.LOWER_BETTER,
        classification_semantics={
            "positive_class": "MUTAGENIC",
            "negative_class": "NON_MUTAGENIC",
            "decision_threshold": 0.50,
        },
        required_metadata=("strains", "s9_activation"),
        comparison_rules={"assay_family_must_match": "Ames"},
        ensemble_compatibility_rules={
            "required_assay": "Ames Reverse Mutation",
            "required_canonical_unit": "probability",
            "forbidden_mix_types": ["mammalian_micronucleus", "chromosomal_aberration"],
        },
    ),

    # 18. DILI CLINICAL LIABILITY
    "DILI clinical liability": EndpointContract(
        endpoint_id="safety_dili_clinical_prob",
        display_name="Drug-Induced Liver Injury (DILI) Clinical Concern",
        category=EndpointCategory.SAFETY,
        scientific_definition="Binary probability of clinical drug-induced liver injury (DILI) concern based on FDA-NCTR Liver Toxicity Knowledge Base (LTKB) annotations. Represents human clinical post-marketing liver injury risk.",
        species="Human",
        assay_type="Clinical Post-Marketing Pharmacovigilance & Physician Assessment (FDA-NCTR LTKB)",
        raw_unit="probability",
        canonical_unit="probability",
        transformation="Identity (probability in [0.0, 1.0])",
        output_type=OutputType.BINARY_CLASSIFICATION,
        directionality=Directionality.LOWER_BETTER,
        classification_semantics={
            "positive_class": "DILI_CONCERN",
            "negative_class": "NO_DILI_CONCERN",
            "decision_threshold": 0.50,
        },
        required_metadata=("annotation_source", "clinical_severity"),
        comparison_rules={"clinical_scope_must_match": "Human Drug-Level DILI"},
        ensemble_compatibility_rules={
            "required_scope": "Clinical DILI",
            "required_canonical_unit": "probability",
            "forbidden_mix_types": ["in_vitro_cytotoxicity", "isolated_mitochondrial_toxicity", "bsep_inhibition"],
        },
    ),

    # 19. METABOLIC SOFT SPOTS (SoM)
    "Metabolic soft spots": EndpointContract(
        endpoint_id="som_metabolic_soft_spots",
        display_name="Metabolic Site-of-Metabolism (SoM) Ranking",
        category=EndpointCategory.METABOLISM,
        scientific_definition="Atom-level ranking of vulnerability to Phase I (oxidation, dealkylation, hydrolysis) and Phase II (glucuronidation, sulfation) metabolic transformations.",
        species="Human",
        assay_type="Rule-Based SMARTS Reaction Library / Atom Reactivity Ranking (SyGMa)",
        raw_unit="atom_rank_and_score",
        canonical_unit="atom_index_ranking",
        transformation="Rank-ordered atom indices with reaction pathway metadata",
        output_type=OutputType.RANKING,
        directionality=Directionality.NEUTRAL,
        required_metadata=("phase_1_reactions", "phase_2_reactions", "atom_indices"),
        comparison_rules={"atom_mapping_must_match": True},
        ensemble_compatibility_rules={
            "ensemble_aggregation_method": "Reciprocal Rank Fusion (RRF) / Borda Count",
            "forbidden_mix_types": ["numeric_average_of_uncalibrated_scores"],
        },
    ),

    # 20. PHARMACOKINETIC CLEARANCE
    "PK Systemic Clearance": EndpointContract(
        endpoint_id="pk_clearance_systemic",
        display_name="Systemic In Vivo Clearance (CL)",
        category=EndpointCategory.PHARMACOKINETICS,
        scientific_definition="Total systemic in vivo clearance of drug from circulating plasma, in mL/min/kg.",
        species="Multi-Species (Mouse, Rat, Dog, Monkey, Human)",
        assay_type="In Vivo Intravenous PK Non-Compartmental Analysis (NCA)",
        raw_unit="mL/min/kg",
        canonical_unit="mL/min/kg",
        transformation="Identity (mL/min/kg)",
        output_type=OutputType.REGRESSION,
        directionality=Directionality.LOWER_BETTER,
        required_metadata=("species", "route", "matrix", "dose_mg_per_kg"),
        experimental_compatibility_rules={
            "compatible_types": ["In Vivo IV NCA Clearance CL"],
            "incompatible_types": ["Apparent Oral Clearance CL/F (unless normalized by F)", "In Vitro Clint (requires IVIVE scaling)"],
        },
        comparison_rules={"species_must_match": True, "route_must_match": "IV"},
        ensemble_compatibility_rules={
            "hierarchy": "EXPERIMENTAL_NCA > HEPATIC_IVIVE > ALLOMETRIC_SCALING > QSAR_PREDICTION",
            "forbidden_mix_types": ["apparent_oral_cl_f_as_iv_cl"],
        },
    ),

    # 21. PHARMACOKINETIC VOLUME OF DISTRIBUTION
    "PK Volume of Distribution": EndpointContract(
        endpoint_id="pk_volume_distribution_vss",
        display_name="Steady-State Volume of Distribution (Vss)",
        category=EndpointCategory.PHARMACOKINETICS,
        scientific_definition="Apparent steady-state volume of distribution reflecting the extent of tissue partitioning, in L/kg.",
        species="Multi-Species (Mouse, Rat, Dog, Monkey, Human)",
        assay_type="In Vivo Intravenous PK NCA / Mechanistic Tissue Ionization Model",
        raw_unit="L/kg",
        canonical_unit="L/kg",
        transformation="Identity (L/kg)",
        output_type=OutputType.REGRESSION,
        directionality=Directionality.NEUTRAL,
        required_metadata=("species", "route", "matrix", "ionization_class"),
        experimental_compatibility_rules={
            "compatible_types": ["In Vivo IV NCA Vss", "In Vivo IV NCA Vz"],
            "incompatible_types": ["Oral Vz/F (unadjusted for F)"],
        },
        comparison_rules={"species_must_match": True},
        ensemble_compatibility_rules={
            "hierarchy": "EXPERIMENTAL_VSS > EXPERIMENTAL_VZ > LOMBARDO_IONIZATION_VD > QSAR_PREDICTION",
        },
    ),

    # 22. PHARMACOKINETIC BIOAVAILABILITY
    "PK Bioavailability": EndpointContract(
        endpoint_id="pk_bioavailability_fraction",
        display_name="Oral Absolute Bioavailability (F)",
        category=EndpointCategory.PHARMACOKINETICS,
        scientific_definition="Fraction of an orally administered dose that reaches the systemic circulation intact, expressed as % (0-100%).",
        species="Multi-Species (Mouse, Rat, Dog, Monkey, Human)",
        assay_type="Matched IV and PO In Vivo PK Study AUC Normalization",
        raw_unit="%",
        canonical_unit="%",
        transformation="F (%) = (AUC_PO * Dose_IV) / (AUC_IV * Dose_PO) * 100",
        output_type=OutputType.REGRESSION,
        directionality=Directionality.HIGHER_BETTER,
        required_metadata=("species", "dose_iv", "dose_po", "auc_iv", "auc_po"),
        experimental_compatibility_rules={
            "compatible_types": ["Matched IV/PO in vivo experimental F"],
            "incompatible_types": ["Unmatched single-dose oral AUC"],
        },
        comparison_rules={"species_must_match": True},
        ensemble_compatibility_rules={
            "hierarchy": "MATCHED_EXPERIMENTAL_F > MECHANISTIC_DECOMPOSITION (F = Fa * Fg * Fh) > UNAVAILABLE",
        },
    ),
}


# ==============================================================================
# ENSEMBLE COMPATIBILITY & VALIDATION ENGINE
# ==============================================================================

def get_endpoint_contract(endpoint_name: str) -> Optional[EndpointContract]:
    """Retrieve the authoritative contract for a given endpoint."""
    return ENDPOINT_CONTRACTS.get(endpoint_name)


def check_ensemble_compatibility(
    contract_a: EndpointContract,
    contract_b: EndpointContract,
) -> Tuple[bool, str]:
    """
    Deterministically check if two endpoint contracts can be ensembled.
    Returns (is_compatible, reason).
    """
    # 1. Endpoint ID must match
    if contract_a.endpoint_id != contract_b.endpoint_id:
        return False, f"Incompatible endpoint IDs: '{contract_a.endpoint_id}' vs '{contract_b.endpoint_id}'"

    # 2. Species must match
    if contract_a.species.lower() != contract_b.species.lower():
        return False, f"Incompatible species: '{contract_a.species}' vs '{contract_b.species}'"

    # 3. Output types must match
    if contract_a.output_type != contract_b.output_type:
        return False, f"Incompatible output types: {contract_a.output_type.value} vs {contract_b.output_type.value}"

    # 4. Canonical units must match
    if contract_a.canonical_unit != contract_b.canonical_unit:
        return False, f"Incompatible canonical units: '{contract_a.canonical_unit}' vs '{contract_b.canonical_unit}'"

    # 5. Role/Isoform compatibility for CYP and Transporters
    if contract_a.category in (EndpointCategory.CYP_PANEL, EndpointCategory.TRANSPORTER):
        sem_a = contract_a.classification_semantics or {}
        sem_b = contract_b.classification_semantics or {}
        if sem_a.get("positive_class") != sem_b.get("positive_class"):
            return False, f"Incompatible classification semantics: {sem_a.get('positive_class')} vs {sem_b.get('positive_class')}"

    return True, "Ensemble compatible: identical endpoint semantics, species, and canonical units."


# ==============================================================================
# FUTURE COMMON MODEL ADAPTER CONTRACT & RESULT SCHEMA
# ==============================================================================

@dataclass
class AdapterPredictionResult:
    """
    Standardized result payload returned by every qualified model adapter.
    """
    model_id: str
    model_version: str
    model_family: str
    endpoint_id: str
    canonical_unit: str
    value: float
    raw_value: Optional[float] = None
    probability: Optional[float] = None
    classification: Optional[str] = None
    applicability_domain: str = "IN_DOMAIN"  # IN_DOMAIN, BORDERLINE, OUT_OF_DOMAIN
    applicability_distance: float = 0.0
    confidence: str = "MEDIUM"               # HIGH, MEDIUM, LOW, NOT_APPLICABLE
    conformal_interval: Optional[Tuple[float, float]] = None
    conformal_set: Optional[List[str]] = None
    runtime_ms: float = 0.0
    execution_tier: ExecutionTier = ExecutionTier.TIER_1_LOCAL_FAST
    warnings: List[str] = field(default_factory=list)
    provenance: Dict[str, Any] = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    def to_dict(self) -> Dict[str, Any]:
        return {
            "model_id": self.model_id,
            "model_version": self.model_version,
            "model_family": self.model_family,
            "endpoint_id": self.endpoint_id,
            "canonical_unit": self.canonical_unit,
            "value": self.value,
            "raw_value": self.raw_value,
            "probability": self.probability,
            "classification": self.classification,
            "applicability_domain": self.applicability_domain,
            "applicability_distance": self.applicability_distance,
            "confidence": self.confidence,
            "conformal_interval": list(self.conformal_interval) if self.conformal_interval else None,
            "conformal_set": self.conformal_set,
            "runtime_ms": self.runtime_ms,
            "execution_tier": self.execution_tier.value,
            "warnings": self.warnings,
            "provenance": self.provenance,
            "timestamp": self.timestamp,
        }


class BaseModelAdapter:
    """
    Abstract base interface for all future multi-model ensemble adapters.
    """
    model_id: str
    model_version: str
    model_family: str
    supported_endpoints: Set[str]
    execution_tier: ExecutionTier
    arm64_status: ARM64Status

    def is_available(self) -> Tuple[bool, str]:
        """Check if model checkpoint and runtime dependencies are available."""
        raise NotImplementedError

    def predict(
        self,
        canonical_smiles: str,
        endpoint_contract: EndpointContract,
        compound_metadata: Optional[Dict[str, Any]] = None,
    ) -> AdapterPredictionResult:
        """
        Execute deterministic inference and return standard AdapterPredictionResult.
        Must enforce endpoint contract validation before execution.
        """
        raise NotImplementedError


# ==============================================================================
# FAILURE ISOLATION & ENSEMBLE AGGREGATION ENGINE
# ==============================================================================

@dataclass
class EnsembleExecutionResult:
    """
    Result of multi-model execution with fault-tolerant member aggregation.
    """
    endpoint_id: str
    successful_predictions: List[AdapterPredictionResult]
    failed_models: List[Dict[str, str]]
    total_models_attempted: int
    is_valid: bool
    status_summary: str

    @property
    def member_count(self) -> int:
        return len(self.successful_predictions)


def execute_fault_tolerant_ensemble(
    adapters: List[BaseModelAdapter],
    canonical_smiles: str,
    contract: EndpointContract,
    compound_metadata: Optional[Dict[str, Any]] = None,
) -> EnsembleExecutionResult:
    """
    Executes multiple model adapters with strict failure isolation:
    If Model 1 succeeds, Model 2 succeeds, and Model 3 fails (raises exception or returns error),
    the endpoint continues cleanly with Model 1 and Model 2 without crashing.
    """
    successful: List[AdapterPredictionResult] = []
    failed: List[Dict[str, str]] = []

    for adapter in adapters:
        if contract.endpoint_id not in adapter.supported_endpoints:
            failed.append({
                "model_id": getattr(adapter, "model_id", "unknown"),
                "reason": f"Model does not support endpoint '{contract.endpoint_id}'",
            })
            continue

        try:
            avail, reason = adapter.is_available()
            if not avail:
                failed.append({
                    "model_id": adapter.model_id,
                    "reason": f"Model unavailable: {reason}",
                })
                continue

            result = adapter.predict(canonical_smiles, contract, compound_metadata)
            successful.append(result)
        except Exception as exc:
            failed.append({
                "model_id": getattr(adapter, "model_id", "unknown"),
                "reason": f"Inference execution failure: {type(exc).__name__}: {str(exc)}",
            })

    is_valid = len(successful) > 0
    status = (
        f"SUCCESS ({len(successful)}/{len(adapters)} models passed)"
        if is_valid else
        f"FAILED (0/{len(adapters)} models passed)"
    )

    return EnsembleExecutionResult(
        endpoint_id=contract.endpoint_id,
        successful_predictions=successful,
        failed_models=failed,
        total_models_attempted=len(adapters),
        is_valid=is_valid,
        status_summary=status,
    )
