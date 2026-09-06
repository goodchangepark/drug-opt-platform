"""Versioned semantic endpoint and unit registry for experiment/prediction joins.

Raw source labels are intentionally not used as comparison keys. This module
is the authoritative source where a persisted observation or prediction is translated
to the canonical scientific endpoint used by the comparison API.
"""
from __future__ import annotations

import math
import re
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Tuple

from backend.unit_normalization import (
    convert_concentration, convert_solubility, convert_ppb,
    convert_cyp_herg_ic50, convert_pic50_to_ic50, convert_clearance,
    convert_volume, convert_time, convert_caco2_papp, convert_auc,
    normalize_scientific_observation, clean_unit_str,
    TYPE_DETERMINISTIC, TYPE_SCALE, TYPE_DIRECT_IDENTITY, TYPE_IVIVE_DERIVED
)


CANONICAL_ENDPOINT_VERSION = "drugopt-canonical-endpoint-v2"
COMPARISON_UNIT_VERSION = "drugopt-comparison-unit-v2"
EXPERIMENTAL_NORMALIZATION_VERSION = "drugopt-experimental-normalization-v2"

# Prediction output provenance
PREDICTION_MODEL = "MODEL"
PREDICTION_MECHANISTIC = "MECHANISTIC_ESTIMATE"
PREDICTION_RULE = "RULE_ESTIMATE"
PREDICTION_DERIVED = "DERIVED_ESTIMATE"
PREDICTION_UNAVAILABLE = "MODEL_UNAVAILABLE"

# Four Authoritative Comparison Status Classes
DIRECT = "DIRECTLY_COMPARABLE"
CONVERTED = "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION"
CONDITIONAL = "CONDITIONALLY_COMPARABLE"
RELATED = "RELATED_NOT_SAME_ENDPOINT"
UNSUPPORTED = "UNSUPPORTED"
NOT_COMPARABLE = "NOT_COMPARABLE"
NO_MODEL = "NO_MATCHING_PREDICTION_MODEL"


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
    domain: str = ""
    allowed_relations: tuple[str, ...] = ("=", "<", ">", "<=", ">=", "~")
    prediction_endpoint_aliases: tuple[str, ...] = ()
    experimental_endpoint_aliases: tuple[str, ...] = ()


def _ep(id_, section, label, definition, value_type, unit, scale, **kwargs):
    return CanonicalEndpoint(id_, section, label, definition, value_type, unit, scale, **kwargs)


REGISTRY: dict[str, CanonicalEndpoint] = {
    # -----------------------------------------------------------------
    # Group 1: Physicochemical Properties (10 endpoints)
    # -----------------------------------------------------------------
    "MW": _ep("MW", "PROPERTIES", "Molecular Weight", "molecular mass of parent structure", "numeric", "g/mol", "LINEAR", domain="physicochemical", prediction_endpoint_aliases=("mw", "molecular weight", "molecular_weight")),
    "CLOGP": _ep("CLOGP", "PROPERTIES", "cLogP", "calculated octanol-water partition coefficient", "numeric", "logP", "LINEAR", domain="physicochemical", prediction_endpoint_aliases=("clogp", "logp", "calc_logp"), experimental_endpoint_aliases=("logp", "clogp")),
    "TPSA": _ep("TPSA", "PROPERTIES", "TPSA", "topological polar surface area", "numeric", "Å²", "LINEAR", domain="physicochemical", prediction_endpoint_aliases=("tpsa", "polar surface area")),
    "HBD": _ep("HBD", "PROPERTIES", "H-Bond Donors", "count of hydrogen bond donors", "integer", "count", "LINEAR", domain="physicochemical", prediction_endpoint_aliases=("hbd", "h_bond_donors")),
    "HBA": _ep("HBA", "PROPERTIES", "H-Bond Acceptors", "count of hydrogen bond acceptors", "integer", "count", "LINEAR", domain="physicochemical", prediction_endpoint_aliases=("hba", "h_bond_acceptors")),
    "ROTB": _ep("ROTB", "PROPERTIES", "Rotatable Bonds", "count of rotatable single bonds", "integer", "count", "LINEAR", domain="physicochemical", prediction_endpoint_aliases=("rotb", "rotatable_bonds")),
    "FSP3": _ep("FSP3", "PROPERTIES", "Fsp3", "fraction of sp3-hybridized carbon atoms", "numeric", "fraction", "LINEAR", domain="physicochemical", prediction_endpoint_aliases=("fsp3", "fraction_sp3")),
    "QED": _ep("QED", "PROPERTIES", "QED Score", "quantitative estimate of drug-likeness", "numeric", "score (0-1)", "LINEAR", domain="physicochemical", prediction_endpoint_aliases=("qed", "drug_likeness")),
    "FORMAL_CHARGE": _ep("FORMAL_CHARGE", "PROPERTIES", "Formal Charge", "net molecular formal charge", "integer", "charge", "LINEAR", domain="physicochemical", prediction_endpoint_aliases=("formal_charge", "net_charge")),
    "HEAVY_ATOM_COUNT": _ep("HEAVY_ATOM_COUNT", "PROPERTIES", "Heavy Atom Count", "total count of non-hydrogen atoms", "integer", "count", "LINEAR", domain="physicochemical", prediction_endpoint_aliases=("heavy_atom_count", "num_heavy_atoms")),

    # -----------------------------------------------------------------
    # Group 2: Solubility & Permeability (7 endpoints)
    # -----------------------------------------------------------------
    "SOLUBILITY_GENERIC": _ep("SOLUBILITY_GENERIC", "ADMET", "Solubility", "aqueous solubility without a validated solid-state subtype", "numeric", "log10(mol/L)", "LOG10", domain="absorption", experimental_endpoint_aliases=("solubility", "aqueous solubility", "water solubility", "aqueous_solubility"), prediction_endpoint_aliases=("solubility", "aqueous solubility")),
    "SOLUBILITY_KINETIC": _ep("SOLUBILITY_KINETIC", "ADMET", "Kinetic solubility", "kinetic aqueous solubility", "numeric", "log10(mol/L)", "LOG10", domain="absorption", assay_requirement="kinetic", experimental_endpoint_aliases=("kinetic solubility",)),
    "SOLUBILITY_THERMODYNAMIC": _ep("SOLUBILITY_THERMODYNAMIC", "ADMET", "Thermodynamic solubility", "thermodynamic aqueous solubility", "numeric", "log10(mol/L)", "LOG10", domain="absorption", assay_requirement="thermodynamic", experimental_endpoint_aliases=("thermodynamic solubility",)),
    "SOLUBILITY_INTRINSIC": _ep("SOLUBILITY_INTRINSIC", "ADMET", "Intrinsic solubility", "intrinsic solubility of the neutral form", "numeric", "log10(mol/L)", "LOG10", domain="absorption", assay_requirement="intrinsic", experimental_endpoint_aliases=("intrinsic solubility",)),
    "CACO2_PAPP_AB": _ep("CACO2_PAPP_AB", "ADMET", "Caco-2 Papp A→B", "Caco-2 apparent permeability from apical to basolateral", "numeric", "log10(cm/s)", "LOG10", domain="absorption", matrix_requirement="Caco-2", direction_requirement="A→B", experimental_endpoint_aliases=("caco-2 permeability", "caco2 permeability", "caco-2 papp a-b", "papp a-b", "caco2", "caco-2"), prediction_endpoint_aliases=("permeability", "caco-2 permeability", "caco2")),
    "CACO2_PAPP_BA": _ep("CACO2_PAPP_BA", "ADMET", "Caco-2 Papp B→A", "Caco-2 apparent permeability from basolateral to apical", "numeric", "log10(cm/s)", "LOG10", domain="absorption", matrix_requirement="Caco-2", direction_requirement="B→A", experimental_endpoint_aliases=("caco-2 papp b-a", "papp b-a")),
    "CACO2_EFFLUX_RATIO": _ep("CACO2_EFFLUX_RATIO", "ADMET", "Caco-2 efflux ratio", "ratio of B→A to A→B permeability", "numeric", "ratio", "LINEAR", domain="absorption", matrix_requirement="Caco-2", experimental_endpoint_aliases=("efflux ratio", "caco2 efflux ratio")),

    # -----------------------------------------------------------------
    # Group 3: Plasma Protein Binding & fu (4 endpoints)
    # -----------------------------------------------------------------
    "HUMAN_PPB": _ep("HUMAN_PPB", "ADMET", "Human plasma protein binding", "fraction of compound bound to human plasma protein", "numeric", "% bound", "PERCENT", domain="distribution", species_requirement="HUMAN", matrix_requirement="plasma", experimental_endpoint_aliases=("ppb", "plasma protein binding", "protein binding", "fraction bound", "fu", "human ppb", "human plasma protein binding"), prediction_endpoint_aliases=("plasma protein binding", "ppb", "human ppb")),
    "RAT_PPB": _ep("RAT_PPB", "ADMET", "Rat plasma protein binding", "fraction bound to rat plasma protein", "numeric", "% bound", "PERCENT", domain="distribution", species_requirement="RAT", matrix_requirement="plasma", experimental_endpoint_aliases=("rat ppb", "rat plasma protein binding")),
    "MOUSE_PPB": _ep("MOUSE_PPB", "ADMET", "Mouse plasma protein binding", "fraction bound to mouse plasma protein", "numeric", "% bound", "PERCENT", domain="distribution", species_requirement="MOUSE", matrix_requirement="plasma", experimental_endpoint_aliases=("mouse ppb", "mouse plasma protein binding")),
    "PPB_UNSPECIFIED": _ep("PPB_UNSPECIFIED", "ADMET", "Plasma protein binding (unspecified)", "plasma protein binding with unresolved species context", "numeric", "% bound", "PERCENT", domain="distribution"),

    # -----------------------------------------------------------------
    # Group 4: Clearance & Microsomal Metabolism (11 endpoints)
    # -----------------------------------------------------------------
    "HLM_CLINT": _ep("HLM_CLINT", "ADMET", "HLM intrinsic clearance", "human liver microsomal intrinsic clearance", "numeric", "log10(mL/min/kg)", "LOG10", domain="clearance", species_requirement="HUMAN", matrix_requirement="microsomes", experimental_endpoint_aliases=("hlm clint", "hlm", "human liver microsomal clearance", "human microsomal intrinsic clearance", "hlm clearance"), prediction_endpoint_aliases=("hlm intrinsic clearance", "hlm clint", "hlm")),
    "RLM_CLINT": _ep("RLM_CLINT", "ADMET", "RLM intrinsic clearance", "rat liver microsomal intrinsic clearance", "numeric", "log10(mL/min/kg)", "LOG10", domain="clearance", species_requirement="RAT", matrix_requirement="microsomes", experimental_endpoint_aliases=("rlm clint", "rlm", "rat liver microsomal clearance", "rat microsomal intrinsic clearance"), prediction_endpoint_aliases=("rlm intrinsic clearance", "rlm clint", "rlm")),
    "MLM_CLINT": _ep("MLM_CLINT", "ADMET", "MLM intrinsic clearance", "mouse liver microsomal intrinsic clearance", "numeric", "log10(mL/min/kg)", "LOG10", domain="clearance", species_requirement="MOUSE", matrix_requirement="microsomes", experimental_endpoint_aliases=("mlm clint", "mlm", "mouse liver microsomal clearance", "mouse microsomal intrinsic clearance"), prediction_endpoint_aliases=("mlm intrinsic clearance", "mlm clint", "mlm")),
    "HEPATOCYTE_CLINT": _ep("HEPATOCYTE_CLINT", "METABOLISM", "Hepatocyte intrinsic clearance", "intrinsic clearance measured in intact hepatocytes", "numeric", "µL/min/10^6 cells", "LINEAR", domain="clearance", matrix_requirement="hepatocytes", experimental_endpoint_aliases=("hepatocyte", "hepatocyte clint", "clh", "hepatocytes")),
    "METABOLIC_SOFT_SPOTS": _ep("METABOLIC_SOFT_SPOTS", "METABOLISM", "Metabolic soft spots", "ranked atom-level metabolic transformation hypotheses", "ranking", "ranked sites", "RANKING", domain="metabolism", prediction_endpoint_aliases=("soft spots", "metabolic soft spots")),
    "METABOLITE_HYPOTHESES": _ep("METABOLITE_HYPOTHESES", "METABOLISM", "Metabolite hypotheses", "rule-generated predicted metabolite structures", "ranking", "hypotheses", "RANKING", domain="metabolism", prediction_endpoint_aliases=("metabolite hypotheses", "predicted metabolites")),
    "METABOLITE_OBSERVATION": _ep("METABOLITE_OBSERVATION", "METABOLISM", "Metabolite observation", "observed metabolite identity or exposure in vivo", "numeric_or_qualitative", "", "MIXED", domain="metabolism", experimental_endpoint_aliases=("metabolite", "metabolites", "metabolite observation")),
    "EXCRETION_FECAL": _ep("EXCRETION_FECAL", "METABOLISM", "Fecal excretion", "dose recovered in feces", "numeric", "% dose", "PERCENT", domain="excretion", experimental_endpoint_aliases=("feces", "fecal excretion")),
    "EXCRETION_URINARY": _ep("EXCRETION_URINARY", "METABOLISM", "Urinary excretion", "dose recovered in urine", "numeric", "% dose", "PERCENT", domain="excretion", experimental_endpoint_aliases=("urine", "urinary excretion")),
    "CYP3A4_METABOLIC_CONTRIBUTION": _ep("CYP3A4_METABOLIC_CONTRIBUTION", "METABOLISM", "CYP3A4 metabolic contribution", "fraction metabolized by CYP3A4 (fm CYP3A4)", "numeric_or_qualitative", "%", "PERCENT", domain="metabolism", experimental_endpoint_aliases=("cyp3a", "cyp3a4", "fm cyp3a", "fm_cyp3a4")),
    "CYP2D6_METABOLIC_CONTRIBUTION": _ep("CYP2D6_METABOLIC_CONTRIBUTION", "METABOLISM", "CYP2D6 metabolic contribution", "fraction metabolized by CYP2D6 (fm CYP2D6)", "numeric_or_qualitative", "%", "PERCENT", domain="metabolism", experimental_endpoint_aliases=("cyp2d6", "fm cyp2d6", "fm_cyp2d6")),

    # -----------------------------------------------------------------
    # Group 5: Ionization & Distribution (4 endpoints)
    # -----------------------------------------------------------------
    "PKA": _ep("PKA", "ADMET", "pKa", "acid/base dissociation constant", "numeric", "pKa", "LINEAR", domain="physicochemical", experimental_endpoint_aliases=("pka",), prediction_endpoint_aliases=("pka", "pka (quantitative ml)")),
    "PKA_ACID": _ep("PKA_ACID", "ADMET", "Acidic pKa", "strongest acidic dissociation constant (macro pKa)", "numeric", "pKa", "LINEAR", domain="physicochemical", experimental_endpoint_aliases=("acidic pka", "pka_acid", "pka acid", "pka (acidic)")),
    "PKA_BASE": _ep("PKA_BASE", "ADMET", "Basic pKa", "strongest basic dissociation constant (macro pKa)", "numeric", "pKa", "LINEAR", domain="physicochemical", experimental_endpoint_aliases=("basic pka", "pka_base", "pka base", "pka (basic)")),
    "LOGD_7_4": _ep("LOGD_7_4", "ADMET", "logD 7.4", "octanol-water distribution coefficient at pH 7.4", "numeric", "logD", "LINEAR", domain="physicochemical", experimental_endpoint_aliases=("logd", "logd7.4", "logd 7.4"), prediction_endpoint_aliases=("logd7.4", "logd7.4 (quantitative ml)")),
    "LOGP_RELATED": _ep("LOGP_RELATED", "ADMET", "logP (related)", "partition coefficient; distinct from logD 7.4", "numeric", "logP", "LINEAR", domain="physicochemical", experimental_endpoint_aliases=("logp",)),
    "VDSS": _ep("VDSS", "ADMET", "Volume of distribution at steady state (Vdss)", "human steady-state volume of distribution", "numeric", "L/kg", "LINEAR", domain="distribution", species_requirement="HUMAN", prediction_endpoint_aliases=("vdss", "vss", "vd")),

    # -----------------------------------------------------------------
    # Group 6: CYP Panel (Quantitative pIC50 + Classifier + Substrates) (13 endpoints)
    # -----------------------------------------------------------------
    "CYP1A2_INHIBITION": _ep("CYP1A2_INHIBITION", "METABOLISM", "CYP1A2 inhibition", "CYP1A2 inhibition potency (pIC50)", "numeric", "pIC50", "PIC50", domain="metabolism", experimental_endpoint_aliases=("cyp1a2 ic50", "cyp1a2 inhibition", "cyp1a2"), prediction_endpoint_aliases=("cyp1a2 inhibitor", "cyp1a2_inhibition")),
    "CYP1A2_INHIBITOR_CLASS": _ep("CYP1A2_INHIBITOR_CLASS", "METABOLISM", "CYP1A2 inhibitor class", "probability of CYP1A2 inhibition (>50% at 10 µM)", "numeric", "probability", "PROBABILITY", domain="metabolism", prediction_endpoint_aliases=("cyp1a2 inhibitor class",)),
    "CYP2C9_INHIBITION": _ep("CYP2C9_INHIBITION", "METABOLISM", "CYP2C9 inhibition", "CYP2C9 inhibition potency (pIC50)", "numeric", "pIC50", "PIC50", domain="metabolism", experimental_endpoint_aliases=("cyp2c9 ic50", "cyp2c9 inhibition", "cyp2c9"), prediction_endpoint_aliases=("cyp2c9 inhibitor", "cyp2c9_inhibition")),
    "CYP2C9_INHIBITOR_CLASS": _ep("CYP2C9_INHIBITOR_CLASS", "METABOLISM", "CYP2C9 inhibitor class", "probability of CYP2C9 inhibition", "numeric", "probability", "PROBABILITY", domain="metabolism", prediction_endpoint_aliases=("cyp2c9 inhibitor class",)),
    "CYP2C19_INHIBITION": _ep("CYP2C19_INHIBITION", "METABOLISM", "CYP2C19 inhibition", "CYP2C19 inhibition potency", "numeric", "pIC50", "PIC50", domain="metabolism", experimental_endpoint_aliases=("cyp2c19 ic50", "cyp2c19 inhibition", "cyp2c19"), prediction_endpoint_aliases=("cyp2c19 inhibitor", "cyp2c19_inhibition")),
    "CYP2C19_INHIBITOR_CLASS": _ep("CYP2C19_INHIBITOR_CLASS", "METABOLISM", "CYP2C19 inhibitor class", "probability of CYP2C19 inhibition", "numeric", "probability", "PROBABILITY", domain="metabolism", prediction_endpoint_aliases=("cyp2c19 inhibitor class",)),
    "CYP2D6_INHIBITION": _ep("CYP2D6_INHIBITION", "METABOLISM", "CYP2D6 inhibition", "CYP2D6 inhibition potency (pIC50)", "numeric", "pIC50", "PIC50", domain="metabolism", experimental_endpoint_aliases=("cyp2d6 ic50", "cyp2d6 inhibition", "cyp2d6"), prediction_endpoint_aliases=("cyp2d6 inhibitor", "cyp2d6_inhibition")),
    "CYP2D6_INHIBITOR_CLASS": _ep("CYP2D6_INHIBITOR_CLASS", "METABOLISM", "CYP2D6 inhibitor class", "probability of CYP2D6 inhibition", "numeric", "probability", "PROBABILITY", domain="metabolism", prediction_endpoint_aliases=("cyp2d6 inhibitor class",)),
    "CYP3A4_INHIBITION": _ep("CYP3A4_INHIBITION", "METABOLISM", "CYP3A4 inhibition", "CYP3A4 inhibition potency (pIC50)", "numeric", "pIC50", "PIC50", domain="metabolism", experimental_endpoint_aliases=("cyp3a4 ic50", "cyp3a4 inhibition", "cyp3a4", "cyp3a", "cyp3a ic50"), prediction_endpoint_aliases=("cyp3a4 inhibitor", "cyp3a4_inhibition")),
    "CYP3A4_INHIBITOR_CLASS": _ep("CYP3A4_INHIBITOR_CLASS", "METABOLISM", "CYP3A4 inhibitor class", "probability of CYP3A4 inhibition", "numeric", "probability", "PROBABILITY", domain="metabolism", prediction_endpoint_aliases=("cyp3a4 inhibitor class",)),
    "CYP2C9_SUBSTRATE": _ep("CYP2C9_SUBSTRATE", "METABOLISM", "CYP2C9 substrate", "probability of being a CYP2C9 substrate", "numeric", "probability", "PROBABILITY", domain="metabolism", prediction_endpoint_aliases=("cyp2c9 substrate",)),
    "CYP2D6_SUBSTRATE": _ep("CYP2D6_SUBSTRATE", "METABOLISM", "CYP2D6 substrate", "probability of being a CYP2D6 substrate", "numeric", "probability", "PROBABILITY", domain="metabolism", prediction_endpoint_aliases=("cyp2d6 substrate",)),
    "CYP3A4_SUBSTRATE": _ep("CYP3A4_SUBSTRATE", "METABOLISM", "CYP3A4 substrate", "probability of being a CYP3A4 substrate", "numeric", "probability", "PROBABILITY", domain="metabolism", prediction_endpoint_aliases=("cyp3a4 substrate",)),

    # -----------------------------------------------------------------
    # Group 7: Transporters & Permeability Classifiers (12 endpoints)
    # -----------------------------------------------------------------
    "PGP_INHIBITION": _ep("PGP_INHIBITION", "METABOLISM", "P-gp interaction / inhibition", "P-glycoprotein inhibition probability", "numeric", "probability", "PROBABILITY", domain="transporters", experimental_endpoint_aliases=("p-gp", "pgp", "p-gp inhibition", "pgp inhibition"), prediction_endpoint_aliases=("p-gp inhibitor", "pgp inhibitor")),
    "PGP_INHIBITION_QUANT": _ep("PGP_INHIBITION_QUANT", "METABOLISM", "P-gp quantitative IC50", "quantitative P-gp inhibition IC50", "numeric", "µM", "LINEAR", domain="transporters", experimental_endpoint_aliases=("p-gp ic50", "pgp ic50")),
    "PGP_SUBSTRATE": _ep("PGP_SUBSTRATE", "METABOLISM", "P-gp substrate", "probability of being a P-gp substrate", "numeric", "probability", "PROBABILITY", domain="transporters", prediction_endpoint_aliases=("p-gp substrate",)),
    "BCRP_INHIBITION": _ep("BCRP_INHIBITION", "METABOLISM", "BCRP interaction / inhibition", "BCRP / ABCG2 inhibition probability", "numeric", "probability", "PROBABILITY", domain="transporters", experimental_endpoint_aliases=("bcrp", "bcrp inhibition", "bcrp ki", "abcg2")),
    "BCRP_INHIBITOR": _ep("BCRP_INHIBITOR", "METABOLISM", "BCRP inhibitor", "BCRP inhibitor classifier", "numeric", "probability", "PROBABILITY", domain="transporters", prediction_endpoint_aliases=("bcrp inhibitor",)),
    "BCRP_INHIBITOR_QUANT": _ep("BCRP_INHIBITOR_QUANT", "METABOLISM", "BCRP quantitative IC50", "quantitative BCRP inhibition IC50", "numeric", "µM", "LINEAR", domain="transporters", experimental_endpoint_aliases=("bcrp ic50",)),
    "BCRP_SUBSTRATE": _ep("BCRP_SUBSTRATE", "METABOLISM", "BCRP substrate", "probability of being a BCRP substrate", "numeric", "probability", "PROBABILITY", domain="transporters", prediction_endpoint_aliases=("bcrp substrate",)),
    "OATP1B1_INHIBITOR": _ep("OATP1B1_INHIBITOR", "METABOLISM", "OATP1B1 inhibitor", "OATP1B1 transporter inhibition probability", "numeric", "probability", "PROBABILITY", domain="transporters", prediction_endpoint_aliases=("oatp1b1 inhibitor",)),
    "OATP1B3_INHIBITOR": _ep("OATP1B3_INHIBITOR", "METABOLISM", "OATP1B3 inhibitor", "OATP1B3 transporter inhibition probability", "numeric", "probability", "PROBABILITY", domain="transporters", prediction_endpoint_aliases=("oatp1b3 inhibitor",)),
    "OCT1_INHIBITOR": _ep("OCT1_INHIBITOR", "METABOLISM", "OCT1 inhibitor", "OCT1 transporter inhibition probability", "numeric", "probability", "PROBABILITY", domain="transporters", prediction_endpoint_aliases=("oct1 inhibitor",)),
    "OCT2_INHIBITOR": _ep("OCT2_INHIBITOR", "METABOLISM", "OCT2 inhibitor", "OCT2 transporter inhibition probability", "numeric", "probability", "PROBABILITY", domain="transporters", prediction_endpoint_aliases=("oct2 inhibitor",)),
    "HIA": _ep("HIA", "ADMET", "Human Intestinal Absorption (HIA)", "probability of high human intestinal absorption (>30%)", "numeric", "probability", "PROBABILITY", domain="absorption", prediction_endpoint_aliases=("hia", "human intestinal absorption")),
    "BBB_PENETRATION": _ep("BBB_PENETRATION", "ADMET", "Blood-Brain Barrier (BBB)", "probability of blood-brain barrier penetration", "numeric", "probability", "PROBABILITY", domain="distribution", prediction_endpoint_aliases=("bbb", "bbb penetration", "blood brain barrier")),

    # -----------------------------------------------------------------
    # Group 8: Safety & Toxicology (4 endpoints)
    # -----------------------------------------------------------------
    "HERG_LIABILITY": _ep("HERG_LIABILITY", "TOXICITY", "hERG liability", "hERG potassium channel inhibition potency (pIC50)", "numeric", "pIC50", "PIC50", domain="toxicity", experimental_endpoint_aliases=("herg", "herg ic50", "herg liability"), prediction_endpoint_aliases=("herg liability", "herg")),
    "HERG_CLASS": _ep("HERG_CLASS", "TOXICITY", "hERG risk class", "probability of significant hERG inhibition (IC50 < 10 µM)", "numeric", "probability", "PROBABILITY", domain="toxicity", prediction_endpoint_aliases=("herg class",)),
    "AMES_MUTAGENICITY": _ep("AMES_MUTAGENICITY", "TOXICITY", "Ames mutagenicity", "probability of bacterial mutagenicity", "numeric", "probability", "PROBABILITY", domain="toxicity", experimental_endpoint_aliases=("ames", "ames mutagenicity"), prediction_endpoint_aliases=("ames mutagenicity", "ames")),
    "DILI_LIABILITY": _ep("DILI_LIABILITY", "TOXICITY", "DILI clinical liability", "probability of drug-induced liver injury liability", "numeric", "probability", "PROBABILITY", domain="toxicity", prediction_endpoint_aliases=("dili clinical liability", "dili")),

    # -----------------------------------------------------------------
    # Group 9: Pharmacokinetics Disposition (Human & Preclinical)
    # -----------------------------------------------------------------
    "HUMAN_PK_CL_IV": _ep("HUMAN_PK_CL_IV", "PK", "Human IV systemic clearance", "human systemic total blood clearance (IV)", "numeric", "mL/min/kg", "LINEAR", domain="pharmacokinetics", species_requirement="HUMAN", route_requirement="IV", experimental_endpoint_aliases=("human cl", "human clearance", "clearance iv"), prediction_endpoint_aliases=("human_pk_cl_iv",)),
    "HUMAN_PK_VD_IV": _ep("HUMAN_PK_VD_IV", "PK", "Human IV volume of distribution", "human volume of distribution (IV)", "numeric", "L/kg", "LINEAR", domain="pharmacokinetics", species_requirement="HUMAN", route_requirement="IV", experimental_endpoint_aliases=("human vd", "volume of distribution iv")),
    "HUMAN_PK_F_IV": _ep("HUMAN_PK_F_IV", "PK", "Human IV bioavailability reference", "IV reference bioavailability (100%)", "numeric", "%", "PERCENT", domain="pharmacokinetics", species_requirement="HUMAN", route_requirement="IV"),
    "HUMAN_PK_T_HALF_IV": _ep("HUMAN_PK_T_HALF_IV", "PK", "Human IV terminal half-life", "elimination half-life following IV administration", "numeric", "hours", "LINEAR", domain="pharmacokinetics", species_requirement="HUMAN", route_requirement="IV"),

    "HUMAN_PK_CLF_ORAL": _ep("HUMAN_PK_CLF_ORAL", "PK", "Human oral clearance CL/F", "apparent oral clearance in humans", "numeric", "L/h", "LINEAR", domain="pharmacokinetics", species_requirement="HUMAN", route_requirement="ORAL", experimental_endpoint_aliases=("cl/f", "cl_f", "oral clearance", "human cl/f"), prediction_endpoint_aliases=("human_pk_clf_oral", "human_pk_clf_oral_oral")),
    "HUMAN_PK_VDF_ORAL": _ep("HUMAN_PK_VDF_ORAL", "PK", "Human oral Vd/F", "apparent oral volume of distribution in humans", "numeric", "L/kg", "LINEAR", domain="pharmacokinetics", species_requirement="HUMAN", route_requirement="ORAL", experimental_endpoint_aliases=("vd/f", "vd_f", "oral volume"), prediction_endpoint_aliases=("human_pk_vdf_oral", "human_pk_vdf_oral_oral")),
    "HUMAN_PK_VSSF_ORAL": _ep("HUMAN_PK_VSSF_ORAL", "PK", "Human oral Vss/F", "apparent volume of distribution at steady state in humans", "numeric", "L", "LINEAR", domain="pharmacokinetics", species_requirement="HUMAN", route_requirement="ORAL", experimental_endpoint_aliases=("vss/f", "vss_f")),
    "HUMAN_PK_F_ORAL": _ep("HUMAN_PK_F_ORAL", "PK", "Human oral bioavailability F", "fraction of oral dose reaching systemic circulation", "numeric", "%", "PERCENT", domain="pharmacokinetics", species_requirement="HUMAN", route_requirement="ORAL", experimental_endpoint_aliases=("bioavailability", "oral bioavailability", "f%"), prediction_endpoint_aliases=("human_pk_f_oral",)),

    "HUMAN_PK_CMAX_ORAL": _ep("HUMAN_PK_CMAX_ORAL", "PK", "Human oral Cmax", "maximum observed plasma concentration after oral dose", "numeric", "ng/mL", "LINEAR", domain="pharmacokinetics", species_requirement="HUMAN", route_requirement="ORAL", experimental_endpoint_aliases=("cmax", "c_max", "peak concentration"), prediction_endpoint_aliases=("human_pk_cmax_oral",)),
    "HUMAN_PK_TMAX_ORAL": _ep("HUMAN_PK_TMAX_ORAL", "PK", "Human oral Tmax", "time to maximum plasma concentration after oral dose", "numeric", "hours", "LINEAR", domain="pharmacokinetics", species_requirement="HUMAN", route_requirement="ORAL", experimental_endpoint_aliases=("tmax", "t_max"), prediction_endpoint_aliases=("human_pk_tmax_oral",)),
    "HUMAN_PK_AUC_ORAL": _ep("HUMAN_PK_AUC_ORAL", "PK", "Human oral AUC", "area under the plasma concentration-time curve after oral dose", "numeric", "ng*h/mL", "LINEAR", domain="pharmacokinetics", species_requirement="HUMAN", route_requirement="ORAL", experimental_endpoint_aliases=("auc", "auc0-t", "auc0-inf", "auctau", "auc_last"), prediction_endpoint_aliases=("human_pk_auc_oral",)),
    "HUMAN_PK_AUC0_INF_ORAL": _ep("HUMAN_PK_AUC0_INF_ORAL", "PK", "Human oral AUC0-inf", "AUC extrapolated to infinity after oral dose", "numeric", "ng*h/mL", "LINEAR", domain="pharmacokinetics", species_requirement="HUMAN", route_requirement="ORAL", experimental_endpoint_aliases=("auc0-inf", "aucinf"), prediction_endpoint_aliases=("human_pk_auc0_inf_oral",)),
    "HUMAN_PK_AUC0_T_ORAL": _ep("HUMAN_PK_AUC0_T_ORAL", "PK", "Human oral AUC0-t", "AUC from time zero to last quantifiable concentration", "numeric", "ng*h/mL", "LINEAR", domain="pharmacokinetics", species_requirement="HUMAN", route_requirement="ORAL", experimental_endpoint_aliases=("auc0-t", "auclast", "auc_last"), prediction_endpoint_aliases=("human_pk_auc0_t_oral",)),
    "HUMAN_PK_T_HALF_ORAL": _ep("HUMAN_PK_T_HALF_ORAL", "PK", "Human oral terminal half-life", "apparent elimination half-life following oral administration", "numeric", "hours", "LINEAR", domain="pharmacokinetics", species_requirement="HUMAN", route_requirement="ORAL", experimental_endpoint_aliases=("half-life", "t1/2", "t_half", "terminal half-life"), prediction_endpoint_aliases=("human_pk_t_half_oral",)),

    "HUMAN_PK_DDI_RELATIVE_RATIO": _ep("HUMAN_PK_DDI_RELATIVE_RATIO", "PK", "DDI relative exposure ratio", "AUC or Cmax geometric mean ratio in drug-drug interaction studies", "numeric", "ratio", "LINEAR", domain="pharmacokinetics", species_requirement="HUMAN", experimental_endpoint_aliases=("ddi ratio", "ratio cmax", "ratio auc")),

    # Preclinical Animal PK
    "RAT_PK_CL_IV": _ep("RAT_PK_CL_IV", "PK", "Rat IV clearance", "systemic clearance in rat", "numeric", "mL/min/kg", "LINEAR", domain="pharmacokinetics", species_requirement="RAT", route_requirement="IV"),
    "RAT_PK_VD_IV": _ep("RAT_PK_VD_IV", "PK", "Rat IV volume of distribution", "volume of distribution in rat", "numeric", "L/kg", "LINEAR", domain="pharmacokinetics", species_requirement="RAT", route_requirement="IV"),
    "RAT_PK_F_ORAL": _ep("RAT_PK_F_ORAL", "PK", "Rat oral bioavailability F", "oral bioavailability in rat", "numeric", "%", "PERCENT", domain="pharmacokinetics", species_requirement="RAT", route_requirement="ORAL"),
    "DOG_PK_CL_IV": _ep("DOG_PK_CL_IV", "PK", "Dog IV clearance", "systemic clearance in dog", "numeric", "mL/min/kg", "LINEAR", domain="pharmacokinetics", species_requirement="DOG", route_requirement="IV"),
    "DOG_PK_VD_IV": _ep("DOG_PK_VD_IV", "PK", "Dog IV volume of distribution", "volume of distribution in dog", "numeric", "L/kg", "LINEAR", domain="pharmacokinetics", species_requirement="DOG", route_requirement="IV"),
    "DOG_PK_F_ORAL": _ep("DOG_PK_F_ORAL", "PK", "Dog oral bioavailability F", "oral bioavailability in dog", "numeric", "%", "PERCENT", domain="pharmacokinetics", species_requirement="DOG", route_requirement="ORAL"),
    "MONKEY_PK_CL_IV": _ep("MONKEY_PK_CL_IV", "PK", "Monkey IV clearance", "systemic clearance in monkey", "numeric", "mL/min/kg", "LINEAR", domain="pharmacokinetics", species_requirement="MONKEY", route_requirement="IV"),
    "MONKEY_PK_F_ORAL": _ep("MONKEY_PK_F_ORAL", "PK", "Monkey oral bioavailability F", "oral bioavailability in monkey", "numeric", "%", "PERCENT", domain="pharmacokinetics", species_requirement="MONKEY", route_requirement="ORAL"),

    # -----------------------------------------------------------------
    # Group 10: Activity & Literature (6 endpoints)
    # -----------------------------------------------------------------
    "ACTIVITY_IC50": _ep("ACTIVITY_IC50", "ACTIVITY", "IC50", "half-maximal inhibitory concentration against target", "numeric", "nM", "LINEAR", domain="pharmacodynamics", experimental_endpoint_aliases=("ic50", "target ic50")),
    "ACTIVITY_EC50": _ep("ACTIVITY_EC50", "ACTIVITY", "EC50", "half-maximal effective concentration against target", "numeric", "nM", "LINEAR", domain="pharmacodynamics", experimental_endpoint_aliases=("ec50", "target ec50")),
    "ACTIVITY_KI": _ep("ACTIVITY_KI", "ACTIVITY", "Ki", "inhibition constant against target", "numeric", "nM", "LINEAR", domain="pharmacodynamics", experimental_endpoint_aliases=("ki", "target ki")),
    "ACTIVITY_KD": _ep("ACTIVITY_KD", "ACTIVITY", "Kd", "dissociation constant against target", "numeric", "nM", "LINEAR", domain="pharmacodynamics", experimental_endpoint_aliases=("kd", "target kd")),
    "ACTIVITY_RATIO": _ep("ACTIVITY_RATIO", "ACTIVITY", "Selectivity ratio", "potency or selectivity ratio between targets/variants", "numeric", "ratio", "LINEAR", domain="pharmacodynamics", experimental_endpoint_aliases=("ratio ic50", "selectivity ratio", "fold selectivity")),
    "LITERATURE_CITATION": _ep("LITERATURE_CITATION", "UNCLASSIFIED", "Literature reference", "bibliographic reference or publication mention", "qualitative", "", "MIXED", domain="literature", experimental_endpoint_aliases=("literature candidate", "citation", "pubmed")),
}

_SPECIES_ALIASES = {
    "human": "HUMAN", "homo sapiens": "HUMAN", "patient": "HUMAN", "patients": "HUMAN",
    "healthy volunteer": "HUMAN", "healthy volunteers": "HUMAN", "clinical": "HUMAN",
    "rat": "RAT", "sd rat": "RAT", "sprague-dawley rat": "RAT", "sprague dawley rat": "RAT", "rattus norvegicus": "RAT",
    "mouse": "MOUSE", "mus musculus": "MOUSE", "mice": "MOUSE",
    "dog": "DOG", "beagle": "DOG", "canine": "DOG",
    "monkey": "MONKEY", "cynomolgus": "MONKEY", "cynomolgus monkey": "MONKEY", "nonhuman primate": "MONKEY", "rhesus": "MONKEY",
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
    if re.search(r"\b(ip|intraperitoneal)\b", text): return "IP"
    return "UNSPECIFIED"


def _number(value: Any) -> float | None:
    if value is None: return None
    match = re.search(r"[-+]?\d+(?:,\d{3})*(?:\.\d+)?(?:[eE][-+]?\d+)?", str(value))
    if not match: return None
    try: return float(match.group(0).replace(",", ""))
    except ValueError: return None


def _pk_key(raw: str, species: str, route: str, context: str) -> tuple[str, str]:
    text = raw.lower().replace("–", "-")
    if "cmax" in text: parameter = "CMAX"
    elif "tmax" in text: parameter = "TMAX"
    elif re.search(r"auc", text):
        if re.search(r"auc\s*(?:0\s*[- ]\s*inf|inf)", text + " " + context): parameter = "AUC0_INF"
        elif "tau" in text + context: parameter = "AUC_TAU"
        elif re.search(r"auc\s*(?:0\s*[- ]\s*t|last|tlast)", text + " " + context): parameter = "AUC0_T"
        else: parameter = "AUC"
    elif "half" in text or "t1/2" in text: parameter = "T_HALF"
    elif "cl/f" in text or ("apparent" in text + " " + context and "clearance" in text + " " + context): parameter = "CLF_ORAL"
    elif "clearance" in text or re.search(r"\bcl\b", text): parameter = "CL"
    elif ("apparent" in text + " " + context and "vss" in text + " " + context): parameter = "VSSF_ORAL"
    elif "vd/f" in text or ("apparent" in text + " " + context and "volume" in text + " " + context): parameter = "VDF_ORAL"
    elif "vss" in text or "volume of distribution at steady state" in text: parameter = "VSS"
    elif "volume" in text or re.search(r"\bvd\b", text): parameter = "VD"
    elif "bioavailability" in text or re.search(r"\bf\b", text): parameter = "F"
    else: parameter = "UNSPECIFIED"

    if parameter in {"CLF_ORAL", "VDF_ORAL", "VSSF_ORAL"}:
        route = "ORAL"
    suffix = route if route != "UNSPECIFIED" else "ORAL"
    return f"{species}_PK_{parameter}_{suffix}", parameter


def normalize_experimental_observation(
    raw_endpoint: Any,
    raw_value: Any = None,
    raw_unit: Any = "",
    *,
    species: Any = "",
    context: Any = "",
    assay_type: Any = "",
    target: Any = "",
    canonical_hint: Any = "",
    mw: Optional[float] = None,
) -> dict:
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

    raw_l = raw.lower()

    # 1. Literature Citation / Bibliographic records
    if raw_l in {"literature candidate", "citation", "pubmed"} or "literature" in raw_l:
        return {
            "canonical_endpoint_id": "LITERATURE_CITATION",
            "section": "UNCLASSIFIED",
            "display_name": "Literature Citation",
            "species": normalized_species,
            "route": route,
            "measurement_subtype": "CITATION",
            "normalized_value": None,
            "normalized_unit": "",
            "comparability_status": UNSUPPORTED,
            "normalization_rule": "non_numeric_citation",
            "reason": "Literature citation record without discrete scalar measurement",
            "comparison_key": "LITERATURE_CITATION",
        }

    # 2. Activity Measurements (IC50, EC50, Ki, Kd, Ratio)
    if re.search(r"\b(ic50|ec50|ki|kd)\b", raw_l) and not re.search(r"(?:cyp\s*[0-9a-z]+|herg|p-gp|pgp)", f"{raw} {target}".lower()):
        if "ratio" in raw_l or "shift" in raw_l:
            return {
                "canonical_endpoint_id": "ACTIVITY_RATIO",
                "section": "ACTIVITY",
                "display_name": "Selectivity Ratio",
                "species": normalized_species,
                "route": route,
                "measurement_subtype": "Ratio",
                "normalized_value": number,
                "normalized_unit": raw_unit or "ratio",
                "comparability_status": RELATED,
                "normalization_rule": "ratio_identity",
                "reason": "Ratio/fold-shift is related, not direct concentration potency",
                "comparison_key": f"ACTIVITY_RATIO|{target or 'UNSPECIFIED'}|{assay_type or 'UNSPECIFIED'}"
            }
        subtype = re.search(r"ic50|ec50|ki|kd", raw_l).group(0).upper()
        endpoint = f"ACTIVITY_{subtype}"
        # Convert activity concentration to nM if in µM or M
        norm_val = number
        norm_u = raw_unit or "nM"
        if number is not None and raw_unit:
            try:
                c_res = convert_concentration(number, str(raw_unit), "nM", mw=mw)
                norm_val = c_res.normalized_value
                norm_u = "nM"
            except Exception:
                pass
        return {
            "canonical_endpoint_id": endpoint,
            "section": "ACTIVITY",
            "display_name": subtype,
            "species": normalized_species,
            "route": route,
            "measurement_subtype": subtype,
            "normalized_value": norm_val,
            "normalized_unit": norm_u,
            "comparability_status": RELATED,
            "normalization_rule": "activity_semantic_group",
            "reason": "Project assay/target mapping is required for a direct activity comparison",
            "comparison_key": f"{endpoint}|{target or 'UNSPECIFIED'}|{assay_type or 'UNSPECIFIED'}"
        }

    # 3. Plasma Protein Binding (PPB) & fu
    if (re.search(r"protein binding|plasma protein|\bppb\b|fraction unbound|\bfu\b", raw_l) or raw_l in {"protein", "bound fraction"}) and not re.search(r"\b(ic50|ec50|ki|kd)\b", raw_l):
        endpoint = {"HUMAN": "HUMAN_PPB", "RAT": "RAT_PPB", "MOUSE": "MOUSE_PPB"}.get(normalized_species, "HUMAN_PPB" if normalized_species == "UNSPECIFIED" and "human" in all_text else "PPB_UNSPECIFIED")
        if endpoint == "PPB_UNSPECIFIED":
            return {
                "canonical_endpoint_id": endpoint,
                "section": "ADMET",
                "display_name": "Plasma protein binding",
                "species": normalized_species,
                "route": route,
                "measurement_subtype": "PPB",
                "normalized_value": None,
                "normalized_unit": "% bound",
                "comparability_status": CONDITIONAL,
                "normalization_rule": "species_required",
                "reason": "Species is required for PPB comparison",
                "comparison_key": f"{endpoint}|{normalized_species}"
            }
        if number is not None:
            conv = convert_ppb(number, str(raw_unit), "% bound")
            status = DIRECT if conv.conversion_type == TYPE_DIRECT_IDENTITY else CONVERTED
            return {
                "canonical_endpoint_id": endpoint,
                "section": "ADMET",
                "display_name": REGISTRY[endpoint].display_name,
                "species": normalized_species,
                "route": route,
                "measurement_subtype": "PPB",
                "normalized_value": conv.normalized_value,
                "normalized_unit": conv.normalized_unit,
                "comparability_status": status,
                "normalization_rule": conv.conversion_formula,
                "reason": "",
                "comparison_key": f"{endpoint}|{normalized_species}"
            }
        return {
            "canonical_endpoint_id": endpoint,
            "section": "ADMET",
            "display_name": REGISTRY[endpoint].display_name,
            "species": normalized_species,
            "route": route,
            "measurement_subtype": "PPB",
            "normalized_value": None,
            "normalized_unit": "% bound",
            "comparability_status": UNSUPPORTED,
            "normalization_rule": "",
            "reason": "PPB value/unit is not safely numeric",
            "comparison_key": f"{endpoint}|{normalized_species}"
        }

    # 4. Caco-2 Permeability
    if re.search(r"caco[- ]?2|caco2|papp", raw_l) or (raw_l in {"permeability", "permeability assay"} and re.search(r"caco[- ]?2|caco2", all_text)):
        if "efflux" in all_text:
            endpoint = "CACO2_EFFLUX_RATIO"
            return {
                "canonical_endpoint_id": endpoint,
                "section": "ADMET",
                "display_name": REGISTRY[endpoint].display_name,
                "species": normalized_species,
                "route": route,
                "measurement_subtype": "EFFLUX_RATIO",
                "normalized_value": number,
                "normalized_unit": "ratio",
                "comparability_status": DIRECT if number is not None else UNSUPPORTED,
                "normalization_rule": "identity",
                "reason": "",
                "comparison_key": endpoint
            }
        elif re.search(r"b\s*[- >]+\s*a|bto a|basolateral", all_text):
            endpoint = "CACO2_PAPP_BA"
        else:
            endpoint = "CACO2_PAPP_AB"

        if endpoint == "CACO2_PAPP_AB" and not re.search(r"a\s*[- >]+\s*b|ato b|apical.to.basolateral", all_text):
            status, reason = CONDITIONAL, "Caco-2 direction is not recorded"
            return {
                "canonical_endpoint_id": endpoint,
                "section": "ADMET",
                "display_name": REGISTRY[endpoint].display_name,
                "species": normalized_species,
                "route": route,
                "measurement_subtype": "PAPP",
                "normalized_value": None,
                "normalized_unit": "log10(cm/s)",
                "comparability_status": status,
                "normalization_rule": "",
                "reason": reason,
                "comparison_key": endpoint
            }
        if number is not None:
            conv = convert_caco2_papp(number, str(raw_unit))
            status = DIRECT if conv.conversion_type == TYPE_DIRECT_IDENTITY else CONVERTED
            return {
                "canonical_endpoint_id": endpoint,
                "section": "ADMET",
                "display_name": REGISTRY[endpoint].display_name,
                "species": normalized_species,
                "route": route,
                "measurement_subtype": "PAPP",
                "normalized_value": conv.normalized_value,
                "normalized_unit": conv.normalized_unit,
                "comparability_status": status,
                "normalization_rule": conv.conversion_formula,
                "reason": "",
                "comparison_key": endpoint
            }
        return {
            "canonical_endpoint_id": endpoint,
            "section": "ADMET",
            "display_name": REGISTRY[endpoint].display_name,
            "species": normalized_species,
            "route": route,
            "measurement_subtype": "PAPP",
            "normalized_value": None,
            "normalized_unit": "log10(cm/s)",
            "comparability_status": UNSUPPORTED,
            "normalization_rule": "",
            "reason": "Caco-2 Papp unit is not safely supported",
            "comparison_key": endpoint
        }

    # 5. Aqueous Solubility
    if "solubility" in raw_l or (raw_l in {"log s", "logs"} and "aqueous" in all_text):
        endpoint = "SOLUBILITY_INTRINSIC" if "intrinsic" in all_text else ("SOLUBILITY_KINETIC" if "kinetic" in all_text else ("SOLUBILITY_THERMODYNAMIC" if "thermodynamic" in all_text else "SOLUBILITY_GENERIC"))
        if number is not None:
            try:
                conv = convert_solubility(number, str(raw_unit), mw=mw)
                status = DIRECT if conv.conversion_type == TYPE_DIRECT_IDENTITY else CONVERTED
                return {
                    "canonical_endpoint_id": endpoint,
                    "section": "ADMET",
                    "display_name": REGISTRY[endpoint].display_name,
                    "species": normalized_species,
                    "route": route,
                    "measurement_subtype": endpoint,
                    "normalized_value": conv.normalized_value,
                    "normalized_unit": conv.normalized_unit,
                    "comparability_status": status,
                    "normalization_rule": conv.conversion_formula,
                    "reason": "",
                    "comparison_key": endpoint
                }
            except Exception as e:
                return {
                    "canonical_endpoint_id": endpoint,
                    "section": "ADMET",
                    "display_name": REGISTRY[endpoint].display_name,
                    "species": normalized_species,
                    "route": route,
                    "measurement_subtype": endpoint,
                    "normalized_value": None,
                    "normalized_unit": "log10(mol/L)",
                    "comparability_status": UNSUPPORTED,
                    "normalization_rule": "",
                    "reason": str(e),
                    "comparison_key": endpoint
                }

    # 6. Microsomal & Hepatocyte Clearance
    if re.search(r"hepatocyte", all_text) or raw_l == "clh":
        endpoint = "HEPATOCYTE_CLINT"
        return {
            "canonical_endpoint_id": endpoint,
            "section": "METABOLISM",
            "display_name": REGISTRY[endpoint].display_name,
            "species": normalized_species,
            "route": route,
            "measurement_subtype": "CLINT",
            "normalized_value": number,
            "normalized_unit": str(raw_unit or "µL/min/10^6 cells"),
            "comparability_status": RELATED,
            "normalization_rule": "identity",
            "reason": "Hepatocyte Clint is cellular clearance; distinct from microsomal HLM Clint",
            "comparison_key": f"{endpoint}|{normalized_species}"
        }

    if re.search(r"\b(hlm|rlm|mlm)\b|microsom", all_text):
        endpoint = "HLM_CLINT" if "human" in all_text or "hlm" in all_text else ("RLM_CLINT" if "rat" in all_text or "rlm" in all_text else ("MLM_CLINT" if "mouse" in all_text or "mlm" in all_text else "HLM_CLINT"))
        if number is not None:
            conv = convert_clearance(number, str(raw_unit), "log10(mL/min/kg)") if "log" not in str(raw_unit).lower() else None
            norm_val = conv.normalized_value if conv else number
            norm_unit = "log10(mL/min/kg)"
            rule = conv.conversion_formula if conv else "identity"
            status = CONVERTED if conv and conv.conversion_type != TYPE_DIRECT_IDENTITY else DIRECT
            return {
                "canonical_endpoint_id": endpoint,
                "section": "ADMET",
                "display_name": REGISTRY[endpoint].display_name,
                "species": normalized_species,
                "route": route,
                "measurement_subtype": "CLINT",
                "normalized_value": norm_val,
                "normalized_unit": norm_unit,
                "comparability_status": status,
                "normalization_rule": rule,
                "reason": "",
                "comparison_key": endpoint
            }

    # 7. Excretion & Mass Balance
    if raw_l in {"feces", "fecal excretion"}:
        endpoint = "EXCRETION_FECAL"
        return {
            "canonical_endpoint_id": endpoint,
            "section": "METABOLISM",
            "display_name": REGISTRY[endpoint].display_name,
            "species": normalized_species,
            "route": route,
            "measurement_subtype": "EXCRETION",
            "normalized_value": number,
            "normalized_unit": "% dose",
            "comparability_status": DIRECT if number is not None else UNSUPPORTED,
            "normalization_rule": "identity",
            "reason": "",
            "comparison_key": f"{endpoint}|{normalized_species}"
        }
    if raw_l in {"urine", "urinary excretion"}:
        endpoint = "EXCRETION_URINARY"
        return {
            "canonical_endpoint_id": endpoint,
            "section": "METABOLISM",
            "display_name": REGISTRY[endpoint].display_name,
            "species": normalized_species,
            "route": route,
            "measurement_subtype": "EXCRETION",
            "normalized_value": number,
            "normalized_unit": "% dose",
            "comparability_status": DIRECT if number is not None else UNSUPPORTED,
            "normalization_rule": "identity",
            "reason": "",
            "comparison_key": f"{endpoint}|{normalized_species}"
        }
    if raw_l == "metabolite":
        endpoint = "METABOLITE_OBSERVATION"
        return {
            "canonical_endpoint_id": endpoint,
            "section": "METABOLISM",
            "display_name": REGISTRY[endpoint].display_name,
            "species": normalized_species,
            "route": route,
            "measurement_subtype": "METABOLITE",
            "normalized_value": number,
            "normalized_unit": str(raw_unit or ""),
            "comparability_status": DIRECT if number is not None else UNSUPPORTED,
            "normalization_rule": "identity",
            "reason": "" if number is not None else "Non-numeric metabolite observation",
            "comparison_key": f"{endpoint}|{normalized_species}"
        }

    # 8. CYP Inhibition (Quantitative IC50 -> pIC50)
    cyp_match = re.search(r"\b(cyp\s*(?:1a2|2c9|2c19|2d6|3a4?))\b", all_text)
    if cyp_match:
        iso_raw = cyp_match.group(1).replace(" ", "").upper()
        iso = "CYP3A4" if iso_raw in {"CYP3A", "CYP3A4"} else iso_raw
        if re.search(r"inhib|ic50|ki|potency", all_text):
            endpoint = f"{iso}_INHIBITION"
            return {
                "canonical_endpoint_id": endpoint,
                "section": "METABOLISM",
                "display_name": REGISTRY.get(endpoint, _ep(endpoint, "METABOLISM", endpoint, "", "numeric", "pIC50", "PIC50")).display_name,
                "species": normalized_species,
                "route": route,
                "measurement_subtype": "INHIBITION",
                "normalized_value": number,
                "normalized_unit": str(raw_unit or ""),
                "comparability_status": RELATED,
                "normalization_rule": "identity",
                "reason": "Related scientific evidence; continuous potency vs classification model",
                "comparison_key": f"{endpoint}|{normalized_species}"
            }
        elif re.search(r"fm|metabol|substrate", all_text):
            endpoint = f"{iso}_METABOLIC_CONTRIBUTION"
            return {
                "canonical_endpoint_id": endpoint,
                "section": "METABOLISM",
                "display_name": f"{iso} metabolic contribution",
                "species": normalized_species,
                "route": route,
                "measurement_subtype": "METABOLIC_CONTRIBUTION",
                "normalized_value": number,
                "normalized_unit": "%" if "%" in str(raw_unit) else str(raw_unit or ""),
                "comparability_status": DIRECT if number is not None else UNSUPPORTED,
                "normalization_rule": "identity",
                "reason": "",
                "comparison_key": f"{endpoint}|{normalized_species}"
            }

    # 9. hERG Liability (IC50 -> pIC50)
    if "herg" in raw_l:
        endpoint = "HERG_LIABILITY"
        if number is not None and any(tok in clean_unit_str(raw_unit) for tok in ("m", "µm", "um", "nm", "pm", "mol")):
            conv = convert_cyp_herg_ic50(number, str(raw_unit))
            return {
                "canonical_endpoint_id": endpoint,
                "section": "TOXICITY",
                "display_name": REGISTRY[endpoint].display_name,
                "species": normalized_species,
                "route": route,
                "measurement_subtype": "HERG_POTENCY",
                "normalized_value": conv.normalized_value,
                "normalized_unit": "pIC50",
                "comparability_status": CONVERTED,
                "normalization_rule": conv.conversion_formula,
                "reason": f"Converted from {number} {raw_unit}",
                "comparison_key": endpoint
            }
        return {
            "canonical_endpoint_id": endpoint,
            "section": "TOXICITY",
            "display_name": REGISTRY[endpoint].display_name,
            "species": normalized_species,
            "route": route,
            "measurement_subtype": "HERG_RISK",
            "normalized_value": number,
            "normalized_unit": str(raw_unit or ""),
            "comparability_status": RELATED,
            "normalization_rule": "identity",
            "reason": "Qualitative hERG channel interaction",
            "comparison_key": endpoint
        }

    # 10. Transporters (P-gp, BCRP)
    if raw_l in {"p-gp", "pgp", "p-gp inhibition"} or "p-gp" in raw_l:
        endpoint = "PGP_INHIBITION"
        return {
            "canonical_endpoint_id": endpoint,
            "section": "METABOLISM",
            "display_name": REGISTRY[endpoint].display_name,
            "species": normalized_species,
            "route": route,
            "measurement_subtype": "PGP",
            "normalized_value": number,
            "normalized_unit": str(raw_unit or ""),
            "comparability_status": RELATED,
            "normalization_rule": "identity",
            "reason": "P-gp interaction observation; quantitative continuous kinetics remain fail-closed",
            "comparison_key": f"{endpoint}|{normalized_species}"
        }
    if raw_l in {"bcrp", "bcrp inhibition"} or "bcrp" in raw_l:
        endpoint = "BCRP_INHIBITION"
        return {
            "canonical_endpoint_id": endpoint,
            "section": "METABOLISM",
            "display_name": REGISTRY[endpoint].display_name,
            "species": normalized_species,
            "route": route,
            "measurement_subtype": "BCRP",
            "normalized_value": number,
            "normalized_unit": str(raw_unit or ""),
            "comparability_status": RELATED,
            "normalization_rule": "identity",
            "reason": "BCRP interaction observation; quantitative continuous kinetics remain fail-closed",
            "comparison_key": f"{endpoint}|{normalized_species}"
        }

    # 11. In Vivo Pharmacokinetics (Cmax, AUC, Tmax, t1/2, CL, CL/F, Vd, Vd/F, Vss, F)
    if re.search(r"\b(?:cmax|tmax|auc[0-9a-z_-]*|half[- ]?life|t1/2|clearance|cl$|cl/f|volume of distribution|vd|vd/f|vss|vss/f|bioavailability|f)\b", raw_l):
        key, parameter = _pk_key(raw, normalized_species, route, context_s)
        norm_val = number
        norm_u = str(raw_unit or "")
        rule = "identity"
        status = DIRECT if number is not None else UNSUPPORTED

        if number is not None:
            try:
                if parameter == "CMAX":
                    c_res = convert_concentration(number, str(raw_unit), "ng/mL", mw=mw)
                    norm_val = c_res.normalized_value
                    norm_u = "ng/mL"
                    rule = c_res.conversion_formula
                    status = DIRECT if c_res.conversion_type == TYPE_DIRECT_IDENTITY else CONVERTED
                elif parameter.startswith("AUC"):
                    a_res = convert_auc(number, str(raw_unit), "ng*h/mL")
                    norm_val = a_res.normalized_value
                    norm_u = "ng*h/mL"
                    rule = a_res.conversion_formula
                    status = DIRECT if a_res.conversion_type == TYPE_DIRECT_IDENTITY else CONVERTED
                elif parameter in {"TMAX", "T_HALF"}:
                    t_res = convert_time(number, str(raw_unit), "hours")
                    norm_val = t_res.normalized_value
                    norm_u = "hours"
                    rule = t_res.conversion_formula
                    status = DIRECT if t_res.conversion_type == TYPE_DIRECT_IDENTITY else CONVERTED
                elif parameter in {"CL", "CLF_ORAL"}:
                    preferred_u = "L/h" if normalized_species == "HUMAN" and "l/h" in clean_unit_str(raw_unit) else "mL/min/kg"
                    cl_res = convert_clearance(number, str(raw_unit), preferred_u)
                    norm_val = cl_res.normalized_value
                    norm_u = preferred_u
                    rule = cl_res.conversion_formula
                    status = DIRECT if cl_res.conversion_type == TYPE_DIRECT_IDENTITY else CONVERTED
                elif parameter in {"VD", "VSS", "VDF_ORAL", "VSSF_ORAL"}:
                    preferred_u = "L" if "l" == clean_unit_str(raw_unit) else "L/kg"
                    v_res = convert_volume(number, str(raw_unit), preferred_u)
                    norm_val = v_res.normalized_value
                    norm_u = preferred_u
                    rule = v_res.conversion_formula
                    status = DIRECT if v_res.conversion_type == TYPE_DIRECT_IDENTITY else CONVERTED
                elif parameter == "F":
                    if "%" in clean_unit_str(raw_unit):
                        norm_val = number
                        norm_u = "%"
                        status = DIRECT
                    elif 0.0 <= number <= 1.0:
                        norm_val = number * 100.0
                        norm_u = "%"
                        rule = f"{number} * 100"
                        status = CONVERTED
            except Exception:
                pass

        return {
            "canonical_endpoint_id": key,
            "section": "PK",
            "display_name": f"{normalized_species.title()} {parameter.replace('_', ' ')}" if normalized_species != "UNSPECIFIED" else parameter.replace('_', ' '),
            "species": normalized_species,
            "route": route,
            "analyte": analyte,
            "measurement_subtype": parameter,
            "normalized_value": norm_val,
            "normalized_unit": norm_u,
            "comparability_status": status,
            "normalization_rule": rule,
            "reason": "" if norm_val is not None else "PK numeric measurement not parsed",
            "comparison_key": f"{key}|{normalized_species}|{route}|{analyte}"
        }

    # 12. Physicochemical Descriptors (pKa, logD, logP)
    if "pka" in raw_l:
        if any(tok in raw_l for tok in ("acid", "acidic")):
            eid = "PKA_ACID"
            dname = "Acidic pKa"
        elif any(tok in raw_l for tok in ("base", "basic")):
            eid = "PKA_BASE"
            dname = "Basic pKa"
        else:
            eid = "PKA"
            dname = "pKa"
        return {
            "canonical_endpoint_id": eid,
            "section": "ADMET",
            "display_name": dname,
            "species": normalized_species,
            "route": route,
            "measurement_subtype": eid,
            "normalized_value": number,
            "normalized_unit": "pKa",
            "comparability_status": DIRECT if number is not None else UNSUPPORTED,
            "normalization_rule": "identity",
            "reason": "",
            "comparison_key": eid
        }
    if "logd" in raw_l:
        return {
            "canonical_endpoint_id": "LOGD_7_4",
            "section": "ADMET",
            "display_name": "logD 7.4",
            "species": normalized_species,
            "route": route,
            "measurement_subtype": "LOGD",
            "normalized_value": number,
            "normalized_unit": "logD",
            "comparability_status": DIRECT if number is not None else UNSUPPORTED,
            "normalization_rule": "identity",
            "reason": "",
            "comparison_key": "LOGD_7_4"
        }
    if "logp" in raw_l:
        return {
            "canonical_endpoint_id": "LOGP_RELATED",
            "section": "ADMET",
            "display_name": "logP (related)",
            "species": normalized_species,
            "route": route,
            "measurement_subtype": "LOGP",
            "normalized_value": number,
            "normalized_unit": "logP",
            "comparability_status": RELATED,
            "normalization_rule": "identity",
            "reason": "logP is related; not identical to pH 7.4 distribution coefficient logD",
            "comparison_key": "LOGP_RELATED"
        }

    # 13. Canonical Hint fallback
    hint = str(canonical_hint or "").upper()
    if hint in REGISTRY:
        ep = REGISTRY[hint]
        return {
            "canonical_endpoint_id": hint,
            "section": ep.section,
            "display_name": ep.display_name,
            "species": normalized_species,
            "route": route,
            "measurement_subtype": hint,
            "normalized_value": number,
            "normalized_unit": ep.canonical_unit,
            "comparability_status": DIRECT if number is not None else UNSUPPORTED,
            "normalization_rule": "canonical_hint",
            "reason": "" if number is not None else "Non-numeric source value",
            "comparison_key": hint
        }

    return {
        "canonical_endpoint_id": "UNRESOLVED",
        "section": "UNCLASSIFIED",
        "display_name": raw or "Unclassified evidence",
        "species": normalized_species,
        "route": route,
        "measurement_subtype": raw,
        "normalized_value": number,
        "normalized_unit": str(raw_unit or ""),
        "comparability_status": UNSUPPORTED,
        "normalization_rule": "",
        "reason": "No canonical endpoint mapping",
        "comparison_key": f"UNRESOLVED|{raw}"
    }


def canonicalize_prediction_endpoint(endpoint: Any, *, species: Any = "", route: Any = "", context: Any = "") -> dict:
    raw = str(endpoint or "").strip()
    text = raw.lower()
    sp = normalize_species(species, context)
    rt = normalize_route(route or context)

    # First check REGISTRY by exact ID
    upper_raw = raw.upper().replace(" ", "_")
    if upper_raw in REGISTRY:
        ep = REGISTRY[upper_raw]
        return {
            "canonical_endpoint_id": upper_raw,
            "species": ep.species_requirement or sp,
            "route": ep.route_requirement or rt,
            "parameter": "",
            "raw_endpoint": raw,
            "comparison_key": f"{upper_raw}|{ep.species_requirement or sp}|{ep.route_requirement or rt}"
        }

    # Check aliases
    for eid, ep in REGISTRY.items():
        if text == ep.display_name.lower() or any(text == alias.lower() for alias in ep.prediction_endpoint_aliases):
            return {
                "canonical_endpoint_id": eid,
                "species": ep.species_requirement or sp,
                "route": ep.route_requirement or rt,
                "parameter": "",
                "raw_endpoint": raw,
                "comparison_key": f"{eid}|{ep.species_requirement or sp}|{ep.route_requirement or rt}"
            }

    # CYP pattern
    isoform = re.search(r"cyp\s*(1a2|2c9|2c19|2d6|3a4)", text)
    if isoform:
        iso = isoform.group(1).upper()
        eid = f"CYP{iso}_INHIBITION" if re.search(r"inhib|block", text) else (f"CYP{iso}_SUBSTRATE" if "substr" in text else f"CYP{iso}_METABOLIC_CONTRIBUTION")
        return {"canonical_endpoint_id": eid, "species": sp, "route": rt, "parameter": "", "raw_endpoint": raw, "comparison_key": f"{eid}|{sp}|{rt}"}

    # PK pattern
    if re.match(r"^(cmax|tmax|auc|half|terminal|clearance|cl$|cl/f|vd|vd/f|vss|volume|apparent|bioavailability|f$)", text):
        eid, parameter = _pk_key(raw, sp, rt, _context_text(context))
        return {"canonical_endpoint_id": eid, "species": sp, "route": rt, "parameter": parameter, "raw_endpoint": raw, "comparison_key": f"{eid}|{sp}|{rt}"}

    return {"canonical_endpoint_id": upper_raw, "species": sp, "route": rt, "parameter": "", "raw_endpoint": raw, "comparison_key": f"{upper_raw}|{sp}|{rt}"}


def prediction_source_type(*, source: Any = "", prediction_type: Any = "", endpoint: Any = "", default: str = PREDICTION_MODEL) -> str:
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
    return {
        "canonical_endpoint_version": CANONICAL_ENDPOINT_VERSION,
        "comparison_unit_version": COMPARISON_UNIT_VERSION,
        "endpoint_count": len(REGISTRY),
        "endpoints": [asdict(item) for item in REGISTRY.values()]
    }


def reindex_persisted_evidence(db, version_id: int | None = None) -> dict:
    """Reclassify existing raw evidence without another source search.

    Only derived routing/normalization columns are changed. Raw endpoint,
    value, unit, provenance and reference fields remain untouched.
    """
    from sqlalchemy import select
    from .models import ExternalExperimentalEvidence, CompoundVersion

    query = select(ExternalExperimentalEvidence)
    if version_id is not None:
        query = query.where(ExternalExperimentalEvidence.compound_version_id == version_id)
    rows = list(db.scalars(query).all())

    # Cache compound molecular weight for accurate molar/mass conversions
    mw_cache: Dict[int, float] = {}
    changed = 0
    for row in rows:
        vid = row.compound_version_id
        if vid not in mw_cache:
            cver = db.get(CompoundVersion, vid)
            if cver and cver.canonical_smiles:
                try:
                    from rdkit import Chem
                    from rdkit.Chem import Descriptors
                    mol = Chem.MolFromSmiles(cver.canonical_smiles)
                    mw_cache[vid] = round(float(Descriptors.MolWt(mol)), 2) if mol else 0.0
                except Exception:
                    mw_cache[vid] = 0.0
            else:
                mw_cache[vid] = 0.0

        mw = mw_cache.get(vid)
        context = row.assay_conditions_json if isinstance(row.assay_conditions_json, dict) else {"conditions": row.assay_conditions_json or ""}
        mapped = normalize_experimental_observation(
            row.raw_endpoint_name,
            row.raw_value,
            row.raw_unit,
            species=row.species,
            context=context,
            assay_type=row.assay_type,
            target=context.get("target", ""),
            canonical_hint=row.canonical_endpoint_id,
            mw=mw if mw and mw > 0 else None,
        )
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
            for key, value in derived.items():
                setattr(row, key, value)
            changed += 1

    return {
        "examined": len(rows),
        "changed": changed,
        "canonical_endpoint_version": CANONICAL_ENDPOINT_VERSION,
        "comparison_unit_version": COMPARISON_UNIT_VERSION,
    }
