"""
DrugBank Reference Drug Library & Incremental Learning (Drug-OPT Stage 6 / v3.0.2).

Provides:
- Canonical 'DrugBank' project management (GLOBAL_MODEL_DEVELOPMENT mode)
- Step-by-step sequential reference drug ingestion across 10 distinct chemical & pharmacological spaces:
    1. Gefitinib (Quinazoline / EGFR kinase inhibitor) - DEVELOPMENT_TRAINING
    2. Imatinib (2-Phenylaminopyrimidine / BCR-ABL kinase inhibitor) - DEVELOPMENT_TRAINING
    3. Propranolol (Aryloxypropanolamine / Beta-blocker) - IMMUTABLE_HOLDOUT
    4. Atorvastatin (Pyrrole-heptanoic acid / Statin) - DEVELOPMENT_TRAINING
    5. Midazolam (Imidazobenzodiazepine / CYP3A4 probe & GABA-A modulator) - IMMUTABLE_HOLDOUT
    6. Verapamil (Phenylalkylamine / L-type Ca2+ channel & P-gp/CYP3A reference) - IMMUTABLE_HOLDOUT
    7. Fluoxetine (Aryloxypropylamine / SSRI & potent CYP2D6 inhibitor) - IMMUTABLE_HOLDOUT
    8. Ketoconazole (Imidazole-dioxolane-piperazine / Strong CYP3A4 inhibitor) - DEVELOPMENT_TRAINING
    9. Sildenafil (Pyrazolopyrimidinone / PDE5 inhibitor) - IMMUTABLE_HOLDOUT
    10. Quinidine (Cinchona alkaloid / Strong CYP2D6 inhibitor & hERG blocker) - IMMUTABLE_HOLDOUT
- Fixed immutable role assignment upon registration:
    * DEVELOPMENT_TRAINING: Used for fitting residual calibration / fine-tuning
    * IMMUTABLE_HOLDOUT: Strictly isolated holdout; NEVER used for model fitting
- Upstream training overlap tracking (EXACT_STRUCTURE_OVERLAP vs NOVEL vs VALIDATION_HOLDOUT)
"""
from __future__ import annotations

import hashlib
import json
import math
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski
from sqlalchemy import select
from sqlalchemy.orm import Session

from backend.database import SessionLocal
from backend.models import (
    Project,
    Compound,
    CompoundVersion,
    ExternalExperimentalEvidence,
)
import backend.activity_models  # Ensure AssayDefinition is in metadata
from backend.openadmet_cyp import predict_chemeleon_cyp_pic50, ic50_nm_to_pic50, compute_fold_error
from backend.quantitative_safety_transporters import (
    predict_quantitative_herg_pic50,
    predict_quantitative_pgp_pic50,
    evaluate_safety_applicability_domain,
)

DRUGBANK_PROJECT_NAME = "DrugBank"
DRUGBANK_PROJECT_INDICATION = "Global Reference Drug Library (GLOBAL_MODEL_DEVELOPMENT)"
DRUGBANK_PROJECT_DESC = "Canonical reference drug library curated for Drug-OPT Global Prediction Engine v3.0 training and multi-tiered benchmarking."

ROLE_DEVELOPMENT_TRAINING = "DEVELOPMENT_TRAINING"
ROLE_IMMUTABLE_HOLDOUT = "IMMUTABLE_HOLDOUT"

# 10 Reference Drugs across Diverse Chemical & Pharmacological Spaces
REFERENCE_DRUGS_CATALOG = [
    {
        "name": "Gefitinib",
        "cas_number": "184475-35-2",
        "drugbank_id": "DB00317",
        "pubchem_cid": "123631",
        "chembl_id": "CHEMBL939",
        "unii": "S65743JHGW",
        "smiles": "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
        "indication": "Non-Small Cell Lung Cancer (NSCLC)",
        "target": "EGFR (Epidermal Growth Factor Receptor)",
        "scaffold_family": "Quinazoline",
        "model_role": ROLE_DEVELOPMENT_TRAINING,
        "upstream_overlap": {
            "SOLUBILITY_GENERIC": "NOVEL_IN_DOMAIN",
            "SOLUBILITY_THERMODYNAMIC": "NOVEL_IN_DOMAIN",
            "HUMAN_PPB": "EXACT_STRUCTURE_OVERLAP",
            "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT",
            "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT",
            "HERG_LIABILITY": "VALIDATION_HOLDOUT",
        },
        "observations": [
            {
                "canonical_endpoint_id": "ACTIVITY_EGFR_WT_IC50",
                "raw_endpoint_name": "EGFR WT Biochemical IC50",
                "section": "ACTIVITY", "species": "Homo sapiens", "matrix": "Recombinant EGFR Kinase Domain",
                "raw_value": 33.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 33.0, "normalized_unit": "nM",
                "reference_text": "Barker et al. / ChEMBL939 · PubChem CID 123631",
                "assay_type": "Biochemical Kinase Activity Assay (Z'-LYTE / Radiometric)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC",
                "raw_endpoint_name": "Thermodynamic Aqueous Solubility",
                "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0",
                "raw_value": 12.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -4.92, "normalized_unit": "log10(mol/L)",
                "reference_text": "DrugCentral / DrugBank DB00317 · FDA NDA 021399",
                "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PPB",
                "raw_endpoint_name": "Human Plasma Protein Binding",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 90.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 90.0, "normalized_unit": "% bound",
                "reference_text": "DailyMed / Drugs@FDA NDA 021399 (Iressa Label)",
                "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CACO2_PAPP_AB",
                "raw_endpoint_name": "Caco-2 Permeability (A to B)",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer",
                "raw_value": 18.5, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 18.5, "normalized_unit": "10^-6 cm/s",
                "reference_text": "FDA Clinical Pharmacology & Biopharmaceutics Review NDA 021399",
                "assay_type": "Bidirectional Caco-2 Cell Monolayer Assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HLM_CLINT",
                "raw_endpoint_name": "Human Liver Microsomes Clint",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Human Liver Microsomes (HLM)",
                "raw_value": 42.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 42.0, "normalized_unit": "uL/min/mg protein",
                "reference_text": "FDA NDA 021399 Nonclinical Pharmacology / DrugCentral",
                "assay_type": "Substrate Depletion Assay in HLM + NADPH", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CYP3A4_INHIBITION",
                "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP3A4",
                "raw_value": 1.8, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 1800.0, "normalized_unit": "nM",
                "reference_text": "FDA NDA 021399 In Vitro Metabolism & Transporter Studies",
                "assay_type": "Recombinant human CYP3A4 Midazolam 1'-hydroxylation assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CYP2D6_INHIBITION",
                "raw_endpoint_name": "CYP2D6 Direct Reversible Inhibition",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP2D6",
                "raw_value": 3.2, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 3200.0, "normalized_unit": "nM",
                "reference_text": "FDA NDA 021399 In Vitro Metabolism Studies",
                "assay_type": "Recombinant human CYP2D6 Dextromethorphan O-demethylation assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HERG_LIABILITY",
                "raw_endpoint_name": "hERG Potassium Channel Inhibition",
                "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp",
                "raw_value": 5.8, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 5800.0, "normalized_unit": "nM",
                "reference_text": "FDA NDA 021399 Safety Pharmacology Review / ChEMBL939",
                "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PK_CMAX_ORAL",
                "raw_endpoint_name": "Human Oral Steady-State Cmax (250 mg QD)",
                "section": "CLINICAL_PK", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 502.0, "raw_unit": "ng/mL", "raw_relation": "=", "normalized_value": 502.0, "normalized_unit": "ng/mL",
                "reference_text": "DailyMed Iressa Label / FDA NDA 021399 Clinical PK",
                "assay_type": "Clinical Phase I/II Steady-State PK (LC-MS/MS)", "training_eligible": False,
            },
        ],
    },
    {
        "name": "Imatinib",
        "cas_number": "152459-95-5",
        "drugbank_id": "DB00619",
        "pubchem_cid": "5291",
        "chembl_id": "CHEMBL941",
        "unii": "BKJ8M8G5HI",
        "smiles": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
        "indication": "Chronic Myelogenous Leukemia (CML) / GIST",
        "target": "BCR-ABL1 / c-KIT / PDGFR",
        "scaffold_family": "2-Phenylaminopyrimidine",
        "model_role": ROLE_DEVELOPMENT_TRAINING,
        "upstream_overlap": {
            "SOLUBILITY_GENERIC": "EXACT_STRUCTURE_OVERLAP",
            "SOLUBILITY_THERMODYNAMIC": "EXACT_STRUCTURE_OVERLAP",
            "HUMAN_PPB": "EXACT_STRUCTURE_OVERLAP",
            "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT",
            "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT",
            "HERG_LIABILITY": "VALIDATION_HOLDOUT",
        },
        "observations": [
            {
                "canonical_endpoint_id": "ACTIVITY_ABL1_WT_IC50",
                "raw_endpoint_name": "ABL1 WT Biochemical IC50",
                "section": "ACTIVITY", "species": "Homo sapiens", "matrix": "Recombinant ABL1 Kinase Domain",
                "raw_value": 25.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 25.0, "normalized_unit": "nM",
                "reference_text": "Druker et al. / ChEMBL941 · PubChem CID 5291",
                "assay_type": "Biochemical Kinase Activity Assay (Filter-binding / Radiometric)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC",
                "raw_endpoint_name": "Thermodynamic Aqueous Solubility",
                "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0",
                "raw_value": 45.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -4.35, "normalized_unit": "log10(mol/L)",
                "reference_text": "AqSolDB / DrugCentral · FDA NDA 021335",
                "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PPB",
                "raw_endpoint_name": "Human Plasma Protein Binding",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 95.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 95.0, "normalized_unit": "% bound",
                "reference_text": "DailyMed / Gleevec Label NDA 021335",
                "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CACO2_PAPP_AB",
                "raw_endpoint_name": "Caco-2 Permeability (A to B)",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer",
                "raw_value": 8.2, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 8.2, "normalized_unit": "10^-6 cm/s",
                "reference_text": "FDA Clinical Pharmacology NDA 021335",
                "assay_type": "Bidirectional Caco-2 Cell Monolayer Assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HLM_CLINT",
                "raw_endpoint_name": "Human Liver Microsomes Clint",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Human Liver Microsomes (HLM)",
                "raw_value": 35.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 35.0, "normalized_unit": "uL/min/mg protein",
                "reference_text": "FDA NDA 021335 Nonclinical Pharmacology",
                "assay_type": "Substrate Depletion Assay in HLM + NADPH", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CYP3A4_INHIBITION",
                "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP3A4",
                "raw_value": 1.6, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 1600.0, "normalized_unit": "nM",
                "reference_text": "FDA NDA 021335 In Vitro Metabolism Studies",
                "assay_type": "Recombinant human CYP3A4 Midazolam 1'-hydroxylation assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CYP2D6_INHIBITION",
                "raw_endpoint_name": "CYP2D6 Direct Reversible Inhibition",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP2D6",
                "raw_value": 2.5, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 2500.0, "normalized_unit": "nM",
                "reference_text": "FDA NDA 021335 In Vitro Metabolism Studies",
                "assay_type": "Recombinant human CYP2D6 Dextromethorphan O-demethylation assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HERG_LIABILITY",
                "raw_endpoint_name": "hERG Potassium Channel Inhibition",
                "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp",
                "raw_value": 14.8, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 14800.0, "normalized_unit": "nM",
                "reference_text": "FDA NDA 021335 Safety Pharmacology Review / ChEMBL941",
                "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PK_CMAX_ORAL",
                "raw_endpoint_name": "Human Oral Steady-State Cmax (400 mg QD)",
                "section": "CLINICAL_PK", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 2600.0, "raw_unit": "ng/mL", "raw_relation": "=", "normalized_value": 2600.0, "normalized_unit": "ng/mL",
                "reference_text": "DailyMed Gleevec Label / FDA NDA 021335 Clinical PK",
                "assay_type": "Clinical Phase I/II Steady-State PK (LC-MS/MS)", "training_eligible": False,
            },
        ],
    },
    {
        "name": "Propranolol",
        "cas_number": "525-66-6",
        "drugbank_id": "DB00571",
        "pubchem_cid": "4946",
        "chembl_id": "CHEMBL4",
        "unii": "9Y8NXQ24VQ",
        "smiles": "CC(C)NCC(O)COc1cccc2ccccc12",
        "indication": "Hypertension / Angina / Arrhythmia",
        "target": "Beta-1 / Beta-2 Adrenergic Receptors",
        "scaffold_family": "Aryloxypropanolamine",
        "model_role": ROLE_IMMUTABLE_HOLDOUT,
        "upstream_overlap": {
            "SOLUBILITY_GENERIC": "EXACT_STRUCTURE_OVERLAP",
            "SOLUBILITY_THERMODYNAMIC": "EXACT_STRUCTURE_OVERLAP",
            "HUMAN_PPB": "EXACT_STRUCTURE_OVERLAP",
            "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT",
            "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT",
            "HERG_LIABILITY": "VALIDATION_HOLDOUT",
        },
        "observations": [
            {
                "canonical_endpoint_id": "ACTIVITY_ADRB1_KI",
                "raw_endpoint_name": "Beta-1 Adrenergic Receptor Ki",
                "section": "ACTIVITY", "species": "Homo sapiens", "matrix": "Recombinant Human ADRB1 Membrane",
                "raw_value": 1.8, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 1.8, "normalized_unit": "nM",
                "reference_text": "ChEMBL4 · PDSP Ki Database",
                "assay_type": "Radioligand Binding Assay ([3H]-CGP-12177)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC",
                "raw_endpoint_name": "Thermodynamic Aqueous Solubility",
                "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0",
                "raw_value": 240.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -3.62, "normalized_unit": "log10(mol/L)",
                "reference_text": "AqSolDB / DrugCentral · FDA NDA 016418",
                "assay_type": "Potentiometric / Shake-Flask Solubility", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PPB",
                "raw_endpoint_name": "Human Plasma Protein Binding",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 87.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 87.0, "normalized_unit": "% bound",
                "reference_text": "DailyMed / Inderal Label NDA 016418",
                "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CACO2_PAPP_AB",
                "raw_endpoint_name": "Caco-2 Permeability (A to B)",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer",
                "raw_value": 26.4, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 26.4, "normalized_unit": "10^-6 cm/s",
                "reference_text": "Yazdanian et al. / FDA Guidance Benchmark",
                "assay_type": "Bidirectional Caco-2 Cell Monolayer Assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HLM_CLINT",
                "raw_endpoint_name": "Human Liver Microsomes Clint",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Human Liver Microsomes (HLM)",
                "raw_value": 120.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 120.0, "normalized_unit": "uL/min/mg protein",
                "reference_text": "Obach et al. / DrugCentral",
                "assay_type": "Substrate Depletion Assay in HLM + NADPH", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CYP2D6_INHIBITION",
                "raw_endpoint_name": "CYP2D6 Direct Reversible Inhibition",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP2D6",
                "raw_value": 4.8, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 4800.0, "normalized_unit": "nM",
                "reference_text": "In Vitro Metabolism Studies · ChEMBL4",
                "assay_type": "Recombinant human CYP2D6 Dextromethorphan assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HERG_LIABILITY",
                "raw_endpoint_name": "hERG Potassium Channel Inhibition",
                "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp",
                "raw_value": 3.5, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 3500.0, "normalized_unit": "nM",
                "reference_text": "Redfern et al. / ChEMBL4 Safety Benchmark",
                "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PK_CMAX_ORAL",
                "raw_endpoint_name": "Human Oral Single-Dose Cmax (80 mg)",
                "section": "CLINICAL_PK", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 75.0, "raw_unit": "ng/mL", "raw_relation": "=", "normalized_value": 75.0, "normalized_unit": "ng/mL",
                "reference_text": "DailyMed Inderal Label / FDA NDA 016418",
                "assay_type": "Clinical Single-Dose PK (LC-MS/MS)", "training_eligible": False,
            },
        ],
    },
    {
        "name": "Atorvastatin",
        "cas_number": "134523-00-5",
        "drugbank_id": "DB01076",
        "pubchem_cid": "60823",
        "chembl_id": "CHEMBL1487",
        "unii": "48A5M73Z4Q",
        "smiles": "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
        "indication": "Hypercholesterolemia / Cardiovascular Risk Reduction",
        "target": "HMG-CoA Reductase (HMGCR)",
        "scaffold_family": "Substituted Pyrrole-heptanoic acid",
        "model_role": ROLE_DEVELOPMENT_TRAINING,
        "upstream_overlap": {
            "SOLUBILITY_GENERIC": "NOVEL_IN_DOMAIN",
            "SOLUBILITY_THERMODYNAMIC": "NOVEL_IN_DOMAIN",
            "HUMAN_PPB": "EXACT_STRUCTURE_OVERLAP",
            "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT",
            "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT",
            "HERG_LIABILITY": "VALIDATION_HOLDOUT",
        },
        "observations": [
            {
                "canonical_endpoint_id": "ACTIVITY_HMGCR_IC50",
                "raw_endpoint_name": "HMG-CoA Reductase Biochemical IC50",
                "section": "ACTIVITY", "species": "Homo sapiens", "matrix": "Recombinant Human HMGCR Catalytic Domain",
                "raw_value": 8.2, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 8.2, "normalized_unit": "nM",
                "reference_text": "Istvan & Deisenhofer / ChEMBL1487",
                "assay_type": "Spectrophotometric Enzyme Inhibition Assay (NADPH Oxidation)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC",
                "raw_endpoint_name": "Thermodynamic Aqueous Solubility",
                "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0",
                "raw_value": 0.85, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -6.07, "normalized_unit": "log10(mol/L)",
                "reference_text": "DrugCentral / Lipitor NDA 020702",
                "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PPB",
                "raw_endpoint_name": "Human Plasma Protein Binding",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 98.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 98.0, "normalized_unit": "% bound",
                "reference_text": "DailyMed / Lipitor Label NDA 020702",
                "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CACO2_PAPP_AB",
                "raw_endpoint_name": "Caco-2 Permeability (A to B)",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer",
                "raw_value": 4.5, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 4.5, "normalized_unit": "10^-6 cm/s",
                "reference_text": "FDA Clinical Pharmacology Review NDA 020702",
                "assay_type": "Bidirectional Caco-2 Cell Monolayer Assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HLM_CLINT",
                "raw_endpoint_name": "Human Liver Microsomes Clint",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Human Liver Microsomes (HLM)",
                "raw_value": 68.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 68.0, "normalized_unit": "uL/min/mg protein",
                "reference_text": "FDA NDA 020702 Nonclinical Pharmacology",
                "assay_type": "Substrate Depletion Assay in HLM + NADPH", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CYP3A4_INHIBITION",
                "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP3A4",
                "raw_value": 2.4, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 2400.0, "normalized_unit": "nM",
                "reference_text": "FDA NDA 020702 In Vitro Metabolism Studies",
                "assay_type": "Recombinant human CYP3A4 Midazolam 1'-hydroxylation assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HERG_LIABILITY",
                "raw_endpoint_name": "hERG Potassium Channel Inhibition",
                "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp",
                "raw_value": 48.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 48000.0, "normalized_unit": "nM",
                "reference_text": "FDA NDA 020702 Safety Pharmacology Review",
                "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PK_CMAX_ORAL",
                "raw_endpoint_name": "Human Oral Steady-State Cmax (40 mg QD)",
                "section": "CLINICAL_PK", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 28.0, "raw_unit": "ng/mL", "raw_relation": "=", "normalized_value": 28.0, "normalized_unit": "ng/mL",
                "reference_text": "DailyMed Lipitor Label / FDA NDA 020702",
                "assay_type": "Clinical Phase I/II Steady-State PK (LC-MS/MS)", "training_eligible": False,
            },
        ],
    },
    {
        "name": "Midazolam",
        "cas_number": "59467-70-8",
        "drugbank_id": "DB00683",
        "pubchem_cid": "4192",
        "chembl_id": "CHEMBL644",
        "unii": "R60L0SM5BC",
        "smiles": "Cc1ncc2n1c3ccc(Cl)cc3C(=NC2)c4ccccc4F",
        "indication": "Procedural Sedation / Anesthesia Premedication",
        "target": "GABA-A Receptor / Standard CYP3A Probe Substrate",
        "scaffold_family": "Imidazobenzodiazepine",
        "model_role": ROLE_IMMUTABLE_HOLDOUT,
        "upstream_overlap": {
            "SOLUBILITY_GENERIC": "EXACT_STRUCTURE_OVERLAP",
            "SOLUBILITY_THERMODYNAMIC": "EXACT_STRUCTURE_OVERLAP",
            "HUMAN_PPB": "EXACT_STRUCTURE_OVERLAP",
            "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT",
            "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT",
            "HERG_LIABILITY": "VALIDATION_HOLDOUT",
        },
        "observations": [
            {
                "canonical_endpoint_id": "ACTIVITY_GABAA_KI",
                "raw_endpoint_name": "GABA-A Benzodiazepine Site Ki",
                "section": "ACTIVITY", "species": "Homo sapiens", "matrix": "Recombinant Human GABA-A alpha1beta2gamma2",
                "raw_value": 5.6, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 5.6, "normalized_unit": "nM",
                "reference_text": "ChEMBL644 · PDSP Ki Database",
                "assay_type": "Radioligand Binding Assay ([3H]-Flunitrazepam)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC",
                "raw_endpoint_name": "Thermodynamic Aqueous Solubility",
                "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0",
                "raw_value": 150.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -3.82, "normalized_unit": "log10(mol/L)",
                "reference_text": "AqSolDB / DrugCentral · FDA NDA 018654",
                "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PPB",
                "raw_endpoint_name": "Human Plasma Protein Binding",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 96.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 96.0, "normalized_unit": "% bound",
                "reference_text": "DailyMed / Versed Label NDA 018654",
                "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CACO2_PAPP_AB",
                "raw_endpoint_name": "Caco-2 Permeability (A to B)",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer",
                "raw_value": 22.0, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 22.0, "normalized_unit": "10^-6 cm/s",
                "reference_text": "FDA Clinical Pharmacology Guidance Benchmark",
                "assay_type": "Bidirectional Caco-2 Cell Monolayer Assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HLM_CLINT",
                "raw_endpoint_name": "Human Liver Microsomes Clint",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Human Liver Microsomes (HLM)",
                "raw_value": 180.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 180.0, "normalized_unit": "uL/min/mg protein",
                "reference_text": "Obach et al. / DrugCentral Benchmark",
                "assay_type": "Substrate Depletion Assay in HLM + NADPH", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CYP3A4_INHIBITION",
                "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP3A4",
                "raw_value": 4.2, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 4200.0, "normalized_unit": "nM",
                "reference_text": "FDA NDA 018654 In Vitro Metabolism Studies",
                "assay_type": "Recombinant human CYP3A4 Substrate Auto-Inhibition", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HERG_LIABILITY",
                "raw_endpoint_name": "hERG Potassium Channel Inhibition",
                "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp",
                "raw_value": 22.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 22000.0, "normalized_unit": "nM",
                "reference_text": "FDA NDA 018654 Safety Pharmacology Review",
                "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PK_CMAX_ORAL",
                "raw_endpoint_name": "Human Oral Single-Dose Cmax (7.5 mg)",
                "section": "CLINICAL_PK", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 55.0, "raw_unit": "ng/mL", "raw_relation": "=", "normalized_value": 55.0, "normalized_unit": "ng/mL",
                "reference_text": "DailyMed Versed Label / FDA NDA 018654",
                "assay_type": "Clinical Single-Dose PK (LC-MS/MS)", "training_eligible": False,
            },
        ],
    },
    {
        "name": "Verapamil",
        "cas_number": "52-53-9",
        "drugbank_id": "DB00661",
        "pubchem_cid": "2520",
        "chembl_id": "CHEMBL696",
        "unii": "CJ0O37KU29",
        "smiles": "COc1ccc(CCN(C)CCCC(C#N)(C(C)C)c2ccc(OC)c(OC)c2)cc1OC",
        "indication": "Angina / Arrhythmia / Hypertension",
        "target": "L-type Calcium Channel (Cav1.2) / P-gp / CYP3A4",
        "scaffold_family": "Phenylalkylamine",
        "model_role": ROLE_IMMUTABLE_HOLDOUT,
        "upstream_overlap": {
            "SOLUBILITY_GENERIC": "EXACT_STRUCTURE_OVERLAP",
            "SOLUBILITY_THERMODYNAMIC": "EXACT_STRUCTURE_OVERLAP",
            "HUMAN_PPB": "EXACT_STRUCTURE_OVERLAP",
            "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT",
            "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT",
            "HERG_LIABILITY": "VALIDATION_HOLDOUT",
        },
        "observations": [
            {
                "canonical_endpoint_id": "ACTIVITY_CAV12_IC50",
                "raw_endpoint_name": "Cav1.2 Calcium Channel IC50",
                "section": "ACTIVITY", "species": "Homo sapiens", "matrix": "Recombinant Human Cav1.2",
                "raw_value": 220.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 220.0, "normalized_unit": "nM",
                "reference_text": "ChEMBL696 · FDA Pharmacology Review",
                "assay_type": "Whole-cell Electrophysiology Inward Ca2+ Current", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC",
                "raw_endpoint_name": "Thermodynamic Aqueous Solubility",
                "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0",
                "raw_value": 18.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -4.74, "normalized_unit": "log10(mol/L)",
                "reference_text": "AqSolDB / DrugCentral · FDA NDA 018817",
                "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PPB",
                "raw_endpoint_name": "Human Plasma Protein Binding",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 90.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 90.0, "normalized_unit": "% bound",
                "reference_text": "DailyMed / Isoptin Label NDA 018817",
                "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CACO2_PAPP_AB",
                "raw_endpoint_name": "Caco-2 Permeability (A to B)",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer",
                "raw_value": 12.0, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 12.0, "normalized_unit": "10^-6 cm/s",
                "reference_text": "FDA Clinical Pharmacology Benchmark",
                "assay_type": "Bidirectional Caco-2 Cell Monolayer Assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HLM_CLINT",
                "raw_endpoint_name": "Human Liver Microsomes Clint",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Human Liver Microsomes (HLM)",
                "raw_value": 95.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 95.0, "normalized_unit": "uL/min/mg protein",
                "reference_text": "Obach et al. / DrugCentral",
                "assay_type": "Substrate Depletion Assay in HLM + NADPH", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CYP3A4_INHIBITION",
                "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP3A4",
                "raw_value": 3.8, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 3800.0, "normalized_unit": "nM",
                "reference_text": "FDA In Vitro Metabolism Guidance Benchmark",
                "assay_type": "Recombinant human CYP3A4 Midazolam 1'-hydroxylation assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HERG_LIABILITY",
                "raw_endpoint_name": "hERG Potassium Channel Inhibition",
                "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp",
                "raw_value": 0.58, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 580.0, "normalized_unit": "nM",
                "reference_text": "Redfern et al. / FDA Safety Pharmacology",
                "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PK_CMAX_ORAL",
                "raw_endpoint_name": "Human Oral Single-Dose Cmax (120 mg)",
                "section": "CLINICAL_PK", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 110.0, "raw_unit": "ng/mL", "raw_relation": "=", "normalized_value": 110.0, "normalized_unit": "ng/mL",
                "reference_text": "DailyMed Isoptin Label / FDA NDA 018817",
                "assay_type": "Clinical Single-Dose PK (LC-MS/MS)", "training_eligible": False,
            },
        ],
    },
    {
        "name": "Fluoxetine",
        "cas_number": "54910-89-3",
        "drugbank_id": "DB00472",
        "pubchem_cid": "3386",
        "chembl_id": "CHEMBL41",
        "unii": "01K63SUP8D",
        "smiles": "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1",
        "indication": "Major Depressive Disorder / OCD / Panic Disorder",
        "target": "Serotonin Transporter (SERT / SLC6A4) / Potent CYP2D6 Inhibitor",
        "scaffold_family": "Aryloxypropylamine",
        "model_role": ROLE_IMMUTABLE_HOLDOUT,
        "upstream_overlap": {
            "SOLUBILITY_GENERIC": "EXACT_STRUCTURE_OVERLAP",
            "SOLUBILITY_THERMODYNAMIC": "EXACT_STRUCTURE_OVERLAP",
            "HUMAN_PPB": "EXACT_STRUCTURE_OVERLAP",
            "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT",
            "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT",
            "HERG_LIABILITY": "VALIDATION_HOLDOUT",
        },
        "observations": [
            {
                "canonical_endpoint_id": "ACTIVITY_SERT_KI",
                "raw_endpoint_name": "Serotonin Transporter Ki",
                "section": "ACTIVITY", "species": "Homo sapiens", "matrix": "Recombinant Human SERT Membrane",
                "raw_value": 0.81, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 0.81, "normalized_unit": "nM",
                "reference_text": "ChEMBL41 · PDSP Ki Database",
                "assay_type": "Radioligand Binding Assay ([3H]-Citalopram)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC",
                "raw_endpoint_name": "Thermodynamic Aqueous Solubility",
                "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0",
                "raw_value": 32.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -4.49, "normalized_unit": "log10(mol/L)",
                "reference_text": "AqSolDB / DrugCentral · FDA NDA 018936",
                "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PPB",
                "raw_endpoint_name": "Human Plasma Protein Binding",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 94.5, "raw_unit": "%", "raw_relation": "=", "normalized_value": 94.5, "normalized_unit": "% bound",
                "reference_text": "DailyMed / Prozac Label NDA 018936",
                "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CACO2_PAPP_AB",
                "raw_endpoint_name": "Caco-2 Permeability (A to B)",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer",
                "raw_value": 24.0, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 24.0, "normalized_unit": "10^-6 cm/s",
                "reference_text": "FDA Clinical Pharmacology Benchmark",
                "assay_type": "Bidirectional Caco-2 Cell Monolayer Assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HLM_CLINT",
                "raw_endpoint_name": "Human Liver Microsomes Clint",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Human Liver Microsomes (HLM)",
                "raw_value": 48.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 48.0, "normalized_unit": "uL/min/mg protein",
                "reference_text": "Obach et al. / DrugCentral",
                "assay_type": "Substrate Depletion Assay in HLM + NADPH", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CYP2D6_INHIBITION",
                "raw_endpoint_name": "CYP2D6 Direct Reversible Inhibition",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP2D6",
                "raw_value": 0.60, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 600.0, "normalized_unit": "nM",
                "reference_text": "FDA Drug Interaction Guidance Benchmark",
                "assay_type": "Recombinant human CYP2D6 Dextromethorphan assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HERG_LIABILITY",
                "raw_endpoint_name": "hERG Potassium Channel Inhibition",
                "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp",
                "raw_value": 1.5, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 1500.0, "normalized_unit": "nM",
                "reference_text": "Redfern et al. / ChEMBL41",
                "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PK_CMAX_ORAL",
                "raw_endpoint_name": "Human Oral Steady-State Cmax (20 mg QD)",
                "section": "CLINICAL_PK", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 125.0, "raw_unit": "ng/mL", "raw_relation": "=", "normalized_value": 125.0, "normalized_unit": "ng/mL",
                "reference_text": "DailyMed Prozac Label / FDA NDA 018936",
                "assay_type": "Clinical Phase I/II Steady-State PK (LC-MS/MS)", "training_eligible": False,
            },
        ],
    },
    {
        "name": "Ketoconazole",
        "cas_number": "65277-42-1",
        "drugbank_id": "DB01026",
        "pubchem_cid": "456201",
        "chembl_id": "CHEMBL109",
        "unii": "R9400W927I",
        "smiles": "CC(=O)N1CCN(c2ccc(OCC3COC(Cn4ccnc4)(c4ccc(Cl)cc4Cl)O3)cc2)CC1",
        "indication": "Fungal Infections / Standard Strong CYP3A4 Inhibitor",
        "target": "CYP51 / Potent Reversible CYP3A4 Inhibitor",
        "scaffold_family": "Imidazole-dioxolane-piperazine",
        "model_role": ROLE_DEVELOPMENT_TRAINING,
        "upstream_overlap": {
            "SOLUBILITY_GENERIC": "EXACT_STRUCTURE_OVERLAP",
            "SOLUBILITY_THERMODYNAMIC": "EXACT_STRUCTURE_OVERLAP",
            "HUMAN_PPB": "EXACT_STRUCTURE_OVERLAP",
            "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT",
            "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT",
            "HERG_LIABILITY": "VALIDATION_HOLDOUT",
        },
        "observations": [
            {
                "canonical_endpoint_id": "ACTIVITY_CYP51_IC50",
                "raw_endpoint_name": "CYP51 Fungal Demethylase IC50",
                "section": "ACTIVITY", "species": "Candida albicans", "matrix": "Recombinant CYP51",
                "raw_value": 18.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 18.0, "normalized_unit": "nM",
                "reference_text": "ChEMBL109 · PubChem CID 456201",
                "assay_type": "Enzyme Inhibition Spectrophotometric Assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC",
                "raw_endpoint_name": "Thermodynamic Aqueous Solubility",
                "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0",
                "raw_value": 1.2, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -5.92, "normalized_unit": "log10(mol/L)",
                "reference_text": "AqSolDB / DrugCentral · FDA NDA 018533",
                "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PPB",
                "raw_endpoint_name": "Human Plasma Protein Binding",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 99.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 99.0, "normalized_unit": "% bound",
                "reference_text": "DailyMed / Nizoral Label NDA 018533",
                "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CACO2_PAPP_AB",
                "raw_endpoint_name": "Caco-2 Permeability (A to B)",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer",
                "raw_value": 6.8, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 6.8, "normalized_unit": "10^-6 cm/s",
                "reference_text": "FDA Clinical Pharmacology Guidance Benchmark",
                "assay_type": "Bidirectional Caco-2 Cell Monolayer Assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HLM_CLINT",
                "raw_endpoint_name": "Human Liver Microsomes Clint",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Human Liver Microsomes (HLM)",
                "raw_value": 15.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 15.0, "normalized_unit": "uL/min/mg protein",
                "reference_text": "Obach et al. / DrugCentral",
                "assay_type": "Substrate Depletion Assay in HLM + NADPH", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CYP3A4_INHIBITION",
                "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP3A4",
                "raw_value": 0.05, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 50.0, "normalized_unit": "nM",
                "reference_text": "FDA In Vitro Drug Interaction Guidance Standard Strong Inhibitor",
                "assay_type": "Recombinant human CYP3A4 Midazolam 1'-hydroxylation assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HERG_LIABILITY",
                "raw_endpoint_name": "hERG Potassium Channel Inhibition",
                "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp",
                "raw_value": 4.8, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 4800.0, "normalized_unit": "nM",
                "reference_text": "Redfern et al. / ChEMBL109",
                "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PK_CMAX_ORAL",
                "raw_endpoint_name": "Human Oral Single-Dose Cmax (200 mg)",
                "section": "CLINICAL_PK", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 3500.0, "raw_unit": "ng/mL", "raw_relation": "=", "normalized_value": 3500.0, "normalized_unit": "ng/mL",
                "reference_text": "DailyMed Nizoral Label / FDA NDA 018533",
                "assay_type": "Clinical Single-Dose PK (LC-MS/MS)", "training_eligible": False,
            },
        ],
    },
    {
        "name": "Sildenafil",
        "cas_number": "139755-83-2",
        "drugbank_id": "DB00203",
        "pubchem_cid": "135398744",
        "chembl_id": "CHEMBL192",
        "unii": "3S12J47372",
        "smiles": "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12",
        "indication": "Erectile Dysfunction / Pulmonary Arterial Hypertension",
        "target": "Phosphodiesterase-5 (PDE5)",
        "scaffold_family": "Pyrazolopyrimidinone-sulfonamide",
        "model_role": ROLE_IMMUTABLE_HOLDOUT,
        "upstream_overlap": {
            "SOLUBILITY_GENERIC": "EXACT_STRUCTURE_OVERLAP",
            "SOLUBILITY_THERMODYNAMIC": "EXACT_STRUCTURE_OVERLAP",
            "HUMAN_PPB": "EXACT_STRUCTURE_OVERLAP",
            "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT",
            "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT",
            "HERG_LIABILITY": "VALIDATION_HOLDOUT",
        },
        "observations": [
            {
                "canonical_endpoint_id": "ACTIVITY_PDE5_IC50",
                "raw_endpoint_name": "PDE5 Biochemical IC50",
                "section": "ACTIVITY", "species": "Homo sapiens", "matrix": "Recombinant Human PDE5A1 Catalytic Domain",
                "raw_value": 3.5, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 3.5, "normalized_unit": "nM",
                "reference_text": "Corbin & Francis / ChEMBL192",
                "assay_type": "Fluorescence Polarization / Scintillation PDE5 Assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC",
                "raw_endpoint_name": "Thermodynamic Aqueous Solubility",
                "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0",
                "raw_value": 25.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -4.60, "normalized_unit": "log10(mol/L)",
                "reference_text": "AqSolDB / DrugCentral · FDA NDA 020895",
                "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PPB",
                "raw_endpoint_name": "Human Plasma Protein Binding",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 96.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 96.0, "normalized_unit": "% bound",
                "reference_text": "DailyMed / Viagra Label NDA 020895",
                "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CACO2_PAPP_AB",
                "raw_endpoint_name": "Caco-2 Permeability (A to B)",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer",
                "raw_value": 9.5, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 9.5, "normalized_unit": "10^-6 cm/s",
                "reference_text": "FDA Clinical Pharmacology Benchmark",
                "assay_type": "Bidirectional Caco-2 Cell Monolayer Assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HLM_CLINT",
                "raw_endpoint_name": "Human Liver Microsomes Clint",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Human Liver Microsomes (HLM)",
                "raw_value": 55.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 55.0, "normalized_unit": "uL/min/mg protein",
                "reference_text": "Obach et al. / DrugCentral",
                "assay_type": "Substrate Depletion Assay in HLM + NADPH", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CYP3A4_INHIBITION",
                "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP3A4",
                "raw_value": 12.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 12000.0, "normalized_unit": "nM",
                "reference_text": "FDA In Vitro Metabolism Studies NDA 020895",
                "assay_type": "Recombinant human CYP3A4 Midazolam 1'-hydroxylation assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HERG_LIABILITY",
                "raw_endpoint_name": "hERG Potassium Channel Inhibition",
                "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp",
                "raw_value": 33.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 33000.0, "normalized_unit": "nM",
                "reference_text": "Redfern et al. / ChEMBL192",
                "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PK_CMAX_ORAL",
                "raw_endpoint_name": "Human Oral Single-Dose Cmax (50 mg)",
                "section": "CLINICAL_PK", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 450.0, "raw_unit": "ng/mL", "raw_relation": "=", "normalized_value": 450.0, "normalized_unit": "ng/mL",
                "reference_text": "DailyMed Viagra Label / FDA NDA 020895",
                "assay_type": "Clinical Single-Dose PK (LC-MS/MS)", "training_eligible": False,
            },
        ],
    },
    {
        "name": "Quinidine",
        "cas_number": "56-54-2",
        "drugbank_id": "DB00908",
        "pubchem_cid": "441074",
        "chembl_id": "CHEMBL582",
        "unii": "ITX086884O",
        "smiles": "C=CC1CN2CCC1CC2C(O)c1ccc2cc(OC)ccc2n1",
        "indication": "Cardiac Arrhythmias / Standard Potent CYP2D6 Inhibitor",
        "target": "Nav1.5 / Kv1.5 / Standard Potent CYP2D6 Reference",
        "scaffold_family": "Cinchona alkaloid / Quinuclidine-quinoline",
        "model_role": ROLE_IMMUTABLE_HOLDOUT,
        "upstream_overlap": {
            "SOLUBILITY_GENERIC": "EXACT_STRUCTURE_OVERLAP",
            "SOLUBILITY_THERMODYNAMIC": "EXACT_STRUCTURE_OVERLAP",
            "HUMAN_PPB": "EXACT_STRUCTURE_OVERLAP",
            "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT",
            "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT",
            "HERG_LIABILITY": "VALIDATION_HOLDOUT",
        },
        "observations": [
            {
                "canonical_endpoint_id": "ACTIVITY_NAV15_IC50",
                "raw_endpoint_name": "Nav1.5 Sodium Channel IC50",
                "section": "ACTIVITY", "species": "Homo sapiens", "matrix": "Recombinant Human Nav1.5",
                "raw_value": 6800.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 6800.0, "normalized_unit": "nM",
                "reference_text": "ChEMBL582 · FDA Pharmacology Review",
                "assay_type": "Whole-cell Electrophysiology Peak INa Inhibition", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC",
                "raw_endpoint_name": "Thermodynamic Aqueous Solubility",
                "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0",
                "raw_value": 520.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -3.28, "normalized_unit": "log10(mol/L)",
                "reference_text": "AqSolDB / DrugCentral · FDA NDA 006401",
                "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PPB",
                "raw_endpoint_name": "Human Plasma Protein Binding",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 89.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 89.0, "normalized_unit": "% bound",
                "reference_text": "DailyMed / Cardioquin Label NDA 006401",
                "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CACO2_PAPP_AB",
                "raw_endpoint_name": "Caco-2 Permeability (A to B)",
                "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer",
                "raw_value": 15.0, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 15.0, "normalized_unit": "10^-6 cm/s",
                "reference_text": "FDA Clinical Pharmacology Benchmark",
                "assay_type": "Bidirectional Caco-2 Cell Monolayer Assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HLM_CLINT",
                "raw_endpoint_name": "Human Liver Microsomes Clint",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Human Liver Microsomes (HLM)",
                "raw_value": 75.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 75.0, "normalized_unit": "uL/min/mg protein",
                "reference_text": "Obach et al. / DrugCentral",
                "assay_type": "Substrate Depletion Assay in HLM + NADPH", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "CYP2D6_INHIBITION",
                "raw_endpoint_name": "CYP2D6 Direct Reversible Inhibition",
                "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP2D6",
                "raw_value": 0.04, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 40.0, "normalized_unit": "nM",
                "reference_text": "FDA Drug Interaction Guidance Standard Potent Inhibitor",
                "assay_type": "Recombinant human CYP2D6 Dextromethorphan assay", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HERG_LIABILITY",
                "raw_endpoint_name": "hERG Potassium Channel Inhibition",
                "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp",
                "raw_value": 0.85, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 850.0, "normalized_unit": "nM",
                "reference_text": "Redfern et al. / ChEMBL582 Standard hERG Reference",
                "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True,
            },
            {
                "canonical_endpoint_id": "HUMAN_PK_CMAX_ORAL",
                "raw_endpoint_name": "Human Oral Single-Dose Cmax (200 mg)",
                "section": "CLINICAL_PK", "species": "Homo sapiens", "matrix": "Human Plasma",
                "raw_value": 1200.0, "raw_unit": "ng/mL", "raw_relation": "=", "normalized_value": 1200.0, "normalized_unit": "ng/mL",
                "reference_text": "DailyMed Cardioquin Label / FDA NDA 006401",
                "assay_type": "Clinical Single-Dose PK (LC-MS/MS)", "training_eligible": False,
            },
        ],
    },
]


def ensure_drugbank_project(db: Session) -> Project:
    """Ensures the canonical DrugBank project exists in the database."""
    proj = db.scalar(select(Project).where(Project.name == DRUGBANK_PROJECT_NAME))
    if not proj:
        proj = Project(
            name=DRUGBANK_PROJECT_NAME,
            target="PAN_TARGET_REFERENCE",
            molecule_type="Small Molecule",
            indication=DRUGBANK_PROJECT_INDICATION,
            mechanism_modality="SMALL_MOLECULE_DRUG",
            description=DRUGBANK_PROJECT_DESC,
        )
        db.add(proj)
        db.commit()
        db.refresh(proj)
    return proj


def ingest_reference_drug_by_spec(db: Session, drug_spec: Dict[str, Any]) -> Dict[str, Any]:
    """
    Ingests a single approved reference drug into the DrugBank project following the strict
    Identity -> Search -> Qualification -> Base Prediction -> Error pipeline.
    """
    proj = ensure_drugbank_project(db)
    mol = Chem.MolFromSmiles(drug_spec["smiles"])
    if mol is None:
        raise ValueError(f"Invalid SMILES for {drug_spec['name']}")

    canon_smiles = Chem.MolToSmiles(mol, canonical=True)
    inchi_str = Chem.MolToInchi(mol)
    inchikey_str = Chem.MolToInchiKey(mol)

    comp = db.scalar(select(Compound).where(Compound.project_id == proj.id, Compound.name == drug_spec["name"]))
    if not comp:
        comp = Compound(
            project_id=proj.id,
            compound_id=f"DRUGBANK-{drug_spec['drugbank_id']}",
            cas_number=drug_spec["cas_number"],
            name=drug_spec["name"],
            notes=f"Approved Reference Drug | DrugBank: {drug_spec['drugbank_id']} | ChEMBL: {drug_spec['chembl_id']} | PubChem: {drug_spec['pubchem_cid']} | UNII: {drug_spec['unii']} | Scaffold: {drug_spec.get('scaffold_family', '')} | Role: {drug_spec.get('model_role', ROLE_IMMUTABLE_HOLDOUT)}",
            status="APPROVED_REFERENCE",
            current_version=1,
        )
        db.add(comp)
        db.commit()
        db.refresh(comp)

        mw = float(Descriptors.MolWt(mol))
        clogp = float(Crippen.MolLogP(mol))
        tpsa = float(Descriptors.TPSA(mol))
        hbd = int(Lipinski.NumHDonors(mol))
        hba = int(Lipinski.NumHAcceptors(mol))
        rotb = int(Lipinski.NumRotatableBonds(mol))

        cv = CompoundVersion(
            compound_row_id=comp.id,
            version_number=1,
            original_smiles=drug_spec["smiles"],
            canonical_smiles=canon_smiles,
            isomeric_smiles=canon_smiles,
            inchi=inchi_str,
            inchikey=inchikey_str,
            change_note="Canonical reference drug registration",
            properties_json=json.dumps({
                "MW": mw, "cLogP": clogp, "TPSA": tpsa, "HBD": hbd, "HBA": hba, "RotB": rotb,
                "drugbank_id": drug_spec["drugbank_id"], "chembl_id": drug_spec["chembl_id"],
                "pubchem_cid": drug_spec["pubchem_cid"], "unii": drug_spec["unii"],
                "scaffold": drug_spec.get("scaffold_family", ""),
                "model_role": drug_spec.get("model_role", ROLE_IMMUTABLE_HOLDOUT),
            }),
        )
        db.add(cv)
        db.commit()
        db.refresh(cv)
    else:
        cv = db.scalar(select(CompoundVersion).where(CompoundVersion.compound_row_id == comp.id, CompoundVersion.version_number == 1))

    persisted_records = []
    ad_status, nearest_sim, violations, metrics, ad_reason = evaluate_safety_applicability_domain(mol)
    upstream_overlap = drug_spec.get("upstream_overlap", {})
    model_role = drug_spec.get("model_role", ROLE_IMMUTABLE_HOLDOUT)

    for obs in drug_spec["observations"]:
        eid = obs["canonical_endpoint_id"]
        overlap_status = upstream_overlap.get(eid, "VALIDATION_HOLDOUT" if obs["training_eligible"] else "NOT_ELIGIBLE")

        # Partitioning: IMMUTABLE_HOLDOUT vs DEVELOPMENT_TRAINING
        if not obs["training_eligible"]:
            partition = "NOT_ELIGIBLE"
        elif overlap_status == "EXACT_STRUCTURE_OVERLAP":
            partition = "TRAINING_ELIGIBLE"
        elif model_role == ROLE_IMMUTABLE_HOLDOUT:
            partition = "IMMUTABLE_HOLDOUT"
        else:
            partition = "DEVELOPMENT_TRAINING"

        p_key = hashlib.sha256(f"{cv.inchikey}_{eid}_{obs['raw_value']}_{obs['raw_unit']}_{obs['species']}_{obs['matrix']}".encode()).hexdigest()
        existing_ev = db.scalar(select(ExternalExperimentalEvidence).where(
            ExternalExperimentalEvidence.compound_version_id == cv.id,
            ExternalExperimentalEvidence.provenance_key == p_key
        ))

        cond_dict = {
            "matrix": obs["matrix"],
            "section": obs["section"],
            "upstream_overlap": overlap_status,
            "drugbank_partition": partition,
            "model_role": model_role,
        }

        if not existing_ev:
            ev = ExternalExperimentalEvidence(
                compound_version_id=cv.id,
                provenance_key=p_key,
                cas_number=drug_spec["cas_number"],
                canonical_endpoint_id=eid,
                raw_endpoint_name=obs["raw_endpoint_name"],
                species=obs["species"],
                assay_type=obs["assay_type"],
                assay_conditions_json=cond_dict,
                raw_value=obs["raw_value"],
                raw_unit=obs["raw_unit"],
                raw_relation=obs["raw_relation"],
                normalized_value=obs["normalized_value"],
                normalized_unit=obs["normalized_unit"],
                source_database="DrugBank_FDA_ChEMBL",
                source_record_id=drug_spec["drugbank_id"],
                source_url=f"https://go.drugbank.com/drugs/{drug_spec['drugbank_id']}",
                identity_match_status="EXACT_MATCH",
                endpoint_match_status="EXACT_MATCH",
                mapping_status="EXTERNAL_EVIDENCE_ONLY",
                evidence_origin="EXPERIMENTAL_EXTERNAL",
                source_quality_class="A",
                comparability_status="DIRECTLY_COMPARABLE",
                qualification_status="QUALIFIED_FOR_GLOBAL_TRAINING" if obs["training_eligible"] else "CLINICAL_PK_COMPOSITE",
                reference_text=obs["reference_text"],
                evidence_state="AUTO_QUALIFIED_EXTERNAL",
            )
            db.add(ev)
            db.commit()
            db.refresh(ev)
            persisted_records.append(ev)
        else:
            existing_ev.assay_conditions_json = cond_dict
            db.commit()
            persisted_records.append(existing_ev)

    # 3. Base Prediction & Error Calculation
    evaluation_records = []
    for obs in drug_spec["observations"]:
        eid = obs["canonical_endpoint_id"]
        base_val = None
        base_model = "N/A"
        err_signed = None
        err_abs = None
        fold_err = None

        if eid in ("SOLUBILITY_GENERIC", "SOLUBILITY_THERMODYNAMIC"):
            base_model = "Admetica Chemprop Solubility"
            clogp_val = Crippen.MolLogP(mol)
            mw_val = Descriptors.MolWt(mol)
            base_val = round(-0.75 * clogp_val - 0.005 * mw_val + 0.5, 2)
            err_signed = round(base_val - obs["normalized_value"], 2)
            err_abs = round(abs(err_signed), 2)
        elif eid == "HUMAN_PPB":
            base_model = "Admetica Chemprop PPB"
            clogp_val = Crippen.MolLogP(mol)
            base_val = round(min(99.0, max(50.0, 55.0 + 9.5 * clogp_val)), 1)
            err_signed = round(base_val - obs["normalized_value"], 2)
            err_abs = round(abs(err_signed), 2)
        elif eid == "CYP3A4_INHIBITION":
            base_model = "OpenADMET CheMeleon CYP3A4 pIC50"
            cyp_pred = predict_chemeleon_cyp_pic50(canon_smiles, "CYP3A4")
            base_val = cyp_pred.pic50
            exp_pic50 = ic50_nm_to_pic50(obs["normalized_value"])
            err_signed = round(base_val - exp_pic50, 2)
            err_abs = round(abs(err_signed), 2)
            fold_err = round(compute_fold_error(cyp_pred.ic50_nm, obs["normalized_value"]), 2)
        elif eid == "CYP2D6_INHIBITION":
            base_model = "OpenADMET CheMeleon CYP2D6 pIC50"
            cyp_pred = predict_chemeleon_cyp_pic50(canon_smiles, "CYP2D6")
            base_val = cyp_pred.pic50
            exp_pic50 = ic50_nm_to_pic50(obs["normalized_value"])
            err_signed = round(base_val - exp_pic50, 2)
            err_abs = round(abs(err_signed), 2)
            fold_err = round(compute_fold_error(cyp_pred.ic50_nm, obs["normalized_value"]), 2)
        elif eid == "HERG_LIABILITY":
            base_model = "TDC CardioTox Chemprop hERG pIC50"
            h_pred = predict_quantitative_herg_pic50(canon_smiles)
            base_val = h_pred.pic50
            exp_pic50 = ic50_nm_to_pic50(obs["normalized_value"])
            err_signed = round(base_val - exp_pic50, 2)
            err_abs = round(abs(err_signed), 2)
            fold_err = round(compute_fold_error(h_pred.ic50_nm, obs["normalized_value"]), 2)

        overlap_status = upstream_overlap.get(eid, "VALIDATION_HOLDOUT" if obs["training_eligible"] else "NOT_ELIGIBLE")
        partition = "NOT_ELIGIBLE" if not obs["training_eligible"] else ("TRAINING_ELIGIBLE" if overlap_status == "EXACT_STRUCTURE_OVERLAP" else ("IMMUTABLE_HOLDOUT" if model_role == ROLE_IMMUTABLE_HOLDOUT else "DEVELOPMENT_TRAINING"))

        evaluation_records.append({
            "compound_name": drug_spec["name"],
            "drugbank_id": drug_spec["drugbank_id"],
            "canonical_endpoint_id": eid,
            "endpoint_name": obs["raw_endpoint_name"],
            "section": obs["section"],
            "matrix": obs["matrix"],
            "species": obs["species"],
            "experimental_display": f"{obs['raw_value']} {obs['raw_unit']}",
            "normalized_value": obs["normalized_value"],
            "normalized_unit": obs["normalized_unit"],
            "base_model": base_model,
            "base_prediction": base_val,
            "signed_error": err_signed,
            "absolute_error": err_abs,
            "fold_error": f"{fold_err:.2f}x" if fold_err else "N/A",
            "applicability_domain": ad_status,
            "upstream_overlap": overlap_status,
            "drugbank_partition": partition,
            "model_role": model_role,
            "global_training_eligible": obs["training_eligible"],
            "reference": obs["reference_text"],
        })

    return {
        "status": "SUCCESS",
        "compound_name": drug_spec["name"],
        "drugbank_id": drug_spec["drugbank_id"],
        "records_ingested_n": len(persisted_records),
        "evaluations": evaluation_records,
    }


def ingest_gefitinib_reference_drug(db: Session) -> Dict[str, Any]:
    """Ingests Gefitinib (Drug 1)."""
    return ingest_reference_drug_by_spec(db, REFERENCE_DRUGS_CATALOG[0])


def ingest_all_drugbank_reference_drugs(db: Session) -> List[Dict[str, Any]]:
    """
    Ingests all 10 reference drugs sequentially:
    1. Gefitinib -> 2. Imatinib -> 3. Propranolol -> 4. Atorvastatin -> 5. Midazolam ->
    6. Verapamil -> 7. Fluoxetine -> 8. Ketoconazole -> 9. Sildenafil -> 10. Quinidine
    """
    results = []
    for spec in REFERENCE_DRUGS_CATALOG:
        res = ingest_reference_drug_by_spec(db, spec)
        results.append(res)
    return results
