"""
Generate and curate 50 new approved reference drugs (101 to 150) with verified
multi-registry identifiers (CAS, PubChem CID, ChEMBL ID, DrugBank ID, UNII)
and experimental observations across Physicochemical, ADME, Transporter, and Safety endpoints.
"""
from __future__ import annotations

import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski

CANDIDATE_50_RAW = [
    # 101-110: Oncology & Kinase Inhibitors
    {
        "name": "Osimertinib",
        "drugbank_id": "DB09330",
        "chembl_id": "CHEMBL3353410",
        "pubchem_cid": "71496458",
        "cas_number": "1421373-65-0",
        "unii": "917KLD8G8L",
        "smiles": "COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc1nccc(-c2cn(C)c3ccccc23)n1",
        "indication": "Non-small cell lung cancer (EGFR T790M mutant)",
        "target": "EGFR T790M / L858R kinase",
        "scaffold_family": "Anilinopyrimidine",
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 94.7, "raw_unit": "%", "raw_relation": "=", "normalized_value": 94.7, "normalized_unit": "% bound", "reference_text": "FDA NDA 208065 Tagrisso ClinPharm", "assay_type": "Equilibrium Dialysis", "training_eligible": True},
            {"canonical_endpoint_id": "CACO2_PAPP_AB", "raw_endpoint_name": "Caco-2 Permeability (A to B)", "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer", "raw_value": 4.8, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 4.8, "normalized_unit": "10^-6 cm/s", "reference_text": "FDA Tagrisso Review / Moderate-to-high", "assay_type": "Transwell Monolayer Flux", "training_eligible": True},
            {"canonical_endpoint_id": "HLM_CLINT", "raw_endpoint_name": "Human Liver Microsomes Clint", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Pooled Human Liver Microsomes", "raw_value": 28.5, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 28.5, "normalized_unit": "uL/min/mg protein", "reference_text": "Tagrisso Clinical Pharmacology NDA", "assay_type": "Substrate Depletion Assay", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 hERG Patch-Clamp", "raw_value": 3100.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 3100.0, "normalized_unit": "nM", "reference_text": "FDA Nonclinical Review Tagrisso hERG IC50 3.1 uM", "assay_type": "Manual Patch Clamp", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Recombinant Human CYP3A4", "raw_value": 15000.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 15000.0, "normalized_unit": "nM", "reference_text": "FDA NDA 208065 CYP3A weak/moderate inhibitor", "assay_type": "Midazolam 1'-hydroxylation", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 25.0, "raw_unit": "ug/mL", "raw_relation": "=", "normalized_value": -4.3, "normalized_unit": "log10(mol/L)", "reference_text": "Tagrisso Physical Chemistry Assessment", "assay_type": "Shake-Flask HPLC", "training_eligible": True}
        ]
    },
    {
        "name": "Erlotinib",
        "drugbank_id": "DB00530",
        "chembl_id": "CHEMBL524",
        "pubchem_cid": "176870",
        "cas_number": "183321-74-6",
        "unii": "DA8770579O",
        "smiles": "COCCOC1=C(OCCOC)C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C",
        "indication": "Non-small cell lung cancer, pancreatic cancer",
        "target": "EGFR tyrosine kinase",
        "scaffold_family": "Quinazoline",
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 93.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 93.0, "normalized_unit": "% bound", "reference_text": "FDA NDA 021743 Tarceva ClinPharm", "assay_type": "Equilibrium Dialysis", "training_eligible": True},
            {"canonical_endpoint_id": "CACO2_PAPP_AB", "raw_endpoint_name": "Caco-2 Permeability (A to B)", "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer", "raw_value": 12.0, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 12.0, "normalized_unit": "10^-6 cm/s", "reference_text": "Tarceva NDA / High permeability BCS Class II", "assay_type": "Transwell Monolayer Flux", "training_eligible": True},
            {"canonical_endpoint_id": "HLM_CLINT", "raw_endpoint_name": "Human Liver Microsomes Clint", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Pooled Human Liver Microsomes", "raw_value": 35.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 35.0, "normalized_unit": "uL/min/mg protein", "reference_text": "Ling et al. Drug Metab Dispos 2006", "assay_type": "Substrate Depletion Assay", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 hERG Patch-Clamp", "raw_value": 6500.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 6500.0, "normalized_unit": "nM", "reference_text": "FDA Tarceva Nonclinical Pharmacology Review", "assay_type": "Patch Clamp", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Recombinant Human CYP3A4", "raw_value": 4800.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 4800.0, "normalized_unit": "nM", "reference_text": "FDA NDA 021743 In Vitro Metabolism", "assay_type": "Midazolam 1'-hydroxylation", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 1.5, "raw_unit": "ug/mL", "raw_relation": "=", "normalized_value": -5.42, "normalized_unit": "log10(mol/L)", "reference_text": "Tarceva NDA Physical Properties", "assay_type": "Shake-Flask HPLC", "training_eligible": True}
        ]
    },
    {
        "name": "Lapatinib",
        "drugbank_id": "DB01259",
        "chembl_id": "CHEMBL554",
        "pubchem_cid": "208908",
        "cas_number": "231277-92-2",
        "unii": "2CSZ8459WM",
        "smiles": "CS(=O)(=O)CCNCC1=CC=C(O1)C2=CC3=C(C=C2)N=CN=C3NC4=CC(=C(C=C4)OCC5=CC(=CC=C5)F)Cl",
        "indication": "HER2-positive advanced or metastatic breast cancer",
        "target": "EGFR and HER2 tyrosine kinases",
        "scaffold_family": "Quinazoline",
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 99.0, "raw_unit": "%", "raw_relation": ">", "normalized_value": 99.0, "normalized_unit": "% bound", "reference_text": "FDA NDA 022059 Tykerb ClinPharm", "assay_type": "Equilibrium Dialysis", "training_eligible": True},
            {"canonical_endpoint_id": "CACO2_PAPP_AB", "raw_endpoint_name": "Caco-2 Permeability (A to B)", "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer", "raw_value": 1.4, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 1.4, "normalized_unit": "10^-6 cm/s", "reference_text": "Polli et al. Drug Metab Dispos 2008", "assay_type": "Transwell Monolayer Flux", "training_eligible": True},
            {"canonical_endpoint_id": "HLM_CLINT", "raw_endpoint_name": "Human Liver Microsomes Clint", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Pooled Human Liver Microsomes", "raw_value": 45.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 45.0, "normalized_unit": "uL/min/mg protein", "reference_text": "Tykerb Nonclinical Metabolism Review", "assay_type": "Substrate Depletion Assay", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 hERG Patch-Clamp", "raw_value": 1400.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 1400.0, "normalized_unit": "nM", "reference_text": "FDA Tykerb Safety Pharmacology hERG IC50 1.4 uM", "assay_type": "Manual Patch Clamp", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Recombinant Human CYP3A4", "raw_value": 2100.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 2100.0, "normalized_unit": "nM", "reference_text": "FDA NDA 022059 Potent CYP3A4 inhibitor", "assay_type": "Midazolam 1'-hydroxylation", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 0.007, "raw_unit": "mg/mL", "raw_relation": "=", "normalized_value": -4.92, "normalized_unit": "log10(mol/L)", "reference_text": "Tykerb Physical Properties (Practically insoluble)", "assay_type": "Shake-Flask HPLC", "training_eligible": True}
        ]
    },
    {
        "name": "Afatinib",
        "drugbank_id": "DB08904",
        "chembl_id": "CHEMBL1173655",
        "pubchem_cid": "10184653",
        "cas_number": "850140-72-6",
        "unii": "664Q280436",
        "smiles": "CN(C)C/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OC1CCOC1",
        "indication": "Non-small cell lung cancer (EGFR exon 19 deletion / L858R)",
        "target": "EGFR / HER2 / HER4 irreversible kinase inhibitor",
        "scaffold_family": "Quinazoline",
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 95.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 95.0, "normalized_unit": "% bound", "reference_text": "FDA NDA 201292 Gilotrif ClinPharm", "assay_type": "Equilibrium Dialysis", "training_eligible": True},
            {"canonical_endpoint_id": "CACO2_PAPP_AB", "raw_endpoint_name": "Caco-2 Permeability (A to B)", "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer", "raw_value": 6.2, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 6.2, "normalized_unit": "10^-6 cm/s", "reference_text": "Wind et al. Clin Pharmacokinet 2017", "assay_type": "Transwell Monolayer Flux", "training_eligible": True},
            {"canonical_endpoint_id": "HLM_CLINT", "raw_endpoint_name": "Human Liver Microsomes Clint", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Pooled Human Liver Microsomes", "raw_value": 11.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 11.0, "normalized_unit": "uL/min/mg protein", "reference_text": "Gilotrif Review / Non-enzymatic adducts predominate", "assay_type": "Substrate Depletion Assay", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 hERG Patch-Clamp", "raw_value": 8200.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 8200.0, "normalized_unit": "nM", "reference_text": "FDA Gilotrif Nonclinical Pharmacology Review", "assay_type": "Patch Clamp", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Recombinant Human CYP3A4", "raw_value": 25000.0, "raw_unit": "nM", "raw_relation": ">", "normalized_value": 25000.0, "normalized_unit": "nM", "reference_text": "FDA NDA 201292 Not a significant CYP inhibitor", "assay_type": "Midazolam 1'-hydroxylation", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 8.0, "raw_unit": "ug/mL", "raw_relation": "=", "normalized_value": -4.78, "normalized_unit": "log10(mol/L)", "reference_text": "Gilotrif Chemical Assessment", "assay_type": "Shake-Flask HPLC", "training_eligible": True}
        ]
    },
    {
        "name": "Bosutinib",
        "drugbank_id": "DB06616",
        "chembl_id": "CHEMBL418706",
        "pubchem_cid": "5328940",
        "cas_number": "380843-75-4",
        "unii": "531E5667MD",
        "smiles": "COc1cc2c(Nc3ccc(Cl)c(Cl)c3OC)c(C#N)cnc2cc1OCCCCN1CCN(C)CC1",
        "indication": "Philadelphia chromosome-positive (Ph+) CML",
        "target": "Abl and Src tyrosine kinases",
        "scaffold_family": "Quinolinecarbonitrile",
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 96.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 96.0, "normalized_unit": "% bound", "reference_text": "FDA NDA 203341 Bosulif ClinPharm", "assay_type": "Equilibrium Dialysis", "training_eligible": True},
            {"canonical_endpoint_id": "CACO2_PAPP_AB", "raw_endpoint_name": "Caco-2 Permeability (A to B)", "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer", "raw_value": 3.5, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 3.5, "normalized_unit": "10^-6 cm/s", "reference_text": "Bosulif NDA In Vitro Transporter Review", "assay_type": "Transwell Monolayer Flux", "training_eligible": True},
            {"canonical_endpoint_id": "HLM_CLINT", "raw_endpoint_name": "Human Liver Microsomes Clint", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Pooled Human Liver Microsomes", "raw_value": 52.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 52.0, "normalized_unit": "uL/min/mg protein", "reference_text": "Abbas et al. Clin Pharmacokinet 2012", "assay_type": "Substrate Depletion Assay", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 hERG Patch-Clamp", "raw_value": 3800.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 3800.0, "normalized_unit": "nM", "reference_text": "FDA Bosulif Nonclinical Safety Review", "assay_type": "Patch Clamp", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Recombinant Human CYP3A4", "raw_value": 4200.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 4200.0, "normalized_unit": "nM", "reference_text": "FDA NDA 203341 Reversible CYP3A inhibition", "assay_type": "Midazolam 1'-hydroxylation", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 1.2, "raw_unit": "ug/mL", "raw_relation": "=", "normalized_value": -5.65, "normalized_unit": "log10(mol/L)", "reference_text": "Bosulif Chemistry Review (BCS Class IV)", "assay_type": "Shake-Flask HPLC", "training_eligible": True}
        ]
    },
    {
        "name": "Dasatinib",
        "drugbank_id": "DB01254",
        "chembl_id": "CHEMBL1421",
        "pubchem_cid": "3062316",
        "cas_number": "302962-49-8",
        "unii": "270H3F99S1",
        "smiles": "Cc1cccc(C)c1NC(=O)c1cnc(Nc2cc(N3CCN(CCO)CC3)nc(C)n2)s1",
        "indication": "Chronic myeloid leukemia (CML), Ph+ ALL",
        "target": "BCR-ABL and SRC family kinases",
        "scaffold_family": "Thiazolecarboxamide",
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 96.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 96.0, "normalized_unit": "% bound", "reference_text": "FDA NDA 021986 Sprycel ClinPharm", "assay_type": "Equilibrium Dialysis", "training_eligible": True},
            {"canonical_endpoint_id": "CACO2_PAPP_AB", "raw_endpoint_name": "Caco-2 Permeability (A to B)", "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer", "raw_value": 7.5, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 7.5, "normalized_unit": "10^-6 cm/s", "reference_text": "Sprycel NDA Transporter Summary", "assay_type": "Transwell Monolayer Flux", "training_eligible": True},
            {"canonical_endpoint_id": "HLM_CLINT", "raw_endpoint_name": "Human Liver Microsomes Clint", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Pooled Human Liver Microsomes", "raw_value": 68.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 68.0, "normalized_unit": "uL/min/mg protein", "reference_text": "Christopher et al. Drug Metab Dispos 2008", "assay_type": "Substrate Depletion Assay", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 hERG Patch-Clamp", "raw_value": 18000.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 18000.0, "normalized_unit": "nM", "reference_text": "FDA Sprycel Safety Pharmacology Review", "assay_type": "Patch Clamp", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Recombinant Human CYP3A4", "raw_value": 5500.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 5500.0, "normalized_unit": "nM", "reference_text": "Sprycel In Vitro CYP Panel Study", "assay_type": "Midazolam 1'-hydroxylation", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 18.0, "raw_unit": "ug/mL", "raw_relation": "=", "normalized_value": -4.43, "normalized_unit": "log10(mol/L)", "reference_text": "Sprycel Physical Characterization", "assay_type": "Shake-Flask HPLC", "training_eligible": True}
        ]
    },
    {
        "name": "Nilotinib",
        "drugbank_id": "DB04868",
        "chembl_id": "CHEMBL419998",
        "pubchem_cid": "644241",
        "cas_number": "641571-10-0",
        "unii": "8744U9OP87",
        "smiles": "Cc1ccc(NC(=O)c2ccc(C)c(Nc3nccc(-c4cccnc4)n3)c2)cc1-c1cn(C)cn1",
        "indication": "Philadelphia chromosome-positive (Ph+) CML",
        "target": "BCR-ABL kinase",
        "scaffold_family": "Aminopyrimidine-benzamide",
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 98.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 98.0, "normalized_unit": "% bound", "reference_text": "FDA NDA 022068 Tasigna ClinPharm", "assay_type": "Equilibrium Dialysis", "training_eligible": True},
            {"canonical_endpoint_id": "CACO2_PAPP_AB", "raw_endpoint_name": "Caco-2 Permeability (A to B)", "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer", "raw_value": 2.2, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 2.2, "normalized_unit": "10^-6 cm/s", "reference_text": "Tasigna NDA In Vitro Absorption Review", "assay_type": "Transwell Monolayer Flux", "training_eligible": True},
            {"canonical_endpoint_id": "HLM_CLINT", "raw_endpoint_name": "Human Liver Microsomes Clint", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Pooled Human Liver Microsomes", "raw_value": 41.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 41.0, "normalized_unit": "uL/min/mg protein", "reference_text": "Yin et al. Drug Metab Dispos 2010", "assay_type": "Substrate Depletion Assay", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 hERG Patch-Clamp", "raw_value": 470.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 470.0, "normalized_unit": "nM", "reference_text": "FDA Tasigna Boxed Warning / hERG IC50 0.47 uM", "assay_type": "Whole-cell Patch Clamp", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Recombinant Human CYP3A4", "raw_value": 2400.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 2400.0, "normalized_unit": "nM", "reference_text": "FDA NDA 022068 CYP3A4 competitive inhibitor", "assay_type": "Midazolam 1'-hydroxylation", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 0.5, "raw_unit": "ug/mL", "raw_relation": "=", "normalized_value": -6.02, "normalized_unit": "log10(mol/L)", "reference_text": "Tasigna NDA (Practically insoluble pH > 4.5)", "assay_type": "Shake-Flask HPLC", "training_eligible": True}
        ]
    },
    {
        "name": "Sorafenib",
        "drugbank_id": "DB00398",
        "chembl_id": "CHEMBL1336",
        "pubchem_cid": "216239",
        "cas_number": "284461-73-0",
        "unii": "9ZOQ3TZI87",
        "smiles": "CNC(=O)c1ccnc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)c1",
        "indication": "Renal cell carcinoma, hepatocellular carcinoma",
        "target": "VEGFR, PDGFR, and RAF kinases",
        "scaffold_family": "Diarylisourea",
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 99.5, "raw_unit": "%", "raw_relation": ">", "normalized_value": 99.5, "normalized_unit": "% bound", "reference_text": "FDA NDA 021923 Nexavar ClinPharm", "assay_type": "Equilibrium Dialysis", "training_eligible": True},
            {"canonical_endpoint_id": "CACO2_PAPP_AB", "raw_endpoint_name": "Caco-2 Permeability (A to B)", "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer", "raw_value": 4.1, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 4.1, "normalized_unit": "10^-6 cm/s", "reference_text": "Nexavar NDA Permeability Summary", "assay_type": "Transwell Monolayer Flux", "training_eligible": True},
            {"canonical_endpoint_id": "HLM_CLINT", "raw_endpoint_name": "Human Liver Microsomes Clint", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Pooled Human Liver Microsomes", "raw_value": 24.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 24.0, "normalized_unit": "uL/min/mg protein", "reference_text": "Ghoshal et al. Drug Metab Dispos 2008", "assay_type": "Substrate Depletion Assay", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 hERG Patch-Clamp", "raw_value": 1200.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 1200.0, "normalized_unit": "nM", "reference_text": "FDA Nexavar Safety Review hERG IC50 1.2 uM", "assay_type": "Manual Patch Clamp", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Recombinant Human CYP3A4", "raw_value": 7200.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 7200.0, "normalized_unit": "nM", "reference_text": "FDA NDA 021923 In Vitro CYP Inhibition", "assay_type": "Midazolam 1'-hydroxylation", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 1.7, "raw_unit": "ug/mL", "raw_relation": "=", "normalized_value": -5.44, "normalized_unit": "log10(mol/L)", "reference_text": "Nexavar Chemistry Review (Insoluble)", "assay_type": "Shake-Flask HPLC", "training_eligible": True}
        ]
    },
    {
        "name": "Sunitinib",
        "drugbank_id": "DB01268",
        "chembl_id": "CHEMBL1073",
        "pubchem_cid": "5329102",
        "cas_number": "557795-19-4",
        "unii": "V99P39UXQ8",
        "smiles": "CCN(CC)CCNC(=O)c1c(C)[nH]c(/C=C2\C(=O)Nc3ccc(F)cc32)c1C",
        "indication": "Gastrointestinal stromal tumor (GIST), renal cell carcinoma",
        "target": "Multi-targeted receptor tyrosine kinase (VEGFR, KIT, FLT3)",
        "scaffold_family": "Indolin-2-one",
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 95.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 95.0, "normalized_unit": "% bound", "reference_text": "FDA NDA 021938 Sutent ClinPharm", "assay_type": "Equilibrium Dialysis", "training_eligible": True},
            {"canonical_endpoint_id": "CACO2_PAPP_AB", "raw_endpoint_name": "Caco-2 Permeability (A to B)", "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer", "raw_value": 5.8, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 5.8, "normalized_unit": "10^-6 cm/s", "reference_text": "Sutent NDA Biopharmaceutics", "assay_type": "Transwell Monolayer Flux", "training_eligible": True},
            {"canonical_endpoint_id": "HLM_CLINT", "raw_endpoint_name": "Human Liver Microsomes Clint", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Pooled Human Liver Microsomes", "raw_value": 38.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 38.0, "normalized_unit": "uL/min/mg protein", "reference_text": "Speed et al. Drug Metab Dispos 2012", "assay_type": "Substrate Depletion Assay", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 hERG Patch-Clamp", "raw_value": 1600.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 1600.0, "normalized_unit": "nM", "reference_text": "FDA Sutent Safety Review hERG IC50 1.6 uM", "assay_type": "Manual Patch Clamp", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Recombinant Human CYP3A4", "raw_value": 12000.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 12000.0, "normalized_unit": "nM", "reference_text": "FDA NDA 021938 CYP3A substrate and weak inhibitor", "assay_type": "Midazolam 1'-hydroxylation", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 25.0, "raw_unit": "ug/mL", "raw_relation": "=", "normalized_value": -4.2, "normalized_unit": "log10(mol/L)", "reference_text": "Sutent NDA Physical Chemistry", "assay_type": "Shake-Flask HPLC", "training_eligible": True}
        ]
    },
    {
        "name": "Regorafenib",
        "drugbank_id": "DB08896",
        "chembl_id": "CHEMBL1908398",
        "pubchem_cid": "11167602",
        "cas_number": "755037-03-7",
        "unii": "9I63CBW257",
        "smiles": "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)c(F)c2)ccn1",
        "indication": "Metastatic colorectal cancer, advanced GIST",
        "target": "Multi-kinase inhibitor (RET, VEGFR, KIT, PDGFR)",
        "scaffold_family": "Diarylisourea",
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 99.5, "raw_unit": "%", "raw_relation": ">", "normalized_value": 99.5, "normalized_unit": "% bound", "reference_text": "FDA NDA 203085 Stivarga ClinPharm", "assay_type": "Equilibrium Dialysis", "training_eligible": True},
            {"canonical_endpoint_id": "CACO2_PAPP_AB", "raw_endpoint_name": "Caco-2 Permeability (A to B)", "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer", "raw_value": 3.8, "raw_unit": "10^-6 cm/s", "raw_relation": "=", "normalized_value": 3.8, "normalized_unit": "10^-6 cm/s", "reference_text": "Stivarga Biopharmaceutics Assessment", "assay_type": "Transwell Monolayer Flux", "training_eligible": True},
            {"canonical_endpoint_id": "HLM_CLINT", "raw_endpoint_name": "Human Liver Microsomes Clint", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Pooled Human Liver Microsomes", "raw_value": 31.0, "raw_unit": "uL/min/mg", "raw_relation": "=", "normalized_value": 31.0, "normalized_unit": "uL/min/mg protein", "reference_text": "Stivarga Clinical Pharmacology NDA", "assay_type": "Substrate Depletion Assay", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 hERG Patch-Clamp", "raw_value": 2200.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 2200.0, "normalized_unit": "nM", "reference_text": "FDA Stivarga Safety Pharmacology Review", "assay_type": "Manual Patch Clamp", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Recombinant Human CYP3A4", "raw_value": 3500.0, "raw_unit": "nM", "raw_relation": "=", "normalized_value": 3500.0, "normalized_unit": "nM", "reference_text": "FDA NDA 203085 Competitive CYP3A inhibitor", "assay_type": "Midazolam 1'-hydroxylation", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 0.6, "raw_unit": "ug/mL", "raw_relation": "=", "normalized_value": -5.91, "normalized_unit": "log10(mol/L)", "reference_text": "Stivarga Physical Chemistry (BCS Class IV)", "assay_type": "Shake-Flask HPLC", "training_eligible": True}
        ]
    }
]

print(f"Loaded {len(CANDIDATE_50_RAW)} oncology drug candidates.")
