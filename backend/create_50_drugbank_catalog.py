"""
Script to create the complete 50 approved reference drugs catalog in DrugBank.
Expands from 40 to 50 drugs for Drug-OPT Global Engine v3.1.
"""
import json
from pathlib import Path

# Load original 40 drugs
with open("/home/xavier/chem/drug-opt-platform/backend/reference_drugs_40.json") as f:
    ORIGINAL_40 = json.load(f)

# Update Drugs 36-40 (Atenolol, Caffeine, Ibuprofen, Lorcaserin, Rosuvastatin)
# to FINAL_TEST_COHORT_2_CONSUMED because they were evaluated in v3.0 locked final test.
for d in ORIGINAL_40:
    if d["name"] in ("Atenolol", "Caffeine", "Ibuprofen", "Lorcaserin", "Rosuvastatin"):
        d["model_role"] = "FINAL_TEST_COHORT_2_CONSUMED"
        d["cohort"] = "FINAL_TEST_COHORT_2_CONSUMED"

# Define 10 new approved reference drugs (41-50)
NEW_10_DRUGS = [
    # 41: Amlodipine (DEVELOPMENT_TRAINING)
    {
        "name": "Amlodipine", "cas_number": "88150-42-9", "drugbank_id": "DB00381", "pubchem_cid": "2162", "chembl_id": "CHEMBL1487", "unii": "1J444QC288",
        "smiles": "CCOC(=O)C1=C(COCCN)NC(C)=C(C(=O)OC)C1c1ccccc1Cl", "indication": "Hypertension / Coronary Artery Disease", "target": "Cav1.2 Calcium Channel / CYP3A4 Substrate", "scaffold_family": "Dihydropyridine-ethanolamine",
        "model_role": "DEVELOPMENT_TRAINING", "cohort": "DEV_TRAINING",
        "upstream_overlap": {"SOLUBILITY_GENERIC": "NOVEL_IN_DOMAIN", "SOLUBILITY_THERMODYNAMIC": "NOVEL_IN_DOMAIN", "HUMAN_PPB": "NOVEL_IN_DOMAIN", "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT", "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT", "HERG_LIABILITY": "VALIDATION_HOLDOUT"},
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 98.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 98.0, "normalized_unit": "% bound", "reference_text": "FDA Norvasc NDA 019787 Clinical Pharmacology Review", "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp", "raw_value": 10.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 10000.0, "normalized_unit": "nM", "reference_text": "Redfern et al. / ChEMBL1487", "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP3A4", "raw_value": 8.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 8000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Guidance Review", "assay_type": "Recombinant human CYP3A4 Midazolam 1'-hydroxylation assay", "training_eligible": True},
            {"canonical_endpoint_id": "CYP2D6_INHIBITION", "raw_endpoint_name": "CYP2D6 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP2D6", "raw_value": 50.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 50000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Review", "assay_type": "Recombinant human CYP2D6 Dextromethorphan assay", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 75.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -4.12, "normalized_unit": "log10(mol/L)", "reference_text": "DrugCentral / Norvasc NDA 019787", "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True},
            {"canonical_endpoint_id": "CACO2_PERMEABILITY", "raw_endpoint_name": "Caco-2 Apparent Permeability", "section": "ADMET", "species": "Homo sapiens", "matrix": "Caco-2 Monolayer pH 7.4/7.4", "raw_value": 1.5e-6, "raw_unit": "cm/s", "raw_relation": "=", "normalized_value": -5.82, "normalized_unit": "log10(cm/s)", "reference_text": "FDA Biopharmaceutics Review / Norvasc", "assay_type": "Bidirectional Caco-2 Cell Permeability (A to B)", "training_eligible": True}
        ]
    },
    # 42: Losartan (DEVELOPMENT_TRAINING)
    {
        "name": "Losartan", "cas_number": "114798-26-4", "drugbank_id": "DB00678", "pubchem_cid": "3961", "chembl_id": "CHEMBL1565", "unii": "JMS50WPO89",
        "smiles": "CCCC1=NC(Cl)=C(CO)N1Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1", "indication": "Hypertension / Diabetic Nephropathy", "target": "Angiotensin II Type 1 (AT1) Receptor / CYP2C9 & CYP3A4 Substrate", "scaffold_family": "Imidazole-biphenyl-tetrazole",
        "model_role": "DEVELOPMENT_TRAINING", "cohort": "DEV_TRAINING",
        "upstream_overlap": {"SOLUBILITY_GENERIC": "NOVEL_IN_DOMAIN", "SOLUBILITY_THERMODYNAMIC": "NOVEL_IN_DOMAIN", "HUMAN_PPB": "NOVEL_IN_DOMAIN", "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT", "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT", "HERG_LIABILITY": "VALIDATION_HOLDOUT"},
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 98.7, "raw_unit": "%", "raw_relation": "=", "normalized_value": 98.7, "normalized_unit": "% bound", "reference_text": "FDA Cozaar NDA 020386 Clinical Pharmacology Review", "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp", "raw_value": 25.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 25000.0, "normalized_unit": "nM", "reference_text": "Redfern et al. / ChEMBL1565", "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP3A4", "raw_value": 15.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 15000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Guidance Review", "assay_type": "Recombinant human CYP3A4 Midazolam 1'-hydroxylation assay", "training_eligible": True},
            {"canonical_endpoint_id": "CYP2D6_INHIBITION", "raw_endpoint_name": "CYP2D6 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP2D6", "raw_value": 80.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 80000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Review Cozaar", "assay_type": "Recombinant human CYP2D6 Dextromethorphan assay", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 33.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -4.48, "normalized_unit": "log10(mol/L)", "reference_text": "DrugCentral / Cozaar NDA 020386", "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True},
            {"canonical_endpoint_id": "HLM_INTRINSIC_CLEARANCE", "raw_endpoint_name": "Human Liver Microsomes Intrinsic Clearance", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "Human Liver Microsomes (HLM)", "raw_value": 22.0, "raw_unit": "µL/min/mg", "raw_relation": "=", "normalized_value": 22.0, "normalized_unit": "µL/min/mg", "reference_text": "FDA Clinical Pharmacology Biopharmaceutics Review", "assay_type": "Substrate Depletion (HLM + NADPH)", "training_eligible": True}
        ]
    },
    # 43: Metronidazole (DEVELOPMENT_TRAINING)
    {
        "name": "Metronidazole", "cas_number": "443-48-1", "drugbank_id": "DB00916", "pubchem_cid": "4173", "chembl_id": "CHEMBL128", "unii": "140QMO861O",
        "smiles": "Cc1ncc([N+](=O)[O-])n1CCO", "indication": "Anaerobic Bacterial & Protozoal Infections", "target": "Bacterial DNA / Moderate CYP2C9 Inhibitor", "scaffold_family": "Nitroimidazole",
        "model_role": "DEVELOPMENT_TRAINING", "cohort": "DEV_TRAINING",
        "upstream_overlap": {"SOLUBILITY_GENERIC": "NOVEL_IN_DOMAIN", "SOLUBILITY_THERMODYNAMIC": "NOVEL_IN_DOMAIN", "HUMAN_PPB": "NOVEL_IN_DOMAIN", "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT", "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT", "HERG_LIABILITY": "VALIDATION_HOLDOUT"},
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 20.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 20.0, "normalized_unit": "% bound", "reference_text": "FDA Flagyl NDA 012623 Clinical Pharmacology Review", "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp", "raw_value": 300.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 300000.0, "normalized_unit": "nM", "reference_text": "Redfern et al. / ChEMBL128", "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP3A4", "raw_value": 100.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 100000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Review Flagyl", "assay_type": "Recombinant human CYP3A4 Midazolam 1'-hydroxylation assay", "training_eligible": True},
            {"canonical_endpoint_id": "CYP2D6_INHIBITION", "raw_endpoint_name": "CYP2D6 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP2D6", "raw_value": 100.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 100000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Review Flagyl", "assay_type": "Recombinant human CYP2D6 Dextromethorphan assay", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 55000.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -1.26, "normalized_unit": "log10(mol/L)", "reference_text": "DrugCentral / Flagyl NDA 012623", "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True}
        ]
    },
    # 44: Montelukast (MODEL_SELECTION_VALIDATION)
    {
        "name": "Montelukast", "cas_number": "158966-92-8", "drugbank_id": "DB01147", "pubchem_cid": "5281040", "chembl_id": "CHEMBL1405", "unii": "MHO6E7N6UU",
        "smiles": "CC(C)(O)c1ccccc1CC/C(=C/c1ccc2ccc(Cl)cc2c1)SCCC1(CC(=O)O)CC1", "indication": "Asthma / Allergic Rhinitis", "target": "CysLT1 Receptor / Very High PPB / CYP2C8 & CYP2C9 Substrate", "scaffold_family": "Quinoline-styryl-cyclopropane",
        "model_role": "MODEL_SELECTION_VALIDATION", "cohort": "VALIDATION_COHORT_2",
        "upstream_overlap": {"SOLUBILITY_GENERIC": "NOVEL_IN_DOMAIN", "SOLUBILITY_THERMODYNAMIC": "NOVEL_IN_DOMAIN", "HUMAN_PPB": "NOVEL_IN_DOMAIN", "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT", "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT", "HERG_LIABILITY": "VALIDATION_HOLDOUT"},
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 99.7, "raw_unit": "%", "raw_relation": "=", "normalized_value": 99.7, "normalized_unit": "% bound", "reference_text": "FDA Singulair NDA 020829 Clinical Pharmacology Review", "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp", "raw_value": 1.2, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 1200.0, "normalized_unit": "nM", "reference_text": "Redfern et al. / ChEMBL1405", "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP3A4", "raw_value": 5.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 5000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Guidance Review", "assay_type": "Recombinant human CYP3A4 Midazolam 1'-hydroxylation assay", "training_eligible": True},
            {"canonical_endpoint_id": "CYP2D6_INHIBITION", "raw_endpoint_name": "CYP2D6 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP2D6", "raw_value": 20.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 20000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Review Singulair", "assay_type": "Recombinant human CYP2D6 Dextromethorphan assay", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 0.35, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -6.46, "normalized_unit": "log10(mol/L)", "reference_text": "DrugCentral / Singulair NDA 020829", "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True}
        ]
    },
    # 45: Pantoprazole (MODEL_SELECTION_VALIDATION)
    {
        "name": "Pantoprazole", "cas_number": "102625-70-7", "drugbank_id": "DB00213", "pubchem_cid": "4679", "chembl_id": "CHEMBL1506", "unii": "D8TST4O562",
        "smiles": "COc1cnc(CS(=O)c2nc3ccccc3[nH]2)c(OC)c1OC(F)F", "indication": "GERD / Zollinger-Ellison Syndrome", "target": "H+/K+ ATPase / CYP2C19 & CYP3A4 Substrate", "scaffold_family": "Difluoromethoxy-benzimidazole-sulfinyl",
        "model_role": "MODEL_SELECTION_VALIDATION", "cohort": "VALIDATION_COHORT_2",
        "upstream_overlap": {"SOLUBILITY_GENERIC": "NOVEL_IN_DOMAIN", "SOLUBILITY_THERMODYNAMIC": "NOVEL_IN_DOMAIN", "HUMAN_PPB": "NOVEL_IN_DOMAIN", "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT", "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT", "HERG_LIABILITY": "VALIDATION_HOLDOUT"},
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 98.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 98.0, "normalized_unit": "% bound", "reference_text": "FDA Protonix NDA 020987 Clinical Pharmacology Review", "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp", "raw_value": 18.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 18000.0, "normalized_unit": "nM", "reference_text": "Redfern et al. / ChEMBL1506", "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP3A4", "raw_value": 6.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 6000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Guidance Review", "assay_type": "Recombinant human CYP3A4 Midazolam 1'-hydroxylation assay", "training_eligible": True},
            {"canonical_endpoint_id": "CYP2D6_INHIBITION", "raw_endpoint_name": "CYP2D6 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP2D6", "raw_value": 35.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 35000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Review Protonix", "assay_type": "Recombinant human CYP2D6 Dextromethorphan assay", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 52.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -4.28, "normalized_unit": "log10(mol/L)", "reference_text": "DrugCentral / Protonix NDA 020987", "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True}
        ]
    },
    # 46: Raloxifene (LOCKED_FINAL_TEST_COHORT_3)
    {
        "name": "Raloxifene", "cas_number": "84449-90-1", "drugbank_id": "DB00479", "pubchem_cid": "5035", "chembl_id": "CHEMBL53", "unii": "4F86W0707Z",
        "smiles": "O=C(c1ccc(OCCN2CCCCC2)cc1)c1c(-c2ccc(O)cc2)sc2cc(O)ccc12", "indication": "Osteoporosis Prophylaxis / Invasive Breast Cancer Risk Reduction", "target": "Estrogen Receptor Alpha/Beta / SERM", "scaffold_family": "Benzothiophene-ketone",
        "model_role": "LOCKED_FINAL_TEST_COHORT_3", "cohort": "LOCKED_FINAL_TEST_COHORT_3",
        "upstream_overlap": {"SOLUBILITY_GENERIC": "NOVEL_IN_DOMAIN", "SOLUBILITY_THERMODYNAMIC": "NOVEL_IN_DOMAIN", "HUMAN_PPB": "NOVEL_IN_DOMAIN", "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT", "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT", "HERG_LIABILITY": "VALIDATION_HOLDOUT"},
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 95.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 95.0, "normalized_unit": "% bound", "reference_text": "FDA Evista NDA 020815 Clinical Pharmacology Review", "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp", "raw_value": 3.5, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 3500.0, "normalized_unit": "nM", "reference_text": "Redfern et al. / ChEMBL53", "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP3A4", "raw_value": 12.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 12000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Guidance Review", "assay_type": "Recombinant human CYP3A4 Midazolam 1'-hydroxylation assay", "training_eligible": True},
            {"canonical_endpoint_id": "CYP2D6_INHIBITION", "raw_endpoint_name": "CYP2D6 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP2D6", "raw_value": 40.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 40000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Review Evista", "assay_type": "Recombinant human CYP2D6 Dextromethorphan assay", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 0.85, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -6.07, "normalized_unit": "log10(mol/L)", "reference_text": "DrugCentral / Evista NDA 020815", "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True}
        ]
    },
    # 47: Tamoxifen (LOCKED_FINAL_TEST_COHORT_3)
    {
        "name": "Tamoxifen", "cas_number": "10540-29-1", "drugbank_id": "DB00675", "pubchem_cid": "5376", "chembl_id": "CHEMBL83", "unii": "094ZI81Y45",
        "smiles": "CCC(=C(c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1", "indication": "Breast Cancer / Ductal Carcinoma In Situ", "target": "Estrogen Receptor / Potent CYP2D6 Substrate & Inhibitor", "scaffold_family": "Triphenylethylene",
        "model_role": "LOCKED_FINAL_TEST_COHORT_3", "cohort": "LOCKED_FINAL_TEST_COHORT_3",
        "upstream_overlap": {"SOLUBILITY_GENERIC": "NOVEL_IN_DOMAIN", "SOLUBILITY_THERMODYNAMIC": "NOVEL_IN_DOMAIN", "HUMAN_PPB": "NOVEL_IN_DOMAIN", "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT", "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT", "HERG_LIABILITY": "VALIDATION_HOLDOUT"},
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 99.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 99.0, "normalized_unit": "% bound", "reference_text": "FDA Nolvadex NDA 017970 Clinical Pharmacology Review", "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp", "raw_value": 0.90, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 900.0, "normalized_unit": "nM", "reference_text": "Redfern et al. / ChEMBL83", "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP3A4", "raw_value": 3.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 3000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Guidance Review", "assay_type": "Recombinant human CYP3A4 Midazolam 1'-hydroxylation assay", "training_eligible": True},
            {"canonical_endpoint_id": "CYP2D6_INHIBITION", "raw_endpoint_name": "CYP2D6 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP2D6", "raw_value": 8.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 8000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Review Nolvadex", "assay_type": "Recombinant human CYP2D6 Dextromethorphan assay", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 0.40, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -6.40, "normalized_unit": "log10(mol/L)", "reference_text": "DrugCentral / Nolvadex NDA 017970", "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True}
        ]
    },
    # 48: Theophylline (LOCKED_FINAL_TEST_COHORT_3)
    {
        "name": "Theophylline", "cas_number": "58-55-9", "drugbank_id": "DB00277", "pubchem_cid": "2153", "chembl_id": "CHEMBL276", "unii": "C137DTR5RG",
        "smiles": "Cn1c(=O)c2[nH]cnc2n(C)c1=O", "indication": "Asthma / COPD Bronchospasm", "target": "Phosphodiesterase / Adenosine Receptor / FDA Standard CYP1A2 Probe", "scaffold_family": "Dimethylxanthine",
        "model_role": "LOCKED_FINAL_TEST_COHORT_3", "cohort": "LOCKED_FINAL_TEST_COHORT_3",
        "upstream_overlap": {"SOLUBILITY_GENERIC": "NOVEL_IN_DOMAIN", "SOLUBILITY_THERMODYNAMIC": "NOVEL_IN_DOMAIN", "HUMAN_PPB": "NOVEL_IN_DOMAIN", "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT", "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT", "HERG_LIABILITY": "VALIDATION_HOLDOUT"},
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 40.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 40.0, "normalized_unit": "% bound", "reference_text": "FDA Theo-Dur NDA 085078 Clinical Pharmacology Review", "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp", "raw_value": 250.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 250000.0, "normalized_unit": "nM", "reference_text": "Redfern et al. / ChEMBL276", "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP3A4", "raw_value": 80.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 80000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Review", "assay_type": "Recombinant human CYP3A4 Midazolam 1'-hydroxylation assay", "training_eligible": True},
            {"canonical_endpoint_id": "CYP2D6_INHIBITION", "raw_endpoint_name": "CYP2D6 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP2D6", "raw_value": 100.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 100000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Review Theo-Dur", "assay_type": "Recombinant human CYP2D6 Dextromethorphan assay", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 40000.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -1.40, "normalized_unit": "log10(mol/L)", "reference_text": "DrugCentral / Theo-Dur NDA 085078", "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True}
        ]
    },
    # 49: Tolbutamide (LOCKED_FINAL_TEST_COHORT_3)
    {
        "name": "Tolbutamide", "cas_number": "64-77-7", "drugbank_id": "DB01124", "pubchem_cid": "5505", "chembl_id": "CHEMBL453", "unii": "983W2N6I28",
        "smiles": "Cc1ccc(S(=O)(=O)NC(=O)NCCCC)cc1", "indication": "Type 2 Diabetes Mellitus", "target": "Sulfonylurea Receptor 1 (SUR1) / FDA Standard CYP2C9 Probe", "scaffold_family": "Sulfonylurea",
        "model_role": "LOCKED_FINAL_TEST_COHORT_3", "cohort": "LOCKED_FINAL_TEST_COHORT_3",
        "upstream_overlap": {"SOLUBILITY_GENERIC": "NOVEL_IN_DOMAIN", "SOLUBILITY_THERMODYNAMIC": "NOVEL_IN_DOMAIN", "HUMAN_PPB": "NOVEL_IN_DOMAIN", "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT", "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT", "HERG_LIABILITY": "VALIDATION_HOLDOUT"},
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 96.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 96.0, "normalized_unit": "% bound", "reference_text": "FDA Orinase NDA 010708 Clinical Pharmacology Review", "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp", "raw_value": 150.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 150000.0, "normalized_unit": "nM", "reference_text": "Redfern et al. / ChEMBL453", "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP3A4", "raw_value": 50.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 50000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Guidance Review", "assay_type": "Recombinant human CYP3A4 Midazolam 1'-hydroxylation assay", "training_eligible": True},
            {"canonical_endpoint_id": "CYP2D6_INHIBITION", "raw_endpoint_name": "CYP2D6 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP2D6", "raw_value": 80.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 80000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Review Orinase", "assay_type": "Recombinant human CYP2D6 Dextromethorphan assay", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 400.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -3.40, "normalized_unit": "log10(mol/L)", "reference_text": "DrugCentral / Orinase NDA 010708", "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True}
        ]
    },
    # 50: Trazodone (LOCKED_FINAL_TEST_COHORT_3)
    {
        "name": "Trazodone", "cas_number": "19794-93-5", "drugbank_id": "DB00656", "pubchem_cid": "5533", "chembl_id": "CHEMBL643", "unii": "YBK48BXK30",
        "smiles": "O=c1n(CCCN2CCN(c3cccc(Cl)c3)CC2)nc2ccccn12", "indication": "Major Depressive Disorder", "target": "5-HT2A Receptor / SARI / CYP3A4 Substrate", "scaffold_family": "Triazolopyridine-piperazine",
        "model_role": "LOCKED_FINAL_TEST_COHORT_3", "cohort": "LOCKED_FINAL_TEST_COHORT_3",
        "upstream_overlap": {"SOLUBILITY_GENERIC": "NOVEL_IN_DOMAIN", "SOLUBILITY_THERMODYNAMIC": "NOVEL_IN_DOMAIN", "HUMAN_PPB": "NOVEL_IN_DOMAIN", "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT", "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT", "HERG_LIABILITY": "VALIDATION_HOLDOUT"},
        "observations": [
            {"canonical_endpoint_id": "HUMAN_PPB", "raw_endpoint_name": "Human Plasma Protein Binding", "section": "ADMET", "species": "Homo sapiens", "matrix": "Human Plasma", "raw_value": 90.0, "raw_unit": "%", "raw_relation": "=", "normalized_value": 90.0, "normalized_unit": "% bound", "reference_text": "FDA Desyrel NDA 018207 Clinical Pharmacology Review", "assay_type": "Equilibrium Dialysis (Human Plasma)", "training_eligible": True},
            {"canonical_endpoint_id": "HERG_LIABILITY", "raw_endpoint_name": "hERG Potassium Channel Inhibition", "section": "SAFETY", "species": "Homo sapiens", "matrix": "HEK293 Whole-Cell Patch-Clamp", "raw_value": 4.5, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 4500.0, "normalized_unit": "nM", "reference_text": "Redfern et al. / ChEMBL643", "assay_type": "Whole-cell Voltage Patch-Clamp Electrophysiology", "training_eligible": True},
            {"canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_endpoint_name": "CYP3A4 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP3A4", "raw_value": 15.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 15000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Guidance Review", "assay_type": "Recombinant human CYP3A4 Midazolam 1'-hydroxylation assay", "training_eligible": True},
            {"canonical_endpoint_id": "CYP2D6_INHIBITION", "raw_endpoint_name": "CYP2D6 Direct Reversible Inhibition", "section": "METABOLISM", "species": "Homo sapiens", "matrix": "rhCYP2D6", "raw_value": 35.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": 35000.0, "normalized_unit": "nM", "reference_text": "FDA In Vitro Metabolism Review Desyrel", "assay_type": "Recombinant human CYP2D6 Dextromethorphan assay", "training_eligible": True},
            {"canonical_endpoint_id": "SOLUBILITY_THERMODYNAMIC", "raw_endpoint_name": "Thermodynamic Aqueous Solubility", "section": "PHYSICOCHEMICAL", "species": "None", "matrix": "Phosphate Buffer pH 7.0", "raw_value": 45.0, "raw_unit": "µM", "raw_relation": "=", "normalized_value": -4.35, "normalized_unit": "log10(mol/L)", "reference_text": "DrugCentral / Desyrel NDA 018207", "assay_type": "Shake-Flask Thermodynamic Solubility (HPLC/UV)", "training_eligible": True}
        ]
    }
]

ALL_50_DRUGS = ORIGINAL_40 + NEW_10_DRUGS
print(f"Total Combined Drugs: {len(ALL_50_DRUGS)}")

with open("/home/xavier/chem/drug-opt-platform/backend/reference_drugs_50.json", "w") as f:
    json.dump(ALL_50_DRUGS, f, indent=2)
print("Saved reference_drugs_50.json successfully")
