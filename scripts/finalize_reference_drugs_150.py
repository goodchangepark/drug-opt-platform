"""
Merges existing 100 drugs and newly curated 50 drugs into backend/reference_drugs_150.json.
Validates all structures, calculates physicochemical properties, partitions cohorts,
and verifies schema consistency.
"""
from __future__ import annotations

import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski

from scripts.build_complete_150_dataset import NEW_DRUGS_50 as BLOCK_1
from scripts.append_remaining_drugs import BLOCK_2
from scripts.curate_blocks_3_and_4 import BLOCK_3, BLOCK_4

NEW_50 = BLOCK_1 + BLOCK_2 + BLOCK_3 + BLOCK_4
print(f"Total new drugs collected: {len(NEW_50)}")
assert len(NEW_50) == 50, f"Expected 50 new drugs, got {len(NEW_50)}"

# Load existing 100 drugs
with open("backend/reference_drugs_100.json", "r") as f:
    drugs_100 = json.load(f)

# Update roles for drugs 81-100 in existing_100
for i, d in enumerate(drugs_100[80:100], 81):
    if i <= 90:
        d["model_role"] = "DEVELOPMENT_TRAINING"
        d["cohort"] = "DEV_TRAINING"
    else:
        d["model_role"] = "MODEL_SELECTION_VALIDATION"
        d["cohort"] = "MODEL_SELECTION_VALIDATION"
    for obs in d.get("observations", []):
        obs["endpoint_role"] = d["model_role"]

# Validate and format the 50 new drugs
for i, d in enumerate(NEW_50, 101):
    mol = Chem.MolFromSmiles(d["smiles"])
    assert mol is not None, f"Invalid SMILES for {d['name']}: {d['smiles']}"
    canon_smi = Chem.MolToSmiles(mol, canonical=True)
    d["smiles"] = canon_smi
    d["molecular_weight"] = round(float(Descriptors.MolWt(mol)), 2)
    d["clogp"] = round(float(Crippen.MolLogP(mol)), 2)
    d["inchikey"] = Chem.MolToInchiKey(mol)

    # Cohort Partitioning:
    # 101-125: DEVELOPMENT_TRAINING (N=25)
    # 126-137: MODEL_SELECTION_VALIDATION (N=12)
    # 138-150: LOCKED_FINAL_TEST_COHORT_6 (N=13)
    if i <= 125:
        role = "DEVELOPMENT_TRAINING"
        cohort = "DEV_TRAINING"
    elif i <= 137:
        role = "MODEL_SELECTION_VALIDATION"
        cohort = "MODEL_SELECTION_VALIDATION"
    else:
        role = "LOCKED_FINAL_TEST_COHORT_6"
        cohort = "LOCKED_FINAL_TEST_COHORT_6"

    d["model_role"] = role
    d["cohort"] = cohort
    d["upstream_overlap"] = {
        "SOLUBILITY_GENERIC": "VALIDATION_HOLDOUT",
        "SOLUBILITY_THERMODYNAMIC": "VALIDATION_HOLDOUT",
        "HUMAN_PPB": "VALIDATION_HOLDOUT",
        "CYP3A4_INHIBITION": "VALIDATION_HOLDOUT",
        "CYP2D6_INHIBITION": "VALIDATION_HOLDOUT",
        "HERG_LIABILITY": "VALIDATION_HOLDOUT",
        "HLM_CLINT": "VALIDATION_HOLDOUT",
        "CACO2_PAPP_AB": "VALIDATION_HOLDOUT",
        "CYP1A2_INHIBITION": "VALIDATION_HOLDOUT",
        "CYP2C9_INHIBITION": "VALIDATION_HOLDOUT"
    }
    for obs in d.get("observations", []):
        obs["endpoint_role"] = role

# Merge 100 + 50
FULL_150 = drugs_100 + NEW_50
print(f"Total merged reference drugs: {len(FULL_150)}")
assert len(FULL_150) == 150

# Check uniqueness
names = [d["name"] for d in FULL_150]
dbids = [d["drugbank_id"] for d in FULL_150]
inchis = [d.get("inchikey") or Chem.MolToInchiKey(Chem.MolFromSmiles(d["smiles"])) for d in FULL_150]
assert len(set(names)) == 150, f"Duplicate names: {len(names)} vs {len(set(names))}"
assert len(set(dbids)) == 150, f"Duplicate DBIDs: {len(dbids)} vs {len(set(dbids))}"
assert len(set(inchis)) == 150, f"Duplicate InChIKeys: {len(inchis)} vs {len(set(inchis))}"

# Save to backend/reference_drugs_150.json
out_path = Path("backend/reference_drugs_150.json")
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(FULL_150, f, indent=2)

print(f"Successfully generated {out_path} ({out_path.stat().st_size} bytes)")

cohorts = {}
for d in FULL_150:
    c = d.get("cohort", "UNKNOWN")
    cohorts[c] = cohorts.get(c, 0) + 1

print("\nFinal DrugBank 150 Cohort Distribution:")
for k, v in sorted(cohorts.items()):
    print(f"  {k}: {v}")
