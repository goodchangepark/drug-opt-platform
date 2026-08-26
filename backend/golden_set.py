"""Golden Reference Structure Set & Verification Gate (Stage 4C-1).

Includes 52 challenging reference structures covering neutrals, acids, bases,
zwitterions, salts, hydrates, tautomers, stereoisomers, isotopes, heterocycles,
quaternary ammoniums, sulfonamides, phosphates, boron, organometallics,
multicomponent mixtures, covalent warheads, peptides, and steroids.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from backend.standardizer import RDKIT_VERSION, STANDARDIZER_NAME, STANDARDIZER_VERSION, standardize_molecule

GOLDEN_FILE_PATH = Path(__file__).parent / "golden_set.json"

RAW_GOLDEN_INPUTS = [
    # 1-4 Neutral, Acid, Base, Zwitterion
    {"id": "G01", "name": "Aspirin (Neutral)", "smiles": "CC(=O)Oc1ccccc1C(=O)O", "category": "Neutral"},
    {"id": "G02", "name": "Ibuprofen (Simple Acid)", "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "category": "Acid"},
    {"id": "G03", "name": "Amphetamine (Simple Base)", "smiles": "CC(Cc1ccccc1)N", "category": "Base"},
    {"id": "G04", "name": "Glycine (Zwitterion)", "smiles": "[NH3+]CC(=O)[O-]", "category": "Zwitterion"},

    # 5-8 Salts & Hydrates
    {"id": "G05", "name": "Metformin Hydrochloride (Salt)", "smiles": "CC(C)N=C(N)N.Cl", "category": "Salt"},
    {"id": "G06", "name": "Naproxen Sodium (Sodium Salt)", "smiles": "CC(c1ccc2cc(oc2c1)OC)C(=O)[O-].[Na+]", "category": "Salt"},
    {"id": "G07", "name": "Zolpidem Tartrate (Organic Salt)", "smiles": "Cc1ccc2n1c(cc2C(=O)N(C)C)c3ccc(C)cc3.OC(C(=O)O)C(O)C(=O)O", "category": "Salt"},
    {"id": "G08", "name": "Amoxicillin Trihydrate (Hydrate)", "smiles": "CC1(C(N2C(S1)C(C2=O)NC(=O)C(c3ccc(cc3)O)N)C(=O)O)C.O.O.O", "category": "Hydrate"},

    # 9-10 Tautomer Pairs
    {"id": "G09", "name": "2-Pyridone (Tautomer A)", "smiles": "O=c1cccc[nH]1", "category": "Tautomer"},
    {"id": "G10", "name": "2-Hydroxypyridine (Tautomer B)", "smiles": "Oc1ccccn1", "category": "Tautomer"},

    # 11-14 Stereoisomers & Alkene E/Z
    {"id": "G11", "name": "(R)-Thalidomide", "smiles": "O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1", "category": "Stereoisomer"},
    {"id": "G12", "name": "(S)-Thalidomide", "smiles": "O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1", "category": "Stereoisomer"},
    {"id": "G13", "name": "E-Stilbene", "smiles": "C(=C/c1ccccc1)\\c2ccccc2", "category": "Alkene E/Z"},
    {"id": "G14", "name": "Z-Stilbene", "smiles": "C(=C\\c1ccccc1)/c2ccccc2", "category": "Alkene E/Z"},

    # 15-16 Isotope Labeled
    {"id": "G15", "name": "Deuterated Water (Isotope)", "smiles": "[2H]O[2H]", "category": "Isotope"},
    {"id": "G16", "name": "Carbon-13 Glucose (Isotope)", "smiles": "C([13C]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O", "category": "Isotope"},

    # 17-20 Heterocycles & Charged Species
    {"id": "G17", "name": "Imidazole", "smiles": "c1c[nH]cn1", "category": "Heterocycle"},
    {"id": "G18", "name": "Sildenafil", "smiles": "CCCC1=NN(C)C2=C1N=C(C)NC2=O", "category": "Heterocycle"},
    {"id": "G19", "name": "Acetylcholine (Quaternary Ammonium)", "smiles": "CC(=O)OCC[N+](C)(C)C", "category": "Quaternary Ammonium"},
    {"id": "G20", "name": "Paraquat (Pyridinium Salt)", "smiles": "C[n+]1ccc(cc1)c2cc[n+](C)cc2.[Cl-].[Cl-]", "category": "Charged Salt"},

    # 21-24 Sulfonamides & Phosphates
    {"id": "G21", "name": "Sulfamethoxazole (Primary Sulfonamide)", "smiles": "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1", "category": "Sulfonamide"},
    {"id": "G22", "name": "Glibenclamide (Secondary Sulfonamide)", "smiles": "COc1ccc(CC(=O)NCCc2ccc(S(=O)(=O)NC(=O)NC3CCCCC3)cc2)cc1Cl", "category": "Sulfonamide"},
    {"id": "G23", "name": "Fostamatinib (Phosphate)", "smiles": "COc1cc2c(cc1OP(=O)(O)O)c(n2)Nc3ccc(F)c(Cl)c3", "category": "Phosphate"},
    {"id": "G24", "name": "Tenofovir (Phosphonate)", "smiles": "CC(CN1C=NC2=C(N)N=CN=C21)OCP(=O)(O)O", "category": "Phosphonate"},

    # 25-27 Boron, Metal & Multi-Component Mixture
    {"id": "G25", "name": "Bortezomib (Boronate)", "smiles": "B(C(CC(C)C)NC(=O)C(CC1=CC=CC=C1)NC(=O)C2=NC=CN=C2)(O)O", "category": "Boron"},
    {"id": "G26", "name": "Cisplatin (Metal Complex)", "smiles": "Cl[Pt](Cl)(N)N", "category": "Metal"},
    {"id": "G27", "name": "Aspirin + Paracetamol (50:50 Mixture)", "smiles": "CC(=O)Oc1ccccc1C(=O)O.CC(=O)Nc1ccc(O)cc1", "category": "Multicomponent"},

    # 28-31 Warheads, Peptides, Macrocycles
    {"id": "G28", "name": "Ibrutinib (Acrylamide Warhead)", "smiles": "C=CC(=O)N1CCC(CC1)N2C3=NC=NC(=C3C(=N2)C4=CC=C(C=C4)OC5=CC=CC=C5)N", "category": "Warhead"},
    {"id": "G29", "name": "Carfilzomib (Epoxide Warhead)", "smiles": "CC(C)CC(C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(C)C)C(=O)C2(CO2)C)NC(=O)C(Cc3ccccc3)NC(=O)CN4CCOCC4", "category": "Warhead"},
    {"id": "G30", "name": "Leu-Enkephalin (Peptide)", "smiles": "CC(C)CC(C(=O)O)NC(=O)C(CC1=CC=CC=C1)NC(=O)CNC(=O)CNC(=O)C(CC2=CC=C(C=C2)O)N", "category": "Peptide"},
    {"id": "G31", "name": "Cyclosporine A (Macrocycle)", "smiles": "CCC1C(=O)N(CC(=O)N(C(C(=O)NC(C(=O)N(C(C(=O)NC(C(=O)NC(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N1)C(C(C)CC=CC)O)C)C(C)C)C)CC(C)C)C)C(C)C)C)CC(C)C)C)C)C)C)C", "category": "Macrocycle"},

    # 32-36 Fluorinated, Halogenated, Sulfoxides
    {"id": "G32", "name": "Atorvastatin (Fluorinated Aromatic)", "smiles": "CC(C)c1c(c(c(n1CCC(CC(CC(=O)O)O)O)c2ccc(cc2)F)c3ccccc3)C(=O)Nc4ccccc4", "category": "Fluorinated"},
    {"id": "G33", "name": "Fluoxetine (Trifluoromethyl)", "smiles": "CNCCC(c1ccccc1)Oc2ccc(cc2)C(F)(F)F", "category": "Fluorinated"},
    {"id": "G34", "name": "Amiodarone (Poly-Iodinated)", "smiles": "CCCCc1oc2ccccc2c1C(=O)c3cc(I)c(OCCN(CC)CC)c(I)c3", "category": "Halogenated"},
    {"id": "G35", "name": "Omeprazole (Sulfoxide)", "smiles": "COc1ccc2[nH]c(nc2c1)S(=O)Cc3ncc(C)c(OC)c3C", "category": "Sulfoxide"},
    {"id": "G36", "name": "Dapsone (Sulfone)", "smiles": "Nc1ccc(cc1)S(=O)(=O)c2ccc(N)cc2", "category": "Sulfone"},

    # 37-41 Thiones, Oximes, N-Oxides, Amidines
    {"id": "G37", "name": "Propylthiouracil (Thione)", "smiles": "CCCC1=CC(=O)NC(=S)N1", "category": "Thione"},
    {"id": "G38", "name": "Pralidoxime (Oxime)", "smiles": "C[n+]1ccccc1C=NO.[Cl-]", "category": "Oxime"},
    {"id": "G39", "name": "Chlordiazepoxide N-Oxide", "smiles": "CNC1=NC2=C(C=C(C=C2)Cl)C(=[N+](C1)O)C3=CC=CC=C3", "category": "N-Oxide"},
    {"id": "G40", "name": "Guanfacine (Guanidine)", "smiles": "NC(=N)NC(=O)Cc1c(Cl)cccc1Cl", "category": "Guanidine"},
    {"id": "G41", "name": "Pentamidine (Amidine)", "smiles": "NC(=N)c1ccc(OCCCCCOc2ccc(cc2)C(=N)N)cc1", "category": "Amidine"},

    # 42-45 Heterocyclic Scaffolds & Tetrazoles
    {"id": "G42", "name": "Valsartan Tetrazole", "smiles": "CCCCC(=O)N(Cc1ccc(cc1)c2ccccc2c3nnn[nh]3)C(C(C)C)C(=O)O", "category": "Tetrazole"},
    {"id": "G43", "name": "Ziprasidone (Isothiazole)", "smiles": "c1ccc2c(c1)c(no2)N3CCN(CC3)CCc4ccc5c(c4)NC(=O)C5", "category": "Isothiazole"},
    {"id": "G44", "name": "Risperidone (Benzisoxazole)", "smiles": "CC1=C(C(=O)N2CCC(CC2)C3=NOC4=C3C=CC(=C4)F)AT1=O", "category": "Benzisoxazole"},
    {"id": "G45", "name": "Caffeine (Purine)", "smiles": "Cn1cnc2n(C)c(=O)n(C)c(=O)c12", "category": "Purine"},

    # 46-49 Steroids, Sugars, Prodrugs
    {"id": "G46", "name": "Fluorouracil (Pyrimidine)", "smiles": "O=c1[nH]cc(F)c(=O)[nH]1", "category": "Pyrimidine"},
    {"id": "G47", "name": "Dexamethasone (Steroid)", "smiles": "CC1CC2C3CCC4CC(=O)C=CC4(C3(C(CC2(C1(C(=O)CO)O)C)O)F)C", "category": "Steroid"},
    {"id": "G48", "name": "Empagliflozin (Sugar Derivative)", "smiles": "Cc1ccc(cc1C(=O)O)C2C(C(C(C(O2)CO)O)O)O", "category": "Sugar"},
    {"id": "G49", "name": "Enalapril (Ester Prodrug)", "smiles": "CCOC(=O)C(CCc1ccccc1)NC(C)C(=O)N2CCCC2C(=O)O", "category": "Prodrug"},

    # 50-52 Edge Cases / Invalid Input
    {"id": "G50", "name": "Octreotide (Unnatural Cyclic Peptide)", "smiles": "CC(C(C(=O)NC(CO)C(=O)O)NC(=O)C(CC1=c2ccccc2[nH]c1)NC(=O)C(Cc3ccccc3)NC(=O)C(CSSCC4NC(=O)C(Cc5ccccc5)NC(=O)C(NC(=O)C4N)C(C)O)O)O", "category": "Peptide"},
    {"id": "G51", "name": "Acetaminophen (Paracetamol)", "smiles": "CC(=O)Nc1ccc(O)cc1", "category": "Neutral"},
    {"id": "G52", "name": "Invalid SMILES String", "smiles": "INVALID_SMILES_STRING_123", "category": "Invalid"},
]


def generate_golden_reference_dataset() -> dict[str, Any]:
    """Generates the canonical golden reference dictionary by standardizing all 52 structures."""
    items = []
    for raw in RAW_GOLDEN_INPUTS:
        res = standardize_molecule(raw["smiles"])
        items.append({
            "id": raw["id"],
            "name": raw["name"],
            "category": raw["category"],
            "input_smiles": raw["smiles"],
            "status": res["status"],
            "canonical_smiles": res["canonical_smiles"],
            "isomeric_smiles": res["isomeric_smiles"],
            "inchikey": res["inchikey"],
            "num_heavy_atoms": res.get("num_heavy_atoms"),
            "salt_extracted": res.get("salt_extracted", False),
            "warnings": res.get("warnings", []),
        })

    return {
        "dataset_name": "CHEM_GOLDEN_REFERENCE_SET_V1",
        "total_items": len(items),
        "standardizer_name": STANDARDIZER_NAME,
        "standardizer_version": STANDARDIZER_VERSION,
        "rdkit_version": RDKIT_VERSION,
        "items": items,
    }


def save_golden_set():
    data = generate_golden_reference_dataset()
    GOLDEN_FILE_PATH.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return GOLDEN_FILE_PATH


def load_golden_set() -> dict[str, Any]:
    if not GOLDEN_FILE_PATH.exists():
        save_golden_set()
    return json.loads(GOLDEN_FILE_PATH.read_text(encoding="utf-8"))


def run_golden_gate_test() -> dict[str, Any]:
    """Runs standardizer against stored golden reference set and produces structured diff report."""
    golden_ref = load_golden_set()
    current_rdkit = RDKIT_VERSION

    passed_count = 0
    failed_count = 0
    diffs = []

    for item in golden_ref["items"]:
        res = standardize_molecule(item["input_smiles"])
        match_status = (res["status"] == item["status"])
        match_cano = (res["canonical_smiles"] == item["canonical_smiles"])
        match_iso = (res["isomeric_smiles"] == item["isomeric_smiles"])
        match_key = (res["inchikey"] == item["inchikey"])

        if match_status and match_cano and match_iso and match_key:
            passed_count += 1
        else:
            failed_count += 1
            diffs.append({
                "id": item["id"],
                "name": item["name"],
                "expected": {
                    "status": item["status"],
                    "canonical_smiles": item["canonical_smiles"],
                    "isomeric_smiles": item["isomeric_smiles"],
                    "inchikey": item["inchikey"],
                },
                "actual": {
                    "status": res["status"],
                    "canonical_smiles": res["canonical_smiles"],
                    "isomeric_smiles": res["isomeric_smiles"],
                    "inchikey": res["inchikey"],
                },
            })

    total = len(golden_ref["items"])
    gate_passed = (failed_count == 0)

    return {
        "gate_passed": gate_passed,
        "dataset_name": golden_ref["dataset_name"],
        "total_items": total,
        "passed_count": passed_count,
        "failed_count": failed_count,
        "reference_rdkit_version": golden_ref["rdkit_version"],
        "current_rdkit_version": current_rdkit,
        "diffs": diffs,
    }


if __name__ == "__main__":
    save_golden_set()
    print(f"Saved golden reference set ({len(RAW_GOLDEN_INPUTS)} items) to {GOLDEN_FILE_PATH}")
