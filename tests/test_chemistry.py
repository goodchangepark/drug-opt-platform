import pytest
from rdkit import Chem

from backend.chemistry import ChemistryError, analyze_smiles

REFERENCE_COMPOUNDS = [
    ("caffeine", "Cn1cnc2n(C)c(=O)n(C)c(=O)c12", 194.19, -0.07, 58.44, 0, 6),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O", 180.16, 1.31, 63.60, 1, 3),
    ("acetaminophen", "CC(=O)Nc1ccc(O)cc1", 151.16, 0.46, 49.33, 2, 2),
    ("ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1", 206.28, 3.97, 37.30, 1, 1),
    ("gefitinib", "COCCCNc1ncnc2cc(Cl)c(Nc3ccc(C#C)c(OC)c3)nc12", 446.90, 3.36, 68.43, 2, 7),
    ("osimertinib", "CS(=O)(=O)CCNc1ncnc2cc(Cl)c(Nc3ccc(C#C)c(OC)c3)nc12", 499.61, 2.77, 78.27, 2, 8),
]


@pytest.mark.parametrize("name,smiles,mw,clogp,tpsa,hbd,hba", REFERENCE_COMPOUNDS)
def test_reference_compound_ranges(name, smiles, mw, clogp, tpsa, hbd, hba):
    result = analyze_smiles(smiles)
    props = result["properties"]
    assert result["identity"]["inchikey"]
    assert abs(props["molecular_weight"] - mw) <= max(1.5, mw * 0.15)
    assert abs(props["clogp"] - clogp) <= 1.2
    assert abs(props["tpsa"] - tpsa) <= max(12, tpsa * 0.45)
    assert props["hbd"] == hbd or abs(props["hbd"] - hbd) <= 2
    assert props["hba"] >= hba or abs(props["hba"] - hba) <= 2
    assert 0 <= props["qed"] <= 1


def test_identity_stereochemistry_and_charge():
    chiral = analyze_smiles("C[C@H](N)C(=O)O")
    assert "@" in chiral["identity"]["isomeric_smiles"]
    charged = analyze_smiles("CC(=O)[O-]")
    assert charged["properties"]["formal_charge"] == -1
    salt = analyze_smiles("[Na+]CC(=O)[O-]")
    assert salt["identity"]["inchikey"] != analyze_smiles("CCC(=O)[O-]")["identity"]["inchikey"]


def test_invalid_and_empty_smiles():
    with pytest.raises(ChemistryError):
        analyze_smiles("not_a_molecule")
    with pytest.raises(ChemistryError):
        analyze_smiles("")


def test_duplicate_detection_and_large_molecule():
    aspirin = analyze_smiles("CC(=O)Oc1ccccc1C(=O)O")
    duplicate = analyze_smiles("OC(=O)c1ccccc1OC(C)=O")
    assert aspirin["identity"]["canonical_smiles"] == duplicate["identity"]["canonical_smiles"]
    large = analyze_smiles("C" * 60 + "O")
    assert large["properties"]["heavy_atom_count"] > 50


def test_drug_likeness_failure_reasons():
    result = analyze_smiles("CCCCCCCCCCCCCCCCCC(=O)O")
    lipinski = result["rules"]["Lipinski Rule of Five"]
    assert lipinski["result"] == "FAIL"
    assert any("cLogP" in reason for reason in lipinski["reasons"])
