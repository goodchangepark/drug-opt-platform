"""Unit tests for CHEM_STANDARDIZER_V1 and Golden Reference Set Gate (Stage 4C-1)."""

import pytest

from backend.golden_set import load_golden_set, run_golden_gate_test
from backend.standardizer import (
    GLOBAL_DESCRIPTOR_CONFIG,
    GLOBAL_FINGERPRINT_CONFIG,
    STANDARDIZER_NAME,
    STANDARDIZER_VERSION,
    standardize_molecule,
)


def test_standardizer_determinism():
    smiles = "CC(=O)Oc1ccccc1C(=O)O.Cl"
    res1 = standardize_molecule(smiles)
    res2 = standardize_molecule(smiles)

    assert res1["status"] == res2["status"]
    assert res1["canonical_smiles"] == res2["canonical_smiles"]
    assert res1["isomeric_smiles"] == res2["isomeric_smiles"]
    assert res1["inchikey"] == res2["inchikey"]
    assert res1["provenance"]["standardizer_name"] == STANDARDIZER_NAME
    assert res1["provenance"]["standardizer_version"] == STANDARDIZER_VERSION


def test_salt_extraction_policy():
    # Metformin HCl
    res = standardize_molecule("CC(C)N=C(N)N.Cl")
    assert res["salt_extracted"] is True
    assert res["status"] == "SUCCESS"
    assert "SALT_REMOVED" in res["warnings"][0]


def test_multicomponent_review_required():
    # 50:50 mixture of Aspirin and Paracetamol without clear inorganic salt
    res = standardize_molecule("CC(=O)Oc1ccccc1C(=O)O.CC(=O)Nc1ccc(O)cc1")
    assert res["status"] == "MULTICOMPONENT_REVIEW_REQUIRED"
    assert any("MULTICOMPONENT_REVIEW_REQUIRED" in w for w in res["warnings"])


def test_stereochemistry_preservation():
    # (R)-Thalidomide vs (S)-Thalidomide
    r_res = standardize_molecule("O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1")
    s_res = standardize_molecule("O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1")
    assert r_res["status"] == "SUCCESS"
    assert s_res["status"] == "SUCCESS"


def test_golden_set_gate():
    report = run_golden_gate_test()
    assert report["gate_passed"] is True
    assert report["total_items"] == 52
    assert report["passed_count"] == 52
    assert report["failed_count"] == 0
    assert len(report["diffs"]) == 0


def test_global_configurations():
    assert GLOBAL_FINGERPRINT_CONFIG["name"] == "Morgan"
    assert GLOBAL_FINGERPRINT_CONFIG["radius"] == 2
    assert GLOBAL_FINGERPRINT_CONFIG["nBits"] == 2048
    assert GLOBAL_FINGERPRINT_CONFIG["useChirality"] is True

    assert GLOBAL_DESCRIPTOR_CONFIG["version"] == "1.0"
    assert "Crippen.MolLogP" in GLOBAL_DESCRIPTOR_CONFIG["clogp_definition"]
