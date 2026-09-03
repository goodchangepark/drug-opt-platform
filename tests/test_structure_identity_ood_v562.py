"""Unit and integration tests for Structure Identity & OOD Correction v5.6.2."""
import pytest
from backend.database import SessionLocal
from backend.models import Compound, CompoundVersion
from backend.openadmet_cyp import (
    evaluate_cyp_applicability_domain,
    predict_chemeleon_cyp_pic50,
    compute_fold_error,
)
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, inchi


def test_hard_acceptance_compound_identities_and_mw():
    """Verify formula, MW, and free-base structure identity for Mobocertinib, Orforglipron, Poziotinib, Sunvozertinib."""
    db = SessionLocal()
    try:
        # 1. Mobocertinib: C32H39N7O4, MW ~ 585.71
        mobo = db.query(Compound).filter(Compound.name.ilike("%Mobocertinib%")).first()
        assert mobo is not None
        cv_mobo = mobo.versions[-1]
        mol_mobo = Chem.MolFromSmiles(cv_mobo.canonical_smiles)
        assert rdMolDescriptors.CalcMolFormula(mol_mobo) == "C32H39N7O4"
        assert Descriptors.MolWt(mol_mobo) == pytest.approx(585.71, abs=0.1)

        # 2. Orforglipron: C48H48F2N10O5, MW ~ 882.97
        orfor = db.query(Compound).filter(Compound.name.ilike("%Orforglipron%")).first()
        assert orfor is not None
        cv_orfor = orfor.versions[-1]
        mol_orfor = Chem.MolFromSmiles(cv_orfor.canonical_smiles)
        assert rdMolDescriptors.CalcMolFormula(mol_orfor) == "C48H48F2N10O5"
        assert Descriptors.MolWt(mol_orfor) == pytest.approx(882.97, abs=0.1)

        # 3. Poziotinib: C23H21Cl2FN4O3, MW ~ 491.35
        pozio = db.query(Compound).filter(Compound.name.ilike("%Poziotinib%")).first()
        assert pozio is not None
        cv_pozio = pozio.versions[-1]
        mol_pozio = Chem.MolFromSmiles(cv_pozio.canonical_smiles)
        assert rdMolDescriptors.CalcMolFormula(mol_pozio) == "C23H21Cl2FN4O3"
        assert Descriptors.MolWt(mol_pozio) == pytest.approx(491.35, abs=0.1)

        # 4. Sunvozertinib: C29H35ClFN7O3, MW ~ 584.10
        sunvo = db.query(Compound).filter(Compound.name.ilike("%Sunvozertinib%")).first()
        assert sunvo is not None
        cv_sunvo = sunvo.versions[-1]
        mol_sunvo = Chem.MolFromSmiles(cv_sunvo.canonical_smiles)
        assert rdMolDescriptors.CalcMolFormula(mol_sunvo) == "C29H35ClFN7O3"
        assert Descriptors.MolWt(mol_sunvo) == pytest.approx(584.10, abs=0.1)
    finally:
        db.close()


def test_applicability_domain_granular_attribution():
    """Verify AD evaluates MW/cLogP envelope + Morgan nearest neighbor similarity with explicit reasons."""
    # Poziotinib should be IN_DOMAIN with high similarity to kinase/CYP scaffolds
    pozio_mol = Chem.MolFromSmiles("C=CC(=O)N1CCC(Oc2cc3c(Nc4ccc(Cl)c(Cl)c4F)ncnc3cc2OC)CC1")
    status_pozio, sim_pozio, viol_pozio, metrics_pozio, reason_pozio = evaluate_cyp_applicability_domain(pozio_mol)
    assert status_pozio == "IN_DOMAIN"
    assert sim_pozio >= 0.30
    assert len(viol_pozio) == 0
    assert metrics_pozio["MW"] == pytest.approx(491.35, abs=0.1)

    # Mobocertinib should be BORDERLINE with 0 violations and MW ~ 585.7
    mobo_mol = Chem.MolFromSmiles("C=CC(=O)Nc1cc(Nc2ncc(C(=O)OC(C)C)c(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C")
    status_mobo, sim_mobo, viol_mobo, metrics_mobo, reason_mobo = evaluate_cyp_applicability_domain(mobo_mol)
    assert status_mobo == "BORDERLINE"
    assert len(viol_mobo) == 0
    assert metrics_mobo["MW"] == pytest.approx(585.71, abs=0.1)

    # Orforglipron should be OUT_OF_DOMAIN due to MW > 800 and cLogP > 6.5
    orfor_mol = Chem.MolFromSmiles("Cc1cc(-n2nc3c(c2-n2ccn(-c4ccc5c(cnn5C)c4F)c2=O)[C@H](C)N(C(=O)c2cc4cc([C@H]5CCOC(C)(C)C5)ccc4n2[C@@]2(c4noc(=O)[nH]4)C[C@@H]2C)CC3)cc(C)c1F")
    status_orfor, sim_orfor, viol_orfor, metrics_orfor, reason_orfor = evaluate_cyp_applicability_domain(orfor_mol)
    assert status_orfor == "OUT_OF_DOMAIN"
    assert any("MW" in v for v in viol_orfor)
    assert metrics_orfor["MW"] == pytest.approx(882.97, abs=0.1)


def test_cyp_quantitative_prediction_with_recomputed_ad():
    """Verify quantitative CYP predictions output correct AD status and numeric values."""
    pred_mobo = predict_chemeleon_cyp_pic50("C=CC(=O)Nc1cc(Nc2ncc(C(=O)OC(C)C)c(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C", "CYP1A2")
    assert pred_mobo.applicability_domain == "BORDERLINE"
    assert pred_mobo.pic50 == pytest.approx(5.96, abs=0.1)
    assert pred_mobo.ic50_um == pytest.approx(1.10, abs=0.2)

    pred_pozio = predict_chemeleon_cyp_pic50("C=CC(=O)N1CCC(Oc2cc3c(Nc4ccc(Cl)c(Cl)c4F)ncnc3cc2OC)CC1", "CYP3A4")
    assert pred_pozio.applicability_domain == "IN_DOMAIN"
    assert pred_pozio.pic50 == pytest.approx(7.12, abs=0.1)
