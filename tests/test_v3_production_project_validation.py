"""
End-to-End Validation for Global Prediction Engine v3.0 Production Release on a New Test Project.
"""
import pytest
from rdkit import Chem
from sqlalchemy import select

from backend.database import SessionLocal
from backend.models import Project, Compound, CompoundVersion, ExternalExperimentalEvidence
from backend.admet import ADMETPredictionRun
from backend.engine_v3_learning import evaluate_global_engine_v3_readiness, predict_global_v3_endpoint


def test_v3_production_release_governance_and_promotion_statuses():
    """Verify that CYP3A4, CYP2D6, Solubility are GLOBAL_V3_PRIMARY; PPB is V3_CANDIDATE; hERG is RETAIN_BASE."""
    db = SessionLocal()
    try:
        readiness = evaluate_global_engine_v3_readiness(db)
        assert readiness["release_status"] == "GLOBAL_ENGINE_V3_PRODUCTION_RELEASE"

        statuses = {ep["endpoint_id"]: ep["promotion_status"] for ep in readiness["endpoints_evaluated"]}
        assert statuses["CYP3A4_INHIBITION"] == "GLOBAL_V3_PRIMARY"
        assert statuses["CYP2D6_INHIBITION"] == "GLOBAL_V3_PRIMARY"
        assert statuses["SOLUBILITY_GENERIC"] == "GLOBAL_V3_PRIMARY"
        assert statuses["HUMAN_PPB"] == "V3_CANDIDATE"
        assert statuses["HERG_LIABILITY"] == "RETAIN_BASE"

        # Verify separated Validation vs Final-Test performance metrics
        for ep in readiness["endpoints_evaluated"]:
            assert "validation_base_error" in ep
            assert "validation_v3_error" in ep
            assert "validation_improvement" in ep
            assert "final_test_base_error" in ep
            assert "final_test_v3_error" in ep
            assert "final_test_improvement" in ep
    finally:
        db.close()


def test_v3_production_routing_on_new_project_compounds():
    """Verify runtime prediction routing and Project Adapter on newly registered compounds in a separate test project."""
    db = SessionLocal()
    try:
        # Create a separate test project
        test_proj_name = "V3_Production_Validation_Project"
        proj = db.scalar(select(Project).where(Project.name == test_proj_name))
        if not proj:
            proj = Project(
                name=test_proj_name,
                target="EGFR_KRAS",
                indication="Oncology Non-Small Cell Lung Cancer",
                description="Test project for Global Engine v3.0 production validation",
            )
            db.add(proj)
            db.commit()
            db.refresh(proj)

        # Register 3 new test compounds in the DB
        new_test_compounds = [
            ("Osimertinib", "CC(=O)Nc1cc(Nc2ncc(Cl)c(Nc3ccccc3)n2)c(OC)cc1N(C)CCN(C)C"),
            ("Sotorasib", "C=CC(=O)N1CCN(c2nc(F)c(F)c(-c3c(C(C)C)n[nH]c3-c3c(F)ccc(F)c3)c2F)CC1"),
            ("Lapatinib", "CS(=O)(=O)CCNCC1=CC=C(O1)C2=CC3=C(C=C2)N=CN=C3NC4=CC(=C(C=C4)OCC5=CC(=CC=C5)F)Cl"),
        ]

        for name, smi in new_test_compounds:
            comp = db.scalar(select(Compound).where(Compound.project_id == proj.id, Compound.name == name))
            if not comp:
                comp = Compound(
                    project_id=proj.id,
                    compound_id=f"CMPD-{name.upper()}",
                    name=name,
                    status="ACTIVE",
                    current_version=1,
                )
                db.add(comp)
                db.commit()
                db.refresh(comp)

                mol = Chem.MolFromSmiles(smi)
                canon_smi = Chem.MolToSmiles(mol, canonical=True) if mol else smi
                inchi_str = Chem.MolToInchi(mol) if mol else ""
                inchikey_str = Chem.MolToInchiKey(mol) if mol else ""

                cv = CompoundVersion(
                    compound_row_id=comp.id,
                    version_number=1,
                    original_smiles=smi,
                    canonical_smiles=canon_smi,
                    isomeric_smiles=canon_smi,
                    inchi=inchi_str,
                    inchikey=inchikey_str,
                    change_note="Initial version for v3 validation",
                )
                db.add(cv)
                db.commit()
                db.refresh(cv)

            # 1. CYP3A4 (GLOBAL_V3_PRIMARY) -> routes to Global v3
            res_cyp3a4 = predict_global_v3_endpoint(db, smi, "CYP3A4_INHIBITION", project_id=proj.id)
            assert res_cyp3a4["model_tier"] == "GLOBAL_V3_PRIMARY"
            assert res_cyp3a4["production_prediction"] == res_cyp3a4["v3_prediction"]
            assert res_cyp3a4["engine_version"] == "global-prediction-engine-v3.0.0"

            # 2. CYP2D6 (GLOBAL_V3_PRIMARY) -> routes to Global v3
            res_cyp2d6 = predict_global_v3_endpoint(db, smi, "CYP2D6_INHIBITION", project_id=proj.id)
            assert res_cyp2d6["model_tier"] == "GLOBAL_V3_PRIMARY"
            assert res_cyp2d6["production_prediction"] == res_cyp2d6["v3_prediction"]

            # 3. Solubility (GLOBAL_V3_PRIMARY) -> routes to Global v3
            res_sol = predict_global_v3_endpoint(db, smi, "SOLUBILITY_GENERIC", project_id=proj.id)
            assert res_sol["model_tier"] == "GLOBAL_V3_PRIMARY"
            assert res_sol["production_prediction"] == res_sol["v3_prediction"]

            # 4. PPB (V3_CANDIDATE) -> routes safely to Base Production
            res_ppb = predict_global_v3_endpoint(db, smi, "HUMAN_PPB", project_id=proj.id)
            assert res_ppb["model_tier"] == "BASE_PRODUCTION"
            assert res_ppb["production_prediction"] == res_ppb["base_prediction"]

            # 5. hERG (RETAIN_BASE) -> routes safely to Base Production
            res_herg = predict_global_v3_endpoint(db, smi, "HERG_LIABILITY", project_id=proj.id)
            assert res_herg["model_tier"] == "BASE_PRODUCTION"
            assert res_herg["production_prediction"] == res_herg["base_prediction"]
    finally:
        db.close()
