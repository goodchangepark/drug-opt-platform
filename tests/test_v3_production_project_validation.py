"""
End-to-End Validation for Global Prediction Engine v3.1 Production Release and Project Adapter Governance.
"""
import pytest
from rdkit import Chem
from sqlalchemy import select, delete

from backend.database import SessionLocal
from backend.models import Project, Compound, CompoundVersion, ExternalExperimentalEvidence
from backend.admet import ADMETPredictionRun
from backend.engine_v3_learning import (
    evaluate_global_engine_v3_readiness,
    predict_global_v3_endpoint,
    build_global_learning_dataset,
)


def test_v3_production_release_governance_and_promotion_statuses():
    """Verify that CYP3A4, CYP2D6, Solubility, hERG, HLM, CYP1A2, and CYP2C9 are GLOBAL_V3_PRIMARY; Caco-2 is V3_CANDIDATE; PPB is RETAIN_BASE."""
    db = SessionLocal()
    try:
        readiness = evaluate_global_engine_v3_readiness(db)
        assert readiness["release_status"] == "GLOBAL_ENGINE_V3_3_PRODUCTION_RELEASE"

        statuses = {ep["endpoint_id"]: ep["promotion_status"] for ep in readiness["endpoints_evaluated"]}
        assert statuses["CYP3A4_INHIBITION"] == "GLOBAL_V3_PRIMARY"
        assert statuses["CYP2D6_INHIBITION"] == "GLOBAL_V3_PRIMARY"
        assert statuses["SOLUBILITY_GENERIC"] == "GLOBAL_V3_PRIMARY"
        assert statuses["HERG_LIABILITY"] == "GLOBAL_V3_PRIMARY"
        assert statuses["HLM_INTRINSIC_CLEARANCE"] == "GLOBAL_V3_PRIMARY"
        assert statuses["CYP1A2_INHIBITION"] == "GLOBAL_V3_PRIMARY"
        assert statuses["CYP2C9_INHIBITION"] == "GLOBAL_V3_PRIMARY"
        assert statuses["CACO2_PERMEABILITY"] == "V3_CANDIDATE"
        assert statuses["HUMAN_PPB"] == "RETAIN_BASE"
        assert statuses["CYP2C19_INHIBITION"] == "MODEL_UNAVAILABLE"

        # Verify separated Validation vs Final-Test performance metrics
        for ep in readiness["endpoints_evaluated"]:
            assert "validation_base_error" in ep
            assert "validation_v3_error" in ep
            assert "validation_improvement" in ep
            assert "final_test_base_error" in ep
            assert "final_test_v3_error" in ep
            assert "final_test_improvement" in ep
            assert "prospective_metrics" in ep
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
                description="Test project for Global Engine v3.3 production validation",
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
            assert res_cyp3a4["engine_version"] == "global-prediction-engine-v3.3.0"
            assert res_cyp3a4["global_prediction"] == res_cyp3a4["v3_prediction"]
            assert res_cyp3a4["project_adjusted_prediction"] is None
            assert res_cyp3a4["project_adapter_status"] in ("INSUFFICIENT_DATA", "OUT_OF_DOMAIN_DISABLED")
            assert "prediction_uncertainty" in res_cyp3a4
            assert "descriptor_envelope" in res_cyp3a4
            assert "ad_extrapolation_guard_applied" in res_cyp3a4

            # 2. CYP2D6 (GLOBAL_V3_PRIMARY) -> routes to Global v3
            res_cyp2d6 = predict_global_v3_endpoint(db, smi, "CYP2D6_INHIBITION", project_id=proj.id)
            assert res_cyp2d6["model_tier"] == "GLOBAL_V3_PRIMARY"
            assert res_cyp2d6["production_prediction"] == res_cyp2d6["v3_prediction"]

            # 3. Solubility (GLOBAL_V3_PRIMARY) -> routes to Global v3
            res_sol = predict_global_v3_endpoint(db, smi, "SOLUBILITY_GENERIC", project_id=proj.id)
            assert res_sol["model_tier"] == "GLOBAL_V3_PRIMARY"
            assert res_sol["production_prediction"] == res_sol["v3_prediction"]

            # 4. PPB (RETAIN_BASE) -> routes safely to Base Production
            res_ppb = predict_global_v3_endpoint(db, smi, "HUMAN_PPB", project_id=proj.id)
            assert res_ppb["model_tier"] == "BASE_PRODUCTION"
            assert res_ppb["production_prediction"] == res_ppb["base_prediction"]

            # 5. hERG (GLOBAL_V3_PRIMARY) -> routes to Global v3
            res_herg = predict_global_v3_endpoint(db, smi, "HERG_LIABILITY", project_id=proj.id)
            assert res_herg["model_tier"] == "GLOBAL_V3_PRIMARY"
            assert res_herg["production_prediction"] == res_herg["global_prediction"]

            # 6. CYP1A2 (GLOBAL_V3_PRIMARY in v3.3) -> routes to Global v3
            res_cyp1a2 = predict_global_v3_endpoint(db, smi, "CYP1A2_INHIBITION", project_id=proj.id)
            assert res_cyp1a2["model_tier"] == "GLOBAL_V3_PRIMARY"
            assert res_cyp1a2["production_prediction"] == res_cyp1a2["v3_prediction"]

            # 7. CYP2C9 (GLOBAL_V3_PRIMARY in v3.3) -> routes to Global v3
            res_cyp2c9 = predict_global_v3_endpoint(db, smi, "CYP2C9_INHIBITION", project_id=proj.id)
            assert res_cyp2c9["model_tier"] == "GLOBAL_V3_PRIMARY"
            assert res_cyp2c9["production_prediction"] == res_cyp2c9["v3_prediction"]
    finally:
        db.close()


def test_project_adapter_independent_compound_governance():
    """
    Verify strict Project Adapter governance:
    1. Independent compound count N < 5 -> INSUFFICIENT_DATA, adapter inactive, Global/Base preserved.
    2. Independent compound count N >= 5 -> LOCO CV evaluated. If CV MAE improves, status = ACTIVE_ADAPTED.
    3. Global and Project-adjusted predictions are returned separately (None when unadapted).
    4. Zero leakage: Project data does NOT pollute the global DrugBank library.
    """
    db = SessionLocal()
    try:
        # Ensure clean state for test project
        gov_proj_name = "Project_Adapter_Governance_Test"
        existing_proj = db.scalar(select(Project).where(Project.name == gov_proj_name))
        if existing_proj:
            for c in db.scalars(select(Compound).where(Compound.project_id == existing_proj.id)).all():
                for cv in db.scalars(select(CompoundVersion).where(CompoundVersion.compound_row_id == c.id)).all():
                    db.execute(delete(ExternalExperimentalEvidence).where(ExternalExperimentalEvidence.compound_version_id == cv.id))
                    db.delete(cv)
                db.delete(c)
            db.delete(existing_proj)
            db.commit()

        proj = Project(
            name=gov_proj_name,
            target="TEST_TARGET",
            indication="Benchmarking Governance",
            description="Test project for Project Adapter K>=5 threshold and LOCO CV verification",
        )
        db.add(proj)
        db.commit()
        db.refresh(proj)

        # 5 distinct compounds with known structures
        adapter_test_specs = [
            ("Adapter_C1", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", 5.2),
            ("Adapter_C2", "CC1=C(C(=O)C2=C(C1=O)N3CC4=C(C=CC=C4)C3=C2)C", 5.8),
            ("Adapter_C3", "Clc1ccc(cc1)C(c2ccccc2)N3CCNCC3", 6.1),
            ("Adapter_C4", "CN1C(=O)CN=C(c2ccccc2)c3cc(Cl)ccc13", 5.5),
            ("Adapter_C5", "CCN(CC)CCNC(=O)c1c(C)[nH]c(C=C2C(=O)Nc3ccc(F)cc23)c1C", 6.4),
        ]

        # Phase 1: Ingest only 3 compounds (N=3 < 5)
        for name, smi, exp_val in adapter_test_specs[:3]:
            comp = db.scalar(select(Compound).where(Compound.project_id == proj.id, Compound.name == name))
            if not comp:
                comp = Compound(project_id=proj.id, compound_id=f"GOV-{name}", name=name, status="ACTIVE", current_version=1)
                db.add(comp)
                db.commit()
                db.refresh(comp)

                mol = Chem.MolFromSmiles(smi)
                canon_smi = Chem.MolToSmiles(mol, canonical=True)
                cv = CompoundVersion(
                    compound_row_id=comp.id,
                    version_number=1,
                    original_smiles=smi,
                    canonical_smiles=canon_smi,
                    isomeric_smiles=canon_smi,
                    inchi=Chem.MolToInchi(mol),
                    inchikey=Chem.MolToInchiKey(mol),
                    change_note="Governance test compound",
                )
                db.add(cv)
                db.commit()
                db.refresh(cv)

                p_key = f"proj_gov_{cv.id}_cyp3a4_{exp_val}"
                ev = db.scalar(select(ExternalExperimentalEvidence).where(ExternalExperimentalEvidence.provenance_key == p_key))
                if not ev:
                    ev = ExternalExperimentalEvidence(
                        compound_version_id=cv.id,
                        provenance_key=p_key,
                        source_database="PROJECT_INHOUSE_ASSAY",
                        source_record_id=f"ASSAY-{name}",
                        canonical_endpoint_id="CYP3A4_INHIBITION",
                        raw_endpoint_name="CYP3A4 pIC50",
                        raw_value=str(exp_val),
                        raw_unit="pIC50",
                        normalized_value=str(exp_val),
                        normalized_unit="pIC50",
                        identity_match_status="EXACT_MATCH",
                        endpoint_match_status="EXACT_MATCH",
                        qualification_status="QUALIFIED_FOR_GLOBAL_TRAINING",
                        reference_text="In-house Project Assay 2026",
                        assay_conditions_json={"assay": "in_house_luminescent_cyp3a4"},
                    )
                    db.add(ev)
                    db.commit()

        # Check Phase 1: N=3 < 5 -> Must be INSUFFICIENT_DATA
        res_phase1 = predict_global_v3_endpoint(db, adapter_test_specs[0][1], "CYP3A4_INHIBITION", project_id=proj.id)
        assert res_phase1["project_adapter_status"] == "INSUFFICIENT_DATA"
        assert res_phase1["project_adapted"] is False
        assert res_phase1["project_adjusted_prediction"] is None
        assert res_phase1["production_prediction"] == res_phase1["global_prediction"]
        assert res_phase1["project_compound_n"] == 3

        # Phase 2: Add compounds 4 and 5 (N=5 >= 5)
        for name, smi, exp_val in adapter_test_specs[3:]:
            comp = db.scalar(select(Compound).where(Compound.project_id == proj.id, Compound.name == name))
            if not comp:
                comp = Compound(project_id=proj.id, compound_id=f"GOV-{name}", name=name, status="ACTIVE", current_version=1)
                db.add(comp)
                db.commit()
                db.refresh(comp)

                mol = Chem.MolFromSmiles(smi)
                canon_smi = Chem.MolToSmiles(mol, canonical=True)
                cv = CompoundVersion(
                    compound_row_id=comp.id,
                    version_number=1,
                    original_smiles=smi,
                    canonical_smiles=canon_smi,
                    isomeric_smiles=canon_smi,
                    inchi=Chem.MolToInchi(mol),
                    inchikey=Chem.MolToInchiKey(mol),
                    change_note="Governance test compound",
                )
                db.add(cv)
                db.commit()
                db.refresh(cv)

                p_key = f"proj_gov_{cv.id}_cyp3a4_{exp_val}"
                ev = db.scalar(select(ExternalExperimentalEvidence).where(ExternalExperimentalEvidence.provenance_key == p_key))
                if not ev:
                    ev = ExternalExperimentalEvidence(
                        compound_version_id=cv.id,
                        provenance_key=p_key,
                        source_database="PROJECT_INHOUSE_ASSAY",
                        source_record_id=f"ASSAY-{name}",
                        canonical_endpoint_id="CYP3A4_INHIBITION",
                        raw_endpoint_name="CYP3A4 pIC50",
                        raw_value=str(exp_val),
                        raw_unit="pIC50",
                        normalized_value=str(exp_val),
                        normalized_unit="pIC50",
                        identity_match_status="EXACT_MATCH",
                        endpoint_match_status="EXACT_MATCH",
                        qualification_status="QUALIFIED_FOR_GLOBAL_TRAINING",
                        reference_text="In-house Project Assay 2026",
                        assay_conditions_json={"assay": "in_house_luminescent_cyp3a4"},
                    )
                    db.add(ev)
                    db.commit()

        # Check Phase 2: N=5 >= 5 -> LOCO CV evaluated
        res_phase2 = predict_global_v3_endpoint(db, adapter_test_specs[0][1], "CYP3A4_INHIBITION", project_id=proj.id)
        assert res_phase2["project_compound_n"] == 5
        assert res_phase2["project_adapter_status"] in ("ACTIVE_ADAPTED", "EVALUATED_NOT_IMPROVED")
        assert "global_prediction" in res_phase2

        if res_phase2["project_adapter_status"] == "ACTIVE_ADAPTED":
            assert res_phase2["project_adapted"] is True
            assert res_phase2["project_adjusted_prediction"] is not None
            assert res_phase2["production_prediction"] == res_phase2["project_adjusted_prediction"]
        else:
            assert res_phase2["project_adapted"] is False
            assert res_phase2["project_adjusted_prediction"] is None
            assert res_phase2["production_prediction"] == res_phase2["global_prediction"]

        # Phase 3: Zero Leakage Verification
        # DrugBank reference dataset must strictly contain only the 80 reference compounds
        db_summary = build_global_learning_dataset(db)
        assert db_summary["total_compounds_registered"] == 80
        assert db_summary["project_name"] == "DrugBank"
        # None of the adapter test compounds appear in DrugBank dataset
        for ep_key, ep_val in db_summary["endpoints"].items():
            for sample in ep_val.get("development_training_samples", []):
                assert not sample["compound_name"].startswith("Adapter_C")
            for sample in ep_val.get("model_selection_validation_samples", []):
                assert not sample["compound_name"].startswith("Adapter_C")
    finally:
        db.close()


def test_v3_3_applicability_domain_extrapolation_guard_and_uncertainty():
    """
    Verify Directive 9:
    1. IN_DOMAIN compound receives full calibration adjustment and lower uncertainty.
    2. OUT_OF_DOMAIN compound triggers AD extrapolation guard, falling back safely to base prediction.
    3. OUT_OF_DOMAIN compound disables Project Adapter automatically with status OUT_OF_DOMAIN_DISABLED.
    4. Descriptor envelope and calibration residual distribution are tracked in output provenance.
    """
    db = SessionLocal()
    try:
        # Standard drug-like compound: Gefitinib (IN_DOMAIN)
        in_domain_smi = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"
        res_in = predict_global_v3_endpoint(db, in_domain_smi, "CYP3A4_INHIBITION")
        assert res_in["applicability_domain"] == "IN_DOMAIN"
        assert res_in["ad_extrapolation_guard_applied"] is False
        assert res_in["prediction_uncertainty"] > 0.0
        assert "molecular_weight" in res_in["descriptor_envelope"]
        assert "residual_std" in res_in["calibration_residual_distribution"]

        # Giant synthetic polymer-like molecule (OUT_OF_DOMAIN: MW > 800)
        ood_smi = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O" # C68H136O2, MW ~993
        res_ood = predict_global_v3_endpoint(db, ood_smi, "CYP3A4_INHIBITION", project_id=1)
        assert res_ood["applicability_domain"] == "OUT_OF_DOMAIN"
        assert res_ood["ad_extrapolation_guard_applied"] is True
        # Guarded fallback: production prediction falls back to base prediction
        assert res_ood["production_prediction"] == res_ood["base_prediction"]
        # Project adapter must be disabled when OOD
        assert res_ood["project_adapter_status"] == "OUT_OF_DOMAIN_DISABLED"
        assert res_ood["project_adapted"] is False
        # Uncertainty is penalized in OOD
        assert res_ood["prediction_uncertainty"] > res_in["prediction_uncertainty"]
    finally:
        db.close()

