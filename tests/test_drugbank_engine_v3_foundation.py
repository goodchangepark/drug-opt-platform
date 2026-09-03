"""Unit and integration tests for DrugBank Project & Global Engine v3.0 Foundation."""
import pytest
from backend.database import SessionLocal
from backend.drugbank_reference import (
    ensure_drugbank_project,
    ingest_gefitinib_reference_drug,
    ingest_all_drugbank_reference_drugs,
    DRUGBANK_PROJECT_NAME,
)
from backend.engine_v3_learning import build_global_learning_dataset, evaluate_global_engine_v3_readiness
from backend.models import Project, Compound, CompoundVersion, ExternalExperimentalEvidence
from sqlalchemy import select


def test_drugbank_project_creation():
    """Verify DrugBank project is created with GLOBAL_MODEL_DEVELOPMENT designation."""
    db = SessionLocal()
    try:
        proj = ensure_drugbank_project(db)
        assert proj.name == DRUGBANK_PROJECT_NAME
        assert "GLOBAL_MODEL_DEVELOPMENT" in proj.indication or "GLOBAL_MODEL_DEVELOPMENT" in proj.description
    finally:
        db.close()


def test_twenty_reference_drugs_sequential_ingestion_and_roles():
    """Verify all 20 reference drugs are ingested with exact identifiers, qualified records, and fixed roles."""
    db = SessionLocal()
    try:
        res_list = ingest_all_drugbank_reference_drugs(db)
        assert len(res_list) == 20
        names = [r["compound_name"] for r in res_list]
        assert "Gefitinib" in names
        assert "Imatinib" in names
        assert "Propranolol" in names
        assert "Atorvastatin" in names
        assert "Midazolam" in names
        assert "Verapamil" in names
        assert "Fluoxetine" in names
        assert "Ketoconazole" in names
        assert "Sildenafil" in names
        assert "Quinidine" in names
        assert "Dextromethorphan" in names
        assert "Amiodarone" in names
        assert "Clarithromycin" in names
        assert "Duloxetine" in names
        assert "Haloperidol" in names
        assert "Paroxetine" in names
        assert "Metoprolol" in names
        assert "Terbinafine" in names
        assert "Ritonavir" in names
        assert "Cimetidine" in names

        for r in res_list:
            assert r["records_ingested_n"] >= 5
            assert r["status"] == "SUCCESS"
    finally:
        db.close()


def test_three_way_partitioning_and_locked_final_test():
    """Verify observations are partitioned into DEVELOPMENT_TRAINING, MODEL_SELECTION_VALIDATION, and LOCKED_FINAL_TEST."""
    db = SessionLocal()
    try:
        dataset = build_global_learning_dataset(db)
        assert dataset["total_compounds_registered"] == 20
        assert dataset["total_eligible_observations"] >= 120
        assert dataset["total_development_observations"] >= 35
        assert dataset["total_validation_observations"] >= 50
        assert dataset["total_final_test_observations"] >= 5

        # Check CYP3A4 3-way split
        cyp3a4_data = dataset["endpoints"]["CYP3A4_INHIBITION"]
        assert len(cyp3a4_data["development_training_samples"]) == 8
        assert len(cyp3a4_data["model_selection_validation_samples"]) == 8
        assert len(cyp3a4_data["locked_final_test_samples"]) == 1
    finally:
        db.close()


def test_engine_v3_candidate_integrity_and_locked_evaluation():
    """Verify Engine v3 candidate models are evaluated on Validation and Locked Final Test cohorts with Primary promotion gated."""
    db = SessionLocal()
    try:
        v3_eval = evaluate_global_engine_v3_readiness(db)
        assert "global-prediction-engine-v3" in v3_eval["engine_version"]
        assert v3_eval["total_compounds"] == 20

        # CYP3A4 validation & final test evaluation
        cyp3a4_eval = next(e for e in v3_eval["endpoints_evaluated"] if e["endpoint_id"] == "CYP3A4_INHIBITION")
        assert cyp3a4_eval["development_training_n"] == 8
        assert cyp3a4_eval["model_selection_validation_n"] == 8
        assert cyp3a4_eval["locked_final_test_n"] == 1
        assert cyp3a4_eval["actual_candidate_mae"] < cyp3a4_eval["actual_base_mae"]
        assert cyp3a4_eval["evolution_status"] == "V3_CANDIDATE_VALIDATED"
        assert "RETAIN_CANDIDATE_STATUS" in cyp3a4_eval["decision"]
        assert "Primary promotion gated" in cyp3a4_eval["decision"]

        # CYP2D6 validation & final test evaluation
        cyp2d6_eval = next(e for e in v3_eval["endpoints_evaluated"] if e["endpoint_id"] == "CYP2D6_INHIBITION")
        assert cyp2d6_eval["development_training_n"] == 6
        assert cyp2d6_eval["model_selection_validation_n"] == 6
        assert cyp2d6_eval["actual_candidate_mae"] < cyp2d6_eval["actual_base_mae"]
        assert cyp2d6_eval["evolution_status"] == "V3_CANDIDATE_VALIDATED"
        assert "Primary promotion gated" in cyp2d6_eval["decision"]

        # PPB validation
        ppb_eval = next(e for e in v3_eval["endpoints_evaluated"] if e["endpoint_id"] == "HUMAN_PPB")
        assert ppb_eval["development_training_n"] == 4
        assert ppb_eval["actual_candidate_mae"] < ppb_eval["actual_base_mae"]

        # hERG calibration audit: retain base model
        herg_eval = next(e for e in v3_eval["endpoints_evaluated"] if e["endpoint_id"] == "HERG_LIABILITY")
        assert herg_eval["evolution_status"] == "CANDIDATE_EVALUATED_RETAIN_BASE"
        assert "RETAIN_BASE_STATUS" in herg_eval["decision"]
    finally:
        db.close()
