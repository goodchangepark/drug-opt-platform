"""
Tests for Global Prediction Engine v3.0 Release Candidate & 40 DrugBank Reference Library (Stage 6 / v3.0.0 Global Completion).
"""
import pytest
from sqlalchemy import select

from backend.database import SessionLocal
from backend.drugbank_reference import (
    ensure_drugbank_project,
    ingest_gefitinib_reference_drug,
    ingest_all_drugbank_reference_drugs,
    DRUGBANK_PROJECT_NAME,
    REFERENCE_DRUGS_CATALOG,
    ROLE_DEVELOPMENT_TRAINING,
    ROLE_MODEL_SELECTION_VALIDATION,
    ROLE_FINAL_TEST_COHORT_1_CONSUMED,
    ROLE_LOCKED_FINAL_TEST_COHORT_2,
)
from backend.engine_v3_learning import (
    build_global_learning_dataset,
    evaluate_global_engine_v3_readiness,
    predict_global_v3_endpoint,
)
from backend.models import Project, Compound, CompoundVersion, ExternalExperimentalEvidence


def test_drugbank_project_creation_and_provenance():
    """Verify DrugBank canonical project exists with GLOBAL_MODEL_DEVELOPMENT indication."""
    db = SessionLocal()
    try:
        proj = ensure_drugbank_project(db)
        assert proj.name == DRUGBANK_PROJECT_NAME
        assert "GLOBAL_MODEL_DEVELOPMENT" in proj.indication
    finally:
        db.close()


def test_forty_reference_drugs_sequential_ingestion_and_roles():
    """Verify all 40 reference drugs are ingested with exact identifiers, qualified records, and fixed roles."""
    db = SessionLocal()
    try:
        res_list = ingest_all_drugbank_reference_drugs(db)
        assert len(res_list) == 40
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
        assert "Bupropion" in names
        assert "Carvedilol" in names
        assert "Clopidogrel" in names
        assert "Diltiazem" in names
        assert "Erythromycin" in names
        assert "Flecainide" in names
        assert "Lansoprazole" in names
        assert "Nifedipine" in names
        assert "Omeprazole" in names
        assert "Simvastatin" in names
        assert "Celecoxib" in names
        assert "Diazepam" in names
        assert "Diclofenac" in names
        assert "Indomethacin" in names
        assert "Warfarin" in names
        assert "Atenolol" in names
        assert "Caffeine" in names
        assert "Ibuprofen" in names
        assert "Lorcaserin" in names
        assert "Rosuvastatin" in names

        for r in res_list:
            assert r["records_ingested_n"] >= 5
            assert r["status"] == "SUCCESS"
    finally:
        db.close()


def test_four_tier_partitioning_and_locked_final_test_cohort():
    """Verify observations are partitioned into Dev (N=18), Validation (N=16), Consumed (N=1), and Locked Final Test (N=5)."""
    db = SessionLocal()
    try:
        dataset = build_global_learning_dataset(db)
        assert dataset["total_compounds_registered"] == 40
        assert dataset["total_eligible_observations"] >= 200
        assert dataset["total_development_observations"] >= 80
        assert dataset["total_validation_observations"] >= 80
        assert dataset["total_consumed_observations"] >= 5
        assert dataset["total_final_test_observations"] >= 20

        # Check CYP3A4 4-tier split
        cyp3a4_data = dataset["endpoints"]["CYP3A4_INHIBITION"]
        assert len(cyp3a4_data["development_training_samples"]) >= 15
        assert len(cyp3a4_data["model_selection_validation_samples"]) >= 10
        assert len(cyp3a4_data["final_test_consumed_samples"]) == 1
        assert len(cyp3a4_data["locked_final_test_samples"]) == 5
    finally:
        db.close()


def test_engine_v3_release_readiness_and_runtime_routing():
    """Verify Engine v3 readiness evaluates candidates, promotes qualified endpoints, and executes runtime routing."""
    db = SessionLocal()
    try:
        v3_eval = evaluate_global_engine_v3_readiness(db)
        assert "global-prediction-engine-v3" in v3_eval["engine_version"]
        assert v3_eval["total_compounds"] == 40

        # CYP3A4 Promotion
        cyp3a4_eval = next(e for e in v3_eval["endpoints_evaluated"] if e["endpoint_id"] == "CYP3A4_INHIBITION")
        assert cyp3a4_eval["promotion_status"] == "GLOBAL_V3_PRIMARY"
        assert cyp3a4_eval["v3_error"] < cyp3a4_eval["base_error"]
        assert cyp3a4_eval["final_test_v3_error"] < cyp3a4_eval["final_test_base_error"]

        # CYP2D6 Promotion
        cyp2d6_eval = next(e for e in v3_eval["endpoints_evaluated"] if e["endpoint_id"] == "CYP2D6_INHIBITION")
        assert cyp2d6_eval["promotion_status"] == "GLOBAL_V3_PRIMARY"
        assert cyp2d6_eval["v3_error"] < cyp2d6_eval["base_error"]

        # hERG Retain Base Production
        herg_eval = next(e for e in v3_eval["endpoints_evaluated"] if e["endpoint_id"] == "HERG_LIABILITY")
        assert herg_eval["promotion_status"] == "RETAIN_BASE"

        # PPB Candidate Status
        ppb_eval = next(e for e in v3_eval["endpoints_evaluated"] if e["endpoint_id"] == "HUMAN_PPB")
        assert ppb_eval["promotion_status"] == "V3_CANDIDATE"

        # Test Runtime Prediction Routing for CYP3A4 (GLOBAL_V3_PRIMARY)
        test_smi = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1" # Gefitinib
        route_cyp3a4 = predict_global_v3_endpoint(db, test_smi, "CYP3A4_INHIBITION")
        assert route_cyp3a4["model_tier"] == "GLOBAL_V3_PRIMARY"
        assert route_cyp3a4["production_prediction"] != route_cyp3a4["base_prediction"]
        assert "v3-" in route_cyp3a4["model_version_hash"]

        # Test Runtime Prediction Routing for hERG (BASE_PRODUCTION)
        route_herg = predict_global_v3_endpoint(db, test_smi, "HERG_LIABILITY")
        assert route_herg["model_tier"] == "BASE_PRODUCTION"
        assert route_herg["production_prediction"] == route_herg["base_prediction"]
    finally:
        db.close()
