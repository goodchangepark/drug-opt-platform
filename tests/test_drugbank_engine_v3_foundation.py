"""
Tests for Global Prediction Engine v3.1 Release & 50 DrugBank Reference Library (Stage 6 / v3.1 Expansion & Project Adapter Governance).
"""
import pytest
from sqlalchemy import select

from backend.database import SessionLocal
from backend.drugbank_reference import (
    ensure_drugbank_project,
    ingest_gefitinib_reference_drug,
    ingest_all_drugbank_reference_drugs,
    ingest_v3_1_expansion_drugs_sequential,
    DRUGBANK_PROJECT_NAME,
    REFERENCE_DRUGS_CATALOG,
    ROLE_DEVELOPMENT_TRAINING,
    ROLE_MODEL_SELECTION_VALIDATION,
    ROLE_FINAL_TEST_COHORT_1_CONSUMED,
    ROLE_FINAL_TEST_COHORT_2_CONSUMED,
    ROLE_LOCKED_FINAL_TEST_COHORT_3,
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


def test_fifty_reference_drugs_sequential_ingestion_and_roles():
    """Verify all 50 reference drugs are ingested with exact identifiers, qualified records, and fixed roles."""
    db = SessionLocal()
    try:
        res_list = ingest_all_drugbank_reference_drugs(db)
        assert len(res_list) == 50
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
        # 10 New approved reference drugs in v3.1
        assert "Amlodipine" in names
        assert "Losartan" in names
        assert "Metronidazole" in names
        assert "Montelukast" in names
        assert "Pantoprazole" in names
        assert "Raloxifene" in names
        assert "Tamoxifen" in names
        assert "Theophylline" in names
        assert "Tolbutamide" in names
        assert "Trazodone" in names

        for r in res_list:
            assert r["records_ingested_n"] >= 5
            assert r["status"] == "SUCCESS"
    finally:
        db.close()


def test_v3_1_ten_drugs_sequential_lifecycle_execution():
    """
    Verify Directive 2:
    10 approved reference drugs (Drugs 41-50) are ingested one-by-one with:
    1. Identity (distinct scaffold from original 40)
    2. Evidence (prioritized: PPB -> hERG -> Caco-2 -> HLM -> CYPs)
    3. Qualification (EXACT_MATCH, QUALIFIED_FOR_GLOBAL_TRAINING)
    4. Prediction (computed for each qualified endpoint)
    5. Error (calculated before advancing to next drug)
    """
    db = SessionLocal()
    try:
        results = ingest_v3_1_expansion_drugs_sequential(db)
        assert len(results) == 10
        expected_names = [
            "Amlodipine", "Losartan", "Metronidazole", "Montelukast", "Pantoprazole",
            "Raloxifene", "Tamoxifen", "Theophylline", "Tolbutamide", "Trazodone"
        ]
        for r, exp_name in zip(results, expected_names):
            assert r["compound_name"] == exp_name
            assert r["identity"]["status"] == "IDENTITY_VERIFIED"
            assert r["identity"]["scaffold_family"] != ""
            assert len(r["evidence"]) >= 5
            assert len(r["qualification"]) >= 5
            assert len(r["prediction"]) >= 5
            assert len(r["error"]) >= 5
            # Verify endpoint priority ordering: PPB appears before hERG, etc.
            endpoint_order = [e["endpoint_id"] for e in r["evidence"]]
            if "HUMAN_PPB" in endpoint_order and "HERG_LIABILITY" in endpoint_order:
                assert endpoint_order.index("HUMAN_PPB") < endpoint_order.index("HERG_LIABILITY")
            # Verify all error calculations are numeric
            for err in r["error"]:
                assert err["absolute_error"] is not None
                assert err["absolute_error"] >= 0.0
    finally:
        db.close()


def test_five_tier_partitioning_and_locked_final_test_cohorts():
    """Verify observations are partitioned into Dev (N=21), Validation (N=18), Consumed Cohort 1 & 2 (N=6), and Locked Final Test Cohort 3 (N=5)."""
    db = SessionLocal()
    try:
        dataset = build_global_learning_dataset(db)
        assert dataset["total_compounds_registered"] == 50
        assert dataset["total_eligible_observations"] >= 250
        assert dataset["total_development_observations"] >= 100
        assert dataset["total_validation_observations"] >= 90
        assert dataset["total_consumed_observations"] >= 25
        assert dataset["total_final_test_observations"] >= 20

        # Check CYP3A4 5-tier split
        cyp3a4_data = dataset["endpoints"]["CYP3A4_INHIBITION"]
        assert len(cyp3a4_data["development_training_samples"]) >= 15
        assert len(cyp3a4_data["model_selection_validation_samples"]) >= 10
        assert len(cyp3a4_data["final_test_consumed_samples"]) >= 5
        assert len(cyp3a4_data["locked_final_test_samples"]) == 5
    finally:
        db.close()


def test_engine_v3_release_readiness_and_runtime_routing():
    """Verify Engine v3.1 readiness evaluates candidates, promotes qualified endpoints, and executes runtime routing."""
    db = SessionLocal()
    try:
        v3_eval = evaluate_global_engine_v3_readiness(db)
        assert "global-prediction-engine-v3.1" in v3_eval["engine_version"]
        assert v3_eval["total_compounds"] == 50

        # CYP3A4 Promotion
        cyp3a4_eval = next(e for e in v3_eval["endpoints_evaluated"] if e["endpoint_id"] == "CYP3A4_INHIBITION")
        assert cyp3a4_eval["promotion_status"] == "GLOBAL_V3_PRIMARY"
        assert cyp3a4_eval["v3_error"] < cyp3a4_eval["base_error"]
        assert cyp3a4_eval["final_test_v3_error"] < cyp3a4_eval["final_test_base_error"]

        # CYP2D6 Promotion
        cyp2d6_eval = next(e for e in v3_eval["endpoints_evaluated"] if e["endpoint_id"] == "CYP2D6_INHIBITION")
        assert cyp2d6_eval["promotion_status"] == "GLOBAL_V3_PRIMARY"
        assert cyp2d6_eval["v3_error"] < cyp2d6_eval["base_error"]

        # hERG Promoted in v3.1 based on empirical holdout replication
        herg_eval = next(e for e in v3_eval["endpoints_evaluated"] if e["endpoint_id"] == "HERG_LIABILITY")
        assert herg_eval["promotion_status"] == "GLOBAL_V3_PRIMARY"

        # PPB Candidate Status
        ppb_eval = next(e for e in v3_eval["endpoints_evaluated"] if e["endpoint_id"] == "HUMAN_PPB")
        assert ppb_eval["promotion_status"] == "V3_CANDIDATE"

        # Test Runtime Prediction Routing for CYP3A4 (GLOBAL_V3_PRIMARY)
        test_smi = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1" # Gefitinib
        route_cyp3a4 = predict_global_v3_endpoint(db, test_smi, "CYP3A4_INHIBITION")
        assert route_cyp3a4["model_tier"] == "GLOBAL_V3_PRIMARY"
        assert route_cyp3a4["production_prediction"] != route_cyp3a4["base_prediction"]
        assert "v3-" in route_cyp3a4["model_version_hash"]
        assert route_cyp3a4["global_prediction"] is not None
        assert route_cyp3a4["project_adjusted_prediction"] is None

        # Test Runtime Prediction Routing for hERG (GLOBAL_V3_PRIMARY)
        route_herg = predict_global_v3_endpoint(db, test_smi, "HERG_LIABILITY")
        assert route_herg["model_tier"] == "GLOBAL_V3_PRIMARY"
        assert route_herg["production_prediction"] == route_herg["global_prediction"]
        assert route_herg["project_adjusted_prediction"] is None

        # Test Runtime Prediction Routing for PPB (BASE_PRODUCTION)
        route_ppb = predict_global_v3_endpoint(db, test_smi, "HUMAN_PPB")
        assert route_ppb["model_tier"] == "BASE_PRODUCTION"
        assert route_ppb["production_prediction"] == route_ppb["base_prediction"]
        assert route_ppb["project_adjusted_prediction"] is None
    finally:
        db.close()
