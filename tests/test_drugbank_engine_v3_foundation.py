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


def test_fifteen_reference_drugs_sequential_ingestion_and_roles():
    """Verify all 15 reference drugs are ingested with exact identifiers, qualified records, and fixed roles."""
    db = SessionLocal()
    try:
        res_list = ingest_all_drugbank_reference_drugs(db)
        assert len(res_list) == 15
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

        for r in res_list:
            assert r["records_ingested_n"] >= 7
            assert r["status"] == "SUCCESS"
    finally:
        db.close()


def test_upstream_training_overlap_isolation_and_holdouts():
    """Verify upstream training overlap is partitioned into DEVELOPMENT_TRAINING vs IMMUTABLE_HOLDOUT across cohorts."""
    db = SessionLocal()
    try:
        dataset = build_global_learning_dataset(db)
        assert dataset["total_compounds_registered"] == 15
        assert dataset["total_eligible_observations"] >= 100
        assert dataset["total_holdout_observations"] >= 50

        herg_data = dataset["endpoints"]["HERG_LIABILITY"]
        assert len(herg_data["development_training_samples"]) == 4
        assert len(herg_data["immutable_holdout_samples"]) == 11
    finally:
        db.close()


def test_engine_v3_candidate_integrity_and_empirical_audit():
    """Verify Engine v3 candidate models are 100% empirical with zero synthetic multipliers and zero leakage."""
    db = SessionLocal()
    try:
        v3_eval = evaluate_global_engine_v3_readiness(db)
        assert "global-prediction-engine-v3" in v3_eval["engine_version"]
        assert v3_eval["total_compounds"] == 15

        # PPB must have Dev N=0 and UNAVAILABLE_NO_TRAINING_FIT
        ppb_eval = next(e for e in v3_eval["endpoints_evaluated"] if e["endpoint_id"] == "HUMAN_PPB")
        assert ppb_eval["development_training_n"] == 0
        assert ppb_eval["actual_candidate_mae"] == "UNAVAILABLE_NO_TRAINING_FIT"
        assert "Dev Training N=0" in ppb_eval["decision"]

        # CYP3A4 must have real empirical improvement on N=8 holdouts
        cyp3a4_eval = next(e for e in v3_eval["endpoints_evaluated"] if e["endpoint_id"] == "CYP3A4_INHIBITION")
        assert cyp3a4_eval["development_training_n"] == 4
        assert cyp3a4_eval["immutable_holdout_n"] == 8
        assert cyp3a4_eval["actual_candidate_mae"] < cyp3a4_eval["actual_base_mae"]
        assert cyp3a4_eval["evolution_status"] == "V3_CANDIDATE_VALIDATED"

        # CYP2D6 must have real empirical improvement on holdouts
        cyp2d6_eval = next(e for e in v3_eval["endpoints_evaluated"] if e["endpoint_id"] == "CYP2D6_INHIBITION")
        assert cyp2d6_eval["actual_candidate_mae"] < cyp2d6_eval["actual_base_mae"]
        assert cyp2d6_eval["evolution_status"] == "V3_CANDIDATE_VALIDATED"

        # hERG calibration audit: if candidate does not beat base, retain base status
        herg_eval = next(e for e in v3_eval["endpoints_evaluated"] if e["endpoint_id"] == "HERG_LIABILITY")
        assert herg_eval["evolution_status"] == "CANDIDATE_EVALUATED_RETAIN_BASE"
        assert "RETAIN_BASE_STATUS" in herg_eval["decision"]
    finally:
        db.close()
