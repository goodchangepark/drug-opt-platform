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


def test_ten_reference_drugs_sequential_ingestion_and_roles():
    """Verify all 10 reference drugs are ingested with exact identifiers, qualified records, and fixed roles."""
    db = SessionLocal()
    try:
        res_list = ingest_all_drugbank_reference_drugs(db)
        assert len(res_list) == 10
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

        for r in res_list:
            assert r["records_ingested_n"] >= 8
            assert r["status"] == "SUCCESS"
    finally:
        db.close()


def test_upstream_training_overlap_isolation_and_holdouts():
    """Verify upstream training overlap is partitioned into DEVELOPMENT_TRAINING vs IMMUTABLE_HOLDOUT."""
    db = SessionLocal()
    try:
        dataset = build_global_learning_dataset(db)
        assert dataset["total_compounds_registered"] == 10
        assert dataset["total_eligible_observations"] >= 70
        assert dataset["total_holdout_observations"] >= 25

        herg_data = dataset["endpoints"]["HERG_LIABILITY"]
        assert len(herg_data["development_training_samples"]) == 4
        assert len(herg_data["immutable_holdout_samples"]) == 6
    finally:
        db.close()


def test_engine_v3_readiness_audit_and_learning_curve():
    """Verify Engine v3 readiness evaluates immutable holdouts and computes learning curve snapshots."""
    db = SessionLocal()
    try:
        v3_eval = evaluate_global_engine_v3_readiness(db)
        assert "global-prediction-engine-v3" in v3_eval["engine_version"]
        assert v3_eval["total_compounds"] == 10
        assert len(v3_eval["herg_learning_curve"]) == 10

        herg_eval = next(e for e in v3_eval["endpoints_evaluated"] if e["endpoint_id"] == "HERG_LIABILITY")
        assert herg_eval["immutable_holdout_n"] == 6
        assert herg_eval["actual_base_mae"] < 0.50
        assert herg_eval["fine_tuned_v3_mae"] < herg_eval["actual_base_mae"]
        assert "RETAIN_CANDIDATE_STATUS" in herg_eval["decision"]
    finally:
        db.close()
