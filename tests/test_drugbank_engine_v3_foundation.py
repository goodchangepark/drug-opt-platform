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


def test_five_reference_drugs_sequential_ingestion():
    """Verify all 5 reference drugs are ingested with exact identifiers and qualified records."""
    db = SessionLocal()
    try:
        res_list = ingest_all_drugbank_reference_drugs(db)
        assert len(res_list) == 5
        names = [r["compound_name"] for r in res_list]
        assert "Gefitinib" in names
        assert "Imatinib" in names
        assert "Propranolol" in names
        assert "Atorvastatin" in names
        assert "Midazolam" in names

        for r in res_list:
            assert r["records_ingested_n"] >= 8
            assert r["status"] == "SUCCESS"
    finally:
        db.close()


def test_upstream_training_overlap_isolation_and_holdouts():
    """Verify upstream training overlap is partitioned into TRAINING_ELIGIBLE vs VALIDATION_HOLDOUT."""
    db = SessionLocal()
    try:
        dataset = build_global_learning_dataset(db)
        assert dataset["total_compounds_registered"] == 5
        assert dataset["total_eligible_observations"] >= 30
        assert dataset["total_holdout_observations"] >= 20

        ppb_data = dataset["endpoints"]["HUMAN_PPB"]
        # PPB compounds have exact structure overlap in upstream datasets
        assert len(ppb_data["training_eligible_samples"]) == 5
        assert len(ppb_data["validation_holdout_samples"]) == 0

        cyp3a4_data = dataset["endpoints"]["CYP3A4_INHIBITION"]
        assert len(cyp3a4_data["validation_holdout_samples"]) >= 3
    finally:
        db.close()


def test_engine_v3_readiness_audit_zero_unproven_claims():
    """Verify Engine v3 readiness computes real holdout MAE and prohibits unproven v3 claims when N < 5."""
    db = SessionLocal()
    try:
        v3_eval = evaluate_global_engine_v3_readiness(db)
        assert "global-prediction-engine-v3" in v3_eval["engine_version"]
        assert v3_eval["total_compounds"] == 5

        for ep in v3_eval["endpoints_evaluated"]:
            if ep["true_holdout_n"] < 5:
                assert ep["fine_tuned_v3_mae"] == "PENDING_SUFFICIENT_HOLDOUT_N"
                assert "Promotion Gated" in ep["decision"]
                assert "NO_IMPROVEMENT_CLAIMED" in ep["projected_improvement"]
            else:
                assert ep["decision"] == "V3_CALIBRATION_READY_FOR_VALIDATION"
                assert isinstance(ep["fine_tuned_v3_mae"], float)
    finally:
        db.close()
