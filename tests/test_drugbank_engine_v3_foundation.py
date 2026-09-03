"""Unit and integration tests for DrugBank Project & Global Engine v3.0 Foundation."""
import pytest
from backend.database import SessionLocal
from backend.drugbank_reference import ensure_drugbank_project, ingest_gefitinib_reference_drug, DRUGBANK_PROJECT_NAME
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


def test_gefitinib_reference_drug_ingestion_and_identifiers():
    """Verify Gefitinib reference drug is ingested with exact identifiers and qualified records."""
    db = SessionLocal()
    try:
        res = ingest_gefitinib_reference_drug(db)
        assert res["status"] == "SUCCESS"
        assert res["compound_name"] == "Gefitinib"
        assert res["drugbank_id"] == "DB00317"
        assert res["records_ingested_n"] >= 10

        comp = db.scalar(select(Compound).where(Compound.id == res["compound_id"]))
        assert comp is not None
        assert "184475-35-2" in comp.cas_number
        assert "CHEMBL939" in comp.notes
        assert "123631" in comp.notes

        cv = db.scalar(select(CompoundVersion).where(CompoundVersion.compound_row_id == comp.id, CompoundVersion.version_number == 1))
        assert cv is not None
        assert cv.inchikey == "XGALLCVXEZPNRQ-UHFFFAOYSA-N"
    finally:
        db.close()


def test_gefitinib_base_predictions_and_error_evaluation():
    """Verify base predictions and signed/absolute errors are computed for Gefitinib across endpoints."""
    db = SessionLocal()
    try:
        res = ingest_gefitinib_reference_drug(db)
        evals = res["evaluations"]

        sol_eval = next(e for e in evals if e["canonical_endpoint_id"] == "SOLUBILITY_GENERIC")
        assert sol_eval["base_prediction"] is not None
        assert sol_eval["absolute_error"] < 0.5
        assert sol_eval["global_training_eligible"] is True

        cyp3a4_eval = next(e for e in evals if e["canonical_endpoint_id"] == "CYP3A4_INHIBITION")
        assert cyp3a4_eval["base_prediction"] is not None
        assert cyp3a4_eval["global_training_eligible"] is True

        pk_eval = next(e for e in evals if e["canonical_endpoint_id"] == "HUMAN_PK_CMAX_ORAL")
        assert pk_eval["global_training_eligible"] is False  # Composite PK parameter
    finally:
        db.close()


def test_engine_v3_readiness_and_promotion_gating():
    """Verify Engine v3 dataset accumulation and promotion gating (N < 5 blocks promotion)."""
    db = SessionLocal()
    try:
        v3_eval = evaluate_global_engine_v3_readiness(db)
        assert "global-prediction-engine-v3" in v3_eval["engine_version"]
        assert v3_eval["total_eligible_observations"] >= 5

        for ep in v3_eval["endpoints_evaluated"]:
            assert ep["drugbank_reference_n"] >= 1
            assert "Promotion Gated" in ep["decision"] or ep["drugbank_reference_n"] >= 5
    finally:
        db.close()
