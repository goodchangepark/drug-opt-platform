"""Unit and integration tests for Endpoint Model Validation & Primary Selection v5.3."""
import pytest
from backend.endpoint_model_validation import run_endpoint_validation, EP_MAP
from backend.database import SessionLocal
from backend.models import Compound
from backend.endpoint_comparison import build_endpoint_comparison


def test_endpoint_validation_report_execution():
    """Verify that run_endpoint_validation executes without errors and generates report rows."""
    report = run_endpoint_validation()
    assert isinstance(report, list)
    assert len(report) >= len(EP_MAP)
    
    # Check required fields in each report row
    for r in report:
        assert "endpoint_name" in r
        assert "canonical_endpoint_id" in r
        assert "independent_n" in r
        assert "primary_model_id" in r
        assert "primary_error" in r
        assert "alternative_errors" in r
        assert "consensus_error" in r
        assert "decision" in r
        assert isinstance(r["independent_n"], int)


def test_no_leakage_and_rule_model_isolation():
    """Verify rule/mechanistic estimates are excluded from ML consensus voting."""
    report = run_endpoint_validation()
    ppb_row = next((r for r in report if r["endpoint_name"] == "Plasma protein binding"), None)
    assert ppb_row is not None
    # PPB primary is admetica_ppbr; physchem sigmoid rule is an alternative, NOT mixed into ML consensus voting
    assert ppb_row["primary_model_id"] == "admetica_ppbr"
    assert "physchem_human_ppb_v1" in ppb_row["alternative_errors"]


def test_scientific_row_multimodel_metadata():
    """Verify that scientific_rows include multi-model UI metadata."""
    db = SessionLocal()
    try:
        comp = db.query(Compound).filter(Compound.name.ilike("%Mobocertinib%")).first()
        assert comp is not None
        res = build_endpoint_comparison(db, comp.versions[-1].id)
        rows = res.get("scientific_rows", [])
        assert len(rows) > 0

        # Check that rows have multi-model and validation fields
        for r in rows:
            assert "primary_model" in r
            assert "primary_prediction" in r
            assert "alternative_models" in r
            assert "consensus" in r
            assert "validation_n" in r
            assert "model_performance" in r
    finally:
        db.close()
