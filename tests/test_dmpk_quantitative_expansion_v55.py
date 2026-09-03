"""Unit and integration tests for Quantitative DMPK Prediction Expansion v5.5."""
import pytest
from backend.endpoint_model_validation import build_dmpk_quantitative_expansion_report, audit_evidence_funnel


def test_dmpk_quantitative_report_schema():
    """Verify DMPK quantitative expansion report adheres to the required schema."""
    report = build_dmpk_quantitative_expansion_report()
    assert isinstance(report, list)
    assert len(report) >= 12

    for row in report:
        assert "endpoint" in row
        assert "n" in row
        assert "existing_classifier" in row
        assert "quantitative_model" in row
        assert "mae_rmse" in row
        assert "coverage" in row
        assert "ood" in row
        assert "status" in row
        assert row["coverage"] == "100.0%"
        assert row["ood"] >= 0


def test_no_fabricated_models_for_unavailable_endpoints():
    """Verify unavailable endpoints remain MODEL_UNAVAILABLE and not fabricated."""
    report = build_dmpk_quantitative_expansion_report()
    cyp2c19_row = next((r for r in report if r["endpoint"] == "CYP2C19 quantitative inhibition"), None)
    assert cyp2c19_row is not None
    assert cyp2c19_row["quantitative_model"] == "MODEL_UNAVAILABLE"

    herg_row = next((r for r in report if r["endpoint"] == "hERG liability"), None)
    assert herg_row is not None
    assert herg_row["quantitative_model"] == "MODEL_UNAVAILABLE"


def test_reconciliation_audit_zero_evidence_loss():
    """Verify 7-record reconciliation confirms 1,364 total qualified records and silent_loss = 0."""
    funnel = audit_evidence_funnel()
    rec = funnel.get("reconciliation_7_records", {})
    assert rec.get("extracted_records") == 1368
    assert rec.get("classified_records") == 1364
    assert rec.get("silent_evidence_loss") == 0
    assert rec.get("reconciliation_status") == "FULL_RECONCILIATION_VERIFIED"
