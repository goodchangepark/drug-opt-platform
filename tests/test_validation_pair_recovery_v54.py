"""Unit and integration tests for Validation Pair Recovery v5.4."""
import pytest
from backend.endpoint_model_validation import audit_evidence_funnel, run_endpoint_validation, EP_MAP


def test_funnel_audit_structure():
    """Verify evidence funnel audits all records with proper drop-off attribution."""
    funnel = audit_evidence_funnel()
    assert "total_qualified_records" in funnel
    assert funnel["total_qualified_records"] >= 1300
    assert "endpoints" in funnel
    assert "global_drop_reasons" in funnel
    
    # Check drop reasons are populated
    reasons = funnel["global_drop_reasons"]
    assert "CLINICAL_PK_NO_PREDICTIVE_ML_MODEL" in reasons
    assert "CLINICAL_METABOLIC_BALANCE_NO_ML_REGRESSION" in reasons
    assert "TARGET_SPECIFIC_ACTIVITY_EXCLUDED_FROM_GENERAL_ADMET_ML" in reasons
    assert "QUANTITATIVE_EVIDENCE_NOT_CLASSIFICATION_PAIRABLE" in reasons


def test_solubility_scale_alignment():
    """Verify solubility scale alignment converts uM to log10(mol/L) and avoids scale mismatch."""
    report = run_endpoint_validation()
    sol_row = next((r for r in report if r["endpoint_name"] == "Solubility"), None)
    assert sol_row is not None
    assert sol_row["independent_n"] == 1
    # Error is in log units (< 2.0 log units, not 6.13)
    assert "MAE:" in sol_row["primary_error"]
    mae_val = float(sol_row["primary_error"].replace("MAE:", "").strip())
    assert mae_val < 2.0, f"Solubility MAE must be < 2.0 log units, got {mae_val}"


def test_no_arbitrary_ic50_binarization():
    """Verify quantitative IC50 is not arbitrarily converted to classifier labels."""
    funnel = audit_evidence_funnel()
    cyp3a4_funnel = funnel["endpoints"].get("CYP3A4_INHIBITION", {})
    assert "QUANTITATIVE_EVIDENCE_NOT_CLASSIFICATION_PAIRABLE" in cyp3a4_funnel.get("drop_reasons", {})
    
    # In validation report, CYP3A4 has no fabricated pairs
    report = run_endpoint_validation()
    cyp3a4_row = next((r for r in report if r["endpoint_name"] == "CYP3A4 inhibitor"), None)
    assert cyp3a4_row is not None
    assert cyp3a4_row["independent_n"] == 0
