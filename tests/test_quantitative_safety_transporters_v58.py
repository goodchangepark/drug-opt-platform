"""Unit and integration tests for Quantitative Safety & Transporter Prediction Framework v5.8."""
import pytest
from backend.quantitative_safety_transporters import (
    predict_quantitative_herg_pic50,
    predict_quantitative_pgp_pic50,
    build_quantitative_safety_transporter_validation_table,
    audit_safety_transporter_quantitative_validation,
    CONTEXT_MATCHED_PATCH_CLAMP,
    CONTEXT_RELATED_SCREENING_LIMIT,
)
from backend.endpoint_model_validation import build_dmpk_quantitative_expansion_report
from backend.database import SessionLocal
from backend.models import Compound, CompoundVersion


def test_herg_quantitative_candidate_model():
    """Verify hERG candidate regression model outputs valid numeric predictions and real AD."""
    # Pruvonertinib
    pruvo_smi = "C=CC(=O)Nc1cc(Nc2nccc(-c3cnc4c(C)cccn34)n2)c(OC)cc1N(C)CCN(C)C"
    pred = predict_quantitative_herg_pic50(pruvo_smi)
    assert pred.status == "CANDIDATE_EXTERNAL_MODEL"
    assert pred.pic50 == pytest.approx(6.07, abs=0.2)
    assert pred.ic50_um > 0.0
    assert pred.ic50_nm > 0.0
    assert pred.applicability_domain in ("IN_DOMAIN", "BORDERLINE", "OUT_OF_DOMAIN")


def test_pgp_quantitative_model_unavailable_governance():
    """Verify P-gp quantitative regression returns MODEL_UNAVAILABLE with zero fabricated models."""
    mobo_smi = "C=CC(=O)Nc1cc(Nc2ncc(C(=O)OC(C)C)c(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C"
    pred = predict_quantitative_pgp_pic50(mobo_smi)
    assert pred.status == "MODEL_UNAVAILABLE_PENDING_PRETRAINED_REGRESSION_CHECKPOINT"
    assert pred.pic50 is None
    assert pred.ic50_um is None
    assert "Broccatelli" in pred.ad_reason


def test_quantitative_safety_validation_table_and_audit():
    """Verify validation table schema and retrospective holdout metrics for safety endpoints."""
    audit = audit_safety_transporter_quantitative_validation()
    assert audit["audit_version"] == "QUANTITATIVE_SAFETY_TRANSPORTER_V58"

    herg_rep = audit["herg_quantitative"]
    assert herg_rep["status"] == "CANDIDATE_EXTERNAL_MODEL"
    assert "RETAIN_CANDIDATE_STATUS" in herg_rep["promotion_decision"]

    holdout = herg_rep["retrospective_external_holdout"]
    assert holdout["independent_n"] >= 3
    assert holdout["exact_overlap_n"] == 0

    rows = audit["table_rows"]
    pruvo_row = next(r for r in rows if r["compound"] == "Pruvonertinib" and r["target"] == "hERG")
    assert pruvo_row["context_match"] == CONTEXT_MATCHED_PATCH_CLAMP
    assert pruvo_row["eligible_for_mae"] is True

    pgp_row = next(r for r in rows if r["target"] == "P-gp")
    assert pgp_row["prediction"].startswith("MODEL_UNAVAILABLE")


def test_dmpk_expansion_report_includes_herg_and_pgp():
    """Verify DMPK quantitative expansion report reflects v5.8 hERG candidate and P-gp unavailable status."""
    report = build_dmpk_quantitative_expansion_report()
    herg_entry = next(r for r in report if r["endpoint"] == "hERG liability")
    assert herg_entry["status"] == "CANDIDATE_EXTERNAL_MODEL_EVALUATED"
    assert herg_entry["quantitative_model"] == "TDC CardioTox Chemprop hERG pIC50"
    assert herg_entry["n"] >= 1

    pgp_entry = next(r for r in report if r["endpoint"] == "P-gp inhibitor")
    assert pgp_entry["status"] == "MODEL_UNAVAILABLE_PENDING_PRETRAINED_REGRESSION_CHECKPOINT"
    assert pgp_entry["quantitative_model"] == "MODEL_UNAVAILABLE"
