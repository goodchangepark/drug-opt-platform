"""Unit and integration tests for CYP Assay-Matched Validation v5.7."""
import pytest
from backend.openadmet_cyp import (
    classify_cyp_assay_context,
    predict_chemeleon_cyp_pic50,
    compute_fold_error,
    CONTEXT_MATCHED_DIRECT,
    CONTEXT_RELATED_TDI,
    CONTEXT_RELATED_SCREENING_LIMIT,
    CONTEXT_RELATED_HEPATOCYTE,
)
from backend.endpoint_model_validation import (
    build_cyp_assay_matched_validation_table,
    audit_cyp_quantitative_validation,
)


def test_assay_context_classifier():
    """Verify assay context classifier distinguishes direct inhibition, TDI, screening limits, and hepatocytes."""
    # 1. Direct Reversible Inhibition (rhCYP)
    ctx_dir, reason_dir, is_rec_d, is_hlm_d = classify_cyp_assay_context(
        raw_endpoint="CYP3A4_INHIBITION",
        raw_value=0.17,
        raw_unit="µM",
        assay_matrix="rhCYP",
        reference_text="Direct substrate inhibition in rhCYP",
    )
    assert ctx_dir == CONTEXT_MATCHED_DIRECT
    assert is_rec_d is True

    # 2. Time-Dependent Inhibition (TDI)
    ctx_tdi, reason_tdi, is_rec_t, is_hlm_t = classify_cyp_assay_context(
        raw_endpoint="CYP3A4_INHIBITION",
        raw_value=0.0073,
        raw_unit="µM",
        reference_text="30-min pre-incubation with NADPH (TDI shift)",
    )
    assert ctx_tdi == CONTEXT_RELATED_TDI
    assert is_rec_t is False

    # 3. Screening Threshold Limit
    ctx_scr, reason_scr, is_rec_s, is_hlm_s = classify_cyp_assay_context(
        raw_endpoint="CYP1A2_INHIBITION",
        raw_value=1.0,
        raw_relation=">",
        raw_unit="µM",
        reference_text="IC50 > 1 uM screening bound",
    )
    assert ctx_scr == CONTEXT_RELATED_SCREENING_LIMIT
    assert is_rec_s is False

    # 4. Hepatocyte Intact Cell Context
    ctx_hep, reason_hep, is_rec_h, is_hlm_h = classify_cyp_assay_context(
        raw_endpoint="CYP3A4_METABOLISM",
        raw_value=5.0,
        raw_unit="µM",
        assay_matrix="Hepatocyte",
        reference_text="Primary human hepatocyte intact clearance",
    )
    assert ctx_hep == CONTEXT_RELATED_HEPATOCYTE
    assert is_rec_h is False


def test_cyp_assay_matched_validation_table():
    """Verify assay-matched validation table for Poziotinib, Mobocertinib, Sunvozertinib, Orforglipron."""
    rows = build_cyp_assay_matched_validation_table()
    assert len(rows) >= 10

    # Poziotinib rhCYP direct assays should be IN_DOMAIN with fold errors < 1.5x
    pozio_3a4 = next(r for r in rows if r["compound"] == "Poziotinib" and "CYP3A4" in r["endpoint"])
    assert pozio_3a4["ad"] == "IN_DOMAIN"
    assert pozio_3a4["context_match"] == CONTEXT_MATCHED_DIRECT
    assert pozio_3a4["fold_error_float"] < 1.5

    # Sunvozertinib rhCYP direct assays should be BORDERLINE with fold errors < 1.5x
    sunvo_3a4 = next(r for r in rows if r["compound"] == "Sunvozertinib" and "CYP3A4" in r["endpoint"])
    assert sunvo_3a4["ad"] == "BORDERLINE"
    assert sunvo_3a4["context_match"] == CONTEXT_MATCHED_DIRECT
    assert sunvo_3a4["fold_error_float"] < 1.5

    # Orforglipron TDI should be tagged RELATED_CONTEXT_TDI and not eligible for reversible MAE
    orf_tdi = next(r for r in rows if r["compound"] == "Orforglipron" and "TDI" in r["assay"])
    assert orf_tdi["context_match"] == CONTEXT_RELATED_TDI
    assert orf_tdi["eligible_for_mae"] is False


def test_audit_reports_assay_matched_metrics():
    """Verify audit separates matched reversible holdout from related contexts and blocks promotion."""
    audit = audit_cyp_quantitative_validation()
    assert "CYP" in audit["audit_version"]

    for iso in ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]:
        rep = audit["isoforms"][iso]
        assert "RETAIN_CANDIDATE_STATUS" in rep["promotion_decision"]
        holdout = rep["assay_matched_holdout"]
        assert holdout["exact_overlap_n"] == 0
