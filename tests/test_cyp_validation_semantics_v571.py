"""Unit and integration tests for CYP Validation Semantics Fix v5.7.1."""
import pytest
from backend.openadmet_cyp import (
    classify_cyp_assay_context,
    predict_chemeleon_cyp_pic50,
    CONTEXT_MATCHED_RECOMBINANT,
    CONTEXT_RELATED_HLM_DIRECT,
    CONTEXT_RELATED_TDI,
    CONTEXT_RELATED_SCREENING_LIMIT,
)
from backend.endpoint_model_validation import (
    build_cyp_assay_matched_validation_table,
    audit_cyp_quantitative_validation,
)


def test_recombinant_vs_hlm_context_classification():
    """Verify classification separates recombinant enzyme from HLM direct and TDI."""
    # 1. Recombinant human enzyme
    ctx_rec, r_rec, is_rec, is_hlm = classify_cyp_assay_context(
        raw_endpoint="CYP3A4_INHIBITION",
        raw_value=0.10,
        raw_unit="µM",
        assay_matrix="rhCYP",
        reference_text="rhCYP recombinant enzyme direct assay",
    )
    assert ctx_rec == CONTEXT_MATCHED_RECOMBINANT
    assert is_rec is True
    assert is_hlm is False

    # 2. HLM direct inhibition
    ctx_hlm, r_hlm, is_rec_h, is_hlm_h = classify_cyp_assay_context(
        raw_endpoint="CYP3A4_INHIBITION",
        raw_value=0.17,
        raw_unit="µM",
        assay_matrix="HLM",
        reference_text="Human liver microsomes direct substrate assay",
    )
    assert ctx_hlm == CONTEXT_RELATED_HLM_DIRECT
    assert is_rec_h is False
    assert is_hlm_h is True

    # 3. TDI (HLM + NADPH pre-incubation)
    ctx_tdi, r_tdi, is_rec_t, is_hlm_t = classify_cyp_assay_context(
        raw_endpoint="CYP3A4_INHIBITION",
        raw_value=0.0073,
        raw_unit="µM",
        assay_matrix="HLM",
        reference_text="30-min pre-incubation with NADPH (TDI shift)",
    )
    assert ctx_tdi == CONTEXT_RELATED_TDI
    assert is_rec_t is False
    assert is_hlm_t is False


def test_two_tier_validation_report_isolation():
    """Verify audit outputs separate Tier 1 Recombinant and Tier 2 HLM metrics without cross-contamination."""
    audit = audit_cyp_quantitative_validation()
    assert audit["audit_version"] == "CYP_VALIDATION_SEMANTICS_V571"
    assert audit["validation_nomenclature"] == "RETROSPECTIVE_EXTERNAL_VALIDATION"

    cyp3a4 = audit["isoforms"]["CYP3A4"]
    tier1 = cyp3a4["tier1_recombinant_validation"]
    tier2 = cyp3a4["tier2_hlm_related_validation"]

    # Tier 1 contains only recombinant observations (Poziotinib, Sunvozertinib)
    assert tier1["independent_n"] == 2
    assert tier1["mae_pic50"] < 0.20
    assert tier1["ad_breakdown"]["out_of_domain"] == 0

    # Tier 2 contains HLM observation (Orforglipron)
    assert tier2["independent_n"] == 1
    assert tier2["ad_breakdown"]["out_of_domain"] == 1

    # Promotion remains strictly prohibited
    assert "Strictly Prohibited" in cyp3a4["promotion_decision"] or "RETAIN_CANDIDATE_STATUS" in cyp3a4["promotion_decision"]
