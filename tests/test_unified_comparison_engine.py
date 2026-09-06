"""
Tests for Unified Experimental <-> Prediction Comparison Engine
Validates:
- Canonical Endpoint Registry (95 endpoints)
- Deterministic Scientific Unit Normalization
- Comparison Matcher (10 evaluation dimensions)
- Zero pairable evidence loss (PAIRABLE_EVIDENCE_WITHOUT_PAIR = 0)
- Sunvozertinib Full Comparison Recovery
"""
import pytest
from backend.canonical_endpoints import (
    REGISTRY, CANONICAL_ENDPOINT_VERSION, DIRECT, CONVERTED, RELATED, UNSUPPORTED,
    canonicalize_prediction_endpoint, normalize_experimental_observation
)
from backend.unit_normalization import (
    convert_concentration, convert_solubility, convert_ppb,
    convert_clearance, convert_volume, convert_time, convert_caco2_papp, convert_auc,
    clean_unit_str, normalize_scientific_observation
)
from backend.comparison_matcher import (
    ScientificComparisonPair, ComparisonPairMatcher, compute_comparison_coverage
)
from backend.database import SessionLocal
from backend.endpoint_comparison import build_endpoint_comparison


def test_canonical_endpoint_registry_completeness():
    """Verify registry has at least 95 canonical endpoints across ADMET, PK, ACTIVITY, TOX."""
    assert len(REGISTRY) >= 95
    sections = {ep.section for ep in REGISTRY.values()}
    assert {"ADMET", "PK", "METABOLISM", "TOXICITY", "PROPERTIES", "ACTIVITY"}.issubset(sections)
    
    # Check specific critical endpoints
    assert "HUMAN_PK_AUC_ORAL" in REGISTRY
    assert "HUMAN_PK_CMAX_ORAL" in REGISTRY
    assert "HUMAN_PPB" in REGISTRY
    assert "CACO2_PAPP_AB" in REGISTRY
    assert "HLM_CLINT" in REGISTRY


def test_unit_normalization_concentration():
    """Test concentration conversions between molar and mass units."""
    mw = 604.7  # Sunvozertinib MW
    # 60 nM to ng/mL
    res = convert_concentration(60.0, "nM", "ng/mL", mw=mw)
    assert res.normalized_unit == "ng/mL"
    expected = 60.0 * 1e-9 * mw * 1e6  # = 36.282 ng/mL
    assert abs(res.normalized_value - expected) < 1e-3

    # Direct identity
    res2 = convert_concentration(100.0, "ng/mL", "ng/mL")
    assert res2.normalized_value == 100.0


def test_unit_normalization_ppb():
    """Test PPB (% bound <-> fu)."""
    res = convert_ppb(91.46, "% bound", "fu")
    assert abs(res.normalized_value - 0.0854) < 1e-4

    res2 = convert_ppb(0.0854, "fu", "% bound")
    assert abs(res2.normalized_value - 91.46) < 1e-4


def test_unit_normalization_clearance_and_volume():
    """Test clearance and volume conversions."""
    # Clearance mL/min/kg -> L/h/kg
    cl_res = convert_clearance(10.0, "mL/min/kg", "L/h/kg")
    assert abs(cl_res.normalized_value - 0.6) < 1e-4

    # Volume L/kg -> mL/kg
    v_res = convert_volume(1.5, "L/kg", "mL/kg")
    assert abs(v_res.normalized_value - 1500.0) < 1e-4


def test_unit_normalization_auc():
    """Test AUC dot/exponent syntax and standard units."""
    # 1783 ng.hr.mL-1 -> ng*h/mL
    auc_res = convert_auc(1783.0, "ng.hr.mL-1", "ng*h/mL")
    assert auc_res.normalized_value == 1783.0
    assert auc_res.normalized_unit == "ng*h/mL"

    # mg*h/L -> ng*h/mL (factor 1000)
    auc_res2 = convert_auc(5.0, "mg*h/L", "ng*h/mL")
    assert abs(auc_res2.normalized_value - 5000.0) < 1e-3


def test_comparison_matcher_direct_and_converted():
    """Test matcher classifying direct vs converted vs related pairs."""
    matcher = ComparisonPairMatcher(mw=500.0)
    
    # Direct pair
    exp_rec = {
        "canonical_endpoint_id": "HUMAN_PK_AUC_ORAL",
        "species": "HUMAN",
        "route": "ORAL",
        "raw_value": 8000.0,
        "raw_unit": "ng*h/mL",
        "dose": 200.0,
        "dose_unit": "mg"
    }
    pred_rec = {
        "available": True,
        "value": 7200.0,
        "unit": "ng*h/mL",
        "species": "HUMAN",
        "route": "ORAL",
        "dose": 200.0,
        "dose_unit": "mg"
    }
    pair = matcher.match_pair(exp_rec, pred_rec)
    assert pair.comparison_status == DIRECT
    assert pair.fold_error is not None
    assert pair.fold_error <= 2.0

    # Cross-species mismatch should be contextually related, never directly compared
    pred_rat = dict(pred_rec, species="RAT")
    pair_cross = matcher.match_pair(exp_rec, pred_rat)
    assert pair_cross.comparison_status == RELATED
    assert pair_cross.absolute_error is None


def test_sunvozertinib_zero_pairable_loss():
    """Verify Sunvozertinib (compound version 13) achieves PAIRABLE_EVIDENCE_WITHOUT_PAIR = 0."""
    db = SessionLocal()
    comp = build_endpoint_comparison(db, 13)
    db.close()
    
    cov = comp.get("summary", {}).get("comparison_coverage", {})
    assert cov.get("pairable_evidence_without_pair") == 0
    assert cov.get("successfully_paired") >= 3
    assert cov.get("direct_pairs") >= 1
    assert cov.get("converted_pairs") >= 2
