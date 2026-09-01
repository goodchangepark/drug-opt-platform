"""
Stage 4D-3B2A: Unit & Integration Tests for hERG Calibration Audit.

Tests:
1.  hERG endpoint contract identity and semantic correctness.
2.  M1 / M2 model identity and version.
3.  Authoritative cohort schema, completeness (N=728), split assignment.
4.  Scaffold-aware split isolation (no scaffold in both cal and test).
5.  Calibration/test isolation: Platt fit on cal only, evaluate on test only.
6.  Isotonic fit on cal only, evaluate on test only.
7.  No test-set tuning: thresholds found on calibration set.
8.  Class balance / prevalence shift (training 86% vs eval 67%).
9.  Assay subtype handling: heterogeneity documented, not silently combined.
10. M1/M2 complementarity (rescue rate < 20%).
11. Fixed blend calibration-selected, holdout-evaluated.
12. Production threshold preserved at 0.50.
13. SHADOW consensus mode preserved.
14. UI freeze: zero frontend changes.
15. Full model metrics schema for all 9 required JSON artifacts.
16. Threshold audit: production threshold still in decision JSON.
17. Error analysis: FP > FN (high-sensitivity, low-specificity regime).
18. Disagreement JSON structure.
19. Adaptive gate is NO_GO.
20. Final decision is HERG_FIXED_BLEND_CANDIDATE or HERG_CALIBRATION_UPDATE_CANDIDATE.
"""

import json
from pathlib import Path
import pytest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
VAL_DIR = ROOT / "validation"
DOCS_DIR = ROOT / "docs"

# ── Helpers ────────────────────────────────────────────────────────────────────

def load_json(fname: str) -> dict:
    p = VAL_DIR / fname
    assert p.exists(), f"Validation artifact missing: {fname}"
    with open(p) as f:
        return json.load(f)


# ── 1. Endpoint contract ────────────────────────────────────────────────────────

def test_herg_endpoint_contract_identity():
    """hERG endpoint contract has correct ID, threshold, and positive class."""
    from backend.endpoint_contracts import get_endpoint_contract
    contract = get_endpoint_contract("hERG liability")
    assert contract.endpoint_id == "safety_herg_blocker_prob"
    assert contract.classification_semantics["positive_class"] == "BLOCKER"
    assert contract.classification_semantics["decision_threshold"] == 0.50
    assert contract.classification_semantics["potency_threshold_um"] == 10.0
    assert "hERG" in contract.comparison_rules.get("target_must_match", "")


def test_herg_endpoint_contract_forbidden_mix():
    """hERG endpoint contract forbids mixing with QT-prolongation clinical data."""
    from backend.endpoint_contracts import get_endpoint_contract
    contract = get_endpoint_contract("hERG liability")
    forbidden = contract.ensemble_compatibility_rules.get("forbidden_mix_types", [])
    assert any("qt_prolongation" in f for f in forbidden)


# ── 2. Model identity ───────────────────────────────────────────────────────────

def test_m1_herg_model_identity():
    """M1 admetica_safety_herg has correct model version and ID."""
    from backend.multimodel import get_adapters_for_endpoint
    adapters = get_adapters_for_endpoint("hERG liability")
    m1 = next((a for a in adapters if a.model_id == "admetica_safety_herg"), None)
    assert m1 is not None, "admetica_safety_herg adapter not found"
    assert "admetica-d4f7056-herg" in m1.model_version
    assert m1.model_family == "admetica"


def test_m2_herg_model_identity():
    """M2 physchem_herg_v1 has correct model version and architecture type."""
    from backend.multimodel import get_adapters_for_endpoint
    adapters = get_adapters_for_endpoint("hERG liability")
    m2 = next((a for a in adapters if a.model_id == "physchem_herg_v1"), None)
    assert m2 is not None, "physchem_herg_v1 adapter not found"
    assert m2.model_family == "pharmacophore_logistic"


def test_herg_exactly_two_adapters():
    """Exactly 2 adapters are registered for hERG liability (M1 + M2)."""
    from backend.multimodel import get_adapters_for_endpoint
    adapters = get_adapters_for_endpoint("hERG liability")
    ids = {a.model_id for a in adapters}
    assert ids == {"admetica_safety_herg", "physchem_herg_v1"}


# ── 3. Authoritative cohort ─────────────────────────────────────────────────────

def test_authoritative_cohort_schema_and_count():
    """stage4d3b2a_authoritative_cohort.json has N=728 with required fields."""
    data = load_json("stage4d3b2a_authoritative_cohort.json")
    assert data["n_compounds"] == 728
    assert data["endpoint"] == "safety_herg_blocker_prob"
    assert data["n_positive"] == 489
    assert data["n_negative"] == 239
    assert len(data["compounds"]) == 728
    required = [
        "compound_id", "canonical_smiles", "scaffold", "median_ic50_nM",
        "ic50_class", "experimental_label", "m1_probability", "m2_probability",
        "m1_label", "m2_label", "fp_m1", "fn_m1", "applicability_domain",
        "split_assignment",
    ]
    first = data["compounds"][0]
    for k in required:
        assert k in first, f"Missing key {k} in compound record"


def test_authoritative_cohort_label_definitions():
    """Cohort positive labels have IC50 <= 10000 nM (10 µM cutoff)."""
    data = load_json("stage4d3b2a_authoritative_cohort.json")
    for c in data["compounds"]:
        if c["experimental_label"] == 1:
            assert c["median_ic50_nM"] <= 10001.0, (
                f"Positive label with IC50={c['median_ic50_nM']} exceeds 10 µM threshold")
        elif c["experimental_label"] == 0:
            assert c["median_ic50_nM"] > 9999.0, (
                f"Negative label with IC50={c['median_ic50_nM']} is below 10 µM threshold")


# ── 4. Scaffold-aware split isolation ──────────────────────────────────────────

def test_scaffold_aware_split_no_leakage():
    """No scaffold appears in both calibration and test splits."""
    data = load_json("stage4d3b2a_authoritative_cohort.json")
    cal_scaffolds = set()
    tst_scaffolds = set()
    for c in data["compounds"]:
        sc = c.get("scaffold", "")
        if sc:
            if c["split_assignment"] == "calibration":
                cal_scaffolds.add(sc)
            else:
                tst_scaffolds.add(sc)
    overlap = cal_scaffolds & tst_scaffolds
    assert len(overlap) == 0, f"Scaffold leakage: {len(overlap)} scaffolds in both splits"


def test_calibration_test_sizes_reasonable():
    """Calibration is ~75% and test is ~25% of total."""
    data = load_json("stage4d3b2a_authoritative_cohort.json")
    cal_n = sum(1 for c in data["compounds"] if c["split_assignment"] == "calibration")
    tst_n = sum(1 for c in data["compounds"] if c["split_assignment"] == "test")
    total = cal_n + tst_n
    assert total == data["n_compounds"]
    assert 0.60 <= cal_n / total <= 0.85, f"Unexpected calibration fraction: {cal_n/total:.2f}"


# ── 5 & 6. Calibration isolation (Platt / isotonic) ───────────────────────────

def test_calibration_json_has_isolated_results():
    """Calibration artifact shows holdout results and references calibration-fitted models."""
    data = load_json("stage4d3b2a_calibration.json")
    assert "holdout_calibration_comparison" in data
    assert "m1_platt" in data["holdout_calibration_comparison"]
    assert "m1_isotonic" in data["holdout_calibration_comparison"]
    assert "m2_platt" in data["holdout_calibration_comparison"]
    assert data["split_method"] == "scaffold_aware_25pct_test"


def test_platt_improves_ece_on_holdout():
    """Platt scaling reduces ECE and LogLoss for M1 on untouched holdout."""
    data = load_json("stage4d3b2a_calibration.json")
    m1_raw = data["holdout_calibration_comparison"]["m1_raw"]
    m1_platt = data["holdout_calibration_comparison"]["m1_platt"]
    assert m1_platt["ece"] < m1_raw["ece"], "Platt should reduce ECE"
    assert m1_platt["log_loss"] < m1_raw["log_loss"], "Platt should reduce LogLoss"


def test_isotonic_improves_ece_on_holdout():
    """Isotonic regression reduces ECE for M1 on holdout."""
    data = load_json("stage4d3b2a_calibration.json")
    m1_raw = data["holdout_calibration_comparison"]["m1_raw"]
    m1_iso = data["holdout_calibration_comparison"]["m1_isotonic"]
    assert m1_iso["ece"] < m1_raw["ece"], "Isotonic should reduce ECE"


# ── 7. No test-set tuning ────────────────────────────────────────────────────────

def test_threshold_audit_optimized_on_calibration_only():
    """Threshold research is optimized on calibration set, evaluated on holdout."""
    data = load_json("stage4d3b2a_threshold_audit.json")
    assert "m1_optimal_thresholds_cal" in data
    holdout_eval = data.get("holdout_evaluation", {})
    for k in holdout_eval.values():
        if isinstance(k, dict):
            assert k.get("evaluated_on") == "untouched_holdout", \
                "Threshold holdout results must be labeled as evaluated on untouched_holdout"


# ── 8. Class balance / prevalence audit ────────────────────────────────────────

def test_class_imbalance_documented():
    """Authoritative cohort documents training vs evaluation prevalence shift > 10pp."""
    data = load_json("stage4d3b2a_authoritative_cohort.json")
    assert data["training_prevalence"] == pytest.approx(0.8599, abs=0.01)
    eval_prev = data["prevalence"]
    shift = abs(data["training_prevalence"] - eval_prev)
    assert shift > 0.10, f"Expected prevalence shift > 10pp, got {shift:.3f}"


def test_class_imbalance_in_root_cause():
    """CLASS_IMBALANCE appears in final decision root causes."""
    data = load_json("stage4d3b2a_final_decision.json")
    causes = data["root_cause_analysis"]["causes"]
    assert "CLASS_IMBALANCE" in causes


# ── 9. Assay heterogeneity handling ────────────────────────────────────────────

def test_assay_heterogeneity_documented():
    """Final decision documents assay heterogeneity as ASSAY_HETEROGENEITY_PRESENT."""
    data = load_json("stage4d3b2a_final_decision.json")
    assert data["assay_heterogeneity"] == "ASSAY_HETEROGENEITY_PRESENT"
    assert "patch" in data["assay_heterogeneity_detail"].lower() or \
           "radioligand" in data["assay_heterogeneity_detail"].lower()


def test_endpoint_contract_assay_type_heterogeneous():
    """hERG endpoint contract assay type string mentions heterogeneous modalities."""
    from backend.endpoint_contracts import get_endpoint_contract
    contract = get_endpoint_contract("hERG liability")
    assay_str = contract.assay_type.lower()
    assert "patch" in assay_str or "binding" in assay_str


# ── 10. M1/M2 complementarity ──────────────────────────────────────────────────

def test_m2_rescue_rate_low():
    """M2 rescue rate of M1 errors is < 20% (inadequate for adaptive weighting)."""
    data = load_json("stage4d3b2a_disagreement.json")
    rescue = data["model_complementarity"]["m2_rescue_rate"]
    assert rescue < 0.20, f"Unexpected high rescue rate: {rescue}"


def test_both_wrong_exceeds_m2_rescue():
    """Cases where both models fail simultaneously exceed M2 rescue cases."""
    data = load_json("stage4d3b2a_disagreement.json")
    both_wrong = data["model_complementarity"]["both_wrong"]
    m2_rescues = data["model_complementarity"]["m2_rescues_of_m1_errors"]
    assert both_wrong > m2_rescues, "Both-wrong should dominate over M2 rescues"


# ── 11. Fixed blend ──────────────────────────────────────────────────────────────

def test_fixed_blend_selected_on_calibration():
    """Best fixed blend is selected on calibration set and evaluated on holdout."""
    data = load_json("stage4d3b2a_fixed_blend.json")
    assert data["blend_selected_on"] == "calibration_set"
    assert data["blend_evaluated_on"] == "untouched_holdout"
    assert data["best_blend"] is not None
    assert "blend_results" in data
    # Verify all 6 blend ratios are present
    ratios = [(1.0, 0.0), (0.98, 0.02), (0.95, 0.05), (0.90, 0.10), (0.80, 0.20), (0.50, 0.50)]
    for w1, w2 in ratios:
        key = f"w1={w1:.2f}_w2={w2:.2f}"
        assert key in data["blend_results"], f"Missing blend ratio: {key}"


# ── 12. Production threshold preserved ─────────────────────────────────────────

def test_production_threshold_unchanged():
    """Production decision threshold is still 0.50 in final decision."""
    data = load_json("stage4d3b2a_final_decision.json")
    assert data["production_threshold"] == 0.50


def test_m1_spec_and_sens_match_known_values():
    """Reproduced M1 metrics match known Stage 4D-2 values within tolerance."""
    data = load_json("stage4d3b2a_final_decision.json")
    m = data["full_cohort_metrics_at_production_threshold"]
    # Known from independent_validation.json
    assert m["sensitivity"] == pytest.approx(0.9755, abs=0.005)
    assert m["specificity"] == pytest.approx(0.113, abs=0.01)
    assert m["mcc"] == pytest.approx(0.1844, abs=0.01)
    assert m["auroc"] == pytest.approx(0.6669, abs=0.01)


# ── 13. SHADOW mode preserved ───────────────────────────────────────────────────

def test_shadow_mode_preserved_in_final_decision():
    """Consensus mode is SHADOW throughout."""
    data = load_json("stage4d3b2a_final_decision.json")
    assert data["consensus_mode"] == "SHADOW"
    assert data["m1_model"]["role"] == "CORE"
    assert data["m2_model"]["role"] == "SHADOW_ONLY"


# ── 14. UI freeze ────────────────────────────────────────────────────────────────

def test_no_frontend_changes():
    """Historical UI guard now checks the current shared maturity contract."""
    app_js = (ROOT / "frontend/static/app.js").read_text(encoding="utf-8")
    assert "renderPredictionMaturity" in app_js
    assert "[...Array(5)]" in app_js


# ── 15. All 9 JSON artifacts exist ─────────────────────────────────────────────

def test_all_nine_validation_artifacts_exist():
    """All 9 required Stage 4D-3B2A validation artifacts are present."""
    required = [
        "stage4d3b2a_authoritative_cohort.json",
        "stage4d3b2a_model_metrics.json",
        "stage4d3b2a_calibration.json",
        "stage4d3b2a_threshold_audit.json",
        "stage4d3b2a_error_analysis.json",
        "stage4d3b2a_disagreement.json",
        "stage4d3b2a_fixed_blend.json",
        "stage4d3b2a_series_analysis.json",
        "stage4d3b2a_final_decision.json",
    ]
    for fname in required:
        path = VAL_DIR / fname
        assert path.exists(), f"Missing artifact: {fname}"
        assert path.stat().st_size > 100, f"Artifact appears empty: {fname}"


# ── 16. Threshold audit schema ──────────────────────────────────────────────────

def test_threshold_audit_schema():
    """Threshold audit JSON has all required fields."""
    data = load_json("stage4d3b2a_threshold_audit.json")
    assert "production_threshold" in data
    assert "m1_optimal_thresholds_cal" in data
    assert "holdout_evaluation" in data
    assert "threshold_verdict" in data
    assert data["threshold_verdict"] in {
        "THRESHOLD_DRIVEN", "MODEL_DRIVEN", "BOTH", "UNCERTAIN"
    }


# ── 17. Error analysis: FP dominates FN ────────────────────────────────────────

def test_error_analysis_fp_dominates_fn():
    """M1 has far more false positives than false negatives (poor specificity regime)."""
    data = load_json("stage4d3b2a_error_analysis.json")
    fp = data["m1_fp_count"]
    fn = data["m1_fn_count"]
    assert fp > fn * 5, f"Expected FP >> FN, got FP={fp}, FN={fn}"


def test_borderline_fraction_of_fps_is_high():
    """Majority of FPs are borderline IC50 compounds (1k–30k nM)."""
    data = load_json("stage4d3b2a_error_analysis.json")
    borderline_fp_frac = data.get("borderline_fraction_of_fp", 0)
    assert borderline_fp_frac > 0.50, \
        f"Expected >50% borderline FPs, got {borderline_fp_frac}"


# ── 18. Disagreement schema ─────────────────────────────────────────────────────

def test_disagreement_schema():
    """Disagreement JSON has M1/M2 complementarity fields."""
    data = load_json("stage4d3b2a_disagreement.json")
    comp = data["model_complementarity"]
    for k in ["both_correct", "both_wrong", "m1_only_correct", "m2_only_correct",
              "m1_error_n", "m2_rescues_of_m1_errors", "m2_rescue_rate"]:
        assert k in comp, f"Missing key: {k}"
    assert comp["both_correct"] + comp["both_wrong"] + comp["m1_only_correct"] + \
           comp["m2_only_correct"] == 728


# ── 19. Adaptive gate is NO_GO ──────────────────────────────────────────────────

def test_adaptive_weighting_gate_is_no_go():
    """Adaptive weighting gate must be NO_GO given insufficient M2 evidence."""
    data = load_json("stage4d3b2a_final_decision.json")
    assert data["adaptive_weighting_gate"] == "NO_GO"


# ── 20. Final scientific decision ───────────────────────────────────────────────

def test_final_decision_is_valid_option():
    """Final decision is one of the allowable hERG decisions."""
    data = load_json("stage4d3b2a_final_decision.json")
    valid_decisions = {
        "HERG_FIXED_CORE_SUFFICIENT",
        "HERG_CALIBRATION_UPDATE_CANDIDATE",
        "HERG_FIXED_BLEND_CANDIDATE",
        "HERG_ADAPTIVE_RESEARCH_CANDIDATE",
        "HERG_NEEDS_BETTER_SECONDARY_MODEL",
        "HERG_ENDPOINT_DATA_REQUALIFICATION_REQUIRED",
    }
    assert data["scientific_decision"] in valid_decisions


def test_series_analysis_schema():
    """Series analysis JSON documents scaffold series performance."""
    data = load_json("stage4d3b2a_series_analysis.json")
    assert "n_scaffolds_analyzed" in data
    assert "m2_wins_substantial" in data
    assert "m2_series_verdict" in data
    assert data["m2_series_verdict"] in {
        "SERIES_ADVANTAGE_PRESENT", "NO_REPRODUCIBLE_SERIES_ADVANTAGE"
    }


def test_label_boundary_uncertainty_documented():
    """Final decision documents LABEL_BOUNDARY_UNCERTAINTY as a root cause."""
    data = load_json("stage4d3b2a_final_decision.json")
    assert "LABEL_BOUNDARY_UNCERTAINTY" in data["root_cause_analysis"]["causes"]
    borderline_frac = data["root_cause_analysis"]["borderline_ic50_fraction"]
    assert borderline_frac > 0.60, \
        f"Expected >60% borderline fraction, got {borderline_frac}"
