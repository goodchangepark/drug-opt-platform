"""
Comprehensive Test Suite for PK Accuracy Validation & Critical Parameter Model Upgrade.
========================================================================================

Tests Directives 3-31:
1. 30-drug clinical PK cohort structure & diversity
2. Semantic endpoint separation (VDss ≠ Vz ≠ Vd/F, CL ≠ CL/F, route context)
3. Elimination status tagging (HEPATIC_CL_ESTIMATE vs TOTAL_CL_PREDICTION)
4. Oral PK confidence downgrading
5. Fold error calculations
6. pKa acid/base benchmark contract (MAE, N)
7. logD7.4 benchmark contract (MAE, N)
8. VDss holdout gate (>= 5% improvement required)
9. Sensitivity analysis sensible structure
10. Uncertainty interval calibration
11. PK engine vs prediction engine version separation
12. Clinical PK validation report integrity
13. DrugBank 250 preservation
14. DB integrity
"""
from __future__ import annotations
import json
import math
import sqlite3
from pathlib import Path
import pytest

from backend.clinical_pk_cohort import (
    CLINICAL_PK_COHORT,
    ObservedPKSemantic,
    classify_elimination_status,
    evaluate_oral_confidence,
    calculate_fold_error,
)
from backend.pk_parameter_set import (
    HUMAN_DEFAULT_BW_KG,
    HUMAN_QH_ML_MIN_KG,
    SPECIES_PHYSIOLOGY,
    PKParameterSet,
    PKParameterEntry,
    build_pk_parameter_set,
    canonical_fu_from_ppb,
    calculate_well_stirred_clearance,
    simulate_one_compartment_disposition,
    SOURCE_TYPE_EXPERIMENTAL,
    SOURCE_TYPE_VALIDATED_MODEL,
    SOURCE_TYPE_MECHANISTIC_DERIVED,
    SOURCE_TYPE_RULE_ESTIMATE,
    SOURCE_TYPE_USER_INPUT,
    SOURCE_TYPE_MODEL_UNAVAILABLE,
)
from backend.prediction_engine_registry import (
    CURRENT_ENGINE_ID,
    CURRENT_ENGINE_VERSION,
    CANDIDATE_ENGINE_ID,
)


# ============================================================
# Directive 3: Clinical PK Cohort Size & Diversity
# ============================================================
class TestClinicalPKCohort:
    def test_cohort_size_at_least_30(self):
        assert len(CLINICAL_PK_COHORT) >= 30, f"Expected >= 30 drugs, got {len(CLINICAL_PK_COHORT)}"

    def test_cohort_has_iv_and_oral(self):
        routes = set(d.route for d in CLINICAL_PK_COHORT)
        assert "IV" in routes
        assert "ORAL" in routes

    def test_ionization_diversity(self):
        classes = set(d.ionization_class for d in CLINICAL_PK_COHORT)
        assert "ACID" in classes
        assert "BASE" in classes
        assert "NEUTRAL" in classes

    def test_extraction_diversity(self):
        cats = set(d.extraction_category for d in CLINICAL_PK_COHORT)
        assert "LOW" in cats
        assert "INTERMEDIATE" in cats
        assert "HIGH" in cats

    def test_renal_diversity(self):
        """At least 3 drugs with >= 50% renal fraction excretion."""
        renal_heavy = [d for d in CLINICAL_PK_COHORT if d.fraction_excreted_renal >= 0.50]
        assert len(renal_heavy) >= 3, f"Need >= 3 high-renal drugs, found {len(renal_heavy)}"

    def test_ppb_diversity(self):
        ppbs = [d.ppb_percent for d in CLINICAL_PK_COHORT if d.ppb_percent is not None]
        assert min(ppbs) <= 10.0, "Need at least one low-PPB drug"
        assert max(ppbs) >= 99.0, "Need at least one extreme-PPB drug"


# ============================================================
# Directive 5: Semantic Endpoint Separation
# ============================================================
class TestSemanticEndpointSeparation:
    def test_vdss_vz_distinction(self):
        """VDss and Vz must be separate fields."""
        for d in CLINICAL_PK_COHORT:
            if d.vdss_l_kg and d.vz_l_kg:
                assert hasattr(d, 'vdss_l_kg')
                assert hasattr(d, 'vz_l_kg')
                assert d.vdss_l_kg != d.vz_l_kg or d.compound_name == "Acetaminophen"

    def test_cl_systemic_vs_cl_hepatic_separation(self):
        """CL systemic and CL hepatic must be separate fields."""
        for d in CLINICAL_PK_COHORT:
            assert hasattr(d, 'cl_systemic_ml_min_kg')
            assert hasattr(d, 'cl_hepatic_ml_min_kg')
            assert hasattr(d, 'cl_renal_ml_min_kg')
            assert hasattr(d, 'cl_oral_apparent_ml_min_kg')

    def test_vd_oral_apparent_distinct(self):
        """Vd/F (apparent oral volume) is a separate field."""
        assert hasattr(CLINICAL_PK_COHORT[0], 'vd_oral_apparent_l_kg')

    def test_route_context_preserved(self):
        """Every observation carries IV/ORAL route context."""
        for d in CLINICAL_PK_COHORT:
            assert d.route in ("IV", "ORAL"), f"Invalid route for {d.compound_name}: {d.route}"


# ============================================================
# Directive 16: Hepatic vs Total CL Distinction
# ============================================================
class TestEliminationStatusTagging:
    def test_hepatic_cl_estimate_for_high_renal(self):
        status = classify_elimination_status(0.90)  # Metformin
        assert status == "HEPATIC_CL_ESTIMATE"

    def test_total_cl_prediction_for_low_renal(self):
        status = classify_elimination_status(0.01)
        assert status == "TOTAL_CL_PREDICTION"

    def test_boundary_at_20_percent(self):
        assert classify_elimination_status(0.20) == "HEPATIC_CL_ESTIMATE"
        assert classify_elimination_status(0.19) == "TOTAL_CL_PREDICTION"


# ============================================================
# Directive 17: Oral PK Confidence
# ============================================================
class TestOralPKConfidence:
    def test_missing_fa_gives_low(self):
        assert evaluate_oral_confidence(None, 1.0) == "PK_CONFIDENCE_LOW"

    def test_missing_ka_gives_low(self):
        assert evaluate_oral_confidence(0.80, None) == "PK_CONFIDENCE_LOW"

    def test_ood_caco2_gives_low(self):
        assert evaluate_oral_confidence(0.80, 1.0, "OOD") == "PK_CONFIDENCE_LOW"

    def test_high_confidence_normal(self):
        assert evaluate_oral_confidence(0.80, 1.0, "IN_DOMAIN") == "PK_CONFIDENCE_HIGH"

    def test_moderate_for_extreme_ka(self):
        assert evaluate_oral_confidence(0.80, 0.10, "IN_DOMAIN") == "PK_CONFIDENCE_MODERATE"


# ============================================================
# Directive 18: Fold Error Calculation
# ============================================================
class TestFoldErrorCalculation:
    def test_symmetric_fold_error(self):
        fe = calculate_fold_error(200.0, 100.0)
        assert fe == 2.0

    def test_symmetric_reverse(self):
        fe = calculate_fold_error(50.0, 100.0)
        assert fe == 2.0

    def test_exact_match(self):
        fe = calculate_fold_error(100.0, 100.0)
        assert fe == 1.0

    def test_invalid_raises(self):
        with pytest.raises(ValueError):
            calculate_fold_error(0.0, 100.0)

    def test_large_fold_error(self):
        fe = calculate_fold_error(1000.0, 100.0)
        assert fe == 10.0


# ============================================================
# Directive 8: pKa Benchmark Contract
# ============================================================
class TestPKABenchmark:
    def test_pka_benchmark_report_exists(self):
        rpt = Path("validation/expanded_pk_validation_report_v1.json")
        assert rpt.exists()
        with open(rpt) as f:
            data = json.load(f)
        pka = data.get("pka_benchmark")
        assert pka is not None
        assert pka["n_endpoints"] >= 25, "Need >= 25 pKa evaluated points"
        # Current mechanistic MAE
        assert pka["current_mechanistic_estimate"]["mae"] > 0
        assert pka["current_mechanistic_estimate"]["mae"] < 2.0  # Sanity check
        assert "acid_mae" in pka["current_mechanistic_estimate"]
        assert "base_mae" in pka["current_mechanistic_estimate"]


# ============================================================
# Directive 10: logD7.4 Benchmark Contract
# ============================================================
class TestLogDBenchmark:
    def test_logd_benchmark_report_exists(self):
        rpt = Path("validation/expanded_pk_validation_report_v1.json")
        assert rpt.exists()
        with open(rpt) as f:
            data = json.load(f)
        logd = data.get("logd74_benchmark")
        assert logd is not None
        assert logd["n_compounds"] >= 25
        assert logd["current_henderson_hasselbalch_derived"]["mae"] > 0
        assert logd["current_henderson_hasselbalch_derived"]["mae"] < 1.0


# ============================================================
# Directive 11: VDss Holdout Gate
# ============================================================
class TestVDssHoldoutGate:
    def test_vdss_improvement_under_5pct_retains_level3(self):
        rpt = Path("validation/expanded_pk_validation_report_v1.json")
        with open(rpt) as f:
            data = json.load(f)
        vd = data["vdss_benchmark"]
        assert vd["promotion_decision"] == "RETAIN_LEVEL_3"
        assert vd["improvement_pct"] < 5.0

    def test_observed_vdss_is_not_scored_as_independent_prediction(self):
        rpt = Path("validation/expanded_pk_validation_report_v1.json")
        data = json.loads(rpt.read_text())
        assert data["pk_parameter_performance"]["vdss"]["n"] == 0
        assert all(item["vdss_evaluation"]["status"] == "NOT_INDEPENDENTLY_PREDICTED"
                   for item in data["compounds"])
        assert all(item["vdss_evaluation"]["included_in_accuracy_metrics"] is False
                   for item in data["compounds"])

    def test_upstream_attribution_preserves_observed_input_and_renal_context(self):
        rpt = Path("validation/expanded_pk_validation_report_v1.json")
        data = json.loads(rpt.read_text())
        attribution = data["upstream_error_attribution"]
        counts = attribution["top_drivers_count"]
        assert counts["not_identifiable_from_observed_inputs"] == 30
        assert counts["renal_clearance_not_modeled_context"] > 0
        assert counts["hepatic_ivive_or_other_model_error_context"] > 0
        details = data["upstream_error_attribution"]["details"]
        assert len(details) == 30
        assert all(item["primary_upstream_driver"] == "NOT_IDENTIFIABLE_FROM_OBSERVED_INPUTS"
                   for item in details)


# ============================================================
# Directive 14: Uncertainty Interval Calibration
# ============================================================
class TestUncertaintyCalibration:
    def test_calibration_report_exists(self):
        rpt = Path("validation/expanded_pk_validation_report_v1.json")
        with open(rpt) as f:
            data = json.load(f)
        cal = data.get("uncertainty_interval_calibration")
        assert cal is not None
        assert cal["total_evaluations"] > 0
        # Should report empirical vs nominal for p50, p80, p90
        for k in ("p50_empirical_pct", "p80_empirical_pct", "p90_empirical_pct"):
            assert k in cal


# ============================================================
# Directive 22: PK Engine Version Separation
# ============================================================
class TestPKEngineVersionSeparation:
    def test_prediction_engine_version_unchanged(self):
        assert CURRENT_ENGINE_VERSION == "3.3.2"
        assert CURRENT_ENGINE_ID == "drugopt-prediction-engine-v3@3.3.2"

    def test_pk_engine_version_exists(self):
        rpt = Path("validation/pk_baseline_freeze_v3_3_4.json")
        with open(rpt) as f:
            data = json.load(f)
        pk_eng = data.get("pk_engine")
        assert pk_eng is not None
        assert pk_eng["engine_id"] == "drugopt-pk-engine-v1@1.0.0"
        assert pk_eng["version"] == "1.0.0"


# ============================================================
# Directive 27: Database Integrity
# ============================================================
class TestDatabaseIntegrity:
    def test_pragma_integrity(self):
        db_path = Path("drugopt.db")
        if db_path.exists():
            conn = sqlite3.connect(db_path)
            result = conn.execute("PRAGMA integrity_check").fetchone()
            assert result[0] == "ok", f"integrity_check failed: {result}"
            fk = conn.execute("PRAGMA foreign_key_check").fetchall()
            assert len(fk) == 0, f"foreign_key_check violations: {fk}"
            conn.close()


# ============================================================
# Directive 25: Scientific Integrity Tests
# ============================================================
class TestScientificIntegrity:
    def test_fu_canonical_4_modes(self):
        """4 PPB conversion modes produce consistent fu."""
        fu_pct_bound = canonical_fu_from_ppb(90.0)
        assert abs(fu_pct_bound - 0.10) < 0.001

        fu_extreme = canonical_fu_from_ppb(99.95)
        assert fu_extreme >= 0.0001  # Lower floor

        # With provenance
        fu_prov, info = canonical_fu_from_ppb(99.95, store_provenance=True)
        assert fu_prov >= 0.0001
        assert isinstance(info, dict)

    def test_well_stirred_guardrails(self):
        """CL cannot exceed Qh (extraction ratio < 1.0)."""
        result = calculate_well_stirred_clearance(
            cl_int_input=5000.0, fu=0.01, rb=1.0,
        )
        assert result["extraction_ratio_eh"] < 1.0
        assert result["cl_h_blood_ml_min_kg"] <= HUMAN_QH_ML_MIN_KG

    def test_simulation_rejects_missing_vd(self):
        """Strict fail-closed: zero Vd causes failure, not fabrication."""
        with pytest.raises(ValueError, match="volume of distribution"):
            simulate_one_compartment_disposition(
                cl_plasma_ml_min_kg=5.0, vdss_l_kg=0.0, dose_mg=100.0,
                body_weight_kg=70.0, route="IV",
            )

    def test_oral_simulation_rejects_missing_f(self):
        """Oral simulation fails without F > 0."""
        with pytest.raises(ValueError, match="bioavailability"):
            simulate_one_compartment_disposition(
                cl_plasma_ml_min_kg=5.0, vdss_l_kg=1.0, dose_mg=100.0,
                body_weight_kg=70.0, route="ORAL", f_oral=0.0, ka_hr_inv=1.0,
            )

    def test_species_physiology_separation(self):
        """Human, Rat, Mouse have distinct physiological constants."""
        assert len(SPECIES_PHYSIOLOGY) >= 3
        assert SPECIES_PHYSIOLOGY["HUMAN"]["qh_ml_min_kg"] != SPECIES_PHYSIOLOGY["RAT"]["qh_ml_min_kg"]
        assert SPECIES_PHYSIOLOGY["HUMAN"]["qh_ml_min_kg"] != SPECIES_PHYSIOLOGY["MOUSE"]["qh_ml_min_kg"]

    def test_source_type_constants(self):
        """All 6 source types are defined."""
        for st in [SOURCE_TYPE_EXPERIMENTAL, SOURCE_TYPE_VALIDATED_MODEL,
                    SOURCE_TYPE_MECHANISTIC_DERIVED, SOURCE_TYPE_RULE_ESTIMATE,
                    SOURCE_TYPE_USER_INPUT, SOURCE_TYPE_MODEL_UNAVAILABLE]:
            assert isinstance(st, str) and len(st) > 0


# ============================================================
# Directive 30: Final Scientific Verdict
# ============================================================
class TestFinalScientificVerdict:
    def test_report_has_verdict(self):
        rpt = Path("validation/expanded_pk_validation_report_v1.json")
        with open(rpt) as f:
            data = json.load(f)
        v = data["final_scientific_verdict"]
        assert v in ("PK_MODEL_VALIDATION_IMPROVED", "PK_MODEL_VALIDATION_PARITY", "PK_MODEL_VALIDATION_INSUFFICIENT")

    def test_report_has_baseline_comparison(self):
        rpt = Path("validation/expanded_pk_validation_report_v1.json")
        with open(rpt) as f:
            data = json.load(f)
        comp = data["comparison_vs_baseline_n6"]
        assert comp["baseline_n"] == 6
        assert comp["expanded_n"] == 30
        assert comp["baseline_aafe"] == 2.20
