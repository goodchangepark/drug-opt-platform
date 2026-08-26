"""Inductive Conformal Prediction, Statistical Coverage Validation & Uncertainty Governance Engine (Stage 4C-3B).

Scientific Rules & Governance:
1. Conformal prediction provides distribution-free, finite-sample valid prediction intervals/sets.
2. Terminology: "90% Conformal Prediction Interval" / "90% Conformal Prediction Set" (NEVER "confidence interval" or "model certainty").
3. Decoupled Data Provenance & Calibration Quality:
   - Calibration Data Provenance: EXTERNAL, INTERNAL, TRAINING_OVERLAP_UNKNOWN, UNAVAILABLE.
   - Calibration Quality: VALIDATED, BORDERLINE, UNDERCOVERED, OVERCOVERED, INSUFFICIENT_N, INVALID, UNAVAILABLE.
   - An interval is NEVER called CALIBRATED merely because its evaluation dataset is external.
4. Statistical Coverage Acceptance:
   - Evaluated via exact binomial tests (Clopper-Pearson 95% confidence interval, binomial SE, deviation, p-value).
   - Arbitrary +-2% thresholds are prohibited.
5. Minimum Sample Size Policy:
   - Minimum calibration and evaluation N >= 30. Datasets with N < 30 (e.g. Caco-2 N=17) are strictly classified INSUFFICIENT_N.
6. Independent Split Conformal Calibration:
   - Quantiles computed strictly on calibration set ONLY. Evaluated strictly on independent evaluation set ONLY.
   - No tuning on evaluation set.
7. Interval Utility (Regression):
   - Evaluates median width, dynamic range relative width, ratio to MAE.
   - Flag UNINFORMATIVE_INTERVAL if width is excessively broad (e.g. Caco-2).
8. Prediction Set Efficiency (Classification):
   - Evaluates empirical coverage, singleton-set rate, ambiguous-set rate ({0, 1}), empty-set rate.
9. Applicability Domain (AD) Independence & Stratified Conditional Coverage:
   - Conformal intervals do NOT replace chemical space AD.
   - If AD = OUT_OF_DOMAIN, output includes explicit warning: "OUT OF DOMAIN — CONFORMAL COVERAGE MAY NOT GENERALIZE".
   - Stratified conditional coverage computed for IN_DOMAIN, BORDERLINE, OUT_OF_DOMAIN (minimum stratum N >= 15).
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.stats import binomtest


class DataProvenance:
    EXTERNAL = "EXTERNAL"
    INTERNAL = "INTERNAL"
    TRAINING_OVERLAP_UNKNOWN = "TRAINING_OVERLAP_UNKNOWN"
    UNAVAILABLE = "UNAVAILABLE"


class CalibrationQuality:
    VALIDATED = "VALIDATED"
    BORDERLINE = "BORDERLINE"
    UNDERCOVERED = "UNDERCOVERED"
    OVERCOVERED = "OVERCOVERED"
    INSUFFICIENT_N = "INSUFFICIENT_N"
    INVALID = "INVALID"
    UNAVAILABLE = "UNAVAILABLE"


class IntervalUtility:
    INFORMATIVE = "INFORMATIVE"
    UNINFORMATIVE_INTERVAL = "UNINFORMATIVE_INTERVAL"
    NOT_APPLICABLE = "NOT_APPLICABLE"


# Pre-computed empirical nonconformity quantiles and exact statistical validation results
# Generated from Stage 4C-3B Conformal Recalibration & Governance Audit (validation/stage4c3b_conformal_audit.json)
CONFORMAL_CALIBRATION_REGISTRY: dict[str, dict[str, Any]] = {
    "Permeability": {
        "data_provenance": DataProvenance.EXTERNAL,
        "calibration_quality": CalibrationQuality.INSUFFICIENT_N,
        "status": CalibrationQuality.INSUFFICIENT_N,
        "is_validated": False,
        "method": "Inductive Conformal Prediction (ICP)",
        "endpoint_type": "REGRESSION",
        "unit": "LogPapp",
        "dataset_name": "Admetica External Caco-2 Benchmark (34 compounds)",
        "calibration_n": 17,
        "evaluation_n": 17,
        "smiles_overlap_with_training": 0,
        "nominal_coverage": 0.90,
        "empirical_coverage": 0.7647,
        "expected_sampling_uncertainty_se": 0.0728,
        "confidence_interval_95": [0.5010, 0.9319],
        "deviation": -0.1353,
        "z_score": -1.859,
        "p_value": 0.082641,
        "quantiles": {
            "0.80": 10.299,
            "0.90": 11.147,
            "0.95": 11.290,
        },
        "empirical_coverage_levels": {
            "0.80": 0.529,
            "0.90": 0.765,
            "0.95": 0.824,
        },
        "interval_utility": {
            "mean_interval_width_90": 22.294,
            "median_interval_width_90": 22.294,
            "dynamic_range": 2.110,
            "relative_interval_width": 10.566,
            "ratio_to_mae": 2.14,
            "utility_status": IntervalUtility.UNINFORMATIVE_INTERVAL,
        },
        "conformal_governance_message": "Evaluation sample size N=17 < 30; cannot statistically validate coverage. Conformal interval width (22.3 LogPapp units) is uninformative.",
    },
    "HLM intrinsic clearance": {
        "data_provenance": DataProvenance.EXTERNAL,
        "calibration_quality": CalibrationQuality.UNDERCOVERED,
        "status": CalibrationQuality.UNDERCOVERED,
        "is_validated": False,
        "method": "Inductive Conformal Prediction (ICP)",
        "endpoint_type": "REGRESSION",
        "unit": "log10(mL/min/kg)",
        "dataset_name": "Biogen Public ADME Prospective Benchmark",
        "calibration_n": 250,
        "evaluation_n": 250,
        "smiles_overlap_with_training": 0,
        "nominal_coverage": 0.90,
        "empirical_coverage": 0.7960,
        "expected_sampling_uncertainty_se": 0.0190,
        "confidence_interval_95": [0.7407, 0.8442],
        "deviation": -0.1040,
        "z_score": -5.481,
        "p_value": 0.000001,
        "quantiles": {
            "0.80": 0.890,
            "0.90": 1.048,
            "0.95": 1.207,
        },
        "empirical_coverage_levels": {
            "0.80": 0.732,
            "0.90": 0.796,
            "0.95": 0.860,
        },
        "interval_utility": {
            "mean_interval_width_90": 2.096,
            "median_interval_width_90": 2.096,
            "dynamic_range": 2.664,
            "relative_interval_width": 0.787,
            "ratio_to_mae": 3.20,
            "utility_status": IntervalUtility.INFORMATIVE,
        },
        "conformal_governance_message": "Empirical coverage (79.6%) is significantly below nominal 90.0% (p<0.0001, z=-5.48). Interval undercovered on prospective benchmark.",
    },
    "hERG liability": {
        "data_provenance": DataProvenance.EXTERNAL,
        "calibration_quality": CalibrationQuality.UNDERCOVERED,
        "status": CalibrationQuality.UNDERCOVERED,
        "is_validated": False,
        "method": "Conformal Classification Prediction Sets",
        "endpoint_type": "CLASSIFICATION",
        "unit": "probability",
        "dataset_name": "ChEMBL37 Non-overlapping hERG IC50 Benchmark",
        "calibration_n": 250,
        "evaluation_n": 250,
        "smiles_overlap_with_training": 0,
        "nominal_coverage": 0.90,
        "empirical_coverage": 0.8320,
        "expected_sampling_uncertainty_se": 0.0190,
        "confidence_interval_95": [0.7798, 0.8762],
        "deviation": -0.0680,
        "z_score": -3.584,
        "p_value": 0.000961,
        "threshold_0.90": 0.999,
        "quantile_90": 0.999,
        "threshold_low": 0.0008,
        "threshold_high": 0.9992,
        "set_efficiency": {
            "singleton_rate": 0.468,
            "ambiguous_rate": 0.532,
            "empty_rate": 0.0,
            "efficiency_status": "HIGH_AMBIGUITY",
        },
        "conformal_governance_message": "Empirical coverage (83.2%) is significantly below nominal 90.0% (p=0.00096, z=-3.58).",
    },
    "CYP2C9 inhibitor": {
        "data_provenance": DataProvenance.INTERNAL,
        "calibration_quality": CalibrationQuality.UNDERCOVERED,
        "status": CalibrationQuality.UNDERCOVERED,
        "is_validated": False,
        "method": "Conformal Classification Prediction Sets",
        "endpoint_type": "CLASSIFICATION",
        "unit": "probability",
        "dataset_name": "ChEMBL30 CYP2C9 inhibitor Benchmark",
        "calibration_n": 232,
        "evaluation_n": 232,
        "smiles_overlap_with_training": 1,
        "nominal_coverage": 0.90,
        "empirical_coverage": 0.7974,
        "expected_sampling_uncertainty_se": 0.0197,
        "confidence_interval_95": [0.7399, 0.8472],
        "deviation": -0.1026,
        "z_score": -5.208,
        "p_value": 0.000003,
        "threshold_0.90": 0.932,
        "quantile_90": 0.932,
        "threshold_low": 0.0683,
        "threshold_high": 0.9317,
        "set_efficiency": {
            "singleton_rate": 0.487,
            "ambiguous_rate": 0.513,
            "empty_rate": 0.0,
            "efficiency_status": "HIGH_AMBIGUITY",
        },
        "conformal_governance_message": "Empirical coverage (79.7%) is significantly below nominal 90.0% (p<0.0001, z=-5.21). Internal calibration data with 1 compound training overlap.",
    },
    "CYP2D6 inhibitor": {
        "data_provenance": DataProvenance.INTERNAL,
        "calibration_quality": CalibrationQuality.VALIDATED,
        "status": CalibrationQuality.VALIDATED,
        "is_validated": True,
        "method": "Conformal Classification Prediction Sets",
        "endpoint_type": "CLASSIFICATION",
        "unit": "probability",
        "dataset_name": "ChEMBL30 CYP2D6 inhibitor Benchmark",
        "calibration_n": 250,
        "evaluation_n": 250,
        "smiles_overlap_with_training": 4,
        "nominal_coverage": 0.90,
        "empirical_coverage": 0.9040,
        "expected_sampling_uncertainty_se": 0.0190,
        "confidence_interval_95": [0.8605, 0.9375],
        "deviation": 0.0040,
        "z_score": 0.211,
        "p_value": 0.916179,
        "threshold_0.90": 0.977,
        "quantile_90": 0.977,
        "threshold_low": 0.0233,
        "threshold_high": 0.9767,
        "set_efficiency": {
            "singleton_rate": 0.404,
            "ambiguous_rate": 0.596,
            "empty_rate": 0.0,
            "efficiency_status": "HIGH_AMBIGUITY",
        },
        "conformal_governance_message": "Empirical coverage (90.4%) statistically validated within 95% CI [86.1%, 93.8%]. Internal calibration data with 4 compounds training overlap.",
    },
    "CYP3A4 inhibitor": {
        "data_provenance": DataProvenance.EXTERNAL,
        "calibration_quality": CalibrationQuality.VALIDATED,
        "status": CalibrationQuality.VALIDATED,
        "is_validated": True,
        "method": "Conformal Classification Prediction Sets",
        "endpoint_type": "CLASSIFICATION",
        "unit": "probability",
        "dataset_name": "ChEMBL30 CYP3A4 inhibitor Benchmark",
        "calibration_n": 250,
        "evaluation_n": 250,
        "smiles_overlap_with_training": 0,
        "nominal_coverage": 0.90,
        "empirical_coverage": 0.8800,
        "expected_sampling_uncertainty_se": 0.0190,
        "confidence_interval_95": [0.8331, 0.9176],
        "deviation": -0.0200,
        "z_score": -1.054,
        "p_value": 0.291112,
        "threshold_0.90": 0.958,
        "quantile_90": 0.958,
        "threshold_low": 0.0415,
        "threshold_high": 0.9585,
        "set_efficiency": {
            "singleton_rate": 0.356,
            "ambiguous_rate": 0.644,
            "empty_rate": 0.0,
            "efficiency_status": "HIGH_AMBIGUITY",
        },
        "conformal_governance_message": "Empirical coverage (88.0%) statistically validated within 95% CI [83.3%, 91.8%]. External validation set with 0 training overlap.",
    },
}


def validate_conformal_coverage(
    k_covered: int,
    n_eval: int,
    nominal_coverage: float = 0.90,
    min_eval_n: int = 30,
    alpha_sig: float = 0.05,
) -> dict[str, Any]:
    """Perform exact binomial coverage validation and calculate statistical sampling uncertainty.
    
    Statistical Rules:
    1. For nominal coverage p0 = 1 - alpha, calculates expected sampling standard error SE = sqrt(p0*(1-p0)/n).
    2. Calculates exact Clopper-Pearson 95% binomial confidence interval [ci_lower, ci_upper].
    3. If n_eval < min_eval_n (e.g. n < 30): quality status is INSUFFICIENT_N.
    4. If empirical coverage < p0 and (p0 > ci_upper or one-sided binomial test p < alpha_sig): UNDERCOVERED.
    5. If empirical coverage > p0 and p0 < ci_lower and empirical coverage > 0.96: OVERCOVERED.
    6. If p0 is within [ci_lower, ci_upper]: VALIDATED.
    """
    if n_eval <= 0:
        return {
            "quality_status": CalibrationQuality.UNAVAILABLE,
            "empirical_coverage": None,
            "nominal_coverage": nominal_coverage,
            "evaluation_n": 0,
            "covered_n": 0,
            "sampling_uncertainty_se": None,
            "confidence_interval_95": None,
            "deviation": None,
            "z_score": None,
            "p_value_two_sided": None,
            "p_value_undercoverage": None,
            "is_validated": False,
            "message": "No evaluation data available.",
        }

    empirical = k_covered / n_eval
    se = math.sqrt(nominal_coverage * (1.0 - nominal_coverage) / n_eval)
    deviation = empirical - nominal_coverage
    z_score = deviation / se if se > 0 else 0.0

    btest = binomtest(k=k_covered, n=n_eval, p=nominal_coverage, alternative="two-sided")
    btest_less = binomtest(k=k_covered, n=n_eval, p=nominal_coverage, alternative="less")
    ci = btest.proportion_ci(confidence_level=1.0 - alpha_sig, method="exact")
    ci_low = round(float(ci.low), 4)
    ci_high = round(float(ci.high), 4)
    p_two_sided = round(float(btest.pvalue), 6)
    p_less = round(float(btest_less.pvalue), 6)

    if n_eval < min_eval_n:
        quality_status = CalibrationQuality.INSUFFICIENT_N
        is_validated = False
        message = f"Evaluation sample size N={n_eval} < {min_eval_n}; cannot statistically validate coverage."
    elif empirical < nominal_coverage and (nominal_coverage > ci_high or p_less < alpha_sig):
        quality_status = CalibrationQuality.UNDERCOVERED
        is_validated = False
        message = f"Empirical coverage ({empirical:.1%}) is significantly below nominal {nominal_coverage:.1%} (p={p_less:.4f}, z={z_score:.2f})."
    elif empirical > nominal_coverage and nominal_coverage < ci_low and empirical > 0.96:
        quality_status = CalibrationQuality.OVERCOVERED
        is_validated = False
        message = f"Empirical coverage ({empirical:.1%}) is significantly overcovered."
    elif ci_low <= nominal_coverage <= ci_high:
        quality_status = CalibrationQuality.VALIDATED
        is_validated = True
        message = f"Empirical coverage ({empirical:.1%}) statistically validated within 95% CI [{ci_low:.1%}, {ci_high:.1%}]."
    else:
        quality_status = CalibrationQuality.BORDERLINE
        is_validated = False
        message = f"Empirical coverage ({empirical:.1%}) is borderline."

    return {
        "quality_status": quality_status,
        "empirical_coverage": round(empirical, 4),
        "nominal_coverage": nominal_coverage,
        "evaluation_n": n_eval,
        "covered_n": k_covered,
        "sampling_uncertainty_se": round(se, 4),
        "confidence_interval_95": [ci_low, ci_high],
        "deviation": round(deviation, 4),
        "z_score": round(z_score, 3),
        "p_value_two_sided": p_two_sided,
        "p_value_undercoverage": p_less,
        "is_validated": is_validated,
        "message": message,
    }


def evaluate_regression_interval_utility(
    quantiles: dict[str, float],
    eval_errors: list[float] | None = None,
    dynamic_range: float | None = None,
    max_relative_width_ratio: float = 1.2,
) -> dict[str, Any]:
    """Evaluate interval utility and informativeness for regression conformal predictions."""
    q90 = quantiles.get("0.90", quantiles.get(0.90, 0.0))
    interval_width_90 = round(2.0 * q90, 3)

    mae = round(float(np.mean(eval_errors)), 3) if eval_errors else None
    ratio_to_mae = round(interval_width_90 / mae, 2) if mae and mae > 0 else None

    relative_width = round(interval_width_90 / dynamic_range, 3) if dynamic_range and dynamic_range > 0 else None
    is_uninformative = bool(relative_width is not None and relative_width > max_relative_width_ratio)

    return {
        "median_interval_width_90": interval_width_90,
        "mean_interval_width_90": interval_width_90,
        "endpoint_dynamic_range": dynamic_range,
        "relative_interval_width": relative_width,
        "ratio_to_mae": ratio_to_mae,
        "is_uninformative": is_uninformative,
        "utility_status": IntervalUtility.UNINFORMATIVE_INTERVAL if is_uninformative else IntervalUtility.INFORMATIVE,
    }


def evaluate_classification_set_efficiency(
    eval_y_true: list[int],
    eval_y_prob: list[float],
    threshold_low: float,
    threshold_high: float,
) -> dict[str, Any]:
    """Evaluate conformal prediction set validity, singleton rate, and ambiguity for classification."""
    n_eval = len(eval_y_true)
    if n_eval == 0:
        return {"evaluation_n": 0, "status": "NO_DATA"}

    covered_count = 0
    singleton_count = 0
    ambiguous_count = 0
    empty_count = 0

    for y, p in zip(eval_y_true, eval_y_prob):
        pred_set = set()
        if p >= threshold_low:
            pred_set.add(1)
        if p <= threshold_high:
            pred_set.add(0)

        if y in pred_set:
            covered_count += 1
        if len(pred_set) == 1:
            singleton_count += 1
        elif len(pred_set) > 1:
            ambiguous_count += 1
        else:
            empty_count += 1

    emp_coverage = covered_count / n_eval
    singleton_rate = singleton_count / n_eval
    ambiguous_rate = ambiguous_count / n_eval
    empty_rate = empty_count / n_eval

    return {
        "evaluation_n": n_eval,
        "covered_n": covered_count,
        "empirical_coverage": round(emp_coverage, 4),
        "singleton_rate": round(singleton_rate, 3),
        "ambiguous_rate": round(ambiguous_rate, 3),
        "empty_rate": round(empty_rate, 3),
        "efficiency_status": "HIGH_AMBIGUITY" if ambiguous_rate > 0.50 else "EFFICIENT",
        "efficiency_summary": f"{singleton_rate:.1%} singletons, {ambiguous_rate:.1%} ambiguous, {empty_rate:.1%} empty",
    }


def evaluate_ad_stratified_coverage(
    y_true: list[float] | list[int],
    y_pred_or_prob: list[float],
    ad_classifications: list[str],
    quantile_or_threshold: float,
    endpoint_type: str = "REGRESSION",
    nominal_coverage: float = 0.90,
    min_stratum_n: int = 15,
) -> dict[str, Any]:
    """Inspect conditional coverage stratified by Applicability Domain (IN_DOMAIN, BORDERLINE, OUT_OF_DOMAIN)."""
    strata: dict[str, list[int]] = {"IN_DOMAIN": [], "BORDERLINE": [], "OUT_OF_DOMAIN": []}
    for idx, ad in enumerate(ad_classifications):
        if ad in strata:
            strata[ad].append(idx)
        else:
            strata["OUT_OF_DOMAIN"].append(idx)

    stratified_results = {}
    for domain, indices in strata.items():
        n_stratum = len(indices)
        if n_stratum == 0:
            stratified_results[domain] = {
                "n": 0,
                "empirical_coverage": None,
                "quality": "NO_DATA",
                "message": "No compounds in this domain stratum.",
            }
            continue

        if endpoint_type == "REGRESSION":
            hits = sum(1 for i in indices if abs(y_true[i] - y_pred_or_prob[i]) <= quantile_or_threshold)
        else:
            t_low = 1.0 - quantile_or_threshold
            t_high = quantile_or_threshold
            hits = 0
            for i in indices:
                p = y_pred_or_prob[i]
                y = y_true[i]
                pset = set()
                if p >= t_low:
                    pset.add(1)
                if p <= t_high:
                    pset.add(0)
                if y in pset:
                    hits += 1

        emp_cov = round(hits / n_stratum, 4)
        if n_stratum < min_stratum_n:
            stratified_results[domain] = {
                "n": n_stratum,
                "covered_n": hits,
                "empirical_coverage": emp_cov,
                "quality": CalibrationQuality.INSUFFICIENT_N,
                "message": f"Stratum N={n_stratum} < {min_stratum_n}; sample size too small for conditional coverage claim.",
            }
        else:
            btest = validate_conformal_coverage(hits, n_stratum, nominal_coverage=nominal_coverage, min_eval_n=min_stratum_n)
            stratified_results[domain] = {
                "n": n_stratum,
                "covered_n": hits,
                "empirical_coverage": emp_cov,
                "quality": btest["quality_status"],
                "sampling_uncertainty_se": btest["sampling_uncertainty_se"],
                "confidence_interval_95": btest["confidence_interval_95"],
                "message": btest["message"],
            }

    return stratified_results


def compute_calibrated_uncertainty(
    endpoint: str,
    predicted_value: float | None,
    applicability_domain: dict[str, Any] | None = None,
    nominal_level: str = "0.90",
) -> dict[str, Any]:
    """Compute calibrated conformal prediction interval or classification prediction set.

    Scientific & Governance Rules:
    1. Independent Provenance & Quality:
       Returns data_provenance (EXTERNAL/INTERNAL/TRAINING_OVERLAP_UNKNOWN/UNAVAILABLE) and
       calibration_quality (VALIDATED/UNDERCOVERED/INSUFFICIENT_N/UNAVAILABLE) separately.
    2. Preserves similarity AD status: If AD = OUT_OF_DOMAIN, adds explicit warning:
       'OUT OF DOMAIN — CONFORMAL COVERAGE MAY NOT GENERALIZE: Structure is out of chemical-space applicability domain.'
    3. Small-N Policy: Small N (Caco-2 N=17) classified INSUFFICIENT_N with warning.
    4. Regression Interval Utility: Uninformative intervals flagged with UNINFORMATIVE_INTERVAL.
    """
    if endpoint not in CONFORMAL_CALIBRATION_REGISTRY:
        prov = (
            DataProvenance.TRAINING_OVERLAP_UNKNOWN
            if endpoint in {
                "Solubility", "Plasma protein binding", "CYP1A2 inhibitor", "CYP2C19 inhibitor",
                "CYP2C9 substrate", "CYP2D6 substrate", "CYP3A4 substrate", "P-gp inhibitor",
                "Ames mutagenicity", "DILI clinical liability"
            }
            else DataProvenance.UNAVAILABLE
        )
        return {
            "data_provenance": prov,
            "calibration_quality": CalibrationQuality.UNAVAILABLE,
            "status": CalibrationQuality.UNAVAILABLE,
            "is_validated": False,
            "reason": f"Independent external calibration set not provided on disk for endpoint '{endpoint}'.",
            "interval": None,
            "prediction_set": None,
            "warnings": ["CONFORMAL_UNAVAILABLE: No qualified independent calibration set available."],
        }

    cal_spec = CONFORMAL_CALIBRATION_REGISTRY[endpoint]
    ad_status = (applicability_domain or {}).get("classification", "IN_DOMAIN")
    warnings = []

    if ad_status == "OUT_OF_DOMAIN":
        warnings.append("OUT OF DOMAIN — CONFORMAL COVERAGE MAY NOT GENERALIZE: Structure is out of chemical-space applicability domain.")

    if cal_spec["calibration_quality"] == CalibrationQuality.UNDERCOVERED:
        warnings.append(f"CALIBRATION_UNDERCOVERED: Observed empirical coverage ({cal_spec['empirical_coverage']:.1%}) fails nominal {cal_spec['nominal_coverage']:.1%} target.")
    elif cal_spec["calibration_quality"] == CalibrationQuality.INSUFFICIENT_N:
        warnings.append(f"CONFORMAL_INSUFFICIENT_N: Evaluation sample size (N={cal_spec['evaluation_n']}) is below statistical reliability threshold.")

    if cal_spec.get("interval_utility", {}).get("utility_status") == IntervalUtility.UNINFORMATIVE_INTERVAL:
        warnings.append("UNINFORMATIVE_INTERVAL: Conformal prediction interval width exceeds endpoint dynamic range.")

    if cal_spec["endpoint_type"] == "REGRESSION":
        if predicted_value is None or math.isnan(predicted_value):
            return {
                "data_provenance": cal_spec["data_provenance"],
                "calibration_quality": CalibrationQuality.INVALID,
                "status": CalibrationQuality.INVALID,
                "is_validated": False,
                "reason": "Predicted value is null.",
                "interval": None,
            }

        q = cal_spec["quantiles"].get(nominal_level, cal_spec["quantiles"]["0.90"])
        lower = round(predicted_value - q, 3)
        upper = round(predicted_value + q, 3)

        # Enforce physical bounds where appropriate
        if cal_spec["unit"] == "% bound":
            lower = max(0.0, lower)
            upper = min(100.0, upper)
        elif cal_spec["unit"] == "mL/min/kg":
            lower = max(0.0, lower)

        pct_label = int(float(nominal_level) * 100)
        display_label = f"{pct_label}% Conformal Prediction Interval"

        return {
            "data_provenance": cal_spec["data_provenance"],
            "calibration_quality": cal_spec["calibration_quality"],
            "status": cal_spec["calibration_quality"],
            "is_validated": cal_spec["is_validated"],
            "display_label": display_label,
            "nominal_coverage": float(nominal_level),
            "empirical_coverage": cal_spec["empirical_coverage_levels"].get(nominal_level, cal_spec["empirical_coverage"]),
            "expected_sampling_uncertainty_se": cal_spec["expected_sampling_uncertainty_se"],
            "confidence_interval_95": cal_spec["confidence_interval_95"],
            "deviation": cal_spec["deviation"],
            "lower_bound": lower,
            "upper_bound": upper,
            "interval_width": round(upper - lower, 3),
            "unit": cal_spec["unit"],
            "calibration_n": cal_spec["calibration_n"],
            "evaluation_n": cal_spec["evaluation_n"],
            "method": cal_spec["method"],
            "interval_utility": cal_spec.get("interval_utility", {}),
            "warnings": warnings,
        }
    else:  # CLASSIFICATION
        prob = predicted_value if (predicted_value is not None and not math.isnan(predicted_value)) else 0.5
        threshold = cal_spec.get("threshold_0.90", 0.95)

        pred_set = []
        if prob >= (1.0 - threshold):
            pred_set.append("POSITIVE")
        if (1.0 - prob) >= (1.0 - threshold):
            pred_set.append("NEGATIVE")

        is_uncertain_set = len(pred_set) > 1
        if is_uncertain_set:
            warnings.append("HIGH_CONFORMAL_UNCERTAINTY: Conformal prediction set contains both classes {POSITIVE, NEGATIVE}.")

        return {
            "data_provenance": cal_spec["data_provenance"],
            "calibration_quality": cal_spec["calibration_quality"],
            "status": cal_spec["calibration_quality"],
            "is_validated": cal_spec["is_validated"],
            "display_label": "90% Conformal Prediction Set",
            "nominal_coverage": cal_spec["nominal_coverage"],
            "empirical_coverage": cal_spec["empirical_coverage"],
            "expected_sampling_uncertainty_se": cal_spec["expected_sampling_uncertainty_se"],
            "confidence_interval_95": cal_spec["confidence_interval_95"],
            "deviation": cal_spec["deviation"],
            "prediction_set": pred_set,
            "is_uncertain_set": is_uncertain_set,
            "calibration_n": cal_spec["calibration_n"],
            "evaluation_n": cal_spec["evaluation_n"],
            "method": cal_spec["method"],
            "set_efficiency": cal_spec.get("set_efficiency", {}),
            "warnings": warnings,
        }


def evaluate_conformal_calibration_coverage(
    y_true: list[float],
    y_pred: list[float],
    quantile: float,
    endpoint_type: str = "REGRESSION",
) -> dict[str, Any]:
    """Calculate empirical coverage and statistical interval validation on evaluation set."""
    if not y_true or not y_pred or len(y_true) != len(y_pred):
        return {"status": "NO_DATA", "n_eval": 0}

    n_eval = len(y_true)

    if endpoint_type == "REGRESSION":
        errors = [abs(t - p) for t, p in zip(y_true, y_pred)]
        hits = sum(1 for e in errors if e <= quantile)
        cov_test = validate_conformal_coverage(hits, n_eval, nominal_coverage=0.90)
        mean_width = round(2 * quantile, 3)

        return {
            "status": "EVALUATED",
            "n_eval": n_eval,
            "quantile": quantile,
            "empirical_coverage": cov_test["empirical_coverage"],
            "sampling_uncertainty_se": cov_test["sampling_uncertainty_se"],
            "confidence_interval_95": cov_test["confidence_interval_95"],
            "quality_status": cov_test["quality_status"],
            "mean_interval_width": mean_width,
            "median_interval_width": mean_width,
            "coverage_error": round(abs(cov_test["empirical_coverage"] - 0.90), 4),
            "is_validated": cov_test["is_validated"],
        }
    else:
        return {"status": "NOT_IMPLEMENTED_FOR_CLASSIFICATION", "n_eval": n_eval}
