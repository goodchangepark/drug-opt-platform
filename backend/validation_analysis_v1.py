"""Engine v1 internal validation analysis engine.

Computes all primary and secondary metrics for paired prediction/experiment
observations in the validation campaign.

Scientific rules:
  - NO model fitting, retraining, recalibration
  - NO threshold modification
  - Metrics are descriptive only
  - Bootstrap uses fixed seed=42 at compound level
  - Small N: report data_insufficiency instead of fake precision
  - Non-positive log values: NON_POSITIVE_EXCLUDED_FROM_LOG_METRIC
  - Censored values: excluded from primary quantitative metrics
  - HISTORICAL_VISIBLE: kept separate from prospective/blinded-retro
  - Shadow: secondary evidence only, never changes production value
"""

from __future__ import annotations

import json
import math
import statistics
import warnings
from typing import Any, Dict, List, Optional, Tuple

# Minimum N constants
MIN_N_REGRESSION = 5
MIN_N_CLASSIFICATION = 5
MIN_N_AUROC = 5
MIN_N_BOOTSTRAP = 10
MIN_N_SERIES = 3

BOOTSTRAP_SEED = 42
BOOTSTRAP_N = 1000

# Log-scale endpoints (prediction values are already in log10)
LOG10_ENDPOINTS = {
    "solubility_aqueous_logs",
    "permeability_caco2_logpapp",
    "hlm_intrinsic_clearance_scaled_log10",
    "rlm_intrinsic_clearance_scaled_log10",
    "mlm_intrinsic_clearance_scaled_log10",
}

# Classification (probability) endpoints
PROB_ENDPOINTS = {
    "safety_herg_blocker_prob",
    "cyp3a4_inhibitor_prob",
    "cyp1a2_inhibitor_prob",
    "cyp2c19_inhibitor_prob",
    "cyp2c9_inhibitor_prob",
    "cyp2c9_substrate_prob",
    "cyp2d6_inhibitor_prob",
    "cyp2d6_substrate_prob",
    "cyp3a4_substrate_prob",
    "safety_ames_mutagenicity_prob",
    "safety_dili_clinical_prob",
    "transporter_pgp_inhibitor_prob",
}

# Endpoints where fold accuracy makes sense (log10 scale)
FOLD_ACCURACY_ENDPOINTS = {
    "solubility_aqueous_logs",
    "permeability_caco2_logpapp",
    "hlm_intrinsic_clearance_scaled_log10",
    "rlm_intrinsic_clearance_scaled_log10",
    "mlm_intrinsic_clearance_scaled_log10",
}

LOG10_2 = math.log10(2)
LOG10_3 = math.log10(3)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


class PairedObservation:
    """One paired prediction/experiment record for analysis."""

    def __init__(
        self,
        compound_id: str,
        endpoint_id: str,
        prediction_value: Optional[float],
        experimental_value: Optional[float],  # Already normalized to prediction units
        experimental_raw_value: Optional[float],  # Original raw value
        experimental_unit: str,
        qualifier: str,
        censor_flag: bool,
        applicability_domain: str,
        reliability: str,
        prospective_evidence_class: str,
        enters_primary_metrics: bool,
        scaffold_hash: str = "",
        series_label: str = "",
        project_label: str = "",
    ):
        self.compound_id = compound_id
        self.endpoint_id = endpoint_id
        self.prediction_value = prediction_value
        self.experimental_value = experimental_value
        self.experimental_raw_value = experimental_raw_value
        self.experimental_unit = experimental_unit
        self.qualifier = qualifier
        self.censor_flag = censor_flag
        self.applicability_domain = applicability_domain
        self.reliability = reliability
        self.prospective_evidence_class = prospective_evidence_class
        self.enters_primary_metrics = enters_primary_metrics
        self.scaffold_hash = scaffold_hash
        self.series_label = series_label
        self.project_label = project_label

    @property
    def signed_error(self) -> Optional[float]:
        if self.prediction_value is None or self.experimental_value is None:
            return None
        if self.censor_flag:
            return None
        return self.prediction_value - self.experimental_value

    @property
    def absolute_error(self) -> Optional[float]:
        se = self.signed_error
        return abs(se) if se is not None else None


# ---------------------------------------------------------------------------
# Metric computation helpers
# ---------------------------------------------------------------------------


def _percentile(values: List[float], p: float) -> float:
    """Simple percentile (linear interpolation)."""
    if not values:
        return float("nan")
    sv = sorted(values)
    idx = (len(sv) - 1) * p / 100.0
    lo = int(idx)
    hi = lo + 1
    if hi >= len(sv):
        return sv[lo]
    frac = idx - lo
    return sv[lo] * (1 - frac) + sv[hi] * frac


def _spearman_rho(xs: List[float], ys: List[float]) -> Optional[float]:
    """Spearman rank correlation. Returns None if N < 2."""
    if len(xs) < 2 or len(xs) != len(ys):
        return None
    try:
        import scipy.stats as stats
        rho, _ = stats.spearmanr(xs, ys)
        return float(rho)
    except ImportError:
        # Manual Spearman
        def rank(v):
            sv = sorted(enumerate(v), key=lambda x: x[1])
            r = [0.0] * len(v)
            i = 0
            while i < len(sv):
                j = i
                while j < len(sv) - 1 and sv[j + 1][1] == sv[j][1]:
                    j += 1
                avg_rank = (i + j) / 2.0 + 1
                for k in range(i, j + 1):
                    r[sv[k][0]] = avg_rank
                i = j + 1
            return r

        rx, ry = rank(xs), rank(ys)
        n = len(rx)
        mx = sum(rx) / n
        my = sum(ry) / n
        num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
        sx = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(n)))
        sy = math.sqrt(sum((ry[i] - my) ** 2 for i in range(n)))
        if sx == 0 or sy == 0:
            return None
        return num / (sx * sy)


def compute_regression_metrics(
    obs: List[PairedObservation],
    endpoint_id: str,
    filter_primary: bool = True,
) -> Dict[str, Any]:
    """Compute regression metrics for a set of observations."""
    if filter_primary:
        pairs = [
            o
            for o in obs
            if o.enters_primary_metrics
            and o.prediction_value is not None
            and o.experimental_value is not None
            and o.qualifier in ("=", "~", "")
            and not o.censor_flag
        ]
    else:
        pairs = [
            o
            for o in obs
            if o.prediction_value is not None and o.experimental_value is not None
        ]

    # Log-scale: check for non-positive experimental values
    excluded_nonpositive = 0
    if endpoint_id in LOG10_ENDPOINTS:
        valid = []
        for o in pairs:
            raw = o.experimental_raw_value
            if raw is not None and raw <= 0:
                excluded_nonpositive += 1
                continue
            valid.append(o)
        pairs = valid

    n = len(pairs)
    if n < MIN_N_REGRESSION:
        return {
            "n": n,
            "data_insufficient": True,
            "reason": f"N={n} < minimum {MIN_N_REGRESSION} for regression metrics",
            "excluded_nonpositive": excluded_nonpositive,
        }

    preds = [o.prediction_value for o in pairs]
    exps = [o.experimental_value for o in pairs]
    errors = [p - e for p, e in zip(preds, exps)]
    abs_errors = [abs(e) for e in errors]

    mae = statistics.mean(abs_errors)
    rmse = math.sqrt(statistics.mean([e ** 2 for e in errors]))
    bias = statistics.mean(errors)
    median_ae = statistics.median(abs_errors)
    p75 = _percentile(abs_errors, 75)
    p90 = _percentile(abs_errors, 90)
    p95 = _percentile(abs_errors, 95)
    spearman = _spearman_rho(preds, exps)

    result = {
        "n": n,
        "data_insufficient": False,
        "MAE": round(mae, 4),
        "RMSE": round(rmse, 4),
        "Bias": round(bias, 4),
        "Median_AE": round(median_ae, 4),
        "Spearman": round(spearman, 4) if spearman is not None else None,
        "P75_AE": round(p75, 4),
        "P90_AE": round(p90, 4),
        "P95_AE": round(p95, 4),
        "excluded_nonpositive": excluded_nonpositive,
    }

    # Fold accuracy for log-scale endpoints
    if endpoint_id in FOLD_ACCURACY_ENDPOINTS:
        w2 = sum(1 for ae in abs_errors if ae <= LOG10_2) / n
        w3 = sum(1 for ae in abs_errors if ae <= LOG10_3) / n
        result["Within_2fold"] = round(w2, 4)
        result["Within_3fold"] = round(w3, 4)

    # Pearson / R2 secondary
    try:
        import scipy.stats as stats
        r, _ = stats.pearsonr(preds, exps)
        result["Pearson_secondary"] = round(float(r), 4)
        result["R2_secondary"] = round(float(r) ** 2, 4)
    except Exception:
        pass

    return result


def compute_classification_metrics(
    obs: List[PairedObservation],
    endpoint_id: str,
    threshold: float = 0.5,
    filter_primary: bool = True,
) -> Dict[str, Any]:
    """Compute classification metrics for binary endpoints."""
    if filter_primary:
        pairs = [
            o
            for o in obs
            if o.enters_primary_metrics
            and o.prediction_value is not None
            and o.experimental_value is not None
        ]
    else:
        pairs = [
            o
            for o in obs
            if o.prediction_value is not None and o.experimental_value is not None
        ]

    n = len(pairs)
    if n < MIN_N_CLASSIFICATION:
        return {
            "n": n,
            "data_insufficient": True,
            "reason": f"N={n} < minimum {MIN_N_CLASSIFICATION}",
        }

    preds = [o.prediction_value for o in pairs]
    # Experimental values: 1 = positive/blocker, 0 = negative
    actuals = [int(round(o.experimental_value)) for o in pairs]

    prevalence = sum(actuals) / n
    pred_binary = [1 if p >= threshold else 0 for p in preds]

    tp = sum(1 for p, a in zip(pred_binary, actuals) if p == 1 and a == 1)
    tn = sum(1 for p, a in zip(pred_binary, actuals) if p == 0 and a == 0)
    fp = sum(1 for p, a in zip(pred_binary, actuals) if p == 1 and a == 0)
    fn = sum(1 for p, a in zip(pred_binary, actuals) if p == 0 and a == 1)

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else None
    specificity = tn / (tn + fp) if (tn + fp) > 0 else None
    balanced_acc = (
        (sensitivity + specificity) / 2
        if sensitivity is not None and specificity is not None
        else None
    )

    # MCC
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn) - (fp * fn)) / denom if denom > 0 else None

    # Brier score
    brier = statistics.mean([(p - a) ** 2 for p, a in zip(preds, actuals)])

    # Log loss
    eps = 1e-15
    logloss = -statistics.mean(
        [
            a * math.log(max(p, eps)) + (1 - a) * math.log(max(1 - p, eps))
            for p, a in zip(preds, actuals)
        ]
    )

    result = {
        "n": n,
        "data_insufficient": False,
        "threshold_used": threshold,
        "Prevalence": round(prevalence, 4),
        "TP": tp,
        "TN": tn,
        "FP": fp,
        "FN": fn,
        "Sensitivity": round(sensitivity, 4) if sensitivity is not None else None,
        "Specificity": round(specificity, 4) if specificity is not None else None,
        "Balanced_Accuracy": round(balanced_acc, 4) if balanced_acc is not None else None,
        "MCC": round(mcc, 4) if mcc is not None else None,
        "Brier": round(brier, 4),
        "LogLoss": round(logloss, 4),
    }

    # AUROC / AUPRC (need scipy or manual trapezoidal)
    if n >= MIN_N_AUROC and len(set(actuals)) > 1:
        try:
            from sklearn.metrics import roc_auc_score, average_precision_score
            result["AUROC"] = round(float(roc_auc_score(actuals, preds)), 4)
            result["AUPRC"] = round(float(average_precision_score(actuals, preds)), 4)
        except ImportError:
            result["AUROC"] = None
            result["AUPRC"] = None
            result["AUROC_note"] = "sklearn not available"

    # ECE (Expected Calibration Error, 10 bins)
    n_bins = 10
    bin_size = 1.0 / n_bins
    ece_sum = 0.0
    for i in range(n_bins):
        lo, hi = i * bin_size, (i + 1) * bin_size
        in_bin = [(p, a) for p, a in zip(preds, actuals) if lo <= p < hi]
        if in_bin:
            frac_pos = sum(a for _, a in in_bin) / len(in_bin)
            avg_pred = sum(p for p, _ in in_bin) / len(in_bin)
            ece_sum += len(in_bin) * abs(frac_pos - avg_pred) / n
    result["ECE"] = round(ece_sum, 4)

    return result


def analyze_ad_stratification(
    obs: List[PairedObservation],
    endpoint_id: str,
) -> Dict[str, Any]:
    """Analyze whether AD status predicts actual prediction error."""
    ad_groups: Dict[str, List[PairedObservation]] = {}
    for o in obs:
        if o.prediction_value is not None and o.experimental_value is not None:
            ad_groups.setdefault(o.applicability_domain, []).append(o)

    result = {
        "question": "Does error actually worsen as AD becomes weaker?",
        "groups": {},
    }

    for ad_class in ["IN_DOMAIN", "BORDERLINE", "OUT_OF_DOMAIN", "UNKNOWN"]:
        group = ad_groups.get(ad_class, [])
        if not group:
            result["groups"][ad_class] = {"n": 0, "note": "no observations"}
            continue

        if endpoint_id in PROB_ENDPOINTS:
            met = compute_classification_metrics(group, endpoint_id, filter_primary=False)
        else:
            met = compute_regression_metrics(group, endpoint_id, filter_primary=False)

        result["groups"][ad_class] = {
            "n": len(group),
            "metrics": met,
        }

    # Qualitative verdict
    in_d = ad_groups.get("IN_DOMAIN", [])
    out_d = ad_groups.get("OUT_OF_DOMAIN", [])
    if (
        in_d
        and out_d
        and endpoint_id not in PROB_ENDPOINTS
    ):
        errors_in = [abs(o.prediction_value - o.experimental_value) for o in in_d
                     if o.prediction_value is not None and o.experimental_value is not None]
        errors_out = [abs(o.prediction_value - o.experimental_value) for o in out_d
                      if o.prediction_value is not None and o.experimental_value is not None]
        if errors_in and errors_out:
            mae_in = statistics.mean(errors_in)
            mae_out = statistics.mean(errors_out)
            result["ad_monotone_verdict"] = (
                "AD_MONOTONE_AS_EXPECTED" if mae_out >= mae_in else "AD_NOT_MONOTONE"
            )
            result["mae_in_domain"] = round(mae_in, 4)
            result["mae_out_of_domain"] = round(mae_out, 4)
    elif len(ad_groups) < 2:
        result["ad_monotone_verdict"] = "INSUFFICIENT_AD_DIVERSITY"

    return result


def analyze_reliability_stratification(
    obs: List[PairedObservation],
    endpoint_id: str,
) -> Dict[str, Any]:
    """Analyze whether reliability classes correlate with observed accuracy."""
    rel_groups: Dict[str, List[PairedObservation]] = {}
    for o in obs:
        if o.prediction_value is not None and o.experimental_value is not None:
            rel_groups.setdefault(o.reliability, []).append(o)

    return {
        "question": "Does higher-reliability predict better accuracy?",
        "groups": {
            rel: {
                "n": len(grp),
                "metrics": (
                    compute_regression_metrics(grp, endpoint_id, filter_primary=False)
                    if endpoint_id not in PROB_ENDPOINTS
                    else compute_classification_metrics(grp, endpoint_id, filter_primary=False)
                ),
            }
            for rel, grp in rel_groups.items()
        },
    }


def analyze_scaffold_series(
    obs: List[PairedObservation],
    endpoint_id: str,
) -> Dict[str, Any]:
    """Scaffold and series analysis."""
    series_groups: Dict[str, List[PairedObservation]] = {}
    for o in obs:
        key = o.series_label or o.scaffold_hash or "UNKNOWN_SERIES"
        series_groups.setdefault(key, []).append(o)

    result = {
        "n_series": len(series_groups),
        "overall_n": len(obs),
        "series_details": {},
    }

    for series_key, grp in series_groups.items():
        n = len(grp)
        within_spearman = None
        if n >= MIN_N_SERIES and endpoint_id not in PROB_ENDPOINTS:
            preds = [o.prediction_value for o in grp if o.prediction_value is not None]
            exps = [o.experimental_value for o in grp if o.experimental_value is not None]
            if len(preds) >= MIN_N_SERIES and len(preds) == len(exps):
                within_spearman = _spearman_rho(preds, exps)

        result["series_details"][series_key] = {
            "n": n,
            "within_series_Spearman": (
                round(within_spearman, 4) if within_spearman is not None else None
            ),
            "note": (
                "Insufficient N for within-series Spearman" if n < MIN_N_SERIES else ""
            ),
        }

    return result


def bootstrap_regression(
    obs: List[PairedObservation],
    endpoint_id: str,
    n_boot: int = BOOTSTRAP_N,
    seed: int = BOOTSTRAP_SEED,
) -> Dict[str, Any]:
    """Bootstrap CI for regression metrics at compound level."""
    import random

    pairs = [
        o
        for o in obs
        if o.enters_primary_metrics
        and o.prediction_value is not None
        and o.experimental_value is not None
        and not o.censor_flag
    ]

    if len(pairs) < MIN_N_BOOTSTRAP:
        return {
            "data_insufficient": True,
            "reason": f"N={len(pairs)} < minimum {MIN_N_BOOTSTRAP} for bootstrap",
            "n": len(pairs),
        }

    rng = random.Random(seed)
    mae_boot = []
    spearman_boot = []

    for _ in range(n_boot):
        sample = [rng.choice(pairs) for _ in range(len(pairs))]
        preds = [o.prediction_value for o in sample]
        exps = [o.experimental_value for o in sample]
        errors = [abs(p - e) for p, e in zip(preds, exps)]
        mae_boot.append(statistics.mean(errors))
        rho = _spearman_rho(preds, exps)
        if rho is not None:
            spearman_boot.append(rho)

    mae_boot.sort()
    spearman_boot.sort()

    def ci(vals):
        lo = _percentile(vals, 2.5)
        hi = _percentile(vals, 97.5)
        return {"CI_2.5": round(lo, 4), "CI_97.5": round(hi, 4), "n_boot": n_boot}

    return {
        "data_insufficient": False,
        "n": len(pairs),
        "seed": seed,
        "n_bootstrap": n_boot,
        "MAE_bootstrap": ci(mae_boot),
        "Spearman_bootstrap": ci(spearman_boot) if spearman_boot else None,
    }


def run_full_analysis(
    observations_by_endpoint: Dict[str, List[PairedObservation]],
) -> Dict[str, Any]:
    """Run all analyses for all endpoints. Returns complete metrics dict."""
    results = {}

    for endpoint_id, obs in observations_by_endpoint.items():
        if not obs:
            results[endpoint_id] = {
                "n_total": 0,
                "data_insufficient": True,
                "reason": "No paired observations",
                "coverage_gap": endpoint_id not in PROB_ENDPOINTS and endpoint_id not in LOG10_ENDPOINTS,
            }
            continue

        n_total = len(obs)
        n_prosp = sum(1 for o in obs if o.prospective_evidence_class == "TRUE_PROSPECTIVE")
        n_blind = sum(1 for o in obs if o.prospective_evidence_class == "BLINDED_RETROSPECTIVE")
        n_hist = sum(1 for o in obs if o.prospective_evidence_class == "HISTORICAL_VISIBLE")
        n_primary = sum(1 for o in obs if o.enters_primary_metrics)

        ep_result = {
            "endpoint_id": endpoint_id,
            "n_total": n_total,
            "n_true_prospective": n_prosp,
            "n_blinded_retrospective": n_blind,
            "n_historical_visible": n_hist,
            "n_enters_primary_metrics": n_primary,
        }

        if endpoint_id in PROB_ENDPOINTS:
            ep_result["analysis_type"] = "CLASSIFICATION"
            ep_result["primary_metrics"] = compute_classification_metrics(obs, endpoint_id)
            ep_result["bootstrap"] = {"note": "Bootstrap for classification requires N >= 10"}
        else:
            ep_result["analysis_type"] = "REGRESSION"
            ep_result["primary_metrics"] = compute_regression_metrics(obs, endpoint_id)
            ep_result["bootstrap"] = bootstrap_regression(obs, endpoint_id)

        ep_result["ad_analysis"] = analyze_ad_stratification(obs, endpoint_id)
        ep_result["reliability_analysis"] = analyze_reliability_stratification(obs, endpoint_id)
        ep_result["scaffold_series_analysis"] = analyze_scaffold_series(obs, endpoint_id)

        results[endpoint_id] = ep_result

    return results
