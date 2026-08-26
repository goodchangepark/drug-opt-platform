"""Stage 4C-3B Conformal Recalibration, Statistical Acceptance & Uncertainty Governance Audit Script.

Scientific Rules & Governance:
1. Complete separation of Calibration Data Provenance (EXTERNAL, INTERNAL, TRAINING_OVERLAP_UNKNOWN, UNAVAILABLE)
   from Conformal Calibration Quality (VALIDATED, BORDERLINE, UNDERCOVERED, OVERCOVERED, INSUFFICIENT_N, INVALID, UNAVAILABLE).
2. Coverage Acceptance Rules: Exact binomial tests (Clopper-Pearson 95% CI, binomial SE, deviation, p-value).
   No arbitrary +-2% rule.
3. Minimum Data Policy: N_eval >= 30 and N_cal >= 30 required for statistical coverage validation. Small N (e.g. Caco-2 N=17) classified INSUFFICIENT_N.
4. Independent Split Conformal: Quantiles computed on calibration data ONLY. Evaluated on independent evaluation data ONLY.
5. Interval Utility: Median width, dynamic range ratio, ratio to MAE. Flag UNINFORMATIVE_INTERVAL if width is excessive.
6. Conditional AD Coverage: Stratified coverage for IN_DOMAIN, BORDERLINE, OUT_OF_DOMAIN.
7. Classification Set Efficiency: Empirical coverage, singleton rate, ambiguous rate, empty rate.
"""

from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import binomtest

sys.path.insert(0, os.path.abspath("."))

from backend.admet_predictor import applicability_domain, predict_endpoint
from backend.standardizer import standardize_molecule


def get_canonical_smiles(smiles: str) -> str:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        std = standardize_molecule(smiles)
        return std.get("canonical_smiles", Chem.MolToSmiles(mol, canonical=True))
    except Exception:
        return ""


def binomial_coverage_acceptance(
    k_covered: int,
    n_eval: int,
    nominal_coverage: float = 0.90,
    min_eval_n: int = 30,
    alpha_sig: float = 0.05,
) -> dict[str, Any]:
    """Perform exact binomial coverage validation and calculate statistical uncertainty."""
    if n_eval <= 0:
        return {
            "quality_status": "UNAVAILABLE",
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
        quality_status = "INSUFFICIENT_N"
        is_validated = False
        message = f"Evaluation sample size N={n_eval} < {min_eval_n}; cannot statistically validate coverage."
    elif empirical < nominal_coverage and (nominal_coverage > ci_high or p_less < alpha_sig):
        quality_status = "UNDERCOVERED"
        is_validated = False
        message = f"Empirical coverage ({empirical:.1%}) is significantly below nominal {nominal_coverage:.1%} (p={p_less:.4f}, z={z_score:.2f})."
    elif empirical > nominal_coverage and nominal_coverage < ci_low and empirical > 0.96:
        quality_status = "OVERCOVERED"
        is_validated = False
        message = f"Empirical coverage ({empirical:.1%}) is significantly overcovered."
    elif ci_low <= nominal_coverage <= ci_high:
        quality_status = "VALIDATED"
        is_validated = True
        message = f"Empirical coverage ({empirical:.1%}) statistically validated within 95% CI [{ci_low:.1%}, {ci_high:.1%}]."
    else:
        quality_status = "BORDERLINE"
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


def evaluate_ad_stratified(
    y_true: np.ndarray,
    y_pred_or_prob: np.ndarray,
    ad_list: list[str],
    q_or_t: float,
    endpoint_type: str = "REGRESSION",
    nominal: float = 0.90,
    min_n: int = 15,
) -> dict[str, Any]:
    stratified = {}
    for domain in ["IN_DOMAIN", "BORDERLINE", "OUT_OF_DOMAIN"]:
        indices = [i for i, ad in enumerate(ad_list) if ad == domain]
        n_strat = len(indices)
        if n_strat == 0:
            stratified[domain] = {"n": 0, "empirical_coverage": None, "quality": "NO_DATA"}
            continue
        if endpoint_type == "REGRESSION":
            hits = sum(1 for i in indices if abs(y_true[i] - y_pred_or_prob[i]) <= q_or_t)
        else:
            t_low = 1.0 - q_or_t
            t_high = q_or_t
            hits = 0
            for i in indices:
                p = y_pred_or_prob[i]
                yt = y_true[i]
                pset = set()
                if p >= t_low:
                    pset.add(1)
                if p <= t_high:
                    pset.add(0)
                if yt in pset:
                    hits += 1
        emp = round(hits / n_strat, 4)
        if n_strat < min_n:
            stratified[domain] = {
                "n": n_strat,
                "covered_n": hits,
                "empirical_coverage": emp,
                "quality": "INSUFFICIENT_N",
                "message": f"Stratum N={n_strat} < {min_n}",
            }
        else:
            test_res = binomial_coverage_acceptance(hits, n_strat, nominal_coverage=nominal, min_eval_n=min_n)
            stratified[domain] = {
                "n": n_strat,
                "covered_n": hits,
                "empirical_coverage": emp,
                "quality": test_res["quality_status"],
                "sampling_uncertainty_se": test_res["sampling_uncertainty_se"],
                "confidence_interval_95": test_res["confidence_interval_95"],
            }
    return stratified


def run_comprehensive_audit():
    print("================================================================")
    print("STAGE 4C-3B CONFORMAL RECALIBRATION & GOVERNANCE AUDIT")
    print("================================================================")

    results = {}
    base_dir = Path(__file__).resolve().parent.parent

    # -------------------------------------------------------------
    # 1. HLM Intrinsic Clearance
    # -------------------------------------------------------------
    print("\n--- 1. Auditing HLM Intrinsic Clearance ---")
    openadmet_train = pd.read_csv(base_dir / "models/openadmet/microsomal_clearance/X_train.csv")
    biogen_val = pd.read_csv(base_dir / "models/openadmet/validation/biogen_public_3521.csv")

    hlm_biogen = biogen_val.dropna(subset=["SMILES", "LOG HLM_CLint (mL/min/kg)"]).copy()
    if len(hlm_biogen) > 500:
        hlm_biogen = hlm_biogen.iloc[:500]

    train_smiles = set(openadmet_train["OPENADMET_CANONICAL_SMILES"].dropna())
    val_smiles = [get_canonical_smiles(s) for s in hlm_biogen["SMILES"]]
    hlm_overlap = len(set(val_smiles).intersection(train_smiles))

    y_true_hlm, y_pred_hlm, ad_hlm = [], [], []
    for _, row in hlm_biogen.iterrows():
        sm = str(row["SMILES"]).strip()
        if not sm or sm == "nan":
            continue
        try:
            true_log_clint = float(row["LOG HLM_CLint (mL/min/kg)"])
            pred = predict_endpoint(sm, "HLM intrinsic clearance")
            if pred and pred.get("status") == "COMPLETE" and pred.get("predicted_value") is not None:
                log_pred = float(pred["predicted_value"])
                y_true_hlm.append(true_log_clint)
                y_pred_hlm.append(log_pred)
                ad_hlm.append(pred["applicability_domain"]["classification"])
        except Exception:
            continue

    y_true_hlm = np.array(y_true_hlm)
    y_pred_hlm = np.array(y_pred_hlm)
    abs_err_hlm = np.abs(y_true_hlm - y_pred_hlm)

    n_cal_hlm = len(abs_err_hlm) // 2
    n_eval_hlm = len(abs_err_hlm) - n_cal_hlm

    cal_err_hlm = abs_err_hlm[:n_cal_hlm]
    eval_err_hlm = abs_err_hlm[n_cal_hlm:]
    eval_ad_hlm = ad_hlm[n_cal_hlm:]
    eval_y_true_hlm = y_true_hlm[n_cal_hlm:]
    eval_y_pred_hlm = y_pred_hlm[n_cal_hlm:]

    q80_hlm = float(np.quantile(cal_err_hlm, 0.80))
    q90_hlm = float(np.quantile(cal_err_hlm, 0.90))
    q95_hlm = float(np.quantile(cal_err_hlm, 0.95))

    hits90_hlm = int(np.sum(eval_err_hlm <= q90_hlm))
    cov_test_hlm = binomial_coverage_acceptance(hits90_hlm, n_eval_hlm, nominal_coverage=0.90)

    dyn_range_hlm = round(float(np.max(eval_y_true_hlm) - np.min(eval_y_true_hlm)), 3)
    width_90_hlm = round(2.0 * q90_hlm, 3)
    mae_hlm = round(float(np.mean(eval_err_hlm)), 3)
    rel_width_hlm = round(width_90_hlm / dyn_range_hlm, 3) if dyn_range_hlm > 0 else None

    strat_hlm = evaluate_ad_stratified(eval_y_true_hlm, eval_y_pred_hlm, eval_ad_hlm, q90_hlm, "REGRESSION")

    results["HLM intrinsic clearance"] = {
        "endpoint": "HLM intrinsic clearance",
        "data_provenance": "EXTERNAL" if hlm_overlap == 0 else "INTERNAL",
        "calibration_quality": cov_test_hlm["quality_status"],
        "dataset_name": "Biogen Public ADME Prospective Benchmark",
        "local_path": "models/openadmet/validation/biogen_public_3521.csv",
        "calibration_n": n_cal_hlm,
        "evaluation_n": n_eval_hlm,
        "smiles_overlap_with_training": hlm_overlap,
        "nominal_coverage": 0.90,
        "empirical_coverage": cov_test_hlm["empirical_coverage"],
        "expected_sampling_uncertainty_se": cov_test_hlm["sampling_uncertainty_se"],
        "confidence_interval_95": cov_test_hlm["confidence_interval_95"],
        "deviation": cov_test_hlm["deviation"],
        "z_score": cov_test_hlm["z_score"],
        "p_value": cov_test_hlm["p_value_two_sided"],
        "quantiles": {"0.80": round(q80_hlm, 3), "0.90": round(q90_hlm, 3), "0.95": round(q95_hlm, 3)},
        "empirical_coverage_levels": {
            "0.80": round(float(np.mean(eval_err_hlm <= q80_hlm)), 3),
            "0.90": round(float(np.mean(eval_err_hlm <= q90_hlm)), 3),
            "0.95": round(float(np.mean(eval_err_hlm <= q95_hlm)), 3),
        },
        "interval_utility": {
            "mean_interval_width_90": width_90_hlm,
            "median_interval_width_90": width_90_hlm,
            "dynamic_range": dyn_range_hlm,
            "relative_interval_width": rel_width_hlm,
            "ratio_to_mae": round(width_90_hlm / mae_hlm, 2) if mae_hlm > 0 else None,
            "utility_status": "INFORMATIVE" if (rel_width_hlm and rel_width_hlm <= 1.2) else "UNINFORMATIVE_INTERVAL",
        },
        "conditional_coverage_by_ad": strat_hlm,
        "mae": mae_hlm,
        "rmse": round(float(np.sqrt(np.mean(eval_err_hlm ** 2))), 3),
        "conformal_governance_message": cov_test_hlm["message"],
    }
    print(f"HLM: Provenance={results['HLM intrinsic clearance']['data_provenance']}, Quality={results['HLM intrinsic clearance']['calibration_quality']}, EmpCov={cov_test_hlm['empirical_coverage']:.1%}, SE=+-{cov_test_hlm['sampling_uncertainty_se']:.1%}")

    # -------------------------------------------------------------
    # 2. Caco-2 Permeability
    # -------------------------------------------------------------
    print("\n--- 2. Auditing Caco-2 Permeability ---")
    caco2_val = pd.read_csv(base_dir / "models/admetica/validation/caco2_external_34.csv").dropna(subset=["SMILES", "LogPapp(derived)"])
    caco2_train = pd.read_csv(base_dir / "models/admetica/caco2/training.csv")

    train_caco2_smiles = set(caco2_train["Drug"].dropna())
    val_caco2_smiles = [get_canonical_smiles(s) for s in caco2_val["SMILES"]]
    caco2_overlap = len(set(val_caco2_smiles).intersection(train_caco2_smiles))

    y_true_caco2, y_pred_caco2, ad_caco2 = [], [], []
    for _, row in caco2_val.iterrows():
        sm = str(row["SMILES"]).strip()
        if not sm or sm == "nan":
            continue
        try:
            true_val = float(row["LogPapp(derived)"])
            pred = predict_endpoint(sm, "Permeability")
            if pred and pred.get("status") == "COMPLETE" and pred.get("predicted_value") is not None:
                y_true_caco2.append(true_val)
                y_pred_caco2.append(float(pred["predicted_value"]))
                ad_caco2.append(pred["applicability_domain"]["classification"])
        except Exception:
            continue

    y_true_caco2 = np.array(y_true_caco2)
    y_pred_caco2 = np.array(y_pred_caco2)
    abs_err_caco2 = np.abs(y_true_caco2 - y_pred_caco2)

    n_cal_caco2 = len(abs_err_caco2) // 2
    n_eval_caco2 = len(abs_err_caco2) - n_cal_caco2

    cal_err_caco2 = abs_err_caco2[:n_cal_caco2]
    eval_err_caco2 = abs_err_caco2[n_cal_caco2:]
    eval_ad_caco2 = ad_caco2[n_cal_caco2:]
    eval_y_true_caco2 = y_true_caco2[n_cal_caco2:]
    eval_y_pred_caco2 = y_pred_caco2[n_cal_caco2:]

    q80_caco2 = float(np.quantile(cal_err_caco2, 0.80))
    q90_caco2 = float(np.quantile(cal_err_caco2, 0.90))
    q95_caco2 = float(np.quantile(cal_err_caco2, 0.95))

    hits90_caco2 = int(np.sum(eval_err_caco2 <= q90_caco2))
    cov_test_caco2 = binomial_coverage_acceptance(hits90_caco2, n_eval_caco2, nominal_coverage=0.90, min_eval_n=30)

    dyn_range_caco2 = round(float(np.max(eval_y_true_caco2) - np.min(eval_y_true_caco2)), 3)
    width_90_caco2 = round(2.0 * q90_caco2, 3)
    mae_caco2 = round(float(np.mean(eval_err_caco2)), 3)
    rel_width_caco2 = round(width_90_caco2 / dyn_range_caco2, 3) if dyn_range_caco2 > 0 else None

    strat_caco2 = evaluate_ad_stratified(eval_y_true_caco2, eval_y_pred_caco2, eval_ad_caco2, q90_caco2, "REGRESSION")

    results["Permeability"] = {
        "endpoint": "Permeability",
        "data_provenance": "EXTERNAL" if caco2_overlap == 0 else "INTERNAL",
        "calibration_quality": cov_test_caco2["quality_status"],
        "dataset_name": "Admetica External Caco-2 Benchmark (34 compounds)",
        "local_path": "models/admetica/validation/caco2_external_34.csv",
        "calibration_n": n_cal_caco2,
        "evaluation_n": n_eval_caco2,
        "smiles_overlap_with_training": caco2_overlap,
        "nominal_coverage": 0.90,
        "empirical_coverage": cov_test_caco2["empirical_coverage"],
        "expected_sampling_uncertainty_se": cov_test_caco2["sampling_uncertainty_se"],
        "confidence_interval_95": cov_test_caco2["confidence_interval_95"],
        "deviation": cov_test_caco2["deviation"],
        "z_score": cov_test_caco2["z_score"],
        "p_value": cov_test_caco2["p_value_two_sided"],
        "quantiles": {"0.80": round(q80_caco2, 3), "0.90": round(q90_caco2, 3), "0.95": round(q95_caco2, 3)},
        "empirical_coverage_levels": {
            "0.80": round(float(np.mean(eval_err_caco2 <= q80_caco2)), 3),
            "0.90": round(float(np.mean(eval_err_caco2 <= q90_caco2)), 3),
            "0.95": round(float(np.mean(eval_err_caco2 <= q95_caco2)), 3),
        },
        "interval_utility": {
            "mean_interval_width_90": width_90_caco2,
            "median_interval_width_90": width_90_caco2,
            "dynamic_range": dyn_range_caco2,
            "relative_interval_width": rel_width_caco2,
            "ratio_to_mae": round(width_90_caco2 / mae_caco2, 2) if mae_caco2 > 0 else None,
            "utility_status": "UNINFORMATIVE_INTERVAL",
        },
        "conditional_coverage_by_ad": strat_caco2,
        "mae": mae_caco2,
        "rmse": round(float(np.sqrt(np.mean(eval_err_caco2 ** 2))), 3),
        "conformal_governance_message": cov_test_caco2["message"],
    }
    print(f"Caco-2: Provenance={results['Permeability']['data_provenance']}, Quality={results['Permeability']['calibration_quality']}, EvalN={n_eval_caco2}, IntervalWidth={width_90_caco2}")

    # -------------------------------------------------------------
    # 3. hERG Liability (Classification)
    # -------------------------------------------------------------
    print("\n--- 3. Auditing hERG Liability ---")
    herg_val = pd.read_csv(base_dir / "models/admetica/validation/safety/chembl37_herg_ic50_no_exact_training_overlap.csv").dropna(subset=["smiles", "label"])
    if len(herg_val) > 500:
        herg_val = herg_val.iloc[:500]
    herg_train = pd.read_csv(base_dir / "models/admetica/safety/herg/training.csv")

    train_herg_smiles = set(herg_train["Smiles"].dropna())
    val_herg_smiles = [get_canonical_smiles(s) for s in herg_val["smiles"]]
    herg_overlap = len(set(val_herg_smiles).intersection(train_herg_smiles))

    y_true_herg, y_prob_herg, ad_herg = [], [], []
    for _, row in herg_val.iterrows():
        sm = str(row["smiles"]).strip()
        if not sm or sm == "nan":
            continue
        try:
            true_label = int(row["label"])
            pred = predict_endpoint(sm, "hERG liability")
            if pred and pred.get("status") == "COMPLETE" and pred.get("probability") is not None:
                y_true_herg.append(true_label)
                y_prob_herg.append(float(pred["probability"]))
                ad_herg.append(pred["applicability_domain"]["classification"])
        except Exception:
            continue

    y_true_herg = np.array(y_true_herg)
    y_prob_herg = np.array(y_prob_herg)

    scores_herg = np.where(y_true_herg == 1, 1.0 - y_prob_herg, y_prob_herg)

    n_cal_herg = len(scores_herg) // 2
    n_eval_herg = len(scores_herg) - n_cal_herg
    cal_scores_herg = scores_herg[:n_cal_herg]
    eval_yt_herg = y_true_herg[n_cal_herg:]
    eval_yp_herg = y_prob_herg[n_cal_herg:]
    eval_ad_herg = ad_herg[n_cal_herg:]

    q90_herg = float(np.quantile(cal_scores_herg, 0.90))
    t_low_herg = 1.0 - q90_herg
    t_high_herg = q90_herg

    covered_herg = 0
    sing_herg = 0
    amb_herg = 0
    emp_herg = 0
    for yt, p in zip(eval_yt_herg, eval_yp_herg):
        pset = set()
        if p >= t_low_herg:
            pset.add(1)
        if p <= t_high_herg:
            pset.add(0)
        if yt in pset:
            covered_herg += 1
        if len(pset) == 1:
            sing_herg += 1
        elif len(pset) > 1:
            amb_herg += 1
        else:
            emp_herg += 1

    cov_test_herg = binomial_coverage_acceptance(covered_herg, n_eval_herg, nominal_coverage=0.90)
    strat_herg = evaluate_ad_stratified(eval_yt_herg, eval_yp_herg, eval_ad_herg, q90_herg, "CLASSIFICATION")

    results["hERG liability"] = {
        "endpoint": "hERG liability",
        "data_provenance": "EXTERNAL" if herg_overlap == 0 else "INTERNAL",
        "calibration_quality": cov_test_herg["quality_status"],
        "dataset_name": "ChEMBL37 Non-overlapping hERG IC50 Benchmark",
        "local_path": "models/admetica/validation/safety/chembl37_herg_ic50_no_exact_training_overlap.csv",
        "calibration_n": n_cal_herg,
        "evaluation_n": n_eval_herg,
        "smiles_overlap_with_training": herg_overlap,
        "nominal_coverage": 0.90,
        "empirical_coverage": cov_test_herg["empirical_coverage"],
        "expected_sampling_uncertainty_se": cov_test_herg["sampling_uncertainty_se"],
        "confidence_interval_95": cov_test_herg["confidence_interval_95"],
        "deviation": cov_test_herg["deviation"],
        "z_score": cov_test_herg["z_score"],
        "p_value": cov_test_herg["p_value_two_sided"],
        "quantile_90": round(q90_herg, 3),
        "threshold_0.90": round(q90_herg, 3),
        "threshold_low": round(t_low_herg, 4),
        "threshold_high": round(t_high_herg, 4),
        "set_efficiency": {
            "singleton_rate": round(sing_herg / n_eval_herg, 3),
            "ambiguous_rate": round(amb_herg / n_eval_herg, 3),
            "empty_rate": round(emp_herg / n_eval_herg, 3),
            "efficiency_status": "HIGH_AMBIGUITY" if (amb_herg / n_eval_herg > 0.50) else "EFFICIENT",
        },
        "conditional_coverage_by_ad": strat_herg,
        "conformal_governance_message": cov_test_herg["message"],
    }
    print(f"hERG: Provenance={results['hERG liability']['data_provenance']}, Quality={results['hERG liability']['calibration_quality']}, EmpCov={cov_test_herg['empirical_coverage']:.1%}")

    # -------------------------------------------------------------
    # 4. CYP Inhibitors (CYP2C9, CYP2D6, CYP3A4)
    # -------------------------------------------------------------
    for cyp_name, val_file, train_file in [
        ("CYP2C9 inhibitor", "models/admetica/validation/cyp/chembl30_2c9_inhibitor.csv", "models/admetica/cyp/cyp2c9-inhibitor/training.csv"),
        ("CYP2D6 inhibitor", "models/admetica/validation/cyp/chembl30_2d6_inhibitor.csv", "models/admetica/cyp/cyp2d6-inhibitor/training.csv"),
        ("CYP3A4 inhibitor", "models/admetica/validation/cyp/chembl30_3a4_inhibitor.csv", "models/admetica/cyp/cyp3a4-inhibitor/training.csv"),
    ]:
        print(f"\n--- 4. Auditing {cyp_name} ---")
        df_val = pd.read_csv(base_dir / val_file).dropna(subset=["smiles", "class"])
        if len(df_val) > 500:
            df_val = df_val.iloc[:500]
        df_train = pd.read_csv(base_dir / train_file)
        tr_smiles = set(df_train["smiles"].dropna())
        val_smiles_list = [get_canonical_smiles(s) for s in df_val["smiles"]]
        cyp_overlap = len(set(val_smiles_list).intersection(tr_smiles))

        y_t, y_p, ad_c = [], [], []
        for _, row in df_val.iterrows():
            sm = str(row["smiles"]).strip()
            if not sm or sm == "nan":
                continue
            try:
                lbl = int(row["class"])
                pred = predict_endpoint(sm, cyp_name)
                if pred and pred.get("status") == "COMPLETE" and pred.get("probability") is not None:
                    y_t.append(lbl)
                    y_p.append(float(pred["probability"]))
                    ad_c.append(pred["applicability_domain"]["classification"])
            except Exception:
                continue

        y_t = np.array(y_t)
        y_p = np.array(y_p)
        scores = np.where(y_t == 1, 1.0 - y_p, y_p)

        n_cal_cyp = len(scores) // 2
        n_eval_cyp = len(scores) - n_cal_cyp
        cal_s = scores[:n_cal_cyp]
        eval_yt = y_t[n_cal_cyp:]
        eval_yp = y_p[n_cal_cyp:]
        eval_ad = ad_c[n_cal_cyp:]

        q90_cyp = float(np.quantile(cal_s, 0.90))
        t_low = 1.0 - q90_cyp
        t_high = q90_cyp

        cov_c = 0
        sing_c = 0
        amb_c = 0
        emp_c = 0
        for yt, p in zip(eval_yt, eval_yp):
            pset = set()
            if p >= t_low:
                pset.add(1)
            if p <= t_high:
                pset.add(0)
            if yt in pset:
                cov_c += 1
            if len(pset) == 1:
                sing_c += 1
            elif len(pset) > 1:
                amb_c += 1
            else:
                emp_c += 1

        cov_test_cyp = binomial_coverage_acceptance(cov_c, n_eval_cyp, nominal_coverage=0.90)
        strat_cyp = evaluate_ad_stratified(eval_yt, eval_yp, eval_ad, q90_cyp, "CLASSIFICATION")

        results[cyp_name] = {
            "endpoint": cyp_name,
            "data_provenance": "EXTERNAL" if cyp_overlap == 0 else "INTERNAL",
            "calibration_quality": cov_test_cyp["quality_status"],
            "dataset_name": f"ChEMBL30 {cyp_name} Benchmark",
            "local_path": val_file,
            "calibration_n": n_cal_cyp,
            "evaluation_n": n_eval_cyp,
            "smiles_overlap_with_training": cyp_overlap,
            "nominal_coverage": 0.90,
            "empirical_coverage": cov_test_cyp["empirical_coverage"],
            "expected_sampling_uncertainty_se": cov_test_cyp["sampling_uncertainty_se"],
            "confidence_interval_95": cov_test_cyp["confidence_interval_95"],
            "deviation": cov_test_cyp["deviation"],
            "z_score": cov_test_cyp["z_score"],
            "p_value": cov_test_cyp["p_value_two_sided"],
            "quantile_90": round(q90_cyp, 3),
            "threshold_0.90": round(q90_cyp, 3),
            "threshold_low": round(t_low, 4),
            "threshold_high": round(t_high, 4),
            "set_efficiency": {
                "singleton_rate": round(sing_c / n_eval_cyp, 3),
                "ambiguous_rate": round(amb_c / n_eval_cyp, 3),
                "empty_rate": round(emp_c / n_eval_cyp, 3),
                "efficiency_status": "HIGH_AMBIGUITY" if (amb_c / n_eval_cyp > 0.50) else "EFFICIENT",
            },
            "conditional_coverage_by_ad": strat_cyp,
            "conformal_governance_message": cov_test_cyp["message"],
        }
        print(f"{cyp_name}: Provenance={results[cyp_name]['data_provenance']}, Quality={results[cyp_name]['calibration_quality']}, EmpCov={cov_test_cyp['empirical_coverage']:.1%}")

    # -------------------------------------------------------------
    # 5. Uncalibrated / Training-Only Endpoints
    # -------------------------------------------------------------
    for uncal_name, path in [
        ("Solubility", "models/admetica/solubility/training.csv"),
        ("Plasma protein binding", "models/admetica/ppbr/training.csv"),
        ("RLM intrinsic clearance", "models/openadmet/microsomal_clearance/X_train.csv"),
        ("MLM intrinsic clearance", "models/openadmet/microsomal_clearance/X_train.csv"),
        ("CYP1A2 inhibitor", "models/admetica/cyp/cyp1a2-inhibitor/training.csv"),
        ("CYP2C19 inhibitor", "models/admetica/cyp/cyp2c19-inhibitor/training.csv"),
        ("CYP2C9 substrate", "models/admetica/cyp/cyp2c9-substrate/training.csv"),
        ("CYP2D6 substrate", "models/admetica/cyp/cyp2d6-substrate/training.csv"),
        ("CYP3A4 substrate", "models/admetica/cyp/cyp3a4-substrate/training.csv"),
        ("P-gp inhibitor", "models/admetica/transporter/pgp-inhibitor/training.csv"),
        ("Ames mutagenicity", "models/admet_ai/training/ames/training.csv"),
        ("DILI clinical liability", "models/admet_ai/training/dili/training.csv"),
    ]:
        prov = "UNAVAILABLE" if "openadmet" in path else "TRAINING_OVERLAP_UNKNOWN"
        results[uncal_name] = {
            "endpoint": uncal_name,
            "data_provenance": prov,
            "calibration_quality": "UNAVAILABLE",
            "dataset_name": "Training data available only; no independent calibration set",
            "source": "Local model training set",
            "local_path": path,
            "calibration_n": 0,
            "evaluation_n": 0,
            "smiles_overlap_with_training": "UNKNOWN",
            "nominal_coverage": 0.90,
            "empirical_coverage": None,
            "conformal_governance_message": f"Independent external calibration set not provided on disk for '{uncal_name}'.",
        }

    out_path = base_dir / "validation/stage4c3b_conformal_audit.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print("\n================================================================")
    print(f"STAGE 4C-3B AUDIT COMPLETE -> Saved to {out_path}")
    print("================================================================")
    return results


if __name__ == "__main__":
    run_comprehensive_audit()
