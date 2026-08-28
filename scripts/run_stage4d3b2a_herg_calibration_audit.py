"""
Stage 4D-3B2A: hERG Calibration, Threshold & Model Quality Audit
=================================================================

MISSION:
Determine why hERG classification showed high sensitivity / poor specificity.
Root cause: threshold, calibration, base-model discrimination, label/assay
heterogeneity, class imbalance, or a combination.

NO adaptive weighting. NO production changes. SHADOW only.
NO UI modifications.

Produces 9 machine-readable JSON artifacts in validation/
"""

from __future__ import annotations

import json
import logging
import sys
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors, rdFingerprintGenerator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
VAL_DIR = ROOT / "validation"
VAL_DIR.mkdir(exist_ok=True)

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger(__name__)

RANDOM_STATE = 42
HERG_THRESHOLD = 10_000.0   # nM (= 10 µM)
PRODUCTION_THRESHOLD = 0.50  # binary classification cutoff
BORDERLINE_FACTOR = 3.0      # borderline ± 1 log unit from threshold

# ──────────────────────────────────────────────────────────────────────────────
# Section 1: Helpers
# ──────────────────────────────────────────────────────────────────────────────

def _safe_div(a: float, b: float, default: float = 0.0) -> float:
    return a / b if b > 0 else default


def bounded_logloss(y_true: np.ndarray, y_prob: np.ndarray, eps: float = 1e-4) -> float:
    p = np.clip(y_prob, eps, 1 - eps)
    return float(-np.mean(y_true * np.log(p) + (1 - y_true) * np.log(1 - p)))


def brier_score(y_true: np.ndarray, y_prob: np.ndarray) -> float:
    return float(np.mean((y_prob - y_true) ** 2))


def compute_metrics(y_true: np.ndarray, y_prob: np.ndarray,
                    threshold: float = 0.5) -> Dict[str, float]:
    y_pred = (y_prob >= threshold).astype(int)
    tp = int(np.sum((y_pred == 1) & (y_true == 1)))
    tn = int(np.sum((y_pred == 0) & (y_true == 0)))
    fp = int(np.sum((y_pred == 1) & (y_true == 0)))
    fn = int(np.sum((y_pred == 0) & (y_true == 1)))
    sens  = _safe_div(tp, tp + fn)
    spec  = _safe_div(tn, tn + fp)
    bacc  = 0.5 * (sens + spec)
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc   = _safe_div(tp * tn - fp * fn, denom)
    auroc = compute_auroc(y_true, y_prob)
    auprc = compute_auprc(y_true, y_prob)
    bs    = brier_score(y_true, y_prob)
    ll    = bounded_logloss(y_true, y_prob)
    ece   = compute_ece(y_true, y_prob)
    return {
        "n": int(len(y_true)),
        "n_pos": int(np.sum(y_true)),
        "n_neg": int(np.sum(1 - y_true)),
        "prevalence": round(float(np.mean(y_true)), 4),
        "threshold": threshold,
        "tp": tp, "tn": tn, "fp": fp, "fn": fn,
        "sensitivity": round(sens, 4),
        "specificity": round(spec, 4),
        "balanced_accuracy": round(bacc, 4),
        "mcc": round(mcc, 4),
        "auroc": round(auroc, 4),
        "auprc": round(auprc, 4),
        "brier_score": round(bs, 4),
        "log_loss": round(ll, 4),
        "ece": round(ece, 4),
    }


def compute_auroc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    thresholds = np.unique(y_score)[::-1]
    tprs, fprs = [0.0], [0.0]
    pos = np.sum(y_true); neg = len(y_true) - pos
    for t in thresholds:
        pred = (y_score >= t)
        tprs.append(float(np.sum(pred & (y_true == 1)) / pos) if pos > 0 else 0.0)
        fprs.append(float(np.sum(pred & (y_true == 0)) / neg) if neg > 0 else 0.0)
    tprs.append(1.0); fprs.append(1.0)
    return float(np.trapezoid(tprs, fprs))


def compute_auprc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Compute AUPRC via trapezoidal integration over precision-recall curve."""
    thresholds = np.sort(np.unique(y_score))[::-1]  # high → low
    pos = float(np.sum(y_true))
    recalls, precs = [0.0], [1.0]
    for t in thresholds:
        pred = (y_score >= t)
        tp = float(np.sum(pred & (y_true == 1)))
        fp = float(np.sum(pred & (y_true == 0)))
        r = tp / pos if pos > 0 else 0.0
        p = tp / (tp + fp) if (tp + fp) > 0 else 1.0
        recalls.append(r)
        precs.append(p)
    recalls.append(1.0)
    precs.append(pos / len(y_true) if len(y_true) > 0 else 0.0)
    rec_arr = np.array(recalls); prec_arr = np.array(precs)
    order = np.argsort(rec_arr)
    return float(np.trapezoid(prec_arr[order], rec_arr[order]))



def compute_ece(y_true: np.ndarray, y_prob: np.ndarray, n_bins: int = 5) -> float:
    bins = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    for i in range(n_bins):
        mask = (y_prob >= bins[i]) & (y_prob < bins[i + 1])
        if i == n_bins - 1:
            mask = (y_prob >= bins[i]) & (y_prob <= bins[i + 1])
        if mask.sum() == 0:
            continue
        acc = float(y_true[mask].mean())
        conf = float(y_prob[mask].mean())
        ece += mask.sum() * abs(acc - conf)
    return ece / len(y_true)


def compute_calibration_bins(y_true: np.ndarray, y_prob: np.ndarray,
                             n_bins: int = 5) -> List[Dict]:
    edges = np.linspace(0, 1, n_bins + 1)
    bins_out = []
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        mask = (y_prob >= lo) & (y_prob < hi)
        if i == n_bins - 1:
            mask = (y_prob >= lo) & (y_prob <= hi)
        n = int(mask.sum())
        if n == 0:
            continue
        yt = y_true[mask]; yp = y_prob[mask]
        obs_frac = float(yt.mean())
        mean_pred = float(yp.mean())
        bins_out.append({
            "range": [round(lo, 2), round(hi, 2)],
            "count": n,
            "mean_predicted_prob": round(mean_pred, 4),
            "observed_fraction": round(obs_frac, 4),
            "calibration_gap": round(abs(mean_pred - obs_frac), 4),
        })
    return bins_out


def find_optimal_thresholds(y_true: np.ndarray, y_prob: np.ndarray) -> Dict[str, Any]:
    """Find MCC-optimal, BAcc-optimal, Youden-optimal thresholds."""
    thresholds = np.unique(y_prob)
    best = {"mcc": -99, "bacc": -99, "youden": -99}
    best_t = {"mcc": 0.5, "bacc": 0.5, "youden": 0.5}
    pos = np.sum(y_true); neg = len(y_true) - pos
    for t in thresholds:
        pred = (y_prob >= t).astype(int)
        tp = np.sum((pred == 1) & (y_true == 1))
        tn = np.sum((pred == 0) & (y_true == 0))
        fp = np.sum((pred == 1) & (y_true == 0))
        fn = np.sum((pred == 0) & (y_true == 1))
        sens = _safe_div(tp, tp + fn)
        spec = _safe_div(tn, tn + fp)
        denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        mcc = _safe_div(tp * tn - fp * fn, denom)
        bacc = 0.5 * (sens + spec)
        youden = sens + spec - 1.0
        if mcc > best["mcc"]:
            best["mcc"] = mcc; best_t["mcc"] = float(t)
        if bacc > best["bacc"]:
            best["bacc"] = bacc; best_t["bacc"] = float(t)
        if youden > best["youden"]:
            best["youden"] = youden; best_t["youden"] = float(t)
    return {"optimal_thresholds": best_t, "optimal_values": {k: round(v, 4) for k, v in best.items()}}


def compute_roc_curve(y_true: np.ndarray, y_prob: np.ndarray, n_points: int = 50) -> List[Dict]:
    thresholds = np.linspace(0, 1, n_points)
    pos = np.sum(y_true); neg = len(y_true) - pos
    curve = []
    for t in thresholds:
        pred = (y_prob >= t).astype(int)
        tp = np.sum((pred == 1) & (y_true == 1))
        tn = np.sum((pred == 0) & (y_true == 0))
        fp = np.sum((pred == 1) & (y_true == 0))
        fn = np.sum((pred == 0) & (y_true == 1))
        curve.append({
            "threshold": round(float(t), 3),
            "sensitivity": round(_safe_div(tp, tp + fn), 4),
            "specificity": round(_safe_div(tn, tn + fp), 4),
            "fpr": round(_safe_div(fp, fp + tn), 4),
        })
    return curve


def get_physcochem(smiles: str) -> Dict[str, Any]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    try:
        mw    = round(float(Descriptors.MolWt(mol)), 1)
        clogp = round(float(Crippen.MolLogP(mol)), 2)
        tpsa  = round(float(Descriptors.TPSA(mol)), 1)
        hba   = int(Descriptors.NumHAcceptors(mol))
        hbd   = int(Descriptors.NumHDonors(mol))
        n_rot = int(Descriptors.NumRotatableBonds(mol))
        n_rings = int(rdMolDescriptors.CalcNumRings(mol))
        n_arom = sum(1 for ring in mol.GetRingInfo().AtomRings()
                     if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring))
        has_basic_n = bool(
            mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O);!$(NS(=O)=O)]")) or
            mol.HasSubstructMatch(Chem.MolFromSmarts("[$([NX3;H2,H1,H0]),$([NX4+])]"))
        )
        formal_charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())
        return {
            "mw": mw, "clogp": clogp, "tpsa": tpsa,
            "hba": hba, "hbd": hbd, "n_rotatable": n_rot,
            "n_rings": n_rings, "n_aromatic_rings": n_arom,
            "has_basic_n": has_basic_n, "formal_charge": formal_charge,
        }
    except Exception:
        return {}


def get_murcko_scaffold(smiles: str) -> str:
    try:
        from rdkit.Chem.Scaffolds import MurckoScaffold
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        sc = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(sc) if sc else ""
    except Exception:
        return ""


def scaffold_aware_split(df: pd.DataFrame, test_frac: float = 0.25,
                         random_state: int = 42) -> Tuple[pd.Index, pd.Index]:
    """Split by scaffold: scaffolds go entirely into calibration or test."""
    np.random.seed(random_state)
    df = df.copy()
    df["_scaffold"] = df["smiles"].apply(get_murcko_scaffold)
    scaffolds = df["_scaffold"].unique()
    np.random.shuffle(scaffolds)
    n_test = max(1, int(len(df) * test_frac))
    test_scaffs, test_count = set(), 0
    for sc in scaffolds:
        if test_count >= n_test:
            break
        mask = df["_scaffold"] == sc
        test_scaffs.add(sc)
        test_count += mask.sum()
    test_idx = df.index[df["_scaffold"].isin(test_scaffs)]
    cal_idx  = df.index[~df["_scaffold"].isin(test_scaffs)]
    return cal_idx, test_idx


def platt_scaling(y_cal: np.ndarray, p_cal: np.ndarray,
                  p_test: np.ndarray) -> np.ndarray:
    """Fit Platt sigmoid (logistic regression on logit) on calibration, apply to test."""
    from scipy.special import logit
    from scipy.optimize import minimize
    logits_cal = logit(np.clip(p_cal, 1e-6, 1 - 1e-6))

    def neg_ll(params):
        a, b = params
        p = 1.0 / (1.0 + np.exp(-(a * logits_cal + b)))
        p = np.clip(p, 1e-7, 1 - 1e-7)
        return -np.mean(y_cal * np.log(p) + (1 - y_cal) * np.log(1 - p))

    res = minimize(neg_ll, x0=[1.0, 0.0], method="Nelder-Mead")
    a, b = res.x
    logits_test = logit(np.clip(p_test, 1e-6, 1 - 1e-6))
    return np.clip(1.0 / (1.0 + np.exp(-(a * logits_test + b))), 1e-4, 1 - 1e-4)


def isotonic_scaling(y_cal: np.ndarray, p_cal: np.ndarray,
                     p_test: np.ndarray) -> np.ndarray:
    """Fit isotonic regression on calibration, apply to test."""
    from sklearn.isotonic import IsotonicRegression
    ir = IsotonicRegression(out_of_bounds="clip")
    ir.fit(p_cal, y_cal)
    return np.clip(ir.predict(p_test), 1e-4, 1 - 1e-4)


def run_model_m1(smiles_list: List[str]) -> List[float]:
    """Run M1 (Admetica D-MPNN) inference."""
    from backend.admet_predictor import predict_endpoint
    probs = []
    for smiles in smiles_list:
        try:
            res = predict_endpoint(smiles, "hERG liability")
            p = res.get("probability", 0.5)
            probs.append(float(p) if p is not None else 0.5)
        except Exception:
            probs.append(0.5)
    return probs


def run_model_m2(smiles_list: List[str]) -> List[float]:
    """Run M2 (Physchem pharmacophore logistic) inference."""
    from backend.multimodel import get_adapters_for_endpoint
    from backend.endpoint_contracts import get_endpoint_contract
    contract = get_endpoint_contract("hERG liability")
    m2_adapter = None
    for a in get_adapters_for_endpoint("hERG liability"):
        if a.model_id == "physchem_herg_v1":
            m2_adapter = a
            break
    if m2_adapter is None:
        raise RuntimeError("M2 physchem_herg_v1 adapter not found")
    probs = []
    for smiles in smiles_list:
        try:
            res = m2_adapter.execute(smiles, contract)
            p = res.probability if res.probability is not None else res.value
            probs.append(float(p) if p is not None else 0.5)
        except Exception:
            probs.append(0.5)
    return probs


# ──────────────────────────────────────────────────────────────────────────────
# Section 2: Load data
# ──────────────────────────────────────────────────────────────────────────────

def load_validation_data() -> pd.DataFrame:
    val_csv = ROOT / "models/admetica/validation/safety/chembl37_herg_ic50_no_exact_training_overlap.csv"
    df = pd.read_csv(val_csv)
    df = df.rename(columns={"smiles": "smiles"})
    df["smiles"] = df["smiles"].astype(str).str.strip()
    # Confirm valid SMILES
    df["valid_mol"] = df["smiles"].apply(lambda s: Chem.MolFromSmiles(s) is not None)
    n_invalid = (~df["valid_mol"]).sum()
    if n_invalid > 0:
        log.warning(f"Dropping {n_invalid} invalid SMILES from validation set")
    df = df[df["valid_mol"]].copy().reset_index(drop=True)
    # Borderline classification
    df["ic50_nM"] = df["median_ic50_nM"]
    df["strong_positive"] = df["ic50_nM"] <= 1000.0
    df["borderline"]      = (df["ic50_nM"] > 1000.0) & (df["ic50_nM"] <= 30000.0)
    df["strong_negative"] = df["ic50_nM"] > 30000.0
    df["ic50_class"] = df.apply(
        lambda r: "STRONG_POS" if r["strong_positive"] else
                  ("BORDERLINE" if r["borderline"] else "STRONG_NEG"), axis=1
    )
    df["log10_ic50"] = np.log10(df["ic50_nM"])
    df["dist_to_cutoff_log10"] = abs(np.log10(df["ic50_nM"] / HERG_THRESHOLD))
    return df


# ──────────────────────────────────────────────────────────────────────────────
# Main audit entry-point
# ──────────────────────────────────────────────────────────────────────────────

def main():
    print("=" * 74)
    print("STAGE 4D-3B2A: hERG CALIBRATION & MODEL QUALITY AUDIT")
    print("=" * 74)

    # ── Load data ──────────────────────────────────────────────────────────
    print("\nLoading validation cohort...")
    df = load_validation_data()
    N = len(df)
    n_pos = int(df["label"].sum())
    n_neg = N - n_pos
    prevalence = round(n_pos / N, 4)
    print(f"Cohort: N={N}, Positives={n_pos} ({100*prevalence:.1f}%), Negatives={n_neg}")
    print(f"IC50 range: {df['ic50_nM'].min():.1f} – {df['ic50_nM'].max():.1f} nM")
    print(f"Borderline (1k-30k nM): {df['borderline'].sum()}, "
          f"Strong pos (≤1k): {df['strong_positive'].sum()}, "
          f"Strong neg (>30k): {df['strong_negative'].sum()}")

    # ── Run M1 and M2 ──────────────────────────────────────────────────────
    print("\nRunning M1 (Admetica D-MPNN)...")
    df["m1_prob"] = run_model_m1(df["smiles"].tolist())
    print("Running M2 (Physchem pharmacophore logistic)...")
    df["m2_prob"] = run_model_m2(df["smiles"].tolist())

    y = df["label"].values.astype(float)
    p1 = df["m1_prob"].values
    p2 = df["m2_prob"].values

    # ── Full cohort raw metrics ────────────────────────────────────────────
    print("\n--- Full cohort raw metrics (production threshold=0.50) ---")
    m1_full = compute_metrics(y, p1, PRODUCTION_THRESHOLD)
    m2_full = compute_metrics(y, p2, PRODUCTION_THRESHOLD)
    static_50 = compute_metrics(y, 0.5 * p1 + 0.5 * p2, PRODUCTION_THRESHOLD)
    print(f"M1 (Admetica D-MPNN): MCC={m1_full['mcc']}, BAcc={m1_full['balanced_accuracy']}, "
          f"Sens={m1_full['sensitivity']}, Spec={m1_full['specificity']}, "
          f"AUROC={m1_full['auroc']}, Brier={m1_full['brier_score']}, LL={m1_full['log_loss']}")
    print(f"M2 (Physchem):        MCC={m2_full['mcc']}, BAcc={m2_full['balanced_accuracy']}, "
          f"Sens={m2_full['sensitivity']}, Spec={m2_full['specificity']}, "
          f"AUROC={m2_full['auroc']}, Brier={m2_full['brier_score']}, LL={m2_full['log_loss']}")
    print(f"50/50 Static:         MCC={static_50['mcc']}, BAcc={static_50['balanced_accuracy']}, "
          f"Sens={static_50['sensitivity']}, Spec={static_50['specificity']}, "
          f"AUROC={static_50['auroc']}, Brier={static_50['brier_score']}, LL={static_50['log_loss']}")

    # ── Applicability domain ───────────────────────────────────────────────
    print("\nComputing physicochemical features and applicability domain...")
    physchems = [get_physcochem(s) for s in df["smiles"].tolist()]
    for k in ["mw", "clogp", "tpsa", "hba", "hbd", "n_rings",
              "n_aromatic_rings", "has_basic_n", "formal_charge"]:
        df[k] = [pc.get(k) for pc in physchems]
    df["scaffold"] = df["smiles"].apply(get_murcko_scaffold)
    # Simple AD: IN_DOMAIN if mw<=800 and -2<=clogp<=7
    df["ad"] = df.apply(
        lambda r: "IN_DOMAIN" if (r["mw"] is not None and r["mw"] <= 800
                                  and r["clogp"] is not None and -2 <= r["clogp"] <= 7)
                  else ("BORDERLINE" if r["mw"] is not None and r["mw"] <= 1000 else "OUT_OF_DOMAIN"), axis=1
    )
    ad_metrics = {}
    for ad_cat in ["IN_DOMAIN", "BORDERLINE", "OUT_OF_DOMAIN"]:
        mask = df["ad"] == ad_cat
        if mask.sum() < 5:
            ad_metrics[ad_cat] = {"n": int(mask.sum()), "skip": "insufficient_n"}
            continue
        ad_metrics[ad_cat] = compute_metrics(y[mask], p1[mask], PRODUCTION_THRESHOLD)
    print(f"AD IN_DOMAIN ({ad_metrics.get('IN_DOMAIN', {}).get('n', 0)}): "
          f"Spec={ad_metrics.get('IN_DOMAIN', {}).get('specificity', 'n/a')}, "
          f"Sens={ad_metrics.get('IN_DOMAIN', {}).get('sensitivity', 'n/a')}")

    # ── Scaffold-aware split (calibration / test) ──────────────────────────
    print("\nBuilding scaffold-aware calibration/test split (75%/25%)...")
    cal_idx, test_idx = scaffold_aware_split(df, test_frac=0.25, random_state=RANDOM_STATE)
    print(f"Calibration: N={len(cal_idx)}, Test (holdout): N={len(test_idx)}")
    y_cal = y[cal_idx]; p1_cal = p1[cal_idx]; p2_cal = p2[cal_idx]
    y_tst = y[test_idx]; p1_tst = p1[test_idx]; p2_tst = p2[test_idx]
    print(f"  Cal prevalence: {y_cal.mean():.3f}, Test prevalence: {y_tst.mean():.3f}")

    # ── ROC / PR curves ────────────────────────────────────────────────────
    print("\nComputing ROC curves...")
    roc_m1 = compute_roc_curve(y, p1)
    roc_m2 = compute_roc_curve(y, p2)

    # ── Threshold optimization (calibration only!) ─────────────────────────
    print("Finding optimal thresholds on calibration set...")
    t_opts_m1 = find_optimal_thresholds(y_cal, p1_cal)
    t_opts_m2 = find_optimal_thresholds(y_cal, p2_cal)
    print(f"M1 optimal thresholds (cal): {t_opts_m1['optimal_thresholds']}")
    print(f"M2 optimal thresholds (cal): {t_opts_m2['optimal_thresholds']}")

    # Evaluate optimal thresholds on UNTOUCHED HOLDOUT ONLY
    threshold_holdout = {}
    for method_name, t_opt in t_opts_m1["optimal_thresholds"].items():
        m = compute_metrics(y_tst, p1_tst, t_opt)
        m["optimal_threshold"] = t_opt
        m["optimized_on"] = "calibration_set"
        m["evaluated_on"] = "untouched_holdout"
        threshold_holdout[f"m1_{method_name}"] = m
    # Also evaluate production threshold on holdout
    threshold_holdout["m1_production_0.50"] = compute_metrics(y_tst, p1_tst, PRODUCTION_THRESHOLD)
    threshold_holdout["m1_production_0.50"]["evaluated_on"] = "untouched_holdout"

    # Determine whether specificity problem is THRESHOLD_DRIVEN or MODEL_DRIVEN
    spec_prod = m1_full["specificity"]
    best_spec_cal = max(
        compute_metrics(y_cal, p1_cal, t)["specificity"]
        for t in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    )
    print(f"\nSpecificity at 0.5: {spec_prod:.4f}")
    print(f"Best specificity achievable (cal grid): {best_spec_cal:.4f}")

    if best_spec_cal > 0.5 and (best_spec_cal - spec_prod) > 0.1:
        threshold_verdict = "THRESHOLD_DRIVEN"
    elif m1_full["auroc"] < 0.65:
        threshold_verdict = "MODEL_DRIVEN"
    elif best_spec_cal > 0.4 and m1_full["auroc"] >= 0.65:
        threshold_verdict = "BOTH"
    else:
        threshold_verdict = "UNCERTAIN"
    print(f"Threshold verdict: {threshold_verdict}")

    # ── Probability calibration (Platt & Isotonic) ─────────────────────────
    print("\nFitting Platt scaling on calibration set...")
    p1_platt_tst  = platt_scaling(y_cal, p1_cal, p1_tst)
    p2_platt_tst  = platt_scaling(y_cal, p2_cal, p2_tst)
    print("Fitting isotonic regression on calibration set...")
    p1_iso_tst    = isotonic_scaling(y_cal, p1_cal, p1_tst)
    p2_iso_tst    = isotonic_scaling(y_cal, p2_cal, p2_tst)

    cal_results_tst = {
        "m1_raw":         compute_metrics(y_tst, p1_tst, PRODUCTION_THRESHOLD),
        "m1_platt":       compute_metrics(y_tst, p1_platt_tst, PRODUCTION_THRESHOLD),
        "m1_isotonic":    compute_metrics(y_tst, p1_iso_tst, PRODUCTION_THRESHOLD),
        "m2_raw":         compute_metrics(y_tst, p2_tst, PRODUCTION_THRESHOLD),
        "m2_platt":       compute_metrics(y_tst, p2_platt_tst, PRODUCTION_THRESHOLD),
        "m2_isotonic":    compute_metrics(y_tst, p2_iso_tst, PRODUCTION_THRESHOLD),
    }
    print("Calibration comparison on holdout:")
    for k, v in cal_results_tst.items():
        print(f"  {k}: Brier={v['brier_score']}, LL={v['log_loss']}, ECE={v['ece']}, "
              f"Spec={v['specificity']}, Sens={v['sensitivity']}")

    # Calibration bins
    cal_bins_m1_raw  = compute_calibration_bins(y, p1)
    cal_bins_m1_plat = compute_calibration_bins(y_tst, p1_platt_tst)
    cal_bins_m1_iso  = compute_calibration_bins(y_tst, p1_iso_tst)
    cal_bins_m2_raw  = compute_calibration_bins(y, p2)

    # ── Fixed blends (calibration-selected, holdout-evaluated) ────────────
    print("\nEvaluating fixed blends on holdout...")
    blend_ratios = [(1.0, 0.0), (0.98, 0.02), (0.95, 0.05),
                    (0.90, 0.10), (0.80, 0.20), (0.50, 0.50)]
    blend_results = {}
    # Select best blend by BAcc on calibration set
    best_blend_key = None; best_blend_bacc = -99.0
    for w1, w2 in blend_ratios:
        p_blend_cal = w1 * p1_cal + w2 * p2_cal
        m_cal = compute_metrics(y_cal, p_blend_cal, PRODUCTION_THRESHOLD)
        p_blend_tst = w1 * p1_tst + w2 * p2_tst
        m_tst = compute_metrics(y_tst, p_blend_tst, PRODUCTION_THRESHOLD)
        key = f"w1={w1:.2f}_w2={w2:.2f}"
        blend_results[key] = {
            "weights": {"m1": w1, "m2": w2},
            "calibration_metrics": m_cal,
            "holdout_metrics": m_tst,
        }
        if m_cal["balanced_accuracy"] > best_blend_bacc:
            best_blend_bacc = m_cal["balanced_accuracy"]
            best_blend_key = key
    print(f"Best blend (cal BAcc selection): {best_blend_key}")
    best_blend_m = blend_results[best_blend_key]["holdout_metrics"]
    print(f"  Holdout: MCC={best_blend_m['mcc']}, BAcc={best_blend_m['balanced_accuracy']}, "
          f"Spec={best_blend_m['specificity']}, Sens={best_blend_m['sensitivity']}")

    # ── Model disagreement analysis ────────────────────────────────────────
    print("\nAnalyzing M1/M2 model disagreement...")
    df["disagreement"] = abs(p1 - p2)
    df["m1_label"]     = (p1 >= PRODUCTION_THRESHOLD).astype(int)
    df["m2_label"]     = (p2 >= PRODUCTION_THRESHOLD).astype(int)
    df["m1_correct"]   = (df["m1_label"] == df["label"]).astype(int)
    df["m2_correct"]   = (df["m2_label"] == df["label"]).astype(int)
    df["m1_error"]     = (df["m1_label"] != df["label"]).astype(int)

    both_correct = int(((df["m1_correct"] == 1) & (df["m2_correct"] == 1)).sum())
    both_wrong   = int(((df["m1_correct"] == 0) & (df["m2_correct"] == 0)).sum())
    m1_only_correct = int(((df["m1_correct"] == 1) & (df["m2_correct"] == 0)).sum())
    m2_only_correct = int(((df["m1_correct"] == 0) & (df["m2_correct"] == 1)).sum())
    print(f"Both correct: {both_correct}, Both wrong: {both_wrong}")
    print(f"M1 correct / M2 wrong: {m1_only_correct}, M2 correct / M1 wrong: {m2_only_correct}")

    # M2 rescue rate: among M1 errors, how often does M2 correct?
    m1_errors_mask = df["m1_error"] == 1
    m2_rescues = int((df["m2_correct"][m1_errors_mask] == 1).sum())
    m1_error_n = int(m1_errors_mask.sum())
    m2_rescue_rate = _safe_div(m2_rescues, m1_error_n)
    print(f"M1 errors: {m1_error_n}, M2 rescues: {m2_rescues}, Rescue rate: {m2_rescue_rate:.3f}")

    # Disagreement vs Brier error
    high_disagree_mask = df["disagreement"] > 0.3
    disagree_m1_brier = brier_score(y[high_disagree_mask], p1[high_disagree_mask]) if high_disagree_mask.sum() > 0 else 0.0
    agree_m1_brier    = brier_score(y[~high_disagree_mask], p1[~high_disagree_mask]) if (~high_disagree_mask).sum() > 0 else 0.0
    print(f"High disagreement (|Δp|>0.3): {high_disagree_mask.sum()} compounds")
    print(f"  M1 Brier in high-disagree: {disagree_m1_brier:.4f} vs agree: {agree_m1_brier:.4f}")

    disagree_result = {
        "both_correct": both_correct,
        "both_wrong": both_wrong,
        "m1_only_correct": m1_only_correct,
        "m2_only_correct": m2_only_correct,
        "m1_error_n": m1_error_n,
        "m2_rescues_of_m1_errors": m2_rescues,
        "m2_rescue_rate": round(m2_rescue_rate, 4),
        "high_disagree_n": int(high_disagree_mask.sum()),
        "high_disagree_m1_brier": round(disagree_m1_brier, 4),
        "agree_m1_brier": round(agree_m1_brier, 4),
        "disagreement_percentiles": {
            "p25": round(float(np.percentile(df["disagreement"], 25)), 4),
            "p50": round(float(np.percentile(df["disagreement"], 50)), 4),
            "p75": round(float(np.percentile(df["disagreement"], 75)), 4),
            "p90": round(float(np.percentile(df["disagreement"], 90)), 4),
        },
    }

    # ── FP/FN error analysis ───────────────────────────────────────────────
    print("\nRunning FP/FN error analysis...")
    m1_pred = (p1 >= PRODUCTION_THRESHOLD).astype(int)
    df["fp_m1"] = ((m1_pred == 1) & (y == 0)).astype(int)
    df["fn_m1"] = ((m1_pred == 0) & (y == 1)).astype(int)
    fp_df = df[df["fp_m1"] == 1].copy()
    fn_df = df[df["fn_m1"] == 1].copy()
    print(f"M1 FP (false positives / inactive called BLOCKER): {len(fp_df)}")
    print(f"M1 FN (false negatives / active missed):           {len(fn_df)}")

    def physchem_summary(subset_df: pd.DataFrame) -> Dict[str, Any]:
        if len(subset_df) == 0:
            return {}
        result = {}
        for col in ["mw", "clogp", "tpsa", "n_aromatic_rings", "n_rings"]:
            if col in subset_df.columns:
                vals = subset_df[col].dropna()
                if len(vals) > 0:
                    result[col] = {
                        "mean": round(float(vals.mean()), 2),
                        "median": round(float(vals.median()), 2),
                        "std": round(float(vals.std()), 2),
                    }
        if "has_basic_n" in subset_df.columns:
            result["pct_basic_n"] = round(float(subset_df["has_basic_n"].mean() * 100), 1)
        if "ic50_class" in subset_df.columns:
            result["ic50_class_dist"] = subset_df["ic50_class"].value_counts().to_dict()
        return result

    fp_summary = physchem_summary(fp_df)
    fn_summary = physchem_summary(fn_df)
    correct_neg_summary = physchem_summary(df[(df["fp_m1"] == 0) & (y == 0)])

    # ── IC50 borderline analysis ───────────────────────────────────────────
    print("\nIC50 borderline boundary analysis...")
    borderline_mask = df["borderline"].values
    strong_pos_mask = df["strong_positive"].values
    strong_neg_mask = df["strong_negative"].values

    def class_metrics(mask: np.ndarray, label: str) -> Dict:
        n = int(mask.sum())
        if n < 5:
            return {"n": n, "skip": "insufficient_n"}
        m = compute_metrics(y[mask], p1[mask], PRODUCTION_THRESHOLD)
        m["ic50_class"] = label
        return m

    ic50_class_metrics = {
        "strong_positive": class_metrics(strong_pos_mask, "STRONG_POS"),
        "borderline": class_metrics(borderline_mask, "BORDERLINE"),
        "strong_negative": class_metrics(strong_neg_mask, "STRONG_NEG"),
    }
    print(f"Strong pos: Spec={ic50_class_metrics['strong_positive'].get('specificity','n/a')}, "
          f"Sens={ic50_class_metrics['strong_positive'].get('sensitivity','n/a')}")
    print(f"Borderline: Spec={ic50_class_metrics['borderline'].get('specificity','n/a')}, "
          f"Sens={ic50_class_metrics['borderline'].get('sensitivity','n/a')}")

    # ── Chemical subgroup analysis ─────────────────────────────────────────
    print("\nChemical subgroup analysis...")
    subgroup_results = {}
    # Basic amine
    for grp, mask_fn in [
        ("Basic_amine_pos", (df["has_basic_n"] == True) & (df["label"] == 1)),
        ("Basic_amine_neg", (df["has_basic_n"] == True) & (df["label"] == 0)),
        ("No_basic_n_pos",  (df["has_basic_n"] == False) & (df["label"] == 1)),
        ("No_basic_n_neg",  (df["has_basic_n"] == False) & (df["label"] == 0)),
        ("High_clogp_pos",  (df["clogp"] >= 4.0) & (df["label"] == 1)),
        ("High_clogp_neg",  (df["clogp"] >= 4.0) & (df["label"] == 0)),
    ]:
        mask = mask_fn.values if hasattr(mask_fn, "values") else np.array(mask_fn)
        if mask.sum() >= 5:
            subgroup_results[grp] = compute_metrics(y[mask], p1[mask], PRODUCTION_THRESHOLD)
            subgroup_results[grp]["n"] = int(mask.sum())

    # ── Series analysis ────────────────────────────────────────────────────
    print("\nScaffold-series analysis...")
    scaffold_counts = df["scaffold"].value_counts()
    populated_scaffolds = scaffold_counts[scaffold_counts >= 5].index.tolist()
    print(f"Scaffolds with N>=5: {len(populated_scaffolds)}")
    series_results = {}
    m2_wins = 0
    for sc in populated_scaffolds:
        smask = df["scaffold"] == sc
        ys  = y[smask]; p1s = p1[smask]; p2s = p2[smask]
        n_s = int(smask.sum())
        m1m = compute_metrics(ys, p1s, PRODUCTION_THRESHOLD)
        m2m = compute_metrics(ys, p2s, PRODUCTION_THRESHOLD)
        m2_wins_here = (m2m["balanced_accuracy"] > m1m["balanced_accuracy"] + 0.05)
        if m2_wins_here:
            m2_wins += 1
        series_results[sc[:60]] = {
            "n": n_s,
            "n_pos": int(ys.sum()),
            "m1_mcc": m1m["mcc"], "m1_bacc": m1m["balanced_accuracy"],
            "m1_spec": m1m["specificity"], "m1_sens": m1m["sensitivity"],
            "m2_mcc": m2m["mcc"], "m2_bacc": m2m["balanced_accuracy"],
            "m2_spec": m2m["specificity"], "m2_sens": m2m["sensitivity"],
            "m2_better_by_5pct_bacc": m2_wins_here,
        }
    print(f"Series where M2 substantially better: {m2_wins} / {len(populated_scaffolds)}")

    # ── Pseudo-project analysis ────────────────────────────────────────────
    print("\nPseudo-project analysis (5 project groups)...")
    np.random.seed(RANDOM_STATE)
    scaffold_list = df["scaffold"].tolist()
    unique_scaffolds = list(dict.fromkeys(scaffold_list))
    np.random.shuffle(unique_scaffolds)
    n_proj = 5
    chunk_size = max(1, len(unique_scaffolds) // n_proj)
    proj_results = {}
    for i in range(n_proj):
        sc_chunk = set(unique_scaffolds[i * chunk_size: (i + 1) * chunk_size])
        pmask = df["scaffold"].isin(sc_chunk)
        yp = y[pmask]; p1p = p1[pmask]; p2p = p2[pmask]
        if len(yp) < 5:
            continue
        m1p = compute_metrics(yp, p1p, PRODUCTION_THRESHOLD)
        m2p = compute_metrics(yp, p2p, PRODUCTION_THRESHOLD)
        proj_results[f"PROJ_{i+1:02d}"] = {
            "n": int(pmask.sum()), "n_pos": int(yp.sum()),
            "m1_mcc": m1p["mcc"], "m1_bacc": m1p["balanced_accuracy"],
            "m1_spec": m1p["specificity"], "m1_sens": m1p["sensitivity"],
            "m2_mcc": m2p["mcc"], "m2_bacc": m2p["balanced_accuracy"],
            "m2_spec": m2p["specificity"], "m2_sens": m2p["sensitivity"],
        }
    for k, v in proj_results.items():
        print(f"  {k}: N={v['n']}, M1 BAcc={v['m1_bacc']}, Spec={v['m1_spec']}")

    # ── M2 role determination ──────────────────────────────────────────────
    print("\nDetermining M2 role...")
    # M2 as calibration supporter: does small blend soften overconfident M1 errors?
    m1_overconf_mask = ((p1 >= 0.9) & (y == 0)) | ((p1 <= 0.1) & (y == 1))
    n_overconf = int(m1_overconf_mask.sum())
    if n_overconf > 0:
        blend_p = 0.95 * p1[m1_overconf_mask] + 0.05 * p2[m1_overconf_mask]
        ll_m1_overconf  = bounded_logloss(y[m1_overconf_mask], p1[m1_overconf_mask])
        ll_blend_overconf = bounded_logloss(y[m1_overconf_mask], blend_p)
        ll_delta = ll_blend_overconf - ll_m1_overconf
        print(f"M1 overconfident errors: {n_overconf}")
        print(f"  LL at 95/5 blend vs raw M1: {ll_blend_overconf:.4f} vs {ll_m1_overconf:.4f} (Δ={ll_delta:.4f})")
    else:
        ll_delta = 0.0

    # ── Final decision ─────────────────────────────────────────────────────
    print("\nDeriving root cause and final decision...")

    # Evidence-based root cause analysis
    root_causes = []
    auroc_m1 = m1_full["auroc"]
    spec_m1  = m1_full["specificity"]

    # Threshold factor
    best_spec_any = 0.0
    for t_name, t_val in t_opts_m1["optimal_thresholds"].items():
        m_tmp = compute_metrics(y_cal, p1_cal, t_val)
        best_spec_any = max(best_spec_any, m_tmp["specificity"])
    if best_spec_any > spec_m1 + 0.10:
        root_causes.append("THRESHOLD")

    # Calibration factor
    if m1_full["ece"] > 0.06:
        root_causes.append("CALIBRATION")

    # Base model discrimination
    if auroc_m1 < 0.75:
        root_causes.append("BASE_MODEL_DISCRIMINATION")

    # Class imbalance
    training_prevalence = 0.8599  # from training data analysis
    if abs(training_prevalence - prevalence) > 0.10:
        root_causes.append("CLASS_IMBALANCE")

    # Label/assay heterogeneity
    if ic50_class_metrics["borderline"].get("n", 0) > 100:
        root_causes.append("LABEL_BOUNDARY_UNCERTAINTY")

    if len(root_causes) == 0:
        root_causes = ["UNCERTAIN"]

    primary_cause = " + ".join(root_causes)
    print(f"Root causes identified: {primary_cause}")

    # Final decision
    m2_rescue_good = m2_rescue_rate > 0.20
    m2_series_wins = m2_wins > 2
    if m1_full["auroc"] < 0.60:
        final_decision = "HERG_NEEDS_BETTER_SECONDARY_MODEL"
    elif not m2_rescue_good and not m2_series_wins:
        if "CALIBRATION" in root_causes and "BASE_MODEL_DISCRIMINATION" not in root_causes:
            final_decision = "HERG_CALIBRATION_UPDATE_CANDIDATE"
        else:
            final_decision = "HERG_FIXED_BLEND_CANDIDATE"
    elif m2_rescue_good and m2_series_wins:
        final_decision = "HERG_ADAPTIVE_RESEARCH_CANDIDATE"
    else:
        final_decision = "HERG_FIXED_BLEND_CANDIDATE"

    print(f"Final hERG decision: {final_decision}")

    # Adaptive gate recommendation
    adaptive_gate = "GO" if final_decision == "HERG_ADAPTIVE_RESEARCH_CANDIDATE" else "NO_GO"
    adaptive_gate_reason = {
        "HERG_NEEDS_BETTER_SECONDARY_MODEL": "M1 AUROC < 0.60; base model discrimination insufficient for adaptive benefit",
        "HERG_CALIBRATION_UPDATE_CANDIDATE": "Primary problem is calibration, not adaptive weighting; fix calibration first",
        "HERG_FIXED_BLEND_CANDIDATE": f"M2 rescue rate {m2_rescue_rate:.3f} and series wins {m2_wins} insufficient to justify adaptive complexity",
        "HERG_ADAPTIVE_RESEARCH_CANDIDATE": "M2 provides reproducible conditional rescue value; adaptive research justified",
        "HERG_ENDPOINT_DATA_REQUALIFICATION_REQUIRED": "Endpoint label quality insufficient",
    }.get(final_decision, "Unknown")

    # ── Build authoritative cohort JSON ───────────────────────────────────
    print("\nBuilding authoritative cohort JSON (N={})...".format(N))
    compounds = []
    for i, row in df.iterrows():
        compounds.append({
            "compound_id": int(i) + 1,
            "canonical_smiles": row["smiles"],
            "scaffold": row["scaffold"],
            "median_ic50_nM": row["ic50_nM"],
            "ic50_class": row["ic50_class"],
            "log10_ic50": round(row["log10_ic50"], 4),
            "dist_to_cutoff_log10": round(row["dist_to_cutoff_log10"], 4),
            "assay_id_count": int(row["assay_id_count"]) if "assay_id_count" in row else None,
            "experimental_label": int(row["label"]),
            "m1_probability": round(float(p1[df.index.get_loc(i)]), 4),
            "m2_probability": round(float(p2[df.index.get_loc(i)]), 4),
            "m1_label": int(p1[df.index.get_loc(i)] >= PRODUCTION_THRESHOLD),
            "m2_label": int(p2[df.index.get_loc(i)] >= PRODUCTION_THRESHOLD),
            "fp_m1": int(row["fp_m1"]),
            "fn_m1": int(row["fn_m1"]),
            "disagreement": round(abs(float(p1[df.index.get_loc(i)]) - float(p2[df.index.get_loc(i)])), 4),
            "applicability_domain": row["ad"],
            "mw": row["mw"], "clogp": row["clogp"], "tpsa": row["tpsa"],
            "has_basic_n": bool(row["has_basic_n"]) if row["has_basic_n"] is not None else None,
            "n_aromatic_rings": int(row["n_aromatic_rings"]) if row["n_aromatic_rings"] is not None else None,
            "formal_charge": int(row["formal_charge"]) if row["formal_charge"] is not None else None,
            "split_assignment": "calibration" if i in cal_idx else "test",
        })

    cohort_json = {
        "endpoint": "safety_herg_blocker_prob",
        "policy_version": "stage4d3b2a-herg-audit-v1",
        "n_compounds": N,
        "n_positive": n_pos,
        "n_negative": n_neg,
        "prevalence": prevalence,
        "training_prevalence": 0.8599,
        "label_definition": "positive when median IC50 <= 10000 nM (pIC50 >= 5.0)",
        "source_dataset": "ChEMBL37 hERG IC50 overlap-excluded",
        "ic50_borderline_n": int(df["borderline"].sum()),
        "ic50_strong_positive_n": int(df["strong_positive"].sum()),
        "ic50_strong_negative_n": int(df["strong_negative"].sum()),
        "calibration_n": len(cal_idx),
        "test_n": len(test_idx),
        "split_method": "scaffold_aware_25pct_test",
        "compounds": compounds,
    }

    # ── Write all 9 JSON artifacts ─────────────────────────────────────────
    print("\nWriting all 9 Stage 4D-3B2A JSON artifacts to validation/...")

    artifacts = {
        "stage4d3b2a_authoritative_cohort.json": cohort_json,
        "stage4d3b2a_model_metrics.json": {
            "endpoint": "safety_herg_blocker_prob",
            "policy_version": "stage4d3b2a-herg-audit-v1",
            "m1_model": {
                "model_id": "admetica_safety_herg",
                "model_name": "Admetica Chemprop human hERG blocker liability",
                "model_version": "admetica-d4f7056-herg-chemprop-v2.1",
                "model_family": "admetica",
                "training_n": 22248,
                "training_prevalence": 0.8599,
                "training_n_pos": 19130,
                "training_n_neg": 3118,
            },
            "m2_model": {
                "model_id": "physchem_herg_v1",
                "model_name": "Physicochemical Basic Center hERG Blocker Classifier",
                "model_version": "physchem-herg-v1.0",
                "model_family": "pharmacophore_logistic",
                "architecture": "logistic regression on cLogP, MW, TPSA, basic_N, aromatic_rings",
                "training_dataset": "Wang et al. hERG Blocker Compilation (N=22,249)",
            },
            "full_cohort_metrics": {
                "m1_core": m1_full,
                "m2_shadow": m2_full,
                "static_50_50": static_50,
            },
            "applicability_domain_metrics": ad_metrics,
            "ic50_class_metrics": ic50_class_metrics,
            "subgroup_metrics": subgroup_results,
        },
        "stage4d3b2a_calibration.json": {
            "endpoint": "safety_herg_blocker_prob",
            "calibration_n": len(cal_idx),
            "test_n": len(test_idx),
            "split_method": "scaffold_aware_25pct_test",
            "full_cohort_ece": {
                "m1_raw": round(compute_ece(y, p1), 4),
                "m2_raw": round(compute_ece(y, p2), 4),
            },
            "holdout_calibration_comparison": cal_results_tst,
            "calibration_bins": {
                "m1_raw": cal_bins_m1_raw,
                "m1_platt": cal_bins_m1_plat,
                "m1_isotonic": cal_bins_m1_iso,
                "m2_raw": cal_bins_m2_raw,
            },
            "calibration_dominance_assessment": {
                "m1_ece": round(compute_ece(y, p1), 4),
                "calibration_improves_specificity": cal_results_tst["m1_platt"]["specificity"] > m1_full["specificity"],
                "calibration_improvement_magnitude_spec": round(
                    cal_results_tst["m1_platt"]["specificity"] - m1_full["specificity"], 4),
                "calibration_improvement_magnitude_brier": round(
                    m1_full["brier_score"] - cal_results_tst["m1_platt"]["brier_score"], 4),
            },
        },
        "stage4d3b2a_threshold_audit.json": {
            "endpoint": "safety_herg_blocker_prob",
            "production_threshold": PRODUCTION_THRESHOLD,
            "m1_optimal_thresholds_cal": t_opts_m1,
            "m2_optimal_thresholds_cal": t_opts_m2,
            "holdout_evaluation": threshold_holdout,
            "roc_curve_m1": roc_m1,
            "roc_curve_m2": roc_m2,
            "specificity_at_production_threshold": spec_m1,
            "best_specificity_achievable_cal_grid": round(best_spec_any, 4),
            "threshold_verdict": threshold_verdict,
        },
        "stage4d3b2a_error_analysis.json": {
            "endpoint": "safety_herg_blocker_prob",
            "m1_fp_count": len(fp_df),
            "m1_fn_count": len(fn_df),
            "fp_physchem_summary": fp_summary,
            "fn_physchem_summary": fn_summary,
            "correct_negative_physchem_summary": correct_neg_summary,
            "ic50_class_error_rates": {
                "strong_positive": ic50_class_metrics["strong_positive"],
                "borderline": ic50_class_metrics["borderline"],
                "strong_negative": ic50_class_metrics["strong_negative"],
            },
            "borderline_fraction_of_fp": round(
                fp_df["borderline"].mean() if len(fp_df) > 0 else 0.0, 4),
            "borderline_fraction_of_fn": round(
                fn_df["borderline"].mean() if len(fn_df) > 0 else 0.0, 4),
        },
        "stage4d3b2a_disagreement.json": {
            "endpoint": "safety_herg_blocker_prob",
            "model_complementarity": disagree_result,
            "m2_rescue_rate": round(m2_rescue_rate, 4),
            "m2_role_assessment": "CALIBRATION_SUPPORTING" if m2_rescue_rate < 0.15 else
                                   ("PARTIALLY_COMPLEMENTARY" if m2_rescue_rate < 0.35 else "COMPLEMENTARY"),
        },
        "stage4d3b2a_fixed_blend.json": {
            "endpoint": "safety_herg_blocker_prob",
            "production_threshold": PRODUCTION_THRESHOLD,
            "blend_selection_metric": "balanced_accuracy",
            "blend_selected_on": "calibration_set",
            "blend_evaluated_on": "untouched_holdout",
            "best_blend": best_blend_key,
            "blend_results": blend_results,
        },
        "stage4d3b2a_series_analysis.json": {
            "endpoint": "safety_herg_blocker_prob",
            "n_scaffolds_analyzed": len(series_results),
            "min_n_threshold": 5,
            "m2_wins_substantial": m2_wins,
            "project_analysis": proj_results,
            "series_analysis": series_results,
            "m2_series_verdict": "SERIES_ADVANTAGE_PRESENT" if m2_wins > 2 else "NO_REPRODUCIBLE_SERIES_ADVANTAGE",
        },
        "stage4d3b2a_final_decision.json": {
            "endpoint": "safety_herg_blocker_prob",
            "policy_version": "stage4d3b2a-herg-audit-v1",
            "m1_model": {
                "model_id": "admetica_safety_herg",
                "role": "CORE",
                "contribution_status": "CORE_PRIMARY",
            },
            "m2_model": {
                "model_id": "physchem_herg_v1",
                "role": "SHADOW_ONLY",
                "contribution_status": "CALIBRATION_SUPPORTING" if m2_rescue_rate < 0.15 else "PARTIALLY_COMPLEMENTARY",
            },
            "full_cohort_metrics_at_production_threshold": m1_full,
            "root_cause_analysis": {
                "causes": root_causes,
                "primary_cause": primary_cause,
                "threshold_verdict": threshold_verdict,
                "auroc_m1": auroc_m1,
                "specificity_at_0.5": spec_m1,
                "best_specificity_achievable_calibration": round(best_spec_any, 4),
                "m2_rescue_rate": round(m2_rescue_rate, 4),
                "m2_series_wins": m2_wins,
                "training_vs_eval_prevalence_shift": round(abs(0.8599 - prevalence), 4),
                "borderline_ic50_fraction": round(df["borderline"].mean(), 4),
            },
            "scientific_decision": final_decision,
            "consensus_mode": "SHADOW",
            "production_threshold": PRODUCTION_THRESHOLD,
            "threshold_recommendation": "KEEP_CURRENT_THRESHOLD" if threshold_verdict in ("MODEL_DRIVEN", "UNCERTAIN")
                                         else "THRESHOLD_REVALIDATION_REQUIRED",
            "adaptive_weighting_gate": adaptive_gate,
            "adaptive_weighting_gate_reason": adaptive_gate_reason,
            "assay_heterogeneity": "ASSAY_HETEROGENEITY_PRESENT",
            "assay_heterogeneity_detail": (
                "Training labels pool patch-clamp and radioligand binding assays (Wang et al.). "
                "Assay type field not retained per row. Labels must be treated as heterogeneous "
                "screening liability, not a single standardized functional assay."
            ),
        },
    }

    for fname, data in artifacts.items():
        out_path = VAL_DIR / fname
        with open(out_path, "w") as f:
            json.dump(data, f, indent=2, default=str)
        print(f"  ✓ {fname}")

    print("\n✓ All 9 Stage 4D-3B2A JSON artifacts written to validation/")
    print("\nFINAL SUMMARY")
    print(f"  M1 AUROC: {auroc_m1}")
    print(f"  M1 Specificity (0.5): {spec_m1}")
    print(f"  Root causes: {primary_cause}")
    print(f"  Final decision: {final_decision}")
    print(f"  Adaptive gate: {adaptive_gate}")


if __name__ == "__main__":
    main()
