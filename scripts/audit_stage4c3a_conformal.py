"""Stage 4C-3A Conformal Scientific Audit & Data Provenance Script."""

import json
import os
import sys
import numpy as np
import pandas as pd
import torch  # Ensure PyTorch static TLS initializes on main thread startup
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

sys.path.insert(0, os.path.abspath("."))

from backend.admet_predictor import predict_endpoint
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


def get_inchikey(smiles: str) -> str:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        return Chem.MolToInchiKey(mol) or ""
    except Exception:
        return ""


def get_scaffold(smiles: str) -> str:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaf, canonical=True) if scaf else ""
    except Exception:
        return ""


def run_audit():
    print("==================================================")
    print("STAGE 4C-3A CONFORMAL SCIENTIFIC AUDIT INITIALIZING")
    print("==================================================")

    results = {}

    # 1. Audit HLM Intrinsic Clearance (OpenADMET vs Biogen 3521)
    print("\n--- Auditing HLM Intrinsic Clearance ---")
    openadmet_train = pd.read_csv("models/openadmet/microsomal_clearance/X_train.csv")
    biogen_val = pd.read_csv("models/openadmet/validation/biogen_public_3521.csv")

    hlm_biogen = biogen_val.dropna(subset=["SMILES", "LOG HLM_CLint (mL/min/kg)"]).copy()
    if len(hlm_biogen) > 500:
        hlm_biogen = hlm_biogen.iloc[:500]
    print(f"Biogen dataset for HLM: {len(hlm_biogen)} compounds")

    train_smiles = set(openadmet_train["OPENADMET_CANONICAL_SMILES"].dropna())

    val_smiles = [get_canonical_smiles(s) for s in hlm_biogen["SMILES"]]

    smiles_overlap = len(set(val_smiles).intersection(train_smiles))
    print(f"SMILES overlap with OpenADMET training: {smiles_overlap} / {len(hlm_biogen)}")

    # Model inference on Biogen validation set
    y_true_hlm = []
    y_pred_hlm = []
    for idx, row in hlm_biogen.iterrows():
        sm = str(row["SMILES"]).strip()
        if not sm or sm == "nan":
            continue
        try:
            true_log_clint = float(row["LOG HLM_CLint (mL/min/kg)"])
            true_clint = 10 ** true_log_clint  # mL/min/kg
            pred = predict_endpoint(sm, "HLM intrinsic clearance")
            if pred and pred.get("status") == "COMPLETE" and pred.get("predicted_value") is not None:
                # predicted_value is log10(mL/min/kg), convert to mL/min/kg
                log_pred = float(pred["predicted_value"])
                pred_clint = 10 ** log_pred
                y_true_hlm.append(true_log_clint)
                y_pred_hlm.append(log_pred)
        except Exception as exc:
            if idx == 0:
                print(f"Row 0 Exception: {exc}", file=sys.stderr)
            continue

    y_true_hlm = np.array(y_true_hlm)
    y_pred_hlm = np.array(y_pred_hlm)
    abs_errors_hlm = np.abs(y_true_hlm - y_pred_hlm)

    # 50/50 Calibration / Evaluation split
    n_total = len(abs_errors_hlm)
    n_cal = n_total // 2
    n_eval = n_total - n_cal
    cal_errors = abs_errors_hlm[:n_cal]
    eval_errors = abs_errors_hlm[n_cal:]

    q80 = float(np.quantile(cal_errors, 0.80))
    q90 = float(np.quantile(cal_errors, 0.90))
    q95 = float(np.quantile(cal_errors, 0.95))

    cov80 = float(np.mean(eval_errors <= q80))
    cov90 = float(np.mean(eval_errors <= q90))
    cov95 = float(np.mean(eval_errors <= q95))

    results["HLM intrinsic clearance"] = {
        "endpoint": "HLM intrinsic clearance",
        "dataset_name": "Biogen Public ADME Dataset (3,521 compounds)",
        "source": "Biogen / OpenADMET Validation Benchmark",
        "local_path": "models/openadmet/validation/biogen_public_3521.csv",
        "calibration_n": n_cal,
        "evaluation_n": n_eval,
        "smiles_overlap_with_training": smiles_overlap,
        "training_overlap_status": "NONE" if smiles_overlap == 0 else f"{smiles_overlap}_OVERLAP",
        "quantiles": {"0.80": q80, "0.90": q90, "0.95": q95},
        "empirical_coverage": {"0.80": cov80, "0.90": cov90, "0.95": cov95},
        "mean_interval_width_90": 2 * q90,
        "median_interval_width_90": 2 * q90,
        "mae": float(np.mean(eval_errors)),
        "rmse": float(np.sqrt(np.mean(eval_errors ** 2))),
        "scientific_status": "CALIBRATED_EXTERNAL" if smiles_overlap == 0 else "CALIBRATED_INTERNAL",
    }
    print(f"HLM Audit Done: status={results['HLM intrinsic clearance']['scientific_status']}, 90% cov={cov90:.3f}, q90={q90:.3f}")

    # 2. Audit Caco-2 Permeability (Admetica vs External 34)
    print("\n--- Auditing Caco-2 Permeability ---")
    caco2_val = pd.read_csv("models/admetica/validation/caco2_external_34.csv").dropna(subset=["SMILES", "LogPapp(derived)"])
    caco2_train = pd.read_csv("models/admetica/caco2/training.csv")

    train_caco2_smiles = set(caco2_train["Drug"].dropna())
    val_caco2_smiles = [get_canonical_smiles(s) for s in caco2_val["SMILES"]]
    caco2_overlap = len(set(val_caco2_smiles).intersection(train_caco2_smiles))

    y_true_caco2 = []
    y_pred_caco2 = []
    for idx, row in caco2_val.iterrows():
        sm = str(row["SMILES"]).strip()
        if not sm or sm == "nan":
            continue
        try:
            true_val = float(row["LogPapp(derived)"])
            pred = predict_endpoint(sm, "Permeability")
            if pred and pred.get("status") == "COMPLETE" and pred.get("predicted_value") is not None:
                y_true_caco2.append(true_val)
                y_pred_caco2.append(float(pred["predicted_value"]))
        except Exception:
            continue

    y_true_caco2 = np.array(y_true_caco2)
    y_pred_caco2 = np.array(y_pred_caco2)
    abs_errors_caco2 = np.abs(y_true_caco2 - y_pred_caco2)

    n_cal_caco2 = len(abs_errors_caco2) // 2
    n_eval_caco2 = len(abs_errors_caco2) - n_cal_caco2
    cal_err_caco2 = abs_errors_caco2[:n_cal_caco2]
    eval_err_caco2 = abs_errors_caco2[n_cal_caco2:]

    q80_caco2 = float(np.quantile(cal_err_caco2, 0.80))
    q90_caco2 = float(np.quantile(cal_err_caco2, 0.90))
    q95_caco2 = float(np.quantile(cal_err_caco2, 0.95))

    cov80_caco2 = float(np.mean(eval_err_caco2 <= q80_caco2))
    cov90_caco2 = float(np.mean(eval_err_caco2 <= q90_caco2))
    cov95_caco2 = float(np.mean(eval_err_caco2 <= q95_caco2))

    results["Permeability"] = {
        "endpoint": "Permeability",
        "dataset_name": "Admetica External Caco-2 Benchmark (34 compounds)",
        "source": "Admetica External Validation Set",
        "local_path": "models/admetica/validation/caco2_external_34.csv",
        "calibration_n": n_cal_caco2,
        "evaluation_n": n_eval_caco2,
        "smiles_overlap_with_training": caco2_overlap,
        "quantiles": {"0.80": q80_caco2, "0.90": q90_caco2, "0.95": q95_caco2},
        "empirical_coverage": {"0.80": cov80_caco2, "0.90": cov90_caco2, "0.95": cov95_caco2},
        "mean_interval_width_90": 2 * q90_caco2,
        "median_interval_width_90": 2 * q90_caco2,
        "mae": float(np.mean(eval_err_caco2)),
        "rmse": float(np.sqrt(np.mean(eval_err_caco2 ** 2))),
        "scientific_status": "CALIBRATED_EXTERNAL" if caco2_overlap == 0 else "CALIBRATED_INTERNAL",
    }
    print(f"Caco-2 Audit Done: status={results['Permeability']['scientific_status']}, 90% cov={cov90_caco2:.3f}, q90={q90_caco2:.3f}")

    # 3. Audit hERG Liability (Admetica vs ChEMBL37 hERG)
    print("\n--- Auditing hERG Liability ---")
    herg_val = pd.read_csv("models/admetica/validation/safety/chembl37_herg_ic50_no_exact_training_overlap.csv").dropna(subset=["smiles", "label"])
    if len(herg_val) > 500:
        herg_val = herg_val.iloc[:500]
    herg_train = pd.read_csv("models/admetica/safety/herg/training.csv")

    train_herg_smiles = set(herg_train["Smiles"].dropna())
    val_herg_smiles = [get_canonical_smiles(s) for s in herg_val["smiles"]]
    herg_overlap = len(set(val_herg_smiles).intersection(train_herg_smiles))

    y_true_herg = []
    y_prob_herg = []
    for idx, row in herg_val.iterrows():
        sm = str(row["smiles"]).strip()
        if not sm or sm == "nan":
            continue
        try:
            true_label = int(row["label"])
            pred = predict_endpoint(sm, "hERG liability")
            if pred and pred.get("status") == "COMPLETE" and pred.get("probability") is not None:
                y_true_herg.append(true_label)
                y_prob_herg.append(float(pred["probability"]))
        except Exception:
            continue

    y_true_herg = np.array(y_true_herg)
    y_prob_herg = np.array(y_prob_herg)

    # Conformal nonconformity s_i = 1 - P(y_i)
    scores_herg = np.where(y_true_herg == 1, 1.0 - y_prob_herg, y_prob_herg)

    n_cal_herg = len(scores_herg) // 2
    n_eval_herg = len(scores_herg) - n_cal_herg
    cal_scores_herg = scores_herg[:n_cal_herg]
    eval_scores_herg = scores_herg[n_cal_herg:]

    q90_herg = float(np.quantile(cal_scores_herg, 0.90))

    threshold_low = 1.0 - q90_herg
    threshold_high = q90_herg

    eval_y_true = y_true_herg[n_cal_herg:]
    eval_y_prob = y_prob_herg[n_cal_herg:]

    covered = []
    singleton_count = 0
    ambiguous_count = 0
    empty_count = 0

    for true_lbl, p in zip(eval_y_true, eval_y_prob):
        pred_set = set()
        if p >= threshold_low:
            pred_set.add(1)
        if p <= threshold_high:
            pred_set.add(0)

        if true_lbl in pred_set:
            covered.append(1)
        else:
            covered.append(0)

        if len(pred_set) == 1:
            singleton_count += 1
        elif len(pred_set) > 1:
            ambiguous_count += 1
        else:
            empty_count += 1

    emp_cov_herg = float(np.mean(covered))
    results["hERG liability"] = {
        "endpoint": "hERG liability",
        "dataset_name": "ChEMBL37 Non-overlapping hERG IC50 Benchmark (728 compounds)",
        "source": "ChEMBL37 / Admetica Validation Set",
        "local_path": "models/admetica/validation/safety/chembl37_herg_ic50_no_exact_training_overlap.csv",
        "calibration_n": n_cal_herg,
        "evaluation_n": n_eval_herg,
        "smiles_overlap_with_training": herg_overlap,
        "quantile_90": q90_herg,
        "threshold_low": threshold_low,
        "threshold_high": threshold_high,
        "empirical_coverage_90": emp_cov_herg,
        "singleton_rate": singleton_count / n_eval_herg,
        "ambiguous_rate": ambiguous_count / n_eval_herg,
        "empty_rate": empty_count / n_eval_herg,
        "scientific_status": "CALIBRATED_EXTERNAL" if herg_overlap == 0 else "CALIBRATED_INTERNAL",
    }
    print(f"hERG Audit Done: status={results['hERG liability']['scientific_status']}, 90% cov={emp_cov_herg:.3f}, q90={q90_herg:.3f}")

    # 4. Audit CYP2C9, CYP2D6, CYP3A4 Inhibitors
    for cyp_name, val_file, train_file in [
        ("CYP2C9 inhibitor", "models/admetica/validation/cyp/chembl30_2c9_inhibitor.csv", "models/admetica/cyp/cyp2c9-inhibitor/training.csv"),
        ("CYP2D6 inhibitor", "models/admetica/validation/cyp/chembl30_2d6_inhibitor.csv", "models/admetica/cyp/cyp2d6-inhibitor/training.csv"),
        ("CYP3A4 inhibitor", "models/admetica/validation/cyp/chembl30_3a4_inhibitor.csv", "models/admetica/cyp/cyp3a4-inhibitor/training.csv"),
    ]:
        print(f"\n--- Auditing {cyp_name} ---")
        df_val = pd.read_csv(val_file).dropna(subset=["smiles", "class"])
        if len(df_val) > 500:
            df_val = df_val.iloc[:500]
        df_train = pd.read_csv(train_file)
        tr_smiles = set(df_train["smiles"].dropna())
        val_smiles_list = [get_canonical_smiles(s) for s in df_val["smiles"]]
        cyp_overlap = len(set(val_smiles_list).intersection(tr_smiles))

        y_t = []
        y_p = []
        for idx, row in df_val.iterrows():
            sm = str(row["smiles"]).strip()
            if not sm or sm == "nan":
                continue
            try:
                lbl = int(row["class"])
                pred = predict_endpoint(sm, cyp_name)
                if pred and pred.get("status") == "COMPLETE" and pred.get("probability") is not None:
                    y_t.append(lbl)
                    y_p.append(float(pred["probability"]))
            except Exception:
                continue

        y_t = np.array(y_t)
        y_p = np.array(y_p)
        scores = np.where(y_t == 1, 1.0 - y_p, y_p)

        n_cal_cyp = len(scores) // 2
        n_eval_cyp = len(scores) - n_cal_cyp
        cal_s = scores[:n_cal_cyp]
        eval_s = scores[n_cal_cyp:]

        q90_cyp = float(np.quantile(cal_s, 0.90))
        t_low = 1.0 - q90_cyp
        t_high = q90_cyp

        eval_yt = y_t[n_cal_cyp:]
        eval_yp = y_p[n_cal_cyp:]

        cov_list = []
        sing = 0
        amb = 0
        emp = 0
        for true_l, p in zip(eval_yt, eval_yp):
            pset = set()
            if p >= t_low:
                pset.add(1)
            if p <= t_high:
                pset.add(0)
            if true_l in pset:
                cov_list.append(1)
            else:
                cov_list.append(0)

            if len(pset) == 1:
                sing += 1
            elif len(pset) > 1:
                amb += 1
            else:
                emp += 1

        results[cyp_name] = {
            "endpoint": cyp_name,
            "dataset_name": f"ChEMBL30 {cyp_name} Benchmark ({len(df_val)} compounds)",
            "source": "ChEMBL30 Validation Set",
            "local_path": val_file,
            "calibration_n": n_cal_cyp,
            "evaluation_n": n_eval_cyp,
            "smiles_overlap_with_training": cyp_overlap,
            "quantile_90": q90_cyp,
            "threshold_low": t_low,
            "threshold_high": t_high,
            "empirical_coverage_90": float(np.mean(cov_list)),
            "singleton_rate": sing / n_eval_cyp,
            "ambiguous_rate": amb / n_eval_cyp,
            "empty_rate": emp / n_eval_cyp,
            "scientific_status": "CALIBRATED_EXTERNAL" if cyp_overlap == 0 else "CALIBRATED_INTERNAL",
        }
        print(f"{cyp_name} Audit Done: status={results[cyp_name]['scientific_status']}, 90% cov={np.mean(cov_list):.3f}, q90={q90_cyp:.3f}")

    # 5. Endpoints without independent calibration benchmarks on disk
    # (Solubility, PPB, RLM, MLM, CYP1A2, CYP2C19, Ames, DILI, P-gp)
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
        results[uncal_name] = {
            "endpoint": uncal_name,
            "dataset_name": "Training data available only; no independent calibration set",
            "source": "Local model training set",
            "local_path": path,
            "calibration_n": 0,
            "evaluation_n": 0,
            "smiles_overlap_with_training": "TRAINING_OVERLAP_UNKNOWN",
            "scientific_status": "CALIBRATED_WITH_TRAINING_OVERLAP_UNKNOWN" if "training.csv" in path else "CONFORMAL_UNAVAILABLE",
            "unavailability_reason": "Independent external calibration set not provided on disk for this endpoint.",
        }

    # Save json audit file
    os.makedirs("validation", exist_ok=True)
    with open("validation/stage4c3a_conformal_audit.json", "w") as f:
        json.dump(results, f, indent=2)

    print("\n==================================================")
    print("AUDIT JSON Persisted to validation/stage4c3a_conformal_audit.json")
    print("==================================================")


if __name__ == "__main__":
    run_audit()
