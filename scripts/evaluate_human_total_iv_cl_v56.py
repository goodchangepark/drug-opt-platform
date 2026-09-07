"""Leakage-safe development and locked evaluation for HUMAN_TOTAL_IV_CL.

PKSmart is used only for the total-IV target.  The locked clinical N=30 is
never used for fitting, feature selection, AD thresholds, or model choice.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import ExtraTreesRegressor, GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
SOURCE = Path("/tmp/pksmart-human-pk.csv")
OUT = ROOT / "validation" / "human_total_iv_cl_v5_6.json"
FINGERPRINT = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)


def canonical(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles or "")
    return Chem.MolToSmiles(mol, canonical=True) if mol else None


def descriptor_vector(mol: Chem.Mol) -> list[float]:
    return [
        Descriptors.MolWt(mol), Crippen.MolLogP(mol), rdMolDescriptors.CalcTPSA(mol),
        Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol), Lipinski.NumRotatableBonds(mol),
        Lipinski.RingCount(mol), rdMolDescriptors.CalcFractionCSP3(mol),
        Descriptors.HeavyAtomCount(mol), Chem.GetFormalCharge(mol),
    ]


def fp_vector(mol: Chem.Mol) -> np.ndarray:
    fp = FINGERPRINT.GetFingerprint(mol)
    arr = np.zeros((1024,), dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def metric(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    folds = np.maximum(y_pred / y_true, y_true / y_pred)
    log_bias = float(np.mean(np.log10(y_pred) - np.log10(y_true)))
    return {
        "n": int(len(y_true)),
        "afe": float(np.mean(y_pred / y_true)),
        "aafe": float(np.mean(folds)),
        "median_fold_error": float(np.median(folds)),
        "within_1_5_fold_pct": float(np.mean(folds <= 1.5) * 100),
        "within_2_fold_pct": float(np.mean(folds <= 2.0) * 100),
        "within_3_fold_pct": float(np.mean(folds <= 3.0) * 100),
        "bias_log10": log_bias,
        "mae_log10": float(mean_absolute_error(np.log10(y_true), np.log10(y_pred))),
    }


def rows_from_source() -> tuple[list[dict], dict]:
    if not SOURCE.exists():
        raise SystemExit(f"Missing {SOURCE}; download the public source to that path first.")
    grouped: dict[str, list[float]] = defaultdict(list)
    mols: dict[str, Chem.Mol] = {}
    raw = 0
    skipped = 0
    with SOURCE.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            raw += 1
            key = canonical(row.get("smiles_r", ""))
            try:
                value = float(row.get("human_CL_mL_min_kg", ""))
            except (TypeError, ValueError):
                skipped += 1
                continue
            mol = Chem.MolFromSmiles(row.get("smiles_r", ""))
            if not key or not mol or not math.isfinite(value) or value <= 0:
                skipped += 1
                continue
            grouped[key].append(value)
            mols[key] = mol
    records = []
    for key, values in sorted(grouped.items()):
        # Compatible repeated source records are represented by a geometric
        # median in log space; spread is retained as QC metadata.
        logs = np.log10(values)
        records.append({
            "canonical_smiles": key,
            "mol": mols[key],
            "y": float(10 ** np.median(logs)),
            "n_measurements": len(values),
            "log_spread": float(np.ptp(logs)) if len(logs) > 1 else 0.0,
        })
    return records, {"raw_rows": raw, "skipped_rows": skipped, "unique_structures": len(records), "duplicate_groups": sum(len(v) > 1 for v in grouped.values())}


def scaffold_groups(records: list[dict]) -> np.ndarray:
    return np.array([
        MurckoScaffold.MurckoScaffoldSmiles(mol=r["mol"], includeChirality=False) or r["canonical_smiles"]
        for r in records
    ])


def build_matrix(records: list[dict], mode: str) -> np.ndarray:
    if mode == "descriptors":
        return np.asarray([descriptor_vector(r["mol"]) for r in records], dtype=float)
    if mode == "fingerprints":
        return np.asarray([fp_vector(r["mol"]) for r in records], dtype=float)
    if mode == "combined":
        return np.hstack((build_matrix(records, "descriptors"), build_matrix(records, "fingerprints")))
    raise ValueError(mode)


def models(seed: int = 56) -> dict:
    return {
        "ridge_descriptors": make_pipeline(StandardScaler(), Ridge(alpha=10.0)),
        "elasticnet_descriptors": make_pipeline(StandardScaler(), ElasticNet(alpha=0.02, l1_ratio=0.2, max_iter=10000, random_state=seed)),
        "random_forest_descriptors": RandomForestRegressor(n_estimators=300, min_samples_leaf=3, max_features=0.8, random_state=seed, n_jobs=-1),
        "extra_trees_descriptors": ExtraTreesRegressor(n_estimators=300, min_samples_leaf=2, max_features=0.8, random_state=seed, n_jobs=-1),
        "gradient_boosting_descriptors": GradientBoostingRegressor(n_estimators=150, max_depth=2, learning_rate=0.03, loss="huber", random_state=seed),
        "ridge_fingerprints": make_pipeline(StandardScaler(with_mean=False), Ridge(alpha=25.0)),
        "ridge_combined": make_pipeline(StandardScaler(with_mean=False), Ridge(alpha=50.0)),
    }


def main() -> None:
    records, source_qc = rows_from_source()
    from backend.clinical_pk_cohort import CLINICAL_PK_COHORT
    n30_keys = {canonical(x.smiles) for x in CLINICAL_PK_COHORT}
    locked_mols = [Chem.MolFromSmiles(x.smiles) for x in CLINICAL_PK_COHORT if x.route == "IV" and x.cl_systemic_ml_min_kg is not None]
    locked_fps = [FINGERPRINT.GetFingerprint(m) for m in locked_mols if m]
    # A stereochemical/salt serialization can differ while the 2-D molecular
    # graph is identical.  Exclude exact Morgan duplicates as an additional
    # anti-leakage guard, while retaining the canonical-SMILES exclusion.
    near_duplicate_keys = {
        r["canonical_smiles"] for r in records
        if any(DataStructs.TanimotoSimilarity(FINGERPRINT.GetFingerprint(r["mol"]), fp) >= 0.999999 for fp in locked_fps)
    }
    development = [r for r in records if r["canonical_smiles"] not in n30_keys and r["canonical_smiles"] not in near_duplicate_keys]
    locked = [
        {"name": x.compound_name, "canonical_smiles": canonical(x.smiles), "mol": Chem.MolFromSmiles(x.smiles), "y": float(x.cl_systemic_ml_min_kg)}
        for x in CLINICAL_PK_COHORT if x.route == "IV" and x.cl_systemic_ml_min_kg is not None
    ]
    xdev = {m: build_matrix(development, m) for m in ("descriptors", "fingerprints", "combined")}
    ydev = np.log10(np.asarray([r["y"] for r in development]))
    groups = scaffold_groups(development)
    # GroupKFold is used only for development/model selection. It never sees N30.
    n_splits = min(5, len(np.unique(groups)))
    cv_results = {}
    for name, model in models().items():
        mode = "fingerprints" if "fingerprints" in name else ("combined" if "combined" in name else "descriptors")
        fold_metrics = []
        for train_idx, test_idx in GroupKFold(n_splits=n_splits).split(xdev[mode], ydev, groups):
            model.fit(xdev[mode][train_idx], ydev[train_idx])
            pred = 10 ** model.predict(xdev[mode][test_idx])
            fold_metrics.append(metric(10 ** ydev[test_idx], pred))
        cv_results[name] = {
            "mode": mode,
            "folds": fold_metrics,
            "mean_aafe": float(np.mean([m["aafe"] for m in fold_metrics])),
            "mean_within_2_fold_pct": float(np.mean([m["within_2_fold_pct"] for m in fold_metrics])),
        }
    best_name = min(cv_results, key=lambda k: cv_results[k]["mean_aafe"])
    best_mode = cv_results[best_name]["mode"]
    best = models()[best_name]
    best.fit(xdev[best_mode], ydev)
    xlocked = build_matrix(locked, best_mode)
    pred_locked = 10 ** best.predict(xlocked)
    ylocked = np.asarray([r["y"] for r in locked])
    final = metric(ylocked, pred_locked)
    # AD is descriptive only: max Morgan similarity to the development set.
    train_fps = [FINGERPRINT.GetFingerprint(r["mol"]) for r in development]
    ad_rows = []
    for r, pred in zip(locked, pred_locked):
        sims = DataStructs.BulkTanimotoSimilarity(FINGERPRINT.GetFingerprint(r["mol"]), train_fps)
        sim = max(sims) if sims else 0.0
        ad_rows.append({"compound": r["name"], "max_train_tanimoto": float(sim), "ad": "IN_DOMAIN" if sim >= 0.45 else ("BORDERLINE" if sim >= 0.30 else "OOD"), "predicted": float(pred), "observed": float(r["y"])})
    by_ad = {}
    for label in ("IN_DOMAIN", "BORDERLINE", "OOD"):
        subset = [r for r in ad_rows if r["ad"] == label]
        if subset:
            by_ad[label] = metric(np.asarray([r["observed"] for r in subset]), np.asarray([r["predicted"] for r in subset]))
    artifact = {
        "artifact": "human_total_iv_cl_v5_6",
        "target": {"id": "HUMAN_TOTAL_IV_CL", "species": "HUMAN", "route": "IV", "unit": "mL/min/kg", "transformation": "log10(CL)"},
        "source": {"name": "PKSmart Human_PK_data.csv", "path": str(SOURCE), "sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(), "qc": source_qc},
        "cohorts": {"development_n": len(development), "model_selection_n": len(development), "frozen_final_n": len(locked), "n30_overlap_excluded": len(records) - len(development), "exact_morgan_duplicate_keys_excluded": len(near_duplicate_keys)},
        "aggregation": "Geometric median in log10(CL) within canonical structure; repeated-measurement count and log spread retained.",
        "cv": {"strategy": "5-fold GroupKFold by Murcko scaffold", "results": cv_results, "selected": best_name},
        "final_locked": {"model": best_name, "features": best_mode, "metrics": final, "ad": by_ad, "rows": ad_rows},
        "comparison_guard": "This is total systemic IV clearance. It is not hepatic CL, CL/F, renal CL, or an assisted IVIVE target.",
        "production_decision": "CANDIDATE_ONLY_UNTIL_REPRODUCIBLE_REVIEW",
    }
    OUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print(json.dumps({"source": source_qc, "cohorts": artifact["cohorts"], "selected": best_name, "final": final, "ad": by_ad}, indent=2))


if __name__ == "__main__":
    main()
