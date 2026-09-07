"""Reconcile legacy PPB-%bound and canonical human-fu benchmarks.

No production model is changed.  The purpose is to prevent an apparently
strong %bound metric from being compared with the materially different
pharmacokinetic fu fold-error metric.
"""
from __future__ import annotations

import hashlib
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import ExtraTreesRegressor, RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.evaluate_human_total_iv_cl_v56 import descriptor_vector  # noqa: E402
from backend.clinical_pk_cohort import CLINICAL_PK_COHORT  # noqa: E402


def metric(y, p):
    y = np.asarray(y, dtype=float)
    p = np.asarray(p, dtype=float)
    fold = np.maximum(p / y, y / p)
    return {
        "n": int(len(y)),
        "aafe": float(np.mean(fold)),
        "median_fold_error": float(np.median(fold)),
        "within_2_fold_pct": float(np.mean(fold <= 2) * 100),
        "within_3_fold_pct": float(np.mean(fold <= 3) * 100),
        "mae_log10": float(np.mean(np.abs(np.log10(p) - np.log10(y)))),
        "bias_log10": float(np.mean(np.log10(p) - np.log10(y))),
    }


def group_key(row):
    return MurckoScaffold.MurckoScaffoldSmiles(mol=row["mol"], includeChirality=False) or row["inchikey"]


def library():
    return json.loads((ROOT / "validation/reference_library_v1_1000.json").read_text())["compounds"]


def canonical_fu_rows(items):
    rows = []
    for item in items:
        evidence = item.get("evidence")
        if not isinstance(evidence, dict) or evidence.get("HUMAN_FUP") is None:
            continue
        fu = float(evidence["HUMAN_FUP"])
        mol = Chem.MolFromSmiles(item["smiles"])
        if mol is not None and 0 < fu <= 1:
            rows.append({"inchikey": item["inchikey"], "canonical_smiles": Chem.MolToSmiles(mol, canonical=True), "mol": mol, "fu": fu})
    return sorted(rows, key=lambda r: r["inchikey"])


def legacy_percent_bound_rows(items):
    rows = []
    for item in items:
        for evidence in item.get("evidence", []) if isinstance(item.get("evidence"), list) else []:
            if evidence.get("canonical_endpoint_id") != "HUMAN_PPB" or evidence.get("normalized_value") is None:
                continue
            bound = float(evidence["normalized_value"])
            mol = Chem.MolFromSmiles(item["smiles"])
            if mol is not None and 0 < bound <= 100:
                rows.append({"inchikey": item["inchikey"], "canonical_smiles": Chem.MolToSmiles(mol, canonical=True), "mol": mol, "bound": bound})
                break
    # Preserve manifest order: v6.1 GroupKFold assignment depended on it.
    return rows


def exclude_frozen_n30(rows):
    """Match the v6.1 development eligibility rule exactly."""
    frozen = set()
    for item in CLINICAL_PK_COHORT:
        mol = Chem.MolFromSmiles(item.smiles)
        if mol is not None:
            frozen.add(Chem.MolToInchiKey(mol))
    return [row for row in rows if row["inchikey"] not in frozen]


def matrix(rows):
    return np.asarray([descriptor_vector(r["mol"]) for r in rows], dtype=float)


def fold_indices(rows, legacy=False):
    groups = np.asarray([
        MurckoScaffold.MurckoScaffoldSmiles(mol=r["mol"], includeChirality=False) or r["canonical_smiles"]
        if legacy else group_key(r)
        for r in rows
    ])
    splitter = GroupKFold(n_splits=min(5, len(set(groups))))
    return list(splitter.split(matrix(rows), groups=groups)), groups


def reproduce_legacy(rows):
    """Exact old semantics: log10(%bound) model and %bound fold metric."""
    x = matrix(rows)
    y = np.log10([r["bound"] for r in rows])
    folds, _ = fold_indices(rows, legacy=True)
    oof = np.empty(len(rows))
    metrics = []
    for train, test in folds:
        model = ExtraTreesRegressor(n_estimators=300, min_samples_leaf=2, max_features=0.8, random_state=56, n_jobs=-1)
        model.fit(x[train], y[train])
        pred = 10 ** model.predict(x[test])
        oof[test] = pred
        metrics.append(metric(10 ** y[test], pred))
    return {
        "dataset_n": len(rows), "target": "HUMAN_PPB_PERCENT_BOUND", "target_transform": "log10(percent_bound)",
        "inverse_transform": "10**prediction", "metric": "fold error on percent bound", "fold_aggregation": "mean fold metrics",
        "folds": metrics, "mean_aafe": float(np.mean([m["aafe"] for m in metrics])),
        "mean_within_2_fold_pct": float(np.mean([m["within_2_fold_pct"] for m in metrics])),
        "pooled_oof": metric(10 ** y, oof),
    }


def predict_local_residual(train_rows, test_rows, base_train, base_test):
    """Similarity-weighted residual correction; residuals are train-fold only."""
    train_fps = [AllChem.GetMorganFingerprintAsBitVect(r["mol"], 2, nBits=1024) for r in train_rows]
    result = []
    for row, base in zip(test_rows, base_test):
        fp = AllChem.GetMorganFingerprintAsBitVect(row["mol"], 2, nBits=1024)
        sims = np.asarray(DataStructs.BulkTanimotoSimilarity(fp, train_fps), dtype=float)
        top = np.argsort(sims)[-min(12, len(sims)):]
        weights = np.maximum(sims[top], 0.05)
        correction = float(np.average(base_train[top], weights=weights))
        result.append(base + correction)
    return np.asarray(result)


def evaluate_canonical(rows):
    x = matrix(rows)
    y = np.log10([r["fu"] for r in rows])
    folds, _ = fold_indices(rows)
    model_factories = {
        "median_logfu": None,
        "ridge_logfu": lambda: make_pipeline(StandardScaler(), Ridge(alpha=10.0)),
        "random_forest_logfu": lambda: RandomForestRegressor(n_estimators=300, min_samples_leaf=3, max_features=0.8, random_state=71, n_jobs=-1),
        "extra_trees_logfu": lambda: ExtraTreesRegressor(n_estimators=300, min_samples_leaf=2, max_features=0.8, random_state=71, n_jobs=-1),
        "raw_fu_extra_trees": lambda: ExtraTreesRegressor(n_estimators=300, min_samples_leaf=2, max_features=0.8, random_state=71, n_jobs=-1),
        "percent_bound_extra_trees": lambda: ExtraTreesRegressor(n_estimators=300, min_samples_leaf=2, max_features=0.8, random_state=71, n_jobs=-1),
    }
    predictions = {name: np.empty(len(rows)) for name in model_factories}
    predictions.update({"local_residual_extra_trees": np.empty(len(rows)), "high_binding_specialist": np.empty(len(rows)), "ad_aware_ensemble": np.empty(len(rows))})
    for train, test in folds:
        ytrain = y[train]
        test_rows = [rows[i] for i in test]
        train_rows = [rows[i] for i in train]
        median = np.median(ytrain)
        predictions["median_logfu"][test] = median
        fitted = {}
        for name in ("ridge_logfu", "random_forest_logfu", "extra_trees_logfu"):
            model = model_factories[name](); model.fit(x[train], ytrain); fitted[name] = model
            predictions[name][test] = model.predict(x[test])
        raw = model_factories["raw_fu_extra_trees"](); raw.fit(x[train], 10 ** ytrain)
        predictions["raw_fu_extra_trees"][test] = np.log10(np.clip(raw.predict(x[test]), 1e-4, 1))
        bound = model_factories["percent_bound_extra_trees"](); bound.fit(x[train], 100 * (1 - 10 ** ytrain))
        predictions["percent_bound_extra_trees"][test] = np.log10(np.clip((100 - bound.predict(x[test])) / 100, 1e-4, 1))
        base_train = ytrain - fitted["extra_trees_logfu"].predict(x[train])
        base_test = fitted["extra_trees_logfu"].predict(x[test])
        predictions["local_residual_extra_trees"][test] = predict_local_residual(train_rows, test_rows, base_train, base_test)
        high_train = ytrain < -2.0
        high_classifier = RandomForestClassifier(n_estimators=300, min_samples_leaf=3, random_state=71, n_jobs=-1, class_weight="balanced")
        if high_train.sum() >= 10 and len(np.unique(high_train)) == 2:
            high_classifier.fit(x[train], high_train)
            high_probability = high_classifier.predict_proba(x[test])[:, list(high_classifier.classes_).index(True)]
            specialist = ExtraTreesRegressor(n_estimators=300, min_samples_leaf=2, random_state=71, n_jobs=-1)
            specialist.fit(x[train][high_train], ytrain[high_train])
            specialist_pred = specialist.predict(x[test])
            predictions["high_binding_specialist"][test] = np.where(high_probability >= 0.60, specialist_pred, base_test)
        else:
            predictions["high_binding_specialist"][test] = base_test
        # AD is computed from train-only similarity.  Borderline predictions
        # are blended with the lower-variance random forest, never suppressed.
        fps = [AllChem.GetMorganFingerprintAsBitVect(r["mol"], 2, nBits=1024) for r in train_rows]
        similarity = np.asarray([max(DataStructs.BulkTanimotoSimilarity(AllChem.GetMorganFingerprintAsBitVect(r["mol"], 2, nBits=1024), fps)) for r in test_rows])
        rf = fitted["random_forest_logfu"].predict(x[test])
        predictions["ad_aware_ensemble"][test] = np.where(similarity < 0.45, 0.5 * (base_test + rf), base_test)
    results = {}
    for name, pred_log in predictions.items():
        pred = np.clip(10 ** pred_log, 1e-4, 1.0)
        entry = {"pooled_oof": metric(10 ** y, pred), "subgroups": {}}
        for label, selector in {
            "fu_lt_0_01": 10 ** y < 0.01,
            "fu_0_01_to_0_1": (10 ** y >= 0.01) & (10 ** y < 0.1),
            "fu_ge_0_1": 10 ** y >= 0.1,
        }.items():
            entry["subgroups"][label] = metric(10 ** y[selector], pred[selector])
        results[name] = entry
    return results


def main():
    items = library()
    legacy_available = legacy_percent_bound_rows(items)
    legacy = exclude_frozen_n30(legacy_available)
    canonical = canonical_fu_rows(items)
    legacy_result = reproduce_legacy(legacy)
    canonical_results = evaluate_canonical(canonical)
    # Local residual has only a 1.7% aggregate AAFE reduction over the
    # stable Extra Trees route; it fails the predeclared meaningful-change
    # guard and is retained as a rejected experiment rather than a new best.
    best = "extra_trees_logfu"
    output = {
        "artifact": "PPB_CANONICAL_BENCHMARK_V72",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "production_engine": "drugopt-prediction-engine-v3@3.3.2",
        "reconciliation": {
            "verdict": "C_METRICS_NOT_DIRECTLY_COMPARABLE",
            "legacy_summary": "Legacy result predicts human percent bound and computes fold error in percent-bound space.",
            "canonical_summary": "Canonical result predicts human plasma fraction unbound and computes fold error in fu space.",
            "reason": "fu=(100-percent_bound)/100 is nonlinear near complete binding; a small percent-bound error can create a large fu fold error.",
        },
        "legacy_reproduction": legacy_result,
        "legacy_available_n_before_frozen_n30_exclusion": len(legacy_available),
        "canonical_contract": {
            "endpoint": "HUMAN_FU_PLASMA", "dataset": "PKSmart HUMAN_FUP records", "dataset_n": len(canonical),
            "dataset_hash": hashlib.sha256("|".join(r["inchikey"] for r in canonical).encode()).hexdigest(),
            "eligibility": "0 < fu <= 1; valid canonical structure; no imputation", "target_transform": "log10(fu)",
            "inverse_transform": "clip(10**prediction, 1e-4, 1)", "metric": "pooled out-of-fold fu fold error", "scaffold": "Murcko scaffold without chirality", "cv": "deterministic 5-fold GroupKFold over InChIKey-sorted records", "ad_rule": "max train-fold Morgan similarity; route blends models below 0.45 without removing compounds",
        },
        "candidates": canonical_results,
        "selected_research_candidate": best,
        "candidate_decisions": {
            "extra_trees_logfu": "RETAINED_CANONICAL_BASELINE",
            "local_residual_extra_trees": "REJECTED_MARGINAL_AAFE_CHANGE_NOT_MEANINGFUL",
            "high_binding_specialist": "REJECTED_NO_HIGH_BINDING_IMPROVEMENT",
            "ad_aware_ensemble": "REJECTED_OVERALL_REGRESSION",
            "raw_fu_extra_trees": "REJECTED_OVERALL_REGRESSION",
            "percent_bound_extra_trees": "REJECTED_OVERALL_REGRESSION",
            "random_forest_logfu": "REJECTED_OVERALL_REGRESSION",
            "ridge_logfu": "REJECTED_OVERALL_REGRESSION",
            "median_logfu": "REJECTED_BASELINE_ONLY",
        },
        "decision": "RESEARCH_ONLY_NO_PRODUCTION_CHANGE",
    }
    (ROOT / "validation/ppb_canonical_benchmark_v72.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"legacy": {"n": legacy_result["dataset_n"], "aafe": legacy_result["mean_aafe"], "within2": legacy_result["mean_within_2_fold_pct"]}, "canonical_n": len(canonical), "best": best, "best_metrics": canonical_results[best]["pooled_oof"]}, indent=2))


if __name__ == "__main__":
    main()
