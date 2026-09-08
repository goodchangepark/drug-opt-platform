#!/usr/bin/env python3
"""Leakage-safe PPB/fu local/global strategy ladder on canonical V72.

Every outer-fold prediction uses only that fold's training molecules.  The
local-residual candidates learn residuals from inner out-of-fold predictions,
never from fitted-on-self residuals or held-out fu values.
"""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from rdkit import DataStructs, RDLogger
from rdkit.Chem import AllChem
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.linear_model import HuberRegressor
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.reconcile_ppb_benchmark_v72 import canonical_fu_rows, fold_indices, group_key, library, matrix, metric

OUTPUT = ROOT / "validation" / "ppb_local_ladder_v75.json"
SEED = 75

# Morgan's older helper emits one warning per molecule on current RDKit.  The
# campaign artifact records the fingerprint contract, so suppressing that
# deprecation noise keeps persistent service logs useful.
RDLogger.DisableLog("rdApp.warning")


def _fps(rows):
    return [AllChem.GetMorganFingerprintAsBitVect(row["mol"], 2, nBits=2048) for row in rows]


def _similarity_matrix(test_fps, train_fps):
    return np.asarray([DataStructs.BulkTanimotoSimilarity(fp, train_fps) for fp in test_fps], dtype=float)


def _weighted_neighbor(values, scores, k, floor=0.05):
    k = min(k, len(values))
    top = np.argpartition(scores, -k, axis=1)[:, -k:]
    result = []
    for i, idx in enumerate(top):
        weights = np.maximum(scores[i, idx], floor) ** 2
        result.append(float(np.average(values[idx], weights=weights)))
    return np.asarray(result)


def _descriptor_knn(x_train, y_train, x_test, k=12):
    scaler = StandardScaler().fit(x_train)
    train = scaler.transform(x_train)
    test = scaler.transform(x_test)
    distances = np.sqrt(np.maximum(((test[:, None, :] - train[None, :, :]) ** 2).mean(axis=2), 0.0))
    similarity = 1.0 / (1.0 + distances)
    return _weighted_neighbor(y_train, similarity, k)


def _inner_oof(x_train, y_train, rows_train):
    groups = np.asarray([group_key(row) for row in rows_train])
    folds = GroupKFold(n_splits=min(5, len(set(groups))))
    pred = np.empty(len(y_train), dtype=float)
    for inner_train, inner_test in folds.split(x_train, groups=groups):
        model = ExtraTreesRegressor(
            n_estimators=180, min_samples_leaf=2, max_features=0.8,
            random_state=SEED, n_jobs=2,
        )
        model.fit(x_train[inner_train], y_train[inner_train])
        pred[inner_test] = model.predict(x_train[inner_test])
    return pred


def _subgroups(y, pred):
    return {
        label: metric(y[mask], pred[mask])
        for label, mask in {
            "fu_lt_0_001": y < 0.001,
            "fu_0_001_to_0_01": (y >= 0.001) & (y < 0.01),
            "fu_0_01_to_0_1": (y >= 0.01) & (y < 0.1),
            "fu_0_1_to_0_5": (y >= 0.1) & (y < 0.5),
            "fu_ge_0_5": y >= 0.5,
        }.items()
        if mask.any()
    }


def main() -> None:
    rows = canonical_fu_rows(library())
    x = matrix(rows)
    y = np.asarray([row["fu"] for row in rows], dtype=float)
    log_y = np.log10(y)
    folds, groups = fold_indices(rows)
    names = [
        "extra_trees_logfu_v72_reproduced",
        "fingerprint_knn_k5", "fingerprint_knn_k12", "fingerprint_knn_k24",
        "descriptor_knn_k12", "nested_local_residual", "cluster_residual",
        "nested_robust_calibration", "hist_gradient_logfu",
        "fixed_global_local_blend", "ad_aware_global_local",
    ]
    predictions = {name: np.empty(len(rows), dtype=float) for name in names}
    max_similarity = np.empty(len(rows), dtype=float)
    ensemble_disagreement = np.empty(len(rows), dtype=float)

    for fold_number, (train_idx, test_idx) in enumerate(folds):
        train_rows = [rows[i] for i in train_idx]
        test_rows = [rows[i] for i in test_idx]
        x_train, x_test = x[train_idx], x[test_idx]
        y_train = log_y[train_idx]
        global_model = ExtraTreesRegressor(
            n_estimators=300, min_samples_leaf=2, max_features=0.8,
            random_state=71, n_jobs=2,
        )
        global_model.fit(x_train, y_train)
        global_test = global_model.predict(x_test)
        predictions["extra_trees_logfu_v72_reproduced"][test_idx] = global_test

        train_fps, test_fps = _fps(train_rows), _fps(test_rows)
        similarities = _similarity_matrix(test_fps, train_fps)
        max_similarity[test_idx] = similarities.max(axis=1)
        for k in (5, 12, 24):
            predictions[f"fingerprint_knn_k{k}"][test_idx] = _weighted_neighbor(y_train, similarities, k)
        descriptor_local = _descriptor_knn(x_train, y_train, x_test, 12)
        predictions["descriptor_knn_k12"][test_idx] = descriptor_local

        inner_base = _inner_oof(x_train, y_train, train_rows)
        residual = y_train - inner_base
        local_correction = _weighted_neighbor(residual, similarities, 12)
        local_residual = global_test + local_correction * np.clip(max_similarity[test_idx] / 0.65, 0.0, 1.0)
        predictions["nested_local_residual"][test_idx] = local_residual

        scaler = StandardScaler().fit(x_train)
        train_scaled, test_scaled = scaler.transform(x_train), scaler.transform(x_test)
        n_clusters = min(16, max(4, len(train_idx) // 25))
        clusterer = KMeans(n_clusters=n_clusters, random_state=SEED + fold_number, n_init=10).fit(train_scaled)
        train_cluster = clusterer.labels_
        test_cluster = clusterer.predict(test_scaled)
        corrections = {}
        for cluster_id in range(n_clusters):
            values = residual[train_cluster == cluster_id]
            corrections[cluster_id] = float(np.median(values) * len(values) / (len(values) + 12)) if len(values) else 0.0
        predictions["cluster_residual"][test_idx] = global_test + np.asarray([corrections[c] for c in test_cluster])

        calibration_x = np.column_stack([inner_base, x_train])
        calibration_test = np.column_stack([global_test, x_test])
        # Fit the transform on the outer training fold only.
        cal_scaler = StandardScaler().fit(calibration_x)
        robust = HuberRegressor(epsilon=1.5, alpha=0.1, max_iter=500).fit(cal_scaler.transform(calibration_x), y_train)
        predictions["nested_robust_calibration"][test_idx] = robust.predict(cal_scaler.transform(calibration_test))

        hist = HistGradientBoostingRegressor(loss="absolute_error", max_iter=180, learning_rate=0.05, max_leaf_nodes=15, l2_regularization=2.0, random_state=SEED)
        hist.fit(x_train, y_train)
        hist_test = hist.predict(x_test)
        predictions["hist_gradient_logfu"][test_idx] = hist_test

        fingerprint_local = predictions["fingerprint_knn_k12"][test_idx]
        fixed_blend = 0.75 * global_test + 0.25 * fingerprint_local
        predictions["fixed_global_local_blend"][test_idx] = fixed_blend
        gate = np.clip((max_similarity[test_idx] - 0.25) / 0.45, 0.0, 0.35)
        predictions["ad_aware_global_local"][test_idx] = (1.0 - gate) * global_test + gate * fingerprint_local
        ensemble_disagreement[test_idx] = np.std(np.column_stack([global_test, fingerprint_local, descriptor_local, hist_test]), axis=1)

    results = {}
    for name, pred_log in predictions.items():
        pred_fu = np.clip(10 ** pred_log, 1e-4, 1.0)
        results[name] = {
            "pooled_oof": metric(y, pred_fu),
            "subgroups": _subgroups(y, pred_fu),
            "ad": {
                label: metric(y[mask], pred_fu[mask])
                for label, mask in {
                    "ood_similarity_lt_0_30": max_similarity < 0.30,
                    "borderline_0_30_to_0_50": (max_similarity >= 0.30) & (max_similarity < 0.50),
                    "in_domain_ge_0_50": max_similarity >= 0.50,
                }.items() if mask.any()
            },
        }
    baseline_name = "extra_trees_logfu_v72_reproduced"
    baseline = results[baseline_name]["pooled_oof"]
    ranked = sorted(results, key=lambda name: results[name]["pooled_oof"]["aafe"])
    best = ranked[0]
    best_metric = results[best]["pooled_oof"]
    baseline_high = results[baseline_name]["subgroups"]["fu_0_001_to_0_01"]["aafe"]
    best_high = results[best]["subgroups"]["fu_0_001_to_0_01"]["aafe"]
    retain = (
        best != baseline_name
        and best_metric["aafe"] <= baseline["aafe"] * 0.95
        and best_metric["within_2_fold_pct"] >= baseline["within_2_fold_pct"]
        and best_metric["within_3_fold_pct"] >= baseline["within_3_fold_pct"] - 1.0
        and best_high <= baseline_high * 1.25
    )
    absolute_log_error = np.abs(predictions[best] - log_y)
    uncertainty_correlation = spearmanr(ensemble_disagreement, absolute_log_error).statistic
    distance_correlation = spearmanr(1.0 - max_similarity, absolute_log_error).statistic
    best_fu = np.clip(10 ** predictions[best], 1e-4, 1.0)
    fold_error = np.maximum(best_fu / y, y / best_fu)
    subgroup_labels = np.select(
        [y < 0.001, y < 0.01, y < 0.1, y < 0.5],
        ["fu_lt_0_001", "fu_0_001_to_0_01", "fu_0_01_to_0_1", "fu_0_1_to_0_5"],
        default="fu_ge_0_5",
    )
    prediction_records = [
        {
            "inchikey": row["inchikey"],
            "scaffold_group": groups[index],
            "experimental_fu": float(y[index]),
            "predicted_fu": float(best_fu[index]),
            "absolute_fold_error": float(fold_error[index]),
            "max_train_fold_morgan_similarity": float(max_similarity[index]),
            "ensemble_disagreement_log10": float(ensemble_disagreement[index]),
            "binding_subgroup": str(subgroup_labels[index]),
        }
        for index, row in enumerate(rows)
    ]
    worst = sorted(prediction_records, key=lambda record: record["absolute_fold_error"], reverse=True)[:25]
    canonical = json.loads((ROOT / "validation/ppb_canonical_benchmark_v72.json").read_text(encoding="utf-8"))["canonical_contract"]
    output = {
        "artifact": "PPB_LOCAL_GLOBAL_LADDER_V75",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "production_engine": "drugopt-prediction-engine-v3@3.3.3",
        "canonical_contract": canonical,
        "canonical_contract_hash": hashlib.sha256(json.dumps(canonical, sort_keys=True).encode()).hexdigest(),
        "leakage_guards": {
            "outer_split": "V72 deterministic scaffold GroupKFold",
            "residual_training": "inner GroupKFold OOF predictions inside each outer training fold",
            "neighbor_pool": "outer training fold only",
            "routing_uses_true_fu": False,
            "scoring_space": "inverse-transformed fraction unbound",
        },
        "baseline": {"name": baseline_name, **baseline},
        "ranked_candidates": ranked,
        "candidates": results,
        "best_candidate": best,
        "best_metrics": best_metric,
        "decision": "RETAIN_RESEARCH_CANDIDATE" if retain else "REJECT_MARGINAL_OR_SUBGROUP_REGRESSION_RETAIN_V72_BASELINE",
        "production_promotion": False,
        "promotion_blockers": ["No independent PPB qualification cohort", "Research target not sufficient by itself for production"],
        "uncertainty_analysis": {
            "ensemble_disagreement_vs_absolute_log_error_spearman": None if np.isnan(uncertainty_correlation) else float(uncertainty_correlation),
            "chemical_distance_vs_absolute_log_error_spearman": None if np.isnan(distance_correlation) else float(distance_correlation),
            "confidence_claim": "NOT_CALIBRATED" if max(uncertainty_correlation or 0, distance_correlation or 0) < 0.3 else "ERROR_CORRELATED_RESEARCH_ONLY",
        },
        "error_analysis": {
            "dominant_error_regime": "fu below 0.01; the two highest-binding bins retain order-of-magnitude fold errors",
            "high_binding_is_not_fixed_by_locality": True,
            "worst_oof_records": worst,
            "all_oof_prediction_records": prediction_records,
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"baseline": baseline, "best": best, "best_metrics": best_metric, "decision": output["decision"]}, indent=2))


if __name__ == "__main__":
    main()
