#!/usr/bin/env python3
"""PPB/fu high-binding weighting, calibration, and frozen-model audit."""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from sklearn.ensemble import ExtraTreesRegressor, GradientBoostingRegressor, RandomForestRegressor

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.admet_predictor import predict_batch_values
from scripts.reconcile_ppb_benchmark_v72 import canonical_fu_rows, fold_indices, library, matrix, metric
from scripts.run_ppb_local_ladder_v75 import _inner_oof, _subgroups

OUTPUT = ROOT / "validation" / "ppb_weighted_ladder_v76.json"


def weights(log_y, high_weight: float, mid_weight: float = 1.0):
    result = np.ones(len(log_y), dtype=float)
    result[log_y < -1.0] = mid_weight
    result[log_y < -2.0] = high_weight
    return result


def main() -> None:
    rows = canonical_fu_rows(library())
    x = matrix(rows)
    y = np.asarray([row["fu"] for row in rows], dtype=float)
    log_y = np.log10(y)
    folds, _ = fold_indices(rows)
    names = [
        "extra_trees_v72", "weighted_high_2", "weighted_high_4", "weighted_high_8",
        "weighted_high_mid", "random_forest_weighted", "huber_gradient",
        "predicted_regime_calibration", "fixed_base_weighted_blend",
    ]
    oof = {name: np.empty(len(rows), dtype=float) for name in names}
    for train, test in folds:
        base = ExtraTreesRegressor(n_estimators=300, min_samples_leaf=2, max_features=.8, random_state=71, n_jobs=2)
        base.fit(x[train], log_y[train])
        base_test = base.predict(x[test])
        oof["extra_trees_v72"][test] = base_test
        weighted_predictions = {}
        for factor in (2, 4, 8):
            model = ExtraTreesRegressor(n_estimators=350, min_samples_leaf=2, max_features=.8, random_state=76, n_jobs=2)
            model.fit(x[train], log_y[train], sample_weight=weights(log_y[train], factor))
            weighted_predictions[factor] = model.predict(x[test])
            oof[f"weighted_high_{factor}"][test] = weighted_predictions[factor]
        model = ExtraTreesRegressor(n_estimators=350, min_samples_leaf=2, max_features=.8, random_state=77, n_jobs=2)
        model.fit(x[train], log_y[train], sample_weight=weights(log_y[train], 5, 2))
        high_mid = model.predict(x[test])
        oof["weighted_high_mid"][test] = high_mid
        rf = RandomForestRegressor(n_estimators=350, min_samples_leaf=2, max_features=.8, random_state=76, n_jobs=2)
        rf.fit(x[train], log_y[train], sample_weight=weights(log_y[train], 5, 1.5))
        oof["random_forest_weighted"][test] = rf.predict(x[test])
        huber = GradientBoostingRegressor(loss="huber", n_estimators=220, learning_rate=.025, max_depth=2, min_samples_leaf=5, random_state=76)
        huber.fit(x[train], log_y[train], sample_weight=weights(log_y[train], 4, 1.5))
        oof["huber_gradient"][test] = huber.predict(x[test])

        inner_base = _inner_oof(x[train], log_y[train], [rows[i] for i in train])
        residual = log_y[train] - inner_base
        edges = (-np.inf, -2.0, -1.0, -0.3, np.inf)
        corrected = base_test.copy()
        for low, high in zip(edges[:-1], edges[1:]):
            train_mask = (inner_base >= low) & (inner_base < high)
            test_mask = (base_test >= low) & (base_test < high)
            if train_mask.sum() >= 8:
                correction = float(np.median(residual[train_mask]) * train_mask.sum() / (train_mask.sum() + 12))
                corrected[test_mask] += correction
        oof["predicted_regime_calibration"][test] = corrected
        oof["fixed_base_weighted_blend"][test] = .7 * base_test + .3 * high_mid

    results = {}
    for name, pred_log in oof.items():
        pred = np.clip(10 ** pred_log, 1e-4, 1.0)
        results[name] = {"pooled_oof": metric(y, pred), "subgroups": _subgroups(y, pred)}

    # This installed public checkpoint is a useful scientific comparator, but
    # it is excluded from model selection because exact source-lineage overlap
    # with PKSmart has not been ruled out.
    raw_bound = np.asarray(predict_batch_values([row["canonical_smiles"] for row in rows], "Plasma protein binding"), dtype=float)
    valid = (raw_bound >= 0) & (raw_bound <= 100)
    pretrained_fu = np.clip(1.0 - raw_bound / 100.0, 1e-4, 1.0)
    pretrained = {
        "checkpoint": "admetica-d4f7056-chemprop-v2.1",
        "training_dataset": "AstraZeneca/ChEMBL CHEMBL3301361",
        "valid_physical_outputs": int(valid.sum()),
        "metrics_on_valid_outputs": metric(y[valid], pretrained_fu[valid]),
        "subgroups_on_valid_outputs": _subgroups(y[valid], pretrained_fu[valid]),
        "selection_eligible": False,
        "reason": "Exact source/assay-lineage independence from PKSmart has not been established.",
    }
    baseline = results["extra_trees_v72"]["pooled_oof"]
    ranked = sorted(results, key=lambda name: results[name]["pooled_oof"]["aafe"])
    best = ranked[0]
    best_metrics = results[best]["pooled_oof"]
    baseline_high = results["extra_trees_v72"]["subgroups"]["fu_0_001_to_0_01"]["aafe"]
    best_high = results[best]["subgroups"]["fu_0_001_to_0_01"]["aafe"]
    meaningful = best != "extra_trees_v72" and best_metrics["aafe"] <= baseline["aafe"] * .95 and best_high <= baseline_high * 1.25
    output = {
        "artifact": "PPB_WEIGHTED_LADDER_V76",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "production_engine": "drugopt-prediction-engine-v3@3.3.3",
        "canonical_benchmark": "validation/ppb_canonical_benchmark_v72.json",
        "n": len(rows),
        "leakage_guards": {
            "outer_split": "V72 scaffold GroupKFold",
            "weighting_uses_training_targets_only": True,
            "regime_calibration_uses_inner_oof_predictions": True,
            "true_held_out_fu_routing": False,
        },
        "baseline": baseline,
        "ranked_candidates": ranked,
        "candidates": results,
        "best_candidate": best,
        "best_metrics": best_metrics,
        "pretrained_checkpoint_diagnostic": pretrained,
        "decision": "RETAIN_RESEARCH_CANDIDATE" if meaningful else "REJECT_NO_MEANINGFUL_ROBUST_GAIN",
        "production_promotion": False,
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"best": best, "best_metrics": best_metrics, "decision": output["decision"], "pretrained": pretrained["metrics_on_valid_outputs"]}, indent=2))


if __name__ == "__main__":
    main()
