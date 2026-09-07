"""Leakage-safe PPB/fu feature-family experiment under the frozen V72 contract.

This is research only: it writes a separate result artifact and never changes
production routing or the canonical V72 baseline.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.evaluate_human_total_iv_cl_v56 import descriptor_vector  # noqa: E402

FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)


def fp_vector(mol):
    arr = np.zeros(1024, dtype=np.float32)
    DataStructs.ConvertToNumpyArray(FP.GetFingerprint(mol), arr)
    return arr


def metrics(y, pred):
    fold = np.maximum(pred / y, y / pred)
    return {"n": int(len(y)), "aafe": float(np.mean(fold)), "median_fold_error": float(np.median(fold)),
            "within_2_fold_pct": float(np.mean(fold <= 2) * 100), "within_3_fold_pct": float(np.mean(fold <= 3) * 100),
            "mae_log10": float(np.mean(np.abs(np.log10(pred) - np.log10(y))))}


def main():
    library = json.loads((ROOT / "validation/reference_library_v1_1000.json").read_text())["compounds"]
    rows = []
    for item in library:
        fu = (item.get("evidence") or {}).get("HUMAN_FUP") if isinstance(item.get("evidence"), dict) else None
        mol = Chem.MolFromSmiles(item.get("smiles", ""))
        if mol is not None and fu is not None and 0 < float(fu) <= 1:
            rows.append((item["inchikey"], mol, float(fu)))
    rows.sort(key=lambda row: row[0])
    y = np.array([r[2] for r in rows]); logy = np.log10(y)
    xd = np.array([descriptor_vector(r[1]) for r in rows], dtype=float)
    xf = np.array([fp_vector(r[1]) for r in rows], dtype=float)
    xh = np.hstack([xd, xf])
    groups = np.array([MurckoScaffold.MurckoScaffoldSmiles(mol=r[1]) or r[0] for r in rows])
    split = list(GroupKFold(n_splits=5).split(xd, logy, groups))
    factories = {
        "extra_trees_descriptor_fp": lambda: ExtraTreesRegressor(n_estimators=400, min_samples_leaf=2, max_features=0.5, random_state=74, n_jobs=-1),
        "random_forest_descriptor_fp": lambda: RandomForestRegressor(n_estimators=400, min_samples_leaf=2, max_features=0.5, random_state=74, n_jobs=-1),
        "hist_gradient_descriptor": lambda: HistGradientBoostingRegressor(max_iter=250, l2_regularization=1.0, random_state=74),
    }
    features = {"extra_trees_descriptor_fp": xh, "random_forest_descriptor_fp": xh, "hist_gradient_descriptor": xd}
    results = {}
    for name, factory in factories.items():
        oof = np.empty(len(rows))
        for train, test in split:
            model = factory(); model.fit(features[name][train], logy[train]); oof[test] = model.predict(features[name][test])
        pred = np.clip(10 ** oof, 1e-4, 1.0)
        results[name] = {"pooled_oof": metrics(y, pred), "subgroups": {
            "fu_lt_0_01": metrics(y[y < .01], pred[y < .01]),
            "fu_0_01_to_0_1": metrics(y[(y >= .01) & (y < .1)], pred[(y >= .01) & (y < .1)]),
            "fu_ge_0_1": metrics(y[y >= .1], pred[y >= .1]),
        }}
    baseline = json.loads((ROOT / "validation/ppb_canonical_benchmark_v72.json").read_text())["candidates"]["extra_trees_logfu"]["pooled_oof"]
    best = min(results, key=lambda k: results[k]["pooled_oof"]["aafe"])
    decision = "RETAIN_RESEARCH_CANDIDATE" if results[best]["pooled_oof"]["aafe"] < baseline["aafe"] * .95 else "REJECT_NO_MEANINGFUL_IMPROVEMENT"
    output = {"artifact": "PPB_FEATURE_HYBRID_V74", "created_at": datetime.now(timezone.utc).isoformat(), "canonical_benchmark": "PPB_CANONICAL_BENCHMARK_V72", "n": len(rows), "baseline": baseline, "candidates": results, "best": best, "decision": decision, "production_engine": "drugopt-prediction-engine-v3@3.3.2"}
    (ROOT / "validation/ppb_feature_hybrid_v74.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"best": best, "decision": decision, "metrics": results[best]["pooled_oof"]}, indent=2))


if __name__ == "__main__":
    main()
