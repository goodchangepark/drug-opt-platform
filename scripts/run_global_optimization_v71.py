"""Run the current-data global optimization checkpoint.

This research-only checkpoint deliberately evaluates PPB/fu alternatives on
qualified existing evidence, records rejected candidates, and updates the
completion board without changing production routing or the runtime DB.
"""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import ExtraTreesRegressor, GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.evaluate_human_total_iv_cl_v56 import descriptor_vector, FINGERPRINT  # noqa: E402


def fold_metric(y_true, y_pred):
    fold = np.maximum(y_pred / y_true, y_true / y_pred)
    return {
        "n": int(len(y_true)),
        "aafe": float(np.mean(fold)),
        "within_2_fold_pct": float(np.mean(fold <= 2) * 100),
        "within_3_fold_pct": float(np.mean(fold <= 3) * 100),
        "mae_log10": float(np.mean(np.abs(np.log10(y_pred) - np.log10(y_true)))),
        "bias_log10": float(np.mean(np.log10(y_pred) - np.log10(y_true))),
    }


def ppb_rows():
    library = json.loads((ROOT / "validation/reference_library_v1_1000.json").read_text())["compounds"]
    rows = []
    for item in library:
        evidence = item.get("evidence")
        if not isinstance(evidence, dict) or evidence.get("HUMAN_FUP") is None:
            continue
        value = float(evidence["HUMAN_FUP"])
        mol = Chem.MolFromSmiles(item["smiles"])
        if mol and 0 < value <= 1:
            rows.append({"key": item["inchikey"], "mol": mol, "fu": value, "source_family": item.get("source_family")})
    return rows


def candidate_models():
    return {
        "ridge_logfu": (make_pipeline(StandardScaler(), Ridge(alpha=10.0)), "logfu"),
        "random_forest_logfu": (RandomForestRegressor(n_estimators=300, min_samples_leaf=3, max_features=0.8, random_state=71, n_jobs=-1), "logfu"),
        "extra_trees_logfu": (ExtraTreesRegressor(n_estimators=300, min_samples_leaf=2, max_features=0.8, random_state=71, n_jobs=-1), "logfu"),
        "gradient_boosting_huber_logfu": (GradientBoostingRegressor(n_estimators=150, max_depth=2, learning_rate=0.03, loss="huber", random_state=71), "logfu"),
        "extra_trees_percent_bound": (ExtraTreesRegressor(n_estimators=300, min_samples_leaf=2, max_features=0.8, random_state=71, n_jobs=-1), "percent_bound"),
    }


def evaluate_ppb(rows):
    x = np.asarray([descriptor_vector(r["mol"]) for r in rows], dtype=float)
    y_log = np.log10([r["fu"] for r in rows])
    groups = np.asarray([MurckoScaffold.MurckoScaffoldSmiles(mol=r["mol"]) or r["key"] for r in rows])
    splitter = GroupKFold(n_splits=min(5, len(np.unique(groups))))
    results = {}
    for name, (model, target) in candidate_models().items():
        folds = []
        target_values = y_log if target == "logfu" else np.asarray([100 * (1 - r["fu"]) for r in rows])
        for train, test in splitter.split(x, target_values, groups):
            model.fit(x[train], target_values[train])
            raw = model.predict(x[test])
            pred = np.clip(10 ** raw, 1e-4, 1.0) if target == "logfu" else np.clip((100 - raw) / 100, 1e-4, 1.0)
            folds.append(fold_metric(np.asarray([rows[i]["fu"] for i in test]), pred))
        results[name] = {
            "endpoint": "HUMAN_FU_PLASMA",
            "features": "RDKit descriptors",
            "target_transform": target,
            "cv_strategy": "5-fold GroupKFold by Murcko scaffold",
            "folds": folds,
            "median_aafe": float(np.median([f["aafe"] for f in folds])),
            "mean_aafe": float(np.mean([f["aafe"] for f in folds])),
            "worst_aafe": float(np.max([f["aafe"] for f in folds])),
            "mean_within_2_fold_pct": float(np.mean([f["within_2_fold_pct"] for f in folds])),
            "mean_mae_log10": float(np.mean([f["mae_log10"] for f in folds])),
            "mean_bias_log10": float(np.mean([f["bias_log10"] for f in folds])),
            "decision": "REJECTED",
        }
    best = min(results, key=lambda k: (results[k]["mean_aafe"], -results[k]["mean_within_2_fold_pct"]))
    results[best]["decision"] = "RESEARCH_CANDIDATE_BEST_CURRENT_DATA"
    return results, best


def update_board(board, ppb_best, ledger_path):
    stable = {"SOLUBILITY", "CACO2_PAPP_AB", "HLM_CLINT", "CYP3A4", "HERG"}
    mechanistic = {"PKA_ACID", "PKA_BASE", "LOGD_7_4"}
    optimizing = {"HUMAN_FU_PLASMA", "HUMAN_HEPATIC_CL"}
    for item in board["rows"]:
        endpoint = item["endpoint"]
        if endpoint in stable:
            item["status"] = "PRODUCTION_STABLE"
        elif endpoint in mechanistic:
            item["status"] = "MECHANISTIC_ONLY"
        elif endpoint in optimizing:
            item["status"] = "MODEL_OPTIMIZING"
        else:
            item["status"] = "CURRENT_DATA_CEILING"
        item["last_experiment"] = ledger_path
        if endpoint == "HUMAN_FU_PLASMA":
            item["current_value"] = ppb_best["mean_aafe"]
            item["development_N"] = 476
            item["source_independence"] = "PKSmart qualified current-data development; independent validation required"
            item["candidate_model"] = ppb_best["name"]
            item["next_action"] = "Seek independent validation when available; continue high-binding robustness analysis without changing production."
        elif endpoint in mechanistic:
            item["next_action"] = "Improve consensus ionization logic and uncertainty using existing structures; do not claim experimental ML."
        elif item["status"] == "CURRENT_DATA_CEILING":
            item["next_action"] = "Revisit after new qualified labels or a validated upstream dependency; keep fail-closed semantics."
    board["status_policy"] = ["DONE", "PRODUCTION_STABLE", "VALIDATED_CANDIDATE", "RESEARCH_CANDIDATE", "MODEL_OPTIMIZING", "MECHANISTIC_ONLY", "CURRENT_DATA_CEILING"]
    board["control_plane_version"] = "V71"
    board["note"] = "CURRENT_DATA_CEILING is resumable, not an endpoint completion claim; external validation remains required before production promotion."


def main():
    rows = ppb_rows()
    results, best = evaluate_ppb(rows)
    ledger = {
        "artifact": "GLOBAL_OPTIMIZATION_LEDGER_V71",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "production_engine": "drugopt-prediction-engine-v3@3.3.2",
        "experiments": [
            {"experiment_id": "V71-PPB-LOGFU-LADDER", "endpoint": "HUMAN_FU_PLASMA", "training_N": len(rows), "validation_N": len(rows), "models": results, "selected": best, "reason": "Research-only scaffold CV; no independent source validation."}
        ],
        "pka": {"qualified_experimental_N": 0, "status": "MECHANISTIC_ONLY", "reason": "No qualified experimental labels in current accessible assets."},
        "logd74": {"qualified_experimental_N": 0, "status": "MECHANISTIC_ONLY", "reason": "No qualified pH-explicit experimental labels in current accessible assets."},
        "next_queue": ["HUMAN_FU_PLASMA", "HUMAN_HEPATIC_CL", "HUMAN_VDSS", "HUMAN_TOTAL_IV_CL", "IV_HALF_LIFE", "ORAL_AUC", "ORAL_CMAX"],
        "production_decision": "UNCHANGED",
    }
    ledger["dataset_hash"] = hashlib.sha256(json.dumps(rows, sort_keys=True, default=str).encode()).hexdigest()
    ledger_path = "validation/global_optimization_ledger_v71.json"
    (ROOT / ledger_path).write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")
    board_path = ROOT / "validation/endpoint_completion_board.json"
    board = json.loads(board_path.read_text())
    update_board(board, {"name": best, **results[best]}, ledger_path)
    board["last_updated"] = datetime.now(timezone.utc).isoformat()
    board_path.write_text(json.dumps(board, indent=2, sort_keys=True) + "\n")
    graph = {
        "artifact": "ENDPOINT_DEPENDENCY_GRAPH_V71",
        "edges": [
            ["PPB_FU", "HUMAN_HEPATIC_CL"], ["HLM_CLINT", "HUMAN_HEPATIC_CL"],
            ["HUMAN_TOTAL_IV_CL", "IV_HALF_LIFE"], ["HUMAN_VDSS", "IV_HALF_LIFE"],
            ["HUMAN_F", "ORAL_AUC"], ["HUMAN_TOTAL_IV_CL", "ORAL_AUC"],
            ["KA", "ORAL_CMAX"], ["HUMAN_F", "ORAL_CMAX"], ["HUMAN_VDSS", "ORAL_CMAX"],
            ["PKA_ACID", "LOGD_7_4"], ["PKA_BASE", "LOGD_7_4"], ["LOGD_7_4", "HUMAN_VDSS"],
        ],
        "routing_rule": "When a parent candidate is retained, schedule dependent research reevaluation; never silently alter production routing.",
    }
    (ROOT / "validation/endpoint_dependency_graph_v71.json").write_text(json.dumps(graph, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"ppb_n": len(rows), "best": best, "best_mean_aafe": results[best]["mean_aafe"], "best_within_2": results[best]["mean_within_2_fold_pct"], "board_status_counts": {s: sum(r["status"] == s for r in board["rows"]) for s in board["status_policy"]}}, indent=2))


if __name__ == "__main__":
    main()
