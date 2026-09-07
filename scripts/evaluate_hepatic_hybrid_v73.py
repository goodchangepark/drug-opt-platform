"""Frozen, leakage-safe comparison of hepatic-IVIVE residual hybrids.

The clinical cohort is an *assisted* benchmark: measured fu and HLM Clint are
inputs to the mechanistic prior, while measured human hepatic CL is the target.
Each held-out scaffold fold receives a residual model fitted only on the other
folds.  This is research evaluation, never a production-routing change.
"""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from backend.clinical_pk_cohort import CLINICAL_PK_COHORT  # noqa: E402
from backend.pk_parameter_set import calculate_well_stirred_clearance, canonical_fu_from_ppb  # noqa: E402
from scripts.evaluate_human_total_iv_cl_v56 import descriptor_vector  # noqa: E402

OUT = ROOT / "validation/hepatic_hybrid_benchmark_v73.json"
LEDGER_OUT = ROOT / "validation/global_optimization_ledger_v73.json"
BOARD_OUT = ROOT / "validation/endpoint_completion_board.json"


def metrics(y, pred):
    fold = np.maximum(pred / y, y / pred)
    return {
        "n": int(len(y)), "aafe": float(np.mean(fold)),
        "median_fold_error": float(np.median(fold)),
        "within_2_fold_pct": float(np.mean(fold <= 2) * 100),
        "within_3_fold_pct": float(np.mean(fold <= 3) * 100),
        "mae_log10": float(np.mean(np.abs(np.log10(pred) - np.log10(y)))),
        "bias_log10": float(np.mean(np.log10(pred) - np.log10(y))),
    }


def build_rows():
    rows = []
    for drug in CLINICAL_PK_COHORT:
        if not (drug.cl_hepatic_ml_min_kg and drug.cl_hepatic_ml_min_kg > 0 and drug.ppb_percent is not None and drug.cl_int_hlm_ml_min_kg and drug.cl_int_hlm_ml_min_kg > 0):
            continue
        mol = Chem.MolFromSmiles(drug.smiles)
        fu = canonical_fu_from_ppb(drug.ppb_percent)
        prior = calculate_well_stirred_clearance(drug.cl_int_hlm_ml_min_kg, fu, rb=drug.blood_to_plasma_rb)["cl_h_plasma_ml_min_kg"]
        rows.append({"name": drug.compound_name, "mol": mol, "y": float(drug.cl_hepatic_ml_min_kg), "prior": max(float(prior), 1e-4), "extraction": drug.extraction_category})
    return sorted(rows, key=lambda r: r["name"])


def factories():
    return {
        "ridge_residual": lambda: make_pipeline(StandardScaler(), Ridge(alpha=10.0)),
        "random_forest_residual": lambda: RandomForestRegressor(n_estimators=500, min_samples_leaf=3, max_features=0.8, random_state=73, n_jobs=-1),
        "extra_trees_residual": lambda: ExtraTreesRegressor(n_estimators=500, min_samples_leaf=3, max_features=0.8, random_state=73, n_jobs=-1),
    }


def persist_control_plane(artifact):
    """Record this distinct route without altering production routing."""
    ledger = {
        "artifact": "GLOBAL_OPTIMIZATION_LEDGER_V73",
        "created_at": artifact["created_at"],
        "production_engine": artifact["production_engine"],
        "previous_ledger": "validation/global_optimization_ledger_v71.json",
        "experiments": [{
            "experiment_id": "V73-HEPATIC-MECHANISTIC-RESIDUAL-HYBRID",
            "endpoint": "HUMAN_HEPATIC_CL",
            "benchmark": "validation/hepatic_hybrid_benchmark_v73.json",
            "canonical_contract": artifact["canonical_contract"],
            "selected_research_candidate": artifact["selected_research_candidate"],
            "decision": artifact["decision"],
        }],
        "next_queue": ["HUMAN_FU_PLASMA", "HUMAN_HEPATIC_CL", "HUMAN_VDSS", "HUMAN_TOTAL_IV_CL", "IV_HALF_LIFE", "ORAL_AUC", "ORAL_CMAX"],
        "production_decision": "UNCHANGED",
    }
    LEDGER_OUT.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")
    board = json.loads(BOARD_OUT.read_text())
    for row in board["rows"]:
        if row["endpoint"] == "HUMAN_HEPATIC_CL":
            row.update({
                "candidate_model": artifact["selected_research_candidate"],
                "current_value": artifact["models"][artifact["selected_research_candidate"]]["pooled_oof"]["aafe"],
                "development_N": artifact["canonical_contract"]["dataset_n"],
                "last_experiment": "validation/hepatic_hybrid_benchmark_v73.json",
                "next_action": "Obtain an independent, prediction-only human hepatic/nonrenal-clearance cohort; the V73 assisted Extra Trees residual hybrid is research-only and must be compared only on its frozen contract.",
                "blocker": "No independent prediction-only hepatic-clearance validation; V73 result is assisted N=30.",
                "source_independence": "assisted/limited; independent validation required",
                "status": "MODEL_OPTIMIZING",
            })
            break
    board["control_plane_version"] = "V73"
    board["last_updated"] = artifact["created_at"]
    BOARD_OUT.write_text(json.dumps(board, indent=2, sort_keys=True) + "\n")


def main():
    rows = build_rows()
    y = np.asarray([r["y"] for r in rows])
    prior = np.asarray([r["prior"] for r in rows])
    x = np.asarray([descriptor_vector(r["mol"]) + [np.log10(r["prior"])] for r in rows], dtype=float)
    groups = np.asarray([MurckoScaffold.MurckoScaffoldSmiles(mol=r["mol"], includeChirality=False) or r["name"] for r in rows])
    splits = list(GroupKFold(n_splits=5).split(x, y, groups))
    predictions = {"well_stirred_observed_fu_hlm": prior.copy()}
    residual = np.log10(y) - np.log10(prior)
    fold_provenance = []
    for name, factory in factories().items():
        pred = np.empty(len(rows))
        for fold, (train, test) in enumerate(splits):
            model = factory()
            model.fit(x[train], residual[train])
            pred[test] = 10 ** (np.log10(prior[test]) + model.predict(x[test]))
            fold_provenance.append({"candidate": name, "fold": fold, "train_n": int(len(train)), "test_n": int(len(test)), "target_fitting_scope": "train_fold_only"})
        predictions[name] = np.clip(pred, 1e-4, None)
    results = {}
    for name, pred in predictions.items():
        results[name] = {"pooled_oof": metrics(y, pred), "by_extraction": {}}
        for extraction in ("LOW", "INTERMEDIATE", "HIGH"):
            ix = [i for i, r in enumerate(rows) if r["extraction"] == extraction]
            results[name]["by_extraction"][extraction] = metrics(y[ix], pred[ix])
    baseline = results["well_stirred_observed_fu_hlm"]["pooled_oof"]
    # A retained research candidate must be materially better and must not
    # regress any prespecified extraction subgroup by more than 10% AAFE.
    retained = []
    for name in factories():
        candidate = results[name]["pooled_oof"]
        subgroup_ok = all(results[name]["by_extraction"][key]["aafe"] <= results["well_stirred_observed_fu_hlm"]["by_extraction"][key]["aafe"] * 1.10 for key in ("LOW", "INTERMEDIATE", "HIGH"))
        if candidate["aafe"] <= baseline["aafe"] * 0.90 and subgroup_ok:
            retained.append(name)
    selected = min(retained, key=lambda n: results[n]["pooled_oof"]["aafe"]) if retained else "well_stirred_observed_fu_hlm"
    artifact = {
        "artifact": "HEPATIC_HYBRID_BENCHMARK_V73", "created_at": datetime.now(timezone.utc).isoformat(),
        "production_engine": "drugopt-prediction-engine-v3@3.3.2", "decision": "RESEARCH_ONLY_NO_PRODUCTION_CHANGE",
        "canonical_contract": {"endpoint": "HUMAN_HEPATIC_CL", "semantics": "human hepatic clearance; not total CL, CL/F, or renal CL", "unit": "mL/min/kg", "dataset": "clinical_pk_cohort.py assisted hepatic-IVIVE cohort", "dataset_hash": hashlib.sha256(json.dumps([(r["name"], r["y"], r["prior"]) for r in rows], sort_keys=True).encode()).hexdigest(), "dataset_n": len(rows), "inputs": "observed PPB/fu and observed HLM Clint permitted only as assisted mechanistic inputs", "target_aggregation": "one curated clinical value per compound", "split": "deterministic 5-fold GroupKFold by Murcko scaffold after compound-name sort", "seed": 73, "target_transform": "log10(CL) residual around well-stirred prior", "inverse_transform": "10**(log10(prior)+residual)", "metric": "pooled out-of-fold fold error", "eligibility": "positive observed hepatic CL, PPB, HLM Clint, valid structure", "ad_evaluation": "all eligible compounds retained; no AD-based exclusions"},
        "models": results, "selection_rule": "retain only >=10% pooled AAFE improvement without >10% extraction-subgroup AAFE regression", "selected_research_candidate": selected, "fold_provenance": fold_provenance,
        "limitations": ["Assisted benchmark, not a prediction-only chain.", "N=30 and no independent validation cohort.", "No candidate is eligible for production promotion from this artifact."],
    }
    OUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    persist_control_plane(artifact)
    print(json.dumps({"n": len(rows), "baseline": baseline, "results": {k: v["pooled_oof"] for k, v in results.items()}, "selected": selected}, indent=2))


if __name__ == "__main__":
    main()
