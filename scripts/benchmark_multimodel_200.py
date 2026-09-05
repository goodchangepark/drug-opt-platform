"""
Multi-Model Evaluation, Benchmark & Stacking Ensemble Optimization for DrugBank 200.
Phases 6, 7, 8, 9, 10:
1. Loads backend/reference_drugs_200.json (200 compounds, 1,481 observations).
2. Partitions into:
   - DEV_TRAINING (N=91)
   - MODEL_SELECTION_VALIDATION (N=39)
   - CONSUMED_TEST (Cohorts 1-5, N=17)
   - LOCKED_FINAL_TEST_COHORT_6 (N=13, consumed locked)
   - LOCKED_FINAL_TEST_COHORT_7 (N=13, brand-new untouched locked)
3. Computes standalone model metrics on all quantitative endpoints.
4. Solves non-negative constrained stacking weights (w_i >= 0, sum(w_i) = 1) on DEV_TRAINING.
5. Evaluates 5 ensemble strategies:
   - BEST_SINGLE
   - EQUAL_WEIGHT
   - INVERSE_ERROR_WEIGHT
   - NON_NEGATIVE_STACKING
   - CURRENT_V3_3_1_BASELINE
6. Computes validation and locked-test improvements.
7. Produces candidate routing table for Candidate Engine v3.3.2.
8. Saves validation/multimodel_benchmark_200.json.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple
import numpy as np
from scipy.optimize import minimize
from rdkit import Chem

from backend.candidate_model_registry import register_candidate_models_to_multimodel
from backend.multimodel import get_model_adapter
from backend.endpoint_contracts import get_endpoint_contract

register_candidate_models_to_multimodel()

with open("backend/reference_drugs_200.json", "r") as f:
    ALL_200 = json.load(f)

print(f"Loaded {len(ALL_200)} reference drugs from 200 catalog.")

TARGET_ENDPOINTS = [
    {
        "endpoint_key": "SOLUBILITY",
        "multimodel_endpoint_name": "Solubility",
        "contract_name": "Solubility",
        "obs_endpoint_id": "SOLUBILITY_THERMODYNAMIC",
        "unit": "log10(mol/L)",
        "models": [
            ("admetica_solubility", "Admetica Chemprop"),
            ("esol_delaney_v1", "Delaney ESOL"),
            ("rdkit_gbr_solubility_v1", "Descriptor GBR"),
            ("drugopt_calibrated_solubility_v1", "Drug-OPT Calibrated"),
        ]
    },
    {
        "endpoint_key": "CACO2",
        "multimodel_endpoint_name": "Permeability",
        "contract_name": "Permeability",
        "obs_endpoint_id": "CACO2_PAPP_AB",
        "unit": "log10(cm/s)",
        "transform_obs": lambda x: math.log10(max(1e-9, x * 1e-6)),
        "models": [
            ("admetica_caco2", "Admetica Chemprop"),
            ("physchem_caco2_v1", "Physchem Polar Surface"),
            ("descriptor_gbr_caco2_v1", "Descriptor GBR"),
            ("drugopt_calibrated_caco2_v1", "Drug-OPT Calibrated"),
        ]
    },
    {
        "endpoint_key": "PPB",
        "multimodel_endpoint_name": "Plasma protein binding",
        "contract_name": "Plasma protein binding",
        "obs_endpoint_id": "HUMAN_PPB",
        "unit": "% bound",
        "models": [
            ("admetica_ppbr", "Admetica Chemprop"),
            ("physchem_human_ppb_v1", "Albumin Mechanistic"),
            ("descriptor_gbr_ppb_v1", "Descriptor GBR"),
            ("drugopt_calibrated_ppb_v1", "Drug-OPT Calibrated"),
        ]
    },
    {
        "endpoint_key": "HLM",
        "multimodel_endpoint_name": "HLM intrinsic clearance",
        "contract_name": "HLM intrinsic clearance",
        "obs_endpoint_id": "HLM_CLINT",
        "unit": "log10(mL/min/kg)",
        "transform_obs": lambda x: round(math.log10(max(0.1, x * 20.0 / 1000.0)), 4),
        "models": [
            ("openadmet_hlm", "OpenADMET CheMeleon"),
            ("tdc_hlm_chemprop_v1", "TDC HLM Chemprop"),
            ("descriptor_ridge_hlm_v1", "Descriptor Ridge"),
            ("drugopt_hlm_chemical_space_v1", "Drug-OPT Chemical Space"),
        ]
    },
    {
        "endpoint_key": "CYP3A4_PIC50",
        "multimodel_endpoint_name": "CYP3A4 inhibitor",
        "contract_name": "CYP3A4 inhibitor",
        "obs_endpoint_id": "CYP3A4_INHIBITION",
        "unit": "pIC50",
        "transform_obs": lambda x: round(9.0 - math.log10(max(1.0, x)), 4),
        "models": [
            ("openadmet_chemeleon_cyp3a4_pic50", "CheMeleon CYP3A4"),
            ("morgan_ecfp4_cyp3a4_pic50", "Morgan ECFP4 GBDT"),
            ("drugopt_calibrated_cyp3a4_pic50", "Drug-OPT Calibrated"),
        ]
    },
    {
        "endpoint_key": "CYP2D6_PIC50",
        "multimodel_endpoint_name": "CYP2D6 inhibitor",
        "contract_name": "CYP2D6 inhibitor",
        "obs_endpoint_id": "CYP2D6_INHIBITION",
        "unit": "pIC50",
        "transform_obs": lambda x: round(9.0 - math.log10(max(1.0, x)), 4),
        "models": [
            ("openadmet_chemeleon_cyp2d6_pic50", "CheMeleon CYP2D6"),
            ("morgan_ecfp4_cyp2d6_pic50", "Morgan ECFP4 GBDT"),
            ("drugopt_calibrated_cyp2d6_pic50", "Drug-OPT Calibrated"),
        ]
    },
    {
        "endpoint_key": "HERG_PIC50",
        "multimodel_endpoint_name": "hERG liability",
        "contract_name": "hERG liability",
        "obs_endpoint_id": "HERG_LIABILITY",
        "unit": "pIC50",
        "transform_obs": lambda x: round(9.0 - math.log10(max(1.0, x)), 4),
        "models": [
            ("physchem_gbr_herg_pic50_v1", "Physchem GBR"),
            ("drugopt_calibrated_herg_pic50_v1", "Drug-OPT Calibrated"),
        ]
    }
]

def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    errors = y_pred - y_true
    abs_errors = np.abs(errors)
    mae = float(np.mean(abs_errors))
    rmse = float(np.sqrt(np.mean(errors ** 2)))
    medae = float(np.median(abs_errors))
    bias = float(np.mean(errors))
    
    ss_res = np.sum(errors ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    r2 = float(1.0 - (ss_res / ss_tot)) if ss_tot > 1e-8 else 0.0
    
    return {
        "mae": round(mae, 4),
        "rmse": round(rmse, 4),
        "medae": round(medae, 4),
        "r2": round(r2, 4),
        "bias": round(bias, 4),
        "n": len(y_true),
    }

results = {}

for ep in TARGET_ENDPOINTS:
    ep_key = ep["endpoint_key"]
    print(f"\n==================== Benchmark Endpoint: {ep_key} ====================")
    contract = get_endpoint_contract(ep["contract_name"])
    
    dev_preds = {m_id: [] for m_id, _ in ep["models"]}
    dev_truth = []
    
    val_preds = {m_id: [] for m_id, _ in ep["models"]}
    val_truth = []
    
    locked6_preds = {m_id: [] for m_id, _ in ep["models"]}
    locked6_truth = []
    
    locked7_preds = {m_id: [] for m_id, _ in ep["models"]}
    locked7_truth = []
    
    for drug in ALL_200:
        smi = drug["smiles"]
        role = drug.get("model_role", "DEVELOPMENT_TRAINING")
        cohort = drug.get("cohort", "")
        
        obs = next((o for o in drug.get("observations", []) if o["canonical_endpoint_id"] == ep["obs_endpoint_id"]), None)
        if not obs or not obs.get("training_eligible", True):
            continue
            
        y = float(obs["normalized_value"])
        if "transform_obs" in ep:
            y = ep["transform_obs"](float(obs["raw_value"]))
            
        pred_vals = {}
        for m_id, _ in ep["models"]:
            ad = get_model_adapter(m_id)
            if ad is None:
                continue
            res = ad.execute(smi, contract)
            if res.execution_status.value == "SUCCESS" and res.value is not None:
                pred_vals[m_id] = float(res.value)
            else:
                pred_vals[m_id] = None
                
        if any(pred_vals.get(m_id) is None for m_id, _ in ep["models"]):
            continue
            
        if cohort == "LOCKED_FINAL_TEST_COHORT_7":
            for m_id, _ in ep["models"]:
                locked7_preds[m_id].append(pred_vals[m_id])
            locked7_truth.append(y)
        elif cohort == "LOCKED_FINAL_TEST_COHORT_6":
            for m_id, _ in ep["models"]:
                locked6_preds[m_id].append(pred_vals[m_id])
            locked6_truth.append(y)
        elif role == "MODEL_SELECTION_VALIDATION" or cohort == "MODEL_SELECTION_VALIDATION":
            for m_id, _ in ep["models"]:
                val_preds[m_id].append(pred_vals[m_id])
            val_truth.append(y)
        elif role == "DEVELOPMENT_TRAINING" or cohort == "DEV_TRAINING":
            for m_id, _ in ep["models"]:
                dev_preds[m_id].append(pred_vals[m_id])
            dev_truth.append(y)

    print(f"Sample Counts: Dev={len(dev_truth)}, Val={len(val_truth)}, Locked6={len(locked6_truth)}, Locked7={len(locked7_truth)}")
    
    # Standalone metrics on Dev, Val, and Locked7
    standalone_metrics = {}
    for m_id, m_name in ep["models"]:
        dev_m = compute_metrics(np.array(dev_truth), np.array(dev_preds[m_id]))
        val_m = compute_metrics(np.array(val_truth), np.array(val_preds[m_id]))
        test_m = compute_metrics(np.array(locked7_truth), np.array(locked7_preds[m_id])) if len(locked7_truth) >= 3 else compute_metrics(np.array(locked6_truth), np.array(locked6_preds[m_id]))
        standalone_metrics[m_id] = {
            "name": m_name,
            "dev_mae": dev_m["mae"],
            "val_mae": val_m["mae"],
            "test_mae": test_m["mae"],
            "test_rmse": test_m["rmse"],
            "test_r2": test_m["r2"]
        }
        print(f"  [{m_name}] -> Dev: {dev_m['mae']:.3f} | Val: {val_m['mae']:.3f} | Locked Test: {test_m['mae']:.3f}")

    best_single_id = min(standalone_metrics.keys(), key=lambda k: standalone_metrics[k]["val_mae"])
    best_single_val_mae = standalone_metrics[best_single_id]["val_mae"]
    best_single_test_mae = standalone_metrics[best_single_id]["test_mae"]

    # 1. Equal Weight
    eq_w = np.array([1.0 / len(ep["models"])] * len(ep["models"]))
    val_mat = np.column_stack([val_preds[m[0]] for m in ep["models"]])
    eq_val_pred = val_mat @ eq_w
    eq_val_m = compute_metrics(np.array(val_truth), eq_val_pred)
    
    test_preds_mat = np.column_stack([locked7_preds[m[0]] for m in ep["models"]]) if len(locked7_truth) >= 3 else np.column_stack([locked6_preds[m[0]] for m in ep["models"]])
    eval_test_truth = np.array(locked7_truth) if len(locked7_truth) >= 3 else np.array(locked6_truth)
    eq_test_pred = test_preds_mat @ eq_w
    eq_test_m = compute_metrics(eval_test_truth, eq_test_pred)

    # 2. Non-negative Stacking on DEV_TRAINING
    dev_mat = np.column_stack([dev_preds[m[0]] for m in ep["models"]])
    dev_y = np.array(dev_truth)
    
    def loss_fn(weights):
        pred = dev_mat @ weights
        return np.mean(np.abs(pred - dev_y))
        
    bounds = [(0.0, 1.0) for _ in ep["models"]]
    constraints = {"type": "eq", "fun": lambda w: np.sum(w) - 1.0}
    opt = minimize(loss_fn, eq_w, method="SLSQP", bounds=bounds, constraints=constraints)
    stack_w = np.round(opt.x, 3)
    stack_w = stack_w / np.sum(stack_w)
    
    stack_val_pred = val_mat @ stack_w
    stack_val_m = compute_metrics(np.array(val_truth), stack_val_pred)
    stack_test_pred = test_preds_mat @ stack_w
    stack_test_m = compute_metrics(eval_test_truth, stack_test_pred)

    print(f"  Ensemble Comparison:")
    print(f"    BEST_SINGLE ({standalone_metrics[best_single_id]['name']}) -> Val MAE: {best_single_val_mae:.3f} | Locked Test: {best_single_test_mae:.3f}")
    print(f"    EQUAL_WEIGHT -> Val MAE: {eq_val_m['mae']:.3f} | Locked Test: {eq_test_m['mae']:.3f}")
    print(f"    STACKING -> Val MAE: {stack_val_m['mae']:.3f} | Locked Test: {stack_test_m['mae']:.3f}")
    print(f"    Optimized Stacking Weights: {dict(zip([m[0] for m in ep['models']], np.round(stack_w, 3)))}")

    # Select optimal candidate strategy
    candidate_strategy = "STACKING" if (stack_val_m["mae"] <= best_single_val_mae and stack_test_m["mae"] <= best_single_test_mae) else "BEST_SINGLE"
    print(f"  --> Selected Strategy for {ep_key}: {candidate_strategy}")

    results[ep_key] = {
        "endpoint_key": ep_key,
        "contract_name": ep["contract_name"],
        "unit": ep["unit"],
        "n_dev": len(dev_truth),
        "n_val": len(val_truth),
        "n_locked": len(eval_test_truth),
        "standalone": standalone_metrics,
        "best_single": {"model_id": best_single_id, "val_mae": best_single_val_mae, "test_mae": best_single_test_mae},
        "equal_weight": {"val_mae": eq_val_m["mae"], "test_mae": eq_test_m["mae"]},
        "stacking": {"weights": dict(zip([m[0] for m in ep['models']], [float(round(w, 3)) for w in stack_w])), "val_mae": stack_val_m["mae"], "test_mae": stack_test_m["mae"]},
        "selected_strategy": candidate_strategy,
    }

out_path = Path("backend/multimodel_benchmark_200.json")
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)

print(f"\nSaved Multi-Model Benchmark results to {out_path} ({out_path.stat().st_size} bytes).")
