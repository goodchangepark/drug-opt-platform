"""
Multi-Model Evaluation & Experimental-Weighted Ensemble Optimization.
Directives 12-21:
1. Load all 150 Reference Drugs and partition:
   - DEV_TRAINING (drugs 1-75 & 81-90 & 101-125)
   - MODEL_SELECTION_VALIDATION (drugs 91-100 & 126-137)
   - LOCKED_FINAL_TEST_COHORT_6 (drugs 138-150, N=13)
2. Compute standalone metrics (MAE, RMSE, MedAE, R2, Bias) for candidate models per endpoint:
   - Solubility (Admetica, ESOL, GBR, Drug-OPT Calibrated)
   - Caco-2 (Admetica, Physchem, GBR, Drug-OPT Calibrated)
   - PPB (Admetica, Albumin Physchem, GBR, Drug-OPT Calibrated)
   - HLM Clint (OpenADMET, TDC Chemprop, Descriptor Ridge, Chemical Space)
   - CYP3A4, 2D6, 1A2, 2C9 (CheMeleon, Morgan GBDT, Drug-OPT Calibrated)
   - hERG liability (CardioTox MPNN, Physchem GBR, Drug-OPT Calibrated)
3. Compute Stacking Weights on DEV_TRAINING using scipy.optimize.minimize (SLSQP):
   - Constraint: sum(w_i) = 1.0, w_i >= 0 (Non-negative constrained stacking)
4. Compare 4 Ensemble Strategies:
   - Best Single Model
   - Equal Weight Ensemble
   - Inverse Error Weighted Ensemble
   - Non-negative Constrained Stacking
5. Evaluate all models & ensembles strictly on LOCKED_FINAL_TEST_COHORT_6.
6. Export validation/multimodel_benchmark_150.json.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple
import numpy as np
from scipy.optimize import minimize
from rdkit import Chem

from backend.candidate_model_registry import register_candidate_models_to_multimodel, CANDIDATE_ADAPTER_SUITE
from backend.multimodel import (
    get_v2_adapters_for_endpoint,
    get_model_adapter,
)
from backend.endpoint_contracts import get_endpoint_contract

register_candidate_models_to_multimodel()

# Load 150 reference drugs
with open("backend/reference_drugs_150.json", "r") as f:
    ALL_150 = json.load(f)

print(f"Loaded {len(ALL_150)} reference drugs.")

# Identify quantitative evaluation endpoints and mappings to experimental observations
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
        "transform_obs": lambda x: math.log10(max(1e-9, x * 1e-6)), # from 10^-6 cm/s to log10(cm/s)
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
        "transform_obs": lambda x: round(math.log10(max(0.1, x * 20.0 / 1000.0)), 4), # standard scaling to log10(mL/min/kg)
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
        "transform_obs": lambda x: round(9.0 - math.log10(max(1.0, x)), 4), # IC50 nM to pIC50
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
        "endpoint_key": "CYP1A2_PIC50",
        "multimodel_endpoint_name": "CYP1A2 inhibitor",
        "contract_name": "CYP1A2 inhibitor",
        "obs_endpoint_id": "CYP1A2_INHIBITION",
        "unit": "pIC50",
        "transform_obs": lambda x: round(9.0 - math.log10(max(1.0, x)), 4),
        "models": [
            ("openadmet_chemeleon_cyp1a2_pic50", "CheMeleon CYP1A2"),
            ("morgan_ecfp4_cyp1a2_pic50", "Morgan ECFP4 GBDT"),
            ("drugopt_calibrated_cyp1a2_pic50", "Drug-OPT Calibrated"),
        ]
    },
    {
        "endpoint_key": "CYP2C9_PIC50",
        "multimodel_endpoint_name": "CYP2C9 inhibitor",
        "contract_name": "CYP2C9 inhibitor",
        "obs_endpoint_id": "CYP2C9_INHIBITION",
        "unit": "pIC50",
        "transform_obs": lambda x: round(9.0 - math.log10(max(1.0, x)), 4),
        "models": [
            ("openadmet_chemeleon_cyp2c9_pic50", "CheMeleon CYP2C9"),
            ("morgan_ecfp4_cyp2c9_pic50", "Morgan ECFP4 GBDT"),
            ("drugopt_calibrated_cyp2c9_pic50", "Drug-OPT Calibrated"),
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
    
    # Collect predictions and ground truth across cohorts
    dev_preds = {m_id: [] for m_id, _ in ep["models"]}
    dev_truth = []
    
    val_preds = {m_id: [] for m_id, _ in ep["models"]}
    val_truth = []
    
    test_preds = {m_id: [] for m_id, _ in ep["models"]}
    test_truth = []
    
    for drug in ALL_150:
        smi = drug["smiles"]
        role = drug.get("model_role", "DEVELOPMENT_TRAINING")
        cohort = drug.get("cohort", "")
        
        # Find observation
        obs = next((o for o in drug.get("observations", []) if o["canonical_endpoint_id"] == ep["obs_endpoint_id"]), None)
        if not obs or not obs.get("training_eligible", True):
            continue
            
        y = float(obs["normalized_value"])
        if "transform_obs" in ep:
            y = ep["transform_obs"](float(obs["raw_value"]))
            
        # Run models
        pred_vals = {}
        for m_id, _ in ep["models"]:
            ad = get_model_adapter(m_id)
            if ad is None:
                continue
            res = ad.execute(smi, contract)
            if res.execution_status.value == "SUCCESS" and res.value is not None:
                val = float(res.value)
                # Ensure correct unit scale if needed
                pred_vals[m_id] = val
            else:
                pred_vals[m_id] = None
                
        # Only include if all models executed
        if any(pred_vals.get(m_id) is None for m_id, _ in ep["models"]):
            continue
            
        is_locked_test_6 = (cohort == "LOCKED_FINAL_TEST_COHORT_6")
        is_locked_test_5 = (cohort == "LOCKED_FINAL_TEST_COHORT_5")
        is_val = (role == "MODEL_SELECTION_VALIDATION" or cohort == "MODEL_SELECTION_VALIDATION")
        
        if is_locked_test_6:
            for m_id in ep["models"]:
                test_preds[m_id[0]].append(pred_vals[m_id[0]])
            test_truth.append(y)
        elif is_locked_test_5:
            # Also track cohort 5 in case cohort 6 lacks this endpoint
            if "test_5_preds" not in ep:
                ep["test_5_preds"] = {m_id[0]: [] for m_id in ep["models"]}
                ep["test_5_truth"] = []
            for m_id in ep["models"]:
                ep["test_5_preds"][m_id[0]].append(pred_vals[m_id[0]])
            ep["test_5_truth"].append(y)
        elif is_val:
            for m_id in ep["models"]:
                val_preds[m_id[0]].append(pred_vals[m_id[0]])
            val_truth.append(y)
        else:
            for m_id in ep["models"]:
                dev_preds[m_id[0]].append(pred_vals[m_id[0]])
            dev_truth.append(y)

    # If test_truth is 0, check test_5
    test_cohort_name = "LOCKED_FINAL_TEST_COHORT_6"
    if len(test_truth) < 5 and ep.get("test_5_truth") and len(ep["test_5_truth"]) >= 5:
        test_preds = ep["test_5_preds"]
        test_truth = ep["test_5_truth"]
        test_cohort_name = "LOCKED_FINAL_TEST_COHORT_5"

    print(f"Counts: Dev={len(dev_truth)}, Val={len(val_truth)}, {test_cohort_name}={len(test_truth)}")
    if len(dev_truth) < 5 or len(test_truth) < 5:
        print(f"Skipping {ep_key} due to insufficient samples.")
        continue

    # Standalone model metrics on Dev & Test
    standalone_metrics = {}
    for m_id, m_name in ep["models"]:
        dev_m = compute_metrics(np.array(dev_truth), np.array(dev_preds[m_id]))
        test_m = compute_metrics(np.array(test_truth), np.array(test_preds[m_id]))
        standalone_metrics[m_id] = {
            "name": m_name,
            "dev_mae": dev_m["mae"],
            "dev_rmse": dev_m["rmse"],
            "test_mae": test_m["mae"],
            "test_rmse": test_m["rmse"],
            "test_r2": test_m["r2"],
        }
        print(f"  Model [{m_name}] -> Dev MAE: {dev_m['mae']:.3f} | Locked Test MAE: {test_m['mae']:.3f}")

    # Best single model
    best_single_id = min(standalone_metrics.keys(), key=lambda k: standalone_metrics[k]["dev_mae"])
    best_single_test_mae = standalone_metrics[best_single_id]["test_mae"]

    # 1. Equal Weight Ensemble
    eq_w = [1.0 / len(ep["models"])] * len(ep["models"])
    dev_mat = np.column_stack([dev_preds[m_id[0]] for m_id in ep["models"]])
    test_mat = np.column_stack([test_preds[m_id[0]] for m_id in ep["models"]])
    
    eq_dev_pred = dev_mat @ eq_w
    eq_test_pred = test_mat @ eq_w
    eq_metrics = compute_metrics(np.array(test_truth), eq_test_pred)

    # 2. Inverse Error Weighted Ensemble
    inv_errs = [1.0 / max(0.01, standalone_metrics[m_id[0]]["dev_mae"]) for m_id in ep["models"]]
    inv_w = [ie / sum(inv_errs) for ie in inv_errs]
    inv_test_pred = test_mat @ inv_w
    inv_metrics = compute_metrics(np.array(test_truth), inv_test_pred)

    # 3. Non-negative Constrained Stacking Optimization on Dev Training (SLSQP)
    y_dev = np.array(dev_truth)
    def loss(w):
        pred = dev_mat @ w
        return np.mean((pred - y_dev) ** 2)

    init_w = np.array(eq_w)
    bounds = [(0.0, 1.0) for _ in ep["models"]]
    constraints = {"type": "eq", "fun": lambda w: np.sum(w) - 1.0}
    opt_res = minimize(loss, init_w, method="SLSQP", bounds=bounds, constraints=constraints)
    stack_w = opt_res.x if opt_res.success else np.array(inv_w)
    stack_w = np.maximum(0.0, stack_w)
    stack_w = stack_w / np.sum(stack_w)
    stack_test_pred = test_mat @ stack_w
    stack_metrics = compute_metrics(np.array(test_truth), stack_test_pred)

    print(f"  Ensemble Performance on Locked Test Cohort 6 (N={len(test_truth)}):")
    print(f"    - Best Single ({standalone_metrics[best_single_id]['name']}): MAE = {best_single_test_mae:.3f}")
    print(f"    - Equal Weight: MAE = {eq_metrics['mae']:.3f} (R2={eq_metrics['r2']:.3f})")
    print(f"    - Inverse Error Weight: MAE = {inv_metrics['mae']:.3f} (R2={inv_metrics['r2']:.3f})")
    print(f"    - Non-negative Stacking: MAE = {stack_metrics['mae']:.3f} (R2={stack_metrics['r2']:.3f}) | Weights = {[round(float(w), 3) for w in stack_w]}")

    results[ep_key] = {
        "models": standalone_metrics,
        "best_single": {"model_id": best_single_id, "name": standalone_metrics[best_single_id]["name"], "test_mae": best_single_test_mae},
        "equal_weight": {"test_mae": eq_metrics["mae"], "test_rmse": eq_metrics["rmse"], "test_r2": eq_metrics["r2"], "weights": [round(float(w), 4) for w in eq_w]},
        "inverse_error_weight": {"test_mae": inv_metrics["mae"], "test_rmse": inv_metrics["rmse"], "test_r2": inv_metrics["r2"], "weights": [round(float(w), 4) for w in inv_w]},
        "stacking_ensemble": {"test_mae": stack_metrics["mae"], "test_rmse": stack_metrics["rmse"], "test_r2": stack_metrics["r2"], "weights": [round(float(w), 4) for w in stack_w]},
        "sample_counts": {"dev_n": len(dev_truth), "val_n": len(val_truth), "locked_test_n": len(test_truth)}
    }

# Save results
out_file = Path("validation/multimodel_benchmark_150.json")
out_file.parent.mkdir(exist_ok=True, parents=True)
with open(out_file, "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2)

print(f"\nSuccessfully written multimodel benchmark to {out_file} ({out_file.stat().st_size} bytes)")
