#!/usr/bin/env python3
"""Independent, endpoint-safe HLM/RLM/MLM model-family audit for v3.3.3.

The public Biogen set is an already-consumed independent diagnostic cohort.
It is never used to tune candidates or promote routing.  HLM, RLM, and MLM
labels remain distinct; the shared candidate uses a multi-output estimator
only on training rows where every head is present.
"""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.admet_predictor import predict_batch_matrix  # noqa: E402
from scripts.evaluate_human_total_iv_cl_v56 import descriptor_vector  # noqa: E402

TRAIN_X = ROOT / "models/openadmet/microsomal_clearance/X_train.csv"
TRAIN_Y = ROOT / "models/openadmet/microsomal_clearance/y_train.csv"
BIOGEN = ROOT / "models/openadmet/validation/biogen_public_3521.csv"
OUT = ROOT / "validation/microsomal_family_audit_v77.json"
FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=512)
RDLogger.DisableLog("rdApp.warning")


def canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(str(smiles))
    return Chem.MolToSmiles(mol, isomericSmiles=True) if mol else ""


def features(smiles: list[str], include_fp: bool) -> np.ndarray:
    rows = []
    for value in smiles:
        mol = Chem.MolFromSmiles(value)
        desc = np.asarray(descriptor_vector(mol), dtype=np.float32)
        if not include_fp:
            rows.append(desc)
            continue
        arr = np.zeros(512, dtype=np.float32)
        DataStructs.ConvertToNumpyArray(FP.GetFingerprint(mol), arr)
        rows.append(np.concatenate((desc, arr)))
    return np.asarray(rows, dtype=np.float32)


def metrics(observed: np.ndarray, predicted: np.ndarray) -> dict:
    error = predicted - observed
    fold = np.power(10.0, np.abs(error))
    return {
        "n": int(len(observed)),
        "mae_log10": float(mean_absolute_error(observed, predicted)),
        "rmse_log10": float(mean_squared_error(observed, predicted) ** 0.5),
        "r2": float(r2_score(observed, predicted)),
        "spearman": float(spearmanr(observed, predicted).statistic),
        "aafe": float(np.mean(fold)),
        "median_fold_error": float(np.median(fold)),
        "within_2_fold_pct": float(np.mean(fold <= 2.0) * 100.0),
        "within_3_fold_pct": float(np.mean(fold <= 3.0) * 100.0),
        "bias_log10": float(np.mean(error)),
    }


def main() -> None:
    train_smiles = pd.read_csv(TRAIN_X).iloc[:, 0].astype(str)
    train_y = pd.read_csv(TRAIN_Y)
    train_keys = set(train_smiles.map(canonical))
    biogen = pd.read_csv(BIOGEN)
    biogen["canonical"] = biogen.SMILES.map(canonical)
    validation = biogen.loc[~biogen.canonical.isin(train_keys) & biogen.canonical.ne("")].copy()

    x_desc = features(train_smiles.tolist(), False)
    x_hybrid = features(train_smiles.tolist(), True)
    v_desc = features(validation.SMILES.astype(str).tolist(), False)
    v_hybrid = features(validation.SMILES.astype(str).tolist(), True)

    # Existing frozen model is evaluated once on the exact same rows.
    production_matrix = np.asarray(predict_batch_matrix(validation.SMILES.astype(str).tolist(), "HLM intrinsic clearance"))
    models = {
        "extra_trees_descriptors": lambda: ExtraTreesRegressor(
            n_estimators=300, min_samples_leaf=3, max_features=0.8,
            random_state=77, n_jobs=2,
        ),
        "extra_trees_descriptor_morgan": lambda: ExtraTreesRegressor(
            n_estimators=240, min_samples_leaf=3, max_features="sqrt",
            random_state=77, n_jobs=2,
        ),
        "hist_gradient_descriptors": lambda: HistGradientBoostingRegressor(
            loss="absolute_error", max_iter=220, max_leaf_nodes=31,
            learning_rate=0.05, l2_regularization=2.0, random_state=77,
        ),
    }
    endpoint_contracts = (
        ("HLM", "LOG_CLint_HLM", "LOG HLM_CLint (mL/min/kg)", 0),
        ("RLM", "LOG_CLint_RLM", "LOG RLM_CLint (mL/min/kg)", 1),
    )
    results: dict[str, dict] = {}
    for endpoint, train_column, validation_column, production_index in endpoint_contracts:
        train_mask = train_y[train_column].notna().to_numpy()
        validation_mask = validation[validation_column].notna().to_numpy()
        observed = validation.loc[validation_mask, validation_column].to_numpy(float)
        endpoint_results = {
            "production_chemeleon": metrics(observed, production_matrix[validation_mask, production_index]),
        }
        for name, factory in models.items():
            use_hybrid = name.endswith("morgan")
            x_train = x_hybrid if use_hybrid else x_desc
            x_test = v_hybrid if use_hybrid else v_desc
            model = factory()
            model.fit(x_train[train_mask], train_y.loc[train_mask, train_column].to_numpy(float))
            endpoint_results[name] = metrics(observed, model.predict(x_test[validation_mask]))
        results[endpoint] = {
            "semantics": f"{endpoint} scaled intrinsic clearance; species heads are not interchangeable",
            "training_n": int(train_mask.sum()),
            "independent_validation_n": int(validation_mask.sum()),
            "models": endpoint_results,
        }

    # Shared representation diagnostic: one estimator, three outputs, with
    # missing labels masked by complete-case training selection.  MLM has no
    # compatible independent column and is therefore not scored or promoted.
    columns = ["LOG_CLint_HLM", "LOG_CLint_RLM", "LOG_CLint_MLM"]
    complete = train_y[columns].notna().all(axis=1).to_numpy()
    shared = ExtraTreesRegressor(
        n_estimators=240, min_samples_leaf=3, max_features="sqrt",
        random_state=78, n_jobs=2,
    )
    shared.fit(x_hybrid[complete], train_y.loc[complete, columns].to_numpy(float))
    shared_prediction = shared.predict(v_hybrid)
    for index, (endpoint, _, validation_column, _) in enumerate(endpoint_contracts):
        mask = validation[validation_column].notna().to_numpy()
        observed = validation.loc[mask, validation_column].to_numpy(float)
        results[endpoint]["models"]["shared_three_head_descriptor_morgan"] = metrics(observed, shared_prediction[mask, index])

    output = {
        "artifact": "MICROSOMAL_FAMILY_AUDIT_V77",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "engine_release": "drugopt-prediction-engine-v3@3.3.3",
        "canonical_contract": {
            "dataset": "Biogen prospective public ADME set; canonical training overlap excluded",
            "dataset_sha256": hashlib.sha256(BIOGEN.read_bytes()).hexdigest(),
            "training_sha256": hashlib.sha256(TRAIN_X.read_bytes() + TRAIN_Y.read_bytes()).hexdigest(),
            "unit": "log10(mL/min/kg)",
            "eligibility": "endpoint-specific non-null label and valid structure",
            "selection_role": "CONSUMED_INDEPENDENT_DIAGNOSTIC_NOT_MODEL_SELECTION",
            "species_isolation": True,
        },
        "results": results,
        "MLM": {
            "training_n": int(train_y.LOG_CLint_MLM.notna().sum()),
            "independent_validation_n": 0,
            "status": "PRODUCTION_STABLE_NO_COMPATIBLE_INDEPENDENT_VALIDATION",
        },
        "decision": "KEEP_EXISTING_PRODUCTION_ROUTING",
        "promotion": False,
        "reason": "Candidates were inspected on an already-consumed independent cohort; a separate locked cohort is required for selection and promotion.",
    }
    OUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({endpoint: value["models"] for endpoint, value in results.items()}, indent=2))


if __name__ == "__main__":
    main()
