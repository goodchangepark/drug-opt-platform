"""Reproduce Stage 3A scientific sanity checks with public reference compounds."""

import csv
import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
RDLogger.DisableLog("rdApp.*")

from backend.admet_predictor import MODEL_ROOT, predict_batch_values  # noqa: E402


def canonical(value):
    mol = Chem.MolFromSmiles(value)
    return Chem.MolToSmiles(mol) if mol else None


with (MODEL_ROOT / "validation" / "caco2_external_34.csv").open(newline="", encoding="utf-8") as stream:
    external = list(csv.DictReader(stream))
caco_smiles = [row["SMILES"] for row in external]
caco_true_log_cm_s = [float(row["Papp(original)a"]) - 6 for row in external]
caco_pred = predict_batch_values(caco_smiles, "Permeability")

with (MODEL_ROOT / "solubility" / "training.csv").open(newline="", encoding="utf-8") as stream:
    solubility_training = {canonical(row["Drug"]) for row in csv.DictReader(stream)}
solubility_references = {
    "ethanol": "CCO",
    "4,4'-dichlorobiphenyl": "Clc1ccc(cc1)c2ccc(Cl)cc2",
}
solubility_pred = dict(zip(solubility_references, predict_batch_values(list(solubility_references.values()), "Solubility")))

report = {
    "caco2_external_34": {
        "source": "Pham-The et al. 2011 supplementary data, filtered by Admetica to remove its training overlap",
        "n": len(caco_true_log_cm_s),
        "unit": "log10(Papp [cm/s])",
        "MAE": mean_absolute_error(caco_true_log_cm_s, caco_pred),
        "RMSE": mean_squared_error(caco_true_log_cm_s, caco_pred) ** 0.5,
        "R2": r2_score(caco_true_log_cm_s, caco_pred),
        "direction_correct_high_vs_low": caco_pred[int(np.argmax(caco_true_log_cm_s))] > caco_pred[int(np.argmin(caco_true_log_cm_s))],
    },
    "solubility_directional_sanity": {
        "unit": "log10(mol/L)",
        "predictions": solubility_pred,
        "expected_direction": "ethanol > 4,4'-dichlorobiphenyl",
        "direction_correct": solubility_pred["ethanol"] > solubility_pred["4,4'-dichlorobiphenyl"],
        "training_overlap": {name: canonical(smiles) in solubility_training for name, smiles in solubility_references.items()},
        "limitation": "Directional sanity check only; AqSolDB is broad and both public reference structures overlap the packaged training set.",
    },
}
print(json.dumps(report, indent=2))
