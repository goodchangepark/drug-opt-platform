"""Reproduce CYP validation with exact canonical training-overlap removal."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from sklearn.metrics import (average_precision_score, balanced_accuracy_score,
                             confusion_matrix, matthews_corrcoef, roc_auc_score)


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.admet_predictor import MODEL_SPECS, predict_batch_values  # noqa: E402


VALIDATION_ROOT = ROOT / "models" / "admetica" / "validation" / "cyp"


def canonical(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(str(smiles))
    return Chem.MolToSmiles(mol, isomericSmiles=True) if mol else None


def training_smiles(endpoint: str) -> set[str]:
    path = ROOT / "models" / "admetica" / MODEL_SPECS[endpoint]["model_key"] / "training.csv"
    return {value for value in pd.read_csv(path).smiles.map(canonical) if value}


def classification_metrics(labels, probabilities) -> dict:
    labels = np.asarray(labels, dtype=int)
    probabilities = np.asarray(probabilities, dtype=float)
    predicted = (probabilities >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(labels, predicted, labels=[0, 1]).ravel()
    return {
        "n": int(len(labels)),
        "both_classes": bool(len(set(labels)) == 2),
        "decision_threshold": 0.5,
        "AUROC": round(float(roc_auc_score(labels, probabilities)), 4),
        "AUPRC": round(float(average_precision_score(labels, probabilities)), 4),
        "balanced_accuracy": round(float(balanced_accuracy_score(labels, predicted)), 4),
        "sensitivity": round(float(tp / (tp + fn)), 4),
        "specificity": round(float(tn / (tn + fp)), 4),
        "MCC": round(float(matthews_corrcoef(labels, predicted)), 4),
        "confusion_matrix": {"TN": int(tn), "FP": int(fp), "FN": int(fn), "TP": int(tp)},
    }


def validate_inhibitor(isoform: str) -> dict:
    endpoint = f"CYP{isoform} inhibitor"
    path = VALIDATION_ROOT / f"chembl30_{isoform.lower()}_inhibitor.csv"
    frame = pd.read_csv(path).dropna(subset=["smiles", "class"])
    frame["canonical"] = frame.smiles.map(canonical)
    frame = frame.dropna(subset=["canonical"]).drop_duplicates("canonical")
    training = training_smiles(endpoint)
    overlap = int(frame.canonical.isin(training).sum())
    independent = frame[~frame.canonical.isin(training)].copy()
    probabilities = predict_batch_values(independent.canonical.tolist(), endpoint)
    return {
        "dataset": f"ChEMBL 30 {endpoint} set distributed with Admetica",
        "source": "https://github.com/datagrok-ai/admetica/tree/master/comparison/novartis/cyp",
        "license": "MIT",
        "canonical_training_overlap_removed": overlap,
        **classification_metrics(independent["class"].astype(int).tolist(), probabilities),
    }


def validate_substrate_sanity() -> dict:
    endpoint = "CYP3A4 substrate"
    frame = pd.read_csv(VALIDATION_ROOT / "fda_tki_cyp3a4_substrate.csv")
    frame["canonical"] = frame.smiles.map(canonical)
    frame = frame.dropna(subset=["canonical"]).drop_duplicates("canonical")
    training = training_smiles(endpoint)
    overlap = int(frame.canonical.isin(training).sum())
    independent = frame[~frame.canonical.isin(training)].copy()
    probabilities = predict_batch_values(independent.canonical.tolist(), endpoint)
    predicted = np.asarray(probabilities) >= 0.5
    return {
        "dataset": "24 FDA-approved tyrosine kinase inhibitors with experimental CYP3A4 substrate labels",
        "source": "https://github.com/datagrok-ai/admetica/tree/master/comparison/predictors",
        "license": "MIT",
        "canonical_training_overlap_removed": overlap,
        "n": int(len(independent)),
        "both_classes": False,
        "positive_only": True,
        "decision_threshold": 0.5,
        "sensitivity": round(float(predicted.mean()), 4),
        "note": "Positive-only directionality sanity set; AUROC, AUPRC, specificity, balanced accuracy and MCC are not identifiable.",
    }


def main():
    RDLogger.DisableLog("rdApp.*")
    output = {
        "method": "Canonical isomeric-SMILES exact overlap removal; probabilities evaluated at fixed threshold 0.5.",
        "inhibitors": {f"CYP{iso}": validate_inhibitor(iso) for iso in ("2C9", "2D6", "3A4")},
        "substrates": {"CYP3A4": validate_substrate_sanity()},
        "unvalidated": {
            "CYP1A2 inhibitor": "No qualified independent public set packaged in this step.",
            "CYP2C19 inhibitor": "No qualified independent public set packaged in this step.",
            "CYP2C9 substrate": "No qualified independent public set packaged in this step.",
            "CYP2D6 substrate": "No qualified independent public set packaged in this step.",
        },
    }
    path = VALIDATION_ROOT / "independent_validation.json"
    path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
