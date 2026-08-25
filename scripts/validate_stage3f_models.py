#!/usr/bin/env python3
"""Reproduce Stage 3F overlap-filtered hERG metrics and public-reference sanity checks."""
import csv
import json
from pathlib import Path
import sys

import numpy as np
from sklearn.metrics import (
    average_precision_score, balanced_accuracy_score, confusion_matrix,
    matthews_corrcoef, roc_auc_score,
)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.admet_predictor import predict_batch_values, predict_endpoint  # noqa: E402

HERG = ROOT / "models/admetica/validation/safety/chembl37_herg_ic50_no_exact_training_overlap.csv"
ACCEPTANCE = ROOT / "validation/stage3_acceptance_dataset.csv"


def main():
    with HERG.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    y = np.asarray([int(row["label"]) for row in rows])
    p = np.asarray(predict_batch_values([row["smiles"] for row in rows], "hERG liability"))
    predicted = (p >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, predicted, labels=[0, 1]).ravel()
    metrics = {
        "n": len(y), "AUROC": roc_auc_score(y, p), "AUPRC": average_precision_score(y, p),
        "balanced_accuracy": balanced_accuracy_score(y, predicted),
        "sensitivity": tp / (tp + fn), "specificity": tn / (tn + fp),
        "MCC": matthews_corrcoef(y, predicted),
        "confusion_matrix": {"tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)},
    }
    sanity = []
    with ACCEPTANCE.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            result = predict_endpoint(row["canonical_smiles"], row["endpoint"])
            expected = int(row["experimental_value"])
            observed = int(result["probability"] >= result["decision_threshold"])
            sanity.append({
                "compound": row["name"], "endpoint": row["endpoint"], "expected_class": expected,
                "probability": result["probability"], "predicted_class": observed,
                "direction_match": expected == observed, "domain": result["applicability_domain"]["classification"],
                "nearest_similarity": result["applicability_domain"]["nearest_training_similarity"],
                "independent_status": row["independent_status"],
            })
    print(json.dumps({"hERG_independent_validation": metrics, "known_compound_sanity": sanity}, indent=2))


if __name__ == "__main__":
    main()
