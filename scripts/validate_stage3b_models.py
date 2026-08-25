"""Reproducible independent and reference-compound validation for Stage 3B."""

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.admet_predictor import predict_batch_matrix, predict_batch_values  # noqa: E402


BIOGEN = ROOT / "models/openadmet/validation/biogen_public_3521.csv"


def canonical(smiles):
    mol = Chem.MolFromSmiles(str(smiles))
    return Chem.MolToSmiles(mol, isomericSmiles=True) if mol else ""


def metrics(observed, predicted):
    observed, predicted = np.asarray(observed), np.asarray(predicted)
    return {
        "n": int(len(observed)),
        "MAE": round(float(mean_absolute_error(observed, predicted)), 4),
        "RMSE": round(float(mean_squared_error(observed, predicted) ** 0.5), 4),
        "R2": round(float(r2_score(observed, predicted)), 4),
        "Spearman": round(float(spearmanr(observed, predicted).statistic), 4),
    }


def main():
    if not BIOGEN.is_file():
        raise SystemExit("Packaged Biogen public validation CSV is missing")
    data = pd.read_csv(BIOGEN)
    data["canonical"] = data.SMILES.map(canonical)

    clearance_train = set(pd.read_csv(ROOT / "models/openadmet/microsomal_clearance/X_train.csv").iloc[:, 0].map(canonical))
    clearance = data.loc[~data.canonical.isin(clearance_train)].copy()
    matrix = np.asarray(predict_batch_matrix(clearance.SMILES.tolist(), "HLM intrinsic clearance"))
    clearance_metrics = {}
    for endpoint, column, index in (
        ("HLM", "LOG HLM_CLint (mL/min/kg)", 0),
        ("RLM", "LOG RLM_CLint (mL/min/kg)", 1),
    ):
        mask = clearance[column].notna().to_numpy()
        clearance_metrics[endpoint] = metrics(clearance.loc[mask, column], matrix[mask, index])

    ppb_train = set(pd.read_csv(ROOT / "models/admetica/ppbr/training.csv").Drug.map(canonical))
    ppb_column = "LOG PLASMA PROTEIN BINDING (HUMAN) (% unbound)"
    ppb = data.loc[data[ppb_column].notna() & ~data.canonical.isin(ppb_train)].copy()
    observed_bound = 100.0 - np.power(10.0, ppb[ppb_column].to_numpy())
    predicted_bound = predict_batch_values(ppb.SMILES.tolist(), "Plasma protein binding")

    references = [
        {"name": "Rifampicin", "smiles": "C[C@H]1/C=C/C=C(\\C(=O)NC2=C(C(=C3C(=C2O)C(=C(C4=C3C(=O)[C@](O4)(O/C=C/[C@@H]([C@H]([C@H]([C@@H]([C@@H]([C@@H]([C@H]1O)C)O)C)OC(=O)C)C)OC)C)C)O)O)/C=N/N5CCN(CC5)C)/C", "experimental_raw_uL_min_mg": 2.84},
        {"name": "Isoniazid", "smiles": "NNC(=O)c1ccncc1", "experimental_raw_uL_min_mg": 13.9},
        {"name": "Ethionamide", "smiles": "CCC1=NC=CC(=C1)C(=S)N", "experimental_raw_uL_min_mg": 77.1},
    ]
    reference_predictions = predict_batch_values([row["smiles"] for row in references], "HLM intrinsic clearance")
    for row, prediction in zip(references, reference_predictions):
        # Same explicit human scaling factors used in the OpenADMET paired source data.
        row["experimental_scaled_log10_mL_min_kg"] = round(float(np.log10(row["experimental_raw_uL_min_mg"] * 45 * 26 / 1000)), 4)
        row["predicted_scaled_log10_mL_min_kg"] = round(prediction, 4)
    reference_spearman = float(spearmanr(
        [row["experimental_scaled_log10_mL_min_kg"] for row in references],
        [row["predicted_scaled_log10_mL_min_kg"] for row in references],
    ).statistic)

    result = {
        "independent_dataset": "Biogen prospective public ADME set (3,521 compounds; canonical training-overlap excluded)",
        "dataset_license": "MIT",
        "clearance": clearance_metrics,
        "human_ppb_percent_bound": metrics(observed_bound, predicted_bound),
        "reference_compounds": references,
        "reference_directionality_spearman": round(reference_spearman, 4),
        "reference_source": "Lakshminarayana et al., J Med Chem 2020 anti-TB profile; raw HLM µL/min/mg. Scaling: 45 mg microsomal protein/g liver and 26 g liver/kg.",
        "limitations": "The Biogen set is independent by canonical-SMILES exclusion, but no publication-level guarantee rules out series or analogue overlap. The three named compounds are a small directionality sanity check, not a performance estimate.",
    }
    output = ROOT / "models/openadmet/microsomal_clearance/independent_validation.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
