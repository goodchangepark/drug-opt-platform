"""Rebuild packaged ADMET applicability-domain indexes from training CSV files."""

import csv
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator


ROOT = Path(__file__).resolve().parents[1]
GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def descriptors(mol):
    return [
        Descriptors.MolWt(mol), Crippen.MolLogP(mol), Descriptors.TPSA(mol),
        Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol), Lipinski.NumRotatableBonds(mol),
    ]


RDLogger.DisableLog("rdApp.*")
datasets = [
    (ROOT / "models/admetica/solubility/training.csv", "Drug", None, ROOT / "models/admetica/solubility/ad_index.npz"),
    (ROOT / "models/admetica/caco2/training.csv", "Drug", None, ROOT / "models/admetica/caco2/ad_index.npz"),
    (ROOT / "models/admetica/ppbr/training.csv", "Drug", None, ROOT / "models/admetica/ppbr/ad_index.npz"),
]
for cyp_key in (
    "cyp1a2-inhibitor", "cyp2c9-inhibitor", "cyp2c19-inhibitor",
    "cyp2d6-inhibitor", "cyp3a4-inhibitor", "cyp2c9-substrate",
    "cyp2d6-substrate", "cyp3a4-substrate",
):
    cyp_root = ROOT / "models" / "admetica" / "cyp" / cyp_key
    datasets.append((cyp_root / "training.csv", "smiles", None, cyp_root / "ad_index.npz"))
clearance_root = ROOT / "models/openadmet/microsomal_clearance"
for species in ("HLM", "RLM", "MLM"):
    datasets.append((
        clearance_root / "X_train.csv", "OPENADMET_CANONICAL_SMILES", f"LOG_CLint_{species}",
        clearance_root / f"ad_index_{species.lower()}.npz",
    ))

for csv_path, smiles_column, target_column, output_path in datasets:
    fingerprints, descriptor_rows = [], []
    target_values = None
    if target_column:
        with (clearance_root / "y_train.csv").open(newline="", encoding="utf-8") as stream:
            target_values = [row.get(target_column, "") for row in csv.DictReader(stream)]
    with csv_path.open(newline="", encoding="utf-8") as stream:
        for index, row in enumerate(csv.DictReader(stream)):
            if target_values is not None and not target_values[index].strip():
                continue
            mol = Chem.MolFromSmiles(row.get(smiles_column, ""))
            if mol is None:
                continue
            array = np.zeros(2048, dtype=np.uint8)
            DataStructs.ConvertToNumpyArray(GENERATOR.GetFingerprint(mol), array)
            fingerprints.append(array)
            descriptor_rows.append(descriptors(mol))
    matrix = np.asarray(fingerprints, dtype=np.uint8)
    values = np.asarray(descriptor_rows, dtype=np.float32)
    np.savez_compressed(
        output_path,
        fingerprints=matrix,
        bit_counts=matrix.sum(axis=1),
        descriptor_min=values.min(axis=0),
        descriptor_max=values.max(axis=0),
    )
    print(output_path.relative_to(ROOT), matrix.shape)
