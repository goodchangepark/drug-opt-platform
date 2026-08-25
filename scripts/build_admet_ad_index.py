"""Rebuild packaged Stage 3A applicability-domain indexes from training CSV files."""

import csv
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator


ROOT = Path(__file__).resolve().parents[1]
MODEL_ROOT = ROOT / "models" / "admetica"
GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def descriptors(mol):
    return [
        Descriptors.MolWt(mol), Crippen.MolLogP(mol), Descriptors.TPSA(mol),
        Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol), Lipinski.NumRotatableBonds(mol),
    ]


RDLogger.DisableLog("rdApp.*")
for endpoint in ("solubility", "caco2"):
    fingerprints, descriptor_rows = [], []
    with (MODEL_ROOT / endpoint / "training.csv").open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            mol = Chem.MolFromSmiles(row.get("Drug", ""))
            if mol is None:
                continue
            array = np.zeros(2048, dtype=np.uint8)
            DataStructs.ConvertToNumpyArray(GENERATOR.GetFingerprint(mol), array)
            fingerprints.append(array)
            descriptor_rows.append(descriptors(mol))
    matrix = np.asarray(fingerprints, dtype=np.uint8)
    values = np.asarray(descriptor_rows, dtype=np.float32)
    np.savez_compressed(
        MODEL_ROOT / endpoint / "ad_index.npz",
        fingerprints=matrix,
        bit_counts=matrix.sum(axis=1),
        descriptor_min=values.min(axis=0),
        descriptor_max=values.max(axis=0),
    )
    print(endpoint, matrix.shape)
