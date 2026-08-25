"""Endpoint-specific Stage 3A ADMET inference and scientific guardrails."""

from __future__ import annotations

import math
import threading
from functools import lru_cache
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator


ROOT = Path(__file__).resolve().parents[1]
MODEL_ROOT = ROOT / "models" / "admetica"
MODEL_VERSION = "admetica-d4f7056-chemprop-v2.1"

MODEL_SPECS = {
    "Solubility": {
        "model_key": "solubility",
        "display_name": "Admetica Chemprop Solubility",
        "endpoint_definition": "Aqueous solubility LogS = log10(S [mol/L]); aggregate aqueous measurements, not a pH-specific or intrinsic-solubility estimate.",
        "unit": "log10(mol/L)",
        "training_dataset": "AqSolDB curated aqueous solubility (9 merged public sources; 9,982 rows)",
        "validation": {"MAE": 0.714, "RMSE": 1.089, "R2": 0.788, "Spearman": 0.897},
        "source": "https://github.com/datagrok-ai/admetica",
        "license": "MIT",
        "limitations": "AqSolDB combines measurements made under heterogeneous pH, temperature, protocol, and solid-state conditions. No pKa is used; this is not a pH-specific prediction.",
    },
    "Permeability": {
        "model_key": "caco2",
        "display_name": "Admetica Chemprop Caco-2",
        "endpoint_definition": "Caco-2 apparent permeability LogPapp = log10(Papp [cm/s]). The aggregated training file does not retain A→B/B→A direction or detailed assay conditions.",
        "unit": "log10(cm/s)",
        "training_dataset": "Wang et al. Caco-2 compiled permeability set as distributed by Admetica (910 rows)",
        "validation": {"MAE": 0.317, "RMSE": 0.415, "R2": 0.701, "Spearman": 0.832, "external_MAE_34_compounds": 0.412},
        "source": "https://github.com/datagrok-ai/admetica",
        "license": "MIT",
        "limitations": "Assay direction and conditions are absent from the aggregate source. The value must not be treated as PAMPA, MDCK, or an efflux ratio.",
    },
}

_MODEL_LOCK = threading.Lock()
_MODELS: dict[str, object] = {}
_TRAINER = None
_FP_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
_DESCRIPTOR_NAMES = ("MW", "cLogP", "TPSA", "HBD", "HBA", "RotB")


def registry_seed(endpoint: str) -> dict:
    spec = MODEL_SPECS[endpoint]
    return {
        "endpoint_name": endpoint,
        "model_name": spec["display_name"],
        "model_version": MODEL_VERSION,
        "implementation_status": "READY",
        "supported_species": ["in vitro Caco-2"] if endpoint == "Permeability" else [],
        "supported_matrix": ["Caco-2"] if endpoint == "Permeability" else ["aqueous"],
        "output_unit": spec["unit"],
        "provenance_json": {key: value for key, value in spec.items() if key not in {"model_key", "display_name", "unit"}},
        "is_active": True,
    }


def model_files_available(endpoint: str) -> tuple[bool, str]:
    key = MODEL_SPECS[endpoint]["model_key"]
    missing = [name for name in ("model_v2_1.pt", "training.csv", "ad_index.npz") if not (MODEL_ROOT / key / name).is_file()]
    if missing:
        return False, "Packaged model asset missing: " + ", ".join(missing)
    try:
        import chemprop  # noqa: F401
        import torch  # noqa: F401
    except ImportError as exc:
        return False, f"Runtime dependency unavailable: {exc.name}"
    return True, ""


def _descriptors(mol) -> dict[str, float]:
    return {
        "MW": float(Descriptors.MolWt(mol)),
        "cLogP": float(Crippen.MolLogP(mol)),
        "TPSA": float(Descriptors.TPSA(mol)),
        "HBD": float(Lipinski.NumHDonors(mol)),
        "HBA": float(Lipinski.NumHAcceptors(mol)),
        "RotB": float(Lipinski.NumRotatableBonds(mol)),
    }


@lru_cache(maxsize=2)
def _training_index(model_key: str):
    with np.load(MODEL_ROOT / model_key / "ad_index.npz") as index:
        fingerprints = index["fingerprints"]
        bit_counts = index["bit_counts"]
        descriptor_min = index["descriptor_min"]
        descriptor_max = index["descriptor_max"]
    ranges = {name: {"min": float(descriptor_min[i]), "max": float(descriptor_max[i])} for i, name in enumerate(_DESCRIPTOR_NAMES)}
    return fingerprints, bit_counts, ranges


def applicability_domain(smiles: str, endpoint: str) -> dict:
    """Heuristic AD: nearest training-set Morgan similarity plus descriptor envelope."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")
    spec = MODEL_SPECS[endpoint]
    training_fps, training_bit_counts, ranges = _training_index(spec["model_key"])
    fingerprint = _FP_GENERATOR.GetFingerprint(mol)
    query = np.zeros(2048, dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fingerprint, query)
    intersections = training_fps @ query
    unions = training_bit_counts + int(query.sum()) - intersections
    nearest = float(np.max(np.divide(intersections, unions, out=np.zeros_like(intersections, dtype=float), where=unions != 0)))
    values = _descriptors(mol)
    outside = [name for name in _DESCRIPTOR_NAMES if not ranges[name]["min"] <= values[name] <= ranges[name]["max"]]
    if nearest >= 0.40 and not outside:
        classification = "IN_DOMAIN"
    elif nearest >= 0.25 and len(outside) <= 1:
        classification = "BORDERLINE"
    else:
        classification = "OUT_OF_DOMAIN"
    return {
        "classification": classification,
        "nearest_training_similarity": round(nearest, 4),
        "chemical_space_distance": round(1.0 - nearest, 4),
        "fingerprint": "Morgan radius 2, 2048 bits, Tanimoto",
        "descriptor_values": {key: round(value, 4) for key, value in values.items()},
        "descriptor_training_ranges": {key: {bound: round(value, 4) for bound, value in limits.items()} for key, limits in ranges.items()},
        "descriptors_outside_range": outside,
        "method": "Heuristic training-set similarity and full descriptor-range envelope; not a model-calibrated uncertainty estimate.",
    }


def _load_model(endpoint: str):
    global _TRAINER
    from chemprop import models
    from lightning import pytorch as pl

    key = MODEL_SPECS[endpoint]["model_key"]
    if key not in _MODELS:
        _MODELS[key] = models.MPNN.load_from_file(MODEL_ROOT / key / "model_v2_1.pt", map_location="cpu")
        _MODELS[key].eval()
    if _TRAINER is None:
        _TRAINER = pl.Trainer(
            logger=False, enable_checkpointing=False, enable_model_summary=False,
            enable_progress_bar=False, accelerator="cpu", devices=1,
        )
    return _MODELS[key], _TRAINER


def predict_endpoint(smiles: str, endpoint: str) -> dict:
    available, reason = model_files_available(endpoint)
    if not available:
        return {"status": "MODEL_UNAVAILABLE", "reason": reason}
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")
    value = predict_batch_values([Chem.MolToSmiles(mol, isomericSmiles=True)], endpoint)[0]
    domain = applicability_domain(smiles, endpoint)
    # A single deterministic checkpoint provides no ensemble uncertainty, so HIGH is never assigned.
    confidence = "MEDIUM" if domain["classification"] == "IN_DOMAIN" else "LOW"
    return {
        "status": "COMPLETE",
        "predicted_value": value,
        "unit": MODEL_SPECS[endpoint]["unit"],
        "confidence": confidence,
        "applicability_domain": domain,
        "uncertainty": None,
        "uncertainty_reason": "Single checkpoint; no calibrated ensemble uncertainty is available.",
    }


def predict_batch_values(smiles: list[str], endpoint: str) -> list[float]:
    """Run one CPU inference batch; used by endpoint inference and reproducible validation."""
    from chemprop import data, featurizers
    import torch

    datapoints = [data.MoleculeDatapoint.from_smi(value) for value in smiles]
    dataset = data.MoleculeDataset(datapoints, featurizer=featurizers.SimpleMoleculeMolGraphFeaturizer())
    loader = data.build_dataloader(dataset, shuffle=False, num_workers=0, drop_last=False)
    with _MODEL_LOCK, torch.inference_mode():
        model, trainer = _load_model(endpoint)
        batches = trainer.predict(model, loader)
    return [float(value) for batch in batches for value in batch.reshape(-1)]


def _normalise_unit(unit: str) -> str:
    return "".join(str(unit).lower().replace("μ", "µ").split())


def comparable_experimental(endpoint: str, measurement, endpoint_name: str) -> tuple[float | None, str]:
    """Convert compatible experimental values to the model's logarithmic output scale."""
    name = endpoint_name.lower()
    matrix_method = f"{measurement.matrix} {measurement.method}".lower()
    value = measurement.mean_value if measurement.mean_value is not None else measurement.value
    if value is None or measurement.qualifier not in ("", "=", "~"):
        return None, "Missing/qualified experimental value"
    value = float(value)
    unit = _normalise_unit(measurement.unit)
    if endpoint == "Solubility":
        if "intrinsic" in name or "ph" in name or "intrinsic" in matrix_method or "ph" in matrix_method:
            return None, "pH-specific/intrinsic solubility is not comparable to aggregate aqueous LogS"
        if "solub" not in name:
            return None, "Different endpoint"
        if unit in {"log10(mol/l)", "logs", "logmol/l", "log(mol/l)"}:
            return value, "Experimental LogS"
        scales = {"mol/l": 1.0, "m": 1.0, "mmol/l": 1e-3, "mm": 1e-3, "µmol/l": 1e-6, "µm": 1e-6, "um": 1e-6}
        if unit in scales and value > 0:
            return math.log10(value * scales[unit]), f"Converted from {measurement.unit} to log10(mol/L)"
        return None, "Incompatible solubility unit"
    if "pampa" in name or "mdck" in name or "pampa" in matrix_method or "mdck" in matrix_method:
        return None, "PAMPA/MDCK is not comparable to Caco-2"
    if "caco" not in name and "caco" not in matrix_method:
        return None, "Experimental permeability is not identified as Caco-2"
    if unit in {"log10(cm/s)", "logpapp", "log(cm/s)"}:
        return value, "Experimental Caco-2 LogPapp"
    if unit in {"cm/s", "cmsec-1", "cms-1"} and value > 0:
        return math.log10(value), "Converted from cm/s to log10(cm/s)"
    if unit in {"10^-6cm/s", "10−6cm/s", "10-6cm/s"} and value > 0:
        return math.log10(value * 1e-6), f"Converted from {measurement.unit} to log10(cm/s)"
    if unit == "µm/s" and value > 0:
        return math.log10(value * 1e-4), "Converted from µm/s to log10(cm/s)"
    return None, "Incompatible Caco-2 unit"


def comparison_for_prediction(endpoint: str, predicted_value: float, measurements, endpoint_names: dict[int, str]) -> list[dict]:
    comparisons = []
    for measurement in measurements:
        experimental_log, note = comparable_experimental(endpoint, measurement, endpoint_names[measurement.endpoint_id])
        if experimental_log is None:
            continue
        predicted_linear, experimental_linear = 10 ** predicted_value, 10 ** experimental_log
        comparisons.append({
            "measurement_id": measurement.id,
            "experimental_value": measurement.mean_value if measurement.mean_value is not None else measurement.value,
            "experimental_unit": measurement.unit,
            "experimental_normalized": round(experimental_log, 6),
            "predicted_normalized": round(predicted_value, 6),
            "normalized_unit": MODEL_SPECS[endpoint]["unit"],
            "absolute_error": round(abs(predicted_value - experimental_log), 6),
            "relative_error_percent_linear_scale": round(abs(predicted_linear - experimental_linear) / abs(experimental_linear) * 100, 3),
            "conversion": note,
        })
    return comparisons
