"""Endpoint-specific ADMET inference and scientific guardrails (Stages 3A-3C)."""

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
OPENADMET_ROOT = ROOT / "models" / "openadmet" / "microsomal_clearance"
MODEL_VERSION = "admetica-d4f7056-chemprop-v2.1"
CYP_MODEL_VERSION = "admetica-d4f7056-cyp-chemprop-v2.1"

MODEL_SPECS = {
    "Solubility": {
        "model_key": "solubility",
        "model_family": "admetica",
        "model_version": "admetica-d4f7056-chemprop-v2.1",
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
        "model_family": "admetica",
        "model_version": "admetica-d4f7056-chemprop-v2.1",
        "display_name": "Admetica Chemprop Caco-2",
        "endpoint_definition": "Caco-2 apparent permeability LogPapp = log10(Papp [cm/s]). The aggregated training file does not retain A→B/B→A direction or detailed assay conditions.",
        "unit": "log10(cm/s)",
        "training_dataset": "Wang et al. Caco-2 compiled permeability set as distributed by Admetica (910 rows)",
        "validation": {"MAE": 0.317, "RMSE": 0.415, "R2": 0.701, "Spearman": 0.832, "external_MAE_34_compounds": 0.412},
        "source": "https://github.com/datagrok-ai/admetica",
        "license": "MIT",
        "limitations": "Assay direction and conditions are absent from the aggregate source. The value must not be treated as PAMPA, MDCK, or an efflux ratio.",
    },
    "Plasma protein binding": {
        "model_key": "ppbr",
        "model_family": "admetica",
        "model_version": "admetica-d4f7056-chemprop-v2.1",
        "display_name": "Admetica Chemprop Human PPB",
        "endpoint_definition": "Human plasma protein binding rate: percent of compound bound in human plasma. Derived fu is (100 - percent bound) / 100 and is not a separate model output.",
        "unit": "% bound",
        "species": "Human",
        "training_dataset": "AstraZeneca/ChEMBL CHEMBL3301361 human plasma protein binding rate (2,790 rows distributed by Admetica)",
        "validation": {"MAE_percent_bound": 6.919, "RMSE_percent_bound": 11.294, "R2": 0.609, "Spearman": 0.762},
        "independent_validation": {"dataset": "Biogen prospective public ADME, canonical training overlap excluded", "n": 185, "MAE_percent_bound": 14.6194, "RMSE_percent_bound": 21.795, "R2": 0.4389, "Spearman": 0.6105},
        "source": "https://github.com/datagrok-ai/admetica",
        "license": "MIT",
        "limitations": "Single deterministic checkpoint. Assay-level conditions are not retained in the model input. Predictions outside the physical 0-100% interval are not clipped or converted to fu.",
    },
}

# Each CYP classifier is an endpoint-specific checkpoint. Inhibition and
# substrate status are intentionally distinct endpoints and never interchanged.
_CYP_INHIBITOR_METRICS = {
    "CYP1A2": (13239, 0.873, 0.866, 0.870, 0.869),
    "CYP2C9": (12881, 0.830, 0.819, 0.826, 0.824),
    "CYP2C19": (13427, 0.819, 0.830, 0.824, 0.825),
    # The upstream overview's model-size cell says 11,127 while its dataset
    # table/released CSV contains about 13.9k records. Preserve both facts.
    "CYP2D6": (13898, 0.866, 0.751, 0.843, 0.808),
    "CYP3A4": (12997, 0.815, 0.842, 0.826, 0.829),
}
_CYP_SUBSTRATE_METRICS = {
    "CYP2C9": (899, 0.728, 0.757, 0.738, 0.742),
    "CYP2D6": (941, 0.749, 0.769, 0.753, 0.759),
    "CYP3A4": (1149, 0.569, 0.779, 0.718, 0.674),
}

for _isoform, (_count, _specificity, _sensitivity, _accuracy, _balanced) in _CYP_INHIBITOR_METRICS.items():
    _slug = f"{_isoform.lower()}-inhibitor"
    MODEL_SPECS[f"{_isoform} inhibitor"] = {
        "model_key": f"cyp/{_slug}",
        "model_family": "admetica",
        "prediction_type": "binary_classification",
        "model_version": CYP_MODEL_VERSION,
        "display_name": f"Admetica Chemprop {_isoform} inhibitor",
        "isoform": _isoform,
        "role": "INHIBITOR",
        "decision_threshold": 0.5,
        "endpoint_definition": (
            f"Binary {_isoform} functional inhibition at the PubChem AID 1851 activity threshold: "
            "active when AC50 <= 10 µM in the human CYP pro-luciferin dealkylation assay. "
            "The output is a class probability, not IC50 or AC50."
        ),
        "assay_definition": (
            "Human recombinant CYP pro-luciferin substrate at Km; compounds tested from 40 µM "
            "to 0.24 nM with NADPH, 60-minute enzyme incubation and luminescence readout."
        ),
        "unit": "probability",
        "training_dataset": f"PubChem AID 1851 CYP panel as curated by Admetica ({_count:,} reported records)",
        "training_n": _count,
        "validation": {
            "specificity": _specificity, "sensitivity": _sensitivity, "accuracy": _accuracy,
            "balanced_accuracy": _balanced, "scope": "Admetica publisher-reported validation",
        },
        "source": "https://github.com/datagrok-ai/admetica",
        "license": "MIT for the Admetica repository/checkpoint; upstream PubChem dataset terms remain source-specific",
        "limitations": (
            "Single deterministic checkpoint; probability calibration was not reported. The functional "
            "luminescence assay may also be reduced by substrate turnover or assay interference. No IC50 is inferred."
        ),
    }

MODEL_SPECS["CYP2C9 inhibitor"]["validation"]["publisher_documentation_note"] = (
    "The Admetica root overview reports balanced accuracy 0.824; its metabolism page reports 0.890. "
    "The overview value is retained conservatively."
)
MODEL_SPECS["CYP2D6 inhibitor"]["validation"]["publisher_documentation_note"] = (
    "The upstream dataset table/released data is approximately 13.9k records, while the model metrics table's size cell says 11,127."
)

for _isoform, (_count, _specificity, _sensitivity, _accuracy, _balanced) in _CYP_SUBSTRATE_METRICS.items():
    _slug = f"{_isoform.lower()}-substrate"
    _sources = (
        "Carbon-Mangels et al. and Zaretzki et al. substrate compilations"
        if _isoform in {"CYP2D6", "CYP3A4"} else
        "Carbon-Mangels et al. CYP2C9 substrate classification compilation"
    )
    MODEL_SPECS[f"{_isoform} substrate"] = {
        "model_key": f"cyp/{_slug}",
        "model_family": "admetica",
        "prediction_type": "binary_classification",
        "model_version": CYP_MODEL_VERSION,
        "display_name": f"Admetica Chemprop {_isoform} substrate",
        "isoform": _isoform,
        "role": "SUBSTRATE",
        "decision_threshold": 0.5,
        "endpoint_definition": f"Binary molecular substrate classification for {_isoform}; probability is not a turnover rate, Km, Vmax, or intrinsic clearance.",
        "assay_definition": "Literature-aggregated substrate/non-substrate labels; heterogeneous source assays without a single harmonized concentration or kinetic protocol.",
        "unit": "probability",
        "training_dataset": f"{_sources} ({_count:,} reported records)",
        "training_n": _count,
        "validation": {
            "specificity": _specificity, "sensitivity": _sensitivity, "accuracy": _accuracy,
            "balanced_accuracy": _balanced, "scope": "Admetica publisher-reported validation",
        },
        "source": "https://github.com/datagrok-ai/admetica",
        "license": "MIT for the Admetica repository/checkpoint; upstream literature dataset terms remain source-specific",
        "limitations": "Single deterministic checkpoint with heterogeneous literature labels and no reported probability calibration; no kinetic quantity is inferred.",
    }

_CYP_INDEPENDENT_INHIBITOR = {
    "CYP2C9": {"n": 464, "both_classes": True, "AUROC": 0.5851, "AUPRC": 0.3813, "balanced_accuracy": 0.5717, "sensitivity": 0.5111, "specificity": 0.6322, "MCC": 0.1324},
    "CYP2D6": {"n": 639, "both_classes": True, "AUROC": 0.5999, "AUPRC": 0.4164, "balanced_accuracy": 0.5657, "sensitivity": 0.5436, "specificity": 0.5878, "MCC": 0.1216},
    "CYP3A4": {"n": 788, "both_classes": True, "AUROC": 0.6533, "AUPRC": 0.4471, "balanced_accuracy": 0.6096, "sensitivity": 0.6527, "specificity": 0.5665, "MCC": 0.2015},
}
for _isoform, _metrics in _CYP_INDEPENDENT_INHIBITOR.items():
    MODEL_SPECS[f"{_isoform} inhibitor"]["independent_validation"] = {
        "dataset": "ChEMBL 30 inhibitor set; canonical training overlap excluded (0 exact overlaps)",
        "decision_threshold": 0.5,
        **_metrics,
    }
MODEL_SPECS["CYP3A4 substrate"]["independent_validation"] = {
    "dataset": "FDA-approved tyrosine kinase inhibitors; 2 canonical training overlaps excluded",
    "n": 22, "both_classes": False, "positive_only": True, "sensitivity": 0.9545,
    "note": "Directionality-only sanity set; AUROC/AUPRC/specificity/balanced accuracy/MCC cannot be calculated from positives alone.",
}

for _species, _task, _count in (("HLM", 0, 5086), ("RLM", 1, 670), ("MLM", 2, 5086)):
    MODEL_SPECS[f"{_species} intrinsic clearance"] = {
        "model_key": "microsomal_clearance",
        "index_key": _species.lower(),
        "model_family": "openadmet_clearance",
        "task_index": _task,
        "model_version": "openadmet-microsomal-clearance-chemeleon-v1-e135493",
        "display_name": f"OpenADMET CheMeleon {_species} intrinsic clearance",
        "endpoint_definition": f"{_species}: species-specific liver microsomal intrinsic clearance scaled to in-vivo clearance and expressed as log10(mL/min/kg). This is not raw µL/min/mg protein or a microsomal half-life.",
        "unit": "log10(mL/min/kg)",
        "species": {"HLM": "Human", "RLM": "Rat", "MLM": "Mouse"}[_species],
        "matrix": f"{_species} liver microsomes",
        "training_dataset": f"OpenADMET curated ChEMBL 35, ASAP-Polaris and ExpansionRx microsomal clearance; {_count:,} non-missing {_species} training labels in the released checkpoint",
        "validation": {"released_checkpoint": "No numeric held-out metric published for this exact all-data checkpoint; model card plots compare an analogous checkpoint that excluded ExpansionRx test data."},
        "source": "https://huggingface.co/openadmet/microsomal-clearance-chemeleon-v1",
        "license": "Apache-2.0",
        "limitations": "Released checkpoint was trained on all available ExpansionRx train and test data, has no ensemble uncertainty, and predicts scaled clearance only. It must not be compared directly with raw µL/min/mg values without explicit scaling parameters.",
    }

MODEL_SPECS["HLM intrinsic clearance"]["independent_validation"] = {
    "dataset": "Biogen prospective public ADME, canonical training overlap excluded", "n": 3078,
    "MAE": 0.6259, "RMSE": 0.7616, "R2": -0.4911, "Spearman": 0.3700,
}
MODEL_SPECS["RLM intrinsic clearance"]["independent_validation"] = {
    "dataset": "Biogen prospective public ADME, canonical training overlap excluded", "n": 3045,
    "MAE": 0.6263, "RMSE": 0.7716, "R2": -0.0577, "Spearman": 0.4248,
}
MODEL_SPECS["MLM intrinsic clearance"]["independent_validation"] = {
    "status": "No compatible independent MLM endpoint in the selected Biogen prospective set"
}

# Operational research summaries derived from the released model's training-label
# quartiles. They are deliberately not presented as universal biological cutoffs.
MICROSOMAL_THRESHOLDS = {
    "HLM intrinsic clearance": {"stable_max": 0.903090, "unstable_min": 1.741998},
    "RLM intrinsic clearance": {"stable_max": 1.301030, "unstable_min": 2.171360},
    "MLM intrinsic clearance": {"stable_max": 1.930057, "unstable_min": 2.818209},
}
for _endpoint, _limits in MICROSOMAL_THRESHOLDS.items():
    MODEL_SPECS[_endpoint]["assessment_thresholds"] = {
        "stable": f"≤ {_limits['stable_max']} log10(mL/min/kg)",
        "moderate": f"> {_limits['stable_max']} and < {_limits['unstable_min']} log10(mL/min/kg)",
        "unstable": f"≥ {_limits['unstable_min']} log10(mL/min/kg)",
        "basis": "25th and 75th percentiles of the species-specific labels used by the released checkpoint; an operational project summary, not a universal assay standard.",
    }


def metabolic_stability_assessment(endpoint: str, value: float) -> dict | None:
    limits = MICROSOMAL_THRESHOLDS.get(endpoint)
    if not limits:
        return None
    if value <= limits["stable_max"]:
        category = "STABLE"
    elif value >= limits["unstable_min"]:
        category = "UNSTABLE"
    else:
        category = "MODERATE"
    return {
        "category": category,
        "metabolic_liability_flag": "METABOLIC STABILITY CONCERN" if category == "UNSTABLE" else None,
        "thresholds": MODEL_SPECS[endpoint]["assessment_thresholds"],
        "evidence_endpoint": endpoint,
        "evidence_value": value,
        "evidence_unit": MODEL_SPECS[endpoint]["unit"],
    }

_MODEL_LOCK = threading.Lock()
_MODELS: dict[str, object] = {}
_TRAINER = None
_FP_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
_DESCRIPTOR_NAMES = ("MW", "cLogP", "TPSA", "HBD", "HBA", "RotB")


def registry_seed(endpoint: str) -> dict:
    spec = MODEL_SPECS[endpoint]
    if spec.get("role"):
        supported_matrix = ["human recombinant CYP enzyme"]
    elif spec.get("matrix"):
        supported_matrix = [spec["matrix"]]
    elif endpoint == "Permeability":
        supported_matrix = ["Caco-2"]
    elif endpoint == "Plasma protein binding":
        supported_matrix = ["human plasma"]
    else:
        supported_matrix = ["aqueous"]
    return {
        "endpoint_name": endpoint,
        "model_name": spec["display_name"],
        "model_version": spec["model_version"],
        "implementation_status": "READY",
        "supported_species": [spec["species"]] if spec.get("species") else (["Human"] if spec.get("role") else (["in vitro Caco-2"] if endpoint == "Permeability" else [])),
        "supported_matrix": supported_matrix,
        "output_unit": spec["unit"],
        "provenance_json": {key: value for key, value in spec.items() if key not in {"model_key", "display_name", "unit", "task_index"}},
        "is_active": True,
    }


def model_files_available(endpoint: str) -> tuple[bool, str]:
    if endpoint not in MODEL_SPECS:
        return False, "No endpoint- and species-specific model installed; cross-species reuse is prohibited."
    spec = MODEL_SPECS[endpoint]
    if spec["model_family"] == "openadmet_clearance":
        names = ("model.pth", "X_train.csv", "y_train.csv", f"ad_index_{spec['index_key']}.npz")
        missing = [name for name in names if not (OPENADMET_ROOT / name).is_file()]
    else:
        key = spec["model_key"]
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


@lru_cache(maxsize=32)
def _training_index(endpoint: str):
    spec = MODEL_SPECS[endpoint]
    path = (OPENADMET_ROOT / f"ad_index_{spec['index_key']}.npz") if spec["model_family"] == "openadmet_clearance" else (MODEL_ROOT / spec["model_key"] / "ad_index.npz")
    with np.load(path) as index:
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
    training_fps, training_bit_counts, ranges = _training_index(endpoint)
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

    spec = MODEL_SPECS[endpoint]
    key = spec["model_key"]
    if key not in _MODELS:
        if spec["model_family"] == "openadmet_clearance":
            import torch
            from chemprop import nn
            model = models.MPNN(
                nn.BondMessagePassing(d_h=2048, depth=3, dropout=0.25), nn.MeanAggregation(),
                nn.RegressionFFN(n_tasks=3, input_dim=2048, hidden_dim=512, n_layers=3, dropout=0.25,
                                 output_transform=nn.UnscaleTransform([0, 0, 0], [1, 1, 1])),
                batch_norm=False, metrics=[nn.metrics.MSE(), nn.metrics.MAE(), nn.metrics.RMSE()],
            )
            state = torch.load(OPENADMET_ROOT / "model.pth", map_location="cpu", weights_only=True)
            model.load_state_dict(state, strict=True)
            _MODELS[key] = model
        else:
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
    # Clearance generalization was weak on the canonical-overlap-excluded Biogen set;
    # domain membership therefore cannot elevate it above LOW.
    if MODEL_SPECS[endpoint].get("prediction_type") == "binary_classification":
        # Confidence is deliberately independent of the predicted probability.
        # It combines publisher performance, independent evidence availability,
        # and the compound's applicability domain.
        balanced = float(MODEL_SPECS[endpoint]["validation"]["balanced_accuracy"])
        independent = MODEL_SPECS[endpoint].get("independent_validation", {})
        independently_supported = (
            independent.get("n", 0) >= 30 and independent.get("both_classes", False)
            and float(independent.get("balanced_accuracy", 0.0)) >= 0.70
        )
        if domain["classification"] == "IN_DOMAIN" and balanced >= 0.80 and independently_supported:
            confidence = "MEDIUM"
        else:
            confidence = "LOW"
    else:
        confidence = "LOW" if endpoint.endswith("intrinsic clearance") else ("MEDIUM" if domain["classification"] == "IN_DOMAIN" else "LOW")
    result = {
        "status": "COMPLETE",
        "predicted_value": value,
        "unit": MODEL_SPECS[endpoint]["unit"],
        "confidence": confidence,
        "applicability_domain": domain,
        "uncertainty": None,
        "uncertainty_reason": "Single checkpoint; no calibrated ensemble uncertainty is available.",
    }
    if endpoint == "Plasma protein binding":
        result["derived_outputs"] = ({
            "fraction_bound": value / 100.0,
            "fu_fraction": (100.0 - value) / 100.0,
            "fu_percent": 100.0 - value,
            "derivation": "fu = 1 - fraction bound",
        } if 0.0 <= value <= 100.0 else {
            "derivation": "Not calculated because the raw model output is outside the physical 0-100% bound range.",
        })
        if not 0.0 <= value <= 100.0:
            result["confidence"] = "LOW"
    elif endpoint.endswith("intrinsic clearance"):
        result["derived_outputs"] = {
            "scaled_clint_mL_min_kg": 10 ** value,
            "transformation": "10 ** predicted log10(mL/min/kg)",
        }
        result["metabolic_stability_assessment"] = metabolic_stability_assessment(endpoint, value)
    elif MODEL_SPECS[endpoint].get("prediction_type") == "binary_classification":
        spec = MODEL_SPECS[endpoint]
        positive = value >= spec["decision_threshold"]
        positive_label = spec["role"]
        result.update({
            "probability": value,
            "classification": positive_label if positive else f"NON_{positive_label}",
            "isoform": spec["isoform"],
            "role": spec["role"],
            "decision_threshold": spec["decision_threshold"],
            "liability_summary": ({
                "flag": f"Potential {spec['isoform']} inhibition concern",
                "rule": f"inhibitor probability >= {spec['decision_threshold']:.2f}",
                "basis": "The model's fixed binary decision threshold; this is a screening flag, not an IC50 claim.",
            } if spec["role"] == "INHIBITOR" and positive else None),
        })
    return result


def predict_batch_values(smiles: list[str], endpoint: str) -> list[float]:
    """Run one CPU inference batch; used by endpoint inference and reproducible validation."""
    matrix = predict_batch_matrix(smiles, endpoint)
    task = MODEL_SPECS[endpoint].get("task_index")
    return [float(row[task] if task is not None else row[0]) for row in matrix]


def predict_batch_matrix(smiles: list[str], endpoint: str) -> list[list[float]]:
    """Run one CPU inference batch and retain all model tasks."""
    from chemprop import data, featurizers
    import torch

    datapoints = [data.MoleculeDatapoint.from_smi(value) for value in smiles]
    dataset = data.MoleculeDataset(datapoints, featurizer=featurizers.SimpleMoleculeMolGraphFeaturizer())
    loader = data.build_dataloader(dataset, shuffle=False, num_workers=0, drop_last=False)
    with _MODEL_LOCK, torch.inference_mode():
        model, trainer = _load_model(endpoint)
        batches = trainer.predict(model, loader)
    return [[float(value) for value in row] for batch in batches for row in batch.reshape(-1, batch.shape[-1] if batch.ndim > 1 else 1)]


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
    if MODEL_SPECS[endpoint].get("prediction_type") == "binary_classification":
        spec = MODEL_SPECS[endpoint]
        identity = f"{endpoint_name} {measurement.matrix} {measurement.method} {measurement.notes}".lower()
        isoform = spec["isoform"].lower()
        role_terms = ("inhibitor", "inhibition") if spec["role"] == "INHIBITOR" else ("substrate",)
        opposite_terms = ("substrate",) if spec["role"] == "INHIBITOR" else ("inhibitor", "inhibition")
        if isoform not in identity or not any(term in identity for term in role_terms):
            return None, "Different CYP isoform or role"
        if any(term in identity for term in opposite_terms):
            return None, "CYP inhibitor and substrate evidence cannot be interchanged"
        if any(term in identity for term in ("rat", "mouse", "dog", "monkey")):
            return None, "Non-human CYP evidence is not comparable to the human model"
        if unit in {"class", "binary", "classification", "0/1"} and value in {0.0, 1.0}:
            return value, "Experimental binary classification"
        return None, "Quantitative CYP evidence is retained but not numerically compared with a classification prediction"
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
    if endpoint == "Plasma protein binding":
        identity = f"{endpoint_name} {measurement.species} {measurement.matrix} {measurement.method}".lower()
        if "rat" in identity or "mouse" in identity or "dog" in identity or "monkey" in identity:
            return None, "Species-specific PPB is not human PPB"
        if "ppb" not in identity and "protein binding" not in identity and "plasma binding" not in identity:
            return None, "Experimental endpoint is not identified as plasma protein binding"
        if unit in {"%bound", "%bound", "percentbound", "%proteinbound"}:
            return (value, "Experimental human percent bound") if 0 <= value <= 100 else (None, "Percent bound outside physical 0-100 range")
        if unit in {"fractionbound", "fraction_bound"}:
            return (value * 100.0, "Converted fraction bound to percent bound") if 0 <= value <= 1 else (None, "Fraction bound outside physical 0-1 range")
        if unit in {"fu", "fractionunbound", "fraction_unbound"}:
            return ((1.0 - value) * 100.0, "Converted fu fraction using percent bound = (1 - fu) * 100") if 0 <= value <= 1 else (None, "fu outside physical 0-1 range")
        if unit in {"%unbound", "percentunbound", "fu%"}:
            return (100.0 - value, "Converted percent unbound using percent bound = 100 - percent unbound") if 0 <= value <= 100 else (None, "Percent unbound outside physical 0-100 range")
        return None, "Incompatible PPB unit; percent bound and fu must be explicit"
    if endpoint.endswith("intrinsic clearance"):
        species = MODEL_SPECS[endpoint]["species"].lower()
        code = endpoint.split()[0].lower()
        identity = f"{endpoint_name} {measurement.species} {measurement.matrix} {measurement.method}".lower()
        species_terms = {"human": ("human", "hlm"), "rat": ("rat", "rlm"), "mouse": ("mouse", "mlm")}[species]
        if not any(term in identity for term in species_terms):
            return None, f"Experimental clearance is not identified as {code.upper()}"
        if any(term in identity for other, terms in {"human": ("human", "hlm"), "rat": ("rat", "rlm"), "mouse": ("mouse", "mlm")}.items() if other != species for term in terms):
            return None, "Microsomal species mismatch"
        if unit in {"log10(ml/min/kg)", "logclint", "log10ml/min/kg"}:
            return value, f"Experimental {code.upper()} scaled log intrinsic clearance"
        if unit in {"ml/min/kg", "mlmin-1kg-1"} and value > 0:
            return math.log10(value), "Converted scaled mL/min/kg to log10(mL/min/kg)"
        if unit in {"µl/min/mg", "ul/min/mg", "ml/min/g", "µlmin-1mg-1", "mlmin-1g-1"} and value > 0:
            provenance = measurement.provenance_json or {}
            mppgl = provenance.get("microsomal_protein_mg_per_g_liver")
            liver = provenance.get("liver_weight_g_per_kg")
            if mppgl is None or liver is None:
                return None, "Raw microsomal clearance requires explicit MPPGL and liver-weight scaling parameters"
            scaled = value * float(mppgl) * float(liver) / 1000.0
            if scaled <= 0:
                return None, "Invalid microsomal scaling parameters"
            return math.log10(scaled), f"Scaled raw clearance using Clint * {mppgl} mg/g * {liver} g/kg / 1000"
        return None, "Incompatible microsomal clearance unit"
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
        if MODEL_SPECS[endpoint].get("prediction_type") == "binary_classification":
            predicted_class = 1.0 if predicted_value >= MODEL_SPECS[endpoint]["decision_threshold"] else 0.0
            comparisons.append({
                "measurement_id": measurement.id,
                "experimental_value": measurement.mean_value if measurement.mean_value is not None else measurement.value,
                "experimental_unit": measurement.unit,
                "experimental_normalized": experimental_log,
                "predicted_normalized": predicted_class,
                "normalized_unit": "class",
                "absolute_error": None,
                "relative_error_percent_linear_scale": None,
                "classification_match": predicted_class == experimental_log,
                "conversion": note,
            })
            continue
        if endpoint == "Plasma protein binding":
            predicted_linear, experimental_linear = predicted_value, experimental_log
        else:
            predicted_linear, experimental_linear = 10 ** predicted_value, 10 ** experimental_log
        relative_error = (abs(predicted_linear - experimental_linear) / abs(experimental_linear) * 100) if experimental_linear != 0 else None
        comparisons.append({
            "measurement_id": measurement.id,
            "experimental_value": measurement.mean_value if measurement.mean_value is not None else measurement.value,
            "experimental_unit": measurement.unit,
            "experimental_normalized": round(experimental_log, 6),
            "predicted_normalized": round(predicted_value, 6),
            "normalized_unit": MODEL_SPECS[endpoint]["unit"],
            "absolute_error": round(abs(predicted_value - experimental_log), 6),
            "relative_error_percent_linear_scale": round(relative_error, 3) if relative_error is not None else None,
            "conversion": note,
        })
    return comparisons


def cyp_experimental_evidence(endpoint: str, predicted_probability: float, measurements, endpoint_names: dict[int, str]) -> list[dict]:
    """Return matching CYP evidence without inventing cross-type numerical errors."""
    spec = MODEL_SPECS.get(endpoint, {})
    if spec.get("prediction_type") != "binary_classification":
        return []
    result = []
    isoform = spec["isoform"].lower()
    role_terms = ("inhibitor", "inhibition") if spec["role"] == "INHIBITOR" else ("substrate",)
    opposite_terms = ("substrate",) if spec["role"] == "INHIBITOR" else ("inhibitor", "inhibition")
    predicted_class = 1 if predicted_probability >= spec["decision_threshold"] else 0
    for measurement in measurements:
        endpoint_name = endpoint_names.get(measurement.endpoint_id, "")
        identity = f"{endpoint_name} {measurement.matrix} {measurement.method} {measurement.notes}".lower()
        if isoform not in identity or not any(term in identity for term in role_terms):
            continue
        if any(term in identity for term in opposite_terms) or any(term in identity for term in ("rat", "mouse", "dog", "monkey")):
            continue
        value = measurement.mean_value if measurement.mean_value is not None else measurement.value
        if value is None:
            continue
        unit = _normalise_unit(measurement.unit)
        classification = unit in {"class", "binary", "classification", "0/1"} and float(value) in {0.0, 1.0}
        result.append({
            "measurement_id": measurement.id,
            "endpoint": endpoint_name,
            "value": value,
            "unit": measurement.unit,
            "evidence_type": "CLASSIFICATION" if classification else "QUANTITATIVE",
            "comparison": ("AGREES" if predicted_class == int(float(value)) else "DISAGREES") if classification else "NOT_NUMERICALLY_COMPARABLE",
            "comparison_note": (
                "Binary experimental class compared with the model class at the documented threshold."
                if classification else
                "Quantitative experimental CYP evidence (for example IC50) is displayed separately; no error is calculated against a classification probability."
            ),
            "absolute_error": None,
            "relative_error_percent": None,
        })
    return result
