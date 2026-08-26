"""Endpoint-specific ADMET inference and scientific guardrails (Stages 3A-3F)."""

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
ADMET_AI_ROOT = ROOT / "models" / "admet_ai"
MODEL_VERSION = "admetica-d4f7056-chemprop-v2.1"
CYP_MODEL_VERSION = "admetica-d4f7056-cyp-chemprop-v2.1"
TRANSPORTER_MODEL_VERSION = "admetica-d4f7056-pgp-inhibitor-chemprop-v2.1"
SAFETY_HERG_MODEL_VERSION = "admetica-d4f7056-herg-chemprop-v2.1"
ADMET_AI_SAFETY_MODEL_VERSION = "admet-ai-v2.0.1-c65bf04-chemprop-v2-ensemble5"

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

# The only transporter checkpoint qualified for local activation in Stage 3E.
# P-gp substrate and every other requested transporter/role remain explicit
# MODEL_UNAVAILABLE registry entries: endpoint family membership is not a
# scientific basis for reusing this inhibitor model.
MODEL_SPECS["P-gp inhibitor"] = {
    "model_key": "transporter/pgp-inhibitor",
    "model_family": "admetica",
    "prediction_type": "binary_classification",
    "model_version": TRANSPORTER_MODEL_VERSION,
    "display_name": "Admetica Chemprop human P-gp/ABCB1 inhibitor",
    "transporter": "P-gp / ABCB1",
    "role": "INHIBITOR",
    "species": "Human",
    "decision_threshold": 0.5,
    "endpoint_definition": (
        "Binary human P-glycoprotein (P-gp/ABCB1) functional inhibitor classification. "
        "The released model output is a binary-model probability score, not Ki, IC50, "
        "efflux ratio, or substrate status."
    ),
    "assay_definition": (
        "Broccatelli literature aggregation from more than 60 sources with heterogeneous "
        "human P-gp functional assays, cell systems, probes, and conditions. The source "
        "assigned inhibitor labels at IC50 <= 15 µM or >25-30% inhibition and non-inhibitor "
        "labels at IC50 >= 100 µM or <10-12% inhibition where such evidence was reported."
    ),
    "unit": "probability",
    "training_dataset": (
        "Broccatelli et al. human P-gp inhibitor compilation as curated by Admetica "
        "(1,275 reported compounds; 666 inhibitors, 609 non-inhibitors; 1,227 valid "
        "structures in the packaged curated file)"
    ),
    "training_n": 1275,
    "validation": {
        "specificity": 0.916, "sensitivity": 0.863, "accuracy": 0.888,
        "balanced_accuracy": 0.889, "scope": "Admetica publisher-reported validation",
        "probability_calibration": "Not reported",
    },
    "independent_validation": {
        "status": "NOT_AVAILABLE",
        "reason": (
            "No rigorously independent public structure/label set was qualified: accessible "
            "alternatives share the Broccatelli/Chen source lineage or do not publish reusable structures."
        ),
    },
    "source": "https://github.com/datagrok-ai/admetica",
    "license": (
        "MIT for the Admetica repository/checkpoint; the upstream literature aggregation "
        "retains source-specific terms"
    ),
    "limitations": (
        "Single deterministic checkpoint; probability calibration and rigorous independent "
        "validation are unavailable. Training assays are heterogeneous, and exact training-set "
        "overlap cannot establish prospective performance. Confidence is therefore capped at LOW."
    ),
}

MODEL_SPECS["hERG liability"] = {
    "model_key": "safety/herg",
    "model_family": "admetica",
    "prediction_type": "binary_classification",
    "model_version": SAFETY_HERG_MODEL_VERSION,
    "display_name": "Admetica Chemprop human hERG blocker liability",
    "safety_endpoint": "hERG",
    "species": "Human",
    "positive_label": "BLOCKER",
    "negative_label": "NON_BLOCKER",
    "decision_threshold": 0.5,
    "endpoint_definition": (
        "Binary human hERG/KCNH2 blocker-liability classification. The aggregated source does "
        "not retain a uniform assay mode, so this is neither a pure binding endpoint nor a pure "
        "functional patch-clamp endpoint and is not an IC50 prediction."
    ),
    "assay_definition": (
        "Wang et al. literature aggregation of heterogeneous hERG blocker evidence. Binding and "
        "functional assay provenance is not retained per packaged row; class labels must therefore "
        "be interpreted as screening liability rather than one standardized assay."
    ),
    "unit": "probability",
    "training_dataset": (
        "Wang et al. hERG blocker compilation as curated by Admetica (22,249 reported records; "
        "22,248 valid structures packaged; 19,130 positive and 3,118 negative)"
    ),
    "training_n": 22249,
    "validation": {
        "specificity": 0.811, "sensitivity": 0.897, "accuracy": 0.885,
        "balanced_accuracy": 0.854, "scope": "Admetica publisher-reported validation",
        "probability_calibration": "Not reported",
    },
    "independent_validation": {
        "status": "COMPLETE_WITH_LIMITATIONS",
        "dataset": "OpenADMET ChEMBL 37 human hERG IC50 aggregate; exact training overlap excluded",
        "n": 728,
        "both_classes": True,
        "positive_definition": "median IC50 <= 10,000 nM",
        "AUROC": 0.6668976906,
        "AUPRC": 0.7854091863,
        "balanced_accuracy": 0.5442154170,
        "sensitivity": 0.9754601227,
        "specificity": 0.1129707113,
        "MCC": 0.1844230442,
        "threshold": 0.5,
        "overlap_policy": "7,249 exact canonical-SMILES overlaps removed from 7,977 ChEMBL aggregates",
        "limitations": "Exact-structure exclusion does not prove source, series, or assay-lineage independence.",
    },
    "source": "https://github.com/datagrok-ai/admetica",
    "license": (
        "MIT for the Admetica repository/checkpoint; upstream dataset licensing is source-specific"
    ),
    "limitations": (
        "Heterogeneous blocker labels, severe class imbalance, no reported calibration, and no "
        "assay-mode field. The overlap-filtered ChEMBL check had very low specificity (0.113), so "
        "confidence is capped at LOW. It does not distinguish binding from functional current "
        "inhibition and must not be converted to IC50 or used as a clinical QT-risk determination."
    ),
}

for _endpoint, _task, _species, _positive, _negative, _count, _auprc, _auroc in (
    ("Ames mutagenicity", 0, "Salmonella typhimurium", "MUTAGENIC", "NON_MUTAGENIC", 7255, 0.8957980387, 0.8815869914),
    ("DILI clinical liability", 13, "Human", "DILI_CONCERN", "NO_DILI_CONCERN", 475, 0.8777147000, 0.8814562005),
):
    _key = "ames" if _endpoint.startswith("Ames") else "dili"
    MODEL_SPECS[_endpoint] = {
        "model_key": "classification",
        "index_key": _key,
        "model_family": "admet_ai_ensemble",
        "prediction_type": "binary_classification",
        "task_index": _task,
        "model_version": ADMET_AI_SAFETY_MODEL_VERSION,
        "display_name": f"ADMET-AI v2 Chemprop ensemble {_endpoint}",
        "safety_endpoint": "Ames" if _key == "ames" else "DILI",
        "species": _species,
        "positive_label": _positive,
        "negative_label": _negative,
        "decision_threshold": 0.5,
        "endpoint_definition": (
            "Binary Ames mutagenicity: bacterial reverse-mutation positive versus negative across "
            "an aggregate of four public studies."
            if _key == "ames" else
            "Binary human clinical drug-induced liver injury association from the FDA National "
            "Center for Toxicological Research compilation; not a mechanistic hepatotoxicity assay."
        ),
        "assay_definition": (
            "Aggregated Salmonella bacterial reverse-mutation labels; strain, metabolic activation, "
            "dose, and protocol are not harmonized in the molecular input."
            if _key == "ames" else
            "Clinical drug-level DILI annotation compiled by FDA/NCTR; it is distinct from in-vitro "
            "cytotoxicity, mitochondrial toxicity, and a quantitative liver injury measurement."
        ),
        "unit": "probability",
        "training_dataset": (
            f"TDC AMES/Xu et al. aggregate ({_count:,} valid training labels)"
            if _key == "ames" else
            f"TDC DILI/Xu et al. FDA-NCTR compilation ({_count:,} drugs)"
        ),
        "training_n": _count,
        "validation": {
            "AUROC": _auroc, "AUPRC": _auprc,
            "scope": "ADMET-AI v2 release-reported five-fold held-out evaluation",
            "probability_calibration": "Not reported",
        },
        "independent_validation": {
            "status": "NOT_AVAILABLE",
            "reason": "No rigorously independent public structure/label set with non-overlapping source lineage was qualified in Stage 3F.",
        },
        "source": "https://github.com/swansonk14/admet_ai",
        "license": (
            "MIT for ADMET-AI code/checkpoints; the TDC endpoint page lists the upstream dataset "
            "license as not specified while also linking CC BY 4.0, so redistribution/commercial "
            "dataset use requires separate review"
        ),
        "limitations": (
            "Five-model ensemble disagreement is available, but calibration and independent "
            "validation are not. Source assays/annotations are heterogeneous; confidence is capped "
            "at LOW and the output is classification only. For Ames, the transparent AD source "
            "index has 7,278 raw public rows while the v2 release reports 7,255 valid labels; exact "
            "fold membership is unavailable, so this AD must be treated as an approximation."
        ),
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
_RUNTIME_LOCK = threading.Lock()
_RUNTIME_ERROR: str | None = None
_RUNTIME_CHECKED = False
_MODELS: dict[str, object] = {}
_TRAINER = None
_FP_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
_DESCRIPTOR_NAMES = ("MW", "cLogP", "TPSA", "HBD", "HBA", "RotB")


def registry_seed(endpoint: str) -> dict:
    spec = MODEL_SPECS[endpoint]
    if spec.get("safety_endpoint") == "hERG":
        supported_matrix = ["heterogeneous human hERG assays"]
    elif spec.get("safety_endpoint") == "Ames":
        supported_matrix = ["bacterial reverse mutation"]
    elif spec.get("safety_endpoint") == "DILI":
        supported_matrix = ["clinical annotation"]
    elif spec.get("transporter"):
        supported_matrix = ["human transporter functional assay"]
    elif spec.get("role"):
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
        "source": str(spec.get("source") or ""),
        "training_dataset": str(spec.get("training_dataset") or ""),
        "validation_json": spec.get("validation") or {},
        "license": str(spec.get("license") or ""),
        "model_priority": int(spec.get("model_priority", 100)),
        "ensemble_eligible": True,
        "species": str(spec.get("species") or ("Human" if spec.get("role") else "")),
        "output_type": str(spec.get("prediction_type") or "regression"),
        "provenance_json": {key: value for key, value in spec.items() if key not in {"model_key", "display_name", "unit", "task_index"}},
        "is_active": True,
    }


def model_files_available(endpoint: str) -> tuple[bool, str]:
    global _RUNTIME_CHECKED, _RUNTIME_ERROR
    if endpoint not in MODEL_SPECS:
        return False, "No endpoint- and species-specific model installed; cross-species reuse is prohibited."
    spec = MODEL_SPECS[endpoint]
    if spec["model_family"] == "openadmet_clearance":
        names = ("model.pth", "X_train.csv", "y_train.csv", f"ad_index_{spec['index_key']}.npz")
        missing = [name for name in names if not (OPENADMET_ROOT / name).is_file()]
    elif spec["model_family"] == "admet_ai_ensemble":
        model_paths = [ADMET_AI_ROOT / "classification" / f"model_{index}.pt" for index in range(5)]
        training_root = ADMET_AI_ROOT / "training" / spec["index_key"]
        paths = model_paths + [training_root / "training.csv", training_root / "ad_index.npz"]
        missing = [str(path.relative_to(ADMET_AI_ROOT)) for path in paths if not path.is_file()]
    else:
        key = spec["model_key"]
        missing = [name for name in ("model_v2_1.pt", "training.csv", "ad_index.npz") if not (MODEL_ROOT / key / name).is_file()]
    if missing:
        return False, "Packaged model asset missing: " + ", ".join(missing)
    # FastAPI sync endpoints may call this concurrently on a fresh process. PyTorch's
    # extension import is not safe to observe half-initialized from another worker thread.
    with _RUNTIME_LOCK:
        if not _RUNTIME_CHECKED:
            try:
                import torch  # noqa: F401
                import chemprop  # noqa: F401
            except ImportError as exc:
                _RUNTIME_ERROR = f"{exc.__class__.__name__}: {exc}"
            _RUNTIME_CHECKED = True
    if _RUNTIME_ERROR:
        return False, f"Runtime dependency unavailable: {_RUNTIME_ERROR}"
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
    if spec["model_family"] == "openadmet_clearance":
        path = OPENADMET_ROOT / f"ad_index_{spec['index_key']}.npz"
    elif spec["model_family"] == "admet_ai_ensemble":
        path = ADMET_AI_ROOT / "training" / spec["index_key"] / "ad_index.npz"
    else:
        path = MODEL_ROOT / spec["model_key"] / "ad_index.npz"
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
        elif spec["model_family"] == "admet_ai_ensemble":
            from chemprop.models import load_model
            _MODELS[key] = [
                load_model(ADMET_AI_ROOT / "classification" / f"model_{index}.pt", multicomponent=False)
                for index in range(5)
            ]
            for member in _MODELS[key]:
                member.eval()
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
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    ensemble_values = None
    if MODEL_SPECS[endpoint]["model_family"] == "admet_ai_ensemble":
        member_matrices = predict_batch_member_matrices([canonical], endpoint)
        task = MODEL_SPECS[endpoint]["task_index"]
        ensemble_values = [float(matrix[0][task]) for matrix in member_matrices]
        value = float(np.mean(ensemble_values))
    else:
        value = predict_batch_values([canonical], endpoint)[0]
    domain = applicability_domain(smiles, endpoint)
    # A single deterministic checkpoint provides no ensemble uncertainty, so HIGH is never assigned.
    # Clearance generalization was weak on the canonical-overlap-excluded Biogen set;
    # domain membership therefore cannot elevate it above LOW.
    if MODEL_SPECS[endpoint].get("prediction_type") == "binary_classification":
        # Confidence is deliberately independent of the predicted probability.
        # It combines publisher performance, independent evidence availability,
        # and the compound's applicability domain.
        validation = MODEL_SPECS[endpoint]["validation"]
        publisher_score = float(validation.get("balanced_accuracy", validation.get("AUROC", 0.0)))
        independent = MODEL_SPECS[endpoint].get("independent_validation", {})
        independently_supported = (
            independent.get("n", 0) >= 30 and independent.get("both_classes", False)
            and float(independent.get("balanced_accuracy", 0.0)) >= 0.70
        )
        if domain["classification"] == "IN_DOMAIN" and publisher_score >= 0.80 and independently_supported:
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
        "uncertainty": round(float(np.std(ensemble_values)), 6) if ensemble_values is not None else None,
        "uncertainty_reason": (
            "Standard deviation across five ADMET-AI v2 checkpoints; useful as model disagreement but not calibrated uncertainty."
            if ensemble_values is not None else
            "Single checkpoint; no calibrated ensemble uncertainty is available."
        ),
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
        positive_label = spec.get("positive_label", spec.get("role"))
        negative_label = spec.get("negative_label", f"NON_{positive_label}")
        target = spec.get("isoform") or spec.get("transporter")
        result.update({
            "probability": value,
            "classification": positive_label if positive else negative_label,
            "decision_threshold": spec["decision_threshold"],
            "liability_summary": ({
                "flag": f"Potential {target} inhibition concern",
                "rule": f"inhibitor probability >= {spec['decision_threshold']:.2f}",
                "basis": "The model's fixed binary decision threshold; this is a screening flag, not an IC50 claim.",
            } if spec.get("role") == "INHIBITOR" and positive else None),
        })
        if spec.get("role"):
            result["role"] = spec["role"]
        if spec.get("isoform"):
            result["isoform"] = spec["isoform"]
        if spec.get("transporter"):
            result["transporter"] = spec["transporter"]
            result["species"] = spec["species"]
        if spec.get("safety_endpoint"):
            result["safety_endpoint"] = spec["safety_endpoint"]
            result["species"] = spec["species"]
            result["ensemble_probabilities"] = ensemble_values
            result["liability_summary"] = ({
                "flag": {
                    "hERG": "Potential hERG blocker liability",
                    "Ames": "Potential Ames mutagenicity concern",
                    "DILI": "Potential clinical DILI association concern",
                }[spec["safety_endpoint"]],
                "rule": f"positive-class probability >= {spec['decision_threshold']:.2f}",
                "basis": "Fixed binary screening threshold; not a quantitative toxicity or clinical causality claim.",
            } if positive else None)
    return result


def predict_batch_values(smiles: list[str], endpoint: str) -> list[float]:
    """Run one CPU inference batch; used by endpoint inference and reproducible validation."""
    matrix = predict_batch_matrix(smiles, endpoint)
    task = MODEL_SPECS[endpoint].get("task_index")
    return [float(row[task] if task is not None else row[0]) for row in matrix]


def predict_batch_matrix(smiles: list[str], endpoint: str) -> list[list[float]]:
    """Run one CPU inference batch and retain all model tasks."""
    members = predict_batch_member_matrices(smiles, endpoint)
    return np.mean(np.asarray(members, dtype=float), axis=0).tolist()


def predict_batch_member_matrices(smiles: list[str], endpoint: str) -> list[list[list[float]]]:
    """Return one prediction matrix per checkpoint for ensemble-aware validation."""
    from chemprop import data, featurizers
    import torch

    datapoints = [data.MoleculeDatapoint.from_smi(value) for value in smiles]
    dataset = data.MoleculeDataset(datapoints, featurizer=featurizers.SimpleMoleculeMolGraphFeaturizer())
    loader = data.build_dataloader(dataset, shuffle=False, num_workers=0, drop_last=False)
    with _MODEL_LOCK, torch.inference_mode():
        model, trainer = _load_model(endpoint)
        models = model if isinstance(model, list) else [model]
        member_batches = [trainer.predict(member, loader) for member in models]
    matrices = []
    for batches in member_batches:
        matrices.append([
            [float(value) for value in row]
            for batch in batches
            for row in batch.reshape(-1, batch.shape[-1] if batch.ndim > 1 else 1)
        ])
    return matrices


def _normalise_unit(unit: str) -> str:
    return "".join(str(unit).lower().replace("μ", "µ").split())


def _classification_target_terms(spec: dict) -> tuple[str, ...]:
    if spec.get("transporter") == "P-gp / ABCB1":
        return ("p-gp", "pgp", "p-glycoprotein", "abcb1")
    if spec.get("transporter"):
        return (spec["transporter"].lower(),)
    if spec.get("safety_endpoint") == "hERG":
        return ("herg", "kcnh2")
    if spec.get("safety_endpoint") == "Ames":
        return ("ames", "mutagen")
    if spec.get("safety_endpoint") == "DILI":
        return ("dili", "drug induced liver injury", "drug-induced liver injury")
    return (spec["isoform"].lower(),)


def _classification_role_terms(spec: dict) -> tuple[tuple[str, ...], tuple[str, ...]]:
    if spec.get("role") == "INHIBITOR":
        return ("inhibitor", "inhibition"), ("substrate",)
    if spec.get("role") == "SUBSTRATE":
        return ("substrate",), ("inhibitor", "inhibition")
    return (), ()


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
        identity = f"{endpoint_name} {measurement.species} {measurement.matrix} {measurement.method} {measurement.notes}".lower()
        target_terms = _classification_target_terms(spec)
        role_terms, opposite_terms = _classification_role_terms(spec)
        family = "safety" if spec.get("safety_endpoint") else ("transporter" if spec.get("transporter") else "CYP")
        if not any(term in identity for term in target_terms) or (role_terms and not any(term in identity for term in role_terms)):
            return None, f"Different {family} target or role"
        if any(term in identity for term in opposite_terms):
            return None, f"{family} inhibitor and substrate evidence cannot be interchanged"
        if spec.get("species") == "Human" and any(term in identity for term in ("rat", "mouse", "dog", "monkey")):
            return None, f"Non-human {family} evidence is not comparable to the human model"
        if unit in {"class", "binary", "classification", "0/1"} and value in {0.0, 1.0}:
            return value, "Experimental binary classification"
        return None, f"Quantitative {family} evidence is retained but not numerically compared with a classification prediction"
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


def classification_experimental_evidence(endpoint: str, predicted_probability: float, measurements, endpoint_names: dict[int, str]) -> list[dict]:
    """Return matching classification evidence without cross-role/target numerical errors."""
    spec = MODEL_SPECS.get(endpoint, {})
    if spec.get("prediction_type") != "binary_classification":
        return []
    result = []
    target_terms = _classification_target_terms(spec)
    role_terms, opposite_terms = _classification_role_terms(spec)
    predicted_class = 1 if predicted_probability >= spec["decision_threshold"] else 0
    for measurement in measurements:
        endpoint_name = endpoint_names.get(measurement.endpoint_id, "")
        identity = f"{endpoint_name} {measurement.species} {measurement.matrix} {measurement.method} {measurement.notes}".lower()
        if not any(term in identity for term in target_terms) or (role_terms and not any(term in identity for term in role_terms)):
            continue
        if any(term in identity for term in opposite_terms):
            continue
        if spec.get("species") == "Human" and any(term in identity for term in ("rat", "mouse", "dog", "monkey")):
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
                "Quantitative experimental evidence (for example IC50) is displayed separately; no error is calculated against a classification probability."
            ),
            "absolute_error": None,
            "relative_error_percent": None,
        })
    return result


def cyp_experimental_evidence(endpoint: str, predicted_probability: float, measurements, endpoint_names: dict[int, str]) -> list[dict]:
    """Backward-compatible Stage 3C name for generic binary evidence matching."""
    return classification_experimental_evidence(endpoint, predicted_probability, measurements, endpoint_names)
