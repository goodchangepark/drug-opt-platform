"""Runtime platform metadata used by the researcher-facing Help page."""

from __future__ import annotations

from functools import lru_cache
import importlib.metadata
import platform
from pathlib import Path


APP_VERSION = "0.6.0-stage5b4-stable"

PACKAGE_SPECS = (
    ("Python", None, "Runtime for the backend and scientific services"),
    ("RDKit", "rdkit", "Structure parsing, standardization, properties, fingerprints and alerts"),
    ("SyGMa", "sygma", "Rule-based metabolic soft spots and metabolite hypotheses"),
    ("mordredcommunity", "mordredcommunity", "Molecular descriptor library"),
    ("descriptastorus", "descriptastorus", "Descriptor generation used by model tooling"),
    ("mhfp", "mhfp", "MinHash molecular fingerprints"),
    ("padelpy", "padelpy", "PaDEL descriptor integration"),
    ("Chemprop", "chemprop", "Message-passing neural-network inference"),
    ("PyTorch", "torch", "CPU tensor and neural-network runtime"),
    ("Lightning", "lightning", "Model runtime support"),
    ("scikit-learn", "scikit-learn", "Classical ML and validation utilities"),
    ("SciPy", "scipy", "Scientific numerical methods"),
    ("FastAPI", "fastapi", "Backend HTTP API"),
    ("Uvicorn", "uvicorn", "ASGI production/local server"),
    ("SQLAlchemy", "sqlalchemy", "Database persistence and project isolation"),
    ("Selenium", "selenium", "Browser acceptance testing"),
    ("NumPy", "numpy", "Numerical arrays and calculations"),
    ("Pandas", "pandas", "Tabular data and CSV workflows"),
)

STRUCTURE_MODULES = (
    ("RDKit", "Structure parsing, CHEM_STANDARDIZER_V1, Crippen cLogP, TPSA, molecular properties, fingerprints and structural handling", "CALCULATED"),
    ("SyGMa", "Rule-based metabolism, metabolic soft spots and metabolite hypotheses", "LIMITED"),
    ("mordredcommunity", "Extended molecular descriptors", "READY"),
    ("descriptastorus", "Model-compatible molecular descriptors", "READY"),
    ("mhfp", "MinHash molecular fingerprints", "READY"),
    ("padelpy", "PaDEL descriptor adapter", "READY"),
    ("Chemprop", "Graph neural-network inference for installed prediction checkpoints", "READY"),
    ("scikit-learn", "Classical machine learning, similarity and validation utilities", "READY"),
    ("PyTorch", "ARM64 CPU neural-network inference runtime", "READY"),
    ("Lightning", "Neural-network runtime support", "READY"),
)

GLOSSARY = (
    ("READY", "The installed model or function is operational."),
    ("LOW CONFIDENCE", "The function is operational, but its scientific predictive confidence is limited."),
    ("MODEL_UNAVAILABLE", "No qualified local prediction model is installed for the endpoint."),
    ("EXPERIMENTAL", "A measured value recorded separately from predictions."),
    ("CALCULATED", "A deterministic calculation from a structure or recorded inputs."),
    ("PREDICTED", "An output produced by a prediction model."),
    ("RULE_ESTIMATE", "An estimate produced by explicit scientific rules."),
    ("DERIVED_ESTIMATE", "A value calculated from another estimate."),
    ("CONFORMAL_UNAVAILABLE", "The model can operate, but calibrated conformal uncertainty is unavailable."),
    ("OUT_OF_DOMAIN", "The compound lies outside the model's defined chemical space."),
)

LIMITATIONS = (
    "No full PBPK or PK/PD model is implemented.",
    "Animal efficacy prediction is not implemented.",
    "Several transporter endpoints remain MODEL_UNAVAILABLE.",
    "No qualified quantitative ML pKa model is installed; ionization support is rule-based.",
    "Conformal calibration is undercovered or unavailable for some operational endpoints.",
    "HLM, RLM and MLM predictions are operational but carry LOW scientific confidence.",
    "Rule-based metabolic soft spots and metabolites are hypotheses, not identified metabolites.",
    "Generated analogs are medicinal chemistry hypotheses and require synthesis and experimental validation.",
)


@lru_cache(maxsize=1)
def package_inventory() -> list[dict]:
    rows = []
    for label, distribution, purpose in PACKAGE_SPECS:
        if distribution is None:
            version, status = platform.python_version(), "READY"
        else:
            try:
                version, status = importlib.metadata.version(distribution), "READY"
            except importlib.metadata.PackageNotFoundError:
                version, status = "Not installed", "MODEL_UNAVAILABLE"
        rows.append({"package": label, "distribution": distribution, "version": version,
                     "purpose": purpose, "status": status})
    return rows


@lru_cache(maxsize=1)
def build_version() -> str:
    """Read the checkout revision without invoking git or requiring it in production."""
    git_dir = Path(__file__).resolve().parents[1] / ".git"
    try:
        head = (git_dir / "HEAD").read_text(encoding="utf-8").strip()
        if head.startswith("ref: "):
            ref = git_dir / head[5:]
            revision = ref.read_text(encoding="utf-8").strip()
        else:
            revision = head
        return revision[:12] if revision else "unavailable"
    except OSError:
        return "unavailable"


def structure_modules(inventory: list[dict]) -> list[dict]:
    versions = {row["package"]: row["version"] for row in inventory}
    return [{"module": name, "version": versions.get(name, "Not installed"),
             "used_for": purpose, "status": status if versions.get(name) != "Not installed" else "MODEL_UNAVAILABLE"}
            for name, purpose, status in STRUCTURE_MODULES]

