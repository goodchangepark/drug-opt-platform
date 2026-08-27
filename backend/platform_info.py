"""Runtime platform metadata used by the researcher-facing Help page."""

from __future__ import annotations

from functools import lru_cache
import importlib.metadata
import platform
from pathlib import Path


from typing import Any

APP_VERSION = "0.6.3-stage5b4-ui"
CURRENT_STAGE = "5B-4"

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

VERSION_HISTORY = [
    {
        "version": "0.1.0",
        "release_date": "2026-08-20",
        "stage": "Stage 1",
        "milestone": "Core Cheminformatics & Compound Workspace",
        "improvements": "Compound registration, CHEM_STANDARDIZER_V1, 2D structure rendering, physicochemical descriptor calculations (MW, cLogP, TPSA, HBD, HBA, RotB, QED).",
    },
    {
        "version": "0.2.0",
        "release_date": "2026-08-22",
        "stage": "Stage 2",
        "milestone": "Interactive Structure Editor & Similarity Search",
        "improvements": "Interactive molecular structure drawing, substructure filtering, Morgan and RDKit fingerprint similarity search, compound versioning.",
    },
    {
        "version": "0.3.0",
        "release_date": "2026-08-24",
        "stage": "Stage 3",
        "milestone": "ADMET Prediction Suite & Model Registry",
        "improvements": "OpenADMET and Chemprop model inference for solubility, permeability, PPB, microsomal stability, CYP450 panel, and cardiac/liver safety (hERG, DILI, Ames).",
    },
    {
        "version": "0.4.0",
        "release_date": "2026-08-25",
        "stage": "Stage 4A",
        "milestone": "Optimization Strategy Engine & MMP Evidence",
        "improvements": "Lead optimization strategy generation, liability identification, protected vs modifiable atom mapping, matched molecular pair (MMP) transforms.",
    },
    {
        "version": "0.4.1",
        "release_date": "2026-08-25",
        "stage": "Stage 4B",
        "milestone": "Analog Generation Engine & Pareto Ranking",
        "improvements": "Rule-based and fragment-based analog generation, multi-parameter optimization (MPO), Pareto frontier scoring, candidate ranking workspace.",
    },
    {
        "version": "0.4.2",
        "release_date": "2026-08-26",
        "stage": "Stage 4C",
        "milestone": "Scientific Hardening & Conformal Uncertainty",
        "improvements": "Conformal prediction intervals, applicability domain verification, pKa ionization engine and pH-dependent Henderson-Hasselbalch profiles.",
    },
    {
        "version": "0.5.0",
        "release_date": "2026-08-26",
        "stage": "Stage 5A",
        "milestone": "Experimental PK, NCA & Mechanistic IVIVE",
        "improvements": "Noncompartmental analysis (NCA) engine for concentration-time curves, physiological scaling, well-stirred hepatic clearance IVIVE.",
    },
    {
        "version": "0.5.1",
        "release_date": "2026-08-27",
        "stage": "Stage 5B-1/2",
        "milestone": "Multi-Route PK Simulation & Absorption Kinetics",
        "improvements": "1-compartment and 2-compartment pharmacokinetic differential simulation for IV, PO, SC, and IP dosing; absorption rate (ka) parameter estimation.",
    },
    {
        "version": "0.5.2",
        "release_date": "2026-08-27",
        "stage": "Stage 5B-3",
        "milestone": "Cross-Species Translation & LOSO Validation",
        "improvements": "Allometric scaling across mouse, rat, dog, monkey; Leave-One-Species-Out (LOSO) cross-validation and fold-error evaluation.",
    },
    {
        "version": "0.6.0-stage5b4-stable",
        "release_date": "2026-08-27",
        "stage": "Stage 5B-4",
        "milestone": "Human Translational PK & Platform Stabilization",
        "improvements": "Human clearance assembly hierarchy (Experimental > Allometry > IVIVE), clinical simulation, prospective snapshot freezing, database audit.",
    },
    {
        "version": "0.6.1-stage5b4-ui",
        "release_date": "2026-08-27",
        "stage": "Stage 5B-4 Refinement",
        "milestone": "Scientific UI/UX Refinement & Multi-Species Visualization",
        "improvements": "Main dashboard simplification, Noto Sans KR typography standard, centralized interpretation registry, multi-species comparative PK summary and interactive curves.",
    },
    {
        "version": "0.6.2-stage5b4-ui",
        "release_date": "2026-08-27",
        "stage": "Stage 5B-4 Refinement 2",
        "milestone": "UI Polish & Unified Prediction Workflow",
        "improvements": "Dashboard card grid redesign, compound overview primary predict workflow, tab-level re-predict actions, PK overview summary, global Noto Sans KR font, and test project cleanup.",
    },
    {
        "version": "0.6.3-stage5b4-ui",
        "release_date": "2026-08-27",
        "stage": "Stage 5B-4 Refinement 3",
        "milestone": "Dashboard Redesign & Compound Save Workflow Restoration",
        "improvements": "Clean 3-column scientific workspace dashboard, streamlined single-card platform overview, restored robust compound save and initial versioning workflow with instant property calculation.",
    },
]



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
        head_file = git_dir / "HEAD"
        if not head_file.is_file():
            return "6e56deed6e85"
        head_content = head_file.read_text(encoding="utf-8").strip()
        if head_content.startswith("ref: "):
            ref_path = git_dir / head_content[5:]
            if ref_path.is_file():
                return ref_path.read_text(encoding="utf-8").strip()[:12]
        return head_content[:12]
    except Exception:
        return "6e56deed6e85"


def structure_modules(inventory: list[dict]) -> list[dict]:
    versions = {row["package"]: row["version"] for row in inventory}
    return [
        {
            "module": name,
            "version": versions.get(name, "Not installed"),
            "used_for": purpose,
            "status": status if versions.get(name) != "Not installed" else "MODEL_UNAVAILABLE",
        }
        for name, purpose, status in STRUCTURE_MODULES
    ]


def version_history() -> list[dict[str, Any]]:
    """Return the curated product evolution history."""
    return list(VERSION_HISTORY)


def latest_release_date() -> str:
    """Return the release date of the latest product version."""
    if VERSION_HISTORY:
        return VERSION_HISTORY[-1]["release_date"]
    return "2026-08-27"
