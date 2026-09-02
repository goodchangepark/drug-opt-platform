"""Runtime platform metadata used by the researcher-facing Help page."""

from __future__ import annotations

from functools import lru_cache
import importlib.metadata
import platform
from pathlib import Path


from typing import Any

APP_VERSION = "1.0.0"
CURRENT_STAGE = "5B-4"
# Internal engineering stage is retained for provenance; this is the
# researcher-facing operational milestone shown by the dashboard and Help.
CURRENT_STAGE_LABEL = "Internal Validation"
CURRENT_STAGE_STATUS = "Prediction Engine v1 Frozen"
CURRENT_STAGE_SUBSTATUS = "Experimental data collection active"

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
    ("MODEL_UNAVAILABLE", "No qualified local prediction model is installed for the endpoint; outputs are not fabricated."),
    ("EXPERIMENTAL", "A measured value recorded separately from predictions; takes precedence over predictions."),
    ("CALCULATED", "A deterministic calculation from a structure or recorded inputs (e.g. Crippen cLogP, TPSA)."),
    ("PREDICTED", "An output produced by a predictive model or endpoint-specific strategy."),
    ("RULE_ESTIMATE", "An estimate produced by explicit scientific rules (e.g. rule-based pKa)."),
    ("DERIVED_ESTIMATE", "A value calculated from another estimate (e.g. derived logD 7.4)."),
    ("CONFORMAL_UNAVAILABLE", "The model can operate, but calibrated conformal uncertainty is unavailable."),
    ("OUT_OF_DOMAIN", "The compound lies outside the model's defined applicability domain."),
    ("SINGLE_CORE_MODEL", "Production prediction driven by a single validated core model."),
    ("RANK_FUSION", "Score integration by rank aggregation rather than direct raw score averaging (e.g. metabolic soft spots)."),
    ("SHADOW", "Research candidate model executing asynchronously or in evaluation mode without altering production predictions."),
    ("MANUAL_PROMOTION_REQUIRED", "Governance policy requiring explicit human scientific signoff before promoting shadow models to active production."),
    ("READY_FOR_INTERNAL_VALIDATION", "Current system readiness: verified for controlled internal research and prospective qualification (does not imply clinical/regulatory claims)."),
)

LIMITATIONS = (
    "Drug-OPT is a developability & DMPK optimization platform; it does NOT provide target docking, primary efficacy prediction, or autonomous AI-scientist optimization.",
    "Quantitative ML pKa and ML logD models are not installed; pKa is a RULE_ESTIMATE and logD7.4 is a DERIVED_ESTIMATE (cLogP must not be equated to logD).",
    "Transporter endpoints (P-gp substrate, BCRP, BSEP, OATP1B1/3, OCT1/2, MATE1/2-K) remain MODEL_UNAVAILABLE.",
    "Microsomal clearance models are species-specific (HLM, RLM, MLM); Dog/Monkey/unspecified microsomal clearances are MODEL_UNAVAILABLE.",
    "hERG liability is single-core raw M1 in production; Platt calibration is research-stage only; secondary models are shadow-only.",
    "CYP3A4 inhibitor is single-core M1; fixed blend and adaptive weighting are research/shadow-only.",
    "Solubility adaptive weighting is research shadow-only; Caco-2 static consensus is shadow-only (INSUFFICIENT_EVIDENCE).",
    "Metabolic soft spots use SyGMa+SMARTCyp rank fusion; metabolite structures are rule-based hypotheses, not identified metabolites.",
    "Generated analogs are medicinal chemistry hypotheses requiring chemical synthesis and experimental validation.",
    "Mechanistic PK/NCA/IVIVE and human translation are mathematical models and do not constitute clinical or regulatory validation.",
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
    {
        "version": "0.6.4-stage4d3b2a",
        "release_date": "2026-08-29",
        "stage": "Stage 4D-3B2A",
        "milestone": "hERG Calibration & Model Quality Audit",
        "improvements": "Independent external evaluation, threshold and Platt calibration research candidate audit; identified secondary model limitations while keeping production hERG single-core strategy strictly unchanged.",
    },
    {
        "version": "0.6.5-stage4d4",
        "release_date": "2026-08-29",
        "stage": "Stage 4D-4",
        "milestone": "Endpoint-Specific Prediction Strategy Governance",
        "improvements": "Finalized 49 endpoint-specific prediction policies across 22 endpoint contracts and 40 runtime endpoints; codified single-core, rank-fusion, rule/derived estimates, and mechanistic PK isolation without modifying production predictions.",
    },
    {
        "version": "0.6.6-stage4d5",
        "release_date": "2026-08-29",
        "stage": "Stage 4D-5",
        "milestone": "Production Qualification & Prospective Validation Governance",
        "improvements": "Prospective qualification lifecycle (SHADOW -> VALIDATED -> PRODUCTION_CANDIDATE -> ACTIVE -> ROLLBACK), promotion gates, drift review policy, immutable qualification evidence, strategy cards, and verified scientific readiness (READY_FOR_INTERNAL_VALIDATION).",
    },
    {
        "version": "0.6.7-stage4d6",
        "release_date": "2026-08-29",
        "stage": "Stage 4D-6",
        "milestone": "Prediction Runtime Integration — Real Multi-Model Execution",
        "improvements": (
            "Connected Stage 4D multimodel scientific governance to the real Save & Predict runtime. "
            "New canonical PredictionOrchestrator executes all endpoint-authorized CORE + SHADOW/secondary models "
            "per EndpointStrategyPolicy (ESOL M2 for Solubility, physchem Caco-2 M2, Morgan CYP3A4 M2, physchem hERG M2). "
            "Authorized secondary outputs execute each prediction run and are recorded in immutable "
            "qualification prediction-freeze provenance without becoming legacy primary registry rows. "
            "Research-only constraints enforced: CYP3A4 fixed blend (0.9578/0.0422) shadow-only, "
            "hERG raw M1 remains production (threshold 0.50), RDKIT-GBR M3 ADAPTIVE_EXCLUDED. "
            "Multimodel-provenance endpoint enhanced with shadow_model_count, core_model_count. "
            "Legacy single-model path preserved as fallback."
        ),
    },
    {
        "version": "0.6.8-stage4d7",
        "release_date": "2026-08-29",
        "stage": "Stage 4D-7",
        "milestone": "Pre-Experimental Prediction Optimization",
        "improvements": (
            "Independently audited endpoint-specific initial-prediction candidates using existing "
            "leakage-safe Stage 4D evidence. Production changes occur only after Stage 4D-5 "
            "qualification gates pass; this review retained current policies where evidence or "
            "configured promotion requirements were insufficient."
        ),
    },
    {
        "version": "0.6.9-stage4e1",
        "release_date": "2026-08-29",
        "stage": "Stage 4E-1",
        "milestone": "Model Landscape & Qualification Planning",
        "improvements": "Prioritized model gaps; audited candidate model and dataset sources, licenses, overlap risk, and ARM64 feasibility; defined a Stage 4E-2 qualification queue without installing models or changing production predictions.",
    },
    {
        "version": "0.7.0-stage4e2",
        "release_date": "2026-08-29",
        "stage": "Stage 4E-2",
        "milestone": "Candidate Acquisition & Runtime Feasibility Qualification",
        "improvements": "Audited candidate licenses, checkpoints, endpoint compatibility, and Xavier ARM64 feasibility; finalized a fail-closed Stage 4E-3 benchmark gate with no production prediction change.",
    },
    {
        "version": "0.7.1-stage4e2r",
        "release_date": "2026-08-29",
        "stage": "Stage 4E-2R",
        "milestone": "Model Acquisition Blocker Resolution",
        "improvements": "Re-evaluated blocked candidates, acquired a pinned licensed Caco-2 external benchmark, and documented unresolved model gates; production prediction engine unchanged.",
    },
    {
        "version": "0.7.2-stage4e3a",
        "release_date": "2026-08-29",
        "stage": "Stage 4E-3A",
        "milestone": "Caco-2 Independent Scientific Benchmark",
        "improvements": "Benchmarked the existing Caco-2 CORE and authorized SHADOW on the pinned ExpansionRx Papp A→B cohort without refitting; completed overlap, censor, duplicate, AD, scaffold, disagreement, and complementarity analyses. No numeric consensus existed or was created, and production policy remains unchanged.",
    },
    {
        "version": "0.7.3-stage4e3b",
        "release_date": "2026-08-29",
        "stage": "Stage 4E-3B",
        "milestone": "Caco-2 Final Model Expansion & Closure",
        "improvements": "Completed the final targeted Caco-2 candidate search under license, checkpoint, endpoint, lineage, and ARM64 gates; no qualified replacement was found, so the current CORE is frozen for Engine-v1. No candidate was installed, registered, benchmarked, or promoted.",
    },
    {
        "version": "0.7.4-stage4e3c",
        "release_date": "2026-08-29",
        "stage": "Stage 4E-3C",
        "milestone": "hERG Final Qualification & Closure",
        "improvements": "Reviewed raw hERG M1 and the historical leakage-safe Platt holdout, separated calibration from discrimination, and completed the final secondary-model acquisition attempt. No qualified replacement or production calibration was found; raw M1 at threshold 0.50 remains frozen for Engine-v1.",
    },
    {
        "version": "0.7.5-stage4e3d",
        "release_date": "2026-08-29",
        "stage": "Stage 4E-3D",
        "milestone": "Clearance Final Qualification & Closure",
        "improvements": "HLM/RLM/MLM contracts and external evidence audited; species-specific Engine-v1 limitations frozen; production unchanged.",
    },
    {
        "version": "0.7.6-stage4e3e",
        "release_date": "2026-08-29",
        "stage": "Stage 4E-3E",
        "milestone": "pKa / logD Final Qualification & Closure",
        "improvements": "Quantitative pKa/logD candidates reviewed with strict checkpoint and endpoint gates; rule/derived Engine-v1 limitations frozen; production unchanged.",
    },
    {
        "version": "0.8.0-engine-v1",
        "release_date": "2026-08-29",
        "stage": "Stage 4E-4",
        "milestone": "Prediction Engine v1 Freeze",
        "improvements": "Finalized 49 endpoint-specific policies, model roles, reliability dimensions, unavailable states, deterministic policy hash, and Engine-v1 prospective-freeze provenance. Ready for controlled internal prospective validation; not clinical or regulatory validation.",
    },
    {
        "version": "0.8.1-validation-framework",
        "release_date": "2026-08-29",
        "stage": "Stage 6 — Internal Validation",
        "milestone": "Engine v1 Internal Prospective Validation Framework",
        "improvements": (
            "Created internal prospective validation framework for Engine v1: "
            "campaign entity (IVC-engine-v1-2026-08-29), 3 GLP-1 compound cohort enrollment, "
            "18 prediction freezes registered, blinding enforcement (prediction-before-experiment), "
            "evidence classification (TRUE_PROSPECTIVE / BLINDED_RETROSPECTIVE / HISTORICAL_VISIBLE), "
            "endpoint compatibility contracts, experiment import pipeline, "
            "analysis engine (regression/classification/bootstrap/AD/reliability/scaffold), "
            "14 validation artifacts, 48 new tests (all passing). "
            "Framework status: READY. Scientific validation: COLLECTING — awaiting internal experimental data."
        ),
    },
    {
        "version": "1.0.0",
        "release_date": "2026-08-30",
        "stage": "Product Release",
        "milestone": "Drug-OPT v1.0 — Prediction Engine v1 Baseline",
        "improvements": (
            "Prediction Engine v1 frozen with 49 endpoint policies; endpoint-specific model, rule, derived, "
            "mechanistic and MODEL_UNAVAILABLE states finalized; immutable prospective prediction freeze and "
            "same-compound leakage protection enabled; internal prospective validation framework ready for data "
            "collection; test/demo projects removed while scientific provenance was preserved."
        ),
    },
    {
        "version": "v3.5",
        "release_date": "2026-09-01",
        "stage": "Development milestone",
        "milestone": "Unified Experimental × Prediction UX",
        "improvements": "Unified experimental, prediction, difference and project-learning presentation with immutable comparison context.",
    },
    {
        "version": "v3.6",
        "release_date": "2026-09-01",
        "stage": "Development milestone",
        "milestone": "Project Learning Curve Validation",
        "improvements": "Leakage-safe repeated holdout learning curves and conservative endpoint-specific adaptation evidence.",
    },
    {
        "version": "v3.7",
        "release_date": "2026-09-01",
        "stage": "Development milestone",
        "milestone": "Continuous Project Learning",
        "improvements": "Persisted experiment feedback, candidate adapter validation, explicit activation, rollback and future-only learning.",
    },
    {
        "version": "v3.8A",
        "release_date": "2026-09-01",
        "stage": "Development milestone",
        "milestone": "Persisted Experimental × Prediction Foundation",
        "improvements": "External search candidates, immutable prediction snapshots and a DB-backed endpoint comparison contract that survives reload and restart.",
    },
    {
        "version": "v3.8B",
        "release_date": "2026-09-01",
        "stage": "Development milestone",
        "milestone": "Canonical Endpoint Harmonization",
        "improvements": "Versioned semantic endpoint and unit registries for ADMET, metabolism and species/route-aware PK comparisons.",
    },
    {
        "version": "v3.9",
        "release_date": "2026-09-01",
        "stage": "Development milestone",
        "milestone": "Prediction Coverage Expansion",
        "improvements": "Prediction inventory, Stage-5 output indexing and explicit MODEL, MECHANISTIC, RULE and DERIVED provenance states.",
    },
    {
        "version": "v4.0",
        "release_date": "2026-09-01",
        "stage": "Development milestone",
        "milestone": "Qualification Contract",
        "improvements": "Layered identity, reference, numeric, endpoint, context, pairability, comparability, import and adaptation qualification stages.",
    },
    {
        "version": "v4.1",
        "release_date": "2026-09-02",
        "stage": "Development milestone",
        "milestone": "Persistence + Unified Scientific Comparison",
        "improvements": "Navigation-safe restoration and canonical Activity, ADMET, Metabolism and species-first PK comparison tables.",
    },
    {
        "version": "v4.2",
        "release_date": "2026-09-02",
        "stage": "Scientific validation milestone",
        "milestone": "Scientific Prediction Validation & PK Accuracy Profiling",
        "improvements": "PK mechanistic provenance and input-completeness audit, endpoint-specific numeric error reporting, prediction history audit and transparent performance limitations.",
    },
    {
        "version": "v4.3",
        "release_date": "2026-09-02",
        "stage": "Scientific results milestone",
        "milestone": "Final Scientific Results & PK Comparison",
        "improvements": "Cleaner project-learning overview, normalized experimental-to-prediction scientific rows, measurement-type-aware metabolism presentation, species/context-specific PK comparison, and regulatory clinical-PK context qualification.",
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
    # The frontend historically used ``date``/``highlights`` while the
    # backend contract uses the more explicit names below.  Return both so
    # old clients and the current Help page render the same audited history.
    return [
        {**row, "date": row["release_date"], "highlights": row["improvements"]}
        for row in VERSION_HISTORY
    ]


def latest_release_date() -> str:
    """Return the release date of the latest product version."""
    if VERSION_HISTORY:
        return VERSION_HISTORY[-1]["release_date"]
    return "2026-08-29"
