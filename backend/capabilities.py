"""Backend-owned capability registry and Dashboard capability aggregation.

Availability is derived from installed API routes and the live ADMET model
registry.  Confidence and conformal calibration remain independent metadata.
"""

from __future__ import annotations

from collections.abc import Iterable


FEATURE_REGISTRY = (
    # Structure and chemistry
    {"group": "structure", "key": "structure_drawing", "label": "Structure Drawing", "routes": ("/api/structure/validate",)},
    {"group": "structure", "key": "smiles_input", "label": "SMILES Input", "routes": ("/api/structure/validate",)},
    {"group": "structure", "key": "compound_versioning", "label": "Compound Versioning", "routes": ("/api/projects/{project_id}/compounds",)},
    {"group": "structure", "key": "structure_validation", "label": "Structure Validation", "routes": ("/api/structure/validate",)},
    {"group": "structure", "key": "physicochemical_properties", "label": "Physicochemical Properties", "routes": ("/api/compounds/{row_id}/calculate",)},
    {"group": "structure", "key": "structural_alerts", "label": "Structural Alerts", "routes": ("/api/compounds/{row_id}/calculate",)},
    # Activity and SAR
    {"group": "activity", "key": "experimental_activity", "label": "Experimental Activity", "routes": ("/api/assays/{assay_id}/measurements",)},
    {"group": "activity", "key": "activity_types", "label": "IC50 / EC50 / Ki / Kd / GI50", "routes": ("/api/assays/{assay_id}/measurements",)},
    {"group": "activity", "key": "project_qsar", "label": "Project QSAR", "routes": ("/api/assays/{assay_id}/models/train",), "availability": "LIMITED"},
    {"group": "activity", "key": "similarity", "label": "Similarity Analysis", "routes": ("/api/projects/{project_id}/compare",)},
    {"group": "activity", "key": "sar_mmp", "label": "SAR / MMP / Activity Cliff", "routes": ("/api/projects/{project_id}/sar", "/api/projects/{project_id}/mmp", "/api/projects/{project_id}/cliffs")},
    # Non-model ADME capabilities
    {"group": "adme", "key": "metabolic_soft_spots", "label": "Metabolic Soft Spots", "routes": ("/api/metabolism/predict/{version_id}",), "availability": "LIMITED"},
    {"group": "adme", "key": "metabolite_hypotheses", "label": "Metabolite Hypotheses", "routes": ("/api/metabolism/predict/{version_id}",), "availability": "LIMITED"},
    # Safety feature independent of predictive models
    {"group": "safety", "key": "structural_alerts", "label": "Structural Alerts", "routes": ("/api/compounds/{row_id}/calculate",)},
    # Optimization
    {"group": "optimization", "key": "liability_analysis", "label": "Liability Analysis", "routes": ("/api/projects/{project_id}/optimization/runs",)},
    {"group": "optimization", "key": "protected_regions", "label": "Protected / Modifiable Regions", "routes": ("/api/optimization/runs/{run_id}",)},
    {"group": "optimization", "key": "transformations", "label": "Medicinal Chemistry Transformations", "routes": ("/api/optimization/config",)},
    {"group": "optimization", "key": "analog_generation", "label": "Analog Generation", "routes": ("/api/proposals/{proposal_id}/execute",)},
    {"group": "optimization", "key": "rescoring", "label": "Re-scoring", "routes": ("/api/proposals/{proposal_id}/execute",), "availability": "LIMITED"},
    {"group": "optimization", "key": "pareto", "label": "Pareto Optimization", "routes": ("/api/proposals/{proposal_id}",)},
    {"group": "optimization", "key": "top_candidates", "label": "Top Candidate Selection", "routes": ("/api/proposals/{proposal_id}",)},
    # PK / DMPK feature registry
    {"group": "pk", "key": "experimental_pk", "label": "Experimental PK Data Management", "routes": ("/api/compounds/{row_id}/pk-studies",)},
    {"group": "pk", "key": "nca", "label": "NCA", "routes": ("/api/pk-studies/{study_id}/run-nca",)},
    {"group": "pk", "key": "ivive", "label": "IVIVE / Hepatic Clearance", "routes": ("/api/compound-versions/{version_id}/ivive/run", "/api/ivive/methods")},
    {"group": "pk", "key": "distribution_absorption", "label": "Vd / Absorption Foundation", "routes": ("/api/compound-versions/{version_id}/pk-foundation",)},
    {"group": "pk", "key": "iv_simulation", "label": "IV Simulation", "routes": ("/api/compound-versions/{version_id}/pk-simulation/run",)},
    {"group": "pk", "key": "extravascular_simulation", "label": "PO / SC / IP Simulation", "routes": ("/api/compound-versions/{version_id}/pk-simulation/run",)},
    {"group": "pk", "key": "cross_species", "label": "Cross-Species Scaling", "routes": ("/api/compound-versions/{version_id}/translational-pk",)},
    {"group": "pk", "key": "human_pk", "label": "Human Translational PK", "routes": ("/api/compound-versions/{version_id}/human-pk/profile", "/api/compound-versions/{version_id}/human-pk/simulation/run")},
    {"group": "pk", "key": "prospective_freeze", "label": "Prospective Prediction Freeze", "routes": ("/api/compound-versions/{version_id}/human-pk/freeze-snapshot",)},
    {"group": "pk", "key": "retrospective_validation", "label": "Retrospective Validation", "routes": ("/api/compound-versions/{version_id}/human-pk/validation",)},
)


GROUPS = (
    ("structure", "Structure & Chemistry", "Version-controlled chemical identity and calculated molecular properties."),
    ("activity", "Activity & SAR", "Project-local experimental activity, modeling, and structure–activity evidence."),
    ("adme", "ADME", "Absorption, distribution, and metabolic-liability evidence."),
    ("cyp_transporters", "CYP & Transporters", "Endpoint- and role-separated metabolism and transporter classifications."),
    ("safety", "Safety / Toxicology", "Classification evidence and calculated structural safety alerts."),
    ("optimization", "Optimization", "Deterministic strategy, analog generation, filtering, and transparent ranking."),
    ("pk", "PK / DMPK", "Experimental, preclinical, translational, and human PK workflows."),
)


ADME_MODELS = (
    ("Solubility", "Solubility"),
    ("Permeability", "Caco-2"),
    ("Plasma protein binding", "PPB"),
    ("Plasma protein binding", "fu", "fu_from_ppb"),
    ("HLM intrinsic clearance", "HLM"),
    ("RLM intrinsic clearance", "RLM"),
    ("MLM intrinsic clearance", "MLM"),
)

CYP_TRANSPORTER_ORDER = (
    "CYP1A2 inhibitor", "CYP2C9 inhibitor", "CYP2C19 inhibitor", "CYP2D6 inhibitor", "CYP3A4 inhibitor",
    "CYP2C9 substrate", "CYP2D6 substrate", "CYP3A4 substrate", "CYP1A2 substrate", "CYP2C19 substrate",
    "P-gp inhibitor", "P-gp substrate", "BCRP inhibitor", "BCRP substrate", "BSEP inhibitor",
    "OATP1B1 inhibitor", "OATP1B3 inhibitor", "OCT1 inhibitor", "OCT2 inhibitor", "MATE1 inhibitor", "MATE2-K inhibitor",
)

SAFETY_ORDER = (
    "hERG liability", "Ames mutagenicity", "DILI clinical liability", "Mitochondrial toxicity",
    "General cytotoxicity", "Skin sensitization", "BBB penetration", "CNS liability",
)

MODEL_LABELS = {
    "hERG liability": "hERG",
    "Ames mutagenicity": "Ames",
    "DILI clinical liability": "DILI",
}


def _group_status(items: list[dict]) -> str:
    states = {item["availability"] for item in items}
    supported = states & {"READY", "LIMITED"}
    if not supported:
        return "MODEL_UNAVAILABLE"
    if "MODEL_UNAVAILABLE" in states:
        return "PARTIAL"
    if "LIMITED" in states:
        return "LIMITED"
    return "READY"


def _feature_item(entry: dict, route_paths: set[str]) -> dict:
    registered = all(path in route_paths for path in entry["routes"])
    availability = entry.get("availability", "READY") if registered else "MODEL_UNAVAILABLE"
    return {
        "key": entry["key"],
        "label": entry["label"],
        "availability": availability,
        "confidence": "NOT_APPLICABLE",
        "conformal_status": "NOT_APPLICABLE",
        "source": "BACKEND_FEATURE_REGISTRY",
        "required_routes": list(entry["routes"]),
    }


def _model_item(model: dict, label: str | None = None, key: str | None = None) -> dict:
    return {
        "key": key or model["endpoint"],
        "label": label or MODEL_LABELS.get(model["endpoint"], model["endpoint"]),
        "endpoint": model["endpoint"],
        "availability": model.get("availability", model.get("status", "MODEL_UNAVAILABLE")),
        "confidence": model.get("confidence", "NOT_APPLICABLE"),
        "conformal_status": model.get("conformal_status", "CONFORMAL_UNAVAILABLE"),
        "source": "ADMET_MODEL_REGISTRY",
        "model_id": model.get("id"),
        "model_name": model.get("model_name"),
        "model_version": model.get("model_version"),
        "unavailable_reason": model.get("unavailable_reason", ""),
    }


def build_capability_summary(model_rows: Iterable[dict], route_paths: Iterable[str], stage: str = "5B-4") -> dict:
    """Compose the Dashboard from live model rows and registered backend routes."""
    models = {row["endpoint"]: row for row in model_rows}
    paths = set(route_paths)
    feature_items: dict[str, list[dict]] = {}
    for entry in FEATURE_REGISTRY:
        feature_items.setdefault(entry["group"], []).append(_feature_item(entry, paths))

    adme_items = []
    for spec in ADME_MODELS:
        endpoint, label, *key_parts = spec
        if endpoint in models:
            adme_items.append(_model_item(models[endpoint], label, key_parts[0] if key_parts else None))
    adme_items.extend(feature_items.get("adme", []))
    cyp_items = [_model_item(models[name]) for name in CYP_TRANSPORTER_ORDER if name in models]
    safety_items = [_model_item(models[name]) for name in SAFETY_ORDER if name in models] + feature_items.get("safety", [])

    by_group = {
        "structure": feature_items.get("structure", []),
        "activity": feature_items.get("activity", []),
        "adme": adme_items,
        "cyp_transporters": cyp_items,
        "safety": safety_items,
        "optimization": feature_items.get("optimization", []),
        "pk": feature_items.get("pk", []),
    }
    groups = [
        {"key": key, "title": title, "description": description, "status": _group_status(by_group[key]), "items": by_group[key]}
        for key, title, description in GROUPS
    ]
    return {
        "stage": stage,
        "source": "BACKEND_CAPABILITY_REGISTRY",
        "groups": groups,
    }
