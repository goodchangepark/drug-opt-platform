"""Prediction-first, read-only compound developability profile.

The profile is a presentation contract over Stable Core.  It deliberately
contains expected rows that have no value so clients never need to infer the
difference between not-yet-run, unavailable, mechanistic, and contextual
endpoints.  Scientific values are sourced only through
``build_scientific_endpoint_rows``; legacy calculation tables are not queried.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from sqlalchemy import select

from .canonical_endpoints import REGISTRY
from .prediction_engine_registry import (
    ROUTE_MODEL_UNAVAILABLE,
    ROUTE_RETAIN_V3_3,
    get_current_production_engine_info,
    get_current_production_routing,
)
from .model_artifact_authority import model_artifact_registration
from .scientific_core_service import build_scientific_endpoint_rows


AVAILABLE_CURRENT = "AVAILABLE_CURRENT"
ON_DEMAND = "ON_DEMAND"
MECHANISTIC_ONLY = "MECHANISTIC_ONLY"
MODEL_UNAVAILABLE = "MODEL_UNAVAILABLE"
MODEL_NOT_REGISTERED = "MODEL_NOT_REGISTERED"
CONTEXT_REQUIRED = "CONTEXT_REQUIRED"
CURRENT_DATA_CEILING = "CURRENT_DATA_CEILING"
CALCULATED_BUT_NOT_ELIGIBLE = "CALCULATED_BUT_NOT_ELIGIBLE"
FAILED = "FAILED"


@dataclass(frozen=True)
class ProfileEndpoint:
    endpoint: str
    group: str
    display_name: str
    species: str = "HUMAN"
    importance: str = "CORE"
    availability: str | None = None
    semantic: str = ""
    aliases: tuple[str, ...] = ()
    source_endpoint: str = ""
    unit: str = ""


def _p(endpoint: str, group: str, display: str, **kwargs: Any) -> ProfileEndpoint:
    return ProfileEndpoint(endpoint, group, display, **kwargs)


# Ordered deliberately: this is also the desktop/mobile display order and the
# stable query catalog for future integration clients.
CORE_PROFILE: tuple[ProfileEndpoint, ...] = (
    _p("MW", "physchem", "MW"),
    _p("CLOGP", "physchem", "cLogP"),
    _p("LOGD_7_4", "physchem", "logD7.4", availability=MECHANISTIC_ONLY),
    _p("PKA", "physchem", "pKa", availability=MECHANISTIC_ONLY),
    _p("TPSA", "physchem", "TPSA"),
    _p("HBD", "physchem", "HBD"),
    _p("HBA", "physchem", "HBA"),
    _p("ROTB", "physchem", "Rotatable Bonds"),
    _p("FSP3", "physchem", "Fsp3"),
    _p("QED", "physchem", "QED", importance="SECONDARY"),
    _p("FORMAL_CHARGE", "physchem", "Formal Charge", importance="SECONDARY"),
    _p("HEAVY_ATOM_COUNT", "physchem", "Heavy Atom Count", importance="SECONDARY"),
    _p("SOLUBILITY_GENERIC", "physchem", "Solubility", aliases=("HUMAN_SOLUBILITY",)),
    _p("CACO2_PAPP_AB", "absorption", "Caco-2 permeability", aliases=("CACO2_PERMEABILITY",)),
    _p("PAMPA_PERMEABILITY", "absorption", "PAMPA permeability", availability=MODEL_UNAVAILABLE, unit="cm/s"),
    _p("HIA", "absorption", "HIA", availability=MODEL_NOT_REGISTERED),
    _p("HUMAN_PPB", "distribution", "Human PPB", aliases=("PPB",)),
    _p("HUMAN_FU", "distribution", "Human fu", aliases=("FU",), source_endpoint="HUMAN_PPB", semantic="FRACTION_UNBOUND", availability=CURRENT_DATA_CEILING),
    _p("BBB_PENETRATION", "distribution", "BBB", availability=MODEL_NOT_REGISTERED),
    _p("VDSS", "distribution", "Vdss", availability=CURRENT_DATA_CEILING),
    _p("HLM_CLINT", "metabolic_stability", "HLM (microsomal stability)"),
    _p("RLM_CLINT", "metabolic_stability", "RLM (microsomal stability)", species="RAT"),
    _p("MLM_CLINT", "metabolic_stability", "MLM (microsomal stability)", species="MOUSE"),
    _p("PLASMA_STABILITY", "metabolic_stability", "Plasma Stability (PS)", availability=MODEL_UNAVAILABLE, unit="min"),
    *tuple(
        _p(f"CYP{isoform}_INHIBITION", "cyp", f"CYP{isoform} quantitative inhibition", semantic="QUANTITATIVE_INHIBITION")
        for isoform in ("1A2", "2C9", "2C19", "2D6", "3A4")
    ),
    *tuple(
        _p(f"CYP{isoform}_INHIBITOR_CLASS", "cyp", f"CYP{isoform} inhibitor classification", semantic="INHIBITOR_CLASSIFICATION")
        for isoform in ("1A2", "2C9", "2C19", "2D6", "3A4")
    ),
    *tuple(
        _p(f"CYP{isoform}_SUBSTRATE", "cyp", f"CYP{isoform} substrate classification", semantic="SUBSTRATE_CLASSIFICATION")
        for isoform in ("2C9", "2D6", "3A4")
    ),
    _p("PGP_INHIBITION", "transporters", "P-gp inhibitor classification", semantic="INHIBITOR_CLASSIFICATION", aliases=("PGP",)),
    _p("PGP_INHIBITION_QUANT", "transporters", "P-gp quantitative inhibition", semantic="QUANTITATIVE_INHIBITION"),
    _p("BCRP_INHIBITOR", "transporters", "BCRP inhibitor classification", semantic="INHIBITOR_CLASSIFICATION", aliases=("BCRP",)),
    _p("BCRP_INHIBITOR_QUANT", "transporters", "BCRP quantitative inhibition", semantic="QUANTITATIVE_INHIBITION"),
    _p("OATP1B1_INHIBITOR", "transporters", "OATP1B1"),
    _p("OATP1B3_INHIBITOR", "transporters", "OATP1B3"),
    _p("OCT1_INHIBITOR", "transporters", "OCT1"),
    _p("OCT2_INHIBITOR", "transporters", "OCT2"),
    _p("HERG_LIABILITY", "safety", "hERG quantitative inhibition", semantic="QUANTITATIVE_INHIBITION", aliases=("HERG_IC50",)),
    _p("HERG_CLASS", "safety", "hERG risk classification", semantic="RISK_CLASSIFICATION"),
    _p("AMES_MUTAGENICITY", "safety", "Ames", semantic="RISK_CLASSIFICATION", aliases=("AMES",), availability=MODEL_UNAVAILABLE),
    _p("DILI_LIABILITY", "safety", "DILI", semantic="RISK_CLASSIFICATION", aliases=("DILI",)),
    _p("METABOLIC_SOFT_SPOTS", "metabolism", "Metabolic soft spots", availability=MECHANISTIC_ONLY),
    _p("METABOLITE_HYPOTHESES", "metabolism", "Predicted metabolites", availability=MECHANISTIC_ONLY),
    _p("HLM_CLINT", "pk", "HLM / Clint (upstream)", semantic="UPSTREAM", source_endpoint="HLM_CLINT"),
    _p("HUMAN_FU", "pk", "Human fu (upstream)", semantic="UPSTREAM_FU", source_endpoint="HUMAN_PPB", availability=CURRENT_DATA_CEILING),
    _p("HUMAN_HEPATIC_CL", "pk", "Human hepatic CL", availability=CURRENT_DATA_CEILING, semantic="PK_PARAMETER", aliases=("HEPATIC_CL",), unit="mL/min/kg"),
    _p("HUMAN_PK_CL_IV", "pk", "Human total/systemic IV CL", availability=CURRENT_DATA_CEILING, semantic="PK_PARAMETER"),
    _p("HUMAN_PK_CL_UNSPECIFIED", "pk", "Human systemic CL (context incomplete)", availability=CONTEXT_REQUIRED, semantic="PK_PARAMETER"),
    _p("HUMAN_PK_VD_IV", "pk", "Human IV volume of distribution", availability=CURRENT_DATA_CEILING, semantic="PK_PARAMETER"),
    _p("VDSS", "pk", "Human Vdss", availability=CURRENT_DATA_CEILING, semantic="PK_PARAMETER", source_endpoint="VDSS"),
    _p("HUMAN_PK_T_HALF_IV", "pk", "Human IV half-life", availability=CURRENT_DATA_CEILING, semantic="PK_PARAMETER"),
    _p("HUMAN_PK_F_ORAL", "pk", "Oral bioavailability (F)", availability=CONTEXT_REQUIRED, semantic="CONTEXTUAL_PK"),
    _p("HUMAN_PK_KA_ORAL", "pk", "Oral absorption rate (ka)", availability=CONTEXT_REQUIRED, semantic="CONTEXTUAL_PK", aliases=("KA",)),
    _p("HUMAN_PK_AUC_ORAL", "pk", "Oral AUC", availability=CONTEXT_REQUIRED, semantic="CONTEXTUAL_PK"),
    _p("HUMAN_PK_CMAX_ORAL", "pk", "Oral Cmax", availability=CONTEXT_REQUIRED, semantic="CONTEXTUAL_PK"),
    _p("HUMAN_PK_CMAX_UNSPECIFIED", "pk", "Human Cmax (context incomplete)", availability=CONTEXT_REQUIRED, semantic="CONTEXTUAL_PK"),
    _p("HUMAN_PK_TMAX_ORAL", "pk", "Oral Tmax", availability=CONTEXT_REQUIRED, semantic="CONTEXTUAL_PK"),
    _p("HUMAN_PK_CLF_ORAL", "pk", "Oral CL/F", availability=CONTEXT_REQUIRED, semantic="CONTEXTUAL_PK"),
    _p("HUMAN_PK_VDF_ORAL", "pk", "Oral Vd/F", availability=CONTEXT_REQUIRED, semantic="CONTEXTUAL_PK"),
)


def _fraction_unbound(value: dict[str, Any] | None) -> dict[str, Any] | None:
    if not value or value.get("value") is None:
        return value
    unit = str(value.get("unit") or "").lower()
    if "%" not in unit and "bound" not in unit:
        return None
    transformed = dict(value)
    transformed["value"] = max(0.0, min(1.0, 1.0 - float(value["value"]) / 100.0))
    transformed["display_value"] = f"{transformed['value']:.4g}"
    transformed["unit"] = "fraction unbound"
    transformed["derived_from"] = "HUMAN_PPB"
    return transformed


def _select_row(rows: list[dict[str, Any]], endpoint: ProfileEndpoint) -> dict[str, Any] | None:
    source = endpoint.source_endpoint or endpoint.endpoint
    matches = [row for row in rows if row["canonical_endpoint"] == source]
    preferred = [row for row in matches if row.get("species") == endpoint.species]
    candidates = preferred or matches
    if not candidates:
        return None
    return next((row for row in candidates if row.get("prediction") and row.get("experimental")),
                next((row for row in candidates if row.get("prediction")), candidates[0]))


def _route_availability(endpoint_id: str, route: dict[str, Any] | None, declared: str | None) -> str:
    if declared:
        return declared
    if not route or route.get("route") == ROUTE_MODEL_UNAVAILABLE:
        return MODEL_UNAVAILABLE
    if route.get("route") == ROUTE_RETAIN_V3_3:
        return MECHANISTIC_ONLY
    return ON_DEMAND if model_artifact_registration(endpoint_id) is not None else MODEL_NOT_REGISTERED


def _profile_row(
    endpoint: ProfileEndpoint,
    scientific_rows: list[dict[str, Any]],
    routes: dict[str, dict[str, Any]],
    execution_outcomes: dict[str, dict[str, str]],
) -> dict[str, Any]:
    source_endpoint = endpoint.source_endpoint or endpoint.endpoint
    stable_row = _select_row(scientific_rows, endpoint)
    route = routes.get(source_endpoint)
    prediction = dict(stable_row["prediction"]) if stable_row and stable_row.get("prediction") else None
    experimental = dict(stable_row["experimental"]) if stable_row and stable_row.get("experimental") else None
    if endpoint.semantic in {"FRACTION_UNBOUND", "UPSTREAM_FU"}:
        prediction = _fraction_unbound(prediction)
        experimental = _fraction_unbound(experimental)
    resolved_availability = _route_availability(source_endpoint, route, endpoint.availability)
    availability = (
        resolved_availability
        if endpoint.availability in {MECHANISTIC_ONLY, MODEL_UNAVAILABLE, MODEL_NOT_REGISTERED, CURRENT_DATA_CEILING, CONTEXT_REQUIRED}
        else AVAILABLE_CURRENT if prediction else resolved_availability
    )
    if prediction and availability == AVAILABLE_CURRENT:
        status = "EXPERIMENTAL_AVAILABLE" if experimental else "PREDICTED"
    elif experimental:
        status = "EXPERIMENTAL_AVAILABLE"
    else:
        status = (
            availability
        )
    outcome = execution_outcomes.get(source_endpoint)
    if not prediction and outcome and outcome.get("status") in {CALCULATED_BUT_NOT_ELIGIBLE, FAILED}:
        status = outcome["status"]
    comparison = stable_row.get("comparison") if stable_row else None
    difference = None
    if endpoint.semantic not in {"FRACTION_UNBOUND", "UPSTREAM_FU"} and comparison and comparison.get("numeric_pairable"):
        difference = {
            "fold_error": comparison.get("fold_error"),
            "prediction_experimental_ratio": comparison.get("prediction_experimental_ratio"),
            "signed_log_error": comparison.get("signed_log_error"),
            "unit": comparison.get("comparison_unit"),
        }
    definition = REGISTRY.get(source_endpoint)
    model_id = prediction.get("model_id") if prediction else (route or {}).get("model_or_ensemble")
    model_version = prediction.get("model_version") if prediction else (route or {}).get("model_version_hash")
    maturity = stable_row.get("maturity") if stable_row else None
    ad = prediction.get("applicability_domain") if prediction else ({"classification": (route or {}).get("ad_status")} if route else None)
    reason = (
        "Current admitted Stable Core prediction is available."
        if prediction else
        (
            "The release registry describes this endpoint, but no exact executable artifact is admitted by Stable Core."
            if status in {"CURRENT_DATA_CEILING", MODEL_NOT_REGISTERED} else
            (route or {}).get("known_limitation")
            or {
                MODEL_UNAVAILABLE: "No qualified executable model is registered for this endpoint.",
                MODEL_NOT_REGISTERED: "A route label exists, but no complete authoritative executable model registration can be proven.",
                CONTEXT_REQUIRED: "Dose, route, regimen, formulation, or other PK context is required.",
                MECHANISTIC_ONLY: "Only a qualified mechanistic route is available; no structure-only value is implied.",
                CURRENT_DATA_CEILING: "The current evidence ceiling does not permit a canonical structure-only prediction.",
                ON_DEMAND: "A qualified model can run after the explicit Predict action.",
            }[availability]
        )
    )
    if not prediction and outcome and outcome.get("reason"):
        reason = outcome["reason"]
    return {
        "canonical_endpoint": source_endpoint,
        "query_endpoint": endpoint.endpoint,
        "aliases": list(endpoint.aliases),
        "display_name": endpoint.display_name,
        "category": endpoint.group,
        "semantic": endpoint.semantic or (definition.domain if definition else endpoint.group),
        "species": endpoint.species,
        "importance": endpoint.importance,
        "prediction": prediction,
        "experimental": experimental,
        "unit": (prediction or experimental or {}).get("unit") or endpoint.unit or (definition.canonical_unit if definition else (route or {}).get("unit", "")),
        "difference": difference,
        "status": status,
        "availability": availability,
        "model_id": model_id,
        "model_version": model_version,
        "prediction_mode": prediction.get("mode") if prediction else None,
        "maturity": maturity,
        "AD": ad,
        "context": stable_row.get("context", {}) if stable_row else {},
        "availability_reason": reason,
        "stable_core_identity": {
            "snapshot_id": prediction.get("snapshot_id") if prediction else None,
            "observation_id": experimental.get("observation_id") if experimental else None,
            "source_endpoint": source_endpoint,
        },
    }


def build_developability_profile(db: Any, version_id: int) -> dict[str, Any]:
    from .models import PredictionRun

    stable = build_scientific_endpoint_rows(db, version_id)
    routes = {row["endpoint_id"]: row for row in get_current_production_routing()}
    latest_workflow = db.scalar(select(PredictionRun).where(
        PredictionRun.version_id == version_id,
        PredictionRun.stage == "prediction_workflow",
    ).order_by(PredictionRun.created_at.desc(), PredictionRun.id.desc()))
    workflow_output = dict(latest_workflow.outputs_json or {}) if latest_workflow else {}
    execution_outcomes: dict[str, dict[str, str]] = {}
    for row in workflow_output.get("current_publication") or []:
        if row.get("status") == CALCULATED_BUT_NOT_ELIGIBLE:
            execution_outcomes[str(row.get("endpoint"))] = {
                "status": CALCULATED_BUT_NOT_ELIGIBLE,
                "reason": str(row.get("reason") or "Stable Core admission rejected the calculated result."),
            }
    for endpoint_id, output in dict(workflow_output.get("v3_predictions") or {}).items():
        if output.get("execution_status") == "EXECUTION_FAILED":
            execution_outcomes[endpoint_id] = {
                "status": FAILED,
                "reason": str(output.get("reason") or "Qualified model execution failed."),
            }
    groups = {name: [] for name in (
        "physchem", "absorption", "distribution", "metabolic_stability",
        "cyp", "transporters", "safety", "metabolism", "pk",
    )}
    for endpoint in CORE_PROFILE:
        groups[endpoint.group].append(_profile_row(endpoint, stable["rows"], routes, execution_outcomes))
    # Groups may intentionally repeat a canonical value as a compact upstream
    # PK reference.  The capability catalog itself stays unique so clients and
    # Predict summaries never double-count that shared scientific identity.
    entries_by_endpoint: dict[str, dict[str, Any]] = {}
    for group in groups.values():
        for row in group:
            entries_by_endpoint.setdefault(row["query_endpoint"], row)
    entries = list(entries_by_endpoint.values())
    counts = {state: sum(row["availability"] == state for row in entries) for state in (
        AVAILABLE_CURRENT, ON_DEMAND, MECHANISTIC_ONLY, MODEL_UNAVAILABLE,
        MODEL_NOT_REGISTERED, CURRENT_DATA_CEILING, CONTEXT_REQUIRED,
    )}
    return {
        "contract": "CompoundDevelopabilityProfile/stable-core-v1.2",
        "authority": stable["contract"],
        "compound_version_id": version_id,
        "prediction_engine": get_current_production_engine_info(),
        "groups": groups,
        "availability_catalog": entries,
        "availability_summary": counts,
    }
