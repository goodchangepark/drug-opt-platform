"""
Stage 4D-6: Canonical Prediction Orchestrator
==============================================

The single authoritative runtime execution path for Save & Predict.

This module replaces the legacy pattern of:
  run_admet_predictions → admet_model_registry (active=True) → predict_endpoint

with a scientifically correct pipeline:

  PredictionOrchestrator.orchestrate(version)
    → standardize structure
    → read EndpointStrategyPolicy for each endpoint
    → build execution plan (core + shadow model IDs)
    → execute CORE models via admet_predictor.predict_endpoint
    → execute SHADOW/RESEARCH models via multimodel adapters
    → persist individual predictions (core + shadow)
    → evaluate AD, calibration, disagreement
    → apply endpoint-specific production strategy (policy-defined, not auto-averaged)
    → store authorized research/shadow outputs
    → freeze pre-experimental prediction evidence
    → persist provenance
    → return backward-compatible production output

Design principles:
  - CORE production value = M1 unless registry explicitly states otherwise
  - Shadow/secondary execution does NOT change the production value
  - Shadow model failure NEVER breaks CORE execution
  - MODEL_UNAVAILABLE preserved; never fabricated
  - Same-compound leakage protection: freeze is stored BEFORE any
    experimental result is linked
  - Pre-experimental freeze is immutable once written
"""

from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

from sqlalchemy import select
from sqlalchemy.orm import Session

from .admet import (
    ADMETEndpoint,
    ADMETConsensusPrediction,
    ADMETModelRegistry,
    ADMETPrediction,
    ADMETPredictionRun,
)
from .admet_predictor import MODEL_SPECS, MODEL_VERSION, model_files_available, predict_endpoint
from .endpoint_contracts import get_endpoint_contract
from .endpoint_strategy_registry import (
    EndpointStrategyPolicy,
    StrategyType,
    get_all_strategies,
    get_endpoint_strategy,
)
from .prediction_engine_v1_policy import (
    ENGINE_V1_POLICY_ID,
    ENGINE_V1_POLICY_VERSION,
    policy_hash as engine_v1_policy_hash,
)
from .multimodel import (
    ExecutionStatus,
    ModelExecutionPayload,
    get_model_adapter,
    initialize_default_adapters,
)

# Ensure adapters are initialized on import
initialize_default_adapters()

# ---------------------------------------------------------------------------
# Stage 4D-6 orchestrator version
# ---------------------------------------------------------------------------

ORCHESTRATOR_VERSION = "stage4d6-prediction-orchestrator-v1"
# The Stage 4D registry remains executable authority; this identifies its
# immutable Stage 4E-4 Engine-v1 snapshot in every new freeze.
POLICY_VERSION = f"{ENGINE_V1_POLICY_ID}@{ENGINE_V1_POLICY_VERSION}"
STANDARDIZER_VERSION = "CHEM_STANDARDIZER_V1"

# ---------------------------------------------------------------------------
# Shadow model adapter mapping
# Maps strategy-registry model_id → multimodel adapter model_id
# These are the models that exist as adapters but NOT in MODEL_SPECS
# ---------------------------------------------------------------------------

SHADOW_ADAPTER_MAP: Dict[str, str] = {
    # Solubility shadow models
    "esol_delaney_v1": "esol_delaney_v1",
    "rdkit_gbr_solubility_v1": "rdkit_gbr_solubility_v1",
    # Permeability (Caco-2) shadow
    "physchem_caco2_v1": "physchem_caco2_v1",
    # CYP3A4 inhibitor M2
    "morgan_cyp3a4_inh_v1": "morgan_cyp3a4_inh_v1",
    # hERG secondary
    "physchem_herg_v1": "physchem_herg_v1",
}

# Fixed blend weights for CYP3A4 (research/shadow only; never production)
CYP3A4_BLEND_W1 = 0.9578  # admetica_cyp_cyp3a4-inhibitor (M1 CORE)
CYP3A4_BLEND_W2 = 0.0422  # morgan_cyp3a4_inh_v1 (M2 SHADOW)


# ---------------------------------------------------------------------------
# Endpoint execution plan
# ---------------------------------------------------------------------------

@dataclass
class EndpointExecutionPlan:
    """Explicit per-endpoint runtime execution plan derived from strategy registry."""
    endpoint_name: str
    endpoint_id: str
    endpoint_contract_version: str
    production_strategy: str
    # Primary (CORE) model — registry DB model_id (integer)
    core_registry_model_id: Optional[int] = None
    core_model_key: str = ""
    core_model_name: str = ""
    core_model_version: str = ""
    core_adapter_key: Optional[str] = None   # key in MODEL_SPECS
    # Shadow/supporting models (adapter IDs from SHADOW_ADAPTER_MAP)
    shadow_adapter_ids: List[str] = field(default_factory=list)
    shadow_model_roles: Dict[str, str] = field(default_factory=dict)
    shadow_registry_model_ids: Dict[str, int] = field(default_factory=dict)  # adapter_id → db model id
    # Calibration metadata
    calibration_status: str = "RAW"
    decision_threshold: Optional[float] = None
    # Execution control
    is_available: bool = True
    unavailable_reason: str = ""
    canonical_unit: str = ""
    # Fixed blend (CYP3A4)
    fixed_blend_weights: Optional[Dict[str, float]] = None
    blend_is_research_only: bool = True


@dataclass
class IndividualModelResult:
    """Result of executing one model for one endpoint."""
    model_key: str        # adapter_id or MODEL_SPECS key
    model_name: str
    model_version: str
    model_role: str       # CORE / SHADOW / RESEARCH_ONLY / CALIBRATION_SUPPORTING / etc.
    db_model_id: Optional[int]
    endpoint_name: str
    endpoint_id: str
    canonical_unit: str
    execution_status: str  # SUCCESS / MODEL_UNAVAILABLE / RUNTIME_ERROR / etc.
    predicted_value: Optional[float] = None
    probability: Optional[float] = None
    predicted_class: Optional[str] = None
    applicability_domain: str = "UNKNOWN"
    confidence: str = "NOT_APPLICABLE"
    uncertainty: Optional[float] = None
    runtime_ms: float = 0.0
    canonical_smiles: str = ""
    error_message: Optional[str] = None
    raw_outputs: Dict[str, Any] = field(default_factory=dict)
    provenance: Dict[str, Any] = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())


@dataclass
class EndpointOrchestrationResult:
    """Complete result for one endpoint after executing all authorized models."""
    endpoint_name: str
    endpoint_id: str
    canonical_unit: str
    production_strategy: str
    # Production (CORE) output
    production_value: Optional[float]
    production_probability: Optional[float]
    production_class: Optional[str]
    production_model_key: str
    production_model_name: str
    production_ad: str
    production_confidence: str
    production_decision_threshold: Optional[float]
    production_execution_status: str
    # All individual model results
    all_results: List[IndividualModelResult] = field(default_factory=list)
    # Shadow/research aggregated results
    shadow_results: List[IndividualModelResult] = field(default_factory=list)
    # Research consensus (e.g. CYP3A4 fixed blend) — shadow only
    research_consensus: Optional[Dict[str, Any]] = None
    # Model disagreement
    model_disagreement: Optional[float] = None
    model_agreement_class: str = "SINGLE_MODEL"
    # Persisted ADMETPrediction row IDs
    persisted_prediction_ids: List[int] = field(default_factory=list)
    # Pre-experimental freeze ID
    freeze_id: Optional[str] = None
    # Timing
    total_runtime_ms: float = 0.0
    errors: List[str] = field(default_factory=list)


@dataclass
class OrchestratorResult:
    """Complete result from the prediction orchestrator for one compound version."""
    version_id: int
    canonical_smiles: str
    standardizer_version: str
    policy_version: str
    orchestrator_version: str
    run_id: int
    endpoint_results: List[EndpointOrchestrationResult] = field(default_factory=list)
    endpoint_statuses: List[Dict[str, Any]] = field(default_factory=list)
    unavailable: List[str] = field(default_factory=list)
    total_runtime_ms: float = 0.0
    status: str = "COMPLETE"
    message: str = ""
    timestamp: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    def as_legacy_dict(
        self,
        available_models: list,
        measurements: list,
        endpoint_names: dict,
        from_admet_prediction_out_fn,
        consensus_fn,
        db: Session,
    ) -> dict:
        """Produce the same dict structure as the legacy run_admet_predictions return."""
        selected_predictions: Dict[int, ADMETPrediction] = {}
        for ep_res in self.endpoint_results:
            for pid in ep_res.persisted_prediction_ids:
                pred = db.get(ADMETPrediction, pid)
                if pred is not None:
                    selected_predictions[pred.model_id] = pred

        predictions = []
        for model in available_models:
            if model.id in selected_predictions:
                predictions.append(from_admet_prediction_out_fn(
                    selected_predictions[model.id], measurements, endpoint_names
                ))

        return {
            "type": "Predicted",
            "run_id": self.run_id,
            "status": self.status,
            "message": self.message,
            "models_available": len(available_models),
            "cache_hit": False,
            "predictions": predictions,
            "consensus_predictions": [consensus_fn(row) for row in []],  # consensus handled separately
            "endpoint_statuses": self.endpoint_statuses,
            "unavailable": self.unavailable,
            "orchestrator": ORCHESTRATOR_VERSION,
        }


def _sha256(payload: str) -> str:
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _model_role_label(model_id: str, policy: EndpointStrategyPolicy) -> str:
    """Determine the explicit runtime role for a model given the policy."""
    if model_id in policy.primary_model_ids:
        return "CORE"
    role = policy.non_primary_model_roles.get(model_id, "")
    if role:
        return role
    if model_id in policy.shadow_model_ids:
        return "SHADOW"
    return "RESEARCH_ONLY"


def is_core_registry_model(model: ADMETModelRegistry) -> bool:
    """Return true only for an exact Stage 4D primary model identity.

    ``is_active`` describes runtime availability.  It cannot decide production
    selection because authorized shadow rows are executable too.
    """
    policy = get_endpoint_strategy(model.endpoint_name)
    if policy is None or policy.primary_strategy == StrategyType.MODEL_UNAVAILABLE:
        return False
    return (
        model.endpoint_name in MODEL_SPECS
        and model.model_version in policy.primary_model_versions
        and (model.provenance_json or {}).get("production_eligible") is not False
    )


def _build_execution_plan(
    endpoint_name: str,
    db: Session,
    policy: EndpointStrategyPolicy,
) -> EndpointExecutionPlan:
    """
    Build a runtime execution plan from the strategy policy.

    The plan specifies:
    - Which DB-registered CORE model to run
    - Which shadow adapter IDs to run
    - Which roles each model has
    - Whether the endpoint is available at all
    """
    plan = EndpointExecutionPlan(
        endpoint_name=endpoint_name,
        endpoint_id=policy.endpoint_id,
        endpoint_contract_version=policy.endpoint_contract_version,
        production_strategy=policy.primary_strategy.value,
        calibration_status=policy.calibration_status.value,
        decision_threshold=policy.decision_threshold,
    )

    contract = get_endpoint_contract(endpoint_name)
    if contract:
        plan.canonical_unit = contract.canonical_unit

    # Strategy MODEL_UNAVAILABLE → no execution
    if policy.primary_strategy == StrategyType.MODEL_UNAVAILABLE:
        plan.is_available = False
        plan.unavailable_reason = f"MODEL_UNAVAILABLE: {endpoint_name}"
        return plan

    # Derived/Rule-based strategies that have no ML model
    if policy.primary_strategy in (
        StrategyType.DERIVED_ESTIMATE,
        StrategyType.RULE_ESTIMATE,
        StrategyType.MECHANISTIC_NO_CONSENSUS,
        StrategyType.RANK_FUSION,
    ):
        plan.is_available = False
        plan.unavailable_reason = f"Non-ML strategy: {policy.primary_strategy.value}"
        return plan

    # Find CORE model in DB registry
    if endpoint_name in MODEL_SPECS:
        # Try to find a registered active model for this endpoint
        core_model = db.scalar(
            select(ADMETModelRegistry).where(
                ADMETModelRegistry.endpoint_name == endpoint_name,
                ADMETModelRegistry.is_active.is_(True),
                ADMETModelRegistry.model_version.in_(policy.primary_model_versions),
            )
        )
        if core_model is None:
            plan.is_available = False
            plan.unavailable_reason = f"No active model in registry for {endpoint_name}"
            return plan

        available, reason = model_files_available(endpoint_name)
        if not available:
            plan.is_available = False
            plan.unavailable_reason = reason
            return plan

        plan.core_registry_model_id = core_model.id
        plan.core_model_key = policy.primary_model_ids[0] if policy.primary_model_ids else endpoint_name
        plan.core_model_name = core_model.model_name
        plan.core_model_version = core_model.model_version
        plan.core_adapter_key = endpoint_name

    # Find shadow model registry IDs (pre-registered shadow model rows)
    for shadow_id in policy.shadow_model_ids:
        role = policy.non_primary_model_roles.get(shadow_id, "SHADOW")
        # M3 remains auditable in the registry but is excluded by the Stage
        # 4D Solubility policy from adaptive/consensus runtime execution.
        if "EXCLUDED" in role.upper():
            continue
        adapter_key = SHADOW_ADAPTER_MAP.get(shadow_id)
        if adapter_key is None:
            continue
        adapter = get_model_adapter(adapter_key)
        if adapter is None:
            continue
        # Look up shadow model in DB
        shadow_reg = db.scalar(
            select(ADMETModelRegistry).where(
                ADMETModelRegistry.endpoint_name == endpoint_name,
                ADMETModelRegistry.model_version == adapter.model_version,
                ADMETModelRegistry.model_name == adapter.model_name,
            )
        )
        plan.shadow_adapter_ids.append(shadow_id)
        plan.shadow_model_roles[shadow_id] = role
        if shadow_reg is not None:
            plan.shadow_registry_model_ids[shadow_id] = shadow_reg.id

    # CYP3A4 fixed blend (research only)
    if endpoint_name == "CYP3A4 inhibitor" and policy.shadow_strategy == StrategyType.FIXED_WEIGHT_BLEND:
        plan.fixed_blend_weights = {
            "admetica_cyp_cyp3a4-inhibitor": CYP3A4_BLEND_W1,
            "morgan_cyp3a4_inh_v1": CYP3A4_BLEND_W2,
        }
        plan.blend_is_research_only = True

    return plan


def _execute_core_model(
    smiles: str,
    plan: EndpointExecutionPlan,
) -> IndividualModelResult:
    """Execute the primary CORE model via admet_predictor.predict_endpoint."""
    t0 = time.perf_counter()
    endpoint_name = plan.endpoint_name
    try:
        result = predict_endpoint(smiles, endpoint_name)
    except Exception as exc:
        return IndividualModelResult(
            model_key=plan.core_model_key or endpoint_name,
            model_name=plan.core_model_name,
            model_version=plan.core_model_version,
            model_role="CORE",
            db_model_id=plan.core_registry_model_id,
            endpoint_name=endpoint_name,
            endpoint_id=plan.endpoint_id,
            canonical_unit=plan.canonical_unit,
            execution_status="RUNTIME_ERROR",
            error_message=str(exc),
            runtime_ms=round((time.perf_counter() - t0) * 1000.0, 2),
            canonical_smiles=smiles,
        )

    if result.get("status") != "COMPLETE":
        return IndividualModelResult(
            model_key=plan.core_model_key or endpoint_name,
            model_name=plan.core_model_name,
            model_version=plan.core_model_version,
            model_role="CORE",
            db_model_id=plan.core_registry_model_id,
            endpoint_name=endpoint_name,
            endpoint_id=plan.endpoint_id,
            canonical_unit=plan.canonical_unit,
            execution_status="MODEL_UNAVAILABLE",
            error_message=result.get("reason", "model unavailable"),
            runtime_ms=round((time.perf_counter() - t0) * 1000.0, 2),
            canonical_smiles=smiles,
        )

    domain = result.get("applicability_domain", {})
    ad_class = domain.get("classification", "UNKNOWN") if isinstance(domain, dict) else str(domain)
    spec = MODEL_SPECS.get(endpoint_name, {})
    prov = {
        "model_source": spec.get("source", ""),
        "endpoint_definition": spec.get("endpoint_definition", ""),
        "training_dataset": spec.get("training_dataset", ""),
        "validation": spec.get("validation", {}),
        "license": spec.get("license", ""),
        "limitations": spec.get("limitations", ""),
        "applicability_domain_details": domain,
        "uncertainty_reason": result.get("uncertainty_reason", ""),
    }
    for key in ("assay_definition", "training_n", "independent_validation"):
        if spec.get(key) is not None:
            prov[key] = spec[key]
    for key in ("probability", "classification", "isoform", "transporter", "safety_endpoint",
                "species", "role", "decision_threshold", "liability_summary", "ensemble_probabilities"):
        if result.get(key) is not None:
            prov[key] = result[key]
    if result.get("derived_outputs") is not None:
        prov["derived_outputs"] = result["derived_outputs"]
    if result.get("metabolic_stability_assessment") is not None:
        prov["metabolic_stability_assessment"] = result["metabolic_stability_assessment"]
    if result.get("calibrated_uncertainty") is not None:
        prov["calibrated_uncertainty"] = result["calibrated_uncertainty"]

    return IndividualModelResult(
        model_key=plan.core_model_key or endpoint_name,
        model_name=plan.core_model_name,
        model_version=plan.core_model_version,
        model_role="CORE",
        db_model_id=plan.core_registry_model_id,
        endpoint_name=endpoint_name,
        endpoint_id=plan.endpoint_id,
        canonical_unit=result.get("unit", plan.canonical_unit),
        execution_status="SUCCESS",
        predicted_value=result.get("predicted_value"),
        probability=result.get("probability"),
        predicted_class=result.get("classification"),
        applicability_domain=ad_class,
        confidence=result.get("confidence", "UNKNOWN"),
        uncertainty=result.get("uncertainty"),
        runtime_ms=result.get("runtime_ms", round((time.perf_counter() - t0) * 1000.0, 2)),
        canonical_smiles=smiles,
        raw_outputs=prov,
    )


def _execute_shadow_model(
    smiles: str,
    shadow_id: str,
    plan: EndpointExecutionPlan,
) -> IndividualModelResult:
    """Execute a shadow/secondary model via its multimodel adapter. Never raises."""
    t0 = time.perf_counter()
    adapter_key = SHADOW_ADAPTER_MAP.get(shadow_id)
    role = plan.shadow_model_roles.get(shadow_id, "SHADOW")
    db_model_id = plan.shadow_registry_model_ids.get(shadow_id)

    if adapter_key is None:
        return IndividualModelResult(
            model_key=shadow_id,
            model_name=shadow_id,
            model_version="unknown",
            model_role=role,
            db_model_id=db_model_id,
            endpoint_name=plan.endpoint_name,
            endpoint_id=plan.endpoint_id,
            canonical_unit=plan.canonical_unit,
            execution_status="SKIPPED",
            error_message=f"No adapter mapping for {shadow_id}",
            runtime_ms=0.0,
            canonical_smiles=smiles,
        )

    adapter = get_model_adapter(adapter_key)
    if adapter is None:
        return IndividualModelResult(
            model_key=shadow_id,
            model_name=shadow_id,
            model_version="unknown",
            model_role=role,
            db_model_id=db_model_id,
            endpoint_name=plan.endpoint_name,
            endpoint_id=plan.endpoint_id,
            canonical_unit=plan.canonical_unit,
            execution_status="MODEL_UNAVAILABLE",
            error_message=f"Adapter {shadow_id} not registered",
            runtime_ms=0.0,
            canonical_smiles=smiles,
        )

    try:
        available, reason = adapter.is_available()
        if not available:
            return IndividualModelResult(
                model_key=shadow_id,
                model_name=adapter.model_name,
                model_version=adapter.model_version,
                model_role=role,
                db_model_id=db_model_id,
                endpoint_name=plan.endpoint_name,
                endpoint_id=plan.endpoint_id,
                canonical_unit=plan.canonical_unit,
                execution_status="MODEL_UNAVAILABLE",
                error_message=reason,
                runtime_ms=0.0,
                canonical_smiles=smiles,
            )

        contract = get_endpoint_contract(plan.endpoint_name)
        payload: ModelExecutionPayload = adapter.execute(smiles, contract)
        runtime_ms = round((time.perf_counter() - t0) * 1000.0, 2)

        if payload.execution_status == ExecutionStatus.SUCCESS:
            return IndividualModelResult(
                model_key=shadow_id,
                model_name=adapter.model_name,
                model_version=adapter.model_version,
                model_role=role,
                db_model_id=db_model_id,
                endpoint_name=plan.endpoint_name,
                endpoint_id=plan.endpoint_id,
                canonical_unit=payload.canonical_unit or plan.canonical_unit,
                execution_status="SUCCESS",
                predicted_value=payload.value,
                probability=payload.probability,
                predicted_class=payload.predicted_class,
                applicability_domain=payload.applicability_domain,
                confidence=payload.confidence,
                runtime_ms=payload.runtime_ms or runtime_ms,
                canonical_smiles=smiles,
                raw_outputs={
                    "model_role": role,
                    "model_id": shadow_id,
                    "model_family": adapter.model_family,
                    "raw_outputs": payload.raw_outputs,
                    "provenance": payload.provenance,
                    "warnings": payload.warnings,
                },
            )
        else:
            return IndividualModelResult(
                model_key=shadow_id,
                model_name=adapter.model_name,
                model_version=adapter.model_version,
                model_role=role,
                db_model_id=db_model_id,
                endpoint_name=plan.endpoint_name,
                endpoint_id=plan.endpoint_id,
                canonical_unit=plan.canonical_unit,
                execution_status=payload.execution_status.value,
                error_message=payload.error_message,
                runtime_ms=payload.runtime_ms or runtime_ms,
                canonical_smiles=smiles,
            )
    except Exception as exc:
        # Shadow model failure MUST NOT break CORE
        return IndividualModelResult(
            model_key=shadow_id,
            model_name=adapter_key,
            model_version="unknown",
            model_role=role,
            db_model_id=db_model_id,
            endpoint_name=plan.endpoint_name,
            endpoint_id=plan.endpoint_id,
            canonical_unit=plan.canonical_unit,
            execution_status="RUNTIME_ERROR",
            error_message=f"{type(exc).__name__}: {exc}",
            runtime_ms=round((time.perf_counter() - t0) * 1000.0, 2),
            canonical_smiles=smiles,
        )


def _compute_disagreement(
    results: List[IndividualModelResult],
    is_classification: bool,
    positive_label: str = "",
    threshold: float = 0.5,
) -> Tuple[Optional[float], str]:
    """Compute model disagreement across successful individual results."""
    values = [r.predicted_value for r in results
              if r.execution_status == "SUCCESS" and r.predicted_value is not None]
    if len(values) < 2:
        return None, "SINGLE_MODEL"
    if is_classification:
        classes = [("POSITIVE" if v >= threshold else "NEGATIVE") for v in values]
        pos_count = sum(1 for c in classes if c == (positive_label or "POSITIVE"))
        agreement = "HIGH_AGREEMENT" if pos_count == 0 or pos_count == len(classes) else "MODERATE_AGREEMENT"
        std = float((max(values) - min(values)) / 2.0)
        return std, agreement
    else:
        import statistics
        if len(values) < 2:
            return 0.0, "SINGLE_MODEL"
        std = statistics.stdev(values)
        if std <= 0.30:
            agreement = "HIGH_AGREEMENT"
        elif std <= 0.60:
            agreement = "MODERATE_AGREEMENT"
        else:
            agreement = "LOW_AGREEMENT"
        return std, agreement


def _persist_prediction(
    db: Session,
    run_id: int,
    version_id: int,
    project_id: int,
    result: IndividualModelResult,
    endpoint_db_id: int,
) -> Optional[int]:
    """
    Persist one model result to admet_predictions.
    Returns the new prediction row id, or None if not persisted.
    Failures are retained as explicit provenance rows when their registry
    identity exists; they can never overwrite a successful CORE cache entry.
    """
    if result.db_model_id is None:
        return None
    # Check for duplicate (cache)
    existing = db.scalar(
        select(ADMETPrediction).where(
            ADMETPrediction.version_id == version_id,
            ADMETPrediction.model_id == result.db_model_id,
            ADMETPrediction.execution_status == result.execution_status,
        ).order_by(ADMETPrediction.created_at.desc())
    )
    if existing is not None:
        # Check if already cached for this version+model
        if existing.predicted_value == result.predicted_value:
            return existing.id

    outputs_json = dict(result.raw_outputs)
    outputs_json["model_role"] = result.model_role
    outputs_json["stage4d6_orchestrator"] = ORCHESTRATOR_VERSION
    outputs_json["record_type"] = "Predicted"
    outputs_json["compound_version_id"] = version_id
    outputs_json["prediction_timestamp"] = result.timestamp
    if result.error_message:
        outputs_json["error_message"] = result.error_message

    pred = ADMETPrediction(
        run_id=run_id,
        endpoint_id=endpoint_db_id,
        version_id=version_id,
        model_id=result.db_model_id,
        model_version=result.model_version,
        execution_status=result.execution_status,
        standardizer_version=STANDARDIZER_VERSION,
        canonical_smiles=result.canonical_smiles,
        runtime_ms=result.runtime_ms,
        predicted_value=result.predicted_value,
        unit=result.canonical_unit,
        confidence=result.confidence,
        applicability_domain=result.applicability_domain,
        uncertainty=result.uncertainty,
        outputs_json=outputs_json,
    )
    db.add(pred)
    db.flush()
    return pred.id


def _get_or_create_endpoint(db: Session, project_id: int, endpoint_name: str) -> ADMETEndpoint:
    """Get or create an ADMETEndpoint row for this project."""
    ep = db.scalar(
        select(ADMETEndpoint).where(
            ADMETEndpoint.project_id == project_id,
            ADMETEndpoint.name == endpoint_name,
        )
    )
    if ep is None:
        ep = ADMETEndpoint(
            project_id=project_id,
            name=endpoint_name,
            category="ADME",
        )
        db.add(ep)
        db.flush()
    return ep


def _compute_research_consensus(
    core_result: IndividualModelResult,
    shadow_results: List[IndividualModelResult],
    plan: EndpointExecutionPlan,
) -> Optional[Dict[str, Any]]:
    """
    Compute research-only consensus where policy defines one.
    For CYP3A4: compute fixed blend (research/shadow only; NOT production).
    """
    if plan.endpoint_name != "CYP3A4 inhibitor":
        return None
    if plan.fixed_blend_weights is None:
        return None
    if core_result.execution_status != "SUCCESS":
        return None
    # Find M2 result
    m2_result = next(
        (r for r in shadow_results
         if r.model_key == "morgan_cyp3a4_inh_v1" and r.execution_status == "SUCCESS"),
        None,
    )
    if m2_result is None:
        return None
    w1 = CYP3A4_BLEND_W1
    w2 = CYP3A4_BLEND_W2
    v1 = core_result.predicted_value or 0.0
    v2 = m2_result.predicted_value or 0.0
    blend = w1 * v1 + w2 * v2
    return {
        "blend_type": "FIXED_WEIGHT_BLEND_RESEARCH_SHADOW_ONLY",
        "weights": {"M1_admetica": w1, "M2_morgan": w2},
        "M1_value": v1,
        "M2_value": v2,
        "blend_value": round(blend, 6),
        "production_status": "RESEARCH_SHADOW_NOT_PRODUCTION",
        "policy": "Stage 4D-3B1A: fixed blend remains research/shadow; dynamic adaptation NO_ADAPTIVE_VALUE",
    }


def _persist_research_consensus(
    db: Session,
    run_id: int,
    version_id: int,
    endpoint_db_id: int,
    core_result: IndividualModelResult,
    plan: EndpointExecutionPlan,
    research_consensus: Optional[Dict[str, Any]],
    disagreement: Optional[float],
    agreement_class: str,
) -> None:
    """Persist a clearly marked research-only blend without touching CORE."""
    if research_consensus is None:
        return
    existing = db.scalar(
        select(ADMETConsensusPrediction).where(
            ADMETConsensusPrediction.version_id == version_id,
            ADMETConsensusPrediction.endpoint_id == endpoint_db_id,
            ADMETConsensusPrediction.consensus_version == ORCHESTRATOR_VERSION,
        )
    )
    if existing is not None:
        return
    value = research_consensus.get("blend_value")
    db.add(ADMETConsensusPrediction(
        run_id=run_id,
        endpoint_id=endpoint_db_id,
        version_id=version_id,
        consensus_version=ORCHESTRATOR_VERSION,
        consensus_mode="SHADOW",
        combined_value=value,
        unit=core_result.canonical_unit,
        classification=("Inhibitor" if value is not None and value >= (plan.decision_threshold or 0.5) else "Non-inhibitor"),
        confidence=core_result.confidence,
        applicability_domain=core_result.applicability_domain,
        model_agreement=agreement_class,
        dispersion_json={"model_disagreement": disagreement},
        vote_pattern="RESEARCH_FIXED_BLEND",
        weights_json=[
            {"model_id": "admetica_cyp_cyp3a4-inhibitor", "weight": CYP3A4_BLEND_W1},
            {"model_id": "morgan_cyp3a4_inh_v1", "weight": CYP3A4_BLEND_W2},
        ],
        provenance_json={
            **research_consensus,
            "stage4d6_research_only": True,
            "production_value_unchanged": core_result.predicted_value,
            "production_model_role": "CORE",
        },
    ))
    db.flush()


def _freeze_prediction(
    db: Session,
    version_id: int,
    project_id: int,
    smiles: str,
    plan: EndpointExecutionPlan,
    core_result: IndividualModelResult,
    shadow_results: List[IndividualModelResult],
    research_consensus: Optional[Dict[str, Any]],
    ep_result: "EndpointOrchestrationResult",
) -> Optional[str]:
    """
    Freeze the pre-experimental prediction evidence.
    Returns the frozen_prediction_id if successful.
    This is idempotent: if a freeze already exists for this
    (version_id, endpoint_id, smiles hash), it is reused.
    """
    from .production_qualification import QualificationPredictionFreezeRow, StrategyType as QualStrategyType

    endpoint_id = plan.endpoint_id
    smiles_hash = _sha256(smiles)
    # Engine-v1 has a distinct namespace so a re-execution cannot collide with
    # or mutate a historical Stage 4D freeze for the same compound.
    frozen_id = f"FREEZE-V1-{version_id}-{endpoint_id}-{smiles_hash[:12]}"

    # Idempotency: don't duplicate
    existing = db.scalar(
        select(QualificationPredictionFreezeRow).where(
            QualificationPredictionFreezeRow.frozen_prediction_id == frozen_id
        )
    )
    if existing is not None:
        return frozen_id

    # The qualification store's models_json is a strict immutable identity
    # schema. Per-execution values remain in provenance, so later linkage can
    # rebuild a real FrozenModelIdentity without guessing a checkpoint.
    model_identities = []
    model_executions = []
    if core_result.execution_status == "SUCCESS":
        model_identities.append({
            "model_id": core_result.model_key,
            "model_version": core_result.model_version,
            "checkpoint_hash": _sha256(
                f"{core_result.model_key}|{core_result.model_version}|{POLICY_VERSION}"
            ),
        })
        model_executions.append({
            "model_id": core_result.model_key,
            "model_name": core_result.model_name,
            "model_version": core_result.model_version,
            "role": "CORE",
            "predicted_value": core_result.predicted_value,
            "probability": core_result.probability,
            "applicability_domain": core_result.applicability_domain,
            "confidence": core_result.confidence,
        })
    for sr in shadow_results:
        model_identities.append({
            "model_id": sr.model_key,
            "model_version": sr.model_version,
            "checkpoint_hash": _sha256(f"{sr.model_key}|{sr.model_version}|{POLICY_VERSION}"),
        })
        model_executions.append({
            "model_id": sr.model_key,
            "model_name": sr.model_name,
            "model_version": sr.model_version,
            "role": sr.model_role,
            "execution_status": sr.execution_status,
            "predicted_value": sr.predicted_value,
            "probability": sr.probability,
            "applicability_domain": sr.applicability_domain,
            "confidence": sr.confidence,
        })

    provenance = {
        "orchestrator": ORCHESTRATOR_VERSION,
        "policy_version": POLICY_VERSION,
        "engine_v1_policy_id": ENGINE_V1_POLICY_ID,
        "engine_v1_policy_hash": engine_v1_policy_hash(),
        "production_strategy": plan.production_strategy,
        "smiles_hash": smiles_hash,
        "research_consensus": research_consensus,
        "model_disagreement": ep_result.model_disagreement,
        "model_agreement_class": ep_result.model_agreement_class,
        "individual_predictions": model_executions,
    }

    # Map strategy string to enum
    strategy_map = {
        "SINGLE_CORE_MODEL": QualStrategyType.SINGLE_CORE_MODEL,
        "FIXED_WEIGHT_BLEND": QualStrategyType.FIXED_WEIGHT_BLEND,
        "STATIC_CONSENSUS": QualStrategyType.STATIC_CONSENSUS,
        "SINGLE_CORE_WITH_CALIBRATION": QualStrategyType.SINGLE_CORE_WITH_CALIBRATION,
    }
    strategy_enum = strategy_map.get(plan.production_strategy, QualStrategyType.SINGLE_CORE_MODEL)

    frozen_at = datetime.now(timezone.utc)
    record_data = {
        "frozen_prediction_id": frozen_id,
        "compound_version_id": str(version_id),
        "project_id": str(project_id),
        "chemical_series_id": "",
        "endpoint_id": endpoint_id,
        "endpoint_contract_version": plan.endpoint_contract_version,
        "candidate_id": f"{endpoint_id}-production",
        "strategy": strategy_enum,
        "models": model_identities,
        "prediction_value": ep_result.production_value,
        "probability": ep_result.production_probability,
        "unit": ep_result.canonical_unit,
        "policy_version": POLICY_VERSION,
        "engine_v1_policy_hash": engine_v1_policy_hash(),
        "standardizer_version": STANDARDIZER_VERSION,
        "applicability_domain": ep_result.production_ad,
        "provenance": provenance,
        "frozen_at": frozen_at,
    }
    # We create the freeze row directly to avoid the full ProspectivePredictionFreeze dataclass dependency
    record_hash_payload = json.dumps(record_data, sort_keys=True, default=str)
    record_hash = _sha256(record_hash_payload)

    row = QualificationPredictionFreezeRow(
        frozen_prediction_id=frozen_id,
        compound_version_id=str(version_id),
        project_id=str(project_id),
        chemical_series_id="",
        endpoint_id=endpoint_id,
        endpoint_contract_version=plan.endpoint_contract_version,
        candidate_id=f"{endpoint_id}-production",
        candidate_specification_hash=_sha256(json.dumps({
            "endpoint_id": endpoint_id,
            "endpoint_contract_version": plan.endpoint_contract_version,
            "production_strategy": plan.production_strategy,
            "models": model_identities,
            "policy_version": POLICY_VERSION,
            "engine_v1_policy_hash": engine_v1_policy_hash(),
            "standardizer_version": STANDARDIZER_VERSION,
        }, sort_keys=True)),
        strategy=strategy_enum.value,
        models_json=model_identities,
        prediction_value=ep_result.production_value,
        probability=ep_result.production_probability,
        unit=ep_result.canonical_unit,
        frozen_at=frozen_at,
        policy_version=POLICY_VERSION,
        standardizer_version=STANDARDIZER_VERSION,
        applicability_domain=ep_result.production_ad,
        provenance_json=provenance,
        record_hash=record_hash,
    )
    try:
        with db.begin_nested():
            db.add(row)
            db.flush()
        return frozen_id
    except Exception:
        # Freeze evidence must never roll back successfully persisted CORE or
        # shadow results. The outer prediction transaction remains intact.
        return None


class PredictionOrchestrator:
    """
    The single authoritative prediction execution path for Stage 4D-6.

    Usage:
        orchestrator = PredictionOrchestrator(db, version, compound)
        result = orchestrator.orchestrate()
        # result.as_legacy_dict(...) for backward compatibility
    """

    def __init__(self, db: Session, version: Any, compound: Any):
        self.db = db
        self.version = version
        self.compound = compound
        self.project_id = compound.project_id
        self.smiles = version.canonical_smiles
        self.version_id = version.id

    def orchestrate(self) -> OrchestratorResult:
        """Execute all authorized models for all endpoints and persist results."""
        t_total = time.perf_counter()
        db = self.db

        # Build the prediction run record
        digest = hashlib.sha256(
            f"{self.version_id}|{self.smiles}|{MODEL_VERSION}|{ORCHESTRATOR_VERSION}".encode()
        ).hexdigest()
        run = ADMETPredictionRun(
            version_id=self.version_id,
            inputs_hash=digest,
            status="RUNNING",
            message="Stage 4D-6 multimodel prediction orchestrator executing.",
        )
        db.add(run)
        db.flush()

        result = OrchestratorResult(
            version_id=self.version_id,
            canonical_smiles=self.smiles,
            standardizer_version=STANDARDIZER_VERSION,
            policy_version=POLICY_VERSION,
            orchestrator_version=ORCHESTRATOR_VERSION,
            run_id=run.id,
        )

        # The endpoint strategy registry is authoritative.  An active DB row
        # only declares a locally executable implementation and cannot add or
        # replace a production CORE path.
        strategy_endpoints = {
            name for name, policy in get_all_strategies().items()
            if name in MODEL_SPECS and policy.primary_strategy != StrategyType.MODEL_UNAVAILABLE
        }

        endpoint_results: List[EndpointOrchestrationResult] = []
        endpoint_statuses: List[Dict[str, Any]] = []
        unavailable: List[str] = []

        for endpoint_name in sorted(strategy_endpoints):
            policy = get_endpoint_strategy(endpoint_name)
            if policy is None:
                continue

            t_ep = time.perf_counter()
            plan = _build_execution_plan(endpoint_name, db, policy)

            if not plan.is_available:
                unavailable.append(f"{endpoint_name}: {plan.unavailable_reason}")
                endpoint_statuses.append({
                    "endpoint": endpoint_name,
                    "status": "MODEL_UNAVAILABLE",
                    "message": plan.unavailable_reason,
                })
                continue

            # Get or create DB endpoint
            ep_obj = _get_or_create_endpoint(db, self.project_id, endpoint_name)

            # ── Execute CORE model ──────────────────────────────────────────
            core_result = _execute_core_model(self.smiles, plan)

            if core_result.execution_status != "SUCCESS":
                unavailable.append(f"{endpoint_name}: {core_result.error_message}")
                endpoint_statuses.append({
                    "endpoint": endpoint_name,
                    "model_id": plan.core_registry_model_id,
                    "status": core_result.execution_status,
                    "message": core_result.error_message,
                    "cache_hit": False,
                })
                continue

            # ── Execute SHADOW/secondary models ──────────────────────────────
            # Shadow model failures are silently isolated; never break CORE
            shadow_results: List[IndividualModelResult] = []
            for shadow_id in plan.shadow_adapter_ids:
                sr = _execute_shadow_model(self.smiles, shadow_id, plan)
                shadow_results.append(sr)

            # ── Compute disagreement ──────────────────────────────────────────
            all_results = [core_result] + shadow_results
            successful_results = [r for r in all_results if r.execution_status == "SUCCESS"]
            spec = MODEL_SPECS.get(endpoint_name, {})
            is_classification = spec.get("prediction_type") == "binary_classification"
            positive_label = spec.get("positive_label", "POSITIVE")
            threshold = plan.decision_threshold or 0.5
            disagreement, agreement_class = _compute_disagreement(
                successful_results, is_classification, positive_label, threshold
            )

            # ── Research consensus (CYP3A4 fixed blend) ──────────────────────
            research_consensus = _compute_research_consensus(core_result, shadow_results, plan)

            # ── Build endpoint result ─────────────────────────────────────────
            ep_result = EndpointOrchestrationResult(
                endpoint_name=endpoint_name,
                endpoint_id=plan.endpoint_id,
                canonical_unit=core_result.canonical_unit,
                production_strategy=plan.production_strategy,
                production_value=core_result.predicted_value,
                production_probability=core_result.probability,
                production_class=core_result.predicted_class,
                production_model_key=core_result.model_key,
                production_model_name=core_result.model_name,
                production_ad=core_result.applicability_domain,
                production_confidence=core_result.confidence,
                production_decision_threshold=plan.decision_threshold,
                production_execution_status=core_result.execution_status,
                all_results=all_results,
                shadow_results=shadow_results,
                research_consensus=research_consensus,
                model_disagreement=disagreement,
                model_agreement_class=agreement_class,
                total_runtime_ms=round((time.perf_counter() - t_ep) * 1000.0, 2),
            )

            # ── Persist CORE prediction ──────────────────────────────────────
            core_pred_id = _persist_prediction(
                db, run.id, self.version_id, self.project_id, core_result, ep_obj.id
            )
            if core_pred_id is not None:
                ep_result.persisted_prediction_ids.append(core_pred_id)

            # Shadow executions are persisted inside the immutable
            # pre-experimental freeze below, rather than as primary
            # ADMETPrediction rows. This keeps legacy primary-model queries
            # unambiguous even when an installed shadow implementation has an
            # active registry row.

            _persist_research_consensus(
                db, run.id, self.version_id, ep_obj.id, core_result, plan,
                research_consensus, disagreement, agreement_class,
            )

            # ── Pre-experimental freeze ──────────────────────────────────────
            try:
                freeze_id = _freeze_prediction(
                    db, self.version_id, self.project_id, self.smiles,
                    plan, core_result, shadow_results, research_consensus, ep_result,
                )
                ep_result.freeze_id = freeze_id
            except Exception:
                pass  # Freeze failure never breaks prediction

            endpoint_results.append(ep_result)
            endpoint_statuses.append({
                "endpoint": endpoint_name,
                "model_id": plan.core_registry_model_id,
                "status": "COMPLETE",
                "cache_hit": False,
                "shadow_models_executed": len([sr for sr in shadow_results if sr.execution_status == "SUCCESS"]),
                "total_models_executed": len(successful_results),
            })

        # ── Finalize run ─────────────────────────────────────────────────────
        total_ms = round((time.perf_counter() - t_total) * 1000.0, 2)
        run.completed_at = datetime.now(timezone.utc)
        if endpoint_results and unavailable:
            run.status = "PARTIAL"
            run.message = f"Stage 4D-6 multimodel: partial. {len(unavailable)} endpoint(s) unavailable."
        elif endpoint_results:
            run.status = "COMPLETE"
            run.message = f"Stage 4D-6 multimodel: {len(endpoint_results)} endpoints completed."
        else:
            run.status = "MODEL_UNAVAILABLE"
            run.message = "; ".join(unavailable[:3]) or "No ADMET models available."
        db.flush()

        result.endpoint_results = endpoint_results
        result.endpoint_statuses = endpoint_statuses
        result.unavailable = unavailable
        result.total_runtime_ms = total_ms
        result.status = run.status
        result.message = run.message
        return result


# ---------------------------------------------------------------------------
# Shadow model registry seed
# Called from admet.py ensure_admet_schema to register shadow model rows
# ---------------------------------------------------------------------------

SHADOW_MODEL_SEEDS = [
    # Solubility M2 — ESOL
    {
        "endpoint_name": "Solubility",
        "model_name": "Delaney ESOL Physicochemical Regressor",
        "model_version": "esol-delaney-2004-v1.0",
        "implementation_status": "READY",
        "supported_species": ["Chemical"],
        "supported_matrix": ["Aqueous"],
        "output_unit": "log10(mol/L)",
        "source": "Delaney JS (2004) J Chem Inf Comput Sci 44:1000-1005; RDKit 2D descriptors",
        "training_dataset": "Delaney (2004) 1,144 aqueous solubility measurements",
        "validation_json": {
            "note": "Physicochemical rule-based regressor; no ML checkpoint. SHADOW_RESEARCH only.",
            "model_role": "SHADOW_RESEARCH",
        },
        "license": "Open Scientific Literature",
        "model_priority": 200,
        "ensemble_eligible": False,
        "species": "Chemical",
        "output_type": "regression_logS",
        "provenance_json": {
            "model_id": "esol_delaney_v1",
            "model_role": "SHADOW_RESEARCH",
            "stage4d_classification": "SHADOW / RESEARCH_SUPPORTING",
            "consensus_permission": "INSUFFICIENT_EVIDENCE",
            "production_eligible": False,
            "reason": "Stage 4D-3A2: adaptive architecture valid but no accuracy gain vs M1 on held-out validation",
        },
        "is_active": True,
    },
    # Solubility M3 — RDKit GBR (ADAPTIVE_EXCLUDED in production)
    {
        "endpoint_name": "Solubility",
        "model_name": "RDKit Descriptor GBR Intrinsic Solubility",
        "model_version": "rdkit-gbr-sol-v1.0",
        "implementation_status": "READY",
        "supported_species": ["Chemical"],
        "supported_matrix": ["Aqueous"],
        "output_unit": "log10(mol/L)",
        "source": "RDKit 2D topological descriptors + Gradient Boosting Regressor",
        "training_dataset": "AqSolDB subset; pure intrinsic solubility only",
        "validation_json": {
            "note": "Physicochemical GBR; no separate held-out validation. ADAPTIVE_EXCLUDED.",
            "model_role": "ADAPTIVE_EXCLUDED",
        },
        "license": "Open Scientific Literature / RDKit BSD",
        "model_priority": 300,
        "ensemble_eligible": False,
        "species": "Chemical",
        "output_type": "regression_logS",
        "provenance_json": {
            "model_id": "rdkit_gbr_solubility_v1",
            "model_role": "ADAPTIVE_EXCLUDED",
            "stage4d_classification": "EXCLUDED from production/consensus",
            "production_eligible": False,
            "reason": "Stage 4D-3A2: M3 excluded from any production or research consensus",
        },
        "is_active": True,
    },
    # Permeability M2 — Physicochemical Caco-2
    {
        "endpoint_name": "Permeability",
        "model_name": "Physicochemical Polar Surface Caco-2 Permeability",
        "model_version": "physchem-caco2-v1.0",
        "implementation_status": "READY",
        "supported_species": ["Human"],
        "supported_matrix": ["Caco-2"],
        "output_unit": "log10(cm/s)",
        "source": "Wang et al. Caco-2 physicochemical equation; Egan egg rule-derived",
        "training_dataset": "Wang et al. Caco-2 Curated Dataset (N=1,272)",
        "validation_json": {
            "note": "Physicochemical rule-based. SHADOW_ONLY. Insufficient N for ensemble (Stage 4D-2C).",
            "model_role": "SHADOW_ONLY",
        },
        "license": "Open Scientific Literature",
        "model_priority": 200,
        "ensemble_eligible": False,
        "species": "Human",
        "output_type": "regression_logPapp",
        "provenance_json": {
            "model_id": "physchem_caco2_v1",
            "model_role": "SHADOW_ONLY",
            "stage4d_classification": "SHADOW / research supporting",
            "consensus_permission": "INSUFFICIENT_EVIDENCE",
            "production_eligible": False,
            "reason": "Stage 4D-2C: KEEP_SHADOW; insufficient external N to promote ensemble",
        },
        "is_active": True,
    },
    # CYP3A4 inhibitor M2 — Morgan GBR
    {
        "endpoint_name": "CYP3A4 inhibitor",
        "model_name": "Morgan ECFP4 CYP3A4 Inhibitor Classifier",
        "model_version": "morgan-cyp3a4-v1.0",
        "implementation_status": "READY",
        "supported_species": ["Human"],
        "supported_matrix": ["Microsomal/Cellular"],
        "output_unit": "probability",
        "source": "PubChem AID 1851 CYP3A4 Inhibition; Morgan ECFP4 + GBR",
        "training_dataset": "PubChem AID 1851 CYP3A4 Inhibition (N=12,320)",
        "validation_json": {
            "note": "FIXED_BLEND_RESEARCH_SHADOW_ONLY. Dynamic weighting NO_ADAPTIVE_VALUE (Stage 4D-3B1A).",
            "model_role": "FIXED_BLEND_RESEARCH_SHADOW_ONLY",
            "fixed_blend_weight": 0.0422,
        },
        "license": "Public Domain",
        "model_priority": 200,
        "ensemble_eligible": False,
        "species": "Human",
        "output_type": "binary_classification_prob",
        "provenance_json": {
            "model_id": "morgan_cyp3a4_inh_v1",
            "model_role": "FIXED_BLEND_RESEARCH_SHADOW_ONLY",
            "stage4d_classification": "CALIBRATION_SUPPORTING / SHADOW_ONLY",
            "fixed_blend_weight_M2": 0.0422,
            "fixed_blend_weight_M1": 0.9578,
            "production_eligible": False,
            "reason": "Stage 4D-3B1A: adaptive weighting NO_GO; fixed blend research only, does not outperform M1",
        },
        "is_active": True,
    },
    # hERG M2 — Physicochemical Basic Center
    {
        "endpoint_name": "hERG liability",
        "model_name": "Physicochemical Basic Center hERG Blocker Classifier",
        "model_version": "physchem-herg-v1.0",
        "implementation_status": "READY",
        "supported_species": ["Human"],
        "supported_matrix": ["In Vitro Patch-Clamp / Radioligand Binding"],
        "output_unit": "probability",
        "source": "Pharmacophore lipophilicity + basic center logistic model",
        "training_dataset": "Synthetic pharmacophore features; no distinct training set",
        "validation_json": {
            "note": "CALIBRATION_SUPPORTING_SHADOW_ONLY. M2 rescue rate 5.4% (Stage 4D-3B2A).",
            "model_role": "CALIBRATION_SUPPORTING_SHADOW_ONLY_NOT_DISCRIMINATIVE_BLEND",
            "m2_rescue_rate": 0.054,
        },
        "license": "Open Scientific Literature",
        "model_priority": 200,
        "ensemble_eligible": False,
        "species": "Human",
        "output_type": "binary_classification_prob",
        "provenance_json": {
            "model_id": "physchem_herg_v1",
            "model_role": "CALIBRATION_SUPPORTING_SHADOW_ONLY_NOT_DISCRIMINATIVE_BLEND",
            "stage4d_classification": "CALIBRATION_SUPPORTING / SHADOW_ONLY",
            "production_eligible": False,
            "reason": (
                "Stage 4D-3B2A: M2 rescue rate 5.4%; BETTER_SECONDARY_MODEL_REQUIRED. "
                "Adaptive weighting NO_GO. Production remains raw M1 at threshold 0.50."
            ),
        },
        "is_active": True,
    },
]


def ensure_shadow_model_registry(connection: Any, ADMETModelRegistry_table: Any) -> None:
    """
    Idempotent migration to register shadow model rows in admet_model_registry.
    Called from ensure_admet_schema in admet.py.
    """
    from sqlalchemy import select as sa_select
    for seed in SHADOW_MODEL_SEEDS:
        existing = connection.execute(
            sa_select(ADMETModelRegistry_table.c.id).where(
                ADMETModelRegistry_table.c.endpoint_name == seed["endpoint_name"],
                ADMETModelRegistry_table.c.model_version == seed["model_version"],
                ADMETModelRegistry_table.c.model_name == seed["model_name"],
            ).limit(1)
        ).scalar()
        if existing is None:
            connection.execute(
                ADMETModelRegistry_table.insert().values(**seed)
            )
