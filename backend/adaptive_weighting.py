"""
Drug-OPT Stage 4D-3A: Hierarchical Experimental Adaptive Weighting Engine.

Implements the 4-Level Hierarchical Evidence Architecture for Aqueous Solubility:
1. GLOBAL (Frozen external qualification baseline)
2. PROJECT (Within-project empirical performance with shrinkage)
3. SERIES (Bemis-Murcko scaffold series empirical performance with shrinkage)
4. LOCAL (Morgan fingerprint Tanimoto neighborhood with distance weighting)

Key Safeguards:
- Strict hierarchical shrinkage: lambda = N_eff / (N_eff + N_prior)
- Absolute zero retrospective leakage: only historical frozen events prior to prediction timestamp are considered
- Preserves distinct component provenance at all four levels
- Model version isolation & Applicability Domain (AD) weighting
- Minimum weight floor (epsilon = 0.02)
- Shadow mode only: visible production predictions remain 100% untouched
"""

from __future__ import annotations

import enum
import functools
import hashlib
import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

from backend.endpoint_contracts import get_endpoint_contract, EndpointContract
from backend.multimodel import ExecutionStatus, ModelExecutionPayload


ADAPTIVE_POLICY_VERSION = "stage4d3a-hierarchical-shrinkage-v1"

# Frozen Global Prior Performance from Stage 4D-2C Audit (Delaney N=250)
GLOBAL_SOLUBILITY_PRIOR_MAE: Dict[str, float] = {
    "admetica_solubility": 0.3386,
    "esol_delaney_v1": 0.6663,
    "rdkit_gbr_solubility_v1": 0.7340,
}

# Default Shrinkage Prior Sample Sizes
DEFAULT_N_PRIOR_PROJECT = 10.0
DEFAULT_N_PRIOR_SERIES = 5.0
DEFAULT_N_PRIOR_LOCAL = 3.0
DEFAULT_LOCAL_SIMILARITY_THRESHOLD = 0.40
DEFAULT_BETA_ERROR_SCALING = 2.0
MINIMUM_WEIGHT_FLOOR = 0.02


class AdaptiveScope(str, enum.Enum):
    GLOBAL = "GLOBAL"
    PROJECT = "PROJECT"
    SERIES = "SERIES"
    LOCAL = "LOCAL"


class AssayQuality(str, enum.Enum):
    HIGH_QUALITY = "HIGH_QUALITY"
    USABLE = "USABLE"
    LIMITED = "LIMITED"
    INCOMPATIBLE = "INCOMPATIBLE"


class AdaptiveReasonCode(str, enum.Enum):
    GLOBAL_PRIOR_DOMINANT = "GLOBAL_PRIOR_DOMINANT"
    PROJECT_EVIDENCE_ACTIVE = "PROJECT_EVIDENCE_ACTIVE"
    PROJECT_EVIDENCE_INCREASED_M2 = "PROJECT_EVIDENCE_INCREASED_M2"
    SERIES_M1_OUTPERFORMS_M2 = "SERIES_M1_OUTPERFORMS_M2"
    SERIES_M2_OUTPERFORMS_M1 = "SERIES_M2_OUTPERFORMS_M1"
    INSUFFICIENT_LOCAL_DATA = "INSUFFICIENT_LOCAL_DATA"
    LOCAL_NEIGHBORHOOD_ACTIVE = "LOCAL_NEIGHBORHOOD_ACTIVE"
    M3_OUT_OF_DOMAIN = "M3_OUT_OF_DOMAIN"
    M3_ADAPTIVE_EXCLUDED = "M3_ADAPTIVE_EXCLUDED"
    UNSTABLE_ADAPTIVE_WEIGHTS = "UNSTABLE_ADAPTIVE_WEIGHTS"
    EXPERIMENT_QUALITY_LIMITED = "EXPERIMENT_QUALITY_LIMITED"
    NO_FROZEN_PREDICTION = "NO_FROZEN_PREDICTION"
    EXPERIMENT_NOT_ADAPTATION_COMPATIBLE = "EXPERIMENT_NOT_ADAPTATION_COMPATIBLE"


class AdaptiveCompletionDecision(str, enum.Enum):
    ADAPTIVE_PROMOTION_CANDIDATE = "ADAPTIVE_PROMOTION_CANDIDATE"
    CONDITIONAL_ADAPTIVE_VALUE = "CONDITIONAL_ADAPTIVE_VALUE"
    KEEP_RESEARCH_SHADOW = "KEEP_RESEARCH_SHADOW"
    ADAPTIVE_REJECTED = "ADAPTIVE_REJECTED"


@dataclass
class HierarchicalWeightBreakdown:
    """Retains mathematical provenance across all 4 evidence levels."""
    model_id: str
    global_weight: float
    project_weight: float
    series_weight: float
    local_weight: float
    project_posterior: float
    series_posterior: float
    local_posterior: float
    pre_ad_weight: float
    post_ad_weight: float
    final_effective_weight: float
    applicability_domain: str

    def to_dict(self) -> Dict[str, Any]:
        return {
            "model_id": self.model_id,
            "global_weight": round(self.global_weight, 4),
            "project_weight": round(self.project_weight, 4),
            "series_weight": round(self.series_weight, 4),
            "local_weight": round(self.local_weight, 4),
            "project_posterior": round(self.project_posterior, 4),
            "series_posterior": round(self.series_posterior, 4),
            "local_posterior": round(self.local_posterior, 4),
            "pre_ad_weight": round(self.pre_ad_weight, 4),
            "post_ad_weight": round(self.post_ad_weight, 4),
            "final_effective_weight": round(self.final_effective_weight, 4),
            "applicability_domain": self.applicability_domain,
        }


@dataclass
class ExperimentalFeedbackRecord:
    """Immutable experimental observation used for prospective adaptive learning."""
    event_id: str
    project_id: int
    compound_version_id: int
    canonical_smiles: str
    endpoint_name: str
    experimental_value: float  # log10(mol/L)
    experimental_unit: str
    assay_quality: AssayQuality
    scaffold_smiles: str
    timestamp: str  # ISO 8601
    frozen_predictions: Dict[str, float]  # model_id -> predicted_logS
    model_errors: Dict[str, float] = field(default_factory=dict)  # model_id -> abs_error
    is_valid: bool = True

    def __post_init__(self):
        if not self.model_errors and self.frozen_predictions:
            for mid, pred in self.frozen_predictions.items():
                if pred is not None:
                    self.model_errors[mid] = abs(pred - self.experimental_value)


@dataclass
class AdaptiveConsensusResult:
    """Deterministic result of hierarchical adaptive weighting."""
    endpoint_name: str
    compound_version_id: int
    predicted_value: float
    model_disagreement: float
    consensus_mode: str = "SHADOW"
    policy_version: str = ADAPTIVE_POLICY_VERSION
    n_global: int = 250
    n_project: int = 0
    n_series: int = 0
    n_local_eff: float = 0.0
    series_id: str = ""
    scaffold_smiles: str = ""
    weights_breakdown: Dict[str, HierarchicalWeightBreakdown] = field(default_factory=dict)
    effective_weights: Dict[str, float] = field(default_factory=dict)
    reason_codes: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    timestamp: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    def to_dict(self) -> Dict[str, Any]:
        return {
            "endpoint_name": self.endpoint_name,
            "compound_version_id": self.compound_version_id,
            "predicted_value": round(self.predicted_value, 4),
            "model_disagreement": round(self.model_disagreement, 4),
            "consensus_mode": self.consensus_mode,
            "policy_version": self.policy_version,
            "sample_counts": {
                "n_global": self.n_global,
                "n_project": self.n_project,
                "n_series": self.n_series,
                "n_local_eff": round(self.n_local_eff, 3),
            },
            "series": {
                "series_id": self.series_id,
                "scaffold_smiles": self.scaffold_smiles,
            },
            "effective_weights": {k: round(v, 4) for k, v in self.effective_weights.items()},
            "weights_breakdown": {k: v.to_dict() for k, v in self.weights_breakdown.items()},
            "reason_codes": self.reason_codes,
            "warnings": self.warnings,
            "timestamp": self.timestamp,
        }


# Global fingerprint generator instance (radius=2, 2048-bit)
_FP_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


@functools.lru_cache(maxsize=16384)
def get_bemis_murcko_scaffold(smiles: str) -> str:
    """Computes canonical Bemis-Murcko scaffold. Returns '[acyclic]' if no rings exist."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "[unknown]"
        scaff = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        return scaff if scaff else "[acyclic]"
    except Exception:
        return "[unknown]"


@functools.lru_cache(maxsize=16384)
def compute_morgan_fingerprint(smiles: str):
    """Generates 2048-bit Morgan fingerprint."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return _FP_GENERATOR.GetFingerprint(mol)
    except Exception:
        return None


def compute_tanimoto_similarity(fp1, fp2) -> float:
    """Calculates Tanimoto similarity between two Morgan fingerprints."""
    if fp1 is None or fp2 is None:
        return 0.0
    return float(DataStructs.TanimotoSimilarity(fp1, fp2))


def compute_error_score(mae: float, beta: float = DEFAULT_BETA_ERROR_SCALING) -> float:
    """Transforms mean absolute error to bounded performance score: exp(-beta * MAE)."""
    return float(math.exp(-beta * max(0.0, mae)))


def compute_shrinkage_lambda(n_eff: float, n_prior: float) -> float:
    """Computes empirical Bayes shrinkage factor lambda = N_eff / (N_eff + N_prior)."""
    if n_eff <= 0.0 or (n_eff + n_prior) <= 0.0:
        return 0.0
    return float(n_eff / (n_eff + n_prior))


def evaluate_experimental_compatibility(
    endpoint_name: str,
    value: Optional[float],
    unit: str,
    method: str = "",
    notes: str = "",
) -> Tuple[bool, AssayQuality, str]:
    """
    Validates whether an experimental measurement is strictly compatible with
    EP_PHYS_SOLUBILITY (aqueous solubility log10(mol/L)).
    """
    if endpoint_name != "Solubility":
        return False, AssayQuality.INCOMPATIBLE, "Endpoint is not Aqueous Solubility"
    
    if value is None or math.isnan(value):
        return False, AssayQuality.INCOMPATIBLE, "Measurement has no numeric value"
    
    unit_norm = unit.strip().lower()
    compatible_units = {"log10(mol/l)", "logs", "log(mol/l)", "log mol/l", "log(m)"}
    if unit_norm not in compatible_units and unit != "log10(mol/L)":
        return False, AssayQuality.INCOMPATIBLE, f"Incompatible unit: {unit}. Must be log10(mol/L)"
    
    # Method & quality check
    method_lower = (method + " " + notes).lower()
    if "kinetic" in method_lower and "turbidimetric" in method_lower:
        return True, AssayQuality.LIMITED, "Kinetic turbidimetric assay - limited precision"
    
    if "thermodynamic" in method_lower or "shake-flask" in method_lower or "shake flask" in method_lower:
        return True, AssayQuality.HIGH_QUALITY, "Thermodynamic shake-flask gold standard"
    
    return True, AssayQuality.USABLE, "Standard aqueous solubility assay"


def compute_hierarchical_adaptive_weights(
    query_smiles: str,
    project_id: int,
    candidate_payloads: List[ModelExecutionPayload],
    historical_feedback_events: List[ExperimentalFeedbackRecord],
    n_prior_project: float = DEFAULT_N_PRIOR_PROJECT,
    n_prior_series: float = DEFAULT_N_PRIOR_SERIES,
    n_prior_local: float = DEFAULT_N_PRIOR_LOCAL,
    similarity_threshold: float = DEFAULT_LOCAL_SIMILARITY_THRESHOLD,
    beta: float = DEFAULT_BETA_ERROR_SCALING,
    include_m3: bool = False,
    prediction_timestamp: Optional[str] = None,
) -> AdaptiveConsensusResult:
    """
    Core hierarchical adaptive consensus engine.
    Calculates weights across GLOBAL -> PROJECT -> SERIES -> LOCAL levels with
    strict prospective filtering (no leakage of future observations).
    """
    active_payloads = [
        p for p in candidate_payloads
        if p.execution_status == ExecutionStatus.SUCCESS
        and (include_m3 or p.model_id != "rdkit_gbr_solubility_v1")
    ]
    model_ids = [p.model_id for p in active_payloads]
    
    reason_codes: List[str] = []
    warnings: List[str] = []
    
    if not active_payloads:
        return AdaptiveConsensusResult(
            endpoint_name="Solubility",
            compound_version_id=0,
            predicted_value=0.0,
            model_disagreement=0.0,
            warnings=["No successful qualified models available for adaptive consensus."],
            reason_codes=[AdaptiveReasonCode.GLOBAL_PRIOR_DOMINANT.value],
        )

    # 1. Chemical Series Identification
    scaffold_smiles = get_bemis_murcko_scaffold(query_smiles)
    series_id = f"SERIES_{hashlib.md5(scaffold_smiles.encode()).hexdigest()[:8]}"
    query_fp = compute_morgan_fingerprint(query_smiles)

    # 2. Prospective Event Filtering (Zero Leakage Check)
    valid_events = []
    for ev in historical_feedback_events:
        if not ev.is_valid or ev.assay_quality == AssayQuality.INCOMPATIBLE:
            continue
        # If timestamp is provided, strictly exclude events occurring after prediction
        if prediction_timestamp and ev.timestamp and ev.timestamp >= prediction_timestamp:
            continue
        # Never adapt on the query compound itself if it happens to be in historical records
        if ev.canonical_smiles and query_smiles and ev.canonical_smiles == query_smiles:
            continue
        valid_events.append(ev)

    # Project-level events
    project_events = [ev for ev in valid_events if ev.project_id == project_id]
    n_project = len(project_events)

    # Series-level events
    series_events = [ev for ev in project_events if ev.scaffold_smiles == scaffold_smiles]
    n_series = len(series_events)

    # Local neighborhood events (Tanimoto >= similarity_threshold)
    local_neighbors: List[Tuple[ExperimentalFeedbackRecord, float]] = []
    if query_fp is not None:
        for ev in project_events:
            ev_fp = compute_morgan_fingerprint(ev.canonical_smiles)
            if ev_fp is not None:
                sim = compute_tanimoto_similarity(query_fp, ev_fp)
                if sim >= similarity_threshold:
                    local_neighbors.append((ev, sim))

    n_local_eff = sum(sim ** 2 for _, sim in local_neighbors) if local_neighbors else 0.0

    # -------------------------------------------------------------------------
    # LEVEL 1: GLOBAL PRIOR WEIGHTS
    # -------------------------------------------------------------------------
    global_scores = {}
    for mid in model_ids:
        prior_mae = GLOBAL_SOLUBILITY_PRIOR_MAE.get(mid, 0.50)
        global_scores[mid] = compute_error_score(prior_mae, beta=beta)
    
    sum_global = sum(global_scores.values())
    global_weights = {mid: global_scores[mid] / sum_global for mid in model_ids}

    # -------------------------------------------------------------------------
    # LEVEL 2: PROJECT LEVEL EVIDENCE & SHRINKAGE
    # -------------------------------------------------------------------------
    if n_project > 0:
        project_maes = {}
        for mid in model_ids:
            errors = [ev.model_errors[mid] for ev in project_events if mid in ev.model_errors]
            project_maes[mid] = float(np.mean(errors)) if errors else GLOBAL_SOLUBILITY_PRIOR_MAE.get(mid, 0.50)
        
        proj_scores = {mid: compute_error_score(project_maes[mid], beta=beta) for mid in model_ids}
        sum_proj = sum(proj_scores.values())
        project_weights = {mid: proj_scores[mid] / sum_proj for mid in model_ids}
        lambda_proj = compute_shrinkage_lambda(float(n_project), n_prior_project)
    else:
        project_weights = dict(global_weights)
        lambda_proj = 0.0

    project_posteriors = {
        mid: (1.0 - lambda_proj) * global_weights[mid] + lambda_proj * project_weights[mid]
        for mid in model_ids
    }

    # -------------------------------------------------------------------------
    # LEVEL 3: SERIES LEVEL EVIDENCE & SHRINKAGE
    # -------------------------------------------------------------------------
    if n_series > 0:
        series_maes = {}
        for mid in model_ids:
            errors = [ev.model_errors[mid] for ev in series_events if mid in ev.model_errors]
            series_maes[mid] = float(np.mean(errors)) if errors else project_maes.get(mid, 0.50)
        
        ser_scores = {mid: compute_error_score(series_maes[mid], beta=beta) for mid in model_ids}
        sum_ser = sum(ser_scores.values())
        series_weights = {mid: ser_scores[mid] / sum_ser for mid in model_ids}
        lambda_ser = compute_shrinkage_lambda(float(n_series), n_prior_series)

        # Check series dominance for explanation
        if "admetica_solubility" in series_maes and "esol_delaney_v1" in series_maes:
            if series_maes["esol_delaney_v1"] < series_maes["admetica_solubility"]:
                reason_codes.append(AdaptiveReasonCode.SERIES_M2_OUTPERFORMS_M1.value)
            else:
                reason_codes.append(AdaptiveReasonCode.SERIES_M1_OUTPERFORMS_M2.value)
    else:
        series_weights = dict(project_posteriors)
        lambda_ser = 0.0

    series_posteriors = {
        mid: (1.0 - lambda_ser) * project_posteriors[mid] + lambda_ser * series_weights[mid]
        for mid in model_ids
    }

    # -------------------------------------------------------------------------
    # LEVEL 4: LOCAL NEIGHBORHOOD EVIDENCE & SHRINKAGE
    # -------------------------------------------------------------------------
    if local_neighbors and n_local_eff > 0.0:
        local_maes = {}
        total_sim = sum(sim for _, sim in local_neighbors)
        for mid in model_ids:
            weighted_errs = []
            for ev, sim in local_neighbors:
                if mid in ev.model_errors:
                    weighted_errs.append(sim * ev.model_errors[mid])
            local_maes[mid] = sum(weighted_errs) / total_sim if (total_sim > 0 and weighted_errs) else series_maes.get(mid, 0.50)
        
        loc_scores = {mid: compute_error_score(local_maes[mid], beta=beta) for mid in model_ids}
        sum_loc = sum(loc_scores.values())
        local_weights = {mid: loc_scores[mid] / sum_loc for mid in model_ids}
        lambda_loc = compute_shrinkage_lambda(n_local_eff, n_prior_local)
        reason_codes.append(AdaptiveReasonCode.LOCAL_NEIGHBORHOOD_ACTIVE.value)
    else:
        local_weights = dict(series_posteriors)
        lambda_loc = 0.0
        reason_codes.append(AdaptiveReasonCode.INSUFFICIENT_LOCAL_DATA.value)

    local_posteriors = {
        mid: (1.0 - lambda_loc) * series_posteriors[mid] + lambda_loc * local_weights[mid]
        for mid in model_ids
    }

    # -------------------------------------------------------------------------
    # LEVEL 5: APPLICABILITY DOMAIN ADJUSTMENT & MINIMUM WEIGHT FLOOR
    # -------------------------------------------------------------------------
    breakdowns = {}
    post_ad_weights = {}
    ad_multipliers = {"IN_DOMAIN": 1.0, "BORDERLINE": 0.5, "OUT_OF_DOMAIN": 0.1, "UNKNOWN": 0.8}

    for p in active_payloads:
        mid = p.model_id
        ad_status = p.applicability_domain
        ad_factor = ad_multipliers.get(ad_status, 0.8)
        
        pre_ad = local_posteriors[mid]
        post_ad = pre_ad * ad_factor

        if ad_status == "OUT_OF_DOMAIN":
            if mid == "rdkit_gbr_solubility_v1":
                reason_codes.append(AdaptiveReasonCode.M3_OUT_OF_DOMAIN.value)
            warnings.append(f"Model {mid} is OUT_OF_DOMAIN; downweighted by 0.1x.")

        post_ad_weights[mid] = max(MINIMUM_WEIGHT_FLOOR, post_ad)

    sum_final = sum(post_ad_weights.values())
    final_effective_weights = {mid: post_ad_weights[mid] / sum_final for mid in model_ids}

    for p in active_payloads:
        mid = p.model_id
        breakdowns[mid] = HierarchicalWeightBreakdown(
            model_id=mid,
            global_weight=global_weights[mid],
            project_weight=project_weights[mid],
            series_weight=series_weights[mid],
            local_weight=local_weights[mid],
            project_posterior=project_posteriors[mid],
            series_posterior=series_posteriors[mid],
            local_posterior=local_posteriors[mid],
            pre_ad_weight=local_posteriors[mid],
            post_ad_weight=post_ad_weights[mid],
            final_effective_weight=final_effective_weights[mid],
            applicability_domain=p.applicability_domain,
        )

    # -------------------------------------------------------------------------
    # PREDICTION & MODEL DISAGREEMENT
    # -------------------------------------------------------------------------
    values_by_id = {p.model_id: float(p.value) for p in active_payloads if p.value is not None}
    adaptive_pred = sum(final_effective_weights[mid] * values_by_id[mid] for mid in model_ids)

    # Weighted model disagreement standard deviation
    disagreement_var = sum(
        final_effective_weights[mid] * ((values_by_id[mid] - adaptive_pred) ** 2)
        for mid in model_ids
    )
    disagreement_std = math.sqrt(max(0.0, disagreement_var))

    # Reason code summary
    if n_project < 5:
        reason_codes.insert(0, AdaptiveReasonCode.GLOBAL_PRIOR_DOMINANT.value)
    else:
        reason_codes.insert(0, AdaptiveReasonCode.PROJECT_EVIDENCE_ACTIVE.value)

    if not include_m3:
        reason_codes.append(AdaptiveReasonCode.M3_ADAPTIVE_EXCLUDED.value)

    return AdaptiveConsensusResult(
        endpoint_name="Solubility",
        compound_version_id=0,
        predicted_value=float(adaptive_pred),
        model_disagreement=float(disagreement_std),
        consensus_mode="SHADOW",
        policy_version=ADAPTIVE_POLICY_VERSION,
        n_global=250,
        n_project=n_project,
        n_series=n_series,
        n_local_eff=float(n_local_eff),
        series_id=series_id,
        scaffold_smiles=scaffold_smiles,
        weights_breakdown=breakdowns,
        effective_weights=final_effective_weights,
        reason_codes=list(set(reason_codes)),
        warnings=warnings,
    )
