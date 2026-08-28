"""
Drug-OPT Stage 4D-1: Static Consensus & Aggregation Engine.

Provides:
- Strict separation of aggregation types:
  - REGRESSION_WEIGHTED (weighted mean + model disagreement standard deviation)
  - CLASSIFICATION_WEIGHTED (weighted probability + voting pattern)
  - RANK_FUSION (Reciprocal Rank Fusion for Site-of-Metabolism)
  - NO_CONSENSUS (PK mechanistic methods)
- Conservative static weighting combining base quality, applicability domain,
  confidence, and dataset/architecture diversity penalties
- Fault-tolerant failure renormalization
- Agreement classification (HIGH_AGREEMENT, MODERATE_AGREEMENT, LOW_AGREEMENT)
- SHADOW mode enforcement: consensus is calculated and logged for audit without
  overwriting visible production primary predictions.
"""

from __future__ import annotations

import enum
import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np

from backend.endpoint_contracts import (
    ENDPOINT_CONTRACTS,
    EndpointContract,
    OutputType,
    Directionality,
    get_endpoint_contract,
    check_ensemble_compatibility,
)
from backend.multimodel import ExecutionStatus, ModelExecutionPayload


class ConsensusMode(str, enum.Enum):
    OFF = "OFF"
    SHADOW = "SHADOW"
    ACTIVE = "ACTIVE"


class AggregationType(str, enum.Enum):
    REGRESSION_WEIGHTED = "REGRESSION_WEIGHTED"
    CLASSIFICATION_WEIGHTED = "CLASSIFICATION_WEIGHTED"
    RANK_FUSION = "RANK_FUSION"
    NO_CONSENSUS = "NO_CONSENSUS"


class AgreementStatus(str, enum.Enum):
    HIGH_AGREEMENT = "HIGH_AGREEMENT"
    MODERATE_AGREEMENT = "MODERATE_AGREEMENT"
    LOW_AGREEMENT = "LOW_AGREEMENT"
    SINGLE_MODEL = "SINGLE_MODEL"
    NO_CONSENSUS = "NO_CONSENSUS"


# Endpoint aggregation type mapping
ENDPOINT_AGGREGATION_MAP: Dict[str, AggregationType] = {
    # Physicochemical & Absorption & Distribution
    "Solubility": AggregationType.REGRESSION_WEIGHTED,
    "Permeability": AggregationType.REGRESSION_WEIGHTED,
    "Plasma protein binding": AggregationType.REGRESSION_WEIGHTED,
    # Microsomal clearance
    "HLM intrinsic clearance": AggregationType.REGRESSION_WEIGHTED,
    "RLM intrinsic clearance": AggregationType.REGRESSION_WEIGHTED,
    "MLM intrinsic clearance": AggregationType.REGRESSION_WEIGHTED,
    # CYP panel
    "CYP1A2 inhibitor": AggregationType.CLASSIFICATION_WEIGHTED,
    "CYP2C9 inhibitor": AggregationType.CLASSIFICATION_WEIGHTED,
    "CYP2C19 inhibitor": AggregationType.CLASSIFICATION_WEIGHTED,
    "CYP2D6 inhibitor": AggregationType.CLASSIFICATION_WEIGHTED,
    "CYP3A4 inhibitor": AggregationType.CLASSIFICATION_WEIGHTED,
    "CYP2C9 substrate": AggregationType.CLASSIFICATION_WEIGHTED,
    "CYP2D6 substrate": AggregationType.CLASSIFICATION_WEIGHTED,
    "CYP3A4 substrate": AggregationType.CLASSIFICATION_WEIGHTED,
    # Transporters & Safety
    "P-gp inhibitor": AggregationType.CLASSIFICATION_WEIGHTED,
    "hERG liability": AggregationType.CLASSIFICATION_WEIGHTED,
    "Ames mutagenicity": AggregationType.CLASSIFICATION_WEIGHTED,
    "DILI clinical liability": AggregationType.CLASSIFICATION_WEIGHTED,
    # Metabolism
    "Metabolic soft spots": AggregationType.RANK_FUSION,
    # Pharmacokinetics (Mechanistic methods: no consensus)
    "PK Systemic Clearance": AggregationType.NO_CONSENSUS,
    "PK Volume of Distribution": AggregationType.NO_CONSENSUS,
    "PK Bioavailability": AggregationType.NO_CONSENSUS,
}

# Empirical architecture/dataset error correlation diversity penalties (Stage 4D-2 Evidence)
# Penalty Factor = max(0.10, 1.0 - r_error^2)
DIVERSITY_PENALTY_PAIRS: Dict[Tuple[str, str], float] = {
    ("admetica", "admet_ai"): 0.55,
    ("admet_ai", "admetica"): 0.55,
    ("admetica_solubility", "esol_delaney_v1"): 0.85,
    ("esol_delaney_v1", "admetica_solubility"): 0.85,
    ("admetica_solubility", "rdkit_gbr_solubility_v1"): 0.80,
    ("rdkit_gbr_solubility_v1", "admetica_solubility"): 0.80,
    ("esol_delaney_v1", "rdkit_gbr_solubility_v1"): 0.25,
    ("rdkit_gbr_solubility_v1", "esol_delaney_v1"): 0.25,
    ("admetica_caco2", "physchem_caco2_v1"): 0.73,
    ("physchem_caco2_v1", "admetica_caco2"): 0.73,
    ("admetica_cyp_cyp3a4-inhibitor", "morgan_cyp3a4_inh_v1"): 0.95,
    ("morgan_cyp3a4_inh_v1", "admetica_cyp_cyp3a4-inhibitor"): 0.95,
    ("admetica_safety_herg", "physchem_herg_v1"): 0.20,
    ("physchem_herg_v1", "admetica_safety_herg"): 0.20,
}


@dataclass
class ConsensusResult:
    """Standardized record of a multi-model consensus prediction."""
    endpoint_id: str
    endpoint_name: str
    compound_version_id: int
    consensus_version: str = "stage4d1-static-v1"
    consensus_mode: ConsensusMode = ConsensusMode.SHADOW
    aggregation_type: AggregationType = AggregationType.REGRESSION_WEIGHTED
    combined_value: Optional[float] = None
    combined_probability: Optional[float] = None
    consensus_classification: Optional[str] = None
    canonical_unit: str = ""
    models_used: List[str] = field(default_factory=list)
    original_weights: Dict[str, float] = field(default_factory=dict)
    effective_weights: Dict[str, float] = field(default_factory=dict)
    model_agreement: AgreementStatus = AgreementStatus.SINGLE_MODEL
    dispersion: Dict[str, Any] = field(default_factory=dict)
    applicability_summary: str = "IN_DOMAIN"
    vote_pattern: Optional[str] = None
    warnings: List[str] = field(default_factory=list)
    provenance: Dict[str, Any] = field(default_factory=dict)
    generated_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    def to_dict(self) -> Dict[str, Any]:
        return {
            "endpoint_id": self.endpoint_id,
            "endpoint_name": self.endpoint_name,
            "compound_version_id": self.compound_version_id,
            "consensus_version": self.consensus_version,
            "consensus_mode": self.consensus_mode.value,
            "aggregation_type": self.aggregation_type.value,
            "combined_value": self.combined_value,
            "combined_probability": self.combined_probability,
            "consensus_classification": self.consensus_classification,
            "canonical_unit": self.canonical_unit,
            "models_used": self.models_used,
            "original_weights": self.original_weights,
            "effective_weights": self.effective_weights,
            "model_agreement": self.model_agreement.value,
            "dispersion": self.dispersion,
            "applicability_summary": self.applicability_summary,
            "vote_pattern": self.vote_pattern,
            "warnings": self.warnings,
            "provenance": self.provenance,
            "generated_at": self.generated_at,
        }


def calculate_static_model_weight(
    payload: ModelExecutionPayload,
    other_payloads: List[ModelExecutionPayload],
) -> Tuple[float, str]:
    """
    Computes static initial weight for an individual model execution payload.
    Factors: Base Quality x Applicability Domain x Confidence x Diversity Factor.
    """
    if payload.execution_status != ExecutionStatus.SUCCESS:
        return 0.0, f"Excluded: execution status is {payload.execution_status.value}"

    # 1. Applicability domain factor
    ad_factor = {
        "IN_DOMAIN": 1.0,
        "BORDERLINE": 0.70,
        "OUT_OF_DOMAIN": 0.10,
    }.get(payload.applicability_domain, 0.50)

    # 2. Confidence factor
    conf_factor = {
        "HIGH": 1.0,
        "MEDIUM": 0.85,
        "LOW": 0.65,
        "NOT_APPLICABLE": 0.70,
    }.get(payload.confidence, 0.70)

    # 3. Diversity factor (downweights models sharing dataset and architecture)
    diversity_factor = 1.0
    for other in other_payloads:
        if other.model_id != payload.model_id and other.execution_status == ExecutionStatus.SUCCESS:
            id_pair = (payload.model_id, other.model_id)
            fam_pair = (payload.model_family, other.model_family)
            if id_pair in DIVERSITY_PENALTY_PAIRS:
                diversity_factor = min(diversity_factor, DIVERSITY_PENALTY_PAIRS[id_pair])
            elif fam_pair in DIVERSITY_PENALTY_PAIRS:
                diversity_factor = min(diversity_factor, DIVERSITY_PENALTY_PAIRS[fam_pair])

    base_weight = 1.0
    weight = base_weight * ad_factor * conf_factor * diversity_factor
    reason = (
        f"Base(1.0) x AD({ad_factor:.2f}) x Conf({conf_factor:.2f}) x Diversity({diversity_factor:.2f})"
    )
    return max(0.001, weight), reason


def compute_endpoint_consensus(
    endpoint_name: str,
    compound_version_id: int,
    model_payloads: List[ModelExecutionPayload],
    mode: ConsensusMode = ConsensusMode.SHADOW,
) -> ConsensusResult:
    """
    Assembles static initial consensus across multiple model execution payloads
    for a specific endpoint according to its authoritative contract.
    """
    contract = get_endpoint_contract(endpoint_name)
    agg_type = ENDPOINT_AGGREGATION_MAP.get(endpoint_name, AggregationType.REGRESSION_WEIGHTED)

    # Base result if no models or PK exclusion
    if agg_type == AggregationType.NO_CONSENSUS or not model_payloads:
        return ConsensusResult(
            endpoint_id=contract.endpoint_id if contract else endpoint_name,
            endpoint_name=endpoint_name,
            compound_version_id=compound_version_id,
            consensus_mode=mode,
            aggregation_type=agg_type,
            canonical_unit=contract.canonical_unit if contract else "",
            model_agreement=AgreementStatus.NO_CONSENSUS,
            warnings=["Endpoint is designated NO_CONSENSUS (mechanistic/NCA method) or no models provided."],
        )

    # Filter successful payloads
    successful_payloads = [p for p in model_payloads if p.execution_status == ExecutionStatus.SUCCESS]
    failed_payloads = [p for p in model_payloads if p.execution_status != ExecutionStatus.SUCCESS]

    if not successful_payloads:
        return ConsensusResult(
            endpoint_id=contract.endpoint_id if contract else endpoint_name,
            endpoint_name=endpoint_name,
            compound_version_id=compound_version_id,
            consensus_mode=mode,
            aggregation_type=agg_type,
            canonical_unit=contract.canonical_unit if contract else "",
            model_agreement=AgreementStatus.NO_CONSENSUS,
            warnings=[f"All attempted models failed: {[f.model_id for f in failed_payloads]}"],
        )

    # Compute raw weights and reasons
    raw_weights = {}
    weight_reasons = {}
    for p in model_payloads:
        w, r = calculate_static_model_weight(p, model_payloads)
        raw_weights[p.model_id] = round(w, 4)
        weight_reasons[p.model_id] = r

    # Renormalize across successful models
    sum_success = sum(raw_weights[p.model_id] for p in successful_payloads)
    effective_weights = {}
    for p in successful_payloads:
        effective_weights[p.model_id] = round(raw_weights[p.model_id] / sum_success, 4) if sum_success > 0 else 1.0 / len(successful_payloads)

    # Applicability summary
    in_domain_weight = sum(effective_weights[p.model_id] for p in successful_payloads if p.applicability_domain == "IN_DOMAIN")
    out_domain_weight = sum(effective_weights[p.model_id] for p in successful_payloads if p.applicability_domain == "OUT_OF_DOMAIN")
    if in_domain_weight >= 0.60:
        overall_domain = "IN_DOMAIN"
    elif out_domain_weight >= 0.50:
        overall_domain = "OUT_OF_DOMAIN"
    else:
        overall_domain = "BORDERLINE"

    warnings = []
    if failed_payloads:
        warnings.append(f"Failure isolation: {len(failed_payloads)} model(s) failed and weights were renormalized ({[f.model_id for f in failed_payloads]}).")
    if overall_domain == "OUT_OF_DOMAIN":
        warnings.append("Consensus molecule is OUT_OF_DOMAIN; prediction consensus has elevated uncertainty.")

    # --------------------------------------------------------------------------
    # 1. REGRESSION AGGREGATION
    # --------------------------------------------------------------------------
    if agg_type == AggregationType.REGRESSION_WEIGHTED:
        valid_values = [(p, p.value, effective_weights[p.model_id]) for p in successful_payloads if p.value is not None]
        if not valid_values:
            return ConsensusResult(
                endpoint_id=contract.endpoint_id,
                endpoint_name=endpoint_name,
                compound_version_id=compound_version_id,
                consensus_mode=mode,
                aggregation_type=agg_type,
                canonical_unit=contract.canonical_unit,
                model_agreement=AgreementStatus.NO_CONSENSUS,
                warnings=["No finite numerical values emitted from successful models."],
            )

        # Weighted Mean
        combined_val = float(sum(val * w for _, val, w in valid_values))
        values_arr = np.array([val for _, val, _ in valid_values])
        weights_arr = np.array([w for _, _, w in valid_values])

        # Weighted Standard Deviation (Model Disagreement)
        if len(valid_values) > 1:
            variance = float(np.sum(weights_arr * (values_arr - combined_val) ** 2))
            model_disagreement_std = float(math.sqrt(max(0.0, variance)))
            val_min = float(np.min(values_arr))
            val_max = float(np.max(values_arr))
            val_range = round(val_max - val_min, 4)

            # Agreement Classification
            if model_disagreement_std <= 0.30:  # <= ~2-fold spread
                agreement = AgreementStatus.HIGH_AGREEMENT
            elif model_disagreement_std <= 0.60:  # ~2 to 4-fold spread
                agreement = AgreementStatus.MODERATE_AGREEMENT
            else:
                agreement = AgreementStatus.LOW_AGREEMENT
        else:
            model_disagreement_std = 0.0
            val_min = combined_val
            val_max = combined_val
            val_range = 0.0
            agreement = AgreementStatus.SINGLE_MODEL

        dispersion = {
            "model_disagreement_std": round(model_disagreement_std, 4),
            "min_value": round(val_min, 4),
            "max_value": round(val_max, 4),
            "range": val_range,
            "interpretation": "MODEL DISAGREEMENT (weighted standard deviation across model predictions; not a confidence interval)",
        }

        return ConsensusResult(
            endpoint_id=contract.endpoint_id,
            endpoint_name=endpoint_name,
            compound_version_id=compound_version_id,
            consensus_mode=mode,
            aggregation_type=agg_type,
            combined_value=round(combined_val, 4),
            canonical_unit=contract.canonical_unit,
            models_used=[p.model_id for p in successful_payloads],
            original_weights=raw_weights,
            effective_weights=effective_weights,
            model_agreement=agreement,
            dispersion=dispersion,
            applicability_summary=overall_domain,
            warnings=warnings,
            provenance={
                "weight_policy": "BaseQuality x ApplicabilityDomain x Confidence x DiversityPenalty",
                "weight_reasons": weight_reasons,
                "input_canonical_smiles": successful_payloads[0].canonical_smiles,
            },
        )

    # --------------------------------------------------------------------------
    # 2. CLASSIFICATION AGGREGATION
    # --------------------------------------------------------------------------
    elif agg_type == AggregationType.CLASSIFICATION_WEIGHTED:
        valid_probs = [(p, p.probability, effective_weights[p.model_id]) for p in successful_payloads if p.probability is not None]
        if not valid_probs:
            return ConsensusResult(
                endpoint_id=contract.endpoint_id,
                endpoint_name=endpoint_name,
                compound_version_id=compound_version_id,
                consensus_mode=mode,
                aggregation_type=agg_type,
                canonical_unit=contract.canonical_unit,
                model_agreement=AgreementStatus.NO_CONSENSUS,
                warnings=["No valid probabilities emitted."],
            )

        combined_prob = float(sum(prob * w for _, prob, w in valid_probs))
        cutoff = 0.50
        if contract.classification_semantics:
            cutoff = float(contract.classification_semantics.get("decision_threshold", 0.50))
            pos_label = contract.classification_semantics.get("positive_class", "POSITIVE")
            neg_label = contract.classification_semantics.get("negative_class", "NEGATIVE")
        else:
            pos_label = "POSITIVE"
            neg_label = "NEGATIVE"

        consensus_class = pos_label if combined_prob >= cutoff else neg_label

        # Model vote pattern
        votes = {}
        for p, prob, _ in valid_probs:
            votes[p.model_id] = pos_label if prob >= cutoff else neg_label
        vote_pattern = ", ".join(f"{mid}:{cls}" for mid, cls in votes.items())

        # Agreement classification
        pos_votes = sum(1 for cls in votes.values() if cls == pos_label)
        total_votes = len(votes)
        if total_votes == 1:
            agreement = AgreementStatus.SINGLE_MODEL
        elif pos_votes == 0 or pos_votes == total_votes:
            agreement = AgreementStatus.HIGH_AGREEMENT
        elif pos_votes / total_votes >= 0.70 or pos_votes / total_votes <= 0.30:
            agreement = AgreementStatus.MODERATE_AGREEMENT
        else:
            agreement = AgreementStatus.LOW_AGREEMENT

        dispersion = {
            "individual_probabilities": {p.model_id: round(p.probability, 4) for p in successful_payloads if p.probability is not None},
            "vote_pattern": votes,
            "decision_threshold": cutoff,
        }

        return ConsensusResult(
            endpoint_id=contract.endpoint_id,
            endpoint_name=endpoint_name,
            compound_version_id=compound_version_id,
            consensus_mode=mode,
            aggregation_type=agg_type,
            combined_value=round(combined_prob, 4),
            combined_probability=round(combined_prob, 4),
            consensus_classification=consensus_class,
            canonical_unit="probability",
            models_used=[p.model_id for p in successful_payloads],
            original_weights=raw_weights,
            effective_weights=effective_weights,
            model_agreement=agreement,
            dispersion=dispersion,
            applicability_summary=overall_domain,
            vote_pattern=vote_pattern,
            warnings=warnings,
            provenance={
                "weight_policy": "BaseQuality x ApplicabilityDomain x Confidence x DiversityPenalty",
                "weight_reasons": weight_reasons,
                "input_canonical_smiles": successful_payloads[0].canonical_smiles,
            },
        )

    # --------------------------------------------------------------------------
    # 3. RANK FUSION AGGREGATION (Site-of-Metabolism)
    # --------------------------------------------------------------------------
    elif agg_type == AggregationType.RANK_FUSION:
        # SyGMa soft spots rank fusion
        sygma_payload = next((p for p in successful_payloads if p.model_family == "rule_based_smarts"), None)
        fused_spots = []
        if sygma_payload and sygma_payload.raw_outputs:
            spots = sygma_payload.raw_outputs.get("spots", [])
            for rank_idx, spot in enumerate(spots, start=1):
                fused_spots.append({
                    "atom_index": spot.get("atom_index"),
                    "atom_environment": spot.get("atom_environment"),
                    "rank": rank_idx,
                    "rrf_score": round(1.0 / (60.0 + rank_idx), 4),
                    "reactions": spot.get("reactions", []),
                })

        return ConsensusResult(
            endpoint_id=contract.endpoint_id,
            endpoint_name=endpoint_name,
            compound_version_id=compound_version_id,
            consensus_mode=mode,
            aggregation_type=agg_type,
            canonical_unit=contract.canonical_unit,
            models_used=[p.model_id for p in successful_payloads],
            original_weights=raw_weights,
            effective_weights=effective_weights,
            model_agreement=AgreementStatus.SINGLE_MODEL,
            dispersion={"fused_soft_spots_count": len(fused_spots), "top_soft_spots": fused_spots[:5]},
            applicability_summary="IN_DOMAIN",
            warnings=warnings,
            provenance={"aggregation_method": "Reciprocal Rank Fusion (RRF)", "formula": "RRF = sum(1 / (60 + rank))"},
        )

    return ConsensusResult(
        endpoint_id=contract.endpoint_id if contract else endpoint_name,
        endpoint_name=endpoint_name,
        compound_version_id=compound_version_id,
        consensus_mode=mode,
        aggregation_type=agg_type,
        canonical_unit=contract.canonical_unit if contract else "",
        model_agreement=AgreementStatus.NO_CONSENSUS,
    )
