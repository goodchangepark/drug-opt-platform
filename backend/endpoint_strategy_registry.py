"""
Stage 4D-4: Endpoint Strategy Registry
=======================================

Versioned, authoritative scientific prediction governance layer.

This module encodes all endpoint-specific prediction policies derived from
Stages 4D-0 through 4D-3B2A. It is the single source of truth for:

  - Primary production strategy per endpoint
  - Shadow research strategy
  - Calibration status
  - Adaptive weighting permission
  - Consensus permission
  - Model uncertainty / disagreement policy
  - Promotion gate requirements
  - Rollback information

IMMUTABLE PRODUCTION CONTRACT:
  No endpoint policy changes visible production outputs without a separate
  Stage 4D-5 promotion validation cycle. This file encodes policy intent,
  not automatic activation.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional

from backend.endpoint_contracts import ENDPOINT_CONTRACTS


REGISTRY_VERSION = "stage4d4-endpoint-strategy-v1"
REGISTRY_POLICY_DATE = "2026-08-29"


# ──────────────────────────────────────────────────────────────────────────────
# Strategy Enums
# ──────────────────────────────────────────────────────────────────────────────

class StrategyType(str, enum.Enum):
    """Scientific prediction strategy, independent of its lifecycle state."""
    SINGLE_CORE_MODEL         = "SINGLE_CORE_MODEL"
    SINGLE_CORE_WITH_CALIBRATION = "SINGLE_CORE_WITH_CALIBRATION"
    FIXED_WEIGHT_BLEND        = "FIXED_WEIGHT_BLEND"
    STATIC_CONSENSUS          = "STATIC_CONSENSUS"
    ADAPTIVE_RESEARCH_SHADOW  = "ADAPTIVE_RESEARCH_SHADOW"
    RANK_FUSION               = "RANK_FUSION"
    RULE_BASED                = "RULE_BASED"
    MECHANISTIC_NO_CONSENSUS  = "MECHANISTIC_NO_CONSENSUS"
    MODEL_UNAVAILABLE         = "MODEL_UNAVAILABLE"
    DERIVED_ESTIMATE          = "DERIVED_ESTIMATE"
    RULE_ESTIMATE             = "RULE_ESTIMATE"


# Compatibility alias for callers written while the AGY draft used this name.
PrimaryStrategy = StrategyType


class AdaptiveStatus(str, enum.Enum):
    """Adaptive weighting permission status."""
    DISABLED                  = "DISABLED"
    NO_ADAPTIVE_VALUE         = "NO_ADAPTIVE_VALUE"
    ARCHITECTURE_VALID_BUT_NO_ACCURACY_GAIN = "ARCHITECTURE_VALID_BUT_NO_ACCURACY_GAIN"
    RESEARCH_CANDIDATE        = "RESEARCH_CANDIDATE"
    NO_GO                     = "NO_GO"


class CalibrationStatus(str, enum.Enum):
    """Probability calibration lifecycle state."""
    RAW                       = "RAW"
    CALIBRATION_RESEARCH      = "CALIBRATION_RESEARCH"
    CALIBRATION_VALIDATED     = "CALIBRATION_VALIDATED"
    CALIBRATION_PRODUCTION    = "CALIBRATION_PRODUCTION"
    NOT_APPLICABLE            = "NOT_APPLICABLE"


class ConsensusPermission(str, enum.Enum):
    """Whether consensus between models is allowed for this endpoint."""
    ALLOWED_PRODUCTION        = "ALLOWED_PRODUCTION"
    ALLOWED_SHADOW            = "ALLOWED_SHADOW"
    DISABLED                  = "DISABLED"
    INSUFFICIENT_EVIDENCE     = "INSUFFICIENT_EVIDENCE"
    MECHANISTICALLY_FORBIDDEN = "MECHANISTICALLY_FORBIDDEN"
    NOT_APPLICABLE            = "NOT_APPLICABLE"


class ValidationStatus(str, enum.Enum):
    """External validation status."""
    PUBLISHER_REPORTED_ONLY   = "PUBLISHER_REPORTED_ONLY"
    INDEPENDENT_PARTIAL       = "INDEPENDENT_PARTIAL"
    INDEPENDENT_COMPLETE      = "INDEPENDENT_COMPLETE"
    INDEPENDENT_COMPLETE_WITH_LIMITATIONS = "INDEPENDENT_COMPLETE_WITH_LIMITATIONS"
    NOT_VALIDATED             = "NOT_VALIDATED"
    MECHANISTIC               = "MECHANISTIC"


class PromotionStatus(str, enum.Enum):
    """Promotion lifecycle position."""
    SHADOW                    = "SHADOW"
    VALIDATED                 = "VALIDATED"
    PRODUCTION_CANDIDATE      = "PRODUCTION_CANDIDATE"
    ACTIVE                    = "ACTIVE"
    FROZEN                    = "FROZEN"
    DEFERRED                  = "DEFERRED"


class DisagreementPolicy(str, enum.Enum):
    """How to handle / report model disagreement."""
    NO_DISAGREEMENT_SINGLE_MODEL = "NO_DISAGREEMENT_SINGLE_MODEL"
    MODEL_DISAGREEMENT_SIGNAL    = "MODEL_DISAGREEMENT_SIGNAL"
    NOT_APPLICABLE               = "NOT_APPLICABLE"


COMMON_PROMOTION_REQUIREMENTS = [
    "Endpoint contract compatibility confirmed",
    "Held-out validation completed without test-set leakage",
    "Meaningful improvement demonstrated across more than one decision metric",
    "Calibration stable on an untouched holdout when calibration applies",
    "Subgroup and scaffold-series robustness demonstrated",
    "Immutable model and version identity recorded",
    "Deterministic rollback metadata and validation artifact recorded",
]


EVIDENCE_ARTIFACTS = {
    "stage3b": "docs/stage3b-model-selection.md",
    "stage3c": "docs/stage3c-model-selection.md",
    "stage3e": "docs/stage3e-model-selection.md",
    "stage4d0": "validation/stage4d0_candidate_model_qualification.json",
    "stage4d2c": "validation/stage4d2c_promotion_decisions.json",
    "stage4d3a2": "validation/stage4d3a2_final_decision.json",
    "stage4d3b1a": "validation/stage4d3b1a_final_decision.json",
    "stage4d3b2a": "validation/stage4d3b2a_final_decision.json",
    "stage3f": "docs/stage3f-model-selection.md",
    "stage4c4": "docs/stage4c4-pka-model-selection.md",
    "stage5a1": "docs/stage5a2a-ivive-hepatic-clearance.md",
}


@dataclass
class RollbackPolicy:
    """Complete deterministic rollback metadata for an active policy."""

    previous_policy_version: str
    rollback_target: str
    rollback_primary_strategy: StrategyType
    rollback_model_ids: List[str]
    rollback_model_versions: List[str]
    promotion_reason: str
    validation_artifact: str
    model_version_provenance: Dict[str, str]


# ──────────────────────────────────────────────────────────────────────────────
# Policy dataclass
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class EndpointStrategyPolicy:
    """Versioned scientific prediction governance record for one endpoint."""

    # Identification
    endpoint_name: str
    endpoint_id: str
    endpoint_contract_version: str = "stage4d4-endpoint-strategy-v1"

    # Primary production strategy
    primary_strategy: StrategyType = StrategyType.SINGLE_CORE_MODEL
    primary_model_ids: List[str] = field(default_factory=list)
    primary_model_versions: List[str] = field(default_factory=list)

    # Calibration
    calibration_status: CalibrationStatus = CalibrationStatus.RAW
    calibration_policy: str = "Raw model output; no post-hoc calibration applied in production."
    decision_threshold: Optional[float] = None

    # Shadow research
    shadow_strategy: Optional[StrategyType] = None
    shadow_model_ids: List[str] = field(default_factory=list)
    shadow_model_versions: List[str] = field(default_factory=list)
    non_primary_model_roles: Dict[str, str] = field(default_factory=dict)

    # Governance flags
    adaptive_status: AdaptiveStatus = AdaptiveStatus.DISABLED
    consensus_permission: ConsensusPermission = ConsensusPermission.DISABLED
    disagreement_policy: DisagreementPolicy = DisagreementPolicy.NO_DISAGREEMENT_SINGLE_MODEL
    calibration_production_enabled: bool = False

    # Applicability / confidence
    applicability_policy: str = "Chemical applicability domain evaluated per model."
    confidence_policy: str = "Model-reported confidence; not calibrated probability."

    # Experimental feedback
    experimental_feedback_policy: str = "No experimental feedback consumed in production."

    # Promotion
    promotion_status: PromotionStatus = PromotionStatus.ACTIVE
    shadow_promotion_status: PromotionStatus = PromotionStatus.DEFERRED
    promotion_requirements: List[str] = field(default_factory=list)
    rollback_target: Optional[str] = None
    rollback_reason: Optional[str] = None
    rollback_policy: Optional[RollbackPolicy] = None
    validation_status: ValidationStatus = ValidationStatus.PUBLISHER_REPORTED_ONLY
    scientific_status: str = "CONSERVATIVE_BASELINE"

    # Limitations
    limitations: List[str] = field(default_factory=list)

    # Audit trail
    scientific_notes: str = ""
    policy_version: str = REGISTRY_VERSION
    evidence_stage: str = "stage4d0"

    # Fallback
    fallback_behavior: str = "Return MODEL_UNAVAILABLE if primary model fails."

    @property
    def consensus_status(self) -> ConsensusPermission:
        """Canonical Stage 4D-4 name; legacy field remains for compatibility."""
        return self.consensus_permission

    @property
    def production_execution_allowed(self) -> bool:
        return self.primary_strategy != StrategyType.MODEL_UNAVAILABLE

    def to_dict(self) -> Dict[str, Any]:
        def serialize(value: Any) -> Any:
            if isinstance(value, enum.Enum):
                return value.value
            if isinstance(value, dict):
                return {key: serialize(item) for key, item in value.items()}
            if isinstance(value, (list, tuple)):
                return [serialize(item) for item in value]
            return value

        d = serialize(asdict(self))
        d["consensus_status"] = self.consensus_status.value
        d.pop("consensus_permission", None)
        d["production_execution_allowed"] = self.production_execution_allowed
        return d


# ──────────────────────────────────────────────────────────────────────────────
# Registry
# ──────────────────────────────────────────────────────────────────────────────

ENDPOINT_STRATEGY_REGISTRY: Dict[str, EndpointStrategyPolicy] = {}


def _register(policy: EndpointStrategyPolicy) -> None:
    contract = ENDPOINT_CONTRACTS.get(policy.endpoint_name)
    if contract is not None:
        policy.endpoint_contract_version = contract.version

    policy.promotion_requirements = list(dict.fromkeys(
        [*policy.promotion_requirements, *COMMON_PROMOTION_REQUIREMENTS]
    ))

    if policy.shadow_strategy is not None:
        policy.shadow_promotion_status = PromotionStatus.SHADOW

    if policy.promotion_status == PromotionStatus.ACTIVE and policy.rollback_policy is None:
        provenance = dict(zip(policy.primary_model_ids, policy.primary_model_versions))
        policy.rollback_policy = RollbackPolicy(
            previous_policy_version=f"{policy.evidence_stage}-production-baseline",
            rollback_target=policy.rollback_target or "MODEL_UNAVAILABLE",
            rollback_primary_strategy=policy.primary_strategy,
            rollback_model_ids=list(policy.primary_model_ids),
            rollback_model_versions=list(policy.primary_model_versions),
            promotion_reason=(
                "Stage 4D-4 codifies the existing production behavior; no production "
                "prediction strategy is activated or changed."
            ),
            validation_artifact=EVIDENCE_ARTIFACTS.get(
                policy.evidence_stage, "docs/current-platform-capabilities.md"
            ),
            model_version_provenance=provenance,
        )

    ENDPOINT_STRATEGY_REGISTRY[policy.endpoint_name] = policy


def get_endpoint_strategy(endpoint_name: str) -> Optional[EndpointStrategyPolicy]:
    direct = ENDPOINT_STRATEGY_REGISTRY.get(endpoint_name)
    if direct is not None:
        return direct
    return next(
        (policy for policy in ENDPOINT_STRATEGY_REGISTRY.values()
         if policy.endpoint_id == endpoint_name),
        None,
    )


def get_all_strategies() -> Dict[str, EndpointStrategyPolicy]:
    return dict(ENDPOINT_STRATEGY_REGISTRY)


# ──────────────────────────────────────────────────────────────────────────────
# Solubility
# ──────────────────────────────────────────────────────────────────────────────
_register(EndpointStrategyPolicy(
    endpoint_name="Solubility",
    endpoint_id="solubility_aqueous_logs",
    primary_strategy=PrimaryStrategy.SINGLE_CORE_MODEL,
    primary_model_ids=["admetica_solubility"],
    primary_model_versions=["admetica-d4f7056-chemprop-v2.1"],
    calibration_status=CalibrationStatus.NOT_APPLICABLE,
    calibration_policy="Regression output in log10(mol/L); no binary calibration.",
    shadow_strategy=StrategyType.ADAPTIVE_RESEARCH_SHADOW,
    shadow_model_ids=["esol_delaney_v1", "rdkit_gbr_solubility_v1"],
    shadow_model_versions=["esol-delaney-2004-v1.0", "rdkit-gbr-sol-v1.0"],
    non_primary_model_roles={
        "esol_delaney_v1": "SHADOW_RESEARCH",
        "rdkit_gbr_solubility_v1": "ADAPTIVE_EXCLUDED",
    },
    adaptive_status=AdaptiveStatus.ARCHITECTURE_VALID_BUT_NO_ACCURACY_GAIN,
    consensus_permission=ConsensusPermission.INSUFFICIENT_EVIDENCE,
    disagreement_policy=DisagreementPolicy.NO_DISAGREEMENT_SINGLE_MODEL,
    promotion_requirements=[
        "Demonstrated accuracy gain on held-out scaffold split vs M1 alone",
        "No assay-condition contamination across calibration/test",
        "Stable calibration across pH/temperature conditions",
    ],
    rollback_target="admetica_solubility (stage4d0)",
    validation_status=ValidationStatus.INDEPENDENT_PARTIAL,
    scientific_status="ARCHITECTURE_VALID_BUT_NO_ACCURACY_GAIN",
    limitations=[
        "Training aggregates heterogeneous pH/temperature/protocol conditions",
        "Not pH-specific; no intrinsic solubility separation",
        "Adaptive weighting architecture valid but no global accuracy gain (Stage 4D-3A2)",
    ],
    scientific_notes=(
        "Stage 4D-3A2 confirmed adaptive architecture is technically correct but "
        "adds no accuracy vs M1 CORE on held-out validation. M2 (ESOL), M3 (GBR) "
        "remain shadow-only. Adaptive remains RESEARCH_SHADOW."
    ),
    evidence_stage="stage4d3a2",
))

# ──────────────────────────────────────────────────────────────────────────────
# Permeability (Caco-2)
# ──────────────────────────────────────────────────────────────────────────────
_register(EndpointStrategyPolicy(
    endpoint_name="Permeability",
    endpoint_id="permeability_caco2_logpapp",
    primary_strategy=PrimaryStrategy.SINGLE_CORE_MODEL,
    primary_model_ids=["admetica_caco2"],
    primary_model_versions=["admetica-d4f7056-chemprop-v2.1"],
    calibration_status=CalibrationStatus.NOT_APPLICABLE,
    calibration_policy="Regression output in log10(cm/s); no binary calibration.",
    shadow_strategy=StrategyType.STATIC_CONSENSUS,
    shadow_model_ids=["physchem_caco2_v1"],
    shadow_model_versions=["physchem-caco2-v1.0"],
    non_primary_model_roles={"physchem_caco2_v1": "SHADOW_ONLY"},
    adaptive_status=AdaptiveStatus.DISABLED,
    consensus_permission=ConsensusPermission.INSUFFICIENT_EVIDENCE,
    disagreement_policy=DisagreementPolicy.NO_DISAGREEMENT_SINGLE_MODEL,
    promotion_requirements=[
        "Independent external validation cohort (N≥200) with A→B direction recorded",
        "Demonstrated ensemble improvement in MAE on holdout",
    ],
    rollback_target="admetica_caco2 (stage4d0)",
    validation_status=ValidationStatus.INDEPENDENT_PARTIAL,
    scientific_status="INSUFFICIENT_EVIDENCE",
    limitations=[
        "Training dataset lacks assay direction (A→B vs B→A) and pH metadata",
        "Not PAMPA, not MDCK, not efflux ratio",
        "Small external validation cohort (N=34) limits ensemble conclusions (Stage 4D-2C)",
    ],
    scientific_notes=(
        "Stage 4D-2C: KEEP_SHADOW for multi-model; insufficient external evidence to "
        "promote ensemble. M1 remains primary."
    ),
    evidence_stage="stage4d2c",
))

# ──────────────────────────────────────────────────────────────────────────────
# PPB
# ──────────────────────────────────────────────────────────────────────────────
_register(EndpointStrategyPolicy(
    endpoint_name="Plasma protein binding",
    endpoint_id="ppb_human_percent_bound",
    primary_strategy=PrimaryStrategy.SINGLE_CORE_MODEL,
    primary_model_ids=["admetica_ppbr"],
    primary_model_versions=["admetica-d4f7056-chemprop-v2.1"],
    calibration_status=CalibrationStatus.NOT_APPLICABLE,
    calibration_policy="Regression output (% bound); no binary calibration.",
    adaptive_status=AdaptiveStatus.DISABLED,
    consensus_permission=ConsensusPermission.DISABLED,
    disagreement_policy=DisagreementPolicy.NO_DISAGREEMENT_SINGLE_MODEL,
    rollback_target="admetica_ppbr (stage4d0)",
    validation_status=ValidationStatus.INDEPENDENT_PARTIAL,
    limitations=[
        "Single deterministic checkpoint; no uncertainty quantification",
        "Independent validation (Biogen N=185): MAE 14.6%; moderate accuracy",
        "Assay conditions not retained; cannot isolate total vs non-specific binding",
    ],
    scientific_notes="No ensemble evidence; SINGLE_CORE_MODEL appropriate.",
))

# ──────────────────────────────────────────────────────────────────────────────
# Deterministic physicochemical property pipeline
# ──────────────────────────────────────────────────────────────────────────────
_register(EndpointStrategyPolicy(
    endpoint_name="Physicochemical properties",
    endpoint_id="physchem_rdkit_2d_descriptors",
    endpoint_contract_version="stage1-rdkit-property-pipeline-v1",
    primary_strategy=StrategyType.DERIVED_ESTIMATE,
    primary_model_ids=["rdkit_2d_descriptors"],
    primary_model_versions=["rdkit-2025.03.1"],
    calibration_status=CalibrationStatus.NOT_APPLICABLE,
    calibration_policy="Deterministic structure descriptors; probability calibration does not apply.",
    adaptive_status=AdaptiveStatus.DISABLED,
    consensus_permission=ConsensusPermission.NOT_APPLICABLE,
    disagreement_policy=DisagreementPolicy.NOT_APPLICABLE,
    applicability_policy="Requires an RDKit-sanitizable standardized molecular structure.",
    confidence_policy="Calculated provenance is reported separately; no ML confidence is assigned.",
    rollback_target="RDKit property pipeline baseline",
    validation_status=ValidationStatus.MECHANISTIC,
    scientific_status="DETERMINISTIC_CALCULATION",
    limitations=[
        "Crippen cLogP is a calculated descriptor and is not logD at any pH",
        "QED and medicinal-chemistry rules are heuristic summaries, not validated clinical outcomes",
        "Descriptor calculations do not represent experimental measurements",
    ],
    scientific_notes="Covers the existing RDKit molecular formula, mass, cLogP, TPSA, HBD/HBA, rings, flexibility, Fsp3, molar refractivity, and QED outputs.",
    fallback_behavior="Return structure validation failure; do not fabricate property values for an invalid molecule.",
    evidence_stage="stage1",
))

# ──────────────────────────────────────────────────────────────────────────────
# Clearance: HLM / RLM / MLM — species-isolated single-core regressions
# ──────────────────────────────────────────────────────────────────────────────
for _ep, _eid, _species, _model_id in [
    ("HLM intrinsic clearance", "hlm_intrinsic_clearance_scaled_log10", "Human", "openadmet_hlm"),
    ("RLM intrinsic clearance", "rlm_intrinsic_clearance_scaled_log10", "Rat", "openadmet_rlm"),
    ("MLM intrinsic clearance", "mlm_intrinsic_clearance_scaled_log10", "Mouse", "openadmet_mlm"),
]:
    _register(EndpointStrategyPolicy(
        endpoint_name=_ep,
        endpoint_id=_eid,
        primary_strategy=PrimaryStrategy.SINGLE_CORE_MODEL,
        primary_model_ids=[_model_id],
        primary_model_versions=["openadmet-microsomal-clearance-chemeleon-v1-e135493"],
        calibration_status=CalibrationStatus.NOT_APPLICABLE,
        calibration_policy="Regression output in log10(mL/min/kg) after species-specific scaling; no binary calibration.",
        adaptive_status=AdaptiveStatus.DISABLED,
        consensus_permission=ConsensusPermission.INSUFFICIENT_EVIDENCE,
        disagreement_policy=DisagreementPolicy.NOT_APPLICABLE,
        applicability_policy="Species isolation enforced. HLM/RLM/MLM outputs must not be averaged.",
        promotion_requirements=[
            "Species-matched prospective validation cohort",
            "No pooling across species without allometric correction",
        ],
        rollback_target="openadmet_clearance (stage4d0)",
        validation_status=ValidationStatus.PUBLISHER_REPORTED_ONLY,
        limitations=[
            f"Species: {_species} microsomes only. Must not be mixed with other species.",
            "Species-specific scaled clearance prediction; averaging with other species is scientifically invalid.",
        ],
        scientific_notes=(
            "Species isolation is strict. HLM/RLM/MLM use the same OpenADMET multi-task checkpoint "
            "but species-specific task outputs. Averaging across species violates biology."
        ),
    ))

# ──────────────────────────────────────────────────────────────────────────────
# CYP Inhibitors — each isoform is a separate endpoint
# ──────────────────────────────────────────────────────────────────────────────

_CYP_INH_NOTES = {
    "CYP1A2 inhibitor": "No secondary model validated. Publisher validation only (Stage 4D-0).",
    "CYP2C9 inhibitor": "No secondary model validated. Publisher validation only (Stage 4D-0).",
    "CYP2C19 inhibitor": "No secondary model validated. Publisher validation only (Stage 4D-0).",
    "CYP2D6 inhibitor": "No secondary model validated. Publisher validation only (Stage 4D-0).",
    "CYP3A4 inhibitor": (
        "Stage 4D-3B1A: Full Adaptive does NOT outperform Fixed Global Prior. "
        "M2 role: CALIBRATION_SUPPORTING / SHADOW_ONLY. Dynamic adaptation DISABLED. "
        "Fixed global blend (w1=0.9578, w2=0.0422) remains FIXED_BLEND_RESEARCH shadow only."
    ),
}

for _ep, _eid, _isoform in [
    ("CYP1A2 inhibitor",  "cyp1a2_inhibitor_prob",  "CYP1A2"),
    ("CYP2C9 inhibitor",  "cyp2c9_inhibitor_prob",  "CYP2C9"),
    ("CYP2C19 inhibitor", "cyp2c19_inhibitor_prob", "CYP2C19"),
    ("CYP2D6 inhibitor",  "cyp2d6_inhibitor_prob",  "CYP2D6"),
    ("CYP3A4 inhibitor",  "cyp3a4_inhibitor_prob",  "CYP3A4"),
]:
    is_cyp3a4 = _isoform == "CYP3A4"
    _register(EndpointStrategyPolicy(
        endpoint_name=_ep,
        endpoint_id=_eid,
        primary_strategy=PrimaryStrategy.SINGLE_CORE_MODEL,
        primary_model_ids=[f"admetica_cyp_{_isoform.lower()}-inhibitor"],
        primary_model_versions=["admetica-d4f7056-cyp-chemprop-v2.1"],
        calibration_status=CalibrationStatus.RAW,
        calibration_policy="Raw binary classification probability; no post-hoc calibration in production.",
        decision_threshold=0.50,
        shadow_strategy=StrategyType.FIXED_WEIGHT_BLEND if is_cyp3a4 else None,
        shadow_model_ids=["morgan_cyp3a4_inh_v1"] if is_cyp3a4 else [],
        shadow_model_versions=["morgan-cyp3a4-v1.0"] if is_cyp3a4 else [],
        non_primary_model_roles=(
            {"morgan_cyp3a4_inh_v1": "FIXED_BLEND_RESEARCH_SHADOW_ONLY"}
            if is_cyp3a4 else {}
        ),
        adaptive_status=AdaptiveStatus.NO_ADAPTIVE_VALUE if is_cyp3a4 else AdaptiveStatus.DISABLED,
        consensus_permission=ConsensusPermission.INSUFFICIENT_EVIDENCE,
        disagreement_policy=(
            DisagreementPolicy.MODEL_DISAGREEMENT_SIGNAL if is_cyp3a4
            else DisagreementPolicy.NO_DISAGREEMENT_SINGLE_MODEL
        ),
        promotion_requirements=[
            "Independent held-out validation cohort (scaffold-split)",
            "Demonstrated calibration improvement on holdout",
            "No test-set leakage or assay contamination",
        ],
        rollback_target=f"admetica_cyp_{_isoform.lower()}-inhibitor (stage4d0)",
        validation_status=ValidationStatus.PUBLISHER_REPORTED_ONLY,
        scientific_status="FIXED_GLOBAL_BLEND_SUFFICIENT" if is_cyp3a4 else "INSUFFICIENT_EVIDENCE",
        limitations=[
            f"Isoform-specific: {_isoform}. Do not mix inhibitor/substrate endpoints.",
            "No probability calibration reported by publisher.",
            "Binary threshold optimized on PubChem AID 1851 training conditions.",
        ] + ([
            "Stage 4D-3B1A: adaptive weighting NO_GO; fixed global blend does not outperform "
            "Fixed Global Prior on Brier/MCC/BAcc; dynamic feedback adds no value beyond w=0.9578/0.0422.",
        ] if is_cyp3a4 else []),
        scientific_notes=_CYP_INH_NOTES[_ep],
        evidence_stage="stage4d3b1a" if is_cyp3a4 else "stage4d0",
    ))

# ──────────────────────────────────────────────────────────────────────────────
# CYP Substrates
# ──────────────────────────────────────────────────────────────────────────────
for _ep, _eid, _isoform in [
    ("CYP2C9 substrate",  "cyp2c9_substrate_prob",  "CYP2C9"),
    ("CYP2D6 substrate",  "cyp2d6_substrate_prob",  "CYP2D6"),
    ("CYP3A4 substrate",  "cyp3a4_substrate_prob",  "CYP3A4"),
]:
    _register(EndpointStrategyPolicy(
        endpoint_name=_ep,
        endpoint_id=_eid,
        primary_strategy=PrimaryStrategy.SINGLE_CORE_MODEL,
        primary_model_ids=[f"admetica_cyp_{_isoform.lower()}-substrate"],
        primary_model_versions=["admetica-d4f7056-cyp-chemprop-v2.1"],
        calibration_status=CalibrationStatus.RAW,
        calibration_policy="Raw binary classification probability; no post-hoc calibration.",
        decision_threshold=0.50,
        adaptive_status=AdaptiveStatus.DISABLED,
        consensus_permission=ConsensusPermission.INSUFFICIENT_EVIDENCE,
        disagreement_policy=DisagreementPolicy.NO_DISAGREEMENT_SINGLE_MODEL,
        promotion_requirements=[
            "Independent external substrate holdout (distinct from inhibitor training)",
            "Isoform-specific assay conditions documented",
        ],
        rollback_target=f"admetica_cyp_{_isoform.lower()}-substrate (stage4d0)",
        validation_status=ValidationStatus.PUBLISHER_REPORTED_ONLY,
        limitations=[
            f"Substrate vs inhibitor are distinct endpoints; must not be interchanged.",
            "Small training cohort; substrate N << inhibitor N.",
        ],
        scientific_notes=f"{_isoform} substrate; strictly distinct from inhibitor endpoint.",
    ))

# ──────────────────────────────────────────────────────────────────────────────
# P-gp Inhibitor
# ──────────────────────────────────────────────────────────────────────────────
_register(EndpointStrategyPolicy(
    endpoint_name="P-gp inhibitor",
    endpoint_id="transporter_pgp_inhibitor_prob",
    primary_strategy=PrimaryStrategy.SINGLE_CORE_MODEL,
    primary_model_ids=["admetica_transporter_pgp-inhibitor"],
    primary_model_versions=["admetica-d4f7056-pgp-inhibitor-chemprop-v2.1"],
    calibration_status=CalibrationStatus.RAW,
    calibration_policy="Raw binary probability; no production calibration.",
    decision_threshold=0.50,
    adaptive_status=AdaptiveStatus.DISABLED,
    consensus_permission=ConsensusPermission.INSUFFICIENT_EVIDENCE,
    disagreement_policy=DisagreementPolicy.NO_DISAGREEMENT_SINGLE_MODEL,
    rollback_target="admetica_transporter_pgp-inhibitor (stage4d0)",
    validation_status=ValidationStatus.NOT_VALIDATED,
    scientific_status="INDEPENDENT_VALIDATION_NOT_AVAILABLE",
    limitations=[
        "Heterogeneous training: inhibitor labels pooled across cell systems, probes, and IC50 cutoffs.",
        "No independent external validation qualified (Stage 4D-0).",
    ],
    scientific_notes="Independent validation not available; confidence capped LOW.",
))

# ──────────────────────────────────────────────────────────────────────────────
# hERG
# ──────────────────────────────────────────────────────────────────────────────
_register(EndpointStrategyPolicy(
    endpoint_name="hERG liability",
    endpoint_id="safety_herg_blocker_prob",
    primary_strategy=PrimaryStrategy.SINGLE_CORE_MODEL,
    primary_model_ids=["admetica_safety_herg"],
    primary_model_versions=["admetica-d4f7056-herg-chemprop-v2.1"],
    calibration_status=CalibrationStatus.CALIBRATION_RESEARCH,
    calibration_policy=(
        "Raw probability in production (threshold=0.50). Platt scaling on calibration "
        "subset shows ECE reduction 0.265→0.089 and specificity improvement 5.2%→31.0% "
        "on scaffold-split holdout. Calibration is RESEARCH stage; not yet activated."
    ),
    decision_threshold=0.50,
    shadow_strategy=StrategyType.SINGLE_CORE_WITH_CALIBRATION,
    shadow_model_ids=["physchem_herg_v1"],
    shadow_model_versions=["physchem-herg-v1.0"],
    non_primary_model_roles={
        "physchem_herg_v1": "CALIBRATION_SUPPORTING_SHADOW_ONLY_NOT_DISCRIMINATIVE_BLEND"
    },
    adaptive_status=AdaptiveStatus.NO_GO,
    consensus_permission=ConsensusPermission.INSUFFICIENT_EVIDENCE,
    disagreement_policy=DisagreementPolicy.MODEL_DISAGREEMENT_SIGNAL,
    experimental_feedback_policy="NO experimental feedback consumed. Adaptive weighting NO_GO (Stage 4D-3B2A).",
    promotion_requirements=[
        "Prospective calibration validation on separate holdout (no leakage from Stage 4D-3B2A data)",
        "Specificity ≥ 0.30 at production threshold on separate holdout",
        "Sensitivity ≥ 0.85 maintained on same holdout",
        "ECE < 0.10 after calibration on separate holdout",
        "New or qualified secondary model demonstrating reproducible M2 rescue rate > 20%",
        "Series heterogeneity analysis showing reproducible M2 advantage on ≥3 scaffold clusters",
    ],
    rollback_target="admetica_safety_herg raw (stage4d0)",
    validation_status=ValidationStatus.INDEPENDENT_COMPLETE_WITH_LIMITATIONS,
    scientific_status="HERG_CALIBRATION_UPDATE_CANDIDATE",
    limitations=[
        "AUROC=0.667 — moderate base discrimination; cannot reach safety-screening specificity targets by calibration alone",
        "Training prevalence 86% vs evaluation 67%: 18.8 pp prior shift (CLASS_IMBALANCE)",
        "72.9% of evaluation compounds are IC50 borderline (1k–30k nM): LABEL_BOUNDARY_UNCERTAINTY",
        "ASSAY_HETEROGENEITY_PRESENT: Wang et al. pools patch-clamp + radioligand binding",
        "M2 rescue rate 5.4%: BETTER_SECONDARY_MODEL_REQUIRED",
        "FP/FN ratio 17.7:1 at threshold 0.50",
    ],
    scientific_notes=(
        "Stage 4D-3B2A corrected interpretation: HERG_CALIBRATION_UPDATE_CANDIDATE. "
        "Cal-selected best blend=100/0 (pure M1). Do NOT retain HERG_FIXED_BLEND_CANDIDATE. "
        "Primary limitations: calibration + prior shift + label boundary uncertainty. "
        "Adaptive weighting: NO_GO. BETTER_SECONDARY_MODEL_REQUIRED. "
        "Production remains unchanged until Stage 4D-5 promotion cycle."
    ),
    evidence_stage="stage4d3b2a",
    fallback_behavior="Return raw M1 probability with confidence=LOW; never suppress output.",
))

# ──────────────────────────────────────────────────────────────────────────────
# Ames mutagenicity
# ──────────────────────────────────────────────────────────────────────────────
_register(EndpointStrategyPolicy(
    endpoint_name="Ames mutagenicity",
    endpoint_id="safety_ames_mutagenicity_prob",
    primary_strategy=PrimaryStrategy.SINGLE_CORE_MODEL,
    primary_model_ids=["admet_ai_ames"],
    primary_model_versions=["admet-ai-v2.0.1-c65bf04-chemprop-v2-ensemble5"],
    calibration_status=CalibrationStatus.RAW,
    calibration_policy="Raw binary probability; no post-hoc calibration in production.",
    decision_threshold=0.50,
    adaptive_status=AdaptiveStatus.DISABLED,
    consensus_permission=ConsensusPermission.INSUFFICIENT_EVIDENCE,
    disagreement_policy=DisagreementPolicy.NO_DISAGREEMENT_SINGLE_MODEL,
    promotion_requirements=[
        "Independent held-out Ames test set (separate strains/conditions from ADMET-AI training)",
        "Calibration evidence on holdout",
    ],
    rollback_target="admet_ai_ames (stage4d0)",
    validation_status=ValidationStatus.PUBLISHER_REPORTED_ONLY,
    limitations=[
        "Training mixes TA98, TA100, TA1535, TA1537, TA102 with/without S9; strains not per-row",
        "Positive and negative cutoffs vary by strain; aggregate label is a simplified binary",
        "Distinct from mammalian micronucleus or chromosomal aberration endpoints",
    ],
    scientific_notes=(
        "ADMET-AI v2 ensemble (5 models). Ames findings must not be extrapolated to "
        "mammalian genotoxicity endpoints. Strains pooled."
    ),
    evidence_stage="stage3f",
))

# ──────────────────────────────────────────────────────────────────────────────
# DILI
# ──────────────────────────────────────────────────────────────────────────────
_register(EndpointStrategyPolicy(
    endpoint_name="DILI clinical liability",
    endpoint_id="safety_dili_clinical_prob",
    primary_strategy=PrimaryStrategy.SINGLE_CORE_MODEL,
    primary_model_ids=["admet_ai_dili"],
    primary_model_versions=["admet-ai-v2.0.1-c65bf04-chemprop-v2-ensemble5"],
    calibration_status=CalibrationStatus.RAW,
    calibration_policy="Raw binary probability; no production calibration.",
    decision_threshold=0.50,
    adaptive_status=AdaptiveStatus.DISABLED,
    consensus_permission=ConsensusPermission.INSUFFICIENT_EVIDENCE,
    disagreement_policy=DisagreementPolicy.NO_DISAGREEMENT_SINGLE_MODEL,
    rollback_target="admet_ai_dili (stage4d0)",
    validation_status=ValidationStatus.PUBLISHER_REPORTED_ONLY,
    limitations=[
        "Small training set (N=475); DILI classification is clinically heterogeneous",
        "Clinical DILI labels aggregate diverse mechanisms; not mechanism-specific",
        "High clinical uncertainty; output is a screening liability, not causal assessment",
    ],
    scientific_notes="ADMET-AI v2 ensemble. DILI is a complex multi-mechanism endpoint; screening use only.",
    evidence_stage="stage3f",
))

# ──────────────────────────────────────────────────────────────────────────────
# Metabolic Soft Spots (SoM) — RANK_FUSION
# ──────────────────────────────────────────────────────────────────────────────
_register(EndpointStrategyPolicy(
    endpoint_name="Metabolic soft spots",
    endpoint_id="som_metabolic_soft_spots",
    primary_strategy=PrimaryStrategy.RANK_FUSION,
    primary_model_ids=["sygma_phase1_2", "smartcyp_dft_v1"],
    primary_model_versions=["sygma-v1.1.0", "smartcyp-v3.0"],
    calibration_status=CalibrationStatus.NOT_APPLICABLE,
    calibration_policy="Rank-based output; not a probability. Reciprocal-rank fusion of atom-level scores.",
    adaptive_status=AdaptiveStatus.DISABLED,
    consensus_permission=ConsensusPermission.ALLOWED_PRODUCTION,
    disagreement_policy=DisagreementPolicy.NOT_APPLICABLE,
    applicability_policy=(
        "SyGMa: SMARTS-based biotransformation rules; SMARTCyp: DFT-derived fragment energies. "
        "Raw model scores are NOT averaged; reciprocal-rank fusion is applied per atom."
    ),
    promotion_requirements=[
        "Per-atom prospective SoM validation cohort",
        "Independent comparison of RRF vs individual model rank accuracy",
    ],
    rollback_target="sygma_phase1_2@sygma-v1.1.0 + smartcyp_dft_v1@smartcyp-v3.0 rank-fusion baseline",
    validation_status=ValidationStatus.MECHANISTIC,
    limitations=[
        "SyGMa applies reaction rules; coverage limited to known biotransformations",
        "SMARTCyp covers CYP-mediated oxidation; other pathways (UGT, SULT) not included",
        "Raw metabolic scores are mechanistically distinct; do NOT average as scalar values",
    ],
    scientific_notes=(
        "RANK_FUSION is scientifically valid here: both models predict atom-level metabolic "
        "susceptibility ranks. Reciprocal-rank fusion preserves per-model ranks and combines "
        "them without mixing raw probability scales."
    ),
    evidence_stage="stage3f",
))

_register(EndpointStrategyPolicy(
    endpoint_name="Metabolite hypotheses",
    endpoint_id="metabolism_sygma_rule_hypotheses",
    endpoint_contract_version="stage3d-sygma-v1",
    primary_strategy=StrategyType.RULE_BASED,
    primary_model_ids=["sygma_phase1_2"],
    primary_model_versions=["sygma-v1.1.0"],
    calibration_status=CalibrationStatus.NOT_APPLICABLE,
    calibration_policy="Empirical biotransformation rule ranking; not a calibrated probability.",
    adaptive_status=AdaptiveStatus.DISABLED,
    consensus_permission=ConsensusPermission.NOT_APPLICABLE,
    disagreement_policy=DisagreementPolicy.NOT_APPLICABLE,
    applicability_policy="Limited to supported SyGMa Phase I/II transformation rules and chemically valid products.",
    confidence_policy="Rule evidence and rank are reported; no atom-level ML confidence is implied.",
    rollback_target="SyGMa stage3d rule baseline",
    validation_status=ValidationStatus.MECHANISTIC,
    scientific_status="RULE_BASED_METABOLITE_HYPOTHESIS",
    limitations=[
        "Outputs are predicted metabolite hypotheses, not experimentally confirmed metabolites",
        "No qualified atom-level ML probability is emitted",
        "CYP compound-level evidence cannot assign an atom or transformation to an isoform",
    ],
    scientific_notes="SMARTCyp participates only in compatible SoM rank fusion; SyGMa rules generate the existing metabolite hypotheses.",
    fallback_behavior="Return no supported hypothesis when no qualified transformation rule matches.",
    evidence_stage="stage3f",
))

# ──────────────────────────────────────────────────────────────────────────────
# PK: Systemic Clearance, Vd, Bioavailability — MECHANISTIC_NO_CONSENSUS
# ──────────────────────────────────────────────────────────────────────────────
for _ep, _eid in [
    ("PK Systemic Clearance",         "pk_clearance_systemic"),
    ("PK Volume of Distribution",     "pk_volume_distribution_vss"),
    ("PK Bioavailability",            "pk_bioavailability_fraction"),
]:
    _register(EndpointStrategyPolicy(
        endpoint_name=_ep,
        endpoint_id=_eid,
        primary_strategy=PrimaryStrategy.MECHANISTIC_NO_CONSENSUS,
        primary_model_ids=["pk_nca", "pk_ivive", "pk_allometric_scaling"],
        primary_model_versions=["5A-1.0", "5A-2A.1.0", "5B-3.1.0"],
        calibration_status=CalibrationStatus.NOT_APPLICABLE,
        calibration_policy="Mechanistic/allometric estimate; not a probabilistic prediction.",
        adaptive_status=AdaptiveStatus.DISABLED,
        consensus_permission=ConsensusPermission.MECHANISTICALLY_FORBIDDEN,
        disagreement_policy=DisagreementPolicy.NOT_APPLICABLE,
        applicability_policy=(
            "Evidence hierarchy: (1) NCA from in vivo data; (2) IVIVE from microsomal clearance; "
            "(3) allometric scaling from species data. Each method has distinct assumptions."
        ),
        rollback_target="mechanistic_pk_baseline (stage5a1)",
        validation_status=ValidationStatus.MECHANISTIC,
        limitations=[
            "NCA, IVIVE, and compartment simulation are mechanistically incompatible; must not be averaged",
            "IVIVE requires microsomal Clint as input; HLM/RLM/MLM predictions carry uncertainty",
            "Allometric scaling assumes geometric scaling law; extrapolation uncertainty not bounded",
        ],
        scientific_notes=(
            "PK endpoints use an evidence hierarchy: NCA (if in vivo data available) > IVIVE > allometry. "
            "These are NOT interchangeable ML outputs and must NOT be averaged as an ensemble."
        ),
        evidence_stage="stage5a1",
    ))

_register(EndpointStrategyPolicy(
    endpoint_name="PK Simulation",
    endpoint_id="pk_concentration_time_simulation",
    endpoint_contract_version="stage5b2-pk-simulation-v1",
    primary_strategy=StrategyType.MECHANISTIC_NO_CONSENSUS,
    primary_model_ids=["pk_simulation"],
    primary_model_versions=["5B-2.0"],
    calibration_status=CalibrationStatus.NOT_APPLICABLE,
    calibration_policy="Numerical mechanistic simulation; probability calibration does not apply.",
    adaptive_status=AdaptiveStatus.DISABLED,
    consensus_permission=ConsensusPermission.MECHANISTICALLY_FORBIDDEN,
    disagreement_policy=DisagreementPolicy.NOT_APPLICABLE,
    applicability_policy="Requires route-compatible PK parameters with explicit experimental or derived provenance.",
    confidence_policy="Simulation input evidence and parameter provenance remain separate from numerical output.",
    rollback_target="PK simulation 5B-2.0 baseline",
    validation_status=ValidationStatus.MECHANISTIC,
    scientific_status="MECHANISTIC_SIMULATION",
    limitations=[
        "Simulation does not create missing CL, V, F, or ka evidence",
        "Experimental, NCA, IVIVE, allometric, and derived inputs retain distinct assumptions",
        "Mechanistic methods must not be averaged as ML predictions",
    ],
    scientific_notes="Stage 4D-4 governs the existing simulation method without altering profiles or parameter selection.",
    fallback_behavior="Return MODEL_UNAVAILABLE when required route and PK parameters are unavailable.",
    evidence_stage="stage5b2",
))

# ──────────────────────────────────────────────────────────────────────────────
# pKa / Ionization — RULE_ESTIMATE
# ──────────────────────────────────────────────────────────────────────────────
_register(EndpointStrategyPolicy(
    endpoint_name="Ionization (pKa)",
    endpoint_id="ionization_pka_estimated",
    endpoint_contract_version="stage4c4-ionization-v1",
    primary_strategy=PrimaryStrategy.RULE_ESTIMATE,
    primary_model_ids=["ionization_smarts_rules_v1"],
    primary_model_versions=["stage4c4-ionization-v1"],
    calibration_status=CalibrationStatus.NOT_APPLICABLE,
    calibration_policy=(
        "SMARTS-based structural classification with literature range pKa estimates. "
        "This is a rule estimate, not a validated ML pKa prediction."
    ),
    adaptive_status=AdaptiveStatus.DISABLED,
    consensus_permission=ConsensusPermission.NOT_APPLICABLE,
    disagreement_policy=DisagreementPolicy.NOT_APPLICABLE,
    promotion_requirements=[
        "Replacement with validated pKa ML model (e.g. ACD/pKa, Epik, or open equivalent)",
        "Per-ionization-center prospective validation cohort",
        "Assay type and ionic strength recorded",
    ],
    rollback_target="structural_ionization_smarts (stage4c4)",
    validation_status=ValidationStatus.NOT_VALIDATED,
    scientific_status="RULE_ESTIMATE_NOT_VALIDATED_ML",
    limitations=[
        "RULE_ESTIMATE: not a quantitatively validated ML model",
        "Estimates derived from literature SMARTS pKa ranges; accuracy ± 1–2 pKa units",
        "Multi-protic and complex ampholyte systems have high uncertainty",
        "Must not be presented to users as a validated quantitative pKa prediction",
    ],
    scientific_notes=(
        "pKa is computed via SMARTS pattern matching and representative literature ranges. "
        "This is scientifically honest: it is a rule estimate, not an ML model. "
        "A validated quantitative pKa model is a future requirement (Stage 4D-5+)."
    ),
    evidence_stage="stage4c4",
    fallback_behavior="Return ionization class with pKa=MODEL_UNAVAILABLE when no rule or experiment applies.",
))

# ──────────────────────────────────────────────────────────────────────────────
# logD(pH 7.4) — explicit derived estimate, never a validated ML prediction
# ──────────────────────────────────────────────────────────────────────────────
_register(EndpointStrategyPolicy(
    endpoint_name="logD pH7.4 derived estimate",
    endpoint_id="physchem_logd_ph74_derived_estimate",
    endpoint_contract_version="stage4c4-ionization-v1",
    primary_strategy=StrategyType.DERIVED_ESTIMATE,
    primary_model_ids=["henderson_hasselbalch_logd_v1"],
    primary_model_versions=["stage4c4-ionization-v1"],
    calibration_status=CalibrationStatus.NOT_APPLICABLE,
    calibration_policy=(
        "Derived from RDKit Crippen cLogP plus rule/experimental pKa under explicit "
        "monoprotic Henderson-Hasselbalch assumptions; not an ML logD prediction."
    ),
    adaptive_status=AdaptiveStatus.DISABLED,
    consensus_permission=ConsensusPermission.NOT_APPLICABLE,
    disagreement_policy=DisagreementPolicy.NOT_APPLICABLE,
    applicability_policy="Only emitted when the ionization assumptions required by the derivation are available.",
    confidence_policy="DERIVED_ESTIMATE provenance is separate from model confidence.",
    rollback_target="stage4c4-ionization-derived-logd",
    validation_status=ValidationStatus.NOT_VALIDATED,
    scientific_status="DERIVED_ESTIMATE_NOT_VALIDATED_ML",
    limitations=[
        "Not a standalone quantitative ML logD7.4 model",
        "Monoprotic assumptions can fail for polyprotic or zwitterionic compounds",
        "Calculated cLogP is never relabeled as experimental or validated predicted logD",
    ],
    scientific_notes="The existing runtime ionization engine exposes this conditional derivation; Stage 4D-4 does not change it.",
    evidence_stage="stage4c4",
    fallback_behavior="Return MODEL_UNAVAILABLE when pKa, cLogP, or derivation assumptions are unavailable.",
))


# These names mirror the inactive entries created by backend.admet.ensure_admet_schema.
UNAVAILABLE_RUNTIME_ENDPOINTS = [
    ("Microsomal clearance", "microsomal_clearance_generic", "stage3b"),
    ("Dog liver microsomal intrinsic clearance", "dlm_intrinsic_clearance", "stage3b"),
    ("Monkey liver microsomal intrinsic clearance", "cylm_intrinsic_clearance", "stage3b"),
    ("CYP1A2 substrate", "cyp1a2_substrate_prob", "stage3c"),
    ("CYP2C19 substrate", "cyp2c19_substrate_prob", "stage3c"),
    ("P-gp substrate", "transporter_pgp_substrate_prob", "stage3e"),
    ("BCRP substrate", "transporter_bcrp_substrate_prob", "stage3e"),
    ("BCRP inhibitor", "transporter_bcrp_inhibitor_prob", "stage3e"),
    ("BSEP inhibitor", "transporter_bsep_inhibitor_prob", "stage3e"),
    ("OATP1B1 inhibitor", "transporter_oatp1b1_inhibitor_prob", "stage3e"),
    ("OATP1B3 inhibitor", "transporter_oatp1b3_inhibitor_prob", "stage3e"),
    ("OCT1 inhibitor", "transporter_oct1_inhibitor_prob", "stage3e"),
    ("OCT2 inhibitor", "transporter_oct2_inhibitor_prob", "stage3e"),
    ("MATE1 inhibitor", "transporter_mate1_inhibitor_prob", "stage3e"),
    ("MATE2-K inhibitor", "transporter_mate2k_inhibitor_prob", "stage3e"),
    ("Mitochondrial toxicity", "safety_mitochondrial_toxicity", "stage3f"),
    ("General cytotoxicity", "safety_general_cytotoxicity", "stage3f"),
    ("Skin sensitization", "safety_skin_sensitization", "stage3f"),
    ("BBB penetration", "distribution_bbb_penetration", "stage3f"),
    ("CNS liability", "safety_cns_liability", "stage3f"),
    ("pKa (quantitative ML)", "physchem_pka_quantitative_ml", "stage4c4"),
    ("logD7.4 (quantitative ML)", "physchem_logd74_quantitative_ml", "stage4c4"),
]


for _ep, _eid, _evidence_stage in UNAVAILABLE_RUNTIME_ENDPOINTS:
    _register(EndpointStrategyPolicy(
        endpoint_name=_ep,
        endpoint_id=_eid,
        endpoint_contract_version=f"unavailable-{_evidence_stage}",
        primary_strategy=StrategyType.MODEL_UNAVAILABLE,
        primary_model_ids=[],
        primary_model_versions=[],
        calibration_status=CalibrationStatus.NOT_APPLICABLE,
        calibration_policy="No qualified model; no calibration applicable.",
        adaptive_status=AdaptiveStatus.DISABLED,
        consensus_permission=ConsensusPermission.NOT_APPLICABLE,
        disagreement_policy=DisagreementPolicy.NOT_APPLICABLE,
        promotion_status=PromotionStatus.DEFERRED,
        promotion_requirements=[
            "Qualified endpoint-specific model with compatible license and runtime",
            "Independent external validation cohort",
            "Endpoint contract defined and compatible with existing governance",
        ],
        validation_status=ValidationStatus.NOT_VALIDATED,
        scientific_status="MODEL_UNAVAILABLE",
        fallback_behavior="Return MODEL_UNAVAILABLE. Do not synthesize values from AI, LLM, cross-endpoint reuse, or heuristics.",
        limitations=["No scientifically qualified endpoint-specific production model is available."],
        scientific_notes=f"{_ep} is an explicit runtime registry slot but is inactive. No production value may be emitted.",
        evidence_stage=_evidence_stage,
    ))


ACTIVE_ADMET_MODEL_ENDPOINTS = (
    "Solubility", "Permeability", "Plasma protein binding",
    "HLM intrinsic clearance", "RLM intrinsic clearance", "MLM intrinsic clearance",
    "CYP1A2 inhibitor", "CYP2C9 inhibitor", "CYP2C19 inhibitor",
    "CYP2D6 inhibitor", "CYP3A4 inhibitor",
    "CYP2C9 substrate", "CYP2D6 substrate", "CYP3A4 substrate",
    "P-gp inhibitor", "hERG liability", "Ames mutagenicity", "DILI clinical liability",
)

RUNTIME_ADMET_MODEL_ENDPOINTS = frozenset(
    [*ACTIVE_ADMET_MODEL_ENDPOINTS, *(_ep for _ep, _, _ in UNAVAILABLE_RUNTIME_ENDPOINTS)]
)


# ──────────────────────────────────────────────────────────────────────────────
# Registry validation
# ──────────────────────────────────────────────────────────────────────────────

def validate_registry() -> List[str]:
    """Run internal consistency checks. Returns list of violations."""
    violations: List[str] = []
    endpoint_ids: Dict[str, str] = {}

    missing_contracts = set(ENDPOINT_CONTRACTS) - set(ENDPOINT_STRATEGY_REGISTRY)
    for name in sorted(missing_contracts):
        violations.append(f"{name}: endpoint contract has no strategy policy")

    missing_runtime = set(RUNTIME_ADMET_MODEL_ENDPOINTS) - set(ENDPOINT_STRATEGY_REGISTRY)
    for name in sorted(missing_runtime):
        violations.append(f"{name}: runtime model registry endpoint has no strategy policy")

    for name, policy in ENDPOINT_STRATEGY_REGISTRY.items():
        if policy.endpoint_id in endpoint_ids:
            violations.append(
                f"{name}: duplicate endpoint_id also used by {endpoint_ids[policy.endpoint_id]}"
            )
        endpoint_ids[policy.endpoint_id] = name

        contract = ENDPOINT_CONTRACTS.get(name)
        if contract is not None:
            if policy.endpoint_id != contract.endpoint_id:
                violations.append(f"{name}: endpoint_id differs from endpoint contract")
            if policy.endpoint_contract_version != contract.version:
                violations.append(f"{name}: endpoint_contract_version differs from endpoint contract")

        if not isinstance(policy.primary_strategy, StrategyType):
            violations.append(f"{name}: invalid primary strategy enum")
        if policy.shadow_strategy is not None and not isinstance(policy.shadow_strategy, StrategyType):
            violations.append(f"{name}: invalid shadow strategy enum")
        if len(policy.primary_model_ids) != len(policy.primary_model_versions):
            violations.append(f"{name}: primary model IDs/versions are not one-to-one")
        if len(policy.shadow_model_ids) != len(policy.shadow_model_versions):
            violations.append(f"{name}: shadow model IDs/versions are not one-to-one")
        if not set(policy.non_primary_model_roles) <= set(policy.shadow_model_ids):
            violations.append(f"{name}: non-primary model role lacks shadow model provenance")

        # MODEL_UNAVAILABLE must have empty model lists
        if policy.primary_strategy == StrategyType.MODEL_UNAVAILABLE:
            if policy.primary_model_ids or policy.primary_model_versions:
                violations.append(f"{name}: MODEL_UNAVAILABLE but has primary_model_ids")
            if policy.production_execution_allowed:
                violations.append(f"{name}: MODEL_UNAVAILABLE permits production execution")
            if policy.promotion_status == PromotionStatus.ACTIVE:
                violations.append(f"{name}: MODEL_UNAVAILABLE cannot be ACTIVE")
            if policy.shadow_strategy is not None:
                violations.append(f"{name}: MODEL_UNAVAILABLE cannot have a shadow strategy")

        # SINGLE_CORE_MODEL must have exactly one model
        if policy.primary_strategy == StrategyType.SINGLE_CORE_MODEL:
            if len(policy.primary_model_ids) != 1:
                violations.append(
                    f"{name}: SINGLE_CORE_MODEL must have exactly 1 primary_model_id, "
                    f"got {len(policy.primary_model_ids)}"
                )

        # MECHANISTIC_NO_CONSENSUS must not have consensus enabled
        if policy.primary_strategy == StrategyType.MECHANISTIC_NO_CONSENSUS:
            if policy.consensus_permission not in (
                ConsensusPermission.MECHANISTICALLY_FORBIDDEN,
                ConsensusPermission.DISABLED,
                ConsensusPermission.NOT_APPLICABLE,
            ):
                violations.append(
                    f"{name}: MECHANISTIC_NO_CONSENSUS must have MECHANISTICALLY_FORBIDDEN consensus"
                )

        # Adaptive production remains impossible when policy is disabled/no-go.
        if policy.adaptive_status in (AdaptiveStatus.NO_GO, AdaptiveStatus.NO_ADAPTIVE_VALUE,
                                       AdaptiveStatus.DISABLED):
            if policy.shadow_strategy == StrategyType.ADAPTIVE_RESEARCH_SHADOW:
                violations.append(
                    f"{name}: adaptive_status={policy.adaptive_status.value} "
                    f"but shadow_strategy=ADAPTIVE_RESEARCH_SHADOW"
                )

        if policy.calibration_status == CalibrationStatus.CALIBRATION_RESEARCH:
            if policy.calibration_production_enabled:
                violations.append(f"{name}: research calibration cannot be production-enabled")

        if policy.shadow_strategy is None:
            if policy.shadow_model_ids or policy.shadow_model_versions:
                violations.append(f"{name}: shadow models exist without a shadow strategy")
        elif policy.shadow_promotion_status != PromotionStatus.SHADOW:
            violations.append(f"{name}: shadow strategy must remain in SHADOW lifecycle state")

        # Calibration status for regression endpoints
        if policy.decision_threshold is None and policy.calibration_status == CalibrationStatus.RAW:
            # For regression endpoints RAW is OK; not a violation
            pass

        # RANK_FUSION must have ≥ 2 model IDs
        if policy.primary_strategy == StrategyType.RANK_FUSION:
            if len(policy.primary_model_ids) < 2:
                violations.append(f"{name}: RANK_FUSION requires ≥ 2 primary_model_ids")

        if policy.promotion_status == PromotionStatus.ACTIVE:
            rollback = policy.rollback_policy
            if rollback is None:
                violations.append(f"{name}: ACTIVE policy lacks rollback_policy")
            elif not all((
                rollback.previous_policy_version,
                rollback.rollback_target,
                rollback.promotion_reason,
                rollback.validation_artifact,
            )):
                violations.append(f"{name}: rollback_policy is incomplete")
            missing_gates = set(COMMON_PROMOTION_REQUIREMENTS) - set(policy.promotion_requirements)
            if missing_gates:
                violations.append(f"{name}: promotion requirements omit mandatory governance gates")

    # Endpoint-specific scientific invariants.
    sol = ENDPOINT_STRATEGY_REGISTRY.get("Solubility")
    if not sol or sol.primary_strategy != StrategyType.SINGLE_CORE_MODEL or sol.shadow_strategy != StrategyType.ADAPTIVE_RESEARCH_SHADOW:
        violations.append("Solubility: production must remain single-core with adaptive research shadow")
    elif sol.non_primary_model_roles.get("rdkit_gbr_solubility_v1") != "ADAPTIVE_EXCLUDED":
        violations.append("Solubility: M3 must remain adaptive-excluded")

    cyp3a4 = ENDPOINT_STRATEGY_REGISTRY.get("CYP3A4 inhibitor")
    if not cyp3a4 or cyp3a4.primary_strategy != StrategyType.SINGLE_CORE_MODEL or cyp3a4.adaptive_status != AdaptiveStatus.NO_ADAPTIVE_VALUE:
        violations.append("CYP3A4 inhibitor: dynamic adaptation must remain disabled")

    herg = ENDPOINT_STRATEGY_REGISTRY.get("hERG liability")
    if not herg or herg.primary_strategy != StrategyType.SINGLE_CORE_MODEL:
        violations.append("hERG liability: production must remain SINGLE_CORE_MODEL")
    elif (herg.calibration_status != CalibrationStatus.CALIBRATION_RESEARCH
          or herg.calibration_production_enabled
          or "physchem_herg_v1" in herg.primary_model_ids):
        violations.append("hERG liability: calibration/M2 must remain research-only")
    elif "NOT_DISCRIMINATIVE_BLEND" not in herg.non_primary_model_roles.get("physchem_herg_v1", ""):
        violations.append("hERG liability: M2 must not be represented as a discriminative blend")

    som = ENDPOINT_STRATEGY_REGISTRY.get("Metabolic soft spots")
    if not som or som.primary_strategy != StrategyType.RANK_FUSION:
        violations.append("Metabolic soft spots: strategy must be RANK_FUSION")

    for pk_name in ("PK Systemic Clearance", "PK Volume of Distribution", "PK Bioavailability"):
        pk_policy = ENDPOINT_STRATEGY_REGISTRY.get(pk_name)
        if not pk_policy or pk_policy.primary_strategy != StrategyType.MECHANISTIC_NO_CONSENSUS:
            violations.append(f"{pk_name}: PK mechanisms must not enter ML consensus")

    pka = ENDPOINT_STRATEGY_REGISTRY.get("Ionization (pKa)")
    logd = ENDPOINT_STRATEGY_REGISTRY.get("logD pH7.4 derived estimate")
    if not pka or pka.primary_strategy != StrategyType.RULE_ESTIMATE:
        violations.append("Ionization (pKa): provenance must remain RULE_ESTIMATE")
    if not logd or logd.primary_strategy != StrategyType.DERIVED_ESTIMATE:
        violations.append("logD pH7.4 derived estimate: provenance must remain DERIVED_ESTIMATE")

    clearance_models = {
        ENDPOINT_STRATEGY_REGISTRY[name].primary_model_ids[0]
        for name in ("HLM intrinsic clearance", "RLM intrinsic clearance", "MLM intrinsic clearance")
    }
    if clearance_models != {"openadmet_hlm", "openadmet_rlm", "openadmet_mlm"}:
        violations.append("Microsomal clearance: species-specific model outputs are not isolated")

    return violations


def get_registry_summary() -> Dict[str, Any]:
    """Return a serializable summary of the full registry."""
    registry = get_all_strategies()
    violations = validate_registry()
    return {
        "registry_version": REGISTRY_VERSION,
        "policy_date": REGISTRY_POLICY_DATE,
        "total_endpoints": len(registry),
        "violations": violations,
        "endpoints": {k: v.to_dict() for k, v in registry.items()},
    }


def get_registry_api_response() -> Dict[str, Any]:
    """Stable read-only API projection; does not participate in prediction execution."""
    summary = get_registry_summary()
    endpoints = []
    for name, policy in sorted(ENDPOINT_STRATEGY_REGISTRY.items()):
        row = policy.to_dict()
        endpoints.append({
            "endpoint": name,
            "endpoint_id": policy.endpoint_id,
            "endpoint_contract_version": policy.endpoint_contract_version,
            "primary_strategy": policy.primary_strategy.value,
            "primary_models": [
                {"model_id": model_id, "model_version": version}
                for model_id, version in zip(policy.primary_model_ids, policy.primary_model_versions)
            ],
            "calibration_status": policy.calibration_status.value,
            "calibration_policy": policy.calibration_policy,
            "calibration_production_enabled": policy.calibration_production_enabled,
            "decision_threshold": policy.decision_threshold,
            "shadow_strategy": policy.shadow_strategy.value if policy.shadow_strategy else None,
            "shadow_models": [
                {"model_id": model_id, "model_version": version}
                for model_id, version in zip(policy.shadow_model_ids, policy.shadow_model_versions)
            ],
            "non_primary_model_roles": dict(policy.non_primary_model_roles),
            "adaptive_status": policy.adaptive_status.value,
            "consensus_status": policy.consensus_status.value,
            "applicability_policy": policy.applicability_policy,
            "confidence_policy": policy.confidence_policy,
            "disagreement_policy": policy.disagreement_policy.value,
            "experimental_feedback_policy": policy.experimental_feedback_policy,
            "validation_status": policy.validation_status.value,
            "scientific_status": policy.scientific_status,
            "limitations": list(policy.limitations),
            "policy_version": policy.policy_version,
            "promotion_status": policy.promotion_status.value,
            "shadow_promotion_status": policy.shadow_promotion_status.value,
            "promotion_requirements": list(policy.promotion_requirements),
            "production_execution_allowed": policy.production_execution_allowed,
            "fallback_behavior": policy.fallback_behavior,
            "rollback_policy": row["rollback_policy"],
            "scientific_notes": policy.scientific_notes,
        })
    return {
        "registry_version": summary["registry_version"],
        "policy_date": summary["policy_date"],
        "total_endpoints": summary["total_endpoints"],
        "violations": summary["violations"],
        "read_only": True,
        "production_behavior_changed": False,
        "endpoints": endpoints,
    }
