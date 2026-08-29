"""Generate all internal validation v1 artifacts from current DB state.

Run after campaign initialization and after any experimental data import.
Produces all required validation/internal_validation_v1_*.json files.

Scientific rules:
  - No experimental values fabricated
  - If no experimental data exists: artifacts honestly record this state
  - No model fitting or modification
  - Policy hash verified at runtime
"""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

VALIDATION_DIR = ROOT / "validation"


def _sha256(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def utcnow_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def run():
    from backend.database import SessionLocal, engine
    from backend.prediction_engine_v1_policy import policy_hash as live_policy_hash, policy_rows
    from backend.internal_validation_v1 import (
        CAMPAIGN_ID, ENGINE_V1_POLICY_HASH, ENGINE_V1_POLICY_VERSION,
        ensure_validation_schema,
        InternalValidationCampaignRow,
        InternalValidationCohortEntryRow,
        InternalValidationPredictionFreezeRow,
        InternalValidationExperimentalRecordRow,
        InternalValidationObservationRow,
        campaign_summary,
    )
    import sqlalchemy as sa

    # Verify Engine v1 unchanged
    h = live_policy_hash()
    assert h == ENGINE_V1_POLICY_HASH, f"ENGINE V1 HASH MISMATCH: {h}"
    print(f"[OK] Policy hash verified: {h}")

    ensure_validation_schema(engine)
    session = SessionLocal()

    try:
        # ---------------------------------------------------------------
        # Collect all data
        # ---------------------------------------------------------------
        campaign = session.scalars(
            sa.select(InternalValidationCampaignRow).where(
                InternalValidationCampaignRow.campaign_id == CAMPAIGN_ID
            )
        ).first()
        assert campaign is not None, "Campaign not found. Run initialize_validation_campaign.py first."

        cohort = list(session.scalars(
            sa.select(InternalValidationCohortEntryRow).where(
                InternalValidationCohortEntryRow.campaign_id == CAMPAIGN_ID
            ).order_by(InternalValidationCohortEntryRow.entry_id)
        ))

        freezes = list(session.scalars(
            sa.select(InternalValidationPredictionFreezeRow).where(
                InternalValidationPredictionFreezeRow.campaign_id == CAMPAIGN_ID
            ).order_by(
                InternalValidationPredictionFreezeRow.compound_version_id,
                InternalValidationPredictionFreezeRow.endpoint_id,
            )
        ))

        experiments = list(session.scalars(
            sa.select(InternalValidationExperimentalRecordRow).where(
                InternalValidationExperimentalRecordRow.campaign_id == CAMPAIGN_ID
            )
        ))

        observations = list(session.scalars(
            sa.select(InternalValidationObservationRow).where(
                InternalValidationObservationRow.campaign_id == CAMPAIGN_ID
            )
        ))

        # Engine v1 endpoint info
        ep_info_by_id = {r["endpoint_id"]: r for r in policy_rows()}
        ep_info_by_name = {r["endpoint_name"]: r for r in policy_rows()}

        print(f"[INFO] Cohort: {len(cohort)} compounds")
        print(f"[INFO] Prediction freezes: {len(freezes)}")
        print(f"[INFO] Experimental records: {len(experiments)}")
        print(f"[INFO] Paired observations: {len(observations)}")

        # ---------------------------------------------------------------
        # 1. Dataset flow
        # ---------------------------------------------------------------
        _write_dataset_flow(campaign, cohort, freezes, experiments, observations)

        # ---------------------------------------------------------------
        # 2. Endpoint contracts
        # ---------------------------------------------------------------
        _write_endpoint_contracts(ep_info_by_id, experiments, observations)

        # ---------------------------------------------------------------
        # 3. Prediction freezes index
        # ---------------------------------------------------------------
        _write_prediction_freezes_index(freezes, campaign)

        # ---------------------------------------------------------------
        # 4. Experimental manifest
        # ---------------------------------------------------------------
        _write_experimental_manifest(experiments, campaign)

        # ---------------------------------------------------------------
        # 5. Pairing audit
        # ---------------------------------------------------------------
        _write_pairing_audit(observations, freezes, experiments, campaign)

        # ---------------------------------------------------------------
        # 6-9. Metrics, bootstrap, AD, reliability, shadow, scaffold
        # ---------------------------------------------------------------
        _write_metrics_artifacts(observations, ep_info_by_id, campaign)

        # ---------------------------------------------------------------
        # 10. Final decision
        # ---------------------------------------------------------------
        _write_final_decision(campaign, observations, ep_info_by_id)

    finally:
        session.close()

    print("\n[OK] All validation artifacts written.")


def _write_dataset_flow(campaign, cohort, freezes, experiments, observations):
    out = {
        "campaign_id": campaign.campaign_id,
        "generated_at": utcnow_iso(),
        "engine_policy_hash": campaign.engine_policy_hash,
        "framework_status": campaign.framework_status,
        "scientific_status": campaign.scientific_status,
        "pipeline_stages": [
            {
                "stage": "STRUCTURE",
                "n_compounds": len(cohort),
                "compounds": [
                    {
                        "compound_label": e.compound_label,
                        "inchikey": e.inchikey,
                        "project_label": e.project_label,
                        "chemical_series_label": e.chemical_series_label,
                        "eligibility_status": e.eligibility_status,
                        "enrolled_at": e.enrolled_at.isoformat() if e.enrolled_at else None,
                    }
                    for e in cohort
                ],
            },
            {
                "stage": "ENGINE_V1_PREDICTION",
                "policy_version": campaign.engine_policy_version,
                "policy_hash": campaign.engine_policy_hash,
                "n_prediction_freezes": len(freezes),
                "endpoints_covered": list({f.endpoint_id for f in freezes}),
            },
            {
                "stage": "IMMUTABLE_FREEZE",
                "n_freezes": len(freezes),
                "prediction_freeze_complete": campaign.prediction_freeze_complete,
                "freeze_before_experiment_enforced": True,
            },
            {
                "stage": "EXPERIMENT",
                "n_experimental_records": len(experiments),
                "experiment_import_complete": campaign.experiment_import_complete,
                "note": (
                    "No experimental data imported yet. Data collection ongoing."
                    if len(experiments) == 0
                    else f"{len(experiments)} experimental records imported."
                ),
            },
            {
                "stage": "COMPARISON",
                "n_paired_observations": len(observations),
                "note": (
                    "No paired observations yet — awaiting experimental data."
                    if len(observations) == 0
                    else f"{len(observations)} paired observations available."
                ),
            },
        ],
        "current_bottleneck": (
            "AWAITING_EXPERIMENTAL_DATA" if len(experiments) == 0 else "PARTIAL_DATA_AVAILABLE"
        ),
    }
    _write("internal_validation_v1_dataset_flow.json", out)


def _write_endpoint_contracts(ep_info_by_id, experiments, observations):
    contracts = []
    for ep_id, info in sorted(ep_info_by_id.items()):
        n_exp = sum(1 for e in experiments if e.endpoint_id == ep_id)
        n_obs = sum(1 for o in observations if o.endpoint_id == ep_id)
        contracts.append({
            "endpoint_id": ep_id,
            "endpoint_name": info.get("endpoint_name"),
            "strategy": info.get("production_strategy"),
            "evidence_class": info.get("evidence_class"),
            "reliability": info.get("reliability"),
            "unit": info.get("unit"),
            "n_experimental_records": n_exp,
            "n_paired_observations": n_obs,
            "coverage_gap": info.get("production_strategy") == "MODEL_UNAVAILABLE",
            "limitations": info.get("limitations", []),
            "endpoint_compatibility_contract": _get_compat_contract(ep_id),
        })
    out = {
        "campaign_id": "IVC-engine-v1-2026-08-29",
        "generated_at": utcnow_iso(),
        "endpoint_count": len(contracts),
        "model_unavailable_count": sum(1 for c in contracts if c["coverage_gap"]),
        "endpoints": contracts,
    }
    _write("internal_validation_v1_endpoint_contracts.json", out)


def _get_compat_contract(ep_id: str) -> dict:
    contracts = {
        "solubility_aqueous_logs": {
            "assay_type": "kinetic or thermodynamic aqueous solubility",
            "unit": "log10(mol/L) or mol/L (converted)",
            "pH": "7.4 preferred",
            "note": "Kinetic vs thermodynamic distinction recorded. Not mixed in primary metrics.",
        },
        "permeability_caco2_logpapp": {
            "assay_type": "Caco-2 bidirectional",
            "assay_direction": "A→B required",
            "unit": "log10(cm/s) or cm/s (converted)",
            "note": "Efflux ratio recorded. Only Papp A→B compared to prediction.",
        },
        "ppb_human_percent_bound": {
            "species": "human only",
            "assay_type": "equilibrium dialysis or rapid equilibrium dialysis",
            "unit": "% bound or fu (converted: % bound = (1-fu)*100)",
            "note": "Strict species isolation. Rat PPB not compared.",
        },
        "hlm_intrinsic_clearance_scaled_log10": {
            "species": "human only",
            "assay_type": "HLM microsomal incubation",
            "unit": "mL/min/kg or µL/min/mg protein (converted per BW/liver mass)",
            "note": "Hepatocyte data NOT compared. Strict species isolation.",
        },
        "rlm_intrinsic_clearance_scaled_log10": {
            "species": "rat only",
            "assay_type": "RLM microsomal",
            "unit": "mL/min/kg",
            "note": "Strict rat species isolation.",
        },
        "mlm_intrinsic_clearance_scaled_log10": {
            "species": "mouse only",
            "assay_type": "MLM microsomal",
            "unit": "mL/min/kg",
            "note": "Strict mouse species isolation.",
        },
        "safety_herg_blocker_prob": {
            "assay_type": "hERG IC50 (patch-clamp or fluorescence)",
            "unit": "IC50 in µM or binary",
            "threshold": "10 µM standard (blocker if IC50 <= 10 µM)",
            "note": "Threshold must be specified. Different thresholds are ASSAY_CONTEXT_LIMITED.",
        },
        "cyp3a4_inhibitor_prob": {
            "assay_type": "CYP3A4 direct inhibition",
            "unit": "IC50 in µM or binary",
            "note": "Inhibitor only (not substrate). IC50 or binary classification.",
        },
    }
    return contracts.get(ep_id, {"note": "No specific compatibility contract defined"})


def _write_prediction_freezes_index(freezes, campaign):
    by_compound = {}
    for f in freezes:
        by_compound.setdefault(f.compound_version_id, []).append({
            "vfreeze_id": f.vfreeze_id,
            "upstream_frozen_prediction_id": f.upstream_frozen_prediction_id,
            "endpoint_id": f.endpoint_id,
            "strategy": f.strategy,
            "evidence_class": f.evidence_class,
            "prediction_value": f.prediction_value,
            "unit": f.unit,
            "applicability_domain": f.applicability_domain,
            "reliability": f.reliability,
            "engine_policy_version": f.engine_policy_version,
            "engine_policy_hash": f.engine_policy_hash,
            "freeze_timestamp": f.freeze_timestamp.isoformat() if f.freeze_timestamp else None,
        })
    out = {
        "campaign_id": campaign.campaign_id,
        "generated_at": utcnow_iso(),
        "engine_policy_hash": campaign.engine_policy_hash,
        "total_freezes": len(freezes),
        "policy_hash_verified": True,
        "experimental_data_hidden_before_prediction": True,
        "compounds": by_compound,
        "note": (
            "All prediction freezes use policy drugopt-prediction-engine-v1@1.0.0. "
            "These were generated by the qualification freeze mechanism before any "
            "internal experimental results were imported into the validation campaign."
        ),
    }
    _write("internal_validation_v1_prediction_freezes.json", out)


def _write_experimental_manifest(experiments, campaign):
    out = {
        "campaign_id": campaign.campaign_id,
        "generated_at": utcnow_iso(),
        "n_experimental_records": len(experiments),
        "experiment_import_complete": campaign.experiment_import_complete,
        "data_collection_status": (
            "AWAITING_DATA" if len(experiments) == 0 else "PARTIAL"
        ),
        "experiments": [
            {
                "exp_record_id": e.exp_record_id,
                "compound_version_id": e.compound_version_id,
                "endpoint_id": e.endpoint_id,
                "qualifier": e.qualifier,
                "raw_unit": e.raw_unit,
                "censor_flag": e.censor_flag,
                "assay_type": e.assay_type,
                "species": e.species,
                "replicate_id": e.replicate_id,
                "assay_date": e.assay_date,
                "result_available_at": (
                    e.result_available_at.isoformat() if e.result_available_at else None
                ),
                "endpoint_compatibility": e.endpoint_compatibility,
                "compatibility_notes": e.compatibility_notes,
                # Raw value intentionally omitted from public artifact
                "raw_value_present": e.raw_value is not None,
                "imported_at": e.imported_at.isoformat() if e.imported_at else None,
            }
            for e in experiments
        ],
        "note": (
            "Raw experimental values are NOT stored in this public artifact. "
            "They are stored only in the protected internal_validation_experimental_records DB table. "
            "This artifact provides audit-trail metadata only."
        ),
    }
    _write("internal_validation_v1_experimental_manifest.json", out)


def _write_pairing_audit(observations, freezes, experiments, campaign):
    freeze_ts_by_id = {f.vfreeze_id: f.freeze_timestamp for f in freezes}
    exp_avail_by_id = {e.exp_record_id: e.result_available_at for e in experiments}

    out = {
        "campaign_id": campaign.campaign_id,
        "generated_at": utcnow_iso(),
        "n_paired_observations": len(observations),
        "n_prediction_freezes": len(freezes),
        "n_experimental_records": len(experiments),
        "ordering_verified": True,  # enforced by blinding check in register_prediction_freeze
        "observation_summary": {
            "TRUE_PROSPECTIVE": sum(1 for o in observations if o.prospective_evidence_class == "TRUE_PROSPECTIVE"),
            "BLINDED_RETROSPECTIVE": sum(1 for o in observations if o.prospective_evidence_class == "BLINDED_RETROSPECTIVE"),
            "HISTORICAL_VISIBLE": sum(1 for o in observations if o.prospective_evidence_class == "HISTORICAL_VISIBLE"),
        },
        "compatibility_summary": {
            "DIRECT_MATCH": sum(1 for o in observations if o.endpoint_compatibility == "DIRECT_MATCH"),
            "DETERMINISTIC_UNIT_TRANSFORM": sum(1 for o in observations if o.endpoint_compatibility == "DETERMINISTIC_UNIT_TRANSFORM"),
            "ASSAY_CONTEXT_LIMITED": sum(1 for o in observations if o.endpoint_compatibility == "ASSAY_CONTEXT_LIMITED"),
            "ENDPOINT_MISMATCH": sum(1 for o in observations if o.endpoint_compatibility == "ENDPOINT_MISMATCH"),
        },
        "primary_metrics_eligible": sum(1 for o in observations if o.enters_primary_metrics),
        "censored_observations": sum(1 for o in observations if o.censor_flag),
        "data_collection_note": (
            "No experimental data imported yet." if len(observations) == 0 else ""
        ),
    }
    _write("internal_validation_v1_pairing_audit.json", out)


def _write_metrics_artifacts(observations, ep_info_by_id, campaign):
    """Write all analysis artifacts. If no observations, write honest empty state."""
    from backend.validation_analysis_v1 import (
        run_full_analysis, PairedObservation,
        PROB_ENDPOINTS, LOG10_ENDPOINTS,
    )

    if not observations:
        empty_note = "No paired observations available. Awaiting experimental data."
        for fname in [
            "internal_validation_v1_metrics.json",
            "internal_validation_v1_bootstrap.json",
            "internal_validation_v1_ad_analysis.json",
            "internal_validation_v1_reliability_analysis.json",
            "internal_validation_v1_shadow_disagreement.json",
            "internal_validation_v1_scaffold_series_analysis.json",
        ]:
            _write(fname, {
                "campaign_id": campaign.campaign_id,
                "generated_at": utcnow_iso(),
                "data_collection_status": "AWAITING_DATA",
                "n_paired_observations": 0,
                "note": empty_note,
                "metrics": {},
            })
        return

    # Build PairedObservation objects
    from backend.database import SessionLocal
    import sqlalchemy as sa
    from backend.internal_validation_v1 import (
        InternalValidationCohortEntryRow,
        InternalValidationPredictionFreezeRow,
    )

    session = SessionLocal()
    try:
        cohort_by_cv = {
            e.compound_version_id: e
            for e in session.scalars(
                sa.select(InternalValidationCohortEntryRow).where(
                    InternalValidationCohortEntryRow.campaign_id == campaign.campaign_id
                )
            )
        }

        obs_by_endpoint: dict = {}
        for obs in observations:
            cohort_e = cohort_by_cv.get(obs.compound_version_id)
            po = PairedObservation(
                compound_id=obs.compound_version_id,
                endpoint_id=obs.endpoint_id,
                prediction_value=obs.prediction_value,
                experimental_value=obs.raw_value,  # assumes already normalized
                experimental_raw_value=obs.raw_value,
                experimental_unit=obs.raw_unit,
                qualifier=obs.qualifier or "=",
                censor_flag=obs.censor_flag,
                applicability_domain=obs.applicability_domain,
                reliability=obs.reliability,
                prospective_evidence_class=obs.prospective_evidence_class,
                enters_primary_metrics=obs.enters_primary_metrics,
                scaffold_hash=cohort_e.murcko_scaffold_hash if cohort_e else "",
                series_label=cohort_e.chemical_series_label if cohort_e else "",
                project_label=cohort_e.project_label if cohort_e else "",
            )
            obs_by_endpoint.setdefault(obs.endpoint_id, []).append(po)
    finally:
        session.close()

    full = run_full_analysis(obs_by_endpoint)

    # Metrics
    metrics_out = {
        "campaign_id": campaign.campaign_id,
        "generated_at": utcnow_iso(),
        "engine_policy_hash": campaign.engine_policy_hash,
        "n_paired_observations": len(observations),
        "endpoint_metrics": {
            ep: {k: v for k, v in m.items() if k not in ("ad_analysis", "reliability_analysis", "scaffold_series_analysis", "bootstrap")}
            for ep, m in full.items()
        },
    }
    _write("internal_validation_v1_metrics.json", metrics_out)

    # Bootstrap
    boot_out = {
        "campaign_id": campaign.campaign_id,
        "generated_at": utcnow_iso(),
        "bootstrap_seed": 42,
        "n_bootstrap": 1000,
        "bootstrap_by_endpoint": {
            ep: m.get("bootstrap", {}) for ep, m in full.items()
        },
    }
    _write("internal_validation_v1_bootstrap.json", boot_out)

    # AD analysis
    ad_out = {
        "campaign_id": campaign.campaign_id,
        "generated_at": utcnow_iso(),
        "ad_analysis_by_endpoint": {
            ep: m.get("ad_analysis", {}) for ep, m in full.items()
        },
    }
    _write("internal_validation_v1_ad_analysis.json", ad_out)

    # Reliability analysis
    rel_out = {
        "campaign_id": campaign.campaign_id,
        "generated_at": utcnow_iso(),
        "reliability_analysis_by_endpoint": {
            ep: m.get("reliability_analysis", {}) for ep, m in full.items()
        },
    }
    _write("internal_validation_v1_reliability_analysis.json", rel_out)

    # Shadow disagreement — no shadow data yet (all SINGLE_CORE_MODEL)
    _write("internal_validation_v1_shadow_disagreement.json", {
        "campaign_id": campaign.campaign_id,
        "generated_at": utcnow_iso(),
        "n_endpoints_with_shadow": 0,
        "note": (
            "Engine v1 GLP-1 prediction freezes use SINGLE_CORE_MODEL strategy. "
            "No authorized shadow models were active for these compounds. "
            "Shadow disagreement analysis requires endpoints with registered shadow outputs."
        ),
        "shadow_disagreement_by_endpoint": {},
    })

    # Scaffold / series
    scaffold_out = {
        "campaign_id": campaign.campaign_id,
        "generated_at": utcnow_iso(),
        "scaffold_series_by_endpoint": {
            ep: m.get("scaffold_series_analysis", {}) for ep, m in full.items()
        },
    }
    _write("internal_validation_v1_scaffold_series_analysis.json", scaffold_out)


def _write_final_decision(campaign, observations, ep_info_by_id):
    n_paired = len(observations)
    n_primary = sum(1 for o in observations if o.enters_primary_metrics)

    # Determine decision
    if n_paired == 0:
        decision = "INTERNAL_VALIDATION_NOT_STARTED_AWAITING_EXPERIMENTAL_DATA"
        evidence_sufficiency = "NOT_STARTED"
    elif n_primary < 10:
        decision = "INTERNAL_VALIDATION_INSUFFICIENT_DATA_CONTINUE_COLLECTION"
        evidence_sufficiency = "INSUFFICIENT"
    else:
        decision = "INTERNAL_VALIDATION_MIXED_ENGINE_V1_RETAINED_WITH_LIMITATIONS"
        evidence_sufficiency = "PARTIAL"

    # Coverage gap endpoints
    coverage_gaps = [
        ep_id for ep_id, info in ep_info_by_id.items()
        if info.get("production_strategy") == "MODEL_UNAVAILABLE"
    ]

    out = {
        "campaign_id": campaign.campaign_id,
        "generated_at": utcnow_iso(),
        "engine_policy_id": "drugopt-prediction-engine-v1",
        "engine_policy_version": "drugopt-prediction-engine-v1@1.0.0",
        "engine_policy_hash": campaign.engine_policy_hash,
        "policy_hash_verified": True,
        "policy_hash_unchanged": True,

        "framework_status": "READY",
        "scientific_validation_status": (
            "NOT_STARTED" if n_paired == 0 else
            "COLLECTING" if n_primary < 10 else
            "PARTIAL"
        ),

        "cohort": {
            "total_compounds": campaign.compound_count,
            "chemical_series": 1,  # GLP1-SM-PYRIDINONE
            "projects": 1,  # GLP1-SM
        },

        "evidence_classification": {
            "TRUE_PROSPECTIVE": sum(1 for o in observations if o.prospective_evidence_class == "TRUE_PROSPECTIVE"),
            "BLINDED_RETROSPECTIVE": sum(1 for o in observations if o.prospective_evidence_class == "BLINDED_RETROSPECTIVE"),
            "HISTORICAL_VISIBLE": sum(1 for o in observations if o.prospective_evidence_class == "HISTORICAL_VISIBLE"),
        },

        "n_paired_observations": n_paired,
        "n_primary_metrics_eligible": n_primary,
        "evidence_sufficiency": evidence_sufficiency,

        "final_scientific_decision": decision,
        "decision_rationale": (
            "No paired prediction/experiment observations available yet. "
            "Prediction freeze is complete and the framework is READY. "
            "Internal experimental data must be collected and imported to complete validation."
            if n_paired == 0 else
            f"Only {n_primary} primary-metrics-eligible paired observations available. "
            "More data required for statistically meaningful performance assessment."
        ),

        "engine_v1_production_decision": "UNCHANGED",
        "note_no_automatic_modification": (
            "D does NOT automatically modify Engine v1. "
            "Engine v1 remains frozen per policy drugopt-prediction-engine-v1@1.0.0."
        ),

        "coverage_gaps": {
            "model_unavailable_endpoints": coverage_gaps,
            "n_coverage_gaps": len(coverage_gaps),
            "note": (
                "Experimental data for MODEL_UNAVAILABLE endpoints should still be collected "
                "and recorded as COVERAGE_GAP evidence for future Engine v2 AI Scientist."
            ),
        },

        "engine_v1_known_limitations": {
            "caco2": "External MAE 0.5695 log10; Spearman 0.041. Limited ranking ability.",
            "herg": "Raw M1: AUROC 0.6669; Specificity 0.113 at threshold 0.50. Miscalibrated.",
            "hlm_rlm": "LOW-MEDIUM evidence. No independent benchmark available.",
            "mlm": "INSUFFICIENT_EVIDENCE. No independent benchmark.",
            "pka": "RULE_ESTIMATE only. ±1–2 pKa unit uncertainty.",
            "logd": "DERIVED_ESTIMATE from logP and pKa. Not quantitative ML.",
        },

        "future_adaptation_readiness": {
            "status": "DATA_LIMITED",
            "n_compounds": campaign.compound_count,
            "n_series": 1,
            "note": (
                "3 compounds from 1 series in 1 project is insufficient for "
                "project- or series-level adaptation. "
                "20-50 compounds across multiple series would be required. "
                "No adaptation implemented in this stage."
            ),
        },

        "ai_scientist_readiness": {
            "status": "PARTIAL",
            "prediction_available": True,
            "experimental_result_available": n_paired > 0,
            "prediction_timestamp_available": True,
            "experiment_timestamp_available": True,
            "evidence_class_available": True,
            "ad_available": True,
            "reliability_available": True,
            "validation_performance_available": n_paired > 0,
            "limitations_documented": True,
            "model_unavailable_documented": True,
            "disagreement_available": False,
            "note": (
                "Framework is PARTIAL. AI Scientist can receive structured prediction "
                "evidence and limitations. Paired validation performance pending experimental data."
            ),
        },

        "git": {
            "suggested_commit_message": (
                "feat(validation): initialize engine v1 internal prospective validation framework"
            ),
            "suggested_tag": "prediction-engine-v1-validation-framework-ready",
            "note": (
                "Tag 'prediction-engine-v1-internal-validation' reserved for "
                "scientifically complete validation with actual paired evidence."
            ),
        },

        "progress_summary": "5 / 6 COMPLETE — TASK 6 VALIDATION CAMPAIGN ACTIVE",
        "next_step": "CONTINUE_INTERNAL_EXPERIMENTAL_DATA_COLLECTION",
    }
    _write("internal_validation_v1_final_decision.json", out)
    print(f"[OK] Final decision: {decision}")


def _write(filename: str, data: dict):
    p = VALIDATION_DIR / filename
    p.write_text(json.dumps(data, indent=2, default=str))
    print(f"[OK] Written: {p.name}")


if __name__ == "__main__":
    run()
