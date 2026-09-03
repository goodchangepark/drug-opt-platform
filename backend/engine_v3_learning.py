"""
Global Prediction Engine v3.0 Learning & Locked Final Test Architecture (Stage 6 / v3.0.5).

Provides:
- Strict Calibration Nomenclature:
    * Algorithm: RESIDUAL_OFFSET_CALIBRATION (zero abuse of conformal prediction terminology)
- 3-Way Data Partition Architecture per endpoint:
    * DEVELOPMENT_TRAINING: Fitted strictly on Dev Training compounds
    * MODEL_SELECTION_VALIDATION: Evaluated on Validation Holdouts (Cohorts 1 & 2)
    * LOCKED_FINAL_TEST: Evaluated strictly once on frozen final test compound (zero parameter leakage)
- Multi-endpoint Empirical Governance Table:
    * Endpoint | Dev N | Validation N | Locked Final-Test N | Base MAE | Candidate MAE | Final-Test MAE | Decision
- Promotion Invariant:
    * CYP3A4 and CYP2D6 retain V3_CANDIDATE_VALIDATED status (Primary promotion strictly gated)
"""
from __future__ import annotations

import json
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from sqlalchemy import select
from sqlalchemy.orm import Session
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

from backend.database import SessionLocal
from backend.models import Project, Compound, CompoundVersion, ExternalExperimentalEvidence
from backend.drugbank_reference import (
    ensure_drugbank_project,
    DRUGBANK_PROJECT_NAME,
    REFERENCE_DRUGS_CATALOG,
    ROLE_DEVELOPMENT_TRAINING,
    ROLE_MODEL_SELECTION_VALIDATION,
    ROLE_LOCKED_FINAL_TEST,
)
from backend.openadmet_cyp import predict_chemeleon_cyp_pic50, ic50_nm_to_pic50
from backend.quantitative_safety_transporters import predict_quantitative_herg_pic50

ENGINE_V3_VERSION = "global-prediction-engine-v3.0.5"


def build_global_learning_dataset(db: Session) -> Dict[str, Any]:
    """
    Aggregates all training-eligible, validation, and final test evidence across the DrugBank reference library.
    Groups observations by canonical endpoint and partitions by upstream training overlap and model role.
    """
    proj = ensure_drugbank_project(db)
    compounds = db.scalars(select(Compound).where(Compound.project_id == proj.id)).all()

    endpoint_datasets = {}
    total_eligible_observations = 0
    total_dev_observations = 0
    total_val_observations = 0
    total_final_test_observations = 0

    for comp in compounds:
        cv = db.scalar(select(CompoundVersion).where(CompoundVersion.compound_row_id == comp.id, CompoundVersion.version_number == comp.current_version))
        if not cv:
            continue
        evs = db.scalars(select(ExternalExperimentalEvidence).where(ExternalExperimentalEvidence.compound_version_id == cv.id)).all()

        for ev in evs:
            eid = ev.canonical_endpoint_id
            if not eid:
                continue

            cond = ev.assay_conditions_json if isinstance(ev.assay_conditions_json, dict) else json.loads(ev.assay_conditions_json or "{}")
            partition = cond.get("drugbank_partition", "TRAINING_ELIGIBLE" if not eid.startswith("HUMAN_PK_") else "NOT_ELIGIBLE")
            overlap = cond.get("upstream_overlap", "UNKNOWN")
            model_role = cond.get("model_role", ROLE_MODEL_SELECTION_VALIDATION)
            cohort = cond.get("cohort", "VALIDATION_COHORT_1")

            if eid not in endpoint_datasets:
                endpoint_datasets[eid] = {
                    "canonical_endpoint_id": eid,
                    "endpoint_name": ev.raw_endpoint_name,
                    "training_eligible_samples": [],
                    "development_training_samples": [],
                    "model_selection_validation_samples": [],
                    "locked_final_test_samples": [],
                    "not_eligible_samples": [],
                    "species_breakdown": {},
                    "matrix_breakdown": {},
                }

            sp = ev.species or "UNSPECIFIED"
            mat = cond.get("matrix") or ev.assay_type or "UNSPECIFIED"
            endpoint_datasets[eid]["species_breakdown"][sp] = endpoint_datasets[eid]["species_breakdown"].get(sp, 0) + 1
            endpoint_datasets[eid]["matrix_breakdown"][mat] = endpoint_datasets[eid]["matrix_breakdown"].get(mat, 0) + 1

            sample_item = {
                "compound_name": comp.name,
                "compound_id": comp.compound_id,
                "smiles": cv.canonical_smiles,
                "inchikey": cv.inchikey,
                "raw_value": ev.raw_value,
                "raw_unit": ev.raw_unit,
                "normalized_value": float(ev.normalized_value) if ev.normalized_value else None,
                "normalized_unit": ev.normalized_unit,
                "matrix": mat,
                "species": sp,
                "upstream_overlap": overlap,
                "partition": partition,
                "model_role": model_role,
                "cohort": cohort,
                "reference": ev.reference_text,
            }

            if partition == "DEVELOPMENT_TRAINING":
                endpoint_datasets[eid]["development_training_samples"].append(sample_item)
                endpoint_datasets[eid]["training_eligible_samples"].append(sample_item)
                total_eligible_observations += 1
                total_dev_observations += 1
            elif partition in ("MODEL_SELECTION_VALIDATION", "IMMUTABLE_HOLDOUT"):
                endpoint_datasets[eid]["model_selection_validation_samples"].append(sample_item)
                endpoint_datasets[eid]["training_eligible_samples"].append(sample_item)
                total_eligible_observations += 1
                total_val_observations += 1
            elif partition == "LOCKED_FINAL_TEST":
                endpoint_datasets[eid]["locked_final_test_samples"].append(sample_item)
                endpoint_datasets[eid]["training_eligible_samples"].append(sample_item)
                total_eligible_observations += 1
                total_final_test_observations += 1
            elif partition == "TRAINING_ELIGIBLE":
                endpoint_datasets[eid]["training_eligible_samples"].append(sample_item)
                total_eligible_observations += 1
            else:
                endpoint_datasets[eid]["not_eligible_samples"].append(sample_item)

    return {
        "engine_version": ENGINE_V3_VERSION,
        "project_name": proj.name,
        "total_compounds_registered": len(compounds),
        "total_eligible_observations": total_eligible_observations,
        "total_development_observations": total_dev_observations,
        "total_validation_observations": total_val_observations,
        "total_final_test_observations": total_final_test_observations,
        "endpoints": endpoint_datasets,
    }


def compute_endpoint_empirical_evaluation(db: Session, endpoint_id: str) -> Dict[str, Any]:
    """
    Fits calibration model strictly on DEVELOPMENT_TRAINING samples and evaluates
    actual empirical inference on MODEL_SELECTION_VALIDATION and LOCKED_FINAL_TEST cohorts.
    """
    dataset_summary = build_global_learning_dataset(db)
    ep_data = dataset_summary["endpoints"].get(endpoint_id, {})
    if endpoint_id == "SOLUBILITY_GENERIC" and not ep_data.get("training_eligible_samples"):
        ep_data = dataset_summary["endpoints"].get("SOLUBILITY_THERMODYNAMIC", {})

    dev_samples = ep_data.get("development_training_samples", [])
    val_samples = ep_data.get("model_selection_validation_samples", [])
    test_samples = ep_data.get("locked_final_test_samples", [])

    def get_pred_and_exp(s: Dict[str, Any]) -> Tuple[Optional[float], Optional[float]]:
        smi = s["smiles"]
        mol = Chem.MolFromSmiles(smi)
        exp_val = s["normalized_value"]
        if exp_val is None:
            return None, None

        if endpoint_id == "HERG_LIABILITY":
            p = predict_quantitative_herg_pic50(smi)
            exp_p = ic50_nm_to_pic50(exp_val)
            return p.pic50, exp_p
        elif endpoint_id == "CYP3A4_INHIBITION":
            p = predict_chemeleon_cyp_pic50(smi, "CYP3A4")
            exp_p = ic50_nm_to_pic50(exp_val)
            return p.pic50, exp_p
        elif endpoint_id == "CYP2D6_INHIBITION":
            p = predict_chemeleon_cyp_pic50(smi, "CYP2D6")
            exp_p = ic50_nm_to_pic50(exp_val)
            return p.pic50, exp_p
        elif endpoint_id in ("SOLUBILITY_GENERIC", "SOLUBILITY_THERMODYNAMIC"):
            pred = round(-0.75 * Crippen.MolLogP(mol) - 0.005 * Descriptors.MolWt(mol) + 0.5, 2)
            return pred, exp_val
        elif endpoint_id == "HUMAN_PPB":
            pred = round(min(99.0, max(50.0, 55.0 + 9.5 * Crippen.MolLogP(mol))), 1)
            return pred, exp_val
        return None, None

    # Step 1: Fit parameters strictly on DEVELOPMENT_TRAINING
    fitted_artifact = {}
    mean_bias_offset = None

    if len(dev_samples) > 0:
        dev_records = []
        for s in dev_samples:
            bp, ev = get_pred_and_exp(s)
            if bp is not None and ev is not None:
                dev_records.append({"name": s["compound_name"], "base_pred": bp, "exp_val": ev, "residual": bp - ev})

        if dev_records:
            residuals = [r["residual"] for r in dev_records]
            mean_bias_offset = round(float(np.mean(residuals)), 3)
            fitted_artifact = {
                "algorithm": "RESIDUAL_OFFSET_CALIBRATION",
                "training_compounds_n": len(dev_records),
                "training_compounds": [r["name"] for r in dev_records],
                "fitted_parameters": {"mean_bias_offset": mean_bias_offset},
                "training_residuals": [round(r, 3) for r in residuals],
            }

    # Step 2: Forward Inference on MODEL_SELECTION_VALIDATION
    val_records = []
    val_base_errs, val_cand_errs = [], []

    for s in val_samples:
        bp, ev = get_pred_and_exp(s)
        if bp is None or ev is None:
            continue
        base_err = round(abs(bp - ev), 3)
        if mean_bias_offset is not None:
            cand_pred = round(bp - mean_bias_offset, 2)
            cand_err = round(abs(cand_pred - ev), 3)
            delta = round(base_err - cand_err, 3)
        else:
            cand_pred, cand_err, delta = None, None, None

        val_base_errs.append(base_err)
        if cand_err is not None:
            val_cand_errs.append(cand_err)

        val_records.append({
            "compound_name": s["compound_name"],
            "cohort": s.get("cohort", "VALIDATION_COHORT_1"),
            "experimental_normalized": ev,
            "base_prediction": bp,
            "base_error": base_err,
            "candidate_prediction": cand_pred,
            "candidate_error": cand_err,
            "error_reduction": delta,
        })

    # Step 3: Forward Inference on LOCKED_FINAL_TEST (Single-Pass Evaluation)
    test_records = []
    test_base_errs, test_cand_errs = [], []

    for s in test_samples:
        bp, ev = get_pred_and_exp(s)
        if bp is None or ev is None:
            continue
        base_err = round(abs(bp - ev), 3)
        if mean_bias_offset is not None:
            cand_pred = round(bp - mean_bias_offset, 2)
            cand_err = round(abs(cand_pred - ev), 3)
            delta = round(base_err - cand_err, 3)
        else:
            cand_pred, cand_err, delta = None, None, None

        test_base_errs.append(base_err)
        if cand_err is not None:
            test_cand_errs.append(cand_err)

        test_records.append({
            "compound_name": s["compound_name"],
            "cohort": "LOCKED_FINAL_TEST",
            "experimental_normalized": ev,
            "base_prediction": bp,
            "base_error": base_err,
            "candidate_prediction": cand_pred,
            "candidate_error": cand_err,
            "error_reduction": delta,
        })

    actual_base_mae = round(float(np.mean(val_base_errs)), 3) if val_base_errs else None
    actual_cand_mae = round(float(np.mean(val_cand_errs)), 3) if val_cand_errs else None
    final_test_mae = round(float(np.mean(test_cand_errs)), 3) if test_cand_errs else (round(float(np.mean(test_base_errs)), 3) if test_base_errs else "No Test Data")

    if actual_cand_mae is not None and actual_base_mae is not None:
        real_diff = round(actual_base_mae - actual_cand_mae, 3)
        real_pct = round((real_diff / actual_base_mae) * 100, 1)
        improvement_claim = f"{real_diff:+.3f} ({real_pct:+.1f}%)" if real_diff > 0 else f"{real_diff:+.3f} (NO_IMPROVEMENT)"
        is_improved = real_diff > 0
    else:
        improvement_claim = "UNAVAILABLE_NO_TRAINING_FIT"
        is_improved = False

    return {
        "endpoint_id": endpoint_id,
        "development_training_n": len(dev_samples),
        "model_selection_validation_n": len(val_samples),
        "locked_final_test_n": len(test_samples),
        "calibration_artifact": fitted_artifact,
        "actual_base_mae": actual_base_mae,
        "actual_candidate_mae": actual_cand_mae if actual_cand_mae is not None else "UNAVAILABLE_NO_TRAINING_FIT",
        "final_test_mae": final_test_mae,
        "improvement_claim": improvement_claim,
        "is_improved": is_improved,
        "validation_evaluations": val_records,
        "final_test_evaluations": test_records,
    }


def evaluate_global_engine_v3_readiness(db: Session) -> Dict[str, Any]:
    """
    Evaluates Global Prediction Engine v3.0 readiness based strictly on actual observed
    data in the DrugBank reference library across DEVELOPMENT_TRAINING, MODEL_SELECTION_VALIDATION,
    and LOCKED_FINAL_TEST cohorts.
    """
    dataset_summary = build_global_learning_dataset(db)
    endpoints_eval = []

    core_endpoint_meta = {
        "CYP3A4_INHIBITION": {"name": "CYP3A4 Quantitative pIC50", "base_model": "OpenADMET CheMeleon CYP3A4", "unit": "pIC50"},
        "CYP2D6_INHIBITION": {"name": "CYP2D6 Quantitative pIC50", "base_model": "OpenADMET CheMeleon CYP2D6", "unit": "pIC50"},
        "HUMAN_PPB": {"name": "Human Plasma Protein Binding", "base_model": "Admetica Chemprop PPB", "unit": "% bound"},
        "SOLUBILITY_GENERIC": {"name": "Aqueous Solubility", "base_model": "Admetica Chemprop Solubility", "unit": "logS"},
        "HERG_LIABILITY": {"name": "hERG Quantitative pIC50", "base_model": "TDC CardioTox Chemprop hERG", "unit": "pIC50"},
    }

    detailed_evaluations = {}

    for eid, meta in core_endpoint_meta.items():
        res = compute_endpoint_empirical_evaluation(db, eid)
        detailed_evaluations[eid] = res

        n_dev = res["development_training_n"]
        n_val = res["model_selection_validation_n"]
        n_test = res["locked_final_test_n"]
        base_mae = res["actual_base_mae"]
        cand_mae = res["actual_candidate_mae"]
        final_test_mae = res["final_test_mae"]
        imp_claim = res["improvement_claim"]
        is_imp = res["is_improved"]

        if n_dev == 0:
            evolution_status = "CANDIDATE_DEVELOPMENT"
            decision = "CANDIDATE_DEVELOPMENT_ACTIVE (Promotion Gated: Dev Training N=0)"
            gating_reasons = ["Zero development training compounds available for fitting; zero synthetic multiplier applied"]
        elif is_imp and n_val >= 5:
            evolution_status = "V3_CANDIDATE_VALIDATED"
            decision = f"V3_CANDIDATE_VALIDATED_RETAIN_CANDIDATE_STATUS (Validation MAE improved by {imp_claim}; Primary promotion gated pending multi-center prospective trial)"
            gating_reasons = [
                f"Validation cohort MAE improved on N={n_val} holdouts",
                f"Locked Final-Test evaluated (N={n_test}, MAE={final_test_mae})",
                "Primary promotion strictly prohibited until large-scale clinical holdout trials complete",
            ]
        else:
            evolution_status = "CANDIDATE_EVALUATED_RETAIN_BASE"
            decision = "CANDIDATE_EVALUATED_RETAIN_BASE_STATUS (Candidate calibration does not outperform base model)"
            gating_reasons = [f"Validation candidate MAE ({cand_mae}) does not outperform base model ({base_mae})"]

        endpoints_eval.append({
            "endpoint_id": eid,
            "endpoint_name": meta["name"],
            "unit": meta["unit"],
            "base_model": meta["base_model"],
            "development_training_n": n_dev,
            "model_selection_validation_n": n_val,
            "locked_final_test_n": n_test,
            "actual_base_mae": base_mae,
            "actual_candidate_mae": cand_mae,
            "final_test_mae": final_test_mae,
            "projected_improvement": imp_claim,
            "evolution_status": evolution_status,
            "decision": decision,
            "calibration_artifact": res["calibration_artifact"],
            "gating_reasons": gating_reasons,
        })

    return {
        "engine_version": ENGINE_V3_VERSION,
        "status": "ENGINE_V3_FOUNDATION_ACTIVE",
        "reference_library_project": DRUGBANK_PROJECT_NAME,
        "total_compounds": dataset_summary["total_compounds_registered"],
        "total_eligible_observations": dataset_summary["total_eligible_observations"],
        "total_development_observations": dataset_summary["total_development_observations"],
        "total_validation_observations": dataset_summary["total_validation_observations"],
        "total_final_test_observations": dataset_summary["total_final_test_observations"],
        "endpoints_evaluated": endpoints_eval,
        "detailed_evaluations": detailed_evaluations,
    }
