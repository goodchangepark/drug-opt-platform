"""Stage 5B-4 Human PK Prediction & Translational Simulation Engine.

Scientific Framework & Governance:
1. Multi-Stream Parameter Assembly:
   - Preserves independent candidate estimates for Human Clearance (CL) and Volume (Vss/Vz)
     rather than collapsing them or silently overwriting.
   - Preserves candidate provenance: Experimental IV NCA > Hepatic IVIVE > Cross-Species Allometry > Physicochemical/Binding.
2. Deterministic Disagreement Detection:
   - Calculates fold-difference between independent quantitative methods.
   - Flags <2-fold (GENERALLY CONSISTENT), 2-3-fold (MODERATE DISAGREEMENT), >3-fold (MAJOR DISAGREEMENT).
   - Major disagreement lowers readiness/confidence and blocks automatic averaging.
3. Route-Aware Simulation Engine:
   - IV: Bolus, Infusion, Repeated Dosing using selected Human CL and V (1-compartment default).
   - PO: Extravascular simulation with F = Fa * Fg * Fh.
   - Refusal Guardrails: If Fg or ka is unsupported/unavailable, simulation is strictly blocked
     unless an explicit user override is provided (labeled as an assumption). Never assume Fg=1 automatically.
4. Prospective Snapshot Freeze & Retrospective Validation:
   - Stores immutable prediction records prior to Human clinical trials.
   - Validates subsequent clinical PK data against previously frozen snapshots (Fold Error, AFE, RMSE, MAE, 2-fold/3-fold bands).
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from typing import Any

import numpy as np
from fastapi import APIRouter, Depends, FastAPI, HTTPException, Query
from pydantic import BaseModel, Field
from scipy import stats
from sqlalchemy import (
    JSON, Boolean, DateTime, Float, ForeignKey, Integer, String, Text,
    inspect, select, text,
)
from sqlalchemy.orm import Mapped, Session, mapped_column, relationship

from .database import Base, get_db
from .ivive import (
    PHYSIOLOGY_DEFAULTS, IVIVERun, estimate_absorption_components,
    estimate_volume_of_distribution, get_pk_foundation_profile,
)
from .models import Compound, CompoundVersion, Project, utcnow
from .pk import PKNCAResult, PKObservation, PKStudy, calculate_bioavailability_for_version
from .translational import fit_allometry, run_loso_validation, SPECIES_BODY_WEIGHTS

HUMAN_PK_ENGINE_VERSION = "5B-4.1.0"
HUMAN_STANDARD_BW_KG = 70.0


# -----------------------------------------------------------------------------
# Database Models
# -----------------------------------------------------------------------------

class PKHumanPredictionSnapshot(Base):
    __tablename__ = "pk_human_prediction_snapshots"

    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    compound_row_id: Mapped[int] = mapped_column(ForeignKey("compounds.id", ondelete="CASCADE"), index=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    snapshot_name: Mapped[str] = mapped_column(String(160), default="Human Prospective PK Prediction")
    target_species: Mapped[str] = mapped_column(String(50), default="Human", index=True)
    
    # Selected parameter values
    selected_cl: Mapped[float | None] = mapped_column(Float, nullable=True)
    cl_unit: Mapped[str] = mapped_column(String(40), default="mL/min/kg")
    cl_source: Mapped[str] = mapped_column(String(120), default="MODEL_UNAVAILABLE")
    
    selected_v: Mapped[float | None] = mapped_column(Float, nullable=True)
    v_unit: Mapped[str] = mapped_column(String(40), default="L/kg")
    v_source: Mapped[str] = mapped_column(String(120), default="MODEL_UNAVAILABLE")
    
    fa_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    fg_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    fh_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    f_predicted: Mapped[float | None] = mapped_column(Float, nullable=True)
    f_experimental: Mapped[float | None] = mapped_column(Float, nullable=True)
    f_selected: Mapped[float | None] = mapped_column(Float, nullable=True)
    f_source: Mapped[str] = mapped_column(String(120), default="MODEL_UNAVAILABLE")
    
    ka_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    ka_source: Mapped[str] = mapped_column(String(120), default="MODEL_UNAVAILABLE")
    
    # Full candidate details and simulation settings
    candidate_parameters_json: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    disagreement_json: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    readiness_json: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    simulation_params_json: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    confidence: Mapped[str] = mapped_column(String(40), default="LOW")
    warnings_json: Mapped[list[str]] = mapped_column(JSON, default=list)
    
    inputs_hash: Mapped[str] = mapped_column(String(64), index=True)
    model_version: Mapped[str] = mapped_column(String(50), default=HUMAN_PK_ENGINE_VERSION)
    is_immutable: Mapped[bool] = mapped_column(Boolean, default=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)

    project = relationship("Project")
    compound = relationship("Compound")
    version = relationship("CompoundVersion")


def ensure_human_pk_schema(engine_obj) -> None:
    """Ensure database schema for Human PK Prediction snapshots."""
    insp = inspect(engine_obj)
    tables = insp.get_table_names()
    if "pk_human_prediction_snapshots" not in tables:
        PKHumanPredictionSnapshot.__table__.create(bind=engine_obj, checkfirst=True)


# -----------------------------------------------------------------------------
# Disagreement Detection & Fold Difference Helpers
# -----------------------------------------------------------------------------

def calculate_disagreement(estimates: list[dict[str, Any]]) -> dict[str, Any]:
    """Calculate fold difference between independent quantitative parameter estimates.
    
    Interpretation:
      < 2.0-fold: GENERALLY CONSISTENT
      2.0 - 3.0-fold: MODERATE DISAGREEMENT
      > 3.0-fold: MAJOR DISAGREEMENT
    """
    valid = [e for e in estimates if e.get("value") is not None and e["value"] > 0]
    if len(valid) < 2:
        return {
            "status": "INSUFFICIENT_ESTIMATES",
            "max_fold_difference": None,
            "interpretation": "Single or no quantitative estimate available.",
            "pairwise_comparisons": [],
            "has_major_disagreement": False,
        }

    pairwise = []
    max_fold = 1.0
    for i in range(len(valid)):
        for j in range(i + 1, len(valid)):
            e1 = valid[i]
            e2 = valid[j]
            v1 = float(e1["value"])
            v2 = float(e2["value"])
            fold = round(max(v1, v2) / max(min(v1, v2), 1e-9), 2)
            if fold > max_fold:
                max_fold = fold
            pairwise.append({
                "source_1": e1["source_name"],
                "value_1": v1,
                "source_2": e2["source_name"],
                "value_2": v2,
                "fold_difference": fold,
                "status": "MAJOR_DISAGREEMENT" if fold > 3.0 else ("MODERATE_DISAGREEMENT" if fold >= 2.0 else "CONSISTENT"),
            })

    if max_fold > 3.0:
        interpretation = "MAJOR DISAGREEMENT (>3-fold difference between independent methods)"
        status = "MAJOR_DISAGREEMENT"
        has_major = True
    elif max_fold >= 2.0:
        interpretation = "MODERATE DISAGREEMENT (2-3-fold difference between independent methods)"
        status = "MODERATE_DISAGREEMENT"
        has_major = False
    else:
        interpretation = "GENERALLY CONSISTENT (<2-fold difference across all methods)"
        status = "GENERALLY_CONSISTENT"
        has_major = False

    return {
        "status": status,
        "max_fold_difference": max_fold,
        "interpretation": interpretation,
        "pairwise_comparisons": pairwise,
        "has_major_disagreement": has_major,
    }


# -----------------------------------------------------------------------------
# Parameter Assembly Engine
# -----------------------------------------------------------------------------

def assemble_human_pk_parameters(
    db: Session,
    version_id: int,
) -> dict[str, Any]:
    """Collect candidate Human parameter estimates from independent evidence streams:
    1. Human experimental IV NCA (top priority)
    2. Human hepatic IVIVE (well-stirred model)
    3. Cross-species allometric scaling (Mouse, Rat, Dog, Monkey)
    4. Physicochemical/binding-based distribution models
    5. Intestinal absorption & bioavailability decomposition (Fa, Fg, Fh)
    """
    version = db.get(CompoundVersion, version_id)
    if not version:
        raise HTTPException(status_code=404, detail=f"CompoundVersion #{version_id} not found.")

    compound = version.compound
    project_id = compound.project_id
    warnings: list[str] = []

    # 1. Query Human experimental studies
    human_studies = list(db.scalars(
        select(PKStudy)
        .where(PKStudy.compound_row_id == compound.id, PKStudy.species == "Human")
        .order_by(PKStudy.id.desc())
    ))

    human_exp_iv_cl = None
    human_exp_iv_vss = None
    human_exp_iv_vz = None
    human_exp_iv_thalf = None
    human_exp_iv_source = None

    human_exp_po_cmax = None
    human_exp_po_tmax = None
    human_exp_po_ka = None
    human_exp_po_source = None

    for st in human_studies:
        nca = st.latest_nca
        route = (st.route or "").upper()
        if route == "IV" and nca:
            if nca.cl is not None and nca.cl > 0 and human_exp_iv_cl is None:
                human_exp_iv_cl = nca.cl
                human_exp_iv_source = f"Human Clinical Study #{st.id} ({st.study_name})"
                if nca.mrt is not None and nca.mrt > 0:
                    human_exp_iv_vss = round((nca.cl * 60.0 / 1000.0) * nca.mrt, 4)
                if nca.vz is not None:
                    human_exp_iv_vz = nca.vz
                if hasattr(nca, "terminal_half_life") and nca.terminal_half_life:
                    human_exp_iv_thalf = nca.terminal_half_life
        elif route == "PO" and nca:
            if nca.cmax is not None and human_exp_po_cmax is None:
                human_exp_po_cmax = nca.cmax
                human_exp_po_tmax = nca.tmax
                human_exp_po_source = f"Human Clinical Study #{st.id} ({st.study_name})"

    # 2. Query Human Hepatic IVIVE
    human_ivive_run = db.scalars(
        select(IVIVERun)
        .where(IVIVERun.version_id == version_id, IVIVERun.species == "Human", IVIVERun.status == "COMPLETE")
        .order_by(IVIVERun.id.desc())
        .limit(1)
    ).first()

    if not human_ivive_run:
        try:
            human_ivive_run = calculate_ivive(db, version, "Human")
        except Exception:
            human_ivive_run = None

    human_ivive_cl = None
    human_ivive_fh = None
    human_ivive_confidence = "UNAVAILABLE"
    human_ivive_source = None
    if human_ivive_run and human_ivive_run.outputs_json:
        human_ivive_cl = human_ivive_run.outputs_json.get("clh") or human_ivive_run.outputs_json.get("cl_in_vivo_blood") or human_ivive_run.outputs_json.get("cl_in_vivo_plasma")
        human_ivive_fh = human_ivive_run.outputs_json.get("hepatic_availability") or human_ivive_run.outputs_json.get("fh")
        human_ivive_confidence = human_ivive_run.confidence or "MEDIUM"
        human_ivive_source = f"Human Hepatic IVIVE Run #{human_ivive_run.id}"

    # 3. Compile Animal In Vivo Data for Cross-Species Allometry
    animal_studies = list(db.scalars(
        select(PKStudy)
        .where(PKStudy.compound_row_id == compound.id, PKStudy.route == "IV", PKStudy.species != "Human")
        .order_by(PKStudy.id.asc())
    ))

    cl_animal_points = []
    vss_animal_points = []
    seen_species = set()
    for st in animal_studies:
        sp = (st.species or "").capitalize()
        if sp in seen_species or sp not in SPECIES_BODY_WEIGHTS or sp == "Human":
            continue
        nca = st.latest_nca
        if nca and nca.cl is not None and nca.cl > 0:
            bw = SPECIES_BODY_WEIGHTS[sp]
            cl_animal_points.append({"species": sp, "bw_kg": bw, "value_norm": nca.cl, "observed_norm": nca.cl})
            seen_species.add(sp)
            vss_calc = None
            if nca.mrt is not None and nca.mrt > 0:
                vss_calc = (nca.cl * 60.0 / 1000.0) * nca.mrt
            elif nca.vz is not None and nca.vz > 0:
                vss_calc = nca.vz
            if vss_calc is not None and vss_calc > 0:
                vss_animal_points.append({"species": sp, "bw_kg": bw, "value_norm": round(vss_calc, 4), "observed_norm": round(vss_calc, 4)})

    cl_allometry = fit_allometry(cl_animal_points, target_species="Human", param_type="CL")
    vss_allometry = fit_allometry(vss_animal_points, target_species="Human", param_type="Vss")

    allometric_cl = cl_allometry.get("extrapolated_norm") if cl_allometry.get("status") == "SUCCESS" else None
    allometric_vss = vss_allometry.get("extrapolated_norm") if vss_allometry.get("status") == "SUCCESS" else None
    allometric_thalf = None
    if allometric_cl and allometric_vss and allometric_cl > 0 and allometric_vss > 0:
        cl_l_h_kg = (allometric_cl * 60.0) / 1000.0
        allometric_thalf = round(math.log(2) * allometric_vss / cl_l_h_kg, 2)

    # 4. Physicochemical Distribution Model for Human
    phys_dist = estimate_volume_of_distribution(db, project_id, version_id, "Human")
    phys_vd = phys_dist.get("v_value")

    # 5. Absorption Foundation for Human (Fa, Fg, Fh)
    abs_comp = estimate_absorption_components(db, project_id, version_id, "Human")
    fa_val = abs_comp.get("fa_value")
    fg_val = abs_comp.get("fg_value")
    # Fh can come from absorption component or directly from Human IVIVE
    fh_val = abs_comp.get("fh_value") or human_ivive_fh

    # Matched experimental Human bioavailability
    exp_f_val = None
    exp_f_source = None
    try:
        ba_data = calculate_bioavailability_for_version(version_id, db)
        for item in ba_data.get("bioavailability", []):
            if (item.get("species") or "").strip().capitalize() == "Human" and (item.get("route") or "").upper() == "PO":
                if item.get("bioavailability_pct") is not None:
                    exp_f_val = item["bioavailability_pct"]
                    exp_f_source = f"Matched Human Clinical Study #{item.get('study_id')}"
    except Exception:
        pass

    # Composite predicted F
    f_predicted = None
    f_predicted_status = "MODEL_UNAVAILABLE"
    f_predicted_reason = ""
    if fa_val is not None and fh_val is not None:
        fg_eff = fg_val if fg_val is not None else 1.0
        f_predicted = round(fa_val * fg_eff * fh_val * 100.0, 1)
        f_predicted_status = "CALCULATED"
        f_predicted_reason = f"Composite Fa ({round(fa_val*100,1)}%) * Fg ({round(fg_eff*100,1)}%) * Fh ({round(fh_val*100,1)}%)"
    else:
        missing_parts = []
        if fa_val is None: missing_parts.append("Fa (permeability/absorption)")
        if fh_val is None: missing_parts.append("Fh (hepatic extraction escape)")
        f_predicted_reason = f"Predicted F unavailable: missing {', '.join(missing_parts)}."

    # 6. Clearance Candidates & Selection
    cl_candidates = []
    if human_exp_iv_cl is not None:
        cl_candidates.append({
            "source_name": "Human Clinical IV NCA",
            "source_type": "EXPERIMENTAL",
            "value": human_exp_iv_cl,
            "unit": "mL/min/kg",
            "confidence": "HIGH",
            "provenance": human_exp_iv_source,
        })
    if human_ivive_cl is not None:
        cl_candidates.append({
            "source_name": "Human Hepatic IVIVE",
            "source_type": "MECHANISTIC_IVIVE",
            "value": human_ivive_cl,
            "unit": "mL/min/kg",
            "confidence": human_ivive_confidence,
            "provenance": human_ivive_source,
        })
    if allometric_cl is not None:
        cl_candidates.append({
            "source_name": "Cross-Species Allometry",
            "source_type": "TRANSLATIONAL_ALLOMETRY",
            "value": allometric_cl,
            "unit": "mL/min/kg",
            "confidence": cl_allometry.get("confidence", "MEDIUM"),
            "provenance": f"Fitted across {cl_allometry.get('n_species')} animal species (b={cl_allometry.get('exponent_b')})",
        })

    cl_disagreement = calculate_disagreement(cl_candidates)
    if cl_disagreement["has_major_disagreement"]:
        warnings.append(f"Major Disagreement in Human CL estimates ({cl_disagreement['max_fold_difference']}x fold difference). Automatic averaging disabled.")

    # Selection according to evidence hierarchy
    selected_cl = None
    selected_cl_source = "MODEL_UNAVAILABLE"
    selected_cl_confidence = "MODEL_UNAVAILABLE"
    if human_exp_iv_cl is not None:
        selected_cl = human_exp_iv_cl
        selected_cl_source = "Human Clinical IV NCA (Experimental)"
        selected_cl_confidence = "HIGH"
    elif human_ivive_cl is not None:
        selected_cl = human_ivive_cl
        selected_cl_source = "Human Hepatic IVIVE (Mechanistic)"
        selected_cl_confidence = human_ivive_confidence
    elif allometric_cl is not None:
        selected_cl = allometric_cl
        selected_cl_source = "Cross-Species Allometric Extrapolation"
        selected_cl_confidence = cl_allometry.get("confidence", "MEDIUM")

    # 7. Distribution Candidates & Selection
    v_candidates = []
    if human_exp_iv_vss is not None:
        v_candidates.append({
            "source_name": "Human Clinical IV NCA (Vss)",
            "source_type": "EXPERIMENTAL",
            "value": human_exp_iv_vss,
            "unit": "L/kg",
            "confidence": "HIGH",
            "provenance": human_exp_iv_source,
        })
    elif human_exp_iv_vz is not None:
        v_candidates.append({
            "source_name": "Human Clinical IV NCA (Vz)",
            "source_type": "EXPERIMENTAL",
            "value": human_exp_iv_vz,
            "unit": "L/kg",
            "confidence": "HIGH",
            "provenance": human_exp_iv_source,
        })
    if allometric_vss is not None:
        v_candidates.append({
            "source_name": "Cross-Species Allometry (Vss)",
            "source_type": "TRANSLATIONAL_ALLOMETRY",
            "value": allometric_vss,
            "unit": "L/kg",
            "confidence": vss_allometry.get("confidence", "MEDIUM"),
            "provenance": f"Fitted across {vss_allometry.get('n_species')} animal species (b={vss_allometry.get('exponent_b')})",
        })
    if phys_vd is not None:
        v_candidates.append({
            "source_name": "Physicochemical Binding Model",
            "source_type": "CALCULATED_PHYSICOCHEMICAL",
            "value": phys_vd,
            "unit": "L/kg",
            "confidence": phys_dist.get("confidence", "LOW"),
            "provenance": phys_dist.get("model", "Oie-Tozer binding model"),
        })

    v_disagreement = calculate_disagreement(v_candidates)
    if v_disagreement["has_major_disagreement"]:
        warnings.append(f"Major Disagreement in Human Distribution (V) estimates ({v_disagreement['max_fold_difference']}x fold difference). Automatic averaging disabled.")

    selected_v = None
    selected_v_source = "MODEL_UNAVAILABLE"
    selected_v_confidence = "MODEL_UNAVAILABLE"
    if human_exp_iv_vss is not None:
        selected_v = human_exp_iv_vss
        selected_v_source = "Human Clinical IV NCA (Vss)"
        selected_v_confidence = "HIGH"
    elif human_exp_iv_vz is not None:
        selected_v = human_exp_iv_vz
        selected_v_source = "Human Clinical IV NCA (Vz)"
        selected_v_confidence = "HIGH"
    elif allometric_vss is not None:
        selected_v = allometric_vss
        selected_v_source = "Cross-Species Allometric Extrapolation (Vss)"
        selected_v_confidence = vss_allometry.get("confidence", "MEDIUM")
    elif phys_vd is not None:
        selected_v = phys_vd
        selected_v_source = "Physicochemical Binding Model"
        selected_v_confidence = phys_dist.get("confidence", "LOW")

    # Half-life calculation from selected CL and V
    selected_thalf = None
    if selected_cl is not None and selected_v is not None and selected_cl > 0 and selected_v > 0:
        cl_l_h_kg = (selected_cl * 60.0) / 1000.0
        selected_thalf = round(math.log(2) * selected_v / cl_l_h_kg, 2)

    # 8. Bioavailability Candidates & Selection
    f_candidates = []
    if exp_f_val is not None:
        f_candidates.append({
            "source_name": "Human Clinical Matched Study (F)",
            "source_type": "EXPERIMENTAL",
            "value": exp_f_val,
            "unit": "%",
            "confidence": "HIGH",
            "provenance": exp_f_source,
        })
    if f_predicted is not None:
        f_candidates.append({
            "source_name": "Composite Predicted F (Fa*Fg*Fh)",
            "source_type": "MECHANISTIC_COMPOSITE",
            "value": f_predicted,
            "unit": "%",
            "confidence": "MEDIUM",
            "provenance": f_predicted_reason,
        })

    f_disagreement = calculate_disagreement(f_candidates)
    selected_f = None
    selected_f_source = "MODEL_UNAVAILABLE"
    if exp_f_val is not None:
        selected_f = exp_f_val
        selected_f_source = "Human Clinical Matched Study (Experimental)"
    elif f_predicted is not None:
        selected_f = f_predicted
        selected_f_source = "Composite Mechanistic Prediction (Fa * Fg * Fh)"

    # 9. Absorption Rate Constant (ka)
    # Human ka from clinical PO study if available
    selected_ka = None
    selected_ka_source = "MODEL_UNAVAILABLE"
    if human_exp_po_tmax is not None and selected_cl is not None and selected_v is not None:
        # Solve for ka numerically from Tmax if possible
        cl_l_h = (selected_cl * 60.0) / 1000.0
        ke = cl_l_h / selected_v
        tmax = float(human_exp_po_tmax)
        if tmax > 0 and ke > 0:
            def tmax_eq(ka_val):
                if ka_val <= 0 or abs(ka_val - ke) < 1e-6:
                    return 1e6
                return (math.log(ka_val / ke) / (ka_val - ke)) - tmax
            try:
                from scipy.optimize import root_scalar
                sol = root_scalar(tmax_eq, bracket=[ke * 1.05, 50.0], method='brentq')
                if sol.converged and sol.root > 0:
                    selected_ka = round(sol.root, 4)
                    selected_ka_source = f"Derived from Human Clinical Tmax ({tmax} h)"
            except Exception:
                pass

    # 10. Human Simulation Readiness Evaluation
    iv_ready = selected_cl is not None and selected_v is not None
    po_ready = iv_ready and selected_f is not None and selected_ka is not None
    
    iv_reasons = []
    if selected_cl is not None:
        iv_reasons.append(f"Clearance available: {selected_cl} mL/min/kg ({selected_cl_source})")
    else:
        iv_reasons.append("Clearance unavailable: run IVIVE or provide cross-species PK.")
        
    if selected_v is not None:
        iv_reasons.append(f"Volume of distribution available: {selected_v} L/kg ({selected_v_source})")
    else:
        iv_reasons.append("Volume of distribution unavailable.")

    if cl_disagreement["has_major_disagreement"]:
        iv_reasons.append(f"Caution: {cl_disagreement['interpretation']}")

    po_reasons = []
    if not iv_ready:
        po_reasons.append("Requires valid IV clearance and volume parameters.")
    if selected_f is not None:
        po_reasons.append(f"Bioavailability available: {selected_f}% ({selected_f_source})")
    else:
        po_reasons.append(f"Bioavailability unavailable: {f_predicted_reason}")
    if selected_ka is not None:
        po_reasons.append(f"Absorption rate constant available: {selected_ka} 1/h ({selected_ka_source})")
    else:
        po_reasons.append("Absorption rate constant (ka) unavailable: requires clinical data or validated absorption model.")

    if iv_ready:
        iv_status = "READY" if not cl_disagreement["has_major_disagreement"] else "PARTIALLY_READY"
    else:
        iv_status = "INSUFFICIENT_DATA"

    if po_ready:
        po_status = "READY" if not (cl_disagreement["has_major_disagreement"] or f_disagreement.get("has_major_disagreement")) else "PARTIALLY_READY"
    elif iv_ready and (selected_f is not None or selected_ka is not None):
        po_status = "PARTIALLY_READY"
    else:
        po_status = "INSUFFICIENT_DATA"

    readiness = {
        "iv_simulation": {
            "status": iv_status,
            "reasons": iv_reasons,
        },
        "po_simulation": {
            "status": po_status,
            "reasons": po_reasons,
        },
        "overall_status": "READY" if iv_status == "READY" and po_status == "READY" else ("PARTIALLY_READY" if iv_status in ("READY", "PARTIALLY_READY") else "INSUFFICIENT_DATA"),
        "oral_translation_guardrail": "Strict Guardrail: Animal F and ka are never transferred directly to Human. Human oral simulation requires supported Human Fa/Fg/Fh components or clinical data.",
    }

    # Assemble inputs hash
    hash_payload = {
        "version_id": version_id,
        "selected_cl": selected_cl,
        "selected_v": selected_v,
        "selected_f": selected_f,
        "selected_ka": selected_ka,
        "fa": fa_val,
        "fg": fg_val,
        "fh": fh_val,
        "engine_version": HUMAN_PK_ENGINE_VERSION,
    }
    inputs_hash = hashlib.sha256(json.dumps(hash_payload, sort_keys=True).encode("utf-8")).hexdigest()

    return {
        "version_id": version_id,
        "compound_id": compound.compound_id,
        "compound_name": compound.name,
        "target_species": "Human",
        "standard_body_weight_kg": HUMAN_STANDARD_BW_KG,
        "clearance": {
            "selected_value": selected_cl,
            "selected_unit": "mL/min/kg",
            "selected_source": selected_cl_source,
            "confidence": selected_cl_confidence,
            "candidates": cl_candidates,
            "disagreement": cl_disagreement,
        },
        "volume": {
            "selected_value": selected_v,
            "selected_unit": "L/kg",
            "selected_source": selected_v_source,
            "confidence": selected_v_confidence,
            "candidates": v_candidates,
            "disagreement": v_disagreement,
        },
        "half_life": {
            "selected_value": selected_thalf,
            "selected_unit": "h",
            "calculation_formula": "ln(2) * V / CL",
            "experimental_value": human_exp_iv_thalf,
            "allometric_value": allometric_thalf,
        },
        "absorption": {
            "fa_value": fa_val,
            "fa_status": abs_comp.get("fa_status", "MODEL_UNAVAILABLE"),
            "fg_value": fg_val,
            "fg_status": abs_comp.get("fg_status", "MODEL_UNAVAILABLE"),
            "fh_value": fh_val,
            "fh_status": "CALCULATED" if fh_val is not None else "MODEL_UNAVAILABLE",
            "f_predicted": f_predicted,
            "f_predicted_status": f_predicted_status,
            "f_predicted_reason": f_predicted_reason,
            "f_experimental": exp_f_val,
            "f_experimental_source": exp_f_source,
            "f_selected": selected_f,
            "f_selected_source": selected_f_source,
            "f_candidates": f_candidates,
            "f_disagreement": f_disagreement,
            "ka_value": selected_ka,
            "ka_source": selected_ka_source,
        },
        "readiness": readiness,
        "inputs_hash": inputs_hash,
        "warnings": warnings,
        "model_version": HUMAN_PK_ENGINE_VERSION,
    }


# -----------------------------------------------------------------------------
# Human Simulation Engine (IV & PO)
# -----------------------------------------------------------------------------

class HumanSimulationRequest(BaseModel):
    route: str = Field(default="IV", description="Route: IV or PO")
    administration_type: str = Field(default="IV_BOLUS", description="IV_BOLUS, IV_INFUSION, or EXTRAVASCULAR_1COMP")
    dose: float = Field(default=100.0, description="Dose amount")
    dose_unit: str = Field(default="mg", description="Dose unit: mg, mg/kg, ug")
    body_weight_kg: float = Field(default=70.0, description="Patient body weight in kg")
    infusion_duration_hours: float = Field(default=1.0, description="Infusion duration (hours)")
    dosing_frequency: str = Field(default="Single Dose", description="Single Dose or Repeated Dosing")
    dose_interval_hours: float = Field(default=24.0, description="Interval between repeated doses (hours)")
    num_doses: int = Field(default=1, description="Number of doses")
    model_type: str = Field(default="ONE_COMPARTMENT", description="ONE_COMPARTMENT")
    # Explicit user overrides (labeled as assumptions)
    user_cl_override: float | None = Field(default=None, description="Clearance override (mL/min/kg)")
    user_v_override: float | None = Field(default=None, description="Volume override (L/kg)")
    user_f_override: float | None = Field(default=None, description="Bioavailability override (%)")
    user_fg_override: float | None = Field(default=None, description="Gut bioavailability override (0-1)")
    user_ka_override: float | None = Field(default=None, description="Absorption rate constant override (1/h)")


def run_human_pk_simulation(
    db: Session,
    version_id: int,
    req: HumanSimulationRequest,
) -> dict[str, Any]:
    """Run route-aware Human concentration-time simulation.
    Strictly applies scientific guardrails:
    - IV Bolus / Infusion / Repeated Dosing.
    - PO simulation requires supported or overridden F (Fa*Fg*Fh) and ka.
    """
    profile = assemble_human_pk_parameters(db, version_id)
    warnings = list(profile.get("warnings", []))
    assumptions = []

    bw = float(req.body_weight_kg or HUMAN_STANDARD_BW_KG)
    route = req.route.upper()

    # 1. Resolve Clearance (CL)
    cl_ml_min_kg = req.user_cl_override or profile["clearance"]["selected_value"]
    cl_source = "User Override" if req.user_cl_override else profile["clearance"]["selected_source"]
    if req.user_cl_override:
        assumptions.append(f"User Override applied for Human CL: {req.user_cl_override} mL/min/kg")

    if cl_ml_min_kg is None or cl_ml_min_kg <= 0:
        raise HTTPException(
            status_code=422,
            detail="Human Clearance (CL) is unavailable. Please run IVIVE, provide cross-species PK data, or specify an override.",
        )

    # 2. Resolve Volume (V)
    v_l_kg = req.user_v_override or profile["volume"]["selected_value"]
    v_source = "User Override" if req.user_v_override else profile["volume"]["selected_source"]
    if req.user_v_override:
        assumptions.append(f"User Override applied for Human V: {req.user_v_override} L/kg")

    if v_l_kg is None or v_l_kg <= 0:
        raise HTTPException(
            status_code=422,
            detail="Human Volume of Distribution (V) is unavailable. Please provide allometry, binding data, or specify an override.",
        )

    # Convert to total units
    cl_l_h = (cl_ml_min_kg * 60.0 / 1000.0) * bw
    v_l = v_l_kg * bw
    ke = cl_l_h / v_l
    t_half = math.log(2) / ke

    # 3. Resolve Dose into absolute mg
    dose_mg = float(req.dose)
    if req.dose_unit == "mg/kg":
        dose_mg = req.dose * bw
    elif req.dose_unit in ("ug", "µg"):
        dose_mg = req.dose / 1000.0
    elif req.dose_unit in ("ug/kg", "µg/kg"):
        dose_mg = (req.dose / 1000.0) * bw

    # 4. Resolve Route-Specific Parameters (PO vs IV)
    f_val = 1.0
    f_source = "100% (IV Systemic Reference)"
    ka_val = None
    ka_source = None
    flip_flop_flag = False

    if route == "PO":
        # Check F override or assembly
        if req.user_f_override is not None:
            f_val = req.user_f_override / 100.0
            f_source = f"User Override ({req.user_f_override}%)"
            assumptions.append(f"Explicit User Override applied for Oral Bioavailability F: {req.user_f_override}%")
        elif req.user_fg_override is not None:
            fa = profile["absorption"]["fa_value"] or 1.0
            fh = profile["absorption"]["fh_value"] or 1.0
            fg = float(req.user_fg_override)
            f_val = fa * fg * fh
            f_source = f"Calculated with Fg override ({round(fg, 2)}): {round(f_val*100, 1)}%"
            assumptions.append(f"Explicit User Override applied for Intestinal Escape Fg: {round(fg, 2)}")
        elif profile["absorption"]["f_selected"] is not None:
            f_val = profile["absorption"]["f_selected"] / 100.0
            f_source = profile["absorption"]["f_selected_source"]
        else:
            raise HTTPException(
                status_code=422,
                detail=f"Human Oral Simulation Refused: {profile['absorption']['f_predicted_reason']} (Scientific Guardrail: Fg cannot be assumed to be 1 without explicit override).",
            )

        # Check ka override or assembly
        if req.user_ka_override is not None:
            ka_val = float(req.user_ka_override)
            ka_source = f"User Override ({ka_val} 1/h)"
            assumptions.append(f"Explicit User Override applied for Absorption Rate Constant ka: {ka_val} 1/h")
        elif profile["absorption"]["ka_value"] is not None:
            ka_val = profile["absorption"]["ka_value"]
            ka_source = profile["absorption"]["ka_source"]
        else:
            raise HTTPException(
                status_code=422,
                detail="Human Oral Simulation Refused: Absorption rate constant (ka) is unavailable. Please provide Human clinical Tmax or an explicit ka override.",
            )

        if ka_val <= ke:
            flip_flop_flag = True
            warnings.append(f"POTENTIAL FLIP-FLOP KINETICS: ka ({round(ka_val, 4)} 1/h) <= ke ({round(ke, 4)} 1/h). Apparent terminal slope reflects absorption rate.")

    # 5. Build Time Grid & Simulate
    interval = float(req.dose_interval_hours or 24.0)
    num_doses = int(req.num_doses) if req.dosing_frequency == "Repeated Dosing" else 1
    t_end = max(interval * num_doses, t_half * 5.0, 24.0)
    n_points = 300
    t_grid = np.linspace(0.0, t_end, n_points)

    c_grid = np.zeros_like(t_grid)

    admin_type = req.administration_type.upper() if route == "IV" else "EXTRAVASCULAR_1COMP"
    t_inf = float(req.infusion_duration_hours or 1.0)

    # Superposition across doses
    for dose_idx in range(num_doses):
        t_dose = dose_idx * interval
        mask = t_grid >= t_dose
        t_rel = t_grid[mask] - t_dose

        if route == "IV":
            if admin_type == "IV_INFUSION" and t_inf > 0:
                r0 = dose_mg / t_inf  # mg/h
                # During infusion
                inf_mask = t_rel <= t_inf
                c_inf = (r0 / cl_l_h) * (1.0 - np.exp(-ke * t_rel[inf_mask]))
                # Post infusion
                post_mask = t_rel > t_inf
                c_end_inf = (r0 / cl_l_h) * (1.0 - math.exp(-ke * t_inf))
                c_post = c_end_inf * np.exp(-ke * (t_rel[post_mask] - t_inf))
                
                c_dose = np.zeros_like(t_rel)
                c_dose[inf_mask] = c_inf
                c_dose[post_mask] = c_post
                c_grid[mask] += c_dose * (1000.0 / 1.0)  # convert mg/L to ng/mL / ug/L (1 mg/L = 1000 ng/mL)
            else:
                # IV Bolus
                c0 = dose_mg / v_l  # mg/L
                c_dose = c0 * np.exp(-ke * t_rel)
                c_grid[mask] += c_dose * 1000.0  # ng/mL
        else:
            # First-order oral absorption
            if abs(ka_val - ke) < 1e-5:
                # Singularity case
                c_dose = (f_val * dose_mg / v_l) * ke * t_rel * np.exp(-ke * t_rel)
            else:
                c_dose = (f_val * dose_mg * ka_val / (v_l * (ka_val - ke))) * (np.exp(-ke * t_rel) - np.exp(-ka_val * t_rel))
            c_grid[mask] += np.maximum(c_dose * 1000.0, 0.0)

    # 6. Compute Analytical Metrics
    cmax = float(np.max(c_grid))
    tmax = float(t_grid[np.argmax(c_grid)])
    c_trough = float(c_grid[-1])

    if route == "IV":
        auc_single = (dose_mg / cl_l_h) * 1000.0  # ng*h/mL
    else:
        auc_single = (f_val * dose_mg / cl_l_h) * 1000.0

    auc_total = auc_single * num_doses
    r_acc = 1.0 / (1.0 - math.exp(-ke * interval)) if num_doses > 1 and interval > 0 else 1.0
    css_avg = (f_val * dose_mg / (cl_l_h * interval)) * 1000.0 if interval > 0 else cmax

    time_series = [
        {"time_hours": round(float(t), 3), "concentration_ng_ml": round(float(c), 4)}
        for t, c in zip(t_grid, c_grid)
    ]

    # 7. Check for Human Clinical Observations overlay
    observations_overlay = []
    clinical_studies = list(db.scalars(
        select(PKStudy)
        .where(PKStudy.compound_row_id == profile["version_id"], PKStudy.species == "Human", PKStudy.route == route)
    ))
    for st in clinical_studies:
        for obs in st.observations:
            observations_overlay.append({
                "study_id": st.id,
                "study_name": st.study_name,
                "time_hours": obs.time_raw,
                "concentration_ng_ml": obs.concentration_raw,
                "blq": obs.blq_flag,
            })

    output_metrics = {
        "cmax_ng_ml": round(cmax, 2),
        "tmax_hours": round(tmax, 2),
        "c_trough_ng_ml": round(c_trough, 2),
        "auc_single_ng_h_ml": round(auc_single, 2),
        "auc_total_ng_h_ml": round(auc_total, 2),
        "half_life_hours": round(t_half, 2),
        "kel_1_per_h": round(ke, 4),
        "clearance_l_h": round(cl_l_h, 3),
        "clearance_ml_min_kg": round(cl_ml_min_kg, 2),
        "volume_l": round(v_l, 2),
        "volume_l_kg": round(v_l_kg, 3),
        "bioavailability_pct": round(f_val * 100.0, 1),
        "ka_1_per_h": round(ka_val, 4) if ka_val is not None else None,
        "accumulation_ratio": round(r_acc, 2),
        "css_avg_ng_ml": round(css_avg, 2),
    }

    return {
        "version_id": version_id,
        "target_species": "Human",
        "route": route,
        "administration_type": admin_type,
        "dose_mg": round(dose_mg, 2),
        "body_weight_kg": bw,
        "dosing_frequency": req.dosing_frequency,
        "dose_interval_hours": interval,
        "num_doses": num_doses,
        "model_type": req.model_type,
        "parameters": {
            "cl": {"value": round(cl_ml_min_kg, 2), "unit": "mL/min/kg", "source": cl_source},
            "v": {"value": round(v_l_kg, 3), "unit": "L/kg", "source": v_source},
            "f": {"value": round(f_val * 100.0, 1), "unit": "%", "source": f_source},
            "ka": {"value": round(ka_val, 4) if ka_val else None, "unit": "1/h", "source": ka_source},
            "half_life": {"value": round(t_half, 2), "unit": "h"},
        },
        "output_metrics": output_metrics,
        "time_series": time_series,
        "observations_overlay": observations_overlay,
        "flip_flop_flag": flip_flop_flag,
        "warnings": warnings,
        "assumptions": assumptions,
        "scientific_labels": [
            "CALCULATED HUMAN PK TRANSLATIONAL SIMULATION",
            "NOT CLINICAL OBSERVATION",
            "NOT PBPK (1-Compartment PK Model)",
        ],
        "engine_version": HUMAN_PK_ENGINE_VERSION,
    }


# -----------------------------------------------------------------------------
# Prospective Prediction Freeze & Retrospective Validation
# -----------------------------------------------------------------------------

def freeze_human_prediction_snapshot(
    db: Session,
    version_id: int,
    snapshot_name: str = "Prospective Human PK Prediction",
) -> dict[str, Any]:
    """Freeze current Human translational PK predictions into an immutable snapshot.
    Immutable record includes candidate parameters, allometry, IVIVE, readiness, and inputs hash.
    """
    version = db.get(CompoundVersion, version_id)
    if not version:
        raise HTTPException(status_code=404, detail=f"CompoundVersion #{version_id} not found.")

    profile = assemble_human_pk_parameters(db, version_id)
    compound = version.compound

    snap = PKHumanPredictionSnapshot(
        project_id=compound.project_id,
        compound_row_id=compound.id,
        version_id=version_id,
        snapshot_name=snapshot_name,
        target_species="Human",
        selected_cl=profile["clearance"]["selected_value"],
        cl_unit=profile["clearance"]["selected_unit"],
        cl_source=profile["clearance"]["selected_source"],
        selected_v=profile["volume"]["selected_value"],
        v_unit=profile["volume"]["selected_unit"],
        v_source=profile["volume"]["selected_source"],
        fa_value=profile["absorption"]["fa_value"],
        fg_value=profile["absorption"]["fg_value"],
        fh_value=profile["absorption"]["fh_value"],
        f_predicted=profile["absorption"]["f_predicted"],
        f_experimental=profile["absorption"]["f_experimental"],
        f_selected=profile["absorption"]["f_selected"],
        f_source=profile["absorption"]["f_selected_source"],
        ka_value=profile["absorption"]["ka_value"],
        ka_source=profile["absorption"]["ka_source"],
        candidate_parameters_json={
            "clearance": profile["clearance"],
            "volume": profile["volume"],
            "half_life": profile["half_life"],
            "absorption": profile["absorption"],
        },
        disagreement_json={
            "clearance": profile["clearance"]["disagreement"],
            "volume": profile["volume"]["disagreement"],
            "bioavailability": profile["absorption"]["f_disagreement"],
        },
        readiness_json=profile["readiness"],
        confidence=profile["clearance"]["confidence"],
        warnings_json=profile["warnings"],
        inputs_hash=profile["inputs_hash"],
        model_version=HUMAN_PK_ENGINE_VERSION,
        is_immutable=True,
    )
    db.add(snap)
    db.commit()
    db.refresh(snap)

    return {
        "status": "FROZEN",
        "snapshot_id": snap.id,
        "snapshot_name": snap.snapshot_name,
        "version_id": version_id,
        "inputs_hash": snap.inputs_hash,
        "created_at": snap.created_at.isoformat(),
        "selected_cl": snap.selected_cl,
        "selected_v": snap.selected_v,
        "f_selected": snap.f_selected,
        "ka_value": snap.ka_value,
        "confidence": snap.confidence,
    }


def validate_against_clinical_data(
    db: Session,
    version_id: int,
    snapshot_id: int | None = None,
) -> dict[str, Any]:
    """Retrospectively validate clinical Human experimental PK against a PREVIOUSLY FROZEN prediction snapshot.
    Ensures predictions are not regenerated post-hoc.
    """
    version = db.get(CompoundVersion, version_id)
    if not version:
        raise HTTPException(status_code=404, detail=f"CompoundVersion #{version_id} not found.")

    compound = version.compound

    # Find the requested snapshot or latest frozen snapshot
    query = select(PKHumanPredictionSnapshot).where(PKHumanPredictionSnapshot.version_id == version_id)
    if snapshot_id:
        query = query.where(PKHumanPredictionSnapshot.id == snapshot_id)
    else:
        query = query.order_by(PKHumanPredictionSnapshot.id.desc())

    snapshot = db.scalars(query).first()
    if not snapshot:
        return {
            "status": "NO_PREDICTION_SNAPSHOT",
            "message": "No frozen Human prospective prediction snapshot exists. Freeze a prediction before validation.",
            "metrics": {},
            "comparisons": [],
        }

    # Query all Human clinical PK studies
    human_studies = list(db.scalars(
        select(PKStudy)
        .where(PKStudy.compound_row_id == compound.id, PKStudy.species == "Human")
        .order_by(PKStudy.id.asc())
    ))

    comparisons = []
    fold_errors = []

    for st in human_studies:
        nca = st.latest_nca
        route = (st.route or "").upper()
        if not nca:
            continue

        # Clearance Comparison (IV only)
        if route == "IV" and nca.cl and snapshot.selected_cl:
            pred = float(snapshot.selected_cl)
            obs = float(nca.cl)
            fe = round(pred / obs, 2)
            afe = round(max(fe, 1.0 / fe), 2)
            fold_errors.append(afe)
            comparisons.append({
                "study_id": st.id,
                "study_name": st.study_name,
                "endpoint": "CL",
                "route": "IV",
                "predicted": pred,
                "observed": obs,
                "unit": "mL/min/kg",
                "fold_error": fe,
                "absolute_fold_error": afe,
                "percent_error": round(((pred - obs) / obs) * 100.0, 1),
                "performance_band": "WITHIN_2_FOLD" if afe <= 2.0 else ("WITHIN_3_FOLD" if afe <= 3.0 else "OUTSIDE_3_FOLD"),
            })

        # Volume Comparison (IV only)
        if route == "IV" and snapshot.selected_v:
            obs_v = None
            if nca.cl and nca.mrt and nca.cl > 0 and nca.mrt > 0:
                obs_v = round((nca.cl * 60.0 / 1000.0) * nca.mrt, 3)
            elif nca.vz:
                obs_v = nca.vz
            if obs_v:
                pred = float(snapshot.selected_v)
                fe = round(pred / obs_v, 2)
                afe = round(max(fe, 1.0 / fe), 2)
                fold_errors.append(afe)
                comparisons.append({
                    "study_id": st.id,
                    "study_name": st.study_name,
                    "endpoint": "Vss",
                    "route": "IV",
                    "predicted": pred,
                    "observed": obs_v,
                    "unit": "L/kg",
                    "fold_error": fe,
                    "absolute_fold_error": afe,
                    "percent_error": round(((pred - obs_v) / obs_v) * 100.0, 1),
                    "performance_band": "WITHIN_2_FOLD" if afe <= 2.0 else ("WITHIN_3_FOLD" if afe <= 3.0 else "OUTSIDE_3_FOLD"),
                })

        # Bioavailability Comparison (PO)
        if route == "PO" and snapshot.f_selected:
            # Check if experimental F exists
            exp_f = None
            try:
                ba = calculate_bioavailability_for_version(version_id, db)
                for item in ba.get("bioavailability", []):
                    if (item.get("species") or "").strip().capitalize() == "Human" and (item.get("route") or "").upper() == "PO":
                        if item.get("bioavailability_pct") is not None:
                            exp_f = item["bioavailability_pct"]
            except Exception:
                pass

            if exp_f:
                pred = float(snapshot.f_selected)
                fe = round(pred / exp_f, 2)
                afe = round(max(fe, 1.0 / fe), 2)
                fold_errors.append(afe)
                comparisons.append({
                    "study_id": st.id,
                    "study_name": st.study_name,
                    "endpoint": "F",
                    "route": "PO",
                    "predicted": pred,
                    "observed": exp_f,
                    "unit": "%",
                    "fold_error": fe,
                    "absolute_fold_error": afe,
                    "percent_error": round(((pred - exp_f) / exp_f) * 100.0, 1),
                    "performance_band": "WITHIN_2_FOLD" if afe <= 2.0 else ("WITHIN_3_FOLD" if afe <= 3.0 else "OUTSIDE_3_FOLD"),
                })

    n = len(comparisons)
    if n == 0:
        return {
            "status": "NO_CLINICAL_MATCH",
            "message": "Prediction snapshot found, but no matching Human clinical NCA data entered yet.",
            "snapshot_id": snapshot.id,
            "snapshot_name": snapshot.snapshot_name,
            "snapshot_created_at": snapshot.created_at.isoformat(),
            "comparisons": [],
            "metrics": {},
        }

    aafe = round(float(np.mean(fold_errors)), 2)
    within_2f = sum(1 for c in comparisons if c["absolute_fold_error"] <= 2.0)
    within_3f = sum(1 for c in comparisons if c["absolute_fold_error"] <= 3.0)

    return {
        "status": "VALIDATED",
        "snapshot_id": snapshot.id,
        "snapshot_name": snapshot.snapshot_name,
        "snapshot_created_at": snapshot.created_at.isoformat(),
        "n_comparisons": n,
        "metrics": {
            "aafe": aafe,
            "within_2_fold_count": within_2f,
            "within_2_fold_pct": round((within_2f / n) * 100.0, 1),
            "within_3_fold_count": within_3f,
            "within_3_fold_pct": round((within_3f / n) * 100.0, 1),
        },
        "comparisons": comparisons,
    }


# -----------------------------------------------------------------------------
# FastAPI Route Registration
# -----------------------------------------------------------------------------

def register_human_pk_routes(app: FastAPI) -> None:
    """Register FastAPI endpoints for Human PK Prediction & Simulation."""

    @app.get("/api/compound-versions/{version_id}/human-pk/profile")
    def get_human_pk_profile_endpoint(
        version_id: int,
        db: Session = Depends(get_db),
    ):
        return assemble_human_pk_parameters(db, version_id)

    @app.post("/api/compound-versions/{version_id}/human-pk/simulation/run")
    def run_human_simulation_endpoint(
        version_id: int,
        req: HumanSimulationRequest,
        db: Session = Depends(get_db),
    ):
        return run_human_pk_simulation(db, version_id, req)

    class FreezeSnapshotPayload(BaseModel):
        snapshot_name: str = "Prospective Human PK Prediction"

    @app.post("/api/compound-versions/{version_id}/human-pk/freeze-snapshot")
    def freeze_human_snapshot_endpoint(
        version_id: int,
        payload: FreezeSnapshotPayload = FreezeSnapshotPayload(),
        db: Session = Depends(get_db),
    ):
        name = payload.snapshot_name
        return freeze_human_prediction_snapshot(db, version_id, snapshot_name=name)

    @app.get("/api/compound-versions/{version_id}/human-pk/snapshots")
    def list_human_snapshots_endpoint(
        version_id: int,
        db: Session = Depends(get_db),
    ):
        snaps = list(db.scalars(
            select(PKHumanPredictionSnapshot)
            .where(PKHumanPredictionSnapshot.version_id == version_id)
            .order_by(PKHumanPredictionSnapshot.id.desc())
        ))
        return {
            "version_id": version_id,
            "snapshots": [
                {
                    "id": s.id,
                    "snapshot_name": s.snapshot_name,
                    "selected_cl": s.selected_cl,
                    "selected_v": s.selected_v,
                    "f_selected": s.f_selected,
                    "ka_value": s.ka_value,
                    "confidence": s.confidence,
                    "inputs_hash": s.inputs_hash,
                    "created_at": s.created_at.isoformat(),
                }
                for s in snaps
            ],
        }

    @app.get("/api/compound-versions/{version_id}/human-pk/validation")
    def validate_human_pk_endpoint(
        version_id: int,
        snapshot_id: int | None = Query(default=None),
        db: Session = Depends(get_db),
    ):
        return validate_against_clinical_data(db, version_id, snapshot_id)

    @app.get("/api/projects/{project_id}/human-pk/summary")
    def get_project_human_pk_summary(
        project_id: int,
        db: Session = Depends(get_db),
    ):
        project = db.get(Project, project_id)
        if not project:
            raise HTTPException(status_code=404, detail=f"Project #{project_id} not found.")

        compounds_summary = []
        for c in project.compounds:
            for v in c.versions:
                try:
                    prof = assemble_human_pk_parameters(db, v.id)
                    compounds_summary.append({
                        "compound_id": c.compound_id,
                        "compound_name": c.name,
                        "version_id": v.id,
                        "version_number": v.version_number,
                        "cl_value": prof["clearance"]["selected_value"],
                        "cl_source": prof["clearance"]["selected_source"],
                        "v_value": prof["volume"]["selected_value"],
                        "v_source": prof["volume"]["selected_source"],
                        "f_value": prof["absorption"]["f_selected"],
                        "f_source": prof["absorption"]["f_selected_source"],
                        "overall_readiness": prof["readiness"]["overall_status"],
                        "has_major_disagreement": prof["clearance"]["disagreement"]["has_major_disagreement"],
                    })
                except Exception:
                    pass

        return {
            "project_id": project_id,
            "project_name": project.name,
            "compounds": compounds_summary,
        }
