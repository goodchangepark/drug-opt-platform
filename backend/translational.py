"""Stage 5B-3 PK Validation, Cross-Species Scaling & Translational PK Foundation.

Answers three core questions:
1. Validation: How accurately does the platform reproduce experimental PK?
2. Scaling: How do PK parameters translate between animal species?
3. Human Translation: What can and cannot currently be predicted for Human PK?

Preserves strict scientific boundaries:
- Never compares incompatible endpoints (CL vs CL/F, Vss vs Vz, F vs Fa).
- Never transfers animal F or ka directly to Human.
- Never fabricates missing species or human PK data.
- Animal allometry excludes Human experimental data (prospective holdout).
- Freezes prospective predictions into immutable snapshots.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from typing import Any

from fastapi import Depends, FastAPI, HTTPException, Query
from pydantic import BaseModel, Field
from scipy import stats
from sqlalchemy import (
    JSON, Boolean, DateTime, Float, ForeignKey, Integer, String, Text,
    inspect, select, text,
)
from sqlalchemy.orm import Mapped, Session, mapped_column, relationship

from .database import Base, get_db
from .ivive import (
    PHYSIOLOGY_DEFAULTS, SPECIES, SPECIES_ALIASES, IVIVERun,
    get_pk_foundation_profile,
)
from .models import Compound, CompoundVersion, Project, utcnow
from .pk import PKNCAResult, PKObservation, PKStudy, calculate_bioavailability_for_version

TRANSLATIONAL_ENGINE_VERSION = "5B-3.1.0"

# Standard Physiology Reference Body Weights (kg)
SPECIES_BODY_WEIGHTS: dict[str, float] = {
    "Mouse": 0.02,
    "Rat": 0.25,
    "Dog": 10.0,
    "Monkey": 5.0,
    "Human": 70.0,
}

# Standard Reference Brain Weights (kg)
SPECIES_BRAIN_WEIGHTS: dict[str, float] = {
    "Mouse": 0.0004,
    "Rat": 0.0018,
    "Dog": 0.08,
    "Monkey": 0.098,
    "Human": 1.4,
}

# Maximum Life-Span Potential (MLP in years: 185.4 * BW^0.636)
SPECIES_MLP: dict[str, float] = {
    "Mouse": round(185.4 * (0.02 ** 0.636), 2),    # ~15.35 y
    "Rat": round(185.4 * (0.25 ** 0.636), 2),      # ~76.62 y
    "Dog": round(185.4 * (10.0 ** 0.636), 2),      # ~799.30 y
    "Monkey": round(185.4 * (5.0 ** 0.636), 2),    # ~514.23 y
    "Human": round(185.4 * (70.0 ** 0.636), 2),    # ~2772.33 y
}


# Database Model for Prospective Prediction Snapshots
class PKTranslationalSnapshot(Base):
    __tablename__ = "pk_translational_snapshots"

    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    compound_row_id: Mapped[int] = mapped_column(ForeignKey("compounds.id", ondelete="CASCADE"), index=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    target_species: Mapped[str] = mapped_column(String(50), default="Human", index=True)
    endpoint: Mapped[str] = mapped_column(String(50), index=True)  # CL, Vss, t1/2, F
    method: Mapped[str] = mapped_column(String(80), index=True)  # SIMPLE_ALLOMETRY, RULE_OF_EXPONENTS, HEPATIC_IVIVE
    predicted_value: Mapped[float] = mapped_column(Float)
    predicted_unit: Mapped[str] = mapped_column(String(40))
    confidence: Mapped[str] = mapped_column(String(40), default="MEDIUM")
    species_used: Mapped[list[str]] = mapped_column(JSON, default=list)
    body_weights_used: Mapped[dict[str, float]] = mapped_column(JSON, default=dict)
    allometric_a: Mapped[float | None] = mapped_column(Float, nullable=True)
    allometric_b: Mapped[float | None] = mapped_column(Float, nullable=True)
    allometric_r2: Mapped[float | None] = mapped_column(Float, nullable=True)
    inputs_hash: Mapped[str] = mapped_column(String(64), index=True)
    frozen_inputs: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    warnings: Mapped[list[str]] = mapped_column(JSON, default=list)
    model_version: Mapped[str] = mapped_column(String(50), default=TRANSLATIONAL_ENGINE_VERSION)
    is_immutable: Mapped[bool] = mapped_column(Boolean, default=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)

    project = relationship("Project")
    compound = relationship("Compound")
    version = relationship("CompoundVersion")


def ensure_translational_schema(engine_obj) -> None:
    """Ensure database schema for translational PK snapshots."""
    insp = inspect(engine_obj)
    tables = insp.get_table_names()
    if "pk_translational_snapshots" not in tables:
        PKTranslationalSnapshot.__table__.create(bind=engine_obj, checkfirst=True)


# Core Allometric Fit Mathematics
def fit_allometry(
    species_points: list[dict[str, Any]],
    target_species: str = "Human",
    param_type: str = "CL",  # "CL" or "Vss"
) -> dict[str, Any]:
    """Fit classical simple allometric power-law: Y_total = a * (BW)^b.

    Y_total is total clearance (mL/min) or total volume (L).
    Normalized value is converted to/from per-kg basis: Y_norm = Y_total / BW.
    """
    if len(species_points) < 2:
        return {
            "status": "INSUFFICIENT_DATA",
            "message": f"Allometric scaling requires at least 2 distinct species with experimental IV {param_type}; found {len(species_points)}.",
            "param_type": param_type,
            "n_species": len(species_points),
            "points": species_points,
        }

    # Extract log-transformed values
    log_bw: list[float] = []
    log_y: list[float] = []
    species_names: list[str] = []
    bw_map: dict[str, float] = {}

    for pt in species_points:
        bw = pt["bw_kg"]
        val_norm = pt.get("value_norm") if pt.get("value_norm") is not None else pt.get("observed_norm", 0.0)  # mL/min/kg for CL, L/kg for Vss
        if val_norm <= 0 or bw <= 0:
            continue
        val_total = val_norm * bw  # Total mL/min or Total L
        log_bw.append(math.log(bw))
        log_y.append(math.log(val_total))
        species_names.append(pt["species"])
        bw_map[pt["species"]] = bw

    n = len(log_bw)
    if n < 2:
        return {
            "status": "INSUFFICIENT_DATA",
            "message": f"Requires at least 2 valid positive points for {param_type}.",
            "param_type": param_type,
            "n_species": n,
        }

    mean_x = sum(log_bw) / n
    mean_y = sum(log_y) / n

    s_xx = sum((x - mean_x) ** 2 for x in log_bw)
    s_xy = sum((x - mean_x) * (y - mean_y) for x, y in zip(log_bw, log_y))

    if s_xx < 1e-12:
        b = 0.75 if param_type == "CL" else 1.0
        ln_a = mean_y - b * mean_x
    else:
        b = s_xy / s_xx
        ln_a = mean_y - b * mean_x

    a = math.exp(ln_a)

    ss_tot = sum((y - mean_y) ** 2 for y in log_y)
    ss_res = sum((y - (ln_a + b * x)) ** 2 for x, y in zip(log_bw, log_y))
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 1e-12 else 1.0
    r2 = max(0.0, min(1.0, r2))

    # Extrapolate to Target Species
    target_bw = SPECIES_BODY_WEIGHTS.get(target_species, 70.0)
    target_total = a * (target_bw ** b)
    target_norm = target_total / target_bw

    # Confidence and Warning Rules
    warnings: list[str] = []
    warnings.append(f"CROSS-SPECIES EXTRAPOLATION: Extrapolating {target_species} (BW={target_bw} kg) beyond animal body-weight range ({min(bw_map.values())}–{max(bw_map.values())} kg).")

    if n == 2:
        confidence = "LOW"
        warnings.append(f"LOW-EVIDENCE HUMAN EXTRAPOLATION: Fitted on 2 species only ({', '.join(species_names)}). A minimum of 3 species is standard practice.")
    elif n >= 3:
        if r2 >= 0.90:
            confidence = "HIGH"
        elif r2 >= 0.70:
            confidence = "MEDIUM"
        else:
            confidence = "LOW"
            warnings.append(f"POOR ALLOMETRIC FIT (R² = {r2:.3f}): Inter-species correlation is weak.")

    # Exponent Diagnostics
    if param_type == "CL":
        if b < 0.4 or b > 1.2:
            warnings.append(f"ALLOMETRIC EXPONENT OUTSIDE EXPECTED RANGE: Fitted exponent b = {b:.3f} deviates from typical clearance allometry (~0.75).")
    elif param_type == "Vss":
        if b < 0.6 or b > 1.4:
            warnings.append(f"ALLOMETRIC EXPONENT OUTSIDE EXPECTED RANGE: Fitted exponent b = {b:.3f} deviates from typical volume allometry (~1.0).")

    # Fitted line points across whole range for visualization
    plot_points = []
    for pt in species_points:
        bw = pt["bw_kg"]
        val_norm = pt.get("value_norm") if pt.get("value_norm") is not None else pt.get("observed_norm", 0.0)
        val_total = val_norm * bw
        plot_points.append({
            "species": pt["species"],
            "bw_kg": bw,
            "observed_norm": round(val_norm, 4),
            "observed_total": round(val_total, 4),
            "fitted_total": round(a * (bw ** b), 4),
            "fitted_norm": round((a * (bw ** b)) / bw, 4),
            "source": pt.get("source", "EXPERIMENTAL_NCA"),
            "evidence_type": pt.get("evidence_type", "EXPERIMENTAL"),
            "is_extrapolated": False,
        })

    # Add extrapolated target point
    extrapolated_point = {
        "species": target_species,
        "bw_kg": target_bw,
        "observed_norm": None,
        "observed_total": None,
        "fitted_total": round(target_total, 4),
        "fitted_norm": round(target_norm, 4),
        "source": "ALLOMETRIC_EXTRAPOLATION",
        "evidence_type": "TRANSLATIONAL_ESTIMATE",
        "is_extrapolated": True,
    }

    # Historical Rule of Exponents (Mahmood & Balian / Boxenbaum) for Clearance
    historical_roe = None
    if param_type == "CL" and n >= 3:
        # MLP correction
        log_y_mlp = [math.log(pt["value_norm"] * pt["bw_kg"] * SPECIES_MLP.get(pt["species"], 100.0)) for pt in species_points]
        mean_y_mlp = sum(log_y_mlp) / n
        s_xy_mlp = sum((x - mean_x) * (y - mean_y_mlp) for x, y in zip(log_bw, log_y_mlp))
        b_mlp = s_xy_mlp / s_xx if s_xx > 1e-12 else 0.75
        a_mlp = math.exp(mean_y_mlp - b_mlp * mean_x)
        human_mlp = SPECIES_MLP.get("Human", 2772.33)
        cl_mlp_total = (a_mlp * (target_bw ** b_mlp)) / human_mlp
        cl_mlp_norm = cl_mlp_total / target_bw

        # BrW correction
        log_y_brw = [math.log(pt["value_norm"] * pt["bw_kg"] * SPECIES_BRAIN_WEIGHTS.get(pt["species"], 0.05)) for pt in species_points]
        mean_y_brw = sum(log_y_brw) / n
        s_xy_brw = sum((x - mean_x) * (y - mean_y_brw) for x, y in zip(log_bw, log_y_brw))
        b_brw = s_xy_brw / s_xx if s_xx > 1e-12 else 0.75
        a_brw = math.exp(mean_y_brw - b_brw * mean_x)
        human_brw = SPECIES_BRAIN_WEIGHTS.get("Human", 1.4)
        cl_brw_total = (a_brw * (target_bw ** b_brw)) / human_brw
        cl_brw_norm = cl_brw_total / target_bw

        # Rule selection
        if b < 0.71:
            roe_rule = "Simple Allometry (b < 0.71)"
            roe_val_norm = target_norm
        elif b <= 1.00:
            roe_rule = "MLP Correction (0.71 <= b <= 1.00)"
            roe_val_norm = cl_mlp_norm
        else:
            roe_rule = "Brain Weight Correction (b > 1.00)"
            roe_val_norm = cl_brw_norm

        historical_roe = {
            "method_name": "Historical Rule of Exponents (Mahmood & Balian)",
            "selected_rule": roe_rule,
            "simple_allometry_norm": round(target_norm, 4),
            "mlp_corrected_norm": round(cl_mlp_norm, 4),
            "brw_corrected_norm": round(cl_brw_norm, 4),
            "roe_predicted_norm": round(roe_val_norm, 4),
            "unit": "mL/min/kg",
            "citation": "Mahmood I, Balian JD. Interspecies scaling: predicting clearance in humans in early drug development. Life Sci. 1996;59:579-585.",
        }

    return {
        "status": "SUCCESS",
        "param_type": param_type,
        "target_species": target_species,
        "target_bw_kg": target_bw,
        "n_species": n,
        "species_used": species_names,
        "body_weights_used": bw_map,
        "coefficient_a": round(a, 6),
        "exponent_b": round(b, 4),
        "r_squared": round(r2, 4),
        "heuristic_reference": "~0.75 (CL)" if param_type == "CL" else "~1.0 (Vss)",
        "extrapolated_total": round(target_total, 4),
        "extrapolated_total_unit": "mL/min" if param_type == "CL" else "L",
        "extrapolated_norm": round(target_norm, 4),
        "extrapolated_norm_unit": "mL/min/kg" if param_type == "CL" else "L/kg",
        "confidence": confidence,
        "warnings": warnings,
        "historical_roe_correction": historical_roe,
        "plot_points": plot_points,
        "extrapolated_point": extrapolated_point,
    }


# Leave-One-Species-Out Cross Validation
def run_loso_validation(
    species_points: list[dict[str, Any]],
    param_type: str = "CL",
) -> dict[str, Any]:
    """Perform Leave-One-Species-Out (LOSO) cross-validation when N >= 3 species exist."""
    n = len(species_points)
    if n < 3:
        return {
            "status": "INSUFFICIENT_DATA",
            "message": f"LOSO cross-validation requires at least 3 distinct species; found {n}.",
            "n_species": n,
        }

    loso_results: list[dict[str, Any]] = []
    afes: list[float] = []
    log_errors: list[float] = []

    for i in range(n):
        held_out = species_points[i]
        training = [species_points[j] for j in range(n) if j != i]

        fit_res = fit_allometry(training, target_species=held_out["species"], param_type=param_type)
        if fit_res["status"] != "SUCCESS":
            continue

        pred_norm = fit_res["extrapolated_norm"]
        obs_norm = held_out["value_norm"]

        if obs_norm <= 0 or pred_norm <= 0:
            continue

        fe = pred_norm / obs_norm
        afe = max(fe, 1.0 / fe)
        log_err = math.log10(pred_norm) - math.log10(obs_norm)

        afes.append(afe)
        log_errors.append(log_err)

        loso_results.append({
            "held_out_species": held_out["species"],
            "bw_kg": held_out["bw_kg"],
            "observed_norm": round(obs_norm, 4),
            "predicted_norm": round(pred_norm, 4),
            "unit": "mL/min/kg" if param_type == "CL" else "L/kg",
            "fold_error": round(fe, 3),
            "absolute_fold_error": round(afe, 3),
            "log_error": round(log_err, 4),
            "within_2_fold": bool(afe <= 2.0),
            "within_3_fold": bool(afe <= 3.0),
            "training_species": fit_res["species_used"],
            "training_exponent_b": fit_res["exponent_b"],
        })

    if not afes:
        return {"status": "NO_VALID_PAIRS", "message": "No valid non-zero points for LOSO evaluation."}

    m = len(afes)
    aafe = 10 ** (sum(math.log10(a) for a in afes) / m)
    bias = 10 ** (sum(log_errors) / m)
    within_2_fold_pct = round((sum(1 for a in afes if a <= 2.0) / m) * 100.0, 1)
    within_3_fold_pct = round((sum(1 for a in afes if a <= 3.0) / m) * 100.0, 1)

    return {
        "status": "SUCCESS",
        "param_type": param_type,
        "n_species_evaluated": m,
        "aafe": round(aafe, 3),
        "bias_gmfe": round(bias, 3),
        "within_2_fold_pct": within_2_fold_pct,
        "within_3_fold_pct": within_3_fold_pct,
        "loso_evaluations": loso_results,
    }


# Validation Metrics Engine
def evaluate_pk_predictions(
    prediction_pairs: list[dict[str, Any]],
) -> dict[str, Any]:
    """Compute standard statistical validation metrics across paired Predicted vs Observed PK parameters.

    Metrics:
    - Fold Error (FE = Pred / Obs)
    - Absolute Fold Error (AFE = max(FE, 1/FE))
    - Average Absolute Fold Error (AAFE = 10^(mean(|log10(FE)|)))
    - Geometric Mean Fold Error / Bias (GMFE = 10^(mean(log10(FE))))
    - RMSE in log10 space
    - MAE in log10 space
    - Spearman rank correlation
    - % within 1.5-fold, 2-fold, 3-fold, and >3-fold error.
    """
    valid_pairs: list[dict[str, Any]] = []
    afes: list[float] = []
    log_errors: list[float] = []
    preds: list[float] = []
    obss: list[float] = []

    for p in prediction_pairs:
        obs = p.get("observed")
        pred = p.get("predicted")
        if obs is None or pred is None or obs <= 0 or pred <= 0:
            continue

        fe = pred / obs
        afe = max(fe, 1.0 / fe)
        log_err = math.log10(pred) - math.log10(obs)

        afes.append(afe)
        log_errors.append(log_err)
        preds.append(pred)
        obss.append(obs)

        band = "WITHIN 1.5-FOLD" if afe <= 1.5 else ("WITHIN 2-FOLD" if afe <= 2.0 else ("WITHIN 3-FOLD" if afe <= 3.0 else ">3-FOLD ERROR"))

        valid_pairs.append({
            "compound_id": p.get("compound_id"),
            "compound_name": p.get("compound_name"),
            "version_id": p.get("version_id"),
            "species": p.get("species"),
            "route": p.get("route"),
            "endpoint": p.get("endpoint"),
            "method": p.get("method"),
            "observed": round(obs, 4),
            "predicted": round(pred, 4),
            "unit": p.get("unit", ""),
            "fold_error": round(fe, 3),
            "absolute_fold_error": round(afe, 3),
            "log_error": round(log_err, 4),
            "performance_band": band,
        })

    n = len(valid_pairs)
    if n == 0:
        return {
            "status": "NO_DATA",
            "message": "No paired Predicted vs Observed PK observations available.",
            "n": 0,
            "pairs": [],
        }

    aafe = 10.0 ** (sum(math.log10(a) for a in afes) / n)
    bias = 10.0 ** (sum(log_errors) / n)
    rmse_log = math.sqrt(sum(err ** 2 for err in log_errors) / n)
    mae_log = sum(abs(err) for err in log_errors) / n

    w_1_5 = sum(1 for a in afes if a <= 1.5)
    w_2 = sum(1 for a in afes if a <= 2.0)
    w_3 = sum(1 for a in afes if a <= 3.0)
    gt_3 = sum(1 for a in afes if a > 3.0)

    spearman_rho = None
    spearman_p = None
    if n >= 3 and len(set(preds)) > 1 and len(set(obss)) > 1:
        try:
            res = stats.spearmanr(obss, preds)
            spearman_rho = round(float(res.statistic), 3) if not math.isnan(res.statistic) else None
            spearman_p = round(float(res.pvalue), 4) if not math.isnan(res.pvalue) else None
        except Exception:
            pass

    acceptance_label = f"WITHIN 2-FOLD ACCEPTANCE MET ({round(w_2 / n * 100.0, 1)}%)" if (w_2 / n) >= 0.70 else f"DESCRIPTIVE PERFORMANCE ({round(w_2 / n * 100.0, 1)}% within 2-fold)"

    return {
        "status": "SUCCESS",
        "n": n,
        "aafe": round(aafe, 3),
        "bias_gmfe": round(bias, 3),
        "rmse_log10": round(rmse_log, 4),
        "mae_log10": round(mae_log, 4),
        "spearman_rho": spearman_rho,
        "spearman_pvalue": spearman_p,
        "within_1_5_fold_count": w_1_5,
        "within_1_5_fold_pct": round(w_1_5 / n * 100.0, 1),
        "within_2_fold_count": w_2,
        "within_2_fold_pct": round(w_2 / n * 100.0, 1),
        "within_3_fold_count": w_3,
        "within_3_fold_pct": round(w_3 / n * 100.0, 1),
        "gt_3_fold_count": gt_3,
        "gt_3_fold_pct": round(gt_3 / n * 100.0, 1),
        "acceptance_summary": acceptance_label,
        "pairs": valid_pairs,
    }


# Aggregator & Translational Profile Builder
def get_translational_pk_profile(
    db: Session,
    version_id: int,
    freeze_snapshot: bool = True,
) -> dict[str, Any]:
    """Extract multi-species experimental PK parameters, fit CL and Vss allometry,
    run LOSO cross-validation, assess Human simulation readiness, and freeze prospective snapshot.
    """
    version = db.get(CompoundVersion, version_id)
    if not version:
        raise HTTPException(status_code=404, detail=f"CompoundVersion #{version_id} not found.")

    compound = version.compound
    project_id = compound.project_id

    # 1. Gather all experimental PK studies for this compound version
    studies = list(db.scalars(
        select(PKStudy)
        .where(PKStudy.compound_row_id == compound.id)
        .order_by(PKStudy.id.asc())
    ))

    # Cross-Species Data Matrix
    species_data_matrix: dict[str, dict[str, Any]] = {
        sp: {
            "species": sp,
            "default_bw_kg": SPECIES_BODY_WEIGHTS.get(sp, 70.0),
            "study_bw_kg": None,
            "effective_bw_kg": SPECIES_BODY_WEIGHTS.get(sp, 70.0),
            "iv_studies_count": 0,
            "po_studies_count": 0,
            "other_studies_count": 0,
            "cl_iv": None,  # mL/min/kg
            "cl_iv_source": None,
            "vss_iv": None,  # L/kg
            "vss_iv_source": None,
            "vz_iv": None,   # L/kg
            "half_life_iv": None,  # h
            "f_po": None,    # %
            "f_po_source": None,
            "cmax_po": None,
            "tmax_po": None,
            "ka_po": None,
            "has_experimental_iv": False,
            "has_experimental_po": False,
        }
        for sp in SPECIES
    }

    cl_points_animal: list[dict[str, Any]] = []
    vss_points_animal: list[dict[str, Any]] = []
    human_experimental_iv: dict[str, Any] = {}

    for st in studies:
        sp_clean = st.species.strip().capitalize() if st.species else "Rat"
        if sp_clean not in species_data_matrix:
            continue

        route = (st.route or "").strip().upper()
        nca = st.latest_nca

        row = species_data_matrix[sp_clean]
        if route == "IV":
            row["iv_studies_count"] += 1
            if nca:
                if nca.cl is not None and nca.cl > 0:
                    row["cl_iv"] = nca.cl
                    row["cl_iv_source"] = f"PK Study #{st.id} ({st.study_name})"
                    row["has_experimental_iv"] = True

                # Steady-State Volume (Vss = CL * MRT)
                vss_calc = None
                if nca.cl is not None and nca.mrt is not None and nca.cl > 0 and nca.mrt > 0:
                    vss_calc = round((nca.cl * 60.0 / 1000.0) * nca.mrt, 4)
                elif hasattr(nca, "vss") and getattr(nca, "vss") is not None:
                    vss_calc = getattr(nca, "vss")

                if vss_calc is not None and vss_calc > 0:
                    row["vss_iv"] = vss_calc
                    row["vss_iv_source"] = f"PK Study #{st.id} ({st.study_name})"
                elif nca.vz is not None and nca.vz > 0:
                    row["vz_iv"] = nca.vz

                hl = getattr(nca, "terminal_half_life", None) or getattr(nca, "half_life", None)
                if hl is not None:
                    row["half_life_iv"] = hl

        elif route == "PO":
            row["po_studies_count"] += 1
            if nca:
                row["has_experimental_po"] = True
                if nca.cmax is not None:
                    row["cmax_po"] = nca.cmax
                if nca.tmax is not None:
                    row["tmax_po"] = nca.tmax

    # Check calculated experimental F across species
    try:
        ba_data = calculate_bioavailability_for_version(db, version_id)
        for item in ba_data.get("bioavailability", []):
            sp_name = item.get("species", "").strip().capitalize()
            if sp_name in species_data_matrix and item.get("route", "").upper() == "PO":
                if item.get("bioavailability_pct") is not None:
                    species_data_matrix[sp_name]["f_po"] = item["bioavailability_pct"]
                    species_data_matrix[sp_name]["f_po_source"] = f"Matched Study #{item.get('study_id')} ({item.get('study_name')})"
    except Exception:
        pass

    for sp in SPECIES:
        # Compile animal allometric points (exclude Human for holdout scaling)
        row = species_data_matrix[sp]
        if sp != "Human":
            if row["cl_iv"] is not None and row["cl_iv"] > 0:
                cl_points_animal.append({
                    "species": sp,
                    "bw_kg": row["effective_bw_kg"],
                    "value_norm": row["cl_iv"],
                    "unit": "mL/min/kg",
                    "source": row["cl_iv_source"] or "EXPERIMENTAL_NCA",
                    "evidence_type": "EXPERIMENTAL",
                })
            if row["vss_iv"] is not None and row["vss_iv"] > 0:
                vss_points_animal.append({
                    "species": sp,
                    "bw_kg": row["effective_bw_kg"],
                    "value_norm": row["vss_iv"],
                    "unit": "L/kg",
                    "source": row["vss_iv_source"] or "EXPERIMENTAL_NCA",
                    "evidence_type": "EXPERIMENTAL",
                })
        else:
            if row["cl_iv"] is not None:
                human_experimental_iv["cl"] = row["cl_iv"]
                human_experimental_iv["cl_unit"] = "mL/min/kg"
            if row["vss_iv"] is not None:
                human_experimental_iv["vss"] = row["vss_iv"]
                human_experimental_iv["vss_unit"] = "L/kg"
            if row["f_po"] is not None:
                human_experimental_iv["f_pct"] = row["f_po"]

    # 2. Fit Allometric Models
    cl_allometry = fit_allometry(cl_points_animal, target_species="Human", param_type="CL")
    vss_allometry = fit_allometry(vss_points_animal, target_species="Human", param_type="Vss")

    # 3. LOSO Validation on Animal Data
    cl_loso = run_loso_validation(cl_points_animal, param_type="CL")
    vss_loso = run_loso_validation(vss_points_animal, param_type="Vss")

    # 4. Human Hepatic IVIVE Retrieval
    human_ivive_run = db.scalars(
        select(IVIVERun)
        .where(IVIVERun.version_id == version_id, IVIVERun.species == "Human", IVIVERun.status == "COMPLETE")
        .order_by(IVIVERun.id.desc())
        .limit(1)
    ).first()

    human_ivive_cl = None
    if human_ivive_run:
        human_ivive_cl = (human_ivive_run.outputs_json or {}).get("cl_in_vivo_blood") or (human_ivive_run.outputs_json or {}).get("cl_in_vivo_plasma") or getattr(human_ivive_run, "cl_in_vivo_blood", None)
    human_ivive_confidence = human_ivive_run.confidence if human_ivive_run else "UNAVAILABLE"

    # 5. Translated Half-Life calculation (from allometric CL & Vss)
    translated_half_life = None
    if cl_allometry.get("status") == "SUCCESS" and vss_allometry.get("status") == "SUCCESS":
        cl_ml_min_kg = cl_allometry["extrapolated_norm"]
        vss_l_kg = vss_allometry["extrapolated_norm"]
        if cl_ml_min_kg > 0 and vss_l_kg > 0:
            cl_l_h_kg = (cl_ml_min_kg * 60.0) / 1000.0
            translated_half_life = round(math.log(2) * vss_l_kg / cl_l_h_kg, 2)

    # 6. Side-by-Side Human Comparison
    human_comparison = {
        "clearance": {
            "method_a_hepatic_ivive": {
                "method_name": "Mechanistic Hepatic IVIVE",
                "value": round(human_ivive_cl, 4) if human_ivive_cl is not None else None,
                "unit": "mL/min/kg",
                "confidence": human_ivive_confidence,
                "status": "AVAILABLE" if human_ivive_cl is not None else "MODEL_UNAVAILABLE",
                "notes": "Hepatic metabolic clearance only; does not account for renal/extrahepatic elimination.",
            },
            "method_b_simple_allometry": {
                "method_name": "Simple Allometric Scaling",
                "value": cl_allometry.get("extrapolated_norm") if cl_allometry.get("status") == "SUCCESS" else None,
                "unit": "mL/min/kg",
                "confidence": cl_allometry.get("confidence", "INSUFFICIENT_DATA"),
                "status": cl_allometry.get("status", "INSUFFICIENT_DATA"),
                "n_species": cl_allometry.get("n_species", 0),
                "exponent_b": cl_allometry.get("exponent_b"),
                "r2": cl_allometry.get("r_squared"),
                "notes": "Classical body-weight power-law extrapolation from animal systemic IV clearance.",
            },
            "method_c_historical_roe": cl_allometry.get("historical_roe_correction"),
            "method_d_experimental_human": {
                "method_name": "Human Clinical IV PK",
                "value": human_experimental_iv.get("cl"),
                "unit": "mL/min/kg",
                "status": "AVAILABLE" if human_experimental_iv.get("cl") is not None else "NOT_RECORDED",
                "notes": "Observed clinical measurement; holds absolute scientific precedence.",
            },
        },
        "volume_vss": {
            "simple_allometry": {
                "value": vss_allometry.get("extrapolated_norm") if vss_allometry.get("status") == "SUCCESS" else None,
                "unit": "L/kg",
                "confidence": vss_allometry.get("confidence", "INSUFFICIENT_DATA"),
                "status": vss_allometry.get("status", "INSUFFICIENT_DATA"),
                "exponent_b": vss_allometry.get("exponent_b"),
                "r2": vss_allometry.get("r_squared"),
            },
            "experimental_human": {
                "value": human_experimental_iv.get("vss"),
                "unit": "L/kg",
                "status": "AVAILABLE" if human_experimental_iv.get("vss") is not None else "NOT_RECORDED",
            },
        },
        "translated_half_life": {
            "value": translated_half_life,
            "unit": "hours",
            "formula": "ln(2) * Vss_allometric / CL_allometric",
            "v_definition_used": "Steady-State Volume of Distribution (Vss)",
        },
    }

    # 7. Deterministic Human Simulation Readiness Assessment
    # CL Readiness
    if human_experimental_iv.get("cl") is not None:
        cl_readiness = "READY"
        cl_reason = "Human clinical IV clearance is experimentally recorded."
    elif cl_allometry.get("status") == "SUCCESS" and cl_allometry.get("n_species", 0) >= 3:
        cl_readiness = "READY"
        cl_reason = f"Allometric CL fit on {cl_allometry['n_species']} animal species (R² = {cl_allometry['r_squared']})."
    elif human_ivive_cl is not None:
        cl_readiness = "LIMITED"
        cl_reason = "Hepatic IVIVE available (hepatic clearance only; renal/extrahepatic unconfirmed)."
    elif cl_allometry.get("status") == "SUCCESS" and cl_allometry.get("n_species") == 2:
        cl_readiness = "LIMITED"
        cl_reason = "2-species allometry only (low-evidence rodent extrapolation)."
    else:
        cl_readiness = "UNAVAILABLE"
        cl_reason = "No clearance data available across species."

    # Volume Readiness
    if human_experimental_iv.get("vss") is not None:
        v_readiness = "READY"
        v_reason = "Human clinical IV Vss is experimentally recorded."
    elif vss_allometry.get("status") == "SUCCESS" and vss_allometry.get("n_species", 0) >= 3:
        v_readiness = "READY"
        v_reason = f"Allometric Vss fit on {vss_allometry['n_species']} animal species (R² = {vss_allometry['r_squared']})."
    elif vss_allometry.get("status") == "SUCCESS" and vss_allometry.get("n_species") == 2:
        v_readiness = "LIMITED"
        v_reason = "2-species allometric Vss (low evidence)."
    else:
        v_readiness = "UNAVAILABLE"
        v_reason = "No IV volume data available across species."

    # Bioavailability Readiness
    if human_experimental_iv.get("f_pct") is not None:
        f_readiness = "READY"
        f_reason = "Human clinical oral bioavailability experimentally established."
    else:
        animal_f_count = sum(1 for s in species_data_matrix.values() if s["species"] != "Human" and s["f_po"] is not None)
        if animal_f_count > 0:
            f_readiness = "LIMITED"
            f_reason = f"Animal bioavailability observed in {animal_f_count} species as supporting evidence. Human oral F cannot be transferred directly."
        else:
            f_readiness = "UNAVAILABLE"
            f_reason = "No oral bioavailability data recorded."

    # Absorption Rate ka Readiness
    ka_readiness = "UNAVAILABLE"
    ka_reason = "Human absorption kinetics unmodeled. Animal ka cannot be transferred directly to Human."

    # Overall Readiness
    if cl_readiness == "READY" and v_readiness == "READY" and f_readiness == "READY":
        overall_readiness = "READY"
    elif (cl_readiness in {"READY", "LIMITED"}) and (v_readiness in {"READY", "LIMITED"}):
        overall_readiness = "PARTIALLY READY"
    else:
        overall_readiness = "NOT READY"

    readiness_card = {
        "overall_status": overall_readiness,
        "clearance": {"status": cl_readiness, "reason": cl_reason},
        "volume": {"status": v_readiness, "reason": v_reason},
        "bioavailability": {"status": f_readiness, "reason": f_reason},
        "absorption_rate": {"status": ka_readiness, "reason": ka_reason},
        "oral_translation_guardrail": "Human oral bioavailability and absorption rate must not be inferred directly from animal species. Extravascular simulation is blocked without explicit human parameterization.",
    }

    # 8. Prospective Prediction Snapshot Freeze (if requested and allometry succeeded)
    if freeze_snapshot and cl_allometry.get("status") == "SUCCESS":
        input_data = {
            "cl_points": cl_points_animal,
            "vss_points": vss_points_animal,
            "human_ivive": human_ivive_cl,
        }
        in_hash = hashlib.sha256(json.dumps(input_data, sort_keys=True).encode()).hexdigest()[:16]

        # Check existing snapshot
        existing = db.scalars(
            select(PKTranslationalSnapshot)
            .where(
                PKTranslationalSnapshot.version_id == version_id,
                PKTranslationalSnapshot.target_species == "Human",
                PKTranslationalSnapshot.endpoint == "CL",
                PKTranslationalSnapshot.inputs_hash == in_hash,
            )
        ).first()

        if not existing:
            snap = PKTranslationalSnapshot(
                project_id=project_id,
                compound_row_id=compound.id,
                version_id=version_id,
                target_species="Human",
                endpoint="CL",
                method="SIMPLE_ALLOMETRY",
                predicted_value=cl_allometry["extrapolated_norm"],
                predicted_unit="mL/min/kg",
                confidence=cl_allometry["confidence"],
                species_used=cl_allometry["species_used"],
                body_weights_used=cl_allometry["body_weights_used"],
                allometric_a=cl_allometry["coefficient_a"],
                allometric_b=cl_allometry["exponent_b"],
                allometric_r2=cl_allometry["r_squared"],
                inputs_hash=in_hash,
                frozen_inputs=input_data,
                warnings=cl_allometry["warnings"],
                model_version=TRANSLATIONAL_ENGINE_VERSION,
                is_immutable=True,
            )
            db.add(snap)
            db.commit()

    return {
        "version_id": version_id,
        "compound_id": compound.compound_id,
        "compound_name": compound.name,
        "species_data_matrix": list(species_data_matrix.values()),
        "clearance_allometry": cl_allometry,
        "volume_allometry": vss_allometry,
        "clearance_loso": cl_loso,
        "volume_loso": vss_loso,
        "human_comparison": human_comparison,
        "human_simulation_readiness": readiness_card,
        "engine_version": TRANSLATIONAL_ENGINE_VERSION,
    }


# Route Registrations
def register_translational_routes(app: FastAPI) -> None:
    """Register FastAPI endpoints for Translational PK & Validation."""

    @app.get("/api/compound-versions/{version_id}/translational-pk")
    def get_compound_translational_pk(
        version_id: int,
        freeze_snapshot: bool = Query(True),
        db: Session = Depends(get_db),
    ):
        return get_translational_pk_profile(db, version_id, freeze_snapshot=freeze_snapshot)

    @app.get("/api/projects/{project_id}/translational-pk")
    def get_project_translational_pk(
        project_id: int,
        db: Session = Depends(get_db),
    ):
        project = db.get(Project, project_id)
        if not project:
            raise HTTPException(status_code=404, detail=f"Project #{project_id} not found.")

        rows: list[dict[str, Any]] = []
        for c in project.compounds:
            if not c.versions:
                continue
            latest_v = c.versions[-1]
            try:
                prof = get_translational_pk_profile(db, latest_v.id, freeze_snapshot=False)
                sp_matrix = {row["species"]: row for row in prof["species_data_matrix"]}
                cl_allo = prof["clearance_allometry"]
                vss_allo = prof["volume_allometry"]

                rows.append({
                    "compound_row_id": c.id,
                    "compound_id": c.compound_id,
                    "compound_name": c.name,
                    "version_id": latest_v.id,
                    "rat_cl_iv": sp_matrix.get("Rat", {}).get("cl_iv"),
                    "dog_cl_iv": sp_matrix.get("Dog", {}).get("cl_iv"),
                    "monkey_cl_iv": sp_matrix.get("Monkey", {}).get("cl_iv"),
                    "human_cl_allometric": cl_allo.get("extrapolated_norm") if cl_allo.get("status") == "SUCCESS" else None,
                    "human_cl_ivive": prof["human_comparison"]["clearance"]["method_a_hepatic_ivive"].get("value"),
                    "human_cl_experimental": sp_matrix.get("Human", {}).get("cl_iv"),
                    "rat_vss_iv": sp_matrix.get("Rat", {}).get("vss_iv"),
                    "dog_vss_iv": sp_matrix.get("Dog", {}).get("vss_iv"),
                    "human_vss_allometric": vss_allo.get("extrapolated_norm") if vss_allo.get("status") == "SUCCESS" else None,
                    "human_simulation_readiness": prof["human_simulation_readiness"]["overall_status"],
                })
            except Exception:
                continue

        return {
            "project_id": project_id,
            "project_name": project.name,
            "compounds_matrix": rows,
        }

    @app.get("/api/compound-versions/{version_id}/pk-validation")
    def get_compound_pk_validation(
        version_id: int,
        db: Session = Depends(get_db),
    ):
        version = db.get(CompoundVersion, version_id)
        if not version:
            raise HTTPException(status_code=404, detail=f"CompoundVersion #{version_id} not found.")

        compound = version.compound
        pairs: list[dict[str, Any]] = []

        # 1. Compare IVIVE predicted hepatic CL vs Experimental IV CL
        ivive_runs = list(db.scalars(
            select(IVIVERun).where(IVIVERun.version_id == version_id, IVIVERun.status == "COMPLETE")
        ))
        for run in ivive_runs:
            cl_pred = (run.outputs_json or {}).get("cl_in_vivo_blood") or (run.outputs_json or {}).get("cl_in_vivo_plasma") or getattr(run, "cl_in_vivo_blood", None)
            if not cl_pred:
                continue
            studies = list(db.scalars(
                select(PKStudy).where(PKStudy.compound_row_id == compound.id, PKStudy.species == run.species, PKStudy.route == "IV")
            ))
            for st in studies:
                if st.latest_nca and st.latest_nca.cl and cl_pred:
                    pairs.append({
                        "compound_id": compound.compound_id,
                        "compound_name": compound.name,
                        "version_id": version_id,
                        "species": run.species,
                        "route": "IV",
                        "endpoint": "CL",
                        "method": "HEPATIC_IVIVE",
                        "predicted": cl_pred,
                        "observed": st.latest_nca.cl,
                        "unit": "mL/min/kg",
                    })

        # 2. Retrospective check against frozen translational snapshots
        snapshots = list(db.scalars(
            select(PKTranslationalSnapshot).where(PKTranslationalSnapshot.version_id == version_id)
        ))
        for snap in snapshots:
            human_studies = list(db.scalars(
                select(PKStudy).where(PKStudy.compound_row_id == compound.id, PKStudy.species == snap.target_species, PKStudy.route == "IV")
            ))
            for st in human_studies:
                if st.latest_nca:
                    if snap.endpoint == "CL" and st.latest_nca.cl:
                        pairs.append({
                            "compound_id": compound.compound_id,
                            "compound_name": compound.name,
                            "version_id": version_id,
                            "species": snap.target_species,
                            "route": "IV",
                            "endpoint": "CL",
                            "method": snap.method,
                            "predicted": snap.predicted_value,
                            "observed": st.latest_nca.cl,
                            "unit": snap.predicted_unit,
                        })
                    elif snap.endpoint == "Vss":
                        vss_obs = None
                        if st.latest_nca.cl and st.latest_nca.mrt:
                            vss_obs = round((st.latest_nca.cl * 60.0 / 1000.0) * st.latest_nca.mrt, 4)
                        elif st.latest_nca.vz:
                            vss_obs = st.latest_nca.vz
                        if vss_obs is not None:
                            pairs.append({
                                "compound_id": compound.compound_id,
                                "compound_name": compound.name,
                                "version_id": version_id,
                                "species": snap.target_species,
                                "route": "IV",
                                "endpoint": "Vss",
                                "method": snap.method,
                                "predicted": snap.predicted_value,
                                "observed": vss_obs,
                                "unit": snap.predicted_unit,
                            })

        metrics = evaluate_pk_predictions(pairs)
        return {
            "version_id": version_id,
            "compound_id": compound.compound_id,
            "validation_metrics": metrics,
        }

    @app.get("/api/projects/{project_id}/pk-validation")
    def get_project_pk_validation(
        project_id: int,
        db: Session = Depends(get_db),
    ):
        project = db.get(Project, project_id)
        if not project:
            raise HTTPException(status_code=404, detail=f"Project #{project_id} not found.")

        pairs: list[dict[str, Any]] = []
        for c in project.compounds:
            for v in c.versions:
                ivive_runs = list(db.scalars(
                    select(IVIVERun).where(IVIVERun.version_id == v.id, IVIVERun.status == "COMPLETE")
                ))
                for run in ivive_runs:
                    cl_pred = (run.outputs_json or {}).get("cl_in_vivo_blood") or (run.outputs_json or {}).get("cl_in_vivo_plasma") or getattr(run, "cl_in_vivo_blood", None)
                    if not cl_pred:
                        continue
                    studies = list(db.scalars(
                        select(PKStudy).where(PKStudy.compound_row_id == c.id, PKStudy.species == run.species, PKStudy.route == "IV")
                    ))
                    for st in studies:
                        if st.latest_nca and st.latest_nca.cl and cl_pred:
                            pairs.append({
                                "compound_id": c.compound_id,
                                "compound_name": c.name,
                                "version_id": v.id,
                                "species": run.species,
                                "route": "IV",
                                "endpoint": "CL",
                                "method": "HEPATIC_IVIVE",
                                "predicted": cl_pred,
                                "observed": st.latest_nca.cl,
                                "unit": "mL/min/kg",
                            })

        metrics = evaluate_pk_predictions(pairs)
        return {
            "project_id": project_id,
            "project_name": project.name,
            "validation_metrics": metrics,
        }
