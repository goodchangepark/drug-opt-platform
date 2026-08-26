"""Deterministic Pharmacokinetic IV Concentration-Time Simulation Engine (Stage 5B-1)."""

from __future__ import annotations

import math
from datetime import datetime, timezone
from typing import Any

import numpy as np
from fastapi import Depends, HTTPException, Query
from pydantic import BaseModel
from scipy.optimize import curve_fit
from sqlalchemy import JSON, Boolean, DateTime, Float, ForeignKey, Integer, String, Text, inspect, select
from sqlalchemy.orm import Mapped, mapped_column, relationship, Session

from .database import Base, get_db
from .ivive import PKParameterSet, estimate_volume_of_distribution, get_pk_foundation_profile
from .models import Compound, CompoundVersion, Project, utcnow
from .pk import PKObservation, PKStudy

SIMULATION_ENGINE_NAME = "Stage 5B-1 Deterministic IV PK Simulation Engine"
SIMULATION_ENGINE_VERSION = "5B-1.0"


def ensure_simulation_schema(engine):
    inspector = inspect(engine)
    if "projects" not in inspector.get_table_names():
        return
    Base.metadata.create_all(
        bind=engine,
        tables=[
            PKSimulationRun.__table__,
        ],
    )


class PKSimulationRun(Base):
    __tablename__ = "pk_simulation_runs"
    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    compound_row_id: Mapped[int] = mapped_column(ForeignKey("compounds.id", ondelete="CASCADE"), index=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    species: Mapped[str] = mapped_column(String(100), default="Rat", index=True)
    route: Mapped[str] = mapped_column(String(40), default="IV", index=True)
    administration_type: Mapped[str] = mapped_column(String(40), default="IV_BOLUS", index=True)
    dose: Mapped[float] = mapped_column(Float, default=5.0)
    dose_unit: Mapped[str] = mapped_column(String(40), default="mg/kg")
    infusion_duration_hours: Mapped[float] = mapped_column(Float, default=0.0)
    dosing_frequency: Mapped[str] = mapped_column(String(60), default="Single Dose")
    dose_interval_hours: Mapped[float] = mapped_column(Float, default=24.0)
    num_doses: Mapped[int] = mapped_column(Integer, default=1)
    model_type: Mapped[str] = mapped_column(String(60), default="ONE_COMPARTMENT", index=True)
    parameter_snapshot: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    parameter_sources: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    simulation_settings: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    output_metrics: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    time_series: Mapped[list[dict[str, Any]]] = mapped_column(JSON, default=list)
    residuals: Mapped[list[dict[str, Any]]] = mapped_column(JSON, default=list)
    warnings: Mapped[list[str]] = mapped_column(JSON, default=list)
    confidence: Mapped[str] = mapped_column(String(40), default="MEDIUM")
    provenance: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)

    project = relationship("Project")
    compound = relationship("Compound")
    version = relationship("CompoundVersion")


# Pydantic Schemas
class PKSimulationRequest(BaseModel):
    species: str = "Rat"
    administration_type: str = "IV_BOLUS"  # "IV_BOLUS" or "IV_INFUSION"
    dose: float = 5.0
    dose_unit: str = "mg/kg"  # "mg/kg" or "µg/kg"
    infusion_duration_hours: float = 0.0
    dosing_frequency: str = "Single Dose"  # "Single Dose" or "Repeated Dosing"
    dose_interval_hours: float = 24.0
    num_doses: int = 1
    model_type: str = "ONE_COMPARTMENT"  # "ONE_COMPARTMENT" or "TWO_COMPARTMENT"
    user_cl_override: float | None = None
    user_v_override: float | None = None
    custom_t_end: float | None = None


# Core Calculations & Dimensional Checks
def canonicalize_units(dose: float, dose_unit: str, cl: float, cl_unit: str, v: float, v_unit: str) -> dict[str, float]:
    """Convert input parameters into canonical SI-derived units:

    Dose: mg/kg
    Clearance: L/h/kg
    Volume: L/kg
    Concentration conversion: 1 mg/L = 1000 ng/mL.
    Raises ValueError on dimensional mismatch or non-positive values.
    """
    if dose <= 0:
        raise ValueError(f"Dose must be strictly positive (>0), got {dose}")
    if cl <= 0:
        raise ValueError(f"Clearance must be strictly positive (>0), got {cl}")
    if v <= 0:
        raise ValueError(f"Volume of distribution must be strictly positive (>0), got {v}")

    # 1. Dose -> mg/kg
    if dose_unit in ("mg/kg", "mg"):
        dose_mg_kg = float(dose)
    elif dose_unit in ("µg/kg", "ug/kg", "µg", "ug"):
        dose_mg_kg = float(dose) / 1000.0
    else:
        raise ValueError(f"Unsupported dose unit '{dose_unit}'. Expected 'mg/kg' or 'µg/kg'.")

    # 2. Clearance -> L/h/kg
    if cl_unit in ("L/h/kg", "L/h"):
        cl_l_h_kg = float(cl)
    elif cl_unit in ("mL/min/kg", "mL/min"):
        cl_l_h_kg = float(cl) * 60.0 / 1000.0  # 1 mL/min = 0.06 L/h
    elif cl_unit in ("mL/h/kg", "mL/h"):
        cl_l_h_kg = float(cl) / 1000.0
    else:
        raise ValueError(f"Unsupported clearance unit '{cl_unit}'. Expected 'mL/min/kg' or 'L/h/kg'.")

    # 3. Volume -> L/kg
    if v_unit in ("L/kg", "L"):
        v_l_kg = float(v)
    elif v_unit in ("mL/kg", "mL"):
        v_l_kg = float(v) / 1000.0
    else:
        raise ValueError(f"Unsupported volume unit '{v_unit}'. Expected 'L/kg' or 'mL/kg'.")

    return {
        "dose_mg_kg": dose_mg_kg,
        "cl_l_h_kg": cl_l_h_kg,
        "v_l_kg": v_l_kg,
    }


def simulate_one_compartment_iv_bolus(
    dose_mg_kg: float,
    cl_l_h_kg: float,
    v_l_kg: float,
    num_doses: int = 1,
    dose_interval_h: float = 24.0,
    t_end_h: float = 24.0,
    num_points: int = 200,
) -> dict[str, Any]:
    """Single or repeated IV bolus 1-compartment linear PK simulation.

    C(t) = Dose / V * exp(-k * t)  [mg/L = 1000 ng/mL]
    Superposition for multiple doses.
    """
    k_elim = cl_l_h_kg / v_l_kg  # 1/h
    half_life_h = math.log(2.0) / k_elim
    c0_ng_ml = (dose_mg_kg / v_l_kg) * 1000.0
    analytical_auc_inf_ng_h_ml = (dose_mg_kg / cl_l_h_kg) * 1000.0

    # Determine time grid
    total_time_h = max(t_end_h, num_doses * dose_interval_h + 5 * half_life_h)
    dt = total_time_h / max(num_points, 100)

    time_series = []
    for i in range(num_points + 1):
        t = i * dt
        conc_ng_ml = 0.0
        for d in range(num_doses):
            t_dose = d * dose_interval_h
            if t >= t_dose:
                conc_ng_ml += c0_ng_ml * math.exp(-k_elim * (t - t_dose))
        time_series.append({"time": round(t, 4), "concentration": round(conc_ng_ml, 4), "unit": "ng/mL"})

    # Trapezoidal Numerical AUC
    auc_last_num = 0.0
    for i in range(len(time_series) - 1):
        t1, c1 = time_series[i]["time"], time_series[i]["concentration"]
        t2, c2 = time_series[i + 1]["time"], time_series[i + 1]["concentration"]
        auc_last_num += 0.5 * (c1 + c2) * (t2 - t1)

    c_last = time_series[-1]["concentration"]
    auc_inf_num = auc_last_num + (c_last / k_elim if k_elim > 0 else 0.0)
    agreement_pct = round(min(auc_inf_num, analytical_auc_inf_ng_h_ml) / max(auc_inf_num, analytical_auc_inf_ng_h_ml) * 100.0, 2) if analytical_auc_inf_ng_h_ml > 0 else 100.0

    return {
        "k_elim": round(k_elim, 6),
        "half_life_hours": round(half_life_h, 4),
        "c0_ng_ml": round(c0_ng_ml, 4),
        "cmax_ng_ml": round(c0_ng_ml if num_doses == 1 else max(p["concentration"] for p in time_series), 4),
        "tmax_hours": 0.0,
        "auc_last_ng_h_ml": round(auc_last_num, 2),
        "auc_inf_analytical_ng_h_ml": round(analytical_auc_inf_ng_h_ml, 2),
        "auc_inf_numerical_ng_h_ml": round(auc_inf_num, 2),
        "auc_agreement_pct": agreement_pct,
        "time_series": time_series,
    }


def simulate_one_compartment_iv_infusion(
    dose_mg_kg: float,
    infusion_duration_h: float,
    cl_l_h_kg: float,
    v_l_kg: float,
    num_doses: int = 1,
    dose_interval_h: float = 24.0,
    t_end_h: float = 24.0,
    num_points: int = 200,
) -> dict[str, Any]:
    """Single or repeated IV infusion 1-compartment linear PK simulation.

    During infusion (0 <= t <= Tinf): C(t) = (R0 / CL) * (1 - exp(-k * t))
    Post infusion (t > Tinf): C(t) = C(Tinf) * exp(-k * (t - Tinf))
    """
    if infusion_duration_h <= 0:
        raise ValueError("Infusion duration must be strictly positive (>0) for IV infusion.")

    k_elim = cl_l_h_kg / v_l_kg
    half_life_h = math.log(2.0) / k_elim
    r0_mg_kg_h = dose_mg_kg / infusion_duration_h  # mg/kg/h
    css_inf_ng_ml = (r0_mg_kg_h / cl_l_h_kg) * 1000.0
    c_tinf_ng_ml = css_inf_ng_ml * (1.0 - math.exp(-k_elim * infusion_duration_h))
    analytical_auc_inf_ng_h_ml = (dose_mg_kg / cl_l_h_kg) * 1000.0

    total_time_h = max(t_end_h, num_doses * dose_interval_h + 5 * half_life_h)
    dt = total_time_h / max(num_points, 100)

    time_grid = [i * dt for i in range(num_points + 1)]
    if infusion_duration_h not in time_grid:
        time_grid.append(infusion_duration_h)
    time_grid = sorted(set(round(t, 4) for t in time_grid if t <= total_time_h))

    time_series = []
    for t in time_grid:
        conc_ng_ml = 0.0
        for d in range(num_doses):
            t_start = d * dose_interval_h
            t_rel = t - t_start
            if t_rel >= 0:
                if t_rel <= infusion_duration_h:
                    conc_ng_ml += css_inf_ng_ml * (1.0 - math.exp(-k_elim * t_rel))
                else:
                    conc_ng_ml += c_tinf_ng_ml * math.exp(-k_elim * (t_rel - infusion_duration_h))
        time_series.append({"time": round(t, 4), "concentration": round(conc_ng_ml, 4), "unit": "ng/mL"})

    auc_last_num = 0.0
    for i in range(len(time_series) - 1):
        t1, c1 = time_series[i]["time"], time_series[i]["concentration"]
        t2, c2 = time_series[i + 1]["time"], time_series[i + 1]["concentration"]
        auc_last_num += 0.5 * (c1 + c2) * (t2 - t1)

    c_last = time_series[-1]["concentration"]
    auc_inf_num = auc_last_num + (c_last / k_elim if k_elim > 0 else 0.0)
    agreement_pct = round(min(auc_inf_num, analytical_auc_inf_ng_h_ml) / max(auc_inf_num, analytical_auc_inf_ng_h_ml) * 100.0, 2) if analytical_auc_inf_ng_h_ml > 0 else 100.0

    if num_doses == 1:
        cmax_val = c_tinf_ng_ml
        tmax_val = infusion_duration_h
    else:
        cmax_val = max(p["concentration"] for p in time_series)
        tmax_val = min((p["time"] for p in time_series if abs(p["concentration"] - cmax_val) < 1e-4), default=infusion_duration_h)

    return {
        "k_elim": round(k_elim, 6),
        "half_life_hours": round(half_life_h, 4),
        "r0_mg_kg_h": round(r0_mg_kg_h, 4),
        "c_tinf_ng_ml": round(c_tinf_ng_ml, 4),
        "cmax_ng_ml": round(cmax_val, 4),
        "tmax_hours": round(tmax_val, 4),
        "auc_last_ng_h_ml": round(auc_last_num, 2),
        "auc_inf_analytical_ng_h_ml": round(analytical_auc_inf_ng_h_ml, 2),
        "auc_inf_numerical_ng_h_ml": round(auc_inf_num, 2),
        "auc_agreement_pct": agreement_pct,
        "time_series": time_series,
    }


def simulate_two_compartment_iv_bolus(
    dose_mg_kg: float,
    cl_l_h_kg: float,
    vc_l_kg: float,
    q_l_h_kg: float,
    vp_l_kg: float,
    t_end_h: float = 24.0,
    num_points: int = 200,
) -> dict[str, Any]:
    """Forward IV bolus 2-compartment linear PK simulation given microconstants.

    C(t) = A * exp(-alpha * t) + B * exp(-beta * t)
    """
    k10 = cl_l_h_kg / vc_l_kg
    k12 = q_l_h_kg / vc_l_kg
    k21 = q_l_h_kg / vp_l_kg

    term1 = k10 + k12 + k21
    discriminant = term1**2 - 4.0 * k10 * k21
    if discriminant < 0:
        raise ValueError("Invalid 2-compartment microconstants (negative discriminant).")

    sqrt_disc = math.sqrt(discriminant)
    alpha = (term1 + sqrt_disc) / 2.0
    beta = (term1 - sqrt_disc) / 2.0

    if alpha == beta:
        raise ValueError("Degenerate 2-compartment rate constants (alpha == beta).")

    # Amplitudes in mg/L
    c0_total_mg_l = dose_mg_kg / vc_l_kg
    a_amp_mg_l = c0_total_mg_l * (alpha - k21) / (alpha - beta)
    b_amp_mg_l = c0_total_mg_l * (k21 - beta) / (alpha - beta)

    a_amp_ng_ml = a_amp_mg_l * 1000.0
    b_amp_ng_ml = b_amp_mg_l * 1000.0
    c0_ng_ml = c0_total_mg_l * 1000.0

    analytical_auc_inf_ng_h_ml = (dose_mg_kg / cl_l_h_kg) * 1000.0
    half_life_alpha = math.log(2.0) / alpha
    half_life_beta = math.log(2.0) / beta

    dt = max(t_end_h, 5 * half_life_beta) / max(num_points, 100)
    time_series = []
    for i in range(num_points + 1):
        t = i * dt
        conc = a_amp_ng_ml * math.exp(-alpha * t) + b_amp_ng_ml * math.exp(-beta * t)
        time_series.append({"time": round(t, 4), "concentration": round(conc, 4), "unit": "ng/mL"})

    auc_last_num = 0.0
    for i in range(len(time_series) - 1):
        t1, c1 = time_series[i]["time"], time_series[i]["concentration"]
        t2, c2 = time_series[i + 1]["time"], time_series[i + 1]["concentration"]
        auc_last_num += 0.5 * (c1 + c2) * (t2 - t1)

    c_last = time_series[-1]["concentration"]
    auc_inf_num = auc_last_num + (c_last / beta if beta > 0 else 0.0)
    agreement_pct = round(min(auc_inf_num, analytical_auc_inf_ng_h_ml) / max(auc_inf_num, analytical_auc_inf_ng_h_ml) * 100.0, 2)

    return {
        "k10": round(k10, 6),
        "k12": round(k12, 6),
        "k21": round(k21, 6),
        "alpha": round(alpha, 6),
        "beta": round(beta, 6),
        "half_life_alpha_hours": round(half_life_alpha, 4),
        "half_life_beta_hours": round(half_life_beta, 4),
        "a_amp_ng_ml": round(a_amp_ng_ml, 4),
        "b_amp_ng_ml": round(b_amp_ng_ml, 4),
        "c0_ng_ml": round(c0_ng_ml, 4),
        "cmax_ng_ml": round(c0_ng_ml, 4),
        "tmax_hours": 0.0,
        "auc_last_ng_h_ml": round(auc_last_num, 2),
        "auc_inf_analytical_ng_h_ml": round(analytical_auc_inf_ng_h_ml, 2),
        "auc_inf_numerical_ng_h_ml": round(auc_inf_num, 2),
        "auc_agreement_pct": agreement_pct,
        "time_series": time_series,
    }


def fit_two_compartment_experimental(observations: list[PKObservation], dose_mg_kg: float) -> dict[str, Any]:
    """Fit a biexponential C(t) = A * exp(-alpha * t) + B * exp(-beta * t) to dense experimental IV points.

    Minimum data requirement: >= 4 non-BLQ points.
    """
    valid_obs = [o for o in observations if not o.blq_flag and o.concentration_normalized_ng_ml is not None and o.concentration_normalized_ng_ml > 0]
    if len(valid_obs) < 4:
        return {
            "status": "MODEL_UNAVAILABLE",
            "message": f"Experimental 2-compartment fitting requires at least 4 non-BLQ points; got {len(valid_obs)}.",
        }

    t_arr = np.array([o.time_hours for o in valid_obs], dtype=float)
    c_arr = np.array([o.concentration_normalized_ng_ml for o in valid_obs], dtype=float)

    # Initial parameter guesses
    c_max = float(np.max(c_arr))
    p0 = [c_max * 0.7, 1.0, c_max * 0.3, 0.1]
    bounds = ([0.0, 0.001, 0.0, 0.0001], [c_max * 10, 100.0, c_max * 10, 50.0])

    def func(t, a, alpha, b, beta):
        return a * np.exp(-alpha * t) + b * np.exp(-beta * t)

    try:
        popt, pcov = curve_fit(func, t_arr, c_arr, p0=p0, bounds=bounds, maxfev=5000)
        a, alpha, b, beta = popt
        if alpha < beta:
            a, b = b, a
            alpha, beta = beta, alpha

        rss = float(np.sum((c_arr - func(t_arr, *popt)) ** 2))
        n = len(c_arr)
        k = 4
        rmse = math.sqrt(rss / n)
        aic = n * math.log(rss / n) + 2 * k + (2 * k * (k + 1)) / max(n - k - 1, 1)

        # Microconstants recovery
        c0_ng_ml = a + b
        c0_mg_l = c0_ng_ml / 1000.0
        vc_l_kg = dose_mg_kg / c0_mg_l if c0_mg_l > 0 else 0.0
        auc_inf = (a / alpha) + (b / beta)
        cl_l_h_kg = (dose_mg_kg / (auc_inf / 1000.0)) if auc_inf > 0 else 0.0
        k21 = (a * beta + b * alpha) / (a + b) if (a + b) > 0 else 0.0
        k10 = (alpha * beta) / k21 if k21 > 0 else 0.0
        k12 = alpha + beta - k10 - k21
        q_l_h_kg = k12 * vc_l_kg
        vp_l_kg = q_l_h_kg / k21 if k21 > 0 else 0.0

        return {
            "status": "FIT_SUCCESS",
            "a_amp_ng_ml": round(a, 4),
            "alpha": round(alpha, 6),
            "b_amp_ng_ml": round(b, 4),
            "beta": round(beta, 6),
            "rss": round(rss, 4),
            "rmse": round(rmse, 4),
            "aic": round(aic, 2),
            "cl_l_h_kg": round(cl_l_h_kg, 4),
            "vc_l_kg": round(vc_l_kg, 4),
            "q_l_h_kg": round(q_l_h_kg, 4),
            "vp_l_kg": round(vp_l_kg, 4),
            "k10": round(k10, 6),
            "k12": round(k12, 6),
            "k21": round(k21, 6),
            "n_points": n,
        }
    except Exception as exc:
        return {
            "status": "FIT_FAILED",
            "message": f"Nonlinear curve fitting did not converge: {exc}",
        }


def compute_goodness_of_fit(
    observations: list[PKObservation],
    time_series: list[dict[str, Any]],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Calculate RMSE, MAE, Fold Errors, and Residual table comparing experimental observations with simulated values."""
    residuals_list = []
    sq_errors = []
    abs_errors = []

    for obs in observations:
        t_obs = obs.time_hours
        c_obs = obs.concentration_normalized_ng_ml

        if obs.blq_flag or c_obs is None:
            residuals_list.append({
                "time_hours": t_obs,
                "observed_ng_ml": "BLQ",
                "simulated_ng_ml": None,
                "residual_ng_ml": None,
                "fold_error": None,
                "status": "BLQ_EXCLUDED",
            })
            continue

        # Find closest simulated point
        sim_match = min(time_series, key=lambda p: abs(p["time"] - t_obs))
        c_sim = sim_match["concentration"]
        res = c_obs - c_sim
        fold_err = round(c_sim / c_obs, 3) if c_obs > 0 else None

        sq_errors.append(res**2)
        abs_errors.append(abs(res))

        residuals_list.append({
            "time_hours": t_obs,
            "observed_ng_ml": round(c_obs, 4),
            "simulated_ng_ml": round(c_sim, 4),
            "residual_ng_ml": round(res, 4),
            "fold_error": fold_err,
            "status": "VALID",
        })

    n = len(sq_errors)
    if n > 0:
        rmse = math.sqrt(sum(sq_errors) / n)
        mae = sum(abs_errors) / n
        rss = sum(sq_errors)
    else:
        rmse = None
        mae = None
        rss = None

    metrics = {
        "n_points_compared": n,
        "rss": round(rss, 4) if rss is not None else None,
        "rmse_ng_ml": round(rmse, 4) if rmse is not None else None,
        "mae_ng_ml": round(mae, 4) if mae is not None else None,
    }
    return metrics, residuals_list


# Master Simulation Orchestrator
def run_pk_simulation(
    db: Session,
    version_id: int,
    request: PKSimulationRequest,
) -> PKSimulationRun:
    version = db.get(CompoundVersion, version_id)
    if not version:
        raise HTTPException(status_code=404, detail=f"CompoundVersion #{version_id} not found.")

    compound = version.compound
    project_id = compound.project_id
    species = request.species

    # 1. Fetch Stage 5A PK Parameter Foundation & Experimental IV Studies
    foundation = get_pk_foundation_profile(db, version_id, species)
    routes = foundation.get("route_parameter_sets", {})
    iv_set = routes.get("IV", {})
    dist = foundation.get("distribution", {})

    # Parameter selection logic with strict priority
    warnings = []
    parameter_sources = {}

    # Clearance selection
    if request.user_cl_override is not None and request.user_cl_override > 0:
        cl_val = request.user_cl_override
        cl_unit = "mL/min/kg"
        cl_conf = "HIGH"
        parameter_sources["CL"] = {"source": "USER_OVERRIDE", "type": "User Specified", "confidence": "HIGH"}
    elif iv_set.get("cl_value") is not None:
        cl_val = iv_set["cl_value"]
        cl_unit = iv_set.get("cl_unit", "mL/min/kg")
        cl_conf = iv_set.get("confidence", "HIGH")
        parameter_sources["CL"] = {
            "source": iv_set.get("cl_source_type", "EXPERIMENTAL_NCA"),
            "type": "Experimental IV Systemic Clearance",
            "confidence": cl_conf,
        }
    else:
        raise HTTPException(
            status_code=400,
            detail=f"Clearance (CL) is unavailable for species {species}. Provide an experimental IV PK study or IVIVE run, or specify user_cl_override.",
        )

    if iv_set.get("cl_source_type") == "PREDICTED_HEPATIC_IVIVE":
        warnings.append("HEPATIC-CLEARANCE-ONLY APPROXIMATION: Simulation uses predicted hepatic CLh as systemic CL fallback.")
        cl_conf = "MEDIUM" if cl_conf == "HIGH" else cl_conf

    # Volume selection
    if request.user_v_override is not None and request.user_v_override > 0:
        v_val = request.user_v_override
        v_unit = "L/kg"
        v_type = "User Specified"
        v_conf = "HIGH"
        parameter_sources["V"] = {"source": "USER_OVERRIDE", "type": "User Specified", "confidence": "HIGH"}
    elif dist.get("v_value") is not None:
        v_val = dist["v_value"]
        v_unit = dist.get("v_unit", "L/kg")
        v_type = dist.get("v_type", "Estimated Vd")
        v_conf = "HIGH" if v_type in ("Vss", "Vz") else "MEDIUM"
        parameter_sources["V"] = {
            "source": "EXPERIMENTAL_NCA" if v_type in ("Vss", "Vz") else "ESTIMATED_VD",
            "type": v_type,
            "confidence": v_conf,
        }
    else:
        raise HTTPException(
            status_code=400,
            detail=f"Volume of distribution (V) is unavailable for species {species}. Provide an experimental IV study or binding/cLogP data, or specify user_v_override.",
        )

    # Convert to canonical units
    try:
        canon = canonicalize_units(
            dose=request.dose,
            dose_unit=request.dose_unit,
            cl=cl_val,
            cl_unit=cl_unit,
            v=v_val,
            v_unit=v_unit,
        )
    except ValueError as err:
        raise HTTPException(status_code=422, detail=f"Dimensional check error: {err}")

    # Confidence ceiling
    confidence_levels = {"HIGH": 3, "MEDIUM": 2, "LOW": 1, "MODEL_UNAVAILABLE": 0}
    rev_levels = {3: "HIGH", 2: "MEDIUM", 1: "LOW", 0: "MODEL_UNAVAILABLE"}
    min_conf_val = min(confidence_levels.get(cl_conf, 1), confidence_levels.get(v_conf, 1))
    sim_confidence = rev_levels[min_conf_val]

    # Fetch matching experimental IV observations for overlay & residuals
    iv_studies = db.scalars(
        select(PKStudy).where(
            PKStudy.compound_row_id == version.compound_row_id,
            PKStudy.version_id == version_id,
            PKStudy.species == species,
            PKStudy.route == "IV",
        )
    ).all()
    exp_obs = []
    if iv_studies:
        exp_obs = list(iv_studies[0].observations)

    # Perform mathematical simulation according to model choice
    if request.model_type == "TWO_COMPARTMENT":
        # Check if experimental 2-compartment fit is supported
        fit_res = fit_two_compartment_experimental(exp_obs, canon["dose_mg_kg"]) if exp_obs else {"status": "MODEL_UNAVAILABLE"}
        if fit_res.get("status") == "FIT_SUCCESS":
            sim_res = simulate_two_compartment_iv_bolus(
                dose_mg_kg=canon["dose_mg_kg"],
                cl_l_h_kg=fit_res["cl_l_h_kg"],
                vc_l_kg=fit_res["vc_l_kg"],
                q_l_h_kg=fit_res["q_l_h_kg"],
                vp_l_kg=fit_res["vp_l_kg"],
                t_end_h=request.custom_t_end or 24.0,
            )
            sim_res["fit_details"] = fit_res
            parameter_sources["Model"] = {"type": "Experimental 2-Compartment Fit", "method": "SciPy Biexponential LSQ"}
        else:
            raise HTTPException(
                status_code=400,
                detail=f"2-Compartment IV model is unavailable: {fit_res.get('message', 'Sufficient microconstants (Vc, Q, Vp) unavailable without dense experimental data.')}. Never fabricate 2-comp parameters from Vz alone.",
            )
    else:
        # 1-Compartment Model
        if request.administration_type == "IV_INFUSION":
            sim_res = simulate_one_compartment_iv_infusion(
                dose_mg_kg=canon["dose_mg_kg"],
                infusion_duration_h=request.infusion_duration_hours,
                cl_l_h_kg=canon["cl_l_h_kg"],
                v_l_kg=canon["v_l_kg"],
                num_doses=request.num_doses if request.dosing_frequency == "Repeated Dosing" else 1,
                dose_interval_h=request.dose_interval_hours,
                t_end_h=request.custom_t_end or 24.0,
            )
            parameter_sources["Model"] = {"type": "1-Compartment IV Infusion", "method": "Analytical Linear PK"}
        else:
            sim_res = simulate_one_compartment_iv_bolus(
                dose_mg_kg=canon["dose_mg_kg"],
                cl_l_h_kg=canon["cl_l_h_kg"],
                v_l_kg=canon["v_l_kg"],
                num_doses=request.num_doses if request.dosing_frequency == "Repeated Dosing" else 1,
                dose_interval_h=request.dose_interval_hours,
                t_end_h=request.custom_t_end or 24.0,
            )
            parameter_sources["Model"] = {"type": "1-Compartment IV Bolus", "method": "Analytical Linear PK"}

    if request.dosing_frequency == "Repeated Dosing":
        warnings.append("LINEAR PK ASSUMPTION: Multi-dose superposition assumes linear kinetics without accumulation saturation.")

    # Calculate goodness-of-fit & residuals if experimental observations exist
    gof_metrics, residual_table = compute_goodness_of_fit(exp_obs, sim_res["time_series"]) if exp_obs else ({}, [])

    # Assemble Parameter Snapshot
    param_snapshot = {
        "dose_normalized_mg_kg": canon["dose_mg_kg"],
        "cl_l_h_kg": canon["cl_l_h_kg"],
        "cl_ml_min_kg": round(canon["cl_l_h_kg"] * 1000.0 / 60.0, 2),
        "v_l_kg": canon["v_l_kg"],
        "v_type": v_type,
        "k_elim_h": sim_res["k_elim"],
        "half_life_hours": sim_res["half_life_hours"],
    }

    # Assemble Output Metrics
    output_metrics = {
        "c0_ng_ml": sim_res.get("c0_ng_ml"),
        "cmax_ng_ml": sim_res["cmax_ng_ml"],
        "tmax_hours": sim_res["tmax_hours"],
        "auc_last_ng_h_ml": sim_res["auc_last_ng_h_ml"],
        "auc_inf_analytical_ng_h_ml": sim_res.get("auc_inf_analytical_ng_h_ml"),
        "auc_inf_numerical_ng_h_ml": sim_res.get("auc_inf_numerical_ng_h_ml"),
        "auc_agreement_pct": sim_res.get("auc_agreement_pct"),
        "half_life_hours": sim_res["half_life_hours"],
        "goodness_of_fit": gof_metrics,
        "uncertainty_status": "UNCERTAINTY NOT QUANTIFIED",
    }

    provenance = {
        "engine_name": SIMULATION_ENGINE_NAME,
        "engine_version": SIMULATION_ENGINE_VERSION,
        "formula": "Analytical Linear PK Integration",
        "units": {"dose": "mg/kg", "clearance": "L/h/kg", "volume": "L/kg", "concentration": "ng/mL", "time": "hours"},
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    run_record = PKSimulationRun(
        project_id=project_id,
        compound_row_id=version.compound_row_id,
        version_id=version_id,
        species=species,
        route="IV",
        administration_type=request.administration_type,
        dose=request.dose,
        dose_unit=request.dose_unit,
        infusion_duration_hours=request.infusion_duration_hours,
        dosing_frequency=request.dosing_frequency,
        dose_interval_hours=request.dose_interval_hours,
        num_doses=request.num_doses,
        model_type=request.model_type,
        parameter_snapshot=param_snapshot,
        parameter_sources=parameter_sources,
        simulation_settings={"custom_t_end": request.custom_t_end, "num_points": len(sim_res["time_series"])},
        output_metrics=output_metrics,
        time_series=sim_res["time_series"],
        residuals=residual_table,
        warnings=warnings,
        confidence=sim_confidence,
        provenance=provenance,
    )
    db.add(run_record)
    db.commit()
    db.refresh(run_record)
    return run_record


def register_simulation_routes(app):
    @app.get("/api/compound-versions/{version_id}/pk-simulation/preview")
    def preview_simulation(
        version_id: int,
        species: str = Query("Rat"),
        db: Session = Depends(get_db),
    ):
        foundation = get_pk_foundation_profile(db, version_id, species)
        routes = foundation.get("route_parameter_sets", {})
        iv_set = routes.get("IV", {})
        dist = foundation.get("distribution", {})

        cl_val = iv_set.get("cl_value")
        cl_source = iv_set.get("cl_source_type", "UNAVAILABLE")
        v_val = dist.get("v_value")
        v_type = dist.get("v_type", "UNAVAILABLE")

        warnings = []
        if cl_source == "PREDICTED_HEPATIC_IVIVE":
            warnings.append("HEPATIC-CLEARANCE-ONLY APPROXIMATION: Simulation will use predicted hepatic CLh as systemic CL fallback.")

        return {
            "version_id": version_id,
            "species": species,
            "route": "IV",
            "clearance": {"value": cl_val, "unit": "mL/min/kg", "source": cl_source},
            "volume": {"value": v_val, "unit": "L/kg", "type": v_type},
            "available_models": [
                {"key": "ONE_COMPARTMENT_BOLUS", "name": "1-Compartment IV Bolus", "status": "AVAILABLE" if cl_val and v_val else "MODEL_UNAVAILABLE"},
                {"key": "ONE_COMPARTMENT_INFUSION", "name": "1-Compartment IV Infusion", "status": "AVAILABLE" if cl_val and v_val else "MODEL_UNAVAILABLE"},
                {"key": "TWO_COMPARTMENT", "name": "2-Compartment IV", "status": "REQUIRES_MICROCONSTANTS_OR_DENSE_DATA"},
            ],
            "warnings": warnings,
            "confidence_ceiling": iv_set.get("confidence", "MEDIUM"),
        }

    @app.post("/api/compound-versions/{version_id}/pk-simulation/run")
    def run_simulation_endpoint(
        version_id: int,
        request: PKSimulationRequest,
        db: Session = Depends(get_db),
    ):
        run = run_pk_simulation(db, version_id, request)
        return run

    @app.get("/api/pk-simulation-runs/{run_id}")
    def get_simulation_run(
        run_id: int,
        db: Session = Depends(get_db),
    ):
        run = db.get(PKSimulationRun, run_id)
        if not run:
            raise HTTPException(status_code=404, detail=f"PKSimulationRun #{run_id} not found.")
        return run

    @app.get("/api/compound-versions/{version_id}/pk-simulation/history")
    def list_simulation_history(
        version_id: int,
        species: str | None = Query(None),
        db: Session = Depends(get_db),
    ):
        stmt = select(PKSimulationRun).where(PKSimulationRun.version_id == version_id)
        if species:
            stmt = stmt.where(PKSimulationRun.species == species)
        stmt = stmt.order_by(PKSimulationRun.id.desc())
        runs = db.scalars(stmt).all()
        return runs

    @app.delete("/api/pk-simulation-runs/{run_id}")
    def delete_simulation_run(
        run_id: int,
        db: Session = Depends(get_db),
    ):
        run = db.get(PKSimulationRun, run_id)
        if not run:
            raise HTTPException(status_code=404, detail=f"PKSimulationRun #{run_id} not found.")
        db.delete(run)
        db.commit()
        return {"status": "DELETED", "id": run_id}
