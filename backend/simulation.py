"""Stage 5B-2: Extravascular & IV Pharmacokinetic Concentration-Time Simulation Engine.

Scientific Framework & Governance:
1. One-Compartment Extravascular First-Order Absorption & Elimination:
   Single Dose:
     C(t) = (F * Dose * ka) / [V * (ka - ke)] * (exp(-ke * t) - exp(-ka * t))
     where ke = CL / V, F is systemic bioavailability fraction, ka is absorption rate constant.
   Singularity handling (ka == ke):
     C(t) = (F * Dose / V) * ke * t * exp(-ke * t)
2. Route-Specific Isolation & Modeling:
   - IV: Bolus (C0 = Dose/V) and Infusion (R0/CL during infusion, post-infusion decay).
   - PO: Oral absorption with gut lumen Fa, intestinal Fg, hepatic Fh decomposition.
   - SC: Subcutaneous absorption with route-specific F_sc and ka_sc (no intestinal Fg applied).
   - IP: Intraperitoneal absorption with route-specific F_ip and ka_ip.
3. Bioavailability Hierarchy:
   - Matched experimental absolute F (exact CompoundVersion, species, route)
   > Experimentally supported component assembly
   > Mechanistic/calculated component assembly (Fa * Fg * Fh)
   > Validated prediction
   > MODEL_UNAVAILABLE.
   If experimental F exists, simulation proceeds using experimental F while preserving
   and displaying mechanistic Fa/Fg/Fh components without overwriting.
4. Absorption Rate Constant (ka) Hierarchy & Tmax Solver:
   - Priority: Experimental fit > Tmax-derived numerical solver > Validated route model > MODEL_UNAVAILABLE.
   - For 1-compartment 1st-order: Tmax = ln(ka/ke) / (ka - ke).
     Numerical root solving via Brent's method when experimental Tmax and ke are known.
     Returns KA_ESTIMATION_UNRELIABLE if solution is non-identifiable or non-positive.
5. Flip-Flop Kinetics Detection:
   - If ka <= ke: Flag POTENTIAL FLIP-FLOP KINETICS.
     Explains that apparent terminal slope reflects absorption rather than systemic elimination.
6. Repeated Extravascular Dosing:
   - Linear multi-dose superposition across N doses with interval tau.
   - Computes accumulation ratio R_acc = 1 / (1 - exp(-ke * tau)), Css,avg = (F * Dose) / (CL * tau),
     peak/trough concentrations, with explicit LINEAR PK ASSUMPTION warning.
7. Parameter Fitting & Observation Overlay:
   - Fits ka (and F if unconstrained) from dense experimental observations with fixed IV CL and V.
   - Computes RMSE, MAE, Residuals, Fold Error with BLQ exclusion.
"""

from __future__ import annotations

import math
from datetime import datetime, timezone
from typing import Any

import numpy as np
from fastapi import Depends, HTTPException, Query
from pydantic import BaseModel
from scipy.optimize import curve_fit, root_scalar
from sqlalchemy import JSON, Boolean, DateTime, Float, ForeignKey, Integer, String, Text, inspect, select, text
from sqlalchemy.orm import Mapped, mapped_column, relationship, Session

from .database import Base, get_db
from .ivive import PKParameterSet, estimate_absorption_components, estimate_volume_of_distribution, get_pk_foundation_profile
from .models import Compound, CompoundVersion, Project, utcnow
from .pk import PKObservation, PKStudy, calculate_bioavailability_for_version

SIMULATION_ENGINE_NAME = "Stage 5B-2 Extravascular & IV PK Simulation Engine"
SIMULATION_ENGINE_VERSION = "5B-2.0"


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
    # Check for newly added columns if table already exists
    with engine.connect() as conn:
        existing_cols = {col["name"] for col in inspector.get_columns("pk_simulation_runs")}
        if "f_value" not in existing_cols:
            conn.execute(text("ALTER TABLE pk_simulation_runs ADD COLUMN f_value FLOAT"))
        if "f_source" not in existing_cols:
            conn.execute(text("ALTER TABLE pk_simulation_runs ADD COLUMN f_source VARCHAR(100)"))
        if "ka_value" not in existing_cols:
            conn.execute(text("ALTER TABLE pk_simulation_runs ADD COLUMN ka_value FLOAT"))
        if "ka_source" not in existing_cols:
            conn.execute(text("ALTER TABLE pk_simulation_runs ADD COLUMN ka_source VARCHAR(100)"))
        if "flip_flop_flag" not in existing_cols:
            conn.execute(text("ALTER TABLE pk_simulation_runs ADD COLUMN flip_flop_flag BOOLEAN DEFAULT 0"))
        if "absorption_components" not in existing_cols:
            conn.execute(text("ALTER TABLE pk_simulation_runs ADD COLUMN absorption_components JSON DEFAULT '{}'"))
        if "steady_state_metrics" not in existing_cols:
            conn.execute(text("ALTER TABLE pk_simulation_runs ADD COLUMN steady_state_metrics JSON DEFAULT '{}'"))
        conn.commit()


class PKSimulationRun(Base):
    __tablename__ = "pk_simulation_runs"
    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    compound_row_id: Mapped[int] = mapped_column(ForeignKey("compounds.id", ondelete="CASCADE"), index=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    species: Mapped[str] = mapped_column(String(100), default="Rat", index=True)
    route: Mapped[str] = mapped_column(String(40), default="IV", index=True)
    administration_type: Mapped[str] = mapped_column(String(60), default="IV_BOLUS", index=True)
    dose: Mapped[float] = mapped_column(Float, default=5.0)
    dose_unit: Mapped[str] = mapped_column(String(40), default="mg/kg")
    infusion_duration_hours: Mapped[float] = mapped_column(Float, default=0.0)
    dosing_frequency: Mapped[str] = mapped_column(String(60), default="Single Dose")
    dose_interval_hours: Mapped[float] = mapped_column(Float, default=24.0)
    num_doses: Mapped[int] = mapped_column(Integer, default=1)
    model_type: Mapped[str] = mapped_column(String(60), default="ONE_COMPARTMENT", index=True)
    f_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    f_source: Mapped[str | None] = mapped_column(String(100), nullable=True)
    ka_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    ka_source: Mapped[str | None] = mapped_column(String(100), nullable=True)
    flip_flop_flag: Mapped[bool] = mapped_column(Boolean, default=False)
    parameter_snapshot: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    parameter_sources: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    absorption_components: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    steady_state_metrics: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
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
    route: str = "PO"  # "IV", "PO", "SC", "IP"
    administration_type: str = "EXTRAVASCULAR_1COMP"  # "IV_BOLUS", "IV_INFUSION", "EXTRAVASCULAR_1COMP"
    dose: float = 5.0
    dose_unit: str = "mg/kg"  # "mg/kg" or "µg/kg"
    infusion_duration_hours: float = 0.0
    dosing_frequency: str = "Single Dose"  # "Single Dose" or "Repeated Dosing"
    dose_interval_hours: float = 24.0
    num_doses: int = 1
    model_type: str = "ONE_COMPARTMENT"  # "ONE_COMPARTMENT" or "TWO_COMPARTMENT"
    user_cl_override: float | None = None
    user_v_override: float | None = None
    user_f_override: float | None = None  # Fraction (0-1) or percentage (0-100)
    user_ka_override: float | None = None  # 1/h
    custom_t_end: float | None = None


class PKFitRequest(BaseModel):
    species: str = "Rat"
    route: str = "PO"
    dose: float = 5.0
    dose_unit: str = "mg/kg"
    fix_cl_v: bool = True
    user_cl_override: float | None = None
    user_v_override: float | None = None


# Core Calculations & Dimensional Checks
def canonicalize_units(dose: float, dose_unit: str, cl: float, cl_unit: str, v: float, v_unit: str) -> dict[str, float]:
    """Convert input parameters into canonical SI-derived units:
    Dose: mg/kg
    Clearance: L/h/kg
    Volume: L/kg
    Concentration: 1 mg/L = 1000 ng/mL.
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


def solve_ka_from_tmax(tmax_obs: float, ke: float) -> dict[str, Any]:
    """Numerically solve for absorption rate constant ka from observed Tmax and known elimination rate ke.

    1-compartment 1st-order analytical relation:
      Tmax = ln(ka / ke) / (ka - ke)

    Scientific monotonicity:
      g(ka) = ln(ka/ke)/(ka-ke) is strictly monotonically decreasing on (0, inf).
      - If tmax_obs < 1/ke: ka > ke (normal absorption)
      - If tmax_obs == 1/ke: ka = ke
      - If tmax_obs > 1/ke: ka < ke (potential flip-flop kinetics)

    Returns dictionary with status CONVERGED or KA_ESTIMATION_UNRELIABLE.
    """
    if tmax_obs <= 0 or ke <= 0:
        return {
            "status": "KA_ESTIMATION_UNRELIABLE",
            "message": f"Non-identifiable parameters: Tmax ({tmax_obs}) and ke ({ke}) must be strictly positive.",
            "ka": None,
        }

    # At ka == ke, limit is 1/ke
    tmax_limit = 1.0 / ke
    if abs(tmax_obs - tmax_limit) < 1e-4:
        return {
            "status": "CONVERGED",
            "ka": round(ke, 6),
            "method": "Analytical Limit",
            "target_tmax": tmax_obs,
            "ke": ke,
            "flip_flop": True,
            "iterations": 1,
        }

    def func(ka_val):
        if ka_val <= 0:
            return 1e9
        if abs(ka_val - ke) < 1e-7:
            return (1.0 / ke) - tmax_obs
        return (math.log(ka_val / ke) / (ka_val - ke)) - tmax_obs

    # Establish root brackets
    if tmax_obs < tmax_limit:
        # ka > ke
        a = ke + 1e-5
        b = max(10.0 * ke, 20.0)
        # Expand upper bracket if needed
        for _ in range(15):
            if func(b) < 0:
                break
            b *= 4.0
    else:
        # ka < ke (flip-flop)
        a = min(0.0001, ke / 1000.0)
        b = ke - 1e-5
        # Shrink lower bracket if needed
        for _ in range(15):
            if func(a) > 0:
                break
            a /= 4.0

    try:
        sol = root_scalar(func, bracket=[a, b], method="brentq", xtol=1e-6, maxiter=200)
        if sol.converged and sol.root > 0:
            solved_ka = float(sol.root)
            return {
                "status": "CONVERGED",
                "ka": round(solved_ka, 6),
                "method": "BrentQ",
                "target_tmax": tmax_obs,
                "ke": ke,
                "flip_flop": bool(solved_ka <= ke),
                "iterations": sol.iterations,
            }
    except Exception as exc:
        pass

    return {
        "status": "KA_ESTIMATION_UNRELIABLE",
        "message": f"Numerical solver could not identify unique stable ka for Tmax={tmax_obs}h and ke={ke:.4f}/h.",
        "ka": None,
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
    """Single or repeated IV bolus 1-compartment linear PK simulation."""
    k_elim = cl_l_h_kg / v_l_kg  # 1/h
    half_life_h = math.log(2.0) / k_elim
    c0_ng_ml = (dose_mg_kg / v_l_kg) * 1000.0
    analytical_auc_inf_ng_h_ml = (dose_mg_kg / cl_l_h_kg) * 1000.0

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

    # Multi-dose steady-state metrics
    steady_state = {}
    if num_doses > 1:
        r_acc = 1.0 / (1.0 - math.exp(-k_elim * dose_interval_h))
        css_avg = (dose_mg_kg / (cl_l_h_kg * dose_interval_h)) * 1000.0
        steady_state = {
            "accumulation_ratio": round(r_acc, 3),
            "css_avg_ng_ml": round(css_avg, 2),
            "tau_hours": dose_interval_h,
            "doses_administered": num_doses,
        }

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
        "steady_state": steady_state,
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
    """Single or repeated IV infusion 1-compartment linear PK simulation."""
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
    """Forward IV bolus 2-compartment linear PK simulation given microconstants."""
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


def simulate_one_compartment_extravascular(
    dose_mg_kg: float,
    cl_l_h_kg: float,
    v_l_kg: float,
    f_fraction: float,
    ka_h: float,
    num_doses: int = 1,
    dose_interval_h: float = 24.0,
    t_end_h: float = 24.0,
    num_points: int = 200,
) -> dict[str, Any]:
    """Single or repeated extravascular (PO/SC/IP) 1-compartment linear PK simulation with 1st-order absorption.

    C(t) = (F * Dose * ka) / [V * (ka - ke)] * (exp(-ke*t) - exp(-ka*t))  [mg/L = 1000 ng/mL]
    Singularity handled when ka == ke via C(t) = (F * Dose / V) * ke * t * exp(-ke * t).
    """
    if f_fraction <= 0 or f_fraction > 1.5:
        raise ValueError(f"Bioavailability fraction F must be positive and realistic (<=1.5), got {f_fraction}")
    if ka_h <= 0:
        raise ValueError(f"Absorption rate constant ka must be strictly positive (>0), got {ka_h}")

    ke = cl_l_h_kg / v_l_kg  # elimination rate constant (1/h)
    t_half_elim = math.log(2.0) / ke
    t_half_abs = math.log(2.0) / ka_h
    is_flip_flop = bool(ka_h <= ke)

    # Analytical Single-Dose Metrics
    analytical_auc_inf = (f_fraction * dose_mg_kg / cl_l_h_kg) * 1000.0  # ng*h/mL

    if abs(ka_h - ke) > 1e-6:
        analytical_tmax = math.log(ka_h / ke) / (ka_h - ke)
        c_scale = (f_fraction * dose_mg_kg * ka_h / (v_l_kg * (ka_h - ke))) * 1000.0
        analytical_cmax = c_scale * (math.exp(-ke * analytical_tmax) - math.exp(-ka_h * analytical_tmax))
    else:
        analytical_tmax = 1.0 / ke
        analytical_cmax = (f_fraction * dose_mg_kg / v_l_kg) * 1000.0 * (1.0 / math.e)

    # Time grid generation
    effective_t_half = max(t_half_elim, t_half_abs)
    total_time_h = max(t_end_h, num_doses * dose_interval_h + 5 * effective_t_half)
    dt = total_time_h / max(num_points, 100)

    # Ensure key peak points and dosing intervals are present in the time grid
    time_points = [i * dt for i in range(num_points + 1)]
    if analytical_tmax <= total_time_h:
        time_points.append(analytical_tmax)
    for d in range(num_doses):
        t_d = d * dose_interval_h
        if t_d <= total_time_h:
            time_points.append(t_d)
            if t_d + analytical_tmax <= total_time_h:
                time_points.append(t_d + analytical_tmax)

    time_grid = sorted(set(round(t, 4) for t in time_points if 0.0 <= t <= total_time_h))

    # Calculate Concentration Curve via Linear Superposition
    time_series = []
    for t in time_grid:
        conc_ng_ml = 0.0
        for d in range(num_doses):
            t_start = d * dose_interval_h
            t_rel = t - t_start
            if t_rel > 0:
                if abs(ka_h - ke) > 1e-6:
                    conc_ng_ml += (f_fraction * dose_mg_kg * ka_h / (v_l_kg * (ka_h - ke))) * (
                        math.exp(-ke * t_rel) - math.exp(-ka_h * t_rel)
                    ) * 1000.0
                else:
                    conc_ng_ml += (f_fraction * dose_mg_kg / v_l_kg) * ke * t_rel * math.exp(-ke * t_rel) * 1000.0
        time_series.append({"time": round(t, 4), "concentration": round(max(0.0, conc_ng_ml), 4), "unit": "ng/mL"})

    # Trapezoidal Numerical AUC
    auc_last_num = 0.0
    for i in range(len(time_series) - 1):
        t1, c1 = time_series[i]["time"], time_series[i]["concentration"]
        t2, c2 = time_series[i + 1]["time"], time_series[i + 1]["concentration"]
        auc_last_num += 0.5 * (c1 + c2) * (t2 - t1)

    c_last = time_series[-1]["concentration"]
    # Tail extrapolation using limiting terminal slope (ke or ka)
    terminal_k = min(ke, ka_h)
    auc_inf_num = auc_last_num + (c_last / terminal_k if terminal_k > 0 else 0.0)
    agreement_pct = round(min(auc_inf_num, analytical_auc_inf) / max(auc_inf_num, analytical_auc_inf) * 100.0, 2) if analytical_auc_inf > 0 else 100.0

    # Max concentration & Tmax from curve
    if num_doses == 1:
        cmax_val = analytical_cmax
        tmax_val = analytical_tmax
    else:
        cmax_val = max(p["concentration"] for p in time_series)
        tmax_val = min((p["time"] for p in time_series if abs(p["concentration"] - cmax_val) < 1e-4), default=analytical_tmax)

    # Multi-dose steady-state metrics
    steady_state = {}
    if num_doses > 1:
        r_acc = 1.0 / (1.0 - math.exp(-ke * dose_interval_h))
        css_avg = (f_fraction * dose_mg_kg / (cl_l_h_kg * dose_interval_h)) * 1000.0
        steady_state = {
            "accumulation_ratio": round(r_acc, 3),
            "css_avg_ng_ml": round(css_avg, 2),
            "tau_hours": dose_interval_h,
            "doses_administered": num_doses,
        }

    return {
        "k_elim": round(ke, 6),
        "k_abs": round(ka_h, 6),
        "half_life_elim_hours": round(t_half_elim, 4),
        "half_life_abs_hours": round(t_half_abs, 4),
        "half_life_hours": round(t_half_abs if is_flip_flop else t_half_elim, 4),
        "is_flip_flop": is_flip_flop,
        "f_fraction": round(f_fraction, 4),
        "cmax_ng_ml": round(cmax_val, 4),
        "tmax_hours": round(tmax_val, 4),
        "auc_last_ng_h_ml": round(auc_last_num, 2),
        "auc_inf_analytical_ng_h_ml": round(analytical_auc_inf, 2),
        "auc_inf_numerical_ng_h_ml": round(auc_inf_num, 2),
        "auc_agreement_pct": agreement_pct,
        "steady_state": steady_state,
        "time_series": time_series,
    }


def fit_two_compartment_experimental(observations: list[PKObservation], dose_mg_kg: float) -> dict[str, Any]:
    """Fit a biexponential C(t) = A * exp(-alpha * t) + B * exp(-beta * t) to dense experimental IV points."""
    valid_obs = [o for o in observations if not o.blq_flag and o.concentration_normalized_ng_ml is not None and o.concentration_normalized_ng_ml > 0]
    if len(valid_obs) < 4:
        return {
            "status": "MODEL_UNAVAILABLE",
            "message": f"Experimental 2-compartment fitting requires at least 4 non-BLQ points; got {len(valid_obs)}.",
        }

    t_arr = np.array([o.time_hours for o in valid_obs], dtype=float)
    c_arr = np.array([o.concentration_normalized_ng_ml for o in valid_obs], dtype=float)

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


def fit_one_compartment_extravascular(
    observations: list[PKObservation],
    dose_mg_kg: float,
    cl_l_h_kg: float,
    v_l_kg: float,
    f_fixed: float | None = None,
) -> dict[str, Any]:
    """Fit absorption rate ka (and F if unconstrained) to experimental extravascular data.

    Preserves IV-derived CL and V to prevent non-identifiability.
    Requires >= 3 non-BLQ points.
    """
    valid_obs = [o for o in observations if not o.blq_flag and o.concentration_normalized_ng_ml is not None and o.concentration_normalized_ng_ml > 0]
    if len(valid_obs) < 3:
        return {
            "status": "MODEL_UNAVAILABLE",
            "message": f"Extravascular fitting requires at least 3 non-BLQ points; got {len(valid_obs)}.",
        }

    t_arr = np.array([o.time_hours for o in valid_obs], dtype=float)
    c_arr = np.array([o.concentration_normalized_ng_ml for o in valid_obs], dtype=float)
    ke = cl_l_h_kg / v_l_kg

    if f_fixed is not None and f_fixed > 0:
        # Fit ka only
        def model_func(t, ka_val):
            return np.where(
                np.abs(ka_val - ke) < 1e-6,
                (f_fixed * dose_mg_kg / v_l_kg) * ke * t * np.exp(-ke * t) * 1000.0,
                (f_fixed * dose_mg_kg * ka_val / (v_l_kg * (ka_val - ke))) * (np.exp(-ke * t) - np.exp(-ka_val * t)) * 1000.0,
            )

        p0 = [1.0]
        bounds = ([0.001], [100.0])
        try:
            popt, _ = curve_fit(model_func, t_arr, c_arr, p0=p0, bounds=bounds, maxfev=5000)
            fitted_ka = float(popt[0])
            fitted_f = f_fixed
            k_params = 1
        except Exception as exc:
            return {"status": "FIT_FAILED", "message": f"ka fitting did not converge: {exc}"}
    else:
        # Fit both ka and F
        def model_func(t, ka_val, f_val):
            return np.where(
                np.abs(ka_val - ke) < 1e-6,
                (f_val * dose_mg_kg / v_l_kg) * ke * t * np.exp(-ke * t) * 1000.0,
                (f_val * dose_mg_kg * ka_val / (v_l_kg * (ka_val - ke))) * (np.exp(-ke * t) - np.exp(-ka_val * t)) * 1000.0,
            )

        p0 = [1.0, 0.5]
        bounds = ([0.001, 0.001], [100.0, 1.2])
        try:
            popt, _ = curve_fit(model_func, t_arr, c_arr, p0=p0, bounds=bounds, maxfev=5000)
            fitted_ka = float(popt[0])
            fitted_f = float(popt[1])
            k_params = 2
        except Exception as exc:
            return {"status": "FIT_FAILED", "message": f"ka/F fitting did not converge: {exc}"}

    pred_c = model_func(t_arr, fitted_ka, *( [fitted_f] if f_fixed is None else [] ))
    rss = float(np.sum((c_arr - pred_c) ** 2))
    n = len(c_arr)
    rmse = math.sqrt(rss / n)
    safe_mse = max(rss / n, 1e-12)
    aic = n * math.log(safe_mse) + 2 * k_params + (2 * k_params * (k_params + 1)) / max(n - k_params - 1, 1)

    return {
        "status": "FIT_SUCCESS",
        "fitted_ka": round(fitted_ka, 6),
        "fitted_f": round(fitted_f, 4),
        "fitted_f_pct": round(fitted_f * 100.0, 2),
        "ke": round(ke, 6),
        "is_flip_flop": bool(fitted_ka <= ke),
        "rss": round(rss, 4),
        "rmse": round(rmse, 4),
        "aic": round(aic, 2),
        "n_points": n,
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
    admin_type = (request.administration_type or "").strip().upper()
    route = (request.route or "").strip().upper()
    if admin_type in {"IV_BOLUS", "IV_INFUSION"} and route in {"", "PO"}:
        route = "IV"
    if not route or route not in {"IV", "PO", "SC", "IP"}:
        route = "IV" if admin_type in {"IV_BOLUS", "IV_INFUSION"} else "PO"

    # 1. Fetch Stage 5A PK Foundation & Route Sets
    foundation = get_pk_foundation_profile(db, version_id, species)
    routes = foundation.get("route_parameter_sets", {})
    iv_set = routes.get("IV", {})
    route_set = routes.get(route, {})
    dist = foundation.get("distribution", {})
    abs_decomp = foundation.get("absorption", {})

    warnings: list[str] = []
    parameter_sources: dict[str, Any] = {}

    # 2. Clearance Selection (Systemic CL)
    if request.user_cl_override is not None and request.user_cl_override > 0:
        cl_val = request.user_cl_override
        cl_unit = "mL/min/kg"
        cl_conf = "HIGH"
        parameter_sources["CL"] = {"source": "USER_OVERRIDE", "type": "User Specified", "evidence_type": "USER_SPECIFIED", "confidence": "HIGH"}
    elif iv_set.get("cl_value") is not None:
        cl_val = iv_set["cl_value"]
        cl_unit = iv_set.get("cl_unit", "mL/min/kg")
        cl_conf = iv_set.get("confidence", "HIGH")
        parameter_sources["CL"] = {
            "source": iv_set.get("cl_source_type", "EXPERIMENTAL_NCA"),
            "type": "Experimental IV Systemic Clearance",
            "evidence_type": "EXPERIMENTAL" if "EXPERIMENTAL" in str(iv_set.get("cl_source_type")) else "DERIVED_ESTIMATE",
            "confidence": cl_conf,
        }
    elif route_set.get("clh_value") is not None:
        cl_val = route_set["clh_value"]
        cl_unit = "mL/min/kg"
        cl_conf = "LOW"
        parameter_sources["CL"] = {
            "source": "PREDICTED_HEPATIC_IVIVE",
            "type": "Hepatic Clearance Fallback",
            "evidence_type": "DERIVED_ESTIMATE",
            "confidence": "LOW",
        }
    else:
        raise HTTPException(
            status_code=400,
            detail=f"Clearance (CL) is unavailable for species {species}. Provide an experimental IV PK study or IVIVE run, or specify user_cl_override.",
        )

    if parameter_sources.get("CL", {}).get("source") == "PREDICTED_HEPATIC_IVIVE":
        warnings.append("HEPATIC-CLEARANCE-ONLY APPROXIMATION: Simulation uses predicted hepatic CLh as systemic CL fallback.")
        cl_conf = "MEDIUM" if cl_conf == "HIGH" else cl_conf

    # 3. Volume Selection (Systemic V)
    if request.user_v_override is not None and request.user_v_override > 0:
        v_val = request.user_v_override
        v_unit = "L/kg"
        v_type = "User Specified"
        v_conf = "HIGH"
        parameter_sources["V"] = {"source": "USER_OVERRIDE", "type": "User Specified", "evidence_type": "USER_SPECIFIED", "confidence": "HIGH"}
    elif dist.get("v_value") is not None:
        v_val = dist["v_value"]
        v_unit = dist.get("v_unit", "L/kg")
        v_type = dist.get("v_type", "Estimated Vd")
        v_conf = "HIGH" if v_type in ("Vss", "Vz") else "MEDIUM"
        parameter_sources["V"] = {
            "source": "EXPERIMENTAL_NCA" if v_type in ("Vss", "Vz") else "ESTIMATED_VD",
            "type": v_type,
            "evidence_type": "EXPERIMENTAL" if v_type in ("Vss", "Vz") else "DERIVED_ESTIMATE",
            "confidence": v_conf,
        }
    else:
        raise HTTPException(
            status_code=400,
            detail=f"Volume of distribution (V) is unavailable for species {species}. Provide an experimental IV study or binding/cLogP data, or specify user_v_override.",
        )

    # 4. Bioavailability (F) Resolution
    f_val: float | None = None
    f_source: str = "MODEL_UNAVAILABLE"
    f_evidence: str = "MODEL_UNAVAILABLE"
    f_conf: str = "MODEL_UNAVAILABLE"

    if route == "IV":
        f_val = 1.0
        f_source = "IV_ROUTE_DEFINITION"
        f_evidence = "THEORETICAL_MAXIMUM"
        f_conf = "HIGH"
        parameter_sources["F"] = {"value": 1.0, "source": f_source, "evidence_type": f_evidence, "confidence": f_conf}
    else:
        # Extravascular (PO, SC, IP)
        # Check matched experimental F directly from calculate_bioavailability_for_version if not present in route_set
        matched_exp_f = route_set.get("f_experimental")
        if matched_exp_f is None:
            ba_data = calculate_bioavailability_for_version(version_id, db)
            for b in ba_data.get("bioavailability", []):
                if b.get("species") == species and b.get("route") == route and b.get("status") == "MATCHED":
                    matched_exp_f = b.get("bioavailability_pct")
                    break

        if request.user_f_override is not None and request.user_f_override > 0:
            raw_f = request.user_f_override
            f_val = raw_f if raw_f <= 1.0 else (raw_f / 100.0)
            f_source = "USER_OVERRIDE"
            f_evidence = "USER_SPECIFIED"
            f_conf = "HIGH"
            parameter_sources["F"] = {"value": round(f_val * 100.0, 2), "source": f_source, "evidence_type": f_evidence, "confidence": f_conf}
        elif matched_exp_f is not None:
            f_val = float(matched_exp_f) / 100.0
            f_source = "MATCHED_EXPERIMENTAL_F"
            f_evidence = "EXPERIMENTAL"
            f_conf = "HIGH"
            parameter_sources["F"] = {"value": float(matched_exp_f), "source": f_source, "evidence_type": f_evidence, "confidence": f_conf}
        elif route_set.get("f_predicted") is not None:
            f_val = float(route_set["f_predicted"]) / 100.0
            f_source = "MECHANISTIC_COMPONENT_ASSEMBLY"
            f_evidence = "DERIVED_ESTIMATE"
            f_conf = "MEDIUM"
            parameter_sources["F"] = {"value": route_set["f_predicted"], "source": f_source, "evidence_type": f_evidence, "confidence": f_conf}
        elif abs_decomp.get("f_predicted") is not None:
            f_val = float(abs_decomp["f_predicted"]) / 100.0
            f_source = "MECHANISTIC_COMPONENT_ASSEMBLY"
            f_evidence = "DERIVED_ESTIMATE"
            f_conf = "MEDIUM"
            parameter_sources["F"] = {"value": abs_decomp["f_predicted"], "source": f_source, "evidence_type": f_evidence, "confidence": f_conf}
        elif abs_decomp.get("fh_value") is not None:
            f_val = float(abs_decomp["fh_value"])
            f_source = "HEPATIC_ESCAPE_FALLBACK"
            f_evidence = "DERIVED_ESTIMATE"
            f_conf = "LOW"
            parameter_sources["F"] = {"value": round(f_val * 100.0, 1), "source": f_source, "evidence_type": f_evidence, "confidence": f_conf}

        if f_val is None:
            warnings.append("MECHANISTIC F INCOMPLETE: Bioavailability could not be fully resolved.")
            raise HTTPException(
                status_code=400,
                detail=f"Bioavailability (F) is unavailable for route {route} in {species}. Provide matched IV/extravascular experimental studies, or specify user_f_override.",
            )

        # Check for incomplete mechanistic decomposition warning when experimental F is used
        if f_source == "MATCHED_EXPERIMENTAL_F":
            if abs_decomp.get("fa_value") is None or abs_decomp.get("fg_value") is None:
                warnings.append(f"MECHANISTIC F INCOMPLETE: Simulation uses matched experimental absolute bioavailability ({round(f_val*100.0, 1)}%), while mechanistic Fa/Fg decomposition remains incomplete.")

    # 5. Absorption Rate Constant (ka) Resolution
    ka_val: float | None = None
    ka_source: str = "MODEL_UNAVAILABLE"
    ka_evidence: str = "MODEL_UNAVAILABLE"
    ka_conf: str = "MODEL_UNAVAILABLE"
    ka_diagnostics: dict[str, Any] = {}

    # Query matching extravascular study for this route and species
    extra_studies = db.scalars(
        select(PKStudy).where(
            PKStudy.compound_row_id == version.compound_row_id,
            PKStudy.version_id == version_id,
            PKStudy.species == species,
            PKStudy.route == route,
        )
    ).all()
    extra_study = extra_studies[0] if extra_studies else None

    if route == "IV":
        ka_val = None
        ka_source = "NOT_APPLICABLE_FOR_IV"
        ka_evidence = "NOT_APPLICABLE"
        ka_conf = "HIGH"
    else:
        if request.user_ka_override is not None and request.user_ka_override > 0:
            ka_val = float(request.user_ka_override)
            ka_source = "USER_OVERRIDE"
            ka_evidence = "USER_SPECIFIED"
            ka_conf = "HIGH"
            parameter_sources["ka"] = {"value": ka_val, "unit": "1/h", "source": ka_source, "evidence_type": ka_evidence, "confidence": ka_conf}
        elif extra_study and extra_study.latest_nca and extra_study.latest_nca.tmax is not None:
            tmax_obs = float(extra_study.latest_nca.tmax)
            ke_est = (cl_val * 60.0 / 1000.0 if cl_unit in ("mL/min/kg", "mL/min") else cl_val) / v_val
            ka_sol = solve_ka_from_tmax(tmax_obs, ke_est)
            if ka_sol.get("status") == "CONVERGED":
                ka_val = ka_sol["ka"]
                ka_source = "EXPERIMENTAL_TMAX_DERIVED"
                ka_evidence = "DERIVED_ESTIMATE"
                ka_conf = "MEDIUM"
                ka_diagnostics = ka_sol
                parameter_sources["ka"] = {
                    "value": ka_val,
                    "unit": "1/h",
                    "source": ka_source,
                    "evidence_type": ka_evidence,
                    "confidence": ka_conf,
                    "target_tmax": tmax_obs,
                }
            else:
                warnings.append("KA ESTIMATION UNRELIABLE: Numerical estimation of ka from observed Tmax failed or is non-identifiable.")
        elif extra_study and extra_study.observations:
            # Attempt experimental fitting
            fit_res = fit_one_compartment_extravascular(list(extra_study.observations), request.dose, (cl_val * 60.0 / 1000.0 if cl_unit in ("mL/min/kg", "mL/min") else cl_val), v_val, f_fixed=f_val)
            if fit_res.get("status") == "FIT_SUCCESS":
                ka_val = fit_res["fitted_ka"]
                ka_source = "FITTED_FROM_EXPERIMENTAL_PK"
                ka_evidence = "EXPERIMENTAL"
                ka_conf = "HIGH"
                ka_diagnostics = fit_res
                parameter_sources["ka"] = {"value": ka_val, "unit": "1/h", "source": ka_source, "evidence_type": ka_evidence, "confidence": ka_conf}
        elif route_set.get("ka_value") is not None:
            ka_val = float(route_set["ka_value"])
            ka_source = route_set.get("ka_source_type", "DERIVED_FROM_PERMEABILITY")
            ka_evidence = "DERIVED_ESTIMATE"
            ka_conf = "MEDIUM"
            parameter_sources["ka"] = {"value": ka_val, "unit": "1/h", "source": ka_source, "evidence_type": ka_evidence, "confidence": ka_conf}
        elif abs_decomp.get("fa_value") is not None or f_val is not None:
            # Standard fast oral absorption rate default (1.0 1/h)
            ka_val = 1.0
            ka_source = "DERIVED_DEFAULT"
            ka_evidence = "DERIVED_ESTIMATE"
            ka_conf = "LOW"
            parameter_sources["ka"] = {"value": ka_val, "unit": "1/h", "source": ka_source, "evidence_type": ka_evidence, "confidence": ka_conf}

        if ka_val is None:
            warnings.append("KA ESTIMATION UNRELIABLE: Absorption rate constant ka could not be identified.")
            raise HTTPException(
                status_code=400,
                detail=f"Absorption rate constant (ka) is unavailable for route {route} in {species}. Provide an experimental extravascular PK study with observed Tmax, fit ka from observations, or specify user_ka_override.",
            )

    # 6. Convert to canonical SI units
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

    # 7. Confidence Ceiling
    confidence_levels = {"HIGH": 3, "MEDIUM": 2, "LOW": 1, "MODEL_UNAVAILABLE": 0}
    rev_levels = {3: "HIGH", 2: "MEDIUM", 1: "LOW", 0: "MODEL_UNAVAILABLE"}
    conf_factors = [cl_conf, v_conf]
    if route != "IV":
        conf_factors.extend([f_conf, ka_conf])
    min_conf_val = min(confidence_levels.get(c, 1) for c in conf_factors)
    sim_confidence = rev_levels[min_conf_val]

    # 8. Route-Specific Observations Fetching (Strict Isolation)
    exp_obs = []
    if route == "IV":
        iv_studies = db.scalars(
            select(PKStudy).where(
                PKStudy.compound_row_id == version.compound_row_id,
                PKStudy.version_id == version_id,
                PKStudy.species == species,
                PKStudy.route == "IV",
            )
        ).all()
        if iv_studies:
            exp_obs = list(iv_studies[0].observations)
    else:
        if extra_study and extra_study.observations:
            exp_obs = list(extra_study.observations)

    # 9. Perform Mathematical Simulation
    flip_flop_detected = False
    if route == "IV":
        if request.model_type == "TWO_COMPARTMENT":
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
                parameter_sources["Model"] = {"type": "Experimental 2-Compartment Fit", "method": "SciPy Biexponential LSQ", "evidence_type": "EXPERIMENTAL"}
            else:
                raise HTTPException(
                    status_code=400,
                    detail=f"2-Compartment IV model is unavailable: {fit_res.get('message', 'Sufficient microconstants unavailable.')}.",
                )
        else:
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
                parameter_sources["Model"] = {"type": "1-Compartment IV Infusion", "method": "Analytical Linear PK", "evidence_type": "DERIVED_ESTIMATE"}
            else:
                sim_res = simulate_one_compartment_iv_bolus(
                    dose_mg_kg=canon["dose_mg_kg"],
                    cl_l_h_kg=canon["cl_l_h_kg"],
                    v_l_kg=canon["v_l_kg"],
                    num_doses=request.num_doses if request.dosing_frequency == "Repeated Dosing" else 1,
                    dose_interval_h=request.dose_interval_hours,
                    t_end_h=request.custom_t_end or 24.0,
                )
                parameter_sources["Model"] = {"type": "1-Compartment IV Bolus", "method": "Analytical Linear PK", "evidence_type": "DERIVED_ESTIMATE"}
    else:
        # Extravascular Simulation (PO, SC, IP)
        sim_res = simulate_one_compartment_extravascular(
            dose_mg_kg=canon["dose_mg_kg"],
            cl_l_h_kg=canon["cl_l_h_kg"],
            v_l_kg=canon["v_l_kg"],
            f_fraction=f_val,
            ka_h=ka_val,
            num_doses=request.num_doses if request.dosing_frequency == "Repeated Dosing" else 1,
            dose_interval_h=request.dose_interval_hours,
            t_end_h=request.custom_t_end or 24.0,
        )
        parameter_sources["Model"] = {
            "type": f"1-Compartment {route} First-Order Absorption",
            "method": "Analytical Linear PK Integration",
            "evidence_type": "DERIVED_ESTIMATE",
        }
        flip_flop_detected = sim_res.get("is_flip_flop", False)
        if flip_flop_detected:
            ke_h = sim_res["k_elim"]
            warnings.append(
                f"POTENTIAL FLIP-FLOP KINETICS: Absorption rate constant ka ({ka_val:.4f} 1/h) <= elimination rate constant ke ({ke_h:.4f} 1/h). "
                f"The observed terminal slope reflects absorption rate-limiting kinetics rather than elimination clearance."
            )

        if route == "SC":
            warnings.append("SIMPLIFIED ABSORPTION MODEL: Subcutaneous administration assumes single first-order depot absorption (no lymphatic transport or local degradation).")
        elif route == "IP":
            warnings.append("SIMPLIFIED ABSORPTION MODEL: Intraperitoneal administration modeled as peritoneal absorption directly into portal/systemic circulation without human oral GI transit.")

    if request.dosing_frequency == "Repeated Dosing":
        warnings.append("LINEAR PK ASSUMPTION: Multi-dose superposition assumes linear kinetics without accumulation saturation.")

    # 10. Residuals & Goodness of Fit
    gof_metrics, residual_table = compute_goodness_of_fit(exp_obs, sim_res["time_series"]) if exp_obs else ({}, [])

    # 11. Parameter Snapshot & Output Metrics
    param_snapshot = {
        "dose_normalized_mg_kg": canon["dose_mg_kg"],
        "cl_l_h_kg": canon["cl_l_h_kg"],
        "cl_ml_min_kg": round(canon["cl_l_h_kg"] * 1000.0 / 60.0, 2),
        "v_l_kg": canon["v_l_kg"],
        "v_type": v_type,
        "f_fraction": f_val,
        "ka_h": ka_val,
        "k_elim_h": sim_res["k_elim"],
        "half_life_hours": sim_res["half_life_hours"],
    }

    output_metrics = {
        "c0_ng_ml": sim_res.get("c0_ng_ml"),
        "cmax_ng_ml": sim_res["cmax_ng_ml"],
        "tmax_hours": sim_res["tmax_hours"],
        "auc_last_ng_h_ml": sim_res["auc_last_ng_h_ml"],
        "auc_inf_analytical_ng_h_ml": sim_res.get("auc_inf_analytical_ng_h_ml"),
        "auc_inf_numerical_ng_h_ml": sim_res.get("auc_inf_numerical_ng_h_ml"),
        "auc_agreement_pct": sim_res.get("auc_agreement_pct"),
        "half_life_hours": sim_res["half_life_hours"],
        "is_flip_flop": flip_flop_detected,
        "goodness_of_fit": gof_metrics,
        "uncertainty_status": "UNCERTAINTY NOT QUANTIFIED",
    }

    absorption_components = {
        "fa": abs_decomp.get("fa_value"),
        "fa_status": abs_decomp.get("fa_status"),
        "fg": abs_decomp.get("fg_value"),
        "fg_status": abs_decomp.get("fg_status"),
        "fh": abs_decomp.get("fh_value"),
        "f_predicted": abs_decomp.get("f_predicted"),
        "f_experimental": route_set.get("f_experimental"),
    }

    provenance = {
        "engine_name": SIMULATION_ENGINE_NAME,
        "engine_version": SIMULATION_ENGINE_VERSION,
        "formula": "Analytical Linear PK Integration (First-Order Absorption & Elimination)",
        "units": {"dose": "mg/kg", "clearance": "L/h/kg", "volume": "L/kg", "concentration": "ng/mL", "time": "hours", "ka": "1/h"},
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    run_record = PKSimulationRun(
        project_id=project_id,
        compound_row_id=version.compound_row_id,
        version_id=version_id,
        species=species,
        route=route,
        administration_type=request.administration_type,
        dose=request.dose,
        dose_unit=request.dose_unit,
        infusion_duration_hours=request.infusion_duration_hours,
        dosing_frequency=request.dosing_frequency,
        dose_interval_hours=request.dose_interval_hours,
        num_doses=request.num_doses,
        model_type=request.model_type,
        f_value=f_val,
        f_source=f_source,
        ka_value=ka_val,
        ka_source=ka_source,
        flip_flop_flag=flip_flop_detected,
        parameter_snapshot=param_snapshot,
        parameter_sources=parameter_sources,
        absorption_components=absorption_components,
        steady_state_metrics=sim_res.get("steady_state", {}),
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
        route: str = Query("PO"),
        db: Session = Depends(get_db),
    ):
        route_clean = route.strip().upper()
        foundation = get_pk_foundation_profile(db, version_id, species)
        routes = foundation.get("route_parameter_sets", {})
        iv_set = routes.get("IV", {})
        target_set = routes.get(route_clean, {})
        dist = foundation.get("distribution", {})
        abs_decomp = foundation.get("absorption", {})

        cl_val = iv_set.get("cl_value") or target_set.get("clh_value")
        cl_source = iv_set.get("cl_source_type") if iv_set.get("cl_value") else ("PREDICTED_HEPATIC_IVIVE" if target_set.get("clh_value") else "UNAVAILABLE")
        cl_evidence = "EXPERIMENTAL" if "EXPERIMENTAL" in cl_source else ("DERIVED_ESTIMATE" if cl_val else "MODEL_UNAVAILABLE")

        v_val = dist.get("v_value")
        v_type = dist.get("v_type", "UNAVAILABLE")
        v_evidence = "EXPERIMENTAL" if v_type in ("Vss", "Vz") else ("DERIVED_ESTIMATE" if v_val else "MODEL_UNAVAILABLE")

        f_val = 100.0 if route_clean == "IV" else (target_set.get("f_experimental") or target_set.get("f_predicted"))
        f_source = "IV_ROUTE" if route_clean == "IV" else ("MATCHED_EXPERIMENTAL_F" if target_set.get("f_experimental") else ("MECHANISTIC_ASSEMBLY" if target_set.get("f_predicted") else "UNAVAILABLE"))
        f_evidence = "THEORETICAL_MAXIMUM" if route_clean == "IV" else ("EXPERIMENTAL" if target_set.get("f_experimental") else ("DERIVED_ESTIMATE" if target_set.get("f_predicted") else "MODEL_UNAVAILABLE"))

        # Look up ka for extravascular route
        extra_studies = db.scalars(
            select(PKStudy).where(
                PKStudy.version_id == version_id,
                PKStudy.species == species,
                PKStudy.route == route_clean,
            )
        ).all()
        extra_study = extra_studies[0] if extra_studies else None

        ka_val = None
        ka_source = "NOT_APPLICABLE" if route_clean == "IV" else "UNAVAILABLE"
        ka_evidence = "NOT_APPLICABLE" if route_clean == "IV" else "MODEL_UNAVAILABLE"
        if route_clean != "IV" and extra_study and extra_study.latest_nca and extra_study.latest_nca.tmax:
            tmax_obs = float(extra_study.latest_nca.tmax)
            if cl_val and v_val:
                ke_est = (cl_val * 60.0 / 1000.0) / v_val
                ka_sol = solve_ka_from_tmax(tmax_obs, ke_est)
                if ka_sol.get("status") == "CONVERGED":
                    ka_val = ka_sol["ka"]
                    ka_source = "EXPERIMENTAL_TMAX_DERIVED"
                    ka_evidence = "DERIVED_ESTIMATE"

        warnings = []
        if cl_source == "PREDICTED_HEPATIC_IVIVE":
            warnings.append("HEPATIC-CLEARANCE-ONLY APPROXIMATION: Simulation will use predicted hepatic CLh as systemic CL fallback.")
        if route_clean != "IV" and f_source == "MATCHED_EXPERIMENTAL_F" and (abs_decomp.get("fa_value") is None or abs_decomp.get("fg_value") is None):
            warnings.append("MECHANISTIC F INCOMPLETE: Matched experimental F will be used while mechanistic components are partially unavailable.")

        available_models = []
        if route_clean == "IV":
            available_models = [
                {"key": "ONE_COMPARTMENT_BOLUS", "name": "1-Compartment IV Bolus", "status": "AVAILABLE" if cl_val and v_val else "MODEL_UNAVAILABLE"},
                {"key": "ONE_COMPARTMENT_INFUSION", "name": "1-Compartment IV Infusion", "status": "AVAILABLE" if cl_val and v_val else "MODEL_UNAVAILABLE"},
                {"key": "TWO_COMPARTMENT", "name": "2-Compartment IV", "status": "REQUIRES_MICROCONSTANTS_OR_DENSE_DATA"},
            ]
        else:
            available_models = [
                {"key": "ONE_COMPARTMENT_EXTRAVASCULAR", "name": f"1-Compartment {route_clean} (1st-Order Absorption)", "status": "AVAILABLE" if cl_val and v_val and f_val and (ka_val or extra_study) else "REQUIRES_KA_OR_OVERRIDE"},
            ]

        confidence_levels = {"HIGH": 3, "MEDIUM": 2, "LOW": 1, "MODEL_UNAVAILABLE": 0}
        rev_levels = {3: "HIGH", 2: "MEDIUM", 1: "LOW", 0: "MODEL_UNAVAILABLE"}
        c_list = [target_set.get("confidence", "MEDIUM")]
        if not ka_val and route_clean != "IV":
            c_list.append("LOW")
        min_c = min(confidence_levels.get(c, 1) for c in c_list)

        return {
            "version_id": version_id,
            "species": species,
            "route": route_clean,
            "clearance": {"value": cl_val, "unit": "mL/min/kg", "source": cl_source, "evidence_type": cl_evidence},
            "volume": {"value": v_val, "unit": "L/kg", "type": v_type, "evidence_type": v_evidence},
            "bioavailability": {"value": f_val, "unit": "%", "source": f_source, "evidence_type": f_evidence},
            "absorption_rate": {"value": ka_val, "unit": "1/h", "source": ka_source, "evidence_type": ka_evidence},
            "mechanistic_components": {
                "fa": abs_decomp.get("fa_value"),
                "fg": abs_decomp.get("fg_value"),
                "fh": abs_decomp.get("fh_value"),
            },
            "available_models": available_models,
            "warnings": warnings,
            "confidence_ceiling": rev_levels[min_c],
        }

    @app.post("/api/compound-versions/{version_id}/pk-simulation/run")
    def run_simulation_endpoint(
        version_id: int,
        request: PKSimulationRequest,
        db: Session = Depends(get_db),
    ):
        run = run_pk_simulation(db, version_id, request)
        return run

    @app.post("/api/compound-versions/{version_id}/pk-simulation/fit-extravascular")
    def fit_extravascular_endpoint(
        version_id: int,
        request: PKFitRequest,
        db: Session = Depends(get_db),
    ):
        version = db.get(CompoundVersion, version_id)
        if not version:
            raise HTTPException(status_code=404, detail="CompoundVersion not found")

        species = request.species
        route = request.route.strip().upper()
        foundation = get_pk_foundation_profile(db, version_id, species)
        routes = foundation.get("route_parameter_sets", {})
        iv_set = routes.get("IV", {})
        dist = foundation.get("distribution", {})

        cl_val = request.user_cl_override or iv_set.get("cl_value")
        v_val = request.user_v_override or dist.get("v_value")

        if not cl_val or not v_val:
            raise HTTPException(status_code=400, detail="Systemic CL and V must be available or specified to fit extravascular absorption parameters.")

        cl_l_h_kg = float(cl_val) * 60.0 / 1000.0
        v_l_kg = float(v_val)

        studies = db.scalars(
            select(PKStudy).where(
                PKStudy.compound_row_id == version.compound_row_id,
                PKStudy.version_id == version_id,
                PKStudy.species == species,
                PKStudy.route == route,
            )
        ).all()
        if not studies or not studies[0].observations:
            raise HTTPException(status_code=400, detail=f"No experimental {route} PK observations available for species {species}.")

        fit_res = fit_one_compartment_extravascular(
            observations=list(studies[0].observations),
            dose_mg_kg=request.dose,
            cl_l_h_kg=cl_l_h_kg,
            v_l_kg=v_l_kg,
        )
        return fit_res

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
        route: str | None = Query(None),
        db: Session = Depends(get_db),
    ):
        stmt = select(PKSimulationRun).where(PKSimulationRun.version_id == version_id)
        if species:
            stmt = stmt.where(PKSimulationRun.species == species)
        if route:
            stmt = stmt.where(PKSimulationRun.route == route.strip().upper())
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
