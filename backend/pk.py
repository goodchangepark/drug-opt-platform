"""Experimental PK Data Management and Noncompartmental Analysis (NCA) Engine (Stage 5A-1)."""

from __future__ import annotations

import csv
import io
import math
from datetime import datetime, timezone
from typing import Any

from fastapi import Depends, HTTPException, Query, UploadFile, File
from pydantic import BaseModel
from sqlalchemy import JSON, Boolean, DateTime, Float, ForeignKey, Integer, String, Text, UniqueConstraint, inspect, select, text
from sqlalchemy.orm import Mapped, mapped_column, relationship, Session

from .database import Base, get_db
from .models import Compound, CompoundVersion, Project, utcnow

NCA_ENGINE_NAME = "Stage 5A-1 PK Noncompartmental Analysis Engine"
NCA_ENGINE_VERSION = "5A-1.0"


def ensure_pk_schema(engine):
    inspector = inspect(engine)
    if "projects" not in inspector.get_table_names():
        return
    Base.metadata.create_all(
        bind=engine,
        tables=[
            PKStudy.__table__,
            PKObservation.__table__,
            PKNCAResult.__table__,
        ],
    )


class PKStudy(Base):
    __tablename__ = "pk_studies"
    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    compound_row_id: Mapped[int] = mapped_column(ForeignKey("compounds.id", ondelete="CASCADE"), index=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    study_name: Mapped[str] = mapped_column(String(200), index=True)
    species: Mapped[str] = mapped_column(String(100), default="Rat", index=True)
    strain: Mapped[str] = mapped_column(String(100), default="")
    sex: Mapped[str] = mapped_column(String(40), default="Unknown")
    route: Mapped[str] = mapped_column(String(40), default="PO", index=True)
    dose: Mapped[float] = mapped_column(Float, default=10.0)
    dose_unit: Mapped[str] = mapped_column(String(40), default="mg/kg")
    dose_normalized_mg_kg: Mapped[float | None] = mapped_column(Float, nullable=True)
    formulation: Mapped[str] = mapped_column(String(200), default="")
    matrix: Mapped[str] = mapped_column(String(100), default="Plasma")
    dosing_frequency: Mapped[str] = mapped_column(String(60), default="Single Dose")
    fed_fasted: Mapped[str] = mapped_column(String(40), default="Fasted")
    lloq: Mapped[float | None] = mapped_column(Float, nullable=True)
    lloq_unit: Mapped[str] = mapped_column(String(40), default="ng/mL")
    study_date: Mapped[str] = mapped_column(String(30), default="")
    source: Mapped[str] = mapped_column(String(200), default="")
    notes: Mapped[str] = mapped_column(Text, default="")
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    updated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow, onupdate=utcnow)

    project = relationship("Project")
    compound = relationship("Compound")
    version = relationship("CompoundVersion")
    observations = relationship("PKObservation", back_populates="study", cascade="all, delete-orphan")
    nca_results = relationship("PKNCAResult", back_populates="study", cascade="all, delete-orphan")

    @property
    def latest_nca(self):
        latest = [n for n in (self.nca_results or []) if getattr(n, "is_latest", True)]
        return latest[0] if latest else (self.nca_results[-1] if self.nca_results else None)


class PKObservation(Base):
    __tablename__ = "pk_observations"
    id: Mapped[int] = mapped_column(primary_key=True)
    pk_study_id: Mapped[int] = mapped_column(ForeignKey("pk_studies.id", ondelete="CASCADE"), index=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    subject_group_id: Mapped[str] = mapped_column(String(100), default="Group Mean")
    time_raw: Mapped[float] = mapped_column(Float)
    time_unit: Mapped[str] = mapped_column(String(30), default="h")
    time_hours: Mapped[float] = mapped_column(Float)
    concentration_raw: Mapped[float | None] = mapped_column(Float, nullable=True)
    concentration_unit: Mapped[str] = mapped_column(String(40), default="ng/mL")
    concentration_normalized_ng_ml: Mapped[float | None] = mapped_column(Float, nullable=True)
    blq_flag: Mapped[bool] = mapped_column(Boolean, default=False)
    replicate: Mapped[str] = mapped_column(String(40), default="R1")
    notes: Mapped[str] = mapped_column(Text, default="")
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)

    study = relationship("PKStudy", back_populates="observations")
    version = relationship("CompoundVersion")


class PKNCAResult(Base):
    __tablename__ = "pk_nca_results"
    id: Mapped[int] = mapped_column(primary_key=True)
    pk_study_id: Mapped[int] = mapped_column(ForeignKey("pk_studies.id", ondelete="CASCADE"), index=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    analysis_version: Mapped[int] = mapped_column(Integer, default=1)
    is_latest: Mapped[bool] = mapped_column(Boolean, default=True)
    selection_mode: Mapped[str] = mapped_column(String(40), default="AUTO")
    subject_group_id: Mapped[str] = mapped_column(String(100), default="Group Mean")
    cmax: Mapped[float | None] = mapped_column(Float, nullable=True)
    cmax_unit: Mapped[str] = mapped_column(String(40), default="ng/mL")
    tmax: Mapped[float | None] = mapped_column(Float, nullable=True)
    tmax_unit: Mapped[str] = mapped_column(String(40), default="h")
    auclast: Mapped[float | None] = mapped_column(Float, nullable=True)
    auclast_unit: Mapped[str] = mapped_column(String(40), default="ng*h/mL")
    aucinf: Mapped[float | None] = mapped_column(Float, nullable=True)
    aucinf_unit: Mapped[str] = mapped_column(String(40), default="ng*h/mL")
    lambda_z: Mapped[float | None] = mapped_column(Float, nullable=True)
    terminal_half_life: Mapped[float | None] = mapped_column(Float, nullable=True)
    mrt: Mapped[float | None] = mapped_column(Float, nullable=True)
    cl: Mapped[float | None] = mapped_column(Float, nullable=True)
    cl_unit: Mapped[str] = mapped_column(String(40), default="mL/min/kg")
    cl_f: Mapped[float | None] = mapped_column(Float, nullable=True)
    cl_f_unit: Mapped[str] = mapped_column(String(40), default="mL/min/kg")
    vz: Mapped[float | None] = mapped_column(Float, nullable=True)
    vz_unit: Mapped[str] = mapped_column(String(40), default="L/kg")
    vz_f: Mapped[float | None] = mapped_column(Float, nullable=True)
    vz_f_unit: Mapped[str] = mapped_column(String(40), default="L/kg")
    aumclast: Mapped[float | None] = mapped_column(Float, nullable=True)
    aumcinf: Mapped[float | None] = mapped_column(Float, nullable=True)
    auc_extrapolated_pct: Mapped[float | None] = mapped_column(Float, nullable=True)
    terminal_point_count: Mapped[int] = mapped_column(Integer, default=0)
    terminal_points_json: Mapped[list] = mapped_column(JSON, default=list)
    r_squared: Mapped[float | None] = mapped_column(Float, nullable=True)
    adjusted_r2: Mapped[float | None] = mapped_column(Float, nullable=True)
    warnings_json: Mapped[list] = mapped_column(JSON, default=list)
    blq_policy_json: Mapped[dict] = mapped_column(JSON, default=dict)
    nca_engine: Mapped[str] = mapped_column(String(100), default=NCA_ENGINE_NAME)
    nca_engine_version: Mapped[str] = mapped_column(String(40), default=NCA_ENGINE_VERSION)
    calculation_method: Mapped[str] = mapped_column(String(100), default="Linear-up / Log-down trapezoidal")
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)
    calculation_timestamp: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)

    study = relationship("PKStudy", back_populates="nca_results")
    version = relationship("CompoundVersion")


# --- UNIT CONVERSION ENGINE ---

TIME_UNITS = {"sec": 1 / 3600.0, "min": 1 / 60.0, "h": 1.0, "hr": 1.0, "hours": 1.0, "day": 24.0, "days": 24.0}
CONC_MASS_FACTOR = {"pg/mL": 0.001, "ng/mL": 1.0, "µg/mL": 1000.0, "ug/mL": 1000.0, "mg/L": 1000.0}
DOSE_MASS_FACTOR = {"mg/kg": 1.0, "µg/kg": 0.001, "ug/kg": 0.001, "mg": 1.0, "µg": 0.001, "ug": 0.001}


def normalize_time_to_hours(value: float, unit: str) -> float:
    unit_clean = unit.strip().lower()
    if unit_clean not in TIME_UNITS:
        raise HTTPException(status_code=400, detail=f"Unsupported time unit: {unit}. Supported: sec, min, h, day")
    return value * TIME_UNITS[unit_clean]


def normalize_conc_to_ng_ml(value: float | None, unit: str, mw: float | None = None) -> float | None:
    if value is None:
        return None
    unit_clean = unit.strip()
    if unit_clean in CONC_MASS_FACTOR:
        return value * CONC_MASS_FACTOR[unit_clean]
    elif unit_clean in {"nM", "nmol/L"}:
        if not mw or mw <= 0:
            raise HTTPException(status_code=400, detail="Compound MW is required to convert molar concentration (nM) to mass units (ng/mL)")
        return value * (mw / 1000.0)
    elif unit_clean in {"µM", "uM", "umol/L"}:
        if not mw or mw <= 0:
            raise HTTPException(status_code=400, detail="Compound MW is required to convert molar concentration (µM) to mass units (ng/mL)")
        return value * mw
    else:
        raise HTTPException(status_code=400, detail=f"Unsupported concentration unit: {unit}. Supported: pg/mL, ng/mL, µg/mL, mg/L, nM, µM")


def normalize_dose_to_mg_kg(value: float, unit: str) -> float:
    unit_clean = unit.strip().lower()
    if unit_clean not in DOSE_MASS_FACTOR:
        raise HTTPException(status_code=400, detail=f"Unsupported dose unit: {unit}. Supported: mg/kg, µg/kg, mg, µg")
    return value * DOSE_MASS_FACTOR[unit_clean]


# --- NONCOMPARTMENTAL ANALYSIS (NCA) ALGORITHM ---

def run_nca_calculation(
    observations: list[PKObservation],
    route: str,
    dose_mg_kg: float,
    dose_unit: str,
    blq_policy: dict | None = None,
    manual_terminal_indices: list[int] | None = None,
    lloq_val: float | None = None,
) -> dict[str, Any]:
    """
    Executes Noncompartmental Analysis (NCA) using Linear-up / Log-down trapezoidal integration.
    Strictly distinguishes IV (CL, Vz) from Extravascular PO/SC/IP (CL/F, Vz/F).
    """
    blq_policy = blq_policy or {"pre_first": "zero", "intermittent": "missing", "terminal": "missing"}
    warnings = []
    
    if not observations:
        return {
            "status": "FAILED",
            "warnings": ["No observations provided for NCA analysis."],
            "cmax": None, "tmax": None, "auclast": None, "aucinf": None,
            "lambda_z": None, "terminal_half_life": None, "mrt": None,
            "cl": None, "cl_f": None, "vz": None, "vz_f": None,
        }

    # Sort observations by time in hours
    sorted_obs = sorted(observations, key=lambda x: x.time_hours)
    
    # Process BLQ points according to policy
    processed_pts = []
    first_quant_found = False
    
    for obs in sorted_obs:
        is_blq = obs.blq_flag or (obs.concentration_normalized_ng_ml is None) or (lloq_val is not None and obs.concentration_normalized_ng_ml is not None and obs.concentration_normalized_ng_ml < lloq_val)
        conc = obs.concentration_normalized_ng_ml
        
        if is_blq:
            if not first_quant_found:
                if blq_policy.get("pre_first") == "zero":
                    processed_pts.append({"time": obs.time_hours, "conc": 0.0, "blq": True, "obs_id": obs.id})
                # if "missing", skip
            else:
                # Intermittent or terminal BLQ
                if blq_policy.get("intermittent") == "zero":
                    processed_pts.append({"time": obs.time_hours, "conc": 0.0, "blq": True, "obs_id": obs.id})
                # if "missing", skip
        else:
            first_quant_found = True
            processed_pts.append({"time": obs.time_hours, "conc": float(conc), "blq": False, "obs_id": obs.id})

    quant_pts = [p for p in processed_pts if p["conc"] > 0]
    
    if not quant_pts:
        return {
            "status": "FAILED",
            "warnings": ["All concentration points are BLQ or zero; cannot compute NCA parameters."],
            "cmax": 0.0, "tmax": 0.0, "auclast": 0.0, "aucinf": None,
            "lambda_z": None, "terminal_half_life": None, "mrt": None,
            "cl": None, "cl_f": None, "vz": None, "vz_f": None,
        }

    # 1. Cmax and Tmax
    cmax_pt = max(quant_pts, key=lambda p: p["conc"])
    cmax = cmax_pt["conc"]
    tmax = cmax_pt["time"]

    # 2. AUC & AUMC Trapezoidal Integration (Linear-up / Log-down)
    auclast = 0.0
    aumclast = 0.0
    
    for i in range(len(processed_pts) - 1):
        t1, c1 = processed_pts[i]["time"], processed_pts[i]["conc"]
        t2, c2 = processed_pts[i + 1]["time"], processed_pts[i + 1]["conc"]
        dt = t2 - t1
        if dt <= 0:
            continue
        
        if c2 >= c1:
            # Linear-up or flat
            d_auc = 0.5 * (c1 + c2) * dt
            d_aumc = 0.5 * (t1 * c1 + t2 * c2) * dt
        else:
            # Log-down (c2 < c1)
            if c2 > 0 and c1 > 0:
                ln_ratio = math.log(c1 / c2)
                d_auc = (c1 - c2) / ln_ratio * dt
                d_aumc = (t1 * c1 - t2 * c2) / ln_ratio * dt + (c1 - c2) / (ln_ratio ** 2) * (dt ** 2)
            else:
                # One of them is 0
                d_auc = 0.5 * (c1 + c2) * dt
                d_aumc = 0.5 * (t1 * c1 + t2 * c2) * dt
                
        auclast += d_auc
        aumclast += d_aumc

    t_last = quant_pts[-1]["time"]
    c_last = quant_pts[-1]["conc"]

    # 3. Terminal Elimination Phase Regression (Lambda_z)
    lambda_z = None
    r2 = None
    adj_r2 = None
    term_pts = []

    if manual_terminal_indices and len(manual_terminal_indices) >= 2:
        # Manual selection by observation IDs or point indices
        selected_pts = [p for p in quant_pts if p["obs_id"] in manual_terminal_indices or quant_pts.index(p) in manual_terminal_indices]
        if len(selected_pts) >= 2:
            term_pts = selected_pts
            selection_mode = "MANUAL_OVERRIDE"
        else:
            warnings.append("Manual terminal points invalid or < 2 points; falling back to automatic terminal selection.")
            selection_mode = "AUTO"
    else:
        selection_mode = "AUTO"

    if selection_mode == "AUTO":
        # Automated terminal point selection: candidate points at or after Tmax
        candidates = [p for p in quant_pts if p["time"] >= tmax]
        if len(candidates) < 3:
            candidates = quant_pts  # fallback if Tmax is late
            
        best_adj_r2 = -999.0
        best_fit = None
        
        # Test all terminal windows of length >= 3
        for start_idx in range(len(candidates) - 2):
            sub_pts = candidates[start_idx:]
            n_sub = len(sub_pts)
            if n_sub < 3:
                continue
            
            times = [p["time"] for p in sub_pts]
            log_concs = [math.log(p["conc"]) for p in sub_pts]
            
            t_mean = sum(times) / n_sub
            lc_mean = sum(log_concs) / n_sub
            
            ss_tt = sum((t - t_mean) ** 2 for t in times)
            ss_lc_lc = sum((lc - lc_mean) ** 2 for lc in log_concs)
            ss_t_lc = sum((t - t_mean) * (lc - lc_mean) for t, lc in zip(times, log_concs))
            
            if ss_tt > 0 and ss_lc_lc > 0:
                slope = ss_t_lc / ss_tt
                if slope < 0:  # Must be negative slope for elimination
                    lz = -slope
                    r2_val = (ss_t_lc ** 2) / (ss_tt * ss_lc_lc)
                    adj_r2_val = 1.0 - (1.0 - r2_val) * (n_sub - 1) / (n_sub - 2)
                    
                    if adj_r2_val > best_adj_r2:
                        best_adj_r2 = adj_r2_val
                        best_fit = (lz, r2_val, adj_r2_val, sub_pts)
                        
        if best_fit:
            lambda_z, r2, adj_r2, term_pts = best_fit

    elif term_pts:
        # Perform manual regression fit
        n_sub = len(term_pts)
        times = [p["time"] for p in term_pts]
        log_concs = [math.log(p["conc"]) for p in term_pts]
        t_mean = sum(times) / n_sub
        lc_mean = sum(log_concs) / n_sub
        ss_tt = sum((t - t_mean) ** 2 for t in times)
        ss_lc_lc = sum((lc - lc_mean) ** 2 for lc in log_concs)
        ss_t_lc = sum((t - t_mean) * (lc - lc_mean) for t, lc in zip(times, log_concs))
        
        if ss_tt > 0 and ss_lc_lc > 0 and ss_t_lc < 0:
            lambda_z = - (ss_t_lc / ss_tt)
            r2 = (ss_t_lc ** 2) / (ss_tt * ss_lc_lc)
            adj_r2 = 1.0 - (1.0 - r2) * (n_sub - 1) / max(1, (n_sub - 2))

    # Extrapolated parameters
    auc_inf = None
    aumc_inf = None
    half_life = None
    auc_extrap_pct = None
    mrt = None
    cl_iv = None
    cl_f_ev = None
    vz_iv = None
    vz_f_ev = None

    if lambda_z and lambda_z > 0:
        half_life = math.log(2.0) / lambda_z
        auc_inf = auclast + (c_last / lambda_z)
        aumc_inf = aumclast + (t_last * c_last / lambda_z) + (c_last / (lambda_z ** 2))
        
        if auc_inf > 0:
            auc_extrap_pct = ((auc_inf - auclast) / auc_inf) * 100.0
            mrt = aumc_inf / auc_inf
            
            # Dose in mg/kg -> CL in L/h/kg: Dose (mg/kg) / AUCinf (ng*h/mL)
            # 1 mg/kg = 1,000,000 ng/kg. AUCinf in ng*h/mL = ng*h / (0.001 L) = 1,000 ng*h/L.
            # CL = (1,000,000 ng/kg) / (AUCinf * 1,000 ng*h/L) = 1,000 / AUCinf  [L/h/kg].
            # Convert to mL/min/kg: L/h/kg * (1000 mL/L) / (60 min/h) = (1,000 / AUCinf) * (1000/60) = 16666.6667 / AUCinf [mL/min/kg].
            if dose_mg_kg and dose_mg_kg > 0:
                raw_cl_l_h_kg = (dose_mg_kg * 1000.0) / auc_inf
                raw_cl_ml_min_kg = raw_cl_l_h_kg * (1000.0 / 60.0)
                
                route_clean = route.strip().upper()
                if route_clean == "IV":
                    cl_iv = raw_cl_ml_min_kg
                    vz_iv = (raw_cl_l_h_kg / lambda_z)  # L/kg
                else:
                    cl_f_ev = raw_cl_ml_min_kg
                    vz_f_ev = (raw_cl_l_h_kg / lambda_z)  # L/kg

        # Check warnings
        if auc_extrap_pct is not None and auc_extrap_pct > 20.0:
            warnings.append(f"AUC extrapolation exceeds 20% ({auc_extrap_pct:.1f}%); terminal clearance/half-life may be uncertain.")
        if r2 is not None and r2 < 0.85:
            warnings.append(f"Terminal elimination R² is low ({r2:.3f}); half-life estimate may be unreliable.")
        if len(term_pts) < 3:
            warnings.append(f"Terminal elimination phase based on only {len(term_pts)} points.")
    else:
        warnings.append("Could not determine valid terminal elimination phase (lambda_z). AUCinf, t1/2, CL, and Vz cannot be calculated.")

    return {
        "status": "COMPLETE",
        "selection_mode": selection_mode,
        "cmax": cmax,
        "cmax_unit": "ng/mL",
        "tmax": tmax,
        "tmax_unit": "h",
        "auclast": auclast,
        "auclast_unit": "ng*h/mL",
        "aucinf": auc_inf,
        "aucinf_unit": "ng*h/mL",
        "lambda_z": lambda_z,
        "terminal_half_life": half_life,
        "mrt": mrt,
        "cl": cl_iv,
        "cl_unit": "mL/min/kg",
        "cl_f": cl_f_ev,
        "cl_f_unit": "mL/min/kg",
        "vz": vz_iv,
        "vz_unit": "L/kg",
        "vz_f": vz_f_ev,
        "vz_f_unit": "L/kg",
        "aumclast": aumclast,
        "aumcinf": aumc_inf,
        "auc_extrapolated_pct": auc_extrap_pct,
        "terminal_point_count": len(term_pts),
        "terminal_points": [p["obs_id"] for p in term_pts],
        "r_squared": r2,
        "adjusted_r2": adj_r2,
        "warnings": warnings,
        "blq_policy": blq_policy,
        "nca_engine": NCA_ENGINE_NAME,
        "nca_engine_version": NCA_ENGINE_VERSION,
        "calculation_method": "Linear-up / Log-down trapezoidal",
    }


# --- CSV PARSER WITH FLEXIBLE MAPPING ---

def parse_pk_csv(content: str, column_mapping: dict[str, str] | None = None) -> tuple[list[dict], list[str], list[str]]:
    """
    Parses CSV text into observation dictionaries.
    Supports flexible column mapping: e.g. {"TIME_HR": "time", "CONC_NG_ML": "concentration"}.
    Returns (valid_rows, errors, detected_columns).
    """
    reader = csv.DictReader(io.StringIO(content))
    fieldnames = list(reader.fieldnames or [])
    
    # Auto-detect mappings if column_mapping is not provided
    mapping = dict(column_mapping or {})
    if not mapping:
        for col in fieldnames:
            col_lower = col.strip().lower()
            if col_lower in {"time", "time_h", "time_hr", "time_hours", "t", "time(h)"}:
                mapping[col] = "time"
            elif col_lower in {"concentration", "conc", "conc_ng_ml", "cp", "val", "value", "conc(ng/ml)"}:
                mapping[col] = "concentration"
            elif col_lower in {"subject", "animal", "group", "subject_id", "animal_id"}:
                mapping[col] = "subject"
            elif col_lower in {"replicate", "rep"}:
                mapping[col] = "replicate"
            elif col_lower in {"blq", "is_blq"}:
                mapping[col] = "blq"

    # Invert mapping: target -> list of CSV cols
    time_col = next((k for k, v in mapping.items() if v == "time"), None)
    conc_col = next((k for k, v in mapping.items() if v == "concentration"), None)
    subj_col = next((k for k, v in mapping.items() if v == "subject"), None)
    rep_col = next((k for k, v in mapping.items() if v == "replicate"), None)
    blq_col = next((k for k, v in mapping.items() if v == "blq"), None)

    valid_rows = []
    errors = []

    for idx, row in enumerate(reader, start=1):
        if not time_col or time_col not in row:
            errors.append(f"Row {idx}: Missing mapped time column")
            continue
        if not conc_col or conc_col not in row:
            errors.append(f"Row {idx}: Missing mapped concentration column")
            continue

        raw_time = row[time_col].strip() if row[time_col] else ""
        raw_conc = row[conc_col].strip() if row[conc_col] else ""

        try:
            val_time = float(raw_time)
        except ValueError:
            errors.append(f"Row {idx}: Invalid numeric time '{raw_time}'")
            continue

        is_blq = False
        val_conc = None
        if blq_col and row.get(blq_col):
            b_val = str(row[blq_col]).strip().lower()
            if b_val in {"1", "true", "yes", "blq"}:
                is_blq = True

        if raw_conc.upper() in {"BLQ", "<LLOQ", "ND", "BQL"} or raw_conc == "":
            is_blq = True
        else:
            try:
                val_conc = float(raw_conc)
                if val_conc <= 0:
                    is_blq = True
            except ValueError:
                is_blq = True

        valid_rows.append({
            "row_number": idx,
            "subject": row.get(subj_col, "Group Mean").strip() if subj_col else "Group Mean",
            "time_raw": val_time,
            "time_unit": "h",
            "concentration_raw": val_conc,
            "concentration_unit": "ng/mL",
            "blq_flag": is_blq,
            "replicate": row.get(rep_col, "R1").strip() if rep_col else "R1",
        })

    return valid_rows, errors, fieldnames


# --- PYDANTIC SCHEMAS FOR API ---

class PKStudyCreate(BaseModel):
    study_name: str
    species: str = "Rat"
    strain: str = ""
    sex: str = "Unknown"
    route: str = "PO"
    dose: float = 10.0
    dose_unit: str = "mg/kg"
    formulation: str = ""
    matrix: str = "Plasma"
    dosing_frequency: str = "Single Dose"
    fed_fasted: str = "Fasted"
    lloq: float | None = None
    lloq_unit: str = "ng/mL"
    study_date: str = ""
    source: str = ""
    notes: str = ""


class PKObservationCreate(BaseModel):
    subject_group_id: str = "Group Mean"
    time_raw: float
    time_unit: str = "h"
    concentration_raw: float | None = None
    concentration_unit: str = "ng/mL"
    blq_flag: bool = False
    replicate: str = "R1"
    notes: str = ""


class PKNCARunOptions(BaseModel):
    blq_policy: dict = {"pre_first": "zero", "intermittent": "missing", "terminal": "missing"}
    manual_terminal_indices: list[int] | None = None
    lloq_val: float | None = None
    subject_group_id: str | None = None


# --- SERIALIZERS ---

def serialize_study(study: PKStudy, latest_nca: PKNCAResult | None = None) -> dict[str, Any]:
    return {
        "id": study.id,
        "project_id": study.project_id,
        "compound_row_id": study.compound_row_id,
        "version_id": study.version_id,
        "study_name": study.study_name,
        "species": study.species,
        "strain": study.strain,
        "sex": study.sex,
        "route": study.route,
        "dose": study.dose,
        "dose_unit": study.dose_unit,
        "dose_normalized_mg_kg": study.dose_normalized_mg_kg,
        "formulation": study.formulation,
        "matrix": study.matrix,
        "dosing_frequency": study.dosing_frequency,
        "fed_fasted": study.fed_fasted,
        "lloq": study.lloq,
        "lloq_unit": study.lloq_unit,
        "study_date": study.study_date,
        "source": study.source,
        "notes": study.notes,
        "observation_count": len(study.observations or []),
        "latest_nca": serialize_nca(latest_nca) if latest_nca else None,
        "created_at": study.created_at.isoformat() if study.created_at else "",
    }


def serialize_observation(obs: PKObservation) -> dict[str, Any]:
    return {
        "id": obs.id,
        "pk_study_id": obs.pk_study_id,
        "version_id": obs.version_id,
        "subject_group_id": obs.subject_group_id,
        "time_raw": obs.time_raw,
        "time_unit": obs.time_unit,
        "time_hours": obs.time_hours,
        "concentration_raw": obs.concentration_raw,
        "concentration_unit": obs.concentration_unit,
        "concentration_normalized_ng_ml": obs.concentration_normalized_ng_ml,
        "blq_flag": obs.blq_flag,
        "replicate": obs.replicate,
        "notes": obs.notes,
        "created_at": obs.created_at.isoformat() if obs.created_at else "",
    }


def serialize_nca(nca: PKNCAResult) -> dict[str, Any]:
    return {
        "id": nca.id,
        "pk_study_id": nca.pk_study_id,
        "version_id": nca.version_id,
        "analysis_version": nca.analysis_version,
        "is_latest": nca.is_latest,
        "selection_mode": nca.selection_mode,
        "subject_group_id": nca.subject_group_id,
        "cmax": nca.cmax,
        "cmax_unit": nca.cmax_unit,
        "tmax": nca.tmax,
        "tmax_unit": nca.tmax_unit,
        "auclast": nca.auclast,
        "auclast_unit": nca.auclast_unit,
        "aucinf": nca.aucinf,
        "aucinf_unit": nca.aucinf_unit,
        "lambda_z": nca.lambda_z,
        "terminal_half_life": nca.terminal_half_life,
        "mrt": nca.mrt,
        "cl": nca.cl,
        "cl_unit": nca.cl_unit,
        "cl_f": nca.cl_f,
        "cl_f_unit": nca.cl_f_unit,
        "vz": nca.vz,
        "vz_unit": nca.vz_unit,
        "vz_f": nca.vz_f,
        "vz_f_unit": nca.vz_f_unit,
        "aumclast": nca.aumclast,
        "aumcinf": nca.aumcinf,
        "auc_extrapolated_pct": nca.auc_extrapolated_pct,
        "terminal_point_count": nca.terminal_point_count,
        "terminal_points": nca.terminal_points_json or [],
        "r_squared": nca.r_squared,
        "adjusted_r2": nca.adjusted_r2,
        "warnings": nca.warnings_json or [],
        "blq_policy": nca.blq_policy_json or {},
        "nca_engine": nca.nca_engine,
        "nca_engine_version": nca.nca_engine_version,
        "calculation_method": nca.calculation_method,
        "provenance": nca.provenance_json or {},
        "calculation_timestamp": nca.calculation_timestamp.isoformat() if nca.calculation_timestamp else "",
    }


def calculate_bioavailability_for_version(version_id: int, db: Session) -> dict[str, Any]:
    """
    Computes absolute bioavailability (F) for extravascular studies (PO, SC, IP)
    when a matched IV study exists for the same CompoundVersion and species.
    """
    version = db.get(CompoundVersion, version_id)
    if not version:
        return {"bioavailability": [], "message": "CompoundVersion not found"}

    studies = list(db.scalars(
        select(PKStudy).where(PKStudy.version_id == version_id)
    ))

    if not studies:
        return {"bioavailability": [], "message": "No PK studies available for this CompoundVersion"}

    # Group by species
    by_species: dict[str, list[PKStudy]] = {}
    for s in studies:
        by_species.setdefault(s.species.strip(), []).append(s)

    results = []

    for species, s_list in by_species.items():
        iv_studies = [s for s in s_list if s.route.strip().upper() == "IV"]
        extra_studies = [s for s in s_list if s.route.strip().upper() in {"PO", "SC", "IP"}]

        if not iv_studies:
            for s in extra_studies:
                results.append({
                    "species": species,
                    "route": s.route,
                    "study_id": s.id,
                    "study_name": s.study_name,
                    "bioavailability_pct": None,
                    "status": "NO_MATCHED_IV",
                    "message": "Absolute bioavailability cannot be calculated without a matched IV study.",
                })
            continue

        # Use the latest valid IV study with NCA result
        best_iv = None
        best_iv_nca = None
        for iv_s in iv_studies:
            latest = db.scalar(
                select(PKNCAResult).where(
                    PKNCAResult.pk_study_id == iv_s.id,
                    PKNCAResult.is_latest == True
                ).order_by(PKNCAResult.analysis_version.desc())
            )
            if latest and (latest.aucinf or latest.auclast):
                best_iv = iv_s
                best_iv_nca = latest
                break

        if not best_iv or not best_iv_nca:
            for s in extra_studies:
                results.append({
                    "species": species,
                    "route": s.route,
                    "study_id": s.id,
                    "study_name": s.study_name,
                    "bioavailability_pct": None,
                    "status": "NO_VALID_IV_NCA",
                    "message": "Absolute bioavailability cannot be calculated without a matched IV study with valid NCA AUC.",
                })
            continue

        iv_dose = best_iv.dose_normalized_mg_kg or normalize_dose_to_mg_kg(best_iv.dose, best_iv.dose_unit)
        iv_auc = best_iv_nca.aucinf if best_iv_nca.aucinf else best_iv_nca.auclast
        iv_auc_type = "AUCinf" if best_iv_nca.aucinf else "AUClast"

        if not iv_dose or iv_dose <= 0 or not iv_auc or iv_auc <= 0:
            for s in extra_studies:
                results.append({
                    "species": species,
                    "route": s.route,
                    "study_id": s.id,
                    "study_name": s.study_name,
                    "bioavailability_pct": None,
                    "status": "INVALID_IV_DOSE_OR_AUC",
                    "message": "Absolute bioavailability cannot be calculated without valid positive IV dose and AUC.",
                })
            continue

        for extra_s in extra_studies:
            extra_nca = db.scalar(
                select(PKNCAResult).where(
                    PKNCAResult.pk_study_id == extra_s.id,
                    PKNCAResult.is_latest == True
                ).order_by(PKNCAResult.analysis_version.desc())
            )
            if not extra_nca or not (extra_nca.aucinf or extra_nca.auclast):
                results.append({
                    "species": species,
                    "route": extra_s.route,
                    "study_id": extra_s.id,
                    "study_name": extra_s.study_name,
                    "bioavailability_pct": None,
                    "status": "NO_EXTRA_NCA",
                    "message": f"NCA analysis required for {extra_s.route} study to calculate bioavailability.",
                })
                continue

            extra_dose = extra_s.dose_normalized_mg_kg or normalize_dose_to_mg_kg(extra_s.dose, extra_s.dose_unit)
            extra_auc = extra_nca.aucinf if extra_nca.aucinf else extra_nca.auclast
            extra_auc_type = "AUCinf" if extra_nca.aucinf else "AUClast"

            if not extra_dose or extra_dose <= 0 or not extra_auc or extra_auc <= 0:
                results.append({
                    "species": species,
                    "route": extra_s.route,
                    "study_id": extra_s.id,
                    "study_name": extra_s.study_name,
                    "bioavailability_pct": None,
                    "status": "INVALID_EXTRA_DOSE_OR_AUC",
                    "message": "Invalid dose or AUC for bioavailability calculation.",
                })
                continue

            # F = (AUC_extra / Dose_extra) / (AUC_IV / Dose_IV) * 100
            f_val = ((extra_auc / extra_dose) / (iv_auc / iv_dose)) * 100.0

            results.append({
                "species": species,
                "route": extra_s.route,
                "label": f"F_{extra_s.route}",
                "study_id": extra_s.id,
                "study_name": extra_s.study_name,
                "matched_iv_study_id": best_iv.id,
                "matched_iv_study_name": best_iv.study_name,
                "bioavailability_pct": round(f_val, 2),
                "status": "MATCHED",
                "message": f"Calculated against matched {species} IV study #{best_iv.id} using {extra_auc_type}/{iv_auc_type}.",
                "provenance": {
                    "extra_dose_mg_kg": extra_dose,
                    "extra_auc": extra_auc,
                    "extra_auc_type": extra_auc_type,
                    "iv_dose_mg_kg": iv_dose,
                    "iv_auc": iv_auc,
                    "iv_auc_type": iv_auc_type,
                }
            })

    return {"bioavailability": results, "compound_version_id": version_id}


# --- FASTAPI ROUTER ATTACHMENT ---

def register_pk_routes(app):
    @app.post("/api/compounds/{row_id}/pk-studies", status_code=201)
    def create_pk_study_endpoint(row_id: int, payload: PKStudyCreate, db: Session = Depends(get_db)):
        compound = db.get(Compound, row_id)
        if not compound:
            raise HTTPException(status_code=404, detail="Compound not found")
        version = next((v for v in compound.versions if v.version_number == compound.current_version), None)
        if not version:
            raise HTTPException(status_code=400, detail="Draw or enter a valid structure version before adding a PK study.")
        
        dose_norm = normalize_dose_to_mg_kg(payload.dose, payload.dose_unit)
        study = PKStudy(
            project_id=compound.project_id,
            compound_row_id=compound.id,
            version_id=version.id,
            study_name=payload.study_name.strip(),
            species=payload.species.strip(),
            strain=payload.strain.strip(),
            sex=payload.sex.strip(),
            route=payload.route.strip().upper(),
            dose=payload.dose,
            dose_unit=payload.dose_unit.strip(),
            dose_normalized_mg_kg=dose_norm,
            formulation=payload.formulation.strip(),
            matrix=payload.matrix.strip(),
            dosing_frequency=payload.dosing_frequency.strip(),
            fed_fasted=payload.fed_fasted.strip(),
            lloq=payload.lloq,
            lloq_unit=payload.lloq_unit.strip(),
            study_date=payload.study_date.strip(),
            source=payload.source.strip(),
            notes=payload.notes.strip(),
        )
        db.add(study)
        db.commit()
        db.refresh(study)
        try:
            from .ivive import refresh_pk_and_ivive_for_version
            refresh_pk_and_ivive_for_version(db, study.version_id, force=True)
        except Exception:
            pass
        return serialize_study(study)

    @app.get("/api/compounds/{row_id}/pk-studies")
    def list_pk_studies_endpoint(row_id: int, version_id: int | None = Query(None), db: Session = Depends(get_db)):
        compound = db.get(Compound, row_id)
        if not compound:
            raise HTTPException(status_code=404, detail="Compound not found")
        target_version_id = version_id
        if not target_version_id:
            version = next((v for v in compound.versions if v.version_number == compound.current_version), None)
            if version:
                target_version_id = version.id

        if not target_version_id:
            return {"studies": [], "bioavailability": []}

        studies = list(db.scalars(
            select(PKStudy).where(PKStudy.compound_row_id == row_id, PKStudy.version_id == target_version_id).order_by(PKStudy.created_at.desc())
        ))

        study_rows = []
        for study in studies:
            latest_nca = db.scalar(
                select(PKNCAResult).where(PKNCAResult.pk_study_id == study.id, PKNCAResult.is_latest == True).order_by(PKNCAResult.analysis_version.desc())
            )
            study_rows.append(serialize_study(study, latest_nca))

        bioavailability = calculate_bioavailability_for_version(target_version_id, db)
        return {
            "compound_id": compound.compound_id,
            "version_id": target_version_id,
            "studies": study_rows,
            "bioavailability": bioavailability.get("bioavailability", []),
        }

    @app.get("/api/pk-studies/{study_id}")
    def get_pk_study_endpoint(study_id: int, db: Session = Depends(get_db)):
        study = db.get(PKStudy, study_id)
        if not study:
            raise HTTPException(status_code=404, detail="PK study not found")

        observations = list(db.scalars(
            select(PKObservation).where(PKObservation.pk_study_id == study_id).order_by(PKObservation.time_hours, PKObservation.subject_group_id)
        ))

        nca_history = list(db.scalars(
            select(PKNCAResult).where(PKNCAResult.pk_study_id == study_id).order_by(PKNCAResult.analysis_version.desc())
        ))
        latest_nca = nca_history[0] if nca_history else None

        bioavailability = calculate_bioavailability_for_version(study.version_id, db)

        return {
            "study": serialize_study(study, latest_nca),
            "observations": [serialize_observation(obs) for obs in observations],
            "nca_history": [serialize_nca(n) for n in nca_history],
            "latest_nca": serialize_nca(latest_nca) if latest_nca else None,
            "bioavailability": bioavailability.get("bioavailability", []),
        }

    @app.patch("/api/pk-studies/{study_id}")
    def update_pk_study_endpoint(study_id: int, payload: dict, db: Session = Depends(get_db)):
        study = db.get(PKStudy, study_id)
        if not study:
            raise HTTPException(status_code=404, detail="PK study not found")

        fields = ["study_name", "species", "strain", "sex", "route", "dose", "dose_unit",
                  "formulation", "matrix", "dosing_frequency", "fed_fasted", "lloq", "lloq_unit",
                  "study_date", "source", "notes"]
        for field in fields:
            if field in payload and payload[field] is not None:
                setattr(study, field, payload[field])

        study.dose_normalized_mg_kg = normalize_dose_to_mg_kg(study.dose, study.dose_unit)
        study.updated_at = utcnow()
        db.commit()
        db.refresh(study)
        return serialize_study(study)

    @app.delete("/api/pk-studies/{study_id}", status_code=204)
    def delete_pk_study_endpoint(study_id: int, db: Session = Depends(get_db)):
        study = db.get(PKStudy, study_id)
        if not study:
            raise HTTPException(status_code=404, detail="PK study not found")
        v_id = study.version_id
        db.delete(study)
        db.commit()
        try:
            from .ivive import refresh_pk_and_ivive_for_version
            refresh_pk_and_ivive_for_version(db, v_id, force=True)
        except Exception:
            pass

    @app.post("/api/pk-studies/{study_id}/observations")
    def add_pk_observations_endpoint(study_id: int, payload: list[PKObservationCreate], db: Session = Depends(get_db)):
        study = db.get(PKStudy, study_id)
        if not study:
            raise HTTPException(status_code=404, detail="PK study not found")
        
        mw = (study.version.properties_json or {}).get("molecular_weight")
        added = []
        for item in payload:
            t_hours = normalize_time_to_hours(item.time_raw, item.time_unit)
            c_norm = normalize_conc_to_ng_ml(item.concentration_raw, item.concentration_unit, mw=mw)
            is_blq = item.blq_flag or (item.concentration_raw is None) or (study.lloq is not None and c_norm is not None and c_norm < study.lloq)
            
            obs = PKObservation(
                pk_study_id=study.id,
                version_id=study.version_id,
                subject_group_id=item.subject_group_id.strip() or "Group Mean",
                time_raw=item.time_raw,
                time_unit=item.time_unit.strip(),
                time_hours=t_hours,
                concentration_raw=item.concentration_raw,
                concentration_unit=item.concentration_unit.strip(),
                concentration_normalized_ng_ml=c_norm,
                blq_flag=is_blq,
                replicate=item.replicate.strip() or "R1",
                notes=item.notes.strip(),
            )
            db.add(obs)
            added.append(obs)
        db.commit()
        return {"added_count": len(added), "pk_study_id": study_id}

    @app.delete("/api/pk-observations/{obs_id}", status_code=204)
    def delete_pk_observation_endpoint(obs_id: int, db: Session = Depends(get_db)):
        obs = db.get(PKObservation, obs_id)
        if not obs:
            raise HTTPException(status_code=404, detail="Observation not found")
        db.delete(obs)
        db.commit()

    @app.post("/api/pk-studies/{study_id}/preview-csv")
    async def preview_pk_csv_endpoint(study_id: int, payload: dict | None = None, file: UploadFile | None = File(None), mapping_json: str | None = Query(None), db: Session = Depends(get_db)):
        study = db.get(PKStudy, study_id)
        if not study:
            raise HTTPException(status_code=404, detail="PK study not found")
        content = ""
        mapping = None
        if file:
            content = (await file.read()).decode("utf-8", errors="replace")
            if mapping_json:
                import json
                mapping = json.loads(mapping_json)
        elif payload:
            content = str(payload.get("csv_text") or "")
            mapping = payload.get("mapping")
        if not content.strip():
            raise HTTPException(status_code=400, detail="CSV text or file is required")
        valid_rows, errors, fieldnames = parse_pk_csv(content, mapping)
        return {
            "total_rows": len(valid_rows) + len(errors),
            "valid_count": len(valid_rows),
            "error_count": len(errors),
            "errors": errors,
            "detected_columns": fieldnames,
            "preview_rows": valid_rows[:10],
        }

    @app.post("/api/pk-studies/{study_id}/import-csv")
    async def import_pk_csv_endpoint(study_id: int, payload: dict | None = None, file: UploadFile | None = File(None), mapping_json: str | None = Query(None), db: Session = Depends(get_db)):
        study = db.get(PKStudy, study_id)
        if not study:
            raise HTTPException(status_code=404, detail="PK study not found")
        content = ""
        mapping = None
        if file:
            content = (await file.read()).decode("utf-8", errors="replace")
            if mapping_json:
                import json
                mapping = json.loads(mapping_json)
        elif payload:
            content = str(payload.get("csv_text") or "")
            mapping = payload.get("mapping")
        if not content.strip():
            raise HTTPException(status_code=400, detail="CSV text or file is required")
        valid_rows, errors, _ = parse_pk_csv(content, mapping)
        
        if errors and len(valid_rows) == 0:
            raise HTTPException(status_code=400, detail=f"CSV import failed with {len(errors)} error(s): {errors[0]}")
        
        mw = (study.version.properties_json or {}).get("molecular_weight")
        added = []
        for row in valid_rows:
            t_hours = normalize_time_to_hours(row["time_raw"], row["time_unit"])
            c_norm = normalize_conc_to_ng_ml(row["concentration_raw"], row["concentration_unit"], mw=mw)
            is_blq = row["blq_flag"] or (row["concentration_raw"] is None) or (study.lloq is not None and c_norm is not None and c_norm < study.lloq)
            
            obs = PKObservation(
                pk_study_id=study.id,
                version_id=study.version_id,
                subject_group_id=row["subject"],
                time_raw=row["time_raw"],
                time_unit=row["time_unit"],
                time_hours=t_hours,
                concentration_raw=row["concentration_raw"],
                concentration_unit=row["concentration_unit"],
                concentration_normalized_ng_ml=c_norm,
                blq_flag=is_blq,
                replicate=row["replicate"],
                notes="",
            )
            db.add(obs)
            added.append(obs)
        db.commit()
        try:
            from .ivive import refresh_pk_and_ivive_for_version
            refresh_pk_and_ivive_for_version(db, study.version_id, force=True)
        except Exception:
            pass
        return {"imported_count": len(added), "errors": errors, "pk_study_id": study.id}

    @app.post("/api/pk-studies/{study_id}/run-nca")
    def run_nca_endpoint(study_id: int, options: PKNCARunOptions | None = None, db: Session = Depends(get_db)):
        study = db.get(PKStudy, study_id)
        if not study:
            raise HTTPException(status_code=404, detail="PK study not found")
        
        options = options or PKNCARunOptions()
        observations = list(db.scalars(
            select(PKObservation).where(PKObservation.pk_study_id == study_id)
        ))
        if options.subject_group_id:
            observations = [obs for obs in observations if obs.subject_group_id == options.subject_group_id]

        if not observations:
            raise HTTPException(status_code=400, detail="No observations found in study to run NCA.")

        dose_mg_kg = study.dose_normalized_mg_kg or normalize_dose_to_mg_kg(study.dose, study.dose_unit)
        lloq_val = options.lloq_val if options.lloq_val is not None else study.lloq
        
        nca_dict = run_nca_calculation(
            observations=observations,
            route=study.route,
            dose_mg_kg=dose_mg_kg,
            dose_unit=study.dose_unit,
            blq_policy=options.blq_policy,
            manual_terminal_indices=options.manual_terminal_indices,
            lloq_val=lloq_val,
        )

        # Increment analysis version
        existing_count = db.scalar(
            select(PKNCAResult.id).where(PKNCAResult.pk_study_id == study_id)
        )
        version_num = 1
        if existing_count:
            max_ver = db.scalar(
                select(PKNCAResult.analysis_version).where(PKNCAResult.pk_study_id == study_id).order_by(PKNCAResult.analysis_version.desc()).limit(1)
            )
            version_num = (max_ver or 0) + 1
            
            # Set older runs is_latest to False
            db.execute(
                text("UPDATE pk_nca_results SET is_latest = 0 WHERE pk_study_id = :sid"),
                {"sid": study_id}
            )

        nca_record = PKNCAResult(
            pk_study_id=study.id,
            version_id=study.version_id,
            analysis_version=version_num,
            is_latest=True,
            selection_mode=nca_dict.get("selection_mode", "AUTO"),
            subject_group_id=options.subject_group_id or "Group Mean",
            cmax=nca_dict.get("cmax"),
            cmax_unit=nca_dict.get("cmax_unit", "ng/mL"),
            tmax=nca_dict.get("tmax"),
            tmax_unit=nca_dict.get("tmax_unit", "h"),
            auclast=nca_dict.get("auclast"),
            auclast_unit=nca_dict.get("auclast_unit", "ng*h/mL"),
            aucinf=nca_dict.get("aucinf"),
            aucinf_unit=nca_dict.get("aucinf_unit", "ng*h/mL"),
            lambda_z=nca_dict.get("lambda_z"),
            terminal_half_life=nca_dict.get("terminal_half_life"),
            mrt=nca_dict.get("mrt"),
            cl=nca_dict.get("cl"),
            cl_unit=nca_dict.get("cl_unit", "mL/min/kg"),
            cl_f=nca_dict.get("cl_f"),
            cl_f_unit=nca_dict.get("cl_f_unit", "mL/min/kg"),
            vz=nca_dict.get("vz"),
            vz_unit=nca_dict.get("vz_unit", "L/kg"),
            vz_f=nca_dict.get("vz_f"),
            vz_f_unit=nca_dict.get("vz_f_unit", "L/kg"),
            aumclast=nca_dict.get("aumclast"),
            aumcinf=nca_dict.get("aumcinf"),
            auc_extrapolated_pct=nca_dict.get("auc_extrapolated_pct"),
            terminal_point_count=nca_dict.get("terminal_point_count", 0),
            terminal_points_json=nca_dict.get("terminal_points", []),
            r_squared=nca_dict.get("r_squared"),
            adjusted_r2=nca_dict.get("adjusted_r2"),
            warnings_json=nca_dict.get("warnings", []),
            blq_policy_json=nca_dict.get("blq_policy", {}),
            nca_engine=nca_dict.get("nca_engine", NCA_ENGINE_NAME),
            nca_engine_version=nca_dict.get("nca_engine_version", NCA_ENGINE_VERSION),
            calculation_method=nca_dict.get("calculation_method", "Linear-up / Log-down trapezoidal"),
            provenance_json={
                "compound_version_id": study.version_id,
                "study_id": study.id,
                "species": study.species,
                "route": study.route,
                "dose_mg_kg": dose_mg_kg,
                "dose_raw": study.dose,
                "dose_unit": study.dose_unit,
                "selection_mode": nca_dict.get("selection_mode", "AUTO"),
                "manual_terminal_indices": options.manual_terminal_indices,
                "timestamp": datetime.now(timezone.utc).isoformat(),
            },
        )
        db.add(nca_record)
        db.commit()
        db.refresh(nca_record)
        try:
            from .ivive import refresh_pk_and_ivive_for_version
            refresh_pk_and_ivive_for_version(db, study.version_id, force=True)
        except Exception:
            pass
        return serialize_nca(nca_record)

    @app.get("/api/compounds/{row_id}/bioavailability")
    def get_bioavailability_endpoint(row_id: int, version_id: int | None = Query(None), db: Session = Depends(get_db)):
        compound = db.get(Compound, row_id)
        if not compound:
            raise HTTPException(status_code=404, detail="Compound not found")
        target_version_id = version_id
        if not target_version_id:
            version = next((v for v in compound.versions if v.version_number == compound.current_version), None)
            if version:
                target_version_id = version.id
        if not target_version_id:
            return {"bioavailability": [], "message": "No compound version found"}
        return calculate_bioavailability_for_version(target_version_id, db)

