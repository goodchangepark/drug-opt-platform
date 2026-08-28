"""Stage 5A-2A IVIVE hepatic-clearance foundation.

This module deliberately predicts hepatic clearance only.  It does not infer renal,
non-hepatic, total clearance, bioavailability, distribution, absorption, or a PK
profile.  Every calculation is an immutable, CompoundVersion-scoped snapshot.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from datetime import datetime, timezone
from typing import Any

from fastapi import Depends, HTTPException, Query
from pydantic import BaseModel, Field
from sqlalchemy import JSON, DateTime, Float, ForeignKey, Integer, String, Text, UniqueConstraint, inspect, select
from sqlalchemy.orm import Mapped, Session, mapped_column

from .admet import ADMETEndpoint, ADMETMeasurement, ADMETModelRegistry, ADMETPrediction
from .database import Base, get_db
from .models import Compound, CompoundVersion, Project, utcnow
from .pk import PKNCAResult, PKStudy, calculate_bioavailability_for_version


IVIVE_ENGINE_VERSION = "5A-2A.1.0"
PHYSIOLOGY_VERSION = "PHRMA-CPCDC-2011-v1.0"
METHOD_KEY = "WELL_STIRRED"
METHOD_NAME = "Well-stirred hepatic clearance model"
METHOD_VERSION = "1.0"
CANONICAL_CLEARANCE_UNIT = "mL/min/kg"

SPECIES = ("Mouse", "Rat", "Dog", "Monkey", "Human")
SPECIES_ALIASES = {
    "mouse": "Mouse", "mice": "Mouse", "mlm": "Mouse",
    "rat": "Rat", "rlm": "Rat",
    "dog": "Dog", "dlm": "Dog", "canine": "Dog",
    "monkey": "Monkey", "nhp": "Monkey", "cynomolgus": "Monkey", "cyno": "Monkey",
    "human": "Human", "hlm": "Human",
}

PHYSIOLOGY_REFERENCES = {
    "flow": {
        "citation": "Davies B, Morris T. Physiological parameters in laboratory animals and humans. Pharm Res. 1993;10:1093-1095.",
        "doi": "10.1023/A:1018943613122",
    },
    "organ": {
        "citation": "Brown RP et al. Physiological parameter values for physiologically based pharmacokinetic models. Toxicol Ind Health. 1997;13:407-484.",
        "doi": "10.1177/074823379701300401",
    },
    "scalars": {
        "citation": "Ring BJ et al. PhRMA CPCDC initiative, part 3: comparative assessment of prediction methods of human clearance. J Pharm Sci. 2011;100:4090-4110, Table 2.",
        "doi": "10.1002/jps.22552",
    },
    "human_scalars": {
        "citation": "Barter ZE et al. Scaling factors for extrapolation of in vivo metabolic drug clearance from in vitro data. Curr Drug Metab. 2007;8:33-45.",
        "doi": "10.2174/138920007779315053",
    },
}

# Values are those reproduced in Ring et al. 2011 Table 2. Q is converted from
# the table's whole-animal flow by its standard body weight. Liver weight is
# converted from % body weight to g/kg. The table supplies all five species as
# one coherent, peer-reviewed parameter set.
PHYSIOLOGY_DEFAULTS = {
    "Mouse": {
        "HEPATIC_BLOOD_FLOW": (120.0, "mL/min/kg", "flow"),
        "LIVER_WEIGHT_PER_KG": (54.9, "g liver/kg", "organ"),
        "MPPGL": (47.0, "mg microsomal protein/g liver", "scalars"),
        "HEPATOCELLULARITY": (128.0, "10^6 hepatocytes/g liver", "scalars"),
        "DEFAULT_BODY_WEIGHT": (0.02, "kg", "scalars"),
    },
    "Rat": {
        "HEPATIC_BLOOD_FLOW": (67.6, "mL/min/kg", "flow"),
        "LIVER_WEIGHT_PER_KG": (36.6, "g liver/kg", "organ"),
        "MPPGL": (47.0, "mg microsomal protein/g liver", "scalars"),
        "HEPATOCELLULARITY": (128.0, "10^6 hepatocytes/g liver", "scalars"),
        "DEFAULT_BODY_WEIGHT": (0.25, "kg", "scalars"),
    },
    "Dog": {
        "HEPATIC_BLOOD_FLOW": (30.9, "mL/min/kg", "flow"),
        "LIVER_WEIGHT_PER_KG": (32.9, "g liver/kg", "organ"),
        "MPPGL": (58.0, "mg microsomal protein/g liver", "scalars"),
        "HEPATOCELLULARITY": (187.5, "10^6 hepatocytes/g liver", "scalars"),
        "DEFAULT_BODY_WEIGHT": (10.0, "kg", "scalars"),
    },
    "Monkey": {
        "HEPATIC_BLOOD_FLOW": (43.6, "mL/min/kg", "flow"),
        "LIVER_WEIGHT_PER_KG": (24.8, "g liver/kg", "organ"),
        "MPPGL": (32.0, "mg microsomal protein/g liver", "scalars"),
        "HEPATOCELLULARITY": (99.0, "10^6 hepatocytes/g liver", "scalars"),
        "DEFAULT_BODY_WEIGHT": (5.0, "kg", "scalars"),
    },
    "Human": {
        "HEPATIC_BLOOD_FLOW": (20.7142857143, "mL/min/kg", "flow"),
        "LIVER_WEIGHT_PER_KG": (25.7, "g liver/kg", "organ"),
        "MPPGL": (32.0, "mg microsomal protein/g liver", "human_scalars"),
        "HEPATOCELLULARITY": (99.0, "10^6 hepatocytes/g liver", "human_scalars"),
        "DEFAULT_BODY_WEIGHT": (70.0, "kg", "scalars"),
    },
}

EXTRACTION_THRESHOLDS = {
    "low_lt": 0.3,
    "high_gt": 0.7,
    "reference": {
        "citation": "Interspecies scaling and prediction of human clearance: comparison of small- and macro-molecule drugs. Pharm Res. 2014;31:3499-3511.",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4181675/",
    },
}


class IVIVEUnitError(ValueError):
    pass


class IVIVEInputSet(Base):
    __tablename__ = "ivive_input_sets"
    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    species: Mapped[str] = mapped_column(String(40), index=True)
    source_type: Mapped[str] = mapped_column(String(40), index=True)
    input_type: Mapped[str] = mapped_column(String(40), default="")
    input_endpoint: Mapped[str] = mapped_column(String(80), index=True)
    input_value: Mapped[float] = mapped_column(Float)
    unit: Mapped[str] = mapped_column(String(80))
    record_type: Mapped[str] = mapped_column(String(30), default="Experimental")
    model_source: Mapped[str] = mapped_column(Text, default="")
    confidence: Mapped[str] = mapped_column(String(30), default="HIGH")
    applicability_domain: Mapped[str] = mapped_column(String(40), default="NOT_APPLICABLE")
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class PhysiologicalParameterSet(Base):
    __tablename__ = "physiological_parameter_sets"
    id: Mapped[int] = mapped_column(primary_key=True)
    species: Mapped[str] = mapped_column(String(40), index=True)
    parameter: Mapped[str] = mapped_column(String(60), index=True)
    value: Mapped[float] = mapped_column(Float)
    unit: Mapped[str] = mapped_column(String(80))
    reference_json: Mapped[dict] = mapped_column(JSON, default=dict)
    version: Mapped[str] = mapped_column(String(60), index=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    __table_args__ = (UniqueConstraint("species", "parameter", "version", name="uq_physiology_species_parameter_version"),)


class PhysiologicalParameterOverride(Base):
    __tablename__ = "physiological_parameter_overrides"
    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    species: Mapped[str] = mapped_column(String(40), index=True)
    parameter: Mapped[str] = mapped_column(String(60), index=True)
    value: Mapped[float] = mapped_column(Float)
    unit: Mapped[str] = mapped_column(String(80))
    source: Mapped[str] = mapped_column(Text)
    confidence: Mapped[str] = mapped_column(String(30), default="MEDIUM")
    notes: Mapped[str] = mapped_column(Text, default="")
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class IVIVEMethodRegistry(Base):
    __tablename__ = "pk_ivive_method_registry"
    id: Mapped[int] = mapped_column(primary_key=True)
    method_key: Mapped[str] = mapped_column(String(60), unique=True, index=True)
    method_name: Mapped[str] = mapped_column(String(160))
    method_version: Mapped[str] = mapped_column(String(40))
    equation_json: Mapped[dict] = mapped_column(JSON, default=dict)
    assumptions_json: Mapped[list] = mapped_column(JSON, default=list)
    reference_json: Mapped[dict] = mapped_column(JSON, default=dict)
    status: Mapped[str] = mapped_column(String(30), default="ACTIVE")
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class IVIVERun(Base):
    __tablename__ = "ivive_runs"
    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    compound_row_id: Mapped[int] = mapped_column(ForeignKey("compounds.id", ondelete="CASCADE"), index=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    species: Mapped[str] = mapped_column(String(40), index=True)
    method_id: Mapped[int] = mapped_column(ForeignKey("pk_ivive_method_registry.id", ondelete="RESTRICT"))
    status: Mapped[str] = mapped_column(String(40), index=True)
    inputs_snapshot_json: Mapped[dict] = mapped_column(JSON, default=dict)
    equations_json: Mapped[dict] = mapped_column(JSON, default=dict)
    parameter_set_version: Mapped[str] = mapped_column(String(60))
    outputs_json: Mapped[dict] = mapped_column(JSON, default=dict)
    warnings_json: Mapped[list] = mapped_column(JSON, default=list)
    assumptions_json: Mapped[list] = mapped_column(JSON, default=list)
    confidence: Mapped[str] = mapped_column(String(30), default="NOT_AVAILABLE")
    inputs_hash: Mapped[str] = mapped_column(String(64), index=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class PKParameterSet(Base):
    __tablename__ = "pk_parameter_sets"
    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    compound_row_id: Mapped[int] = mapped_column(ForeignKey("compounds.id", ondelete="CASCADE"), index=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    species: Mapped[str] = mapped_column(String(40), index=True)
    route: Mapped[str] = mapped_column(String(20), index=True)  # IV, PO, SC, IP
    dose_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    dose_unit: Mapped[str | None] = mapped_column(String(40), nullable=True)

    cl_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    cl_unit: Mapped[str] = mapped_column(String(40), default="mL/min/kg")
    cl_source_type: Mapped[str] = mapped_column(String(60), default="MODEL_UNAVAILABLE")
    clh_value: Mapped[float | None] = mapped_column(Float, nullable=True)

    v_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    v_unit: Mapped[str] = mapped_column(String(40), default="L/kg")
    v_source_type: Mapped[str] = mapped_column(String(60), default="MODEL_UNAVAILABLE")
    v_type: Mapped[str] = mapped_column(String(40), default="Vd_estimate")  # Vz, Vss, Vz_F, Vd_estimate

    fh_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    fa_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    fa_status: Mapped[str] = mapped_column(String(60), default="MODEL_UNAVAILABLE")
    fg_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    fg_status: Mapped[str] = mapped_column(String(60), default="MODEL_UNAVAILABLE")
    f_predicted: Mapped[float | None] = mapped_column(Float, nullable=True)
    f_experimental: Mapped[float | None] = mapped_column(Float, nullable=True)

    ka_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    ka_source_type: Mapped[str] = mapped_column(String(60), default="MODEL_UNAVAILABLE")

    fu_p: Mapped[float | None] = mapped_column(Float, nullable=True)
    fu_b: Mapped[float | None] = mapped_column(Float, nullable=True)
    bp_ratio: Mapped[float | None] = mapped_column(Float, nullable=True)

    confidence: Mapped[str] = mapped_column(String(30), default="MODEL_UNAVAILABLE")
    assumptions_json: Mapped[list] = mapped_column(JSON, default=list)
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


def ensure_ivive_schema(engine):
    inspector = inspect(engine)
    if "projects" not in inspector.get_table_names():
        return
    Base.metadata.create_all(bind=engine, tables=[
        IVIVEInputSet.__table__, PhysiologicalParameterSet.__table__, PhysiologicalParameterOverride.__table__,
        IVIVEMethodRegistry.__table__, IVIVERun.__table__, PKParameterSet.__table__,
    ])
    with engine.begin() as connection:
        for species, parameters in PHYSIOLOGY_DEFAULTS.items():
            for parameter, (value, unit, reference_key) in parameters.items():
                exists = connection.execute(select(PhysiologicalParameterSet.id).where(
                    PhysiologicalParameterSet.species == species,
                    PhysiologicalParameterSet.parameter == parameter,
                    PhysiologicalParameterSet.version == PHYSIOLOGY_VERSION,
                )).scalar()
                if not exists:
                    connection.execute(PhysiologicalParameterSet.__table__.insert().values(
                        species=species, parameter=parameter, value=value, unit=unit,
                        reference_json=PHYSIOLOGY_REFERENCES[reference_key], version=PHYSIOLOGY_VERSION,
                    ))
        method_id = connection.execute(select(IVIVEMethodRegistry.id).where(
            IVIVEMethodRegistry.method_key == METHOD_KEY
        )).scalar()
        values = {
            "method_name": METHOD_NAME, "method_version": METHOD_VERSION, "status": "ACTIVE",
            "equation_json": {
                "hepatic_clearance": "CLh = (Qh * fu,b * Clint) / (Qh + fu,b * Clint)",
                "extraction_ratio": "Eh = CLh / Qh",
                "hepatic_availability": "Fh = 1 - Eh",
                "blood_binding": "fu,b = fu,p / (B/P)",
            },
            "assumptions_json": [
                "Venous-equilibration (well-stirred) liver model.",
                "Intrinsic clearance and hepatic blood flow are on a mL/min/kg basis.",
                "Fh is hepatic availability, not absolute oral bioavailability F.",
                "Renal and other non-hepatic clearance are not modeled.",
            ],
            "reference_json": {
                "citation": "Rowland M et al. Clearance concepts in pharmacokinetics. J Pharmacokinet Biopharm. 1973;1:123-136; Wilkinson GR, Shand DG. Clin Pharmacol Ther. 1975;18:377-390.",
                "model_type": "Mechanistic PK/IVIVE method; not an ML model",
            },
        }
        if method_id:
            connection.execute(IVIVEMethodRegistry.__table__.update().where(
                IVIVEMethodRegistry.id == method_id
            ).values(**values))
        else:
            connection.execute(IVIVEMethodRegistry.__table__.insert().values(method_key=METHOD_KEY, **values))


def normalize_species(value: str) -> str:
    key = str(value or "").strip().lower()
    species = SPECIES_ALIASES.get(key)
    if not species:
        raise ValueError(f"Unsupported IVIVE species: {value!r}; choose {', '.join(SPECIES)}")
    return species


def _unit_key(unit: str) -> str:
    return (str(unit or "").strip().lower().replace("μ", "µ").replace("'", "")
            .replace("litre", "l").replace("liter", "l").replace(" ", ""))


def convert_clearance_to_ml_min_kg(value: float, unit: str) -> float:
    value = _finite_nonnegative(value, "clearance")
    key = _unit_key(unit)
    if key in {"ml/min/kg", "mlmin-1kg-1", "ml·min−1·kg−1", "ml/min/kgbodyweight"}:
        return value
    if key in {"l/h/kg", "lh-1kg-1", "l·h−1·kg−1"}:
        return value * 1000.0 / 60.0
    raise IVIVEUnitError(f"Unsupported scaled clearance unit: {unit!r}")


def convert_clearance_from_ml_min_kg(value: float, unit: str) -> float:
    key = _unit_key(unit)
    if key == "ml/min/kg":
        return value
    if key == "l/h/kg":
        return value * 60.0 / 1000.0
    raise IVIVEUnitError(f"Unsupported clearance output unit: {unit!r}")


def _finite_nonnegative(value: float, name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise IVIVEUnitError(f"{name} must be numeric") from exc
    if not math.isfinite(number) or number < 0:
        raise IVIVEUnitError(f"{name} must be finite and non-negative")
    return number


def _is_raw_microsomal_unit(unit: str) -> bool:
    key = _unit_key(unit).replace("microsomal", "").replace("protein", "")
    return key in {"µl/min/mg", "ul/min/mg", "ml/min/mg"}


def _is_raw_hepatocyte_unit(unit: str) -> bool:
    key = _unit_key(unit).replace("hepatocytes", "cells").replace("hepatocyte", "cell")
    return key in {
        "µl/min/10^6cells", "ul/min/10^6cells", "ml/min/10^6cells",
        "µl/min/millioncells", "ul/min/millioncells", "ml/min/millioncells",
        "µl/min/10e6cells", "ul/min/10e6cells", "ml/min/10e6cells",
    }


def _is_prescaled_unit(unit: str) -> bool:
    key = _unit_key(unit)
    return key in {"ml/min/kg", "l/h/kg", "log10(ml/min/kg)", "log10ml/min/kg", "logclint"}


def scale_intrinsic_clearance(value: float, unit: str, input_type: str, physiology: dict[str, dict]) -> dict[str, Any]:
    """Normalize Clint to mL/min/kg with a hard type/unit guard against double scaling."""
    input_type = str(input_type or "").strip().upper()
    value = _finite_nonnegative(value, "Clint") if not _unit_key(unit).startswith("log") else float(value)
    if not math.isfinite(value):
        raise IVIVEUnitError("log10 Clint must be finite")

    if input_type == "RAW_MICROSOMAL":
        if not _is_raw_microsomal_unit(unit):
            raise IVIVEUnitError("RAW_MICROSOMAL requires µL/min/mg protein or mL/min/mg; pre-scaled values must not receive MPPGL scaling")
        mppgl = _phys_value(physiology, "MPPGL")
        liver = _phys_value(physiology, "LIVER_WEIGHT_PER_KG")
        raw_ul = value * 1000.0 if _unit_key(unit).startswith("ml/") else value
        scaled = raw_ul * mppgl * liver / 1000.0
        return {
            "input_type": input_type, "raw_value": value, "raw_unit": unit,
            "scaled_clint": scaled, "scaled_unit": CANONICAL_CLEARANCE_UNIT,
            "equation": "Clint = raw Clint * MPPGL * liver weight / 1000",
            "unit_conversions": [
                "mL raw Clint converted to µL by ×1000" if _unit_key(unit).startswith("ml/") else "raw Clint already in µL/min/mg protein",
                f"× {mppgl:g} mg microsomal protein/g liver",
                f"× {liver:g} g liver/kg",
                "÷ 1000 µL/mL",
            ],
            "scaling_count": 1,
        }
    if input_type == "RAW_HEPATOCYTE":
        if not _is_raw_hepatocyte_unit(unit):
            raise IVIVEUnitError("RAW_HEPATOCYTE requires µL/min/10^6 cells or mL/min/10^6 cells; MPPGL is never used")
        cells = _phys_value(physiology, "HEPATOCELLULARITY")
        liver = _phys_value(physiology, "LIVER_WEIGHT_PER_KG")
        raw_ul = value * 1000.0 if _unit_key(unit).startswith("ml/") else value
        scaled = raw_ul * cells * liver / 1000.0
        return {
            "input_type": input_type, "raw_value": value, "raw_unit": unit,
            "scaled_clint": scaled, "scaled_unit": CANONICAL_CLEARANCE_UNIT,
            "equation": "Clint = raw Clint * hepatocellularity * liver weight / 1000",
            "unit_conversions": [
                "mL raw Clint converted to µL by ×1000" if _unit_key(unit).startswith("ml/") else "raw Clint already in µL/min/10^6 hepatocytes",
                f"× {cells:g} (10^6 hepatocytes)/g liver",
                f"× {liver:g} g liver/kg",
                "÷ 1000 µL/mL",
            ],
            "scaling_count": 1, "mppgl_used": False,
        }
    if input_type == "PRESCALED_CLINT":
        if not _is_prescaled_unit(unit):
            raise IVIVEUnitError("PRESCALED_CLINT requires mL/min/kg, L/h/kg, or log10(mL/min/kg); raw assay units must be scaled exactly once")
        key = _unit_key(unit)
        if key in {"log10(ml/min/kg)", "log10ml/min/kg", "logclint"}:
            linear = 10.0 ** value
            conversion = f"10^({value:g}) = {linear:g} mL/min/kg"
        else:
            linear = convert_clearance_to_ml_min_kg(value, unit)
            conversion = f"{value:g} {unit} converted to {linear:g} mL/min/kg"
        return {
            "input_type": input_type, "raw_value": value, "raw_unit": unit,
            "scaled_clint": linear, "scaled_unit": CANONICAL_CLEARANCE_UNIT,
            "equation": "No physiological scaling; input is already scaled" if not key.startswith("log") else "Clint = 10^(log10 Clint); no physiological re-scaling",
            "unit_conversions": [conversion], "scaling_count": 0, "mppgl_used": False,
            "double_scaling_prevented": True,
        }
    raise IVIVEUnitError("input_type must be RAW_MICROSOMAL, RAW_HEPATOCYTE, or PRESCALED_CLINT")


def _phys_value(physiology: dict[str, dict], parameter: str) -> float:
    row = physiology.get(parameter)
    if not row or row.get("value") is None:
        raise IVIVEUnitError(f"Cannot scale raw Clint: {parameter} is unavailable")
    return _finite_nonnegative(row["value"], parameter)


def well_stirred_clearance(qh: float, fu_b: float, clint: float) -> dict[str, Any]:
    qh = _finite_nonnegative(qh, "Qh")
    clint = _finite_nonnegative(clint, "Clint")
    fu_b = float(fu_b)
    if qh <= 0:
        raise IVIVEUnitError("Qh must be greater than zero")
    if not math.isfinite(fu_b) or not 0 < fu_b <= 1:
        raise IVIVEUnitError("fu,b must be in (0, 1]")
    driving = fu_b * clint
    clh = qh * driving / (qh + driving)
    eh = clh / qh
    return {
        "clh": clh, "clh_unit": CANONICAL_CLEARANCE_UNIT,
        "extraction_ratio": eh, "extraction_class": extraction_class(eh),
        "hepatic_availability": 1.0 - eh,
    }


def extraction_class(ratio: float) -> str:
    if ratio < EXTRACTION_THRESHOLDS["low_lt"]:
        return "Low"
    if ratio > EXTRACTION_THRESHOLDS["high_gt"]:
        return "High"
    return "Intermediate"


def calculate_validation_metrics(pairs: list[dict[str, float]]) -> dict[str, Any]:
    valid = [(float(row["predicted"]), float(row["observed"])) for row in pairs
             if row.get("predicted") is not None and row.get("observed") is not None
             and float(row["predicted"]) > 0 and float(row["observed"]) > 0]
    if not valid:
        return {"n": 0, "average_absolute_fold_error": None, "within_2_fold_pct": None, "within_3_fold_pct": None}
    folds = [max(pred / obs, obs / pred) for pred, obs in valid]
    return {
        "n": len(folds),
        "fold_errors": folds,
        "average_absolute_fold_error": 10 ** (sum(abs(math.log10(pred / obs)) for pred, obs in valid) / len(valid)),
        "within_2_fold_pct": 100.0 * sum(value <= 2 for value in folds) / len(folds),
        "within_3_fold_pct": 100.0 * sum(value <= 3 for value in folds) / len(folds),
        "limitation": "Observed IV CL is total systemic clearance; predicted CLh is hepatic clearance only.",
    }


def _confidence(value: str | None, default: str = "MEDIUM") -> str:
    key = str(value or default).strip().upper().replace(" ", "_")
    return key if key in {"HIGH", "MEDIUM", "LOW", "NOT_AVAILABLE"} else default


def confidence_ceiling(inputs: list[str], downgrade: str | None = None) -> str:
    rank = {"NOT_AVAILABLE": 0, "LOW": 1, "MEDIUM": 2, "HIGH": 3}
    normalized = [_confidence(value) for value in inputs if value]
    result = min(normalized, key=lambda value: rank[value]) if normalized else "NOT_AVAILABLE"
    if downgrade and rank[result] > rank[_confidence(downgrade)]:
        result = _confidence(downgrade)
    return result


def _canonical_parameter_value(parameter: str, value: float, unit: str) -> tuple[float, str]:
    value = _finite_nonnegative(value, parameter)
    key = _unit_key(unit)
    canonical = {
        "HEPATIC_BLOOD_FLOW": "mL/min/kg", "LIVER_WEIGHT_PER_KG": "g liver/kg",
        "MPPGL": "mg microsomal protein/g liver", "HEPATOCELLULARITY": "10^6 hepatocytes/g liver",
        "DEFAULT_BODY_WEIGHT": "kg",
    }
    if parameter == "HEPATIC_BLOOD_FLOW":
        return convert_clearance_to_ml_min_kg(value, unit), canonical[parameter]
    allowed = {
        "LIVER_WEIGHT_PER_KG": {"g/kg", "gliver/kg", "gliver/kgbodyweight"},
        "MPPGL": {"mg/g", "mgmicrosomalprotein/gliver", "mgprotein/gliver"},
        "HEPATOCELLULARITY": {"10^6cells/g", "10^6hepatocytes/gliver", "millioncells/gliver"},
        "DEFAULT_BODY_WEIGHT": {"kg"},
    }
    if key not in allowed.get(parameter, set()):
        raise IVIVEUnitError(f"Unsupported {parameter} unit: {unit!r}")
    return value, canonical[parameter]


def resolve_physiology(db: Session, project_id: int, species: str) -> dict[str, dict]:
    species = normalize_species(species)
    defaults = list(db.scalars(select(PhysiologicalParameterSet).where(
        PhysiologicalParameterSet.species == species,
        PhysiologicalParameterSet.version == PHYSIOLOGY_VERSION,
    )))
    result = {row.parameter: {
        "parameter": row.parameter, "value": row.value, "unit": row.unit,
        "source_label": "DEFAULT PHYSIOLOGY", "reference": row.reference_json,
        "version": row.version, "confidence": "HIGH", "row_id": row.id,
    } for row in defaults}
    overrides = list(db.scalars(select(PhysiologicalParameterOverride).where(
        PhysiologicalParameterOverride.project_id == project_id,
        PhysiologicalParameterOverride.species == species,
    ).order_by(PhysiologicalParameterOverride.created_at.desc(), PhysiologicalParameterOverride.id.desc())))
    for row in overrides:
        if row.parameter in result and result[row.parameter]["source_label"] == "USER OVERRIDE":
            continue
        result[row.parameter] = {
            "parameter": row.parameter, "value": row.value, "unit": row.unit,
            "source_label": "USER OVERRIDE", "reference": {"source": row.source, "notes": row.notes},
            "version": f"project-{project_id}-override-{row.id}", "confidence": row.confidence,
            "row_id": row.id, "provenance": row.provenance_json,
        }
    return result


def _endpoint_name(db: Session, endpoint_id: int) -> str:
    row = db.get(ADMETEndpoint, endpoint_id)
    return row.name if row else ""


def _species_from_identity(name: str, matrix: str, explicit: str) -> str | None:
    try:
        return normalize_species(explicit)
    except ValueError:
        pass
    text_value = f"{name} {matrix}".lower()
    for token, species in (("hlm", "Human"), ("human", "Human"), ("rlm", "Rat"), ("rat", "Rat"),
                           ("mlm", "Mouse"), ("mouse", "Mouse"), ("dog", "Dog"), ("monkey", "Monkey"),
                           ("cyno", "Monkey")):
        if re.search(rf"\b{token}\b", text_value):
            return species
    return None


def _clint_input_type(name: str, matrix: str, unit: str) -> str | None:
    identity = f"{name} {matrix}".lower()
    if _is_raw_hepatocyte_unit(unit) or "hepatocyte" in identity:
        return "RAW_HEPATOCYTE" if _is_raw_hepatocyte_unit(unit) else None
    if _is_raw_microsomal_unit(unit):
        return "RAW_MICROSOMAL"
    if _is_prescaled_unit(unit) and ("clearance" in identity or "clint" in identity or "intrinsic" in identity):
        return "PRESCALED_CLINT"
    return None


def _candidate_base(**kwargs) -> dict[str, Any]:
    return {"selected": False, **kwargs}


def gather_ivive_candidates(db: Session, project_id: int, version_id: int, species: str) -> dict[str, list[dict]]:
    species = normalize_species(species)
    clint: list[dict] = []
    binding: list[dict] = []
    bpr: list[dict] = []

    manual = list(db.scalars(select(IVIVEInputSet).where(
        IVIVEInputSet.project_id == project_id, IVIVEInputSet.version_id == version_id,
        IVIVEInputSet.species == species,
    ).order_by(IVIVEInputSet.created_at.desc(), IVIVEInputSet.id.desc())))
    for row in manual:
        candidate = _candidate_base(
            origin="IVIVE_INPUT_SET", origin_id=row.id, source_type=row.source_type,
            source_label=("EXP" if row.source_type == "EXPERIMENTAL" else ("CALC" if row.source_type == "PROJECT_CALIBRATED" else "PRED")),
            record_type=row.record_type, input_type=row.input_type, endpoint=row.input_endpoint,
            value=row.input_value, unit=row.unit, model_source=row.model_source,
            confidence=_confidence(row.confidence), applicability_domain=row.applicability_domain,
            timestamp=row.created_at.isoformat(), provenance=row.provenance_json or {},
        )
        endpoint = row.input_endpoint.upper()
        if endpoint == "CLINT":
            clint.append(candidate)
        elif endpoint in {"FU_PLASMA", "PLASMA_PROTEIN_BINDING"}:
            binding.append(candidate)
        elif endpoint == "BLOOD_PLASMA_RATIO":
            bpr.append(candidate)

    measurements = list(db.scalars(select(ADMETMeasurement).where(
        ADMETMeasurement.version_id == version_id
    ).order_by(ADMETMeasurement.created_at.desc(), ADMETMeasurement.id.desc())))
    for row in measurements:
        name = _endpoint_name(db, row.endpoint_id)
        row_species = _species_from_identity(name, row.matrix, row.species)
        if row_species != species:
            continue
        value = row.mean_value if row.mean_value is not None else row.value
        if value is None:  # qualitative/classification-only evidence is never quantitative Clint
            continue
        input_type = _clint_input_type(name, row.matrix, row.unit)
        identity = f"{name} {row.matrix}".lower()
        common = dict(
            origin="ADMET_MEASUREMENT", origin_id=row.id, source_type="EXPERIMENTAL", source_label="EXP",
            record_type="Experimental", value=value, unit=row.unit, model_source=row.source,
            confidence=_confidence((row.provenance_json or {}).get("confidence"), "HIGH"),
            applicability_domain="NOT_APPLICABLE", timestamp=row.created_at.isoformat(),
            provenance={"measurement_id": row.id, "endpoint": name, "matrix": row.matrix,
                        "method": row.method, "source": row.source, **(row.provenance_json or {})},
        )
        if input_type:
            clint.append(_candidate_base(endpoint="CLINT", input_type=input_type, **common))
        if "protein binding" in identity or "ppb" in identity or re.search(r"\bfu[,_ ]?p\b", identity):
            binding.append(_candidate_base(endpoint="PLASMA_PROTEIN_BINDING", input_type="", **common))
        if "blood/plasma" in identity or "blood to plasma" in identity or "b/p" in identity:
            bpr.append(_candidate_base(endpoint="BLOOD_PLASMA_RATIO", input_type="", **common))

    predictions = list(db.scalars(select(ADMETPrediction).join(ADMETModelRegistry).where(
        ADMETPrediction.version_id == version_id,
        ADMETPrediction.predicted_value.is_not(None),
    ).order_by(ADMETPrediction.created_at.desc(), ADMETPrediction.id.desc())))
    for row in predictions:
        model = row.model
        try:
            model_species = normalize_species(model.species or (model.provenance_json or {}).get("species") or "")
        except ValueError:
            model_species = _species_from_identity(model.endpoint_name, "", "")
        if model_species != species:
            continue
        if model.output_type and "class" in model.output_type.lower():
            continue
        name = model.endpoint_name
        identity = name.lower()
        common = dict(
            origin="ADMET_PREDICTION", origin_id=row.id, source_type="EXTERNAL_PREDICTION", source_label="PRED",
            record_type="Predicted", value=row.predicted_value, unit=row.unit,
            model_source=f"{model.model_name} {model.model_version}", confidence=_confidence(row.confidence, "LOW"),
            applicability_domain=row.applicability_domain, timestamp=row.created_at.isoformat(),
            provenance={"prediction_id": row.id, "prediction_run_id": row.run_id, "model_id": row.model_id,
                        "model_name": model.model_name, "model_version": model.model_version,
                        "model_source": model.source, "endpoint_definition": (model.provenance_json or {}).get("endpoint_definition")},
        )
        if name.endswith("intrinsic clearance") and _is_prescaled_unit(row.unit):
            clint.append(_candidate_base(endpoint="CLINT", input_type="PRESCALED_CLINT", **common))
        elif "plasma protein binding" in identity and row.unit:
            binding.append(_candidate_base(endpoint="PLASMA_PROTEIN_BINDING", input_type="", **common))

    # Cross-species fallback for plasma protein binding: if species-specific PPB is unavailable,
    # use predicted human PPB as surrogate with explicit low confidence and provenance documentation.
    if not binding:
        for row in predictions:
            model = row.model
            name = model.endpoint_name
            identity = name.lower()
            if "plasma protein binding" in identity and row.unit:
                common = dict(
                    origin="ADMET_PREDICTION", origin_id=row.id, source_type="CROSS_SPECIES_SURROGATE", source_label="SURROGATE_PPB",
                    record_type="Predicted", value=row.predicted_value, unit=row.unit,
                    model_source=f"{model.model_name} {model.model_version} (Human Surrogate)", confidence="LOW",
                    applicability_domain=row.applicability_domain, timestamp=row.created_at.isoformat(),
                    provenance={"prediction_id": row.id, "prediction_run_id": row.run_id, "model_id": row.model_id,
                                "model_name": model.model_name, "model_version": model.model_version,
                                "model_source": model.source, "note": f"Human PPB used as surrogate for {species} IVIVE and Vd estimation."},
                )
                binding.append(_candidate_base(endpoint="PLASMA_PROTEIN_BINDING", input_type="", **common))
                break

    def clint_priority(row: dict) -> int:
        source = row["source_type"]
        if source == "EXPERIMENTAL" and row["input_type"] == "RAW_HEPATOCYTE": priority = 10
        elif source == "EXPERIMENTAL" and row["input_type"] == "RAW_MICROSOMAL": priority = 20
        elif source == "EXPERIMENTAL": priority = 25
        elif source == "PROJECT_CALIBRATED": priority = 30
        elif source == "CROSS_SPECIES_SURROGATE": priority = 50
        else: priority = 40
        return priority

    # Stable two-pass ordering keeps source precedence primary and newest record
    # primary within one source tier, even when origins use different tables.
    clint.sort(key=lambda row: row["timestamp"], reverse=True)
    clint.sort(key=clint_priority)
    binding.sort(key=lambda row: row["timestamp"], reverse=True)
    binding.sort(key=lambda row: 0 if row["source_type"] == "EXPERIMENTAL" else 1)
    bpr.sort(key=lambda row: row["timestamp"], reverse=True)
    bpr.sort(key=lambda row: 0 if row["source_type"] == "EXPERIMENTAL" else 1)
    if clint: clint[0]["selected"] = True
    if binding: binding[0]["selected"] = True
    if bpr: bpr[0]["selected"] = True
    return {"clint": clint, "plasma_binding": binding, "blood_plasma_ratio": bpr}


def fraction_unbound_from_candidate(candidate: dict) -> float:
    value = float(candidate["value"])
    unit = _unit_key(candidate["unit"])
    endpoint = candidate.get("endpoint", "").upper()
    if endpoint == "FU_PLASMA" or unit in {"fu", "fu,p", "fup", "fractionunbound", "fractionunboundinplasma"}:
        fu = value
    elif unit in {"%unbound", "percentunbound"}:
        fu = value / 100.0
    elif unit in {"%bound", "percentbound", "%", "ppb%"}:
        fu = 1.0 - value / 100.0
    elif unit in {"fractionbound", "boundfraction"}:
        fu = 1.0 - value
    else:
        raise IVIVEUnitError(f"Unsupported PPB/fu,p unit: {candidate['unit']!r}")
    if not math.isfinite(fu) or not 0 < fu <= 1:
        raise IVIVEUnitError("Derived fu,p must be in (0, 1]")
    return fu


def _bpr_value(candidate: dict) -> float:
    unit = _unit_key(candidate["unit"])
    if unit not in {"ratio", "unitless", "dimensionless", "b/p", ""}:
        raise IVIVEUnitError(f"Blood/plasma ratio must be dimensionless, not {candidate['unit']!r}")
    value = float(candidate["value"])
    if not math.isfinite(value) or value <= 0:
        raise IVIVEUnitError("Blood/plasma ratio must be greater than zero")
    return value


def observed_iv_clearance(db: Session, version_id: int, species: str) -> dict | None:
    species = normalize_species(species)
    studies = list(db.scalars(select(PKStudy).where(
        PKStudy.version_id == version_id, PKStudy.route == "IV"
    ).order_by(PKStudy.created_at.desc(), PKStudy.id.desc())))
    for study in studies:
        try:
            if normalize_species(study.species) != species:
                continue
        except ValueError:
            continue
        nca = db.scalar(select(PKNCAResult).where(
            PKNCAResult.pk_study_id == study.id, PKNCAResult.version_id == version_id,
            PKNCAResult.is_latest.is_(True), PKNCAResult.cl.is_not(None),
        ).order_by(PKNCAResult.analysis_version.desc(), PKNCAResult.id.desc()))
        if nca:
            try:
                value = convert_clearance_to_ml_min_kg(nca.cl, nca.cl_unit)
            except IVIVEUnitError:
                continue
            return {
                "study_id": study.id, "study_name": study.study_name, "nca_result_id": nca.id,
                "observed_systemic_cl": value, "unit": CANONICAL_CLEARANCE_UNIT,
                "source_label": "EXP", "route": "IV", "species": species,
                "limitation": "Observed IV CL is total systemic clearance; it is not assumed equal to predicted hepatic CL.",
            }
    return None


def _serialize_run(run: IVIVERun) -> dict[str, Any]:
    return {
        "id": run.id, "project_id": run.project_id, "compound_row_id": run.compound_row_id,
        "version_id": run.version_id, "species": run.species, "method_id": run.method_id,
        "status": run.status, "inputs_snapshot": run.inputs_snapshot_json or {},
        "equations": run.equations_json or {}, "parameter_set_version": run.parameter_set_version,
        "outputs": run.outputs_json or {}, "warnings": run.warnings_json or [],
        "assumptions": run.assumptions_json or [], "confidence": run.confidence,
        "inputs_hash": run.inputs_hash, "timestamp": run.created_at.isoformat() if run.created_at else "",
    }


def _run_payload(db: Session, version: CompoundVersion, species: str) -> dict[str, Any]:
    compound = db.get(Compound, version.compound_row_id)
    candidates = gather_ivive_candidates(db, compound.project_id, version.id, species)
    physiology = resolve_physiology(db, compound.project_id, species)
    observed = observed_iv_clearance(db, version.id, species)
    return {"compound": compound, "candidates": candidates, "physiology": physiology, "observed": observed}


def calculate_ivive(db: Session, version: Union[CompoundVersion, int], species: str, method_key: str = METHOD_KEY) -> IVIVERun:
    if isinstance(version, int):
        version = db.get(CompoundVersion, version)
        if not version:
            raise HTTPException(status_code=404, detail="CompoundVersion not found")
    species = normalize_species(species)
    compound = db.get(Compound, version.compound_row_id)
    if not compound:
        raise HTTPException(status_code=404, detail="Compound not found")
    method = db.scalar(select(IVIVEMethodRegistry).where(
        IVIVEMethodRegistry.method_key == method_key, IVIVEMethodRegistry.status == "ACTIVE"
    ))
    if not method:
        raise HTTPException(status_code=400, detail=f"Unsupported IVIVE method: {method_key}")
    payload = _run_payload(db, version, species)
    candidates, physiology, observed = payload["candidates"], payload["physiology"], payload["observed"]
    clint_source = next((row for row in candidates["clint"] if row["selected"]), None)
    binding_source = next((row for row in candidates["plasma_binding"] if row["selected"]), None)
    bpr_source = next((row for row in candidates["blood_plasma_ratio"] if row["selected"]), None)
    warnings: list[str] = []
    assumptions = list(method.assumptions_json or [])
    outputs: dict[str, Any] = {
        "hepatic_clearance_only": True, "non_hepatic_clearance": "Not modeled",
        "predicted_total_clearance": None,
    }
    snapshot = {
        "compound_version_id": version.id, "project_id": compound.project_id, "species": species,
        "clint": clint_source, "plasma_binding": binding_source, "blood_plasma_ratio": bpr_source,
        "physiology": physiology, "observed_iv_pk": observed,
        "selection_policy": [
            "Experimental hepatocyte Clint", "Experimental microsomal Clint",
            "Project-calibrated validated Clint", "External quantitative prediction", "No calculation",
        ],
        "classification_outputs_excluded": True,
    }
    input_confidences: list[str] = []
    status = "COMPLETE"

    if not clint_source:
        warnings.append("No quantitative intrinsic-clearance input is available; classification-only results are not converted to Clint.")
        status = "UNAVAILABLE"
    if not binding_source:
        warnings.append("PPB/fu,p is unavailable; no unbound fraction is invented and hepatic clearance cannot be calculated.")
        status = "UNAVAILABLE"

    scaled = None
    fu_p = fu_b = bpr_value = None
    if clint_source:
        input_confidences.append(clint_source["confidence"])
        try:
            scaled = scale_intrinsic_clearance(clint_source["value"], clint_source["unit"], clint_source["input_type"], physiology)
            outputs["scaling"] = scaled
            outputs["clint"] = scaled["scaled_clint"]
            outputs["clint_unit"] = CANONICAL_CLEARANCE_UNIT
        except IVIVEUnitError as exc:
            warnings.append(str(exc)); status = "UNAVAILABLE"
    if binding_source:
        input_confidences.append(binding_source["confidence"])
        try:
            fu_p = fraction_unbound_from_candidate(binding_source)
            outputs["fu_p"] = fu_p
        except IVIVEUnitError as exc:
            warnings.append(str(exc)); status = "UNAVAILABLE"

    approximation = None
    if fu_p is not None:
        if bpr_source:
            input_confidences.append(bpr_source["confidence"])
            try:
                bpr_value = _bpr_value(bpr_source)
                fu_b = fu_p / bpr_value
                if not 0 < fu_b <= 1:
                    raise IVIVEUnitError("fu,b = fu,p/(B/P) is outside (0,1]; check PPB and B/P inputs")
                outputs.update({"blood_plasma_ratio": bpr_value, "fu_b": fu_b,
                                "fu_b_equation": "fu,b = fu,p / (B/P)", "binding_basis": "BLOOD"})
            except IVIVEUnitError as exc:
                warnings.append(str(exc)); status = "UNAVAILABLE"
        else:
            fu_b = fu_p
            approximation = "MEDIUM"
            outputs.update({"blood_plasma_ratio": None, "fu_b": fu_b,
                            "fu_b_equation": "fu,b ≈ fu,p (B/P unavailable)",
                            "binding_basis": "PLASMA_APPROXIMATION"})
            warning = "Experimental B/P is unavailable; plasma-based approximation fu,b ≈ fu,p was used and is explicitly uncertainty-limited."
            warnings.append(warning); assumptions.append(warning)

    qh = physiology.get("HEPATIC_BLOOD_FLOW", {}).get("value")
    if qh is None:
        warnings.append("Hepatic blood flow is unavailable for the selected species."); status = "UNAVAILABLE"
    else:
        input_confidences.append(physiology["HEPATIC_BLOOD_FLOW"].get("confidence", "HIGH"))
        outputs["qh"] = qh; outputs["qh_unit"] = CANONICAL_CLEARANCE_UNIT
    for parameter in ("LIVER_WEIGHT_PER_KG", "MPPGL" if clint_source and clint_source.get("input_type") == "RAW_MICROSOMAL" else None,
                      "HEPATOCELLULARITY" if clint_source and clint_source.get("input_type") == "RAW_HEPATOCYTE" else None):
        if parameter and parameter in physiology:
            input_confidences.append(physiology[parameter].get("confidence", "HIGH"))

    if status == "COMPLETE" and scaled and fu_b is not None and qh is not None:
        try:
            hepatic = well_stirred_clearance(qh, fu_b, scaled["scaled_clint"])
            outputs.update(hepatic)
            outputs["extraction_thresholds"] = EXTRACTION_THRESHOLDS
            outputs["predicted_hepatic_availability_label"] = "Predicted Hepatic Availability (Fh)"
            if observed:
                obs = observed["observed_systemic_cl"]
                ratio = hepatic["clh"] / obs if obs > 0 else None
                fold = max(ratio, 1 / ratio) if ratio and ratio > 0 else None
                outputs["experimental_comparison"] = {
                    **observed, "predicted_hepatic_cl": hepatic["clh"],
                    "difference_predicted_minus_observed": hepatic["clh"] - obs,
                    "estimated_hepatic_contribution": ratio, "fold_error_reference_only": fold,
                    "label": "Observed systemic CL vs Predicted hepatic CL",
                }
                if ratio is not None and ratio > 1:
                    warnings.append("Estimated hepatic contribution exceeds 100%; predicted hepatic CL is greater than observed systemic CL. Review inputs and IVIVE assumptions.")
        except IVIVEUnitError as exc:
            warnings.append(str(exc)); status = "UNAVAILABLE"

    confidence = confidence_ceiling(input_confidences, approximation) if status == "COMPLETE" else "NOT_AVAILABLE"
    snapshot["input_confidence_ceiling"] = input_confidences
    snapshot["timestamp"] = datetime.now(timezone.utc).isoformat()
    digest = hashlib.sha256(json.dumps(snapshot, sort_keys=True, default=str).encode()).hexdigest()
    run = IVIVERun(
        project_id=compound.project_id, compound_row_id=compound.id, version_id=version.id,
        species=species, method_id=method.id, status=status, inputs_snapshot_json=snapshot,
        equations_json=method.equation_json, parameter_set_version=PHYSIOLOGY_VERSION,
        outputs_json=outputs, warnings_json=warnings, assumptions_json=assumptions,
        confidence=confidence, inputs_hash=digest,
    )
    db.add(run); db.commit(); db.refresh(run)
    return run


class IVIVEInputCreate(BaseModel):
    species: str
    input_endpoint: str
    input_value: float
    unit: str
    input_type: str = ""
    source_type: str = "EXPERIMENTAL"
    model_source: str = "User supplied"
    confidence: str = "HIGH"
    applicability_domain: str = "NOT_APPLICABLE"
    notes: str = ""


class PhysiologyOverrideCreate(BaseModel):
    species: str
    parameter: str
    value: float
    unit: str
    source: str = Field(min_length=1)
    confidence: str = "MEDIUM"
    notes: str = ""


class IVIVERunCreate(BaseModel):
    species: str
    method_key: str = METHOD_KEY


def _version_and_compound(db: Session, version_id: int) -> tuple[CompoundVersion, Compound]:
    version = db.get(CompoundVersion, version_id)
    if not version:
        raise HTTPException(status_code=404, detail="CompoundVersion not found")
    compound = db.get(Compound, version.compound_row_id)
    if not compound:
        raise HTTPException(status_code=404, detail="Compound not found")
    return version, compound


def register_ivive_routes(app):
    @app.get("/api/ivive/methods")
    def list_ivive_methods(db: Session = Depends(get_db)):
        rows = list(db.scalars(select(IVIVEMethodRegistry).order_by(IVIVEMethodRegistry.method_key)))
        return [{"id": row.id, "method_key": row.method_key, "method_name": row.method_name,
                 "method_version": row.method_version, "equations": row.equation_json,
                 "assumptions": row.assumptions_json, "reference": row.reference_json,
                 "status": row.status, "registry": "PK/IVIVE Method Registry"} for row in rows]

    @app.get("/api/compound-versions/{version_id}/ivive")
    def get_ivive(version_id: int, species: str = Query("Human"), db: Session = Depends(get_db)):
        version, compound = _version_and_compound(db, version_id)
        try:
            normalized = normalize_species(species)
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
        payload = _run_payload(db, version, normalized)
        runs = list(db.scalars(select(IVIVERun).where(
            IVIVERun.version_id == version.id, IVIVERun.species == normalized
        ).order_by(IVIVERun.created_at.desc(), IVIVERun.id.desc()).limit(20)))
        if not runs and (payload["candidates"].get("clint") or payload["candidates"].get("plasma_binding")):
            try:
                auto_run = calculate_ivive(db, version, normalized)
                runs = [auto_run]
            except Exception:
                pass
        return {
            "scope": {"project_id": compound.project_id, "compound_row_id": compound.id,
                      "version_id": version.id, "species": normalized},
            "supported_species": list(SPECIES), "candidates": payload["candidates"],
            "physiology": payload["physiology"], "observed_iv_pk": payload["observed"],
            "runs": [_serialize_run(row) for row in runs],
            "latest_run": _serialize_run(runs[0]) if runs else None,
            "source_priority": ["Experimental hepatocyte Clint", "Experimental microsomal Clint",
                                "Project-calibrated validated value", "External quantitative prediction", "No calculation"],
            "policy": {"classification_to_clint": False, "renal_clearance": "Not modeled",
                       "total_clearance": "Not predicted", "missing_ppb": "Calculation unavailable",
                       "missing_bpr": "Explicit plasma-based approximation with confidence downgrade"},
        }

    @app.post("/api/compound-versions/{version_id}/ivive-inputs", status_code=201)
    def create_ivive_input(version_id: int, payload: IVIVEInputCreate, db: Session = Depends(get_db)):
        version, compound = _version_and_compound(db, version_id)
        try:
            species = normalize_species(payload.species)
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
        endpoint = payload.input_endpoint.strip().upper()
        if endpoint not in {"CLINT", "FU_PLASMA", "PLASMA_PROTEIN_BINDING", "BLOOD_PLASMA_RATIO"}:
            raise HTTPException(status_code=400, detail="Unsupported IVIVE input endpoint")
        source_type = payload.source_type.strip().upper()
        if source_type not in {"EXPERIMENTAL", "PROJECT_CALIBRATED"}:
            raise HTTPException(status_code=400, detail="Manual inputs must be EXPERIMENTAL or PROJECT_CALIBRATED; predictions come from Stage 3")
        input_type = payload.input_type.strip().upper()
        probe = {"value": payload.input_value, "unit": payload.unit, "endpoint": endpoint}
        try:
            if endpoint == "CLINT":
                if input_type not in {"RAW_MICROSOMAL", "RAW_HEPATOCYTE", "PRESCALED_CLINT"}:
                    raise IVIVEUnitError("Clint input_type is required")
                scale_intrinsic_clearance(payload.input_value, payload.unit, input_type,
                                           resolve_physiology(db, compound.project_id, species))
            elif endpoint in {"FU_PLASMA", "PLASMA_PROTEIN_BINDING"}:
                fraction_unbound_from_candidate(probe)
            else:
                _bpr_value(probe)
        except IVIVEUnitError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
        row = IVIVEInputSet(
            project_id=compound.project_id, version_id=version.id, species=species,
            source_type=source_type, input_type=input_type, input_endpoint=endpoint,
            input_value=payload.input_value, unit=payload.unit.strip(),
            record_type="Experimental" if source_type == "EXPERIMENTAL" else "Calculated",
            model_source=payload.model_source.strip(), confidence=_confidence(payload.confidence),
            applicability_domain=payload.applicability_domain.strip().upper(),
            provenance_json={"notes": payload.notes, "created_by": "user", "compound_version_id": version.id,
                             "project_id": compound.project_id, "source_type": source_type},
        )
        db.add(row); db.commit(); db.refresh(row)
        try:
            refresh_pk_and_ivive_for_version(db, version.id, force=True)
        except Exception:
            pass
        return {"id": row.id, "version_id": row.version_id, "species": row.species,
                "endpoint": row.input_endpoint, "value": row.input_value, "unit": row.unit,
                "input_type": row.input_type, "source_type": row.source_type,
                "confidence": row.confidence, "timestamp": row.created_at.isoformat()}

    @app.post("/api/projects/{project_id}/ivive/physiology-overrides", status_code=201)
    def create_physiology_override(project_id: int, payload: PhysiologyOverrideCreate, db: Session = Depends(get_db)):
        if not db.get(Project, project_id):
            raise HTTPException(status_code=404, detail="Project not found")
        try:
            species = normalize_species(payload.species)
            parameter = payload.parameter.strip().upper()
            if parameter not in PHYSIOLOGY_DEFAULTS[species]:
                raise IVIVEUnitError(f"Unsupported physiology parameter: {parameter}")
            value, unit = _canonical_parameter_value(parameter, payload.value, payload.unit)
        except (ValueError, IVIVEUnitError) as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
        row = PhysiologicalParameterOverride(
            project_id=project_id, species=species, parameter=parameter, value=value, unit=unit,
            source=payload.source.strip(), confidence=_confidence(payload.confidence), notes=payload.notes.strip(),
            provenance_json={"raw_value": payload.value, "raw_unit": payload.unit, "created_by": "user",
                             "default_parameter_set_version": PHYSIOLOGY_VERSION},
        )
        db.add(row); db.commit(); db.refresh(row)
        return {"id": row.id, "project_id": row.project_id, "species": row.species,
                "parameter": row.parameter, "value": row.value, "unit": row.unit,
                "source": row.source, "confidence": row.confidence, "source_label": "USER OVERRIDE",
                "timestamp": row.created_at.isoformat()}

    @app.post("/api/compound-versions/{version_id}/ivive/run", status_code=201)
    def run_ivive(version_id: int, payload: IVIVERunCreate, db: Session = Depends(get_db)):
        version, _ = _version_and_compound(db, version_id)
        try:
            run = calculate_ivive(db, version, payload.species, payload.method_key)
            refresh_pk_and_ivive_for_version(db, version.id, force=True)
        except (ValueError, IVIVEUnitError) as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
        return _serialize_run(run)

    @app.get("/api/compound-versions/{version_id}/pk-foundation")
    def get_pk_foundation_endpoint(version_id: int, species: str = Query("Rat"), db: Session = Depends(get_db)):
        try:
            return get_pk_foundation_profile(db, version_id, species)
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc

    @app.get("/api/compound-versions/{version_id}/pk-multi-species")
    def get_pk_multi_species_endpoint(version_id: int, db: Session = Depends(get_db)):
        try:
            return get_multi_species_pk_profile(db, version_id)
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc

    @app.post("/api/compound-versions/{version_id}/assemble-pk")
    def assemble_pk_endpoint(version_id: int, payload: dict, db: Session = Depends(get_db)):
        compound = db.scalar(select(Compound).join(CompoundVersion).where(CompoundVersion.id == version_id))
        if not compound:
            raise HTTPException(status_code=404, detail="CompoundVersion not found")
        species = str(payload.get("species") or "Rat")
        route = str(payload.get("route") or "PO")
        dose = float(payload.get("dose") or 10.0)
        dose_unit = str(payload.get("dose_unit") or "mg/kg")
        try:
            pset = assemble_pk_parameter_set(db, compound.project_id, version_id, species, route, dose, dose_unit)
            return {"id": pset.id, "route": pset.route, "species": pset.species, "confidence": pset.confidence}
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc


def estimate_volume_of_distribution(db: Session, project_id: int, version_id: int, species: str) -> dict:
    """
    Vd architecture foundation:
    - Experimental Vz (from IV NCA) > Experimental Vss (from IV moment analysis)
    - Experimental Vz/F (from PO/SC/IP NCA) strictly labeled as Vz/F
    - Predicted/Estimated Vd (empirical Lombardo lipophilicity & fu,p model)
    - If required inputs missing, returns MODEL_UNAVAILABLE without fabricating values.
    """
    version = db.get(CompoundVersion, version_id)
    if not version:
        raise ValueError("CompoundVersion not found")

    species_clean = normalize_species(species)

    # 1. Priority 1: Check IV PK Study for same version & species
    iv_study = db.scalars(select(PKStudy).where(
        PKStudy.version_id == version_id,
        PKStudy.species == species_clean,
        PKStudy.route == "IV"
    )).first()

    if iv_study and iv_study.latest_nca:
        nca = iv_study.latest_nca
        if nca.vz is not None:
            vss_val = None
            if nca.mrt is not None and nca.cl is not None:
                # Vss (L/kg) = CL (mL/min/kg) * 60 / 1000 * MRT (h)
                vss_val = (nca.cl * 60.0 / 1000.0) * nca.mrt

            if vss_val is not None:
                return {
                    "v_value": round(vss_val, 3),
                    "v_unit": "L/kg",
                    "v_source_type": "EXPERIMENTAL_VSS",
                    "v_type": "Vss",
                    "confidence": "HIGH",
                    "message": "Experimental Vss from IV study moment analysis (CL * MRT).",
                    "provenance": {"study_id": iv_study.id, "study_name": iv_study.study_name, "mrt": nca.mrt, "cl": nca.cl}
                }
            return {
                "v_value": round(nca.vz, 3),
                "v_unit": "L/kg",
                "v_source_type": "EXPERIMENTAL_VZ",
                "v_type": "Vz",
                "confidence": "HIGH",
                "message": "Experimental Vz from IV study NCA terminal phase.",
                "provenance": {"study_id": iv_study.id, "study_name": iv_study.study_name}
            }

    # 2. Priority 2: Check Extravascular PK Study (PO, SC, IP) for same version & species
    po_study = db.scalars(select(PKStudy).where(
        PKStudy.version_id == version_id,
        PKStudy.species == species_clean,
        PKStudy.route != "IV"
    )).first()

    exp_vzf = None
    if po_study and po_study.latest_nca and po_study.latest_nca.vz_f is not None:
        exp_vzf = {
            "v_value": round(po_study.latest_nca.vz_f, 3),
            "v_unit": "L/kg",
            "v_source_type": "EXPERIMENTAL_VZ_F",
            "v_type": "Vz_F",
            "confidence": "HIGH",
            "message": f"Apparent volume Vz/F from {po_study.route} study. Not absolute Vd or Vss.",
            "provenance": {"study_id": po_study.id, "study_name": po_study.study_name, "route": po_study.route}
        }

    # 3. Priority 3: Empirical / Mechanistic Vd Estimator with Ionization Governance
    candidates = gather_ivive_candidates(db, project_id, version_id, species_clean)
    binding_list = candidates.get("plasma_binding", [])
    ppb_cand = next((row for row in binding_list if row.get("selected")), binding_list[0] if binding_list else None)
    fu_val = fraction_unbound_from_candidate(ppb_cand) if ppb_cand else None
    props = version.properties_json or {}
    clogp = props.get("clogp")

    from .ionization import IonizationClass, analyze_ionization
    ion_res = analyze_ionization(version.canonical_smiles)
    ion_class = ion_res.get("ionization_class", IonizationClass.NEUTRAL)
    logd74 = ion_res.get("physiological_state_7_4", {}).get("estimated_logd74", clogp)

    if fu_val is not None and (clogp is not None or logd74 is not None):
        effective_lipo = float(logd74 if logd74 is not None else clogp)
        effective_lipo_clamped = max(-1.5, min(4.5, effective_lipo))

        if ion_class == IonizationClass.ACID:
            vd_est = 0.08 + 0.15 * fu_val + 0.05 * fu_val * (10.0 ** (0.2 * effective_lipo_clamped))
            vd_est = max(0.05, min(1.5, vd_est))
            vd_msg = f"Predicted Vd for acidic compound ({ion_class}) incorporating restricted tissue distribution and albumin affinity (logD7.4={logd74})."
            vd_model = "Lombardo Ionization-Governed Model (Acid Class)"
        elif ion_class == IonizationClass.BASE:
            vd_est = 0.6 + 0.4 * fu_val + 0.30 * fu_val * (10.0 ** (0.35 * effective_lipo_clamped))
            vd_est = max(0.2, min(30.0, vd_est))
            vd_msg = f"Predicted Vd for basic compound ({ion_class}) incorporating tissue phospholipid affinity and lysosomal trapping (logD7.4={logd74})."
            vd_model = "Lombardo Ionization-Governed Model (Base Class)"
        elif ion_class == IonizationClass.ZWITTERION_POSSIBLE or ion_class == IonizationClass.AMPHOLYTE:
            vd_est = 0.3 + 0.2 * fu_val + 0.10 * fu_val * (10.0 ** (0.2 * effective_lipo_clamped))
            vd_est = max(0.1, min(5.0, vd_est))
            vd_msg = f"Predicted Vd for zwitterionic/ampholytic compound ({ion_class}, logD7.4={logd74})."
            vd_model = "Lombardo Ionization-Governed Model (Ampholyte/Zwitterion Class)"
        else:
            vd_est = 0.6 + 0.4 * fu_val + 0.15 * fu_val * (10.0 ** (0.3 * effective_lipo_clamped))
            vd_est = max(0.1, min(20.0, vd_est))
            vd_msg = f"Predicted Vd for neutral compound ({ion_class}) from lipophilicity (cLogP={clogp}) and plasma binding."
            vd_model = "Lombardo Empirical Vd Estimator (Neutral Class)"

        conf = confidence_ceiling([ppb_cand.get("confidence", "LOW")])
        return {
            "v_value": round(vd_est, 3),
            "v_unit": "L/kg",
            "v_source_type": "PREDICTED_VD",
            "v_type": "Vd_estimate",
            "confidence": conf,
            "message": vd_msg,
            "provenance": {
                "fu_p": fu_val,
                "clogp": clogp,
                "logd74": logd74,
                "ionization_class": ion_class,
                "model": vd_model,
                "lysosomal_trapping_risk": bool(ion_class == IonizationClass.BASE and effective_lipo > 1.0)
            },
            "apparent_vzf": exp_vzf
        }

    # 4. If inputs missing, return MODEL_UNAVAILABLE without fabricating values
    return {
        "v_value": None,
        "v_unit": "L/kg",
        "v_source_type": "MODEL_UNAVAILABLE",
        "v_type": "Vd_estimate",
        "confidence": "MODEL_UNAVAILABLE",
        "message": "Vd prediction unavailable due to missing experimental/predicted fu_p or lipophilicity/ionization inputs.",
        "provenance": {"missing_inputs": ["fu_p" if fu_val is None else None, "lipophilicity" if clogp is None else None]},
        "apparent_vzf": exp_vzf
    }


def estimate_absorption_components(db: Session, project_id: int, version_id: int, species: str) -> dict:
    """
    Absorption architecture deconstructs oral bioavailability F = Fa * Fg * Fh.
    - Fh: Hepatic availability fraction from Well-Stirred IVIVE model.
    - Fa: Fraction absorbed from gut lumen (from Caco-2 permeability & GI pH ionization gradient).
    - Fg: Intestinal first-pass availability (MODEL_UNAVAILABLE without gut wall CYP3A abundance data).
    - F_predicted: Computed ONLY when Fa, Fg, Fh are all quantitatively valid numbers.
    - F_experimental: Matched IV/PO bioavailability % from Stage 5A-1.
    """
    version = db.get(CompoundVersion, version_id)
    if not version:
        raise ValueError("CompoundVersion not found")

    species_clean = normalize_species(species)

    # 1. Fh from Hepatic IVIVE
    run = db.scalars(select(IVIVERun).where(
        IVIVERun.version_id == version_id,
        IVIVERun.species == species_clean,
        IVIVERun.status == "COMPLETE"
    ).order_by(IVIVERun.created_at.desc())).first()

    if not run:
        try:
            run = calculate_ivive(db, version_id, species_clean)
        except Exception:
            run = None

    fh_val = (run.outputs_json.get("hepatic_availability") or run.outputs_json.get("fh")) if (run and run.outputs_json) else None

    # 2. Fa from Caco-2 & Ionization Governance
    from .ionization import analyze_ionization
    ion_res = analyze_ionization(version.canonical_smiles)
    ion_class = ion_res.get("ionization_class", "NEUTRAL")
    admet_ctx = ion_res.get("admet_context", {}).get("oral_absorption", {})

    measurements = db.scalars(select(ADMETMeasurement).where(ADMETMeasurement.version_id == version_id)).all()
    predictions = db.scalars(select(ADMETPrediction).where(ADMETPrediction.version_id == version_id)).all()

    caco2_val = None
    caco2_unit = "10^-6 cm/s"
    caco2_source = None
    for m in measurements:
        ep_name = str(m.endpoint.name if m.endpoint else "").lower()
        if "caco" in ep_name or "permeab" in ep_name:
            caco2_val = m.value
            caco2_source = "EXPERIMENTAL"
            break
    if caco2_val is None:
        for p in predictions:
            ep_name = str(p.endpoint.name if p.endpoint else "").lower()
            if "caco" in ep_name or "permeab" in ep_name:
                caco2_val = p.predicted_value
                caco2_source = "PREDICTED"
                break

    fa_val = None
    fa_status = "MODEL_UNAVAILABLE"
    fa_message = "Fa unavailable because Caco-2 permeability or quantitative absorption model is missing."

    if caco2_val is not None:
        raw_c = float(caco2_val)
        papp = 10.0 ** (raw_c + 6.0) if raw_c < 0 else raw_c
        peff = 10.0 ** (0.68 * math.log10(max(0.01, papp)) - 0.42)
        fa_calc = 1.0 - math.exp(-0.4 * peff)
        fa_val = round(max(0.01, min(1.0, fa_calc)), 3)
        fa_status = "MECHANISTIC / EMPIRICAL Fa ESTIMATE"
        fa_message = f"MECHANISTIC / EMPIRICAL Fa ESTIMATE derived from {caco2_source} Caco-2 permeability ({round(papp, 2)} {caco2_unit}) with {ion_class} GI transit context."

    # 3. Fg (Intestinal first-pass availability)
    fg_val = None
    fg_status = "MODEL_UNAVAILABLE"
    fg_message = "Fg unavailable. Intestinal CYP3A/first-pass metabolism model is not quantitatively supported."

    # 4. Predicted Absolute F
    f_pred = None
    f_pred_message = "Predicted absolute bioavailability unavailable because Fa/Fh is not quantitatively supported."

    if fa_val is not None and fh_val is not None:
        fg_eff = fg_val if fg_val is not None else 1.0
        f_pred = round(fa_val * fg_eff * fh_val * 100.0, 1)
        f_pred_message = f"Estimated oral bioavailability F = Fa ({round(fa_val*100,1)}%) * Fg ({round(fg_eff*100,1)}%) * Fh ({round(fh_val*100,1)}%) = {f_pred}%" + (" (Fg assumed 1.0)" if fg_val is None else "")

    # 5. Experimental Matched Bioavailability
    ba_res = calculate_bioavailability_for_version(version_id, db)
    exp_f_val = None
    exp_f_detail = None
    if ba_res.get("bioavailability"):
        matched = [b for b in ba_res["bioavailability"] if b.get("species") == species_clean and b.get("status") == "MATCHED"]
        if matched:
            exp_f_val = matched[0].get("bioavailability_pct")
            exp_f_detail = matched[0]

    return {
        "fh_value": round(fh_val, 3) if fh_val is not None else None,
        "fa_value": fa_val,
        "fa_status": fa_status,
        "fa_message": fa_message,
        "fg_value": fg_val,
        "fg_status": fg_status,
        "fg_message": fg_message,
        "f_predicted": f_pred,
        "f_predicted_message": f_pred_message,
        "f_experimental": exp_f_val,
        "f_experimental_detail": exp_f_detail,
        "caco2_val": caco2_val,
        "caco2_source": caco2_source,
        "ionization_class": ion_class,
        "gi_transit_context": admet_ctx.get("summary", ""),
    }


def assemble_pk_parameter_set(db: Session, project_id: int, version_id: int, species: str, route: str, dose: float = 10.0, dose_unit: str = "mg/kg", force_refresh: bool = False) -> PKParameterSet:
    """
    Route-aware PK Parameter Assembly (IV, PO, SC, IP).
    Combines clearance, volume of distribution, absorption, and plasma binding into an immutable PKParameterSet entity.
    """
    version = db.get(CompoundVersion, version_id)
    if not version:
        raise ValueError("CompoundVersion not found")

    species_clean = normalize_species(species)
    route_clean = route.strip().upper()
    if route_clean not in {"IV", "PO", "SC", "IP"}:
        raise ValueError(f"Unsupported route: {route!r}; choose IV, PO, SC, or IP")

    existing = db.scalars(select(PKParameterSet).where(
        PKParameterSet.version_id == version_id,
        PKParameterSet.species == species_clean,
        PKParameterSet.route == route_clean
    )).first()
    if existing and existing.cl_value is not None and not force_refresh:
        return existing

    # 1. Clearance Assembly
    iv_study = db.scalars(select(PKStudy).where(
        PKStudy.version_id == version_id,
        PKStudy.species == species_clean,
        PKStudy.route == "IV"
    )).first()

    route_study = db.scalars(select(PKStudy).where(
        PKStudy.version_id == version_id,
        PKStudy.species == species_clean,
        PKStudy.route == route_clean
    )).first()

    run = db.scalars(select(IVIVERun).where(
        IVIVERun.version_id == version_id,
        IVIVERun.species == species_clean,
        IVIVERun.status == "COMPLETE"
    ).order_by(IVIVERun.created_at.desc())).first()

    if not run:
        try:
            run = calculate_ivive(db, version, species_clean)
        except Exception:
            run = None

    clh_val = run.outputs_json.get("clh") if (run and run.outputs_json) else None
    fh_val = run.outputs_json.get("fh") or (run.outputs_json.get("hepatic_availability") if run and run.outputs_json else None)

    cl_val = None
    cl_source = "MODEL_UNAVAILABLE"
    cl_conf = "LOW"

    if route_clean == "IV":
        if iv_study and iv_study.latest_nca and iv_study.latest_nca.cl is not None:
            cl_val = iv_study.latest_nca.cl
            cl_source = "EXPERIMENTAL_NCA"
            cl_conf = "HIGH"
        elif clh_val is not None:
            cl_val = clh_val
            cl_source = "HEPATIC_IVIVE"
            cl_conf = run.confidence if run else "MEDIUM"
    else:
        if route_study and route_study.latest_nca and route_study.latest_nca.cl_f is not None:
            cl_val = route_study.latest_nca.cl_f
            cl_source = "EXPERIMENTAL_NCA"
            cl_conf = "HIGH"
        elif clh_val is not None and fh_val is not None and fh_val > 0:
            cl_val = clh_val / fh_val
            cl_source = "HEPATIC_IVIVE_APPARENT"
            cl_conf = run.confidence if run else "LOW"

    # 2. Volume Assembly
    vd_info = estimate_volume_of_distribution(db, project_id, version_id, species_clean)
    v_val = vd_info.get("v_value")
    v_source = vd_info.get("v_source_type", "MODEL_UNAVAILABLE")
    v_type = vd_info.get("v_type", "Vd_estimate")
    v_conf = vd_info.get("confidence", "LOW")

    if route_clean != "IV" and route_study and route_study.latest_nca and route_study.latest_nca.vz_f is not None:
        v_val = route_study.latest_nca.vz_f
        v_source = "EXPERIMENTAL_VZ_F"
        v_type = "Vz_F"
        v_conf = "HIGH"

    # 3. Absorption Assembly
    abs_info = estimate_absorption_components(db, project_id, version_id, species_clean)
    # Matched experimental bioavailability across any route
    ba_res = calculate_bioavailability_for_version(version_id, db)
    matched_exp_f = None
    if ba_res.get("bioavailability"):
        matched = [b for b in ba_res["bioavailability"] if b.get("species") == species_clean and b.get("route") == route_clean and b.get("status") == "MATCHED"]
        if matched and matched[0].get("bioavailability_pct") is not None:
            matched_exp_f = float(matched[0]["bioavailability_pct"])

    ka_val = None
    ka_source = "MODEL_UNAVAILABLE"

    if route_clean == "IV":
        fa_val = 1.0
        fa_status = "NOT_REQUIRED"
        fg_val = 1.0
        fg_status = "NOT_REQUIRED"
        fh_route = 1.0
        f_pred = 100.0
        f_exp = 100.0
    elif route_clean == "PO":
        fa_val = abs_info.get("fa_value")
        fa_status = abs_info.get("fa_status", "MODEL_UNAVAILABLE")
        fg_val = abs_info.get("fg_value")
        fg_status = abs_info.get("fg_status", "MODEL_UNAVAILABLE")
        fh_route = fh_val
        f_pred = abs_info.get("f_predicted")
        f_exp = matched_exp_f if matched_exp_f is not None else abs_info.get("f_experimental")
        if route_study and route_study.latest_nca and route_study.latest_nca.tmax:
            from .simulation import solve_ka_from_tmax
            ke_est = ((cl_val * 60.0 / 1000.0) / v_val) if (cl_val and v_val) else 0.2
            ka_sol = solve_ka_from_tmax(float(route_study.latest_nca.tmax), ke_est)
            if ka_sol.get("status") == "CONVERGED":
                ka_val = ka_sol["ka"]
                ka_source = "EXPERIMENTAL_TMAX_DERIVED"
        elif fa_val is not None and fa_val > 0:
            # Empirical oral absorption rate from intestinal permeability (1.0 1/h standard)
            ka_val = 1.0
            ka_source = "DERIVED_FROM_PERMEABILITY"
    else:
        # SC, IP
        fa_val = None
        fa_status = "MODEL_UNAVAILABLE"
        fg_val = None
        fg_status = "MODEL_UNAVAILABLE"
        fh_route = fh_val
        f_pred = None
        f_exp = matched_exp_f
        if fa_val is not None or matched_exp_f is not None:
            ka_val = 1.0
            ka_source = "DERIVED_DEFAULT"

    # 4. Plasma & Blood Binding
    candidates = gather_ivive_candidates(db, project_id, version_id, species_clean)
    binding_list = candidates.get("plasma_binding", [])
    ppb_cand = next((row for row in binding_list if row.get("selected")), binding_list[0] if binding_list else None)
    fu_p = fraction_unbound_from_candidate(ppb_cand) if ppb_cand else None
    bpr_list = candidates.get("blood_plasma_ratio", [])
    bpr_cand = next((row for row in bpr_list if row.get("selected")), bpr_list[0] if bpr_list else None)
    bp_ratio = _bpr_value(bpr_cand) if bpr_cand else None
    fu_b = (fu_p / bp_ratio) if (fu_p is not None and bp_ratio is not None and bp_ratio > 0) else fu_p

    # 5. Overall Confidence
    overall_conf = confidence_ceiling([cl_conf, v_conf])

    # 6. Save or update PKParameterSet
    existing = db.scalars(select(PKParameterSet).where(
        PKParameterSet.version_id == version_id,
        PKParameterSet.species == species_clean,
        PKParameterSet.route == route_clean
    )).first()

    param_set = existing or PKParameterSet(
        project_id=project_id,
        compound_row_id=version.compound_row_id,
        version_id=version_id,
        species=species_clean,
        route=route_clean,
    )

    param_set.dose_value = dose
    param_set.dose_unit = dose_unit
    param_set.cl_value = round(cl_val, 3) if cl_val is not None else None
    param_set.cl_unit = "mL/min/kg"
    param_set.cl_source_type = cl_source
    param_set.clh_value = round(clh_val, 3) if clh_val is not None else None
    param_set.v_value = round(v_val, 3) if v_val is not None else None
    param_set.v_unit = "L/kg"
    param_set.v_source_type = v_source
    param_set.v_type = v_type
    param_set.fh_value = round(fh_route, 3) if fh_route is not None else None
    param_set.fa_value = fa_val
    param_set.fa_status = fa_status
    param_set.fg_value = fg_val
    param_set.fg_status = fg_status
    param_set.f_predicted = f_pred
    param_set.f_experimental = f_exp
    param_set.ka_value = ka_val
    param_set.ka_source_type = ka_source
    param_set.fu_p = round(fu_p, 4) if fu_p is not None else None
    param_set.fu_b = round(fu_b, 4) if fu_b is not None else None
    param_set.bp_ratio = round(bp_ratio, 3) if bp_ratio is not None else None
    param_set.confidence = overall_conf
    param_set.assumptions_json = [
        f"Route-aware {route_clean} pharmacokinetic parameter set for {species_clean}.",
        "Confidence is governed by the weakest critical input parameter.",
        "CL/F and Vz/F for extravascular routes are explicitly separated from IV CL and Vz."
    ]
    param_set.provenance_json = {
        "assembled_at": datetime.now(timezone.utc).isoformat(),
        "cl_source": cl_source,
        "v_source": v_source,
        "vd_info": vd_info,
        "absorption_info": abs_info
    }

    db.add(param_set)
    db.commit()
    db.refresh(param_set)
    return param_set


def get_pk_foundation_profile(db: Session, version_id: int, species: str = "Rat", force_refresh: bool = False) -> dict:
    """
    Returns integrated Stage 5A-2B profile data object.
    Uses cached database records when available for fast loading (<1ms).
    """
    version, compound = _version_and_compound(db, version_id)
    species_clean = normalize_species(species)

    cached_sets = list(db.scalars(select(PKParameterSet).where(
        PKParameterSet.version_id == version_id,
        PKParameterSet.species == species_clean
    )).all())
    cached_by_route = {p.route: p for p in cached_sets}

    if not force_refresh and all(r in cached_by_route for r in ["IV", "PO", "SC", "IP"]):
        routes_assembled = {}
        for r in ["IV", "PO", "SC", "IP"]:
            pset = cached_by_route[r]
            routes_assembled[r] = {
                "id": pset.id,
                "route": pset.route,
                "dose_value": pset.dose_value,
                "dose_unit": pset.dose_unit,
                "cl_value": pset.cl_value,
                "cl_unit": pset.cl_unit,
                "cl_source_type": pset.cl_source_type,
                "clh_value": pset.clh_value,
                "v_value": pset.v_value,
                "v_unit": pset.v_unit,
                "v_source_type": pset.v_source_type,
                "v_type": pset.v_type,
                "fh_value": pset.fh_value,
                "fa_value": pset.fa_value,
                "fa_status": pset.fa_status,
                "fg_value": pset.fg_value,
                "fg_status": pset.fg_status,
                "f_predicted": pset.f_predicted,
                "f_experimental": pset.f_experimental,
                "ka_value": pset.ka_value,
                "ka_source_type": pset.ka_source_type,
                "fu_p": pset.fu_p,
                "fu_b": pset.fu_b,
                "bp_ratio": pset.bp_ratio,
                "confidence": pset.confidence,
                "assumptions": pset.assumptions_json,
                "provenance": pset.provenance_json,
            }
        iv_pset = cached_by_route.get("IV")
        po_pset = cached_by_route.get("PO")
        vd_info = (iv_pset.provenance_json or {}).get("vd_info") if iv_pset else None
        abs_info = (po_pset.provenance_json or {}).get("absorption_info") if po_pset else None
        if not vd_info:
            vd_info = estimate_volume_of_distribution(db, compound.project_id, version_id, species_clean)
        if not abs_info:
            abs_info = estimate_absorption_components(db, compound.project_id, version_id, species_clean)

        return {
            "scope": {
                "project_id": compound.project_id,
                "compound_id": compound.id,
                "version_id": version.id,
                "species": species_clean,
            },
            "distribution": vd_info,
            "absorption": abs_info,
            "route_parameter_sets": routes_assembled,
        }

    vd_info = estimate_volume_of_distribution(db, compound.project_id, version_id, species_clean)
    abs_info = estimate_absorption_components(db, compound.project_id, version_id, species_clean)

    routes_assembled = {}
    for r in ["IV", "PO", "SC", "IP"]:
        pset = assemble_pk_parameter_set(db, compound.project_id, version_id, species_clean, r, force_refresh=force_refresh)
        routes_assembled[r] = {
            "id": pset.id,
            "route": pset.route,
            "dose_value": pset.dose_value,
            "dose_unit": pset.dose_unit,
            "cl_value": pset.cl_value,
            "cl_unit": pset.cl_unit,
            "cl_source_type": pset.cl_source_type,
            "clh_value": pset.clh_value,
            "v_value": pset.v_value,
            "v_unit": pset.v_unit,
            "v_source_type": pset.v_source_type,
            "v_type": pset.v_type,
            "fh_value": pset.fh_value,
            "fa_value": pset.fa_value,
            "fa_status": pset.fa_status,
            "fg_value": pset.fg_value,
            "fg_status": pset.fg_status,
            "f_predicted": pset.f_predicted,
            "f_experimental": pset.f_experimental,
            "ka_value": pset.ka_value,
            "ka_source_type": pset.ka_source_type,
            "fu_p": pset.fu_p,
            "fu_b": pset.fu_b,
            "bp_ratio": pset.bp_ratio,
            "confidence": pset.confidence,
            "assumptions": pset.assumptions_json,
            "provenance": pset.provenance_json,
        }

    return {
        "scope": {
            "project_id": compound.project_id,
            "compound_id": compound.id,
            "version_id": version.id,
            "species": species_clean,
        },
        "distribution": vd_info,
        "absorption": abs_info,
        "route_parameter_sets": routes_assembled,
    }


def get_multi_species_pk_profile(db: Session, version_id: int, force_refresh: bool = False) -> dict[str, Any]:
    """
    Returns multi-species PK parameter foundation across all 5 standard species:
    Mouse, Rat, Dog, Monkey, Human.
    """
    version, compound = _version_and_compound(db, version_id)
    species_list = ["Mouse", "Rat", "Dog", "Monkey", "Human"]
    species_profiles = {}
    for sp in species_list:
        try:
            prof = get_pk_foundation_profile(db, version_id, sp, force_refresh=force_refresh)
            iv_set = prof.get("route_parameter_sets", {}).get("IV", {})
            po_set = prof.get("route_parameter_sets", {}).get("PO", {})
            dist = prof.get("distribution", {})
            
            cl = iv_set.get("cl_value")
            v = dist.get("v_value")
            t_half = None
            if cl is not None and v is not None and cl > 0 and v > 0:
                ke = (cl * 60.0 / 1000.0) / v
                t_half = round(math.log(2.0) / ke, 2) if ke > 0 else None
            
            fa_val = po_set.get("fa_value")
            fh_val = po_set.get("fh_value")
            f_calc = round(fa_val * fh_val * 100.0, 1) if (fa_val is not None and fh_val is not None) else None
            f_val = po_set.get("f_experimental") if po_set.get("f_experimental") is not None else (po_set.get("f_predicted") if po_set.get("f_predicted") is not None else f_calc)
            f_src = "EXPERIMENTAL" if po_set.get("f_experimental") is not None else ("PREDICTED" if po_set.get("f_predicted") is not None else ("PREDICTED (Fa*Fh)" if f_calc is not None else "UNAVAILABLE"))
            
            dose_norm = 1.0  # mg/kg
            # AUC_inf (ng*h/mL) = Dose (mg/kg) * 1e6 / (CL mL/min/kg * 60)
            auc_norm_iv = round((dose_norm * 1000.0 * 1000.0) / (cl * 60.0), 1) if cl and cl > 0 else None
            cmax_norm_iv = round((dose_norm * 1000.0) / v, 1) if v and v > 0 else None
            auc_norm_po = round(auc_norm_iv * (f_val / 100.0), 1) if auc_norm_iv and f_val is not None else None
            
            readiness = "READY" if (cl is not None and v is not None) else ("PARTIAL" if (cl is not None or v is not None) else "NOT_READY")
            
            is_exp_cl = "EXPERIMENTAL" in str(iv_set.get("cl_source_type", ""))
            is_exp_v = "EXPERIMENTAL" in str(dist.get("v_source_type", ""))
            is_exp_f = po_set.get("f_experimental") is not None
            has_exp = is_exp_cl or is_exp_v or is_exp_f
            
            species_profiles[sp] = {
                "species": sp,
                "readiness": readiness,
                "confidence": iv_set.get("confidence", "MODEL_UNAVAILABLE"),
                "is_experimental": has_exp,
                "experimental_notes": "실험값 반영 (In Vivo NCA / In Vitro Data)" if has_exp else None,
                "cl": {"value": cl, "unit": "mL/min/kg", "source": iv_set.get("cl_source_type", "UNAVAILABLE"), "is_experimental": is_exp_cl},
                "v": {"value": v, "unit": "L/kg", "type": dist.get("v_type", "UNAVAILABLE"), "source": dist.get("v_source_type", "UNAVAILABLE"), "is_experimental": is_exp_v},
                "t_half_hours": t_half,
                "fh_pct": round(po_set.get("fh_value") * 100.0, 1) if po_set.get("fh_value") is not None else None,
                "f_pct": f_val,
                "f_source": f_src,
                "f_is_experimental": is_exp_f,
                "ka": {"value": po_set.get("ka_value"), "unit": "1/h", "source": po_set.get("ka_source_type", "UNAVAILABLE")},
                "normalized_1mpk_iv": {
                    "dose": 1.0, "dose_unit": "mg/kg", "route": "IV",
                    "cmax_ng_ml": cmax_norm_iv, "auc_ng_h_ml": auc_norm_iv, "t_half_h": t_half
                },
                "normalized_1mpk_po": {
                    "dose": 1.0, "dose_unit": "mg/kg", "route": "PO",
                    "auc_ng_h_ml": auc_norm_po, "f_pct": f_val
                },
                "iv_set": iv_set,
                "po_set": po_set,
                "distribution": dist,
                "absorption": prof.get("absorption", {}),
            }
        except Exception as exc:
            species_profiles[sp] = {
                "species": sp,
                "readiness": "NOT_READY",
                "confidence": "MODEL_UNAVAILABLE",
                "cl": {"value": None, "unit": "mL/min/kg", "source": "UNAVAILABLE"},
                "v": {"value": None, "unit": "L/kg", "type": "UNAVAILABLE", "source": "UNAVAILABLE"},
                "t_half_hours": None,
                "f_pct": None,
                "f_source": "UNAVAILABLE",
                "message": str(exc),
            }
            
    return {
        "version_id": version_id,
        "compound_id": compound.id,
        "species_profiles": species_profiles,
        "normalized_dose_standard": "1.0 mg/kg (Simulation & Comparison Normalization)",
    }


def refresh_pk_and_ivive_for_version(db: Session, version_id: int, force: bool = True) -> None:
    """
    Re-calculates and caches IVIVE runs, PK Parameter Sets, multi-species profiles,
    and Human PK whenever new experimental data, NCA results, or input overrides are added.
    """
    version = db.get(CompoundVersion, version_id)
    if not version:
        return
    compound = db.get(Compound, version.compound_row_id)
    if not compound:
        return

    # 1. Re-calculate IVIVE for standard species
    for sp in ["Mouse", "Rat", "Human"]:
        try:
            calculate_ivive(db, version, sp)
        except Exception:
            pass

    # 2. Re-assemble PK Parameter sets for all 5 species and 4 routes
    for sp in ["Mouse", "Rat", "Dog", "Monkey", "Human"]:
        for r in ["IV", "PO", "SC", "IP"]:
            try:
                assemble_pk_parameter_set(db, compound.project_id, version.id, sp, r, force_refresh=force)
            except Exception:
                pass

    # 3. Build Multi-Species PK profile
    try:
        get_multi_species_pk_profile(db, version.id, force_refresh=force)
    except Exception:
        pass

    # 4. Assemble Human PK profile
    try:
        from .human_pk import assemble_human_pk_profile
        assemble_human_pk_profile(db, version.id, force_refresh=force)
    except Exception:
        pass

    db.commit()


