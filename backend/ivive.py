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
from .pk import PKNCAResult, PKStudy


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


def ensure_ivive_schema(engine):
    inspector = inspect(engine)
    if "projects" not in inspector.get_table_names():
        return
    Base.metadata.create_all(bind=engine, tables=[
        IVIVEInputSet.__table__, PhysiologicalParameterSet.__table__, PhysiologicalParameterOverride.__table__,
        IVIVEMethodRegistry.__table__, IVIVERun.__table__,
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

    def clint_priority(row: dict) -> int:
        source = row["source_type"]
        if source == "EXPERIMENTAL" and row["input_type"] == "RAW_HEPATOCYTE": priority = 10
        elif source == "EXPERIMENTAL" and row["input_type"] == "RAW_MICROSOMAL": priority = 20
        elif source == "EXPERIMENTAL": priority = 25
        elif source == "PROJECT_CALIBRATED": priority = 30
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


def calculate_ivive(db: Session, version: CompoundVersion, species: str, method_key: str = METHOD_KEY) -> IVIVERun:
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
        except (ValueError, IVIVEUnitError) as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
        return _serialize_run(run)
