import csv
import hashlib
import io
from datetime import datetime, timezone

from fastapi import Depends, HTTPException
from fastapi.responses import PlainTextResponse
from sqlalchemy import JSON, DateTime, Float, ForeignKey, Integer, String, Text, UniqueConstraint, select
from sqlalchemy.orm import Mapped, mapped_column, relationship

from .database import Base
from .models import CompoundVersion, utcnow


def ensure_admet_schema(engine):
    from sqlalchemy import inspect
    inspector = inspect(engine)
    if "projects" not in inspector.get_table_names():
        return
    existing = set(inspector.get_table_names())
    if not {"admet_endpoints", "admet_model_registry", "admet_prediction_runs"}.issubset(existing):
        Base.metadata.create_all(
            bind=engine,
            tables=[
                ADMETEndpoint.__table__, ADMETAssayDefinition.__table__, ADMETMeasurement.__table__,
                ADMETModelRegistry.__table__, ADMETPredictionRun.__table__, ADMETPrediction.__table__,
            ],
        )
    with engine.begin() as connection:
        registered = set(connection.execute(select(ADMETModelRegistry.endpoint_name)).scalars())
        for name in ("Solubility", "Permeability", "Microsomal clearance", "Plasma protein binding"):
            if name not in registered:
                connection.execute(
                    ADMETModelRegistry.__table__.insert().values(
                        endpoint_name=name,
                        model_name=f"{name} baseline registry entry",
                    )
                )


ADMET_CSV_COLUMNS = [
    "compound_id", "version_number", "endpoint", "species", "matrix", "value", "unit",
    "qualifier", "replicate", "mean", "sd", "n", "method", "source", "date", "notes",
]


class ADMETEndpoint(Base):
    __tablename__ = "admet_endpoints"
    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    name: Mapped[str] = mapped_column(String(120), index=True)
    category: Mapped[str] = mapped_column(String(80), default="ADME")
    description: Mapped[str] = mapped_column(Text, default="")
    preferred_unit: Mapped[str] = mapped_column(String(40), default="")
    direction: Mapped[str] = mapped_column(String(20), default="lower_better")
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    measurements = relationship("ADMETMeasurement", back_populates="endpoint", cascade="all, delete-orphan")
    predictions = relationship("ADMETPrediction", back_populates="endpoint", cascade="all, delete-orphan")
    __table_args__ = (UniqueConstraint("project_id", "name", name="uq_admet_endpoint_project_name"),)


class ADMETAssayDefinition(Base):
    """Reusable experimental protocol metadata; measurements can override species/matrix."""
    __tablename__ = "admet_assay_definitions"
    id: Mapped[int] = mapped_column(primary_key=True)
    endpoint_id: Mapped[int] = mapped_column(ForeignKey("admet_endpoints.id", ondelete="CASCADE"), index=True)
    assay_uid: Mapped[str] = mapped_column(String(64), unique=True, default=lambda: f"AM-{datetime.now(timezone.utc).timestamp():.0f}")
    method: Mapped[str] = mapped_column(String(200), default="")
    species: Mapped[str] = mapped_column(String(100), default="")
    matrix: Mapped[str] = mapped_column(String(120), default="")
    protocol: Mapped[str] = mapped_column(Text, default="")
    source: Mapped[str] = mapped_column(String(200), default="")
    notes: Mapped[str] = mapped_column(Text, default="")
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    endpoint = relationship("ADMETEndpoint")


class ADMETMeasurement(Base):
    __tablename__ = "admet_measurements"
    id: Mapped[int] = mapped_column(primary_key=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    endpoint_id: Mapped[int] = mapped_column(ForeignKey("admet_endpoints.id", ondelete="CASCADE"), index=True)
    assay_definition_id: Mapped[int | None] = mapped_column(ForeignKey("admet_assay_definitions.id", ondelete="SET NULL"))
    species: Mapped[str] = mapped_column(String(100), default="")
    matrix: Mapped[str] = mapped_column(String(120), default="")
    value: Mapped[float | None] = mapped_column(Float, nullable=True)
    unit: Mapped[str] = mapped_column(String(40))
    qualifier: Mapped[str] = mapped_column(String(5), default="=")
    replicate: Mapped[str] = mapped_column(String(60), default="R1")
    mean_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    standard_deviation: Mapped[float | None] = mapped_column(Float, nullable=True)
    sample_size: Mapped[int | None] = mapped_column(Integer, nullable=True)
    method: Mapped[str] = mapped_column(String(200), default="")
    source: Mapped[str] = mapped_column(String(200), default="User experimental")
    experiment_date: Mapped[str] = mapped_column(String(30), default="")
    notes: Mapped[str] = mapped_column(Text, default="")
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    endpoint = relationship("ADMETEndpoint", back_populates="measurements")
    version = relationship("CompoundVersion")
    assay_definition = relationship("ADMETAssayDefinition")


class ADMETModelRegistry(Base):
    __tablename__ = "admet_model_registry"
    id: Mapped[int] = mapped_column(primary_key=True)
    endpoint_name: Mapped[str] = mapped_column(String(120), index=True)
    model_name: Mapped[str] = mapped_column(String(160))
    model_version: Mapped[str] = mapped_column(String(60), default="0")
    implementation_status: Mapped[str] = mapped_column(String(40), default="NOT_INSTALLED")
    supported_species: Mapped[list] = mapped_column(JSON, default=list)
    supported_matrix: Mapped[list] = mapped_column(JSON, default=list)
    output_unit: Mapped[str] = mapped_column(String(40), default="")
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)
    is_active: Mapped[bool] = mapped_column(default=False)
    registered_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    predictions = relationship("ADMETPrediction", back_populates="model")


class ADMETPrediction(Base):
    __tablename__ = "admet_predictions"
    id: Mapped[int] = mapped_column(primary_key=True)
    run_id: Mapped[int] = mapped_column(ForeignKey("admet_prediction_runs.id", ondelete="CASCADE"), index=True)
    endpoint_id: Mapped[int] = mapped_column(ForeignKey("admet_endpoints.id", ondelete="CASCADE"), index=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    model_id: Mapped[int] = mapped_column(ForeignKey("admet_model_registry.id", ondelete="RESTRICT"))
    predicted_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    unit: Mapped[str] = mapped_column(String(40))
    confidence: Mapped[str] = mapped_column(String(30), default="NOT_AVAILABLE")
    applicability_domain: Mapped[str] = mapped_column(String(40), default="UNKNOWN")
    uncertainty: Mapped[float | None] = mapped_column(Float, nullable=True)
    outputs_json: Mapped[dict] = mapped_column(JSON, default=dict)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    endpoint = relationship("ADMETEndpoint", back_populates="predictions")
    version = relationship("CompoundVersion")
    model = relationship("ADMETModelRegistry", back_populates="predictions")
    run = relationship("ADMETPredictionRun", back_populates="predictions")


class ADMETPredictionRun(Base):
    __tablename__ = "admet_prediction_runs"
    id: Mapped[int] = mapped_column(primary_key=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    requested_by: Mapped[str] = mapped_column(String(160), default="user")
    inputs_hash: Mapped[str] = mapped_column(String(64))
    status: Mapped[str] = mapped_column(String(40), default="NOT_INSTALLED")
    message: Mapped[str] = mapped_column(Text, default="No real ADMET AI/QSAR models are installed in Stage 3.")
    started_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    completed_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    predictions = relationship("ADMETPrediction", back_populates="run", cascade="all, delete-orphan")


def measurement_out(row: ADMETMeasurement):
    return {
        "id": row.id, "version_id": row.version_id, "endpoint_id": row.endpoint_id,
        "species": row.species, "matrix": row.matrix, "value": row.value, "unit": row.unit,
        "qualifier": row.qualifier, "replicate": row.replicate, "mean": row.mean_value,
        "sd": row.standard_deviation, "n": row.sample_size, "method": row.method,
        "source": row.source, "date": row.experiment_date, "notes": row.notes,
        "type": "Experimental", "created_at": row.created_at.isoformat(),
    }


def validate_measurement(payload: dict):
    if not str(payload.get("endpoint") or "").strip():
        raise HTTPException(status_code=400, detail="endpoint is required")
    if not str(payload.get("unit") or "").strip():
        raise HTTPException(status_code=400, detail="unit is required")
    try:
        value = float(payload["value"]) if payload.get("value") not in (None, "") else None
        mean = float(payload.get("mean")) if payload.get("mean") not in (None, "") else None
        sd = float(payload.get("sd")) if payload.get("sd") not in (None, "") else None
        sample_size = int(payload.get("n")) if payload.get("n") not in (None, "") else None
    except (TypeError, ValueError) as exc:
        raise HTTPException(status_code=400, detail=f"Numeric field invalid: {exc}")
    if value is None and mean is None:
        raise HTTPException(status_code=400, detail="Either value or mean is required")
    if sd is not None and sd < 0:
        raise HTTPException(status_code=400, detail="SD cannot be negative")
    if sample_size is not None and sample_size < 1:
        raise HTTPException(status_code=400, detail="Sample size must be at least 1")
    if payload.get("qualifier") not in (None, "", "=", "<", "<=", ">", ">=", "~"):
        raise HTTPException(status_code=400, detail="Qualifier must be one of = < <= > >= ~")
    return value, mean, sd


def parse_csv(text: str):
    reader = csv.DictReader(io.StringIO(text))
    missing = [column for column in ("compound_id", "endpoint", "value", "unit") if column not in (reader.fieldnames or [])]
    if missing:
        raise HTTPException(status_code=400, detail=f"Missing CSV columns: {', '.join(missing)}")
    return list(reader), ADMET_CSV_COLUMNS


def csv_export(rows, labels_by_version=None) -> PlainTextResponse:
    labels_by_version = labels_by_version or {}
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=ADMET_CSV_COLUMNS)
    writer.writeheader()
    for row in rows:
        label, number = labels_by_version.get(row.version_id, ("", ""))
        writer.writerow({
            "compound_id": label, "version_number": number, "endpoint": row.endpoint.name,
            "species": row.species, "matrix": row.matrix, "value": row.value, "unit": row.unit,
            "qualifier": row.qualifier, "replicate": row.replicate, "mean": row.mean_value,
            "sd": row.standard_deviation, "n": row.sample_size, "method": row.method,
            "source": row.source, "date": row.experiment_date, "notes": row.notes.replace("\n", " "),
        })
    return PlainTextResponse(stream.getvalue(), media_type="text/csv")


def inputs_hash(version_ids: list[int]) -> str:
    payload = ",".join(str(value) for value in sorted(version_ids)).encode()
    return hashlib.sha256(payload).hexdigest()
