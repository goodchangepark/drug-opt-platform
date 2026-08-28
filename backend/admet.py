import csv
import hashlib
import io
from datetime import datetime, timezone

from fastapi import Depends, HTTPException
from fastapi.responses import PlainTextResponse
from sqlalchemy import JSON, DateTime, Float, ForeignKey, Integer, String, Text, UniqueConstraint, inspect, select, text
from sqlalchemy.orm import Mapped, mapped_column, relationship

from .database import Base
from .models import CompoundVersion, utcnow
from .admet_predictor import MODEL_SPECS, registry_seed


TRANSPORTER_UNAVAILABLE = {
    "P-gp substrate": {
        "transporter": "P-gp / ABCB1", "role": "SUBSTRATE", "species": "Human",
        "reason": "Distinct public substrate models were identified, but no checkpoint with sufficiently documented validation, reusable local weights, and clear redistribution terms qualified; the P-gp inhibitor checkpoint is never reused.",
    },
    "BCRP substrate": {
        "transporter": "BCRP / ABCG2", "role": "SUBSTRATE", "species": "Human",
        "reason": "No scientifically qualified public, locally deployable BCRP substrate checkpoint with documented assay provenance and validation was found.",
    },
    "BCRP inhibitor": {
        "transporter": "BCRP / ABCG2", "role": "INHIBITOR", "species": "Human",
        "reason": "Web predictors and datasets exist, but no scientifically qualified public checkpoint with clear local redistribution terms was found.",
    },
    "BSEP inhibitor": {
        "transporter": "BSEP / ABCB11", "role": "INHIBITOR", "species": "Human",
        "reason": "Public web models/datasets exist, but no scientifically qualified reusable local checkpoint with clear license and validation provenance was found.",
    },
    "OATP1B1 inhibitor": {
        "transporter": "OATP1B1 / SLCO1B1", "role": "INHIBITOR", "species": "Human",
        "reason": "No scientifically qualified public, locally deployable OATP1B1 inhibitor checkpoint with clear redistribution terms was found.",
    },
    "OATP1B3 inhibitor": {
        "transporter": "OATP1B3 / SLCO1B3", "role": "INHIBITOR", "species": "Human",
        "reason": "No scientifically qualified public, locally deployable OATP1B3 inhibitor checkpoint with clear redistribution terms was found.",
    },
    "OCT1 inhibitor": {
        "transporter": "OCT1 / SLC22A1", "role": "INHIBITOR", "species": "Human",
        "reason": "The available public data are insufficient to qualify a licensed, validated local OCT1 inhibitor checkpoint.",
    },
    "OCT2 inhibitor": {
        "transporter": "OCT2 / SLC22A2", "role": "INHIBITOR", "species": "Human",
        "reason": "The available public data are insufficient to qualify a licensed, validated local OCT2 inhibitor checkpoint.",
    },
    "MATE1 inhibitor": {
        "transporter": "MATE1 / SLC47A1", "role": "INHIBITOR", "species": "Human",
        "reason": "No scientifically qualified public, locally deployable MATE1 inhibitor checkpoint with clear redistribution terms was found.",
    },
    "MATE2-K inhibitor": {
        "transporter": "MATE2-K / SLC47A2", "role": "INHIBITOR", "species": "Human",
        "reason": "No scientifically qualified public, locally deployable MATE2-K inhibitor checkpoint with clear redistribution terms was found.",
    },
}

SAFETY_UNAVAILABLE = {
    "Mitochondrial toxicity": {
        "safety_endpoint": "Mitochondrial toxicity", "species": "Human",
        "reason": "No public checkpoint with a sufficiently specific assay definition, qualified validation, and clear local redistribution terms was identified in Stage 3F.",
    },
    "General cytotoxicity": {
        "safety_endpoint": "General cytotoxicity", "species": "Not standardized",
        "reason": "Cell-line, exposure-time, and assay heterogeneity prevents qualification as one deployable endpoint; no generic cytotoxicity prediction is emitted.",
    },
    "Skin sensitization": {
        "safety_endpoint": "Skin sensitization", "species": "Not standardized",
        "reason": "No validated public checkpoint with sufficiently clear assay and license provenance was qualified for local deployment.",
    },
    "BBB penetration": {
        "safety_endpoint": "BBB penetration", "species": "Human",
        "reason": "BBB is a distribution endpoint rather than a toxicity endpoint, and no model was added merely to increase endpoint count.",
    },
    "CNS liability": {
        "safety_endpoint": "CNS liability", "species": "Human",
        "reason": "CNS liability is not a single assay-defined endpoint; no composite or unsupported prediction is emitted.",
    },
}

PHYSICOCHEM_UNAVAILABLE = {
    "pKa (quantitative ML)": {
        "endpoint": "pKa (quantitative ML)", "species": "Chemical / In Vitro",
        "reason": "Evaluated candidate models (MolGpKa, Epik, Chemprop) were rejected due to ARM64 binary build failures or proprietary licensing. Deterministic structural classification (IonizationEngine_v1) is active; quantitative ML pKa is marked MODEL_UNAVAILABLE to prevent fabricated values.",
    },
    "logD7.4 (quantitative ML)": {
        "endpoint": "logD7.4 (quantitative ML)", "species": "Chemical / In Vitro",
        "reason": "No qualified standalone open-source logD7.4 checkpoint was identified in Stage 4C-4. Calculated cLogP is strictly distinguished and never equated to logD7.4 without explicit ionization assumptions.",
    },
}


def ensure_admet_schema(engine):
    inspector = inspect(engine)
    if "projects" not in inspector.get_table_names():
        return
    Base.metadata.create_all(
        bind=engine,
        tables=[
            ADMETEndpoint.__table__, ADMETAssayDefinition.__table__, ADMETMeasurement.__table__,
            ADMETModelRegistry.__table__, ADMETPredictionRun.__table__, ADMETPrediction.__table__,
            ADMETConsensusPrediction.__table__, ADMETModelComparison.__table__, ADMETModelPerformance.__table__,
        ],
    )
    with engine.begin() as connection:
        measurement_columns = {row["name"] for row in inspect(engine).get_columns("admet_measurements")}
        if "qualitative_value" not in measurement_columns:
            connection.execute(text("ALTER TABLE admet_measurements ADD COLUMN qualitative_value VARCHAR(120) NOT NULL DEFAULT ''"))
        registry_columns = {row["name"] for row in inspect(engine).get_columns("admet_model_registry")}
        registry_additions = {
            "source": "TEXT NOT NULL DEFAULT ''", "training_dataset": "TEXT NOT NULL DEFAULT ''",
            "validation_json": "JSON NOT NULL DEFAULT '{}'", "license": "TEXT NOT NULL DEFAULT ''",
            "model_priority": "INTEGER NOT NULL DEFAULT 100", "ensemble_eligible": "BOOLEAN NOT NULL DEFAULT 1",
            "species": "VARCHAR(100) NOT NULL DEFAULT ''", "output_type": "VARCHAR(60) NOT NULL DEFAULT ''",
        }
        prediction_columns = {row["name"] for row in inspect(engine).get_columns("admet_predictions")}
        prediction_additions = {
            "model_version": "VARCHAR(60) NOT NULL DEFAULT '0'",
            "execution_status": "VARCHAR(40) NOT NULL DEFAULT 'SUCCESS'",
            "standardizer_version": "VARCHAR(40) NOT NULL DEFAULT 'CHEM_STANDARDIZER_V1'",
            "canonical_smiles": "TEXT NOT NULL DEFAULT ''",
            "runtime_ms": "FLOAT NOT NULL DEFAULT 0.0",
        }
        for column, definition in prediction_additions.items():
            if column not in prediction_columns:
                connection.execute(text(f"ALTER TABLE admet_predictions ADD COLUMN {column} {definition}"))

        consensus_columns = {row["name"] for row in inspect(engine).get_columns("admet_consensus_predictions")}
        consensus_additions = {
            "consensus_version": "VARCHAR(60) NOT NULL DEFAULT 'stage4d1-static-v1'",
            "consensus_mode": "VARCHAR(30) NOT NULL DEFAULT 'SHADOW'",
            "model_agreement": "VARCHAR(40) NOT NULL DEFAULT 'HIGH_AGREEMENT'",
            "dispersion_json": "JSON NOT NULL DEFAULT '{}'",
            "vote_pattern": "VARCHAR(120) NOT NULL DEFAULT ''",
        }
        for column, definition in consensus_additions.items():
            if column not in consensus_columns:
                connection.execute(text(f"ALTER TABLE admet_consensus_predictions ADD COLUMN {column} {definition}"))

        registered = set(connection.execute(select(ADMETModelRegistry.endpoint_name)).scalars())
        registry_names = list(MODEL_SPECS) + [
            "Microsomal clearance", "Dog liver microsomal intrinsic clearance",
            "Monkey liver microsomal intrinsic clearance",
            "CYP1A2 substrate", "CYP2C19 substrate",
        ] + list(TRANSPORTER_UNAVAILABLE) + list(SAFETY_UNAVAILABLE) + list(PHYSICOCHEM_UNAVAILABLE)
        for name in registry_names:
            if name not in registered:
                values = registry_seed(name) if name in MODEL_SPECS else {
                    "endpoint_name": name,
                    "model_name": f"{name} — no scientifically qualified model installed",
                    "model_version": (
                        "unavailable-stage4c4" if name in PHYSICOCHEM_UNAVAILABLE else
                        ("unavailable-stage3f" if name in SAFETY_UNAVAILABLE else
                        ("unavailable-stage3e" if name in TRANSPORTER_UNAVAILABLE else
                        ("unavailable-stage3c" if name.startswith("CYP") else "unavailable-stage3b")
                        ))), "implementation_status": "MODEL_UNAVAILABLE",
                    "is_active": False,
                    "provenance_json": (
                        {**PHYSICOCHEM_UNAVAILABLE[name], "status": "MODEL_UNAVAILABLE", "checkpoint_available": False,
                         "model_source": "Stage 4C-4 pKa / logD model qualification audit",
                         "license": "No qualified model/checkpoint license"}
                        if name in PHYSICOCHEM_UNAVAILABLE else
                        ({**SAFETY_UNAVAILABLE[name], "status": "MODEL_UNAVAILABLE", "checkpoint_available": False,
                          "model_source": "Stage 3F public-model qualification audit",
                          "license": "No qualified model/checkpoint license"}
                        if name in SAFETY_UNAVAILABLE else
                        ({**TRANSPORTER_UNAVAILABLE[name], "status": "MODEL_UNAVAILABLE", "checkpoint_available": False}
                        if name in TRANSPORTER_UNAVAILABLE else
                        {"reason": (
                            "A released upstream checkpoint exists, but CYP1A2/CYP2C19 substrate endpoints lack publisher-reported validation and sufficiently clear dataset/assay provenance; disabled rather than emitting an unsupported prediction."
                            if name in {"CYP1A2 substrate", "CYP2C19 substrate"} else
                            "No endpoint- and species-specific public pretrained model qualified in Stage 3B; no cross-species reuse or fake prediction."
                        )}))
                    ),
                }
                connection.execute(
                    ADMETModelRegistry.__table__.insert().values(**values)
                )
            elif name in MODEL_SPECS:
                values = registry_seed(name)
                canonical = connection.execute(
                    select(ADMETModelRegistry.id).where(
                        ADMETModelRegistry.endpoint_name == name,
                        ADMETModelRegistry.model_name == values["model_name"],
                        ADMETModelRegistry.model_version == values["model_version"],
                    ).limit(1)
                ).scalar()
                if canonical:
                    connection.execute(
                        ADMETModelRegistry.__table__.update().where(ADMETModelRegistry.id == canonical)
                        .values(**{key: value for key, value in values.items() if key != "endpoint_name"})
                    )


ADMET_CSV_COLUMNS = [
    "compound_id", "version_number", "endpoint", "species", "matrix", "value", "qualitative_value", "unit",
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
    qualitative_value: Mapped[str] = mapped_column(String(120), default="")
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
    source: Mapped[str] = mapped_column(Text, default="")
    training_dataset: Mapped[str] = mapped_column(Text, default="")
    validation_json: Mapped[dict] = mapped_column(JSON, default=dict)
    license: Mapped[str] = mapped_column(Text, default="")
    model_priority: Mapped[int] = mapped_column(Integer, default=100)
    ensemble_eligible: Mapped[bool] = mapped_column(default=True)
    species: Mapped[str] = mapped_column(String(100), default="")
    output_type: Mapped[str] = mapped_column(String(60), default="")
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
    model_version: Mapped[str] = mapped_column(String(60), default="0")
    execution_status: Mapped[str] = mapped_column(String(40), default="SUCCESS")
    standardizer_version: Mapped[str] = mapped_column(String(40), default="CHEM_STANDARDIZER_V1")
    canonical_smiles: Mapped[str] = mapped_column(Text, default="")
    runtime_ms: Mapped[float] = mapped_column(Float, default=0.0)
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


class ADMETConsensusPrediction(Base):
    """One reproducible endpoint consensus assembled from model-specific predictions."""
    __tablename__ = "admet_consensus_predictions"
    id: Mapped[int] = mapped_column(primary_key=True)
    run_id: Mapped[int] = mapped_column(ForeignKey("admet_prediction_runs.id", ondelete="CASCADE"), index=True)
    endpoint_id: Mapped[int] = mapped_column(ForeignKey("admet_endpoints.id", ondelete="CASCADE"), index=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    consensus_version: Mapped[str] = mapped_column(String(60), default="stage4d1-static-v1")
    consensus_mode: Mapped[str] = mapped_column(String(30), default="SHADOW")
    combined_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    unit: Mapped[str] = mapped_column(String(40), default="")
    classification: Mapped[str] = mapped_column(String(120), default="")
    confidence: Mapped[str] = mapped_column(String(30), default="UNKNOWN")
    applicability_domain: Mapped[str] = mapped_column(String(40), default="UNKNOWN")
    model_agreement: Mapped[str] = mapped_column(String(40), default="HIGH_AGREEMENT")
    dispersion_json: Mapped[dict] = mapped_column(JSON, default=dict)
    vote_pattern: Mapped[str] = mapped_column(String(120), default="")
    weights_json: Mapped[list] = mapped_column(JSON, default=list)
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    endpoint = relationship("ADMETEndpoint")


class ADMETModelComparison(Base):
    """Experimental feedback retained per model without overwriting predictions."""
    __tablename__ = "admet_model_comparisons"
    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    endpoint_id: Mapped[int] = mapped_column(ForeignKey("admet_endpoints.id", ondelete="CASCADE"), index=True)
    model_id: Mapped[int] = mapped_column(ForeignKey("admet_model_registry.id", ondelete="RESTRICT"), index=True)
    prediction_id: Mapped[int] = mapped_column(ForeignKey("admet_predictions.id", ondelete="CASCADE"), index=True)
    measurement_id: Mapped[int] = mapped_column(ForeignKey("admet_measurements.id", ondelete="CASCADE"), index=True)
    task_type: Mapped[str] = mapped_column(String(40), default="regression")
    predicted_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    experimental_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    absolute_error: Mapped[float | None] = mapped_column(Float, nullable=True)
    squared_error: Mapped[float | None] = mapped_column(Float, nullable=True)
    fold_error: Mapped[float | None] = mapped_column(Float, nullable=True)
    predicted_class: Mapped[str] = mapped_column(String(120), default="")
    experimental_class: Mapped[str] = mapped_column(String(120), default="")
    correct: Mapped[bool | None] = mapped_column(nullable=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    __table_args__ = (UniqueConstraint("prediction_id", "measurement_id", name="uq_admet_model_prediction_measurement"),)


class ADMETModelPerformance(Base):
    """Global and project-scoped metrics used by conservative consensus weighting."""
    __tablename__ = "admet_model_performance"
    id: Mapped[int] = mapped_column(primary_key=True)
    scope_key: Mapped[str] = mapped_column(String(80), index=True)
    project_id: Mapped[int | None] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), nullable=True, index=True)
    endpoint_name: Mapped[str] = mapped_column(String(120), index=True)
    model_id: Mapped[int] = mapped_column(ForeignKey("admet_model_registry.id", ondelete="RESTRICT"), index=True)
    task_type: Mapped[str] = mapped_column(String(40), default="regression")
    sample_size: Mapped[int] = mapped_column(Integer, default=0)
    metrics_json: Mapped[dict] = mapped_column(JSON, default=dict)
    performance_factor: Mapped[float] = mapped_column(Float, default=1.0)
    updated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    __table_args__ = (UniqueConstraint("scope_key", "model_id", name="uq_admet_model_performance_scope"),)


def measurement_out(row: ADMETMeasurement):
    return {
        "id": row.id, "version_id": row.version_id, "endpoint_id": row.endpoint_id,
        "species": row.species, "matrix": row.matrix, "value": row.value,
        "qualitative_value": row.qualitative_value, "unit": row.unit,
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
    if value is None and mean is None and not str(payload.get("qualitative_value") or "").strip():
        raise HTTPException(status_code=400, detail="A numeric value, mean, or qualitative value is required")
    if sd is not None and sd < 0:
        raise HTTPException(status_code=400, detail="SD cannot be negative")
    if sample_size is not None and sample_size < 1:
        raise HTTPException(status_code=400, detail="Sample size must be at least 1")
    if payload.get("qualifier") not in (None, "", "=", "<", "<=", ">", ">=", "~"):
        raise HTTPException(status_code=400, detail="Qualifier must be one of = < <= > >= ~")
    endpoint_str = str(payload.get("endpoint") or "").strip()
    if endpoint_str.lower() in {"logd", "logd(ph)", "distribution coefficient", "logd (ph)"} or (endpoint_str.lower().startswith("logd") and "7.4" not in endpoint_str and "ph" not in endpoint_str.lower()):
        ph_val = payload.get("ph") or payload.get("assay_ph") or (payload.get("provenance") or {}).get("ph")
        notes_str = str(payload.get("notes") or "").lower()
        if ph_val is None and "ph=" not in notes_str and "ph " not in notes_str:
            raise HTTPException(status_code=400, detail="Assay pH is mandatory for experimental logD measurements (e.g. pH 7.4).")
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
            "species": row.species, "matrix": row.matrix, "value": row.value,
            "qualitative_value": row.qualitative_value, "unit": row.unit,
            "qualifier": row.qualifier, "replicate": row.replicate, "mean": row.mean_value,
            "sd": row.standard_deviation, "n": row.sample_size, "method": row.method,
            "source": row.source, "date": row.experiment_date, "notes": row.notes.replace("\n", " "),
        })
    return PlainTextResponse(stream.getvalue(), media_type="text/csv")


def inputs_hash(version_ids: list[int]) -> str:
    payload = ",".join(str(value) for value in sorted(version_ids)).encode()
    return hashlib.sha256(payload).hexdigest()
