from datetime import datetime, timezone

from sqlalchemy import JSON, DateTime, Float, ForeignKey, Integer, String, Text
from sqlalchemy.orm import Mapped, mapped_column, relationship

from .database import Base
from .models import utcnow


class AssayDefinition(Base):
    __tablename__ = "assay_definitions"
    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    assay_uid: Mapped[str] = mapped_column(String(64), unique=True, default=lambda: f"AS-{datetime.now(timezone.utc).timestamp():.0f}")
    version_number: Mapped[int] = mapped_column(Integer, default=1)
    supersedes_id: Mapped[int | None] = mapped_column(ForeignKey("assay_definitions.id", ondelete="SET NULL"))
    active: Mapped[bool] = mapped_column(default=True)
    name: Mapped[str] = mapped_column(String(200), index=True)
    target: Mapped[str] = mapped_column(String(300), default="")
    target_type: Mapped[str] = mapped_column(String(120), default="")
    assay_category: Mapped[str] = mapped_column(String(120), default="")
    measurement_type: Mapped[str] = mapped_column(String(80), default="IC50")
    custom_measurement_name: Mapped[str] = mapped_column(String(120), default="")
    unit: Mapped[str] = mapped_column(String(20), default="nM")
    species: Mapped[str] = mapped_column(String(100), default="")
    cell_line: Mapped[str] = mapped_column(String(200), default="")
    mutation_variant: Mapped[str] = mapped_column(String(300), default="")
    protein_construct: Mapped[str] = mapped_column(String(300), default="")
    substrate: Mapped[str] = mapped_column(String(300), default="")
    atp_concentration: Mapped[str] = mapped_column(String(100), default="")
    incubation_time: Mapped[str] = mapped_column(String(100), default="")
    detection_method: Mapped[str] = mapped_column(String(200), default="")
    experimental_conditions: Mapped[str] = mapped_column(Text, default="")
    protocol: Mapped[str] = mapped_column(Text, default="")
    reference_compound: Mapped[str] = mapped_column(String(300), default="")
    reference_structure_smiles: Mapped[str] = mapped_column(Text, default="")
    reference_activity: Mapped[float | None] = mapped_column(Float, nullable=True)
    reference_source: Mapped[str] = mapped_column(Text, default="")
    reference_provenance_url: Mapped[str] = mapped_column(Text, default="")
    notes: Mapped[str] = mapped_column(Text, default="")
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)

    project_relationship = relationship("Project")
    measurements: Mapped[list["ActivityMeasurement"]] = relationship(back_populates="assay", cascade="all, delete-orphan")
    predictions: Mapped[list["ActivityPrediction"]] = relationship(back_populates="assay", cascade="all, delete-orphan")
    models: Mapped[list["QSARModel"]] = relationship(back_populates="assay", cascade="all, delete-orphan")


class ActivityMeasurement(Base):
    __tablename__ = "activity_measurements"
    id: Mapped[int] = mapped_column(primary_key=True)
    assay_id: Mapped[int] = mapped_column(ForeignKey("assay_definitions.id", ondelete="CASCADE"), index=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    raw_value: Mapped[float] = mapped_column(Float)
    original_unit: Mapped[str] = mapped_column(String(20))
    normalized_value_nm: Mapped[float] = mapped_column(Float, index=True)
    qualifier: Mapped[str] = mapped_column(String(3), default="=")
    replicate_label: Mapped[str] = mapped_column(String(50), default="R1")
    experiment_date: Mapped[str] = mapped_column(String(30), default="")
    source: Mapped[str] = mapped_column(String(120), default="User experimental")
    notes: Mapped[str] = mapped_column(Text, default="")
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)

    assay: Mapped[AssayDefinition] = relationship(back_populates="measurements")


class QSARModel(Base):
    __tablename__ = "qsar_models"
    id: Mapped[int] = mapped_column(primary_key=True)
    assay_id: Mapped[int] = mapped_column(ForeignKey("assay_definitions.id", ondelete="CASCADE"), index=True)
    model_uid: Mapped[str] = mapped_column(String(64), unique=True, default=lambda: f"QSR-{datetime.now(timezone.utc).timestamp():.0f}")
    algorithm: Mapped[str] = mapped_column(String(100))
    sklearn_version: Mapped[str] = mapped_column(String(40))
    rdkit_version: Mapped[str] = mapped_column(String(40))
    random_seed: Mapped[int] = mapped_column(Integer, default=42)
    fingerprint_config: Mapped[dict] = mapped_column(JSON, default=dict)
    descriptor_config: Mapped[list] = mapped_column(JSON, default=list)
    training_n: Mapped[int]
    validation_method: Mapped[str] = mapped_column(String(80), default="5-fold CV + scaffold split")
    metrics: Mapped[dict] = mapped_column(JSON, default=dict)
    selection_reason: Mapped[str] = mapped_column(Text, default="")
    pickle_data: Mapped[str] = mapped_column(Text)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)

    assay: Mapped[AssayDefinition] = relationship(back_populates="models")


class ActivityPrediction(Base):
    __tablename__ = "activity_predictions"
    id: Mapped[int] = mapped_column(primary_key=True)
    assay_id: Mapped[int] = mapped_column(ForeignKey("assay_definitions.id", ondelete="CASCADE"), index=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    model_id: Mapped[int | None] = mapped_column(ForeignKey("qsar_models.id", ondelete="SET NULL"))
    prediction_type: Mapped[str] = mapped_column(String(60), default="similarity")
    predicted_pactivity: Mapped[float] = mapped_column(Float)
    predicted_value_nm: Mapped[float] = mapped_column(Float)
    confidence: Mapped[str] = mapped_column(String(30))
    applicability_domain: Mapped[str] = mapped_column(String(30))
    nearest_neighbors: Mapped[list] = mapped_column(JSON, default=list)
    uncertainty: Mapped[float | None] = mapped_column(Float, nullable=True)
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)

    assay: Mapped[AssayDefinition] = relationship(back_populates="predictions")
