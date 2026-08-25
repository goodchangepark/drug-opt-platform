from datetime import datetime, timezone

from sqlalchemy import JSON, DateTime, ForeignKey, Integer, String, Text, UniqueConstraint
from sqlalchemy.orm import Mapped, mapped_column, relationship

from .database import Base


def utcnow():
    return datetime.now(timezone.utc)


class Project(Base):
    __tablename__ = "projects"
    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[str] = mapped_column(String(200), unique=True, index=True)
    target: Mapped[str] = mapped_column(String(300), default="")
    indication: Mapped[str] = mapped_column(String(300), default="")
    mechanism_modality: Mapped[str] = mapped_column(String(300), default="")
    description: Mapped[str] = mapped_column(Text, default="")
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    updated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow, onupdate=utcnow)

    compounds: Mapped[list["Compound"]] = relationship(back_populates="project", cascade="all, delete-orphan")


class Compound(Base):
    __tablename__ = "compounds"
    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    compound_id: Mapped[str] = mapped_column(String(50), index=True)
    name: Mapped[str] = mapped_column(String(200), default="")
    notes: Mapped[str] = mapped_column(Text, default="")
    current_version: Mapped[int] = mapped_column(Integer, default=1)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    updated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow, onupdate=utcnow)

    project: Mapped[Project] = relationship(back_populates="compounds")
    versions: Mapped[list["CompoundVersion"]] = relationship(
        back_populates="compound", cascade="all, delete-orphan", order_by="CompoundVersion.version_number"
    )
    __table_args__ = (UniqueConstraint("project_id", "compound_id", name="uq_compound_project_label"),)


class CompoundVersion(Base):
    __tablename__ = "compound_versions"
    id: Mapped[int] = mapped_column(primary_key=True)
    compound_row_id: Mapped[int] = mapped_column(ForeignKey("compounds.id", ondelete="CASCADE"), index=True)
    version_number: Mapped[int]
    original_smiles: Mapped[str] = mapped_column(Text)
    canonical_smiles: Mapped[str] = mapped_column(String(1000), index=True)
    isomeric_smiles: Mapped[str] = mapped_column(String(1000))
    inchi: Mapped[str] = mapped_column(Text, default="")
    inchikey: Mapped[str] = mapped_column(String(32), index=True)
    change_note: Mapped[str] = mapped_column(String(500), default="")
    properties_json: Mapped[dict | None] = mapped_column(JSON, nullable=True)
    alerts_json: Mapped[list | None] = mapped_column(JSON, nullable=True)
    assessment_json: Mapped[dict | None] = mapped_column(JSON, nullable=True)
    calculation_json: Mapped[dict | None] = mapped_column(JSON, nullable=True)
    svg: Mapped[str] = mapped_column(Text, default="")
    highlighted_svg: Mapped[str] = mapped_column(Text, default="")
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)

    compound: Mapped[Compound] = relationship(back_populates="versions")
    property_runs: Mapped[list["PropertyCalculation"]] = relationship(back_populates="version", cascade="all, delete-orphan")
    structural_alerts: Mapped[list["StructuralAlert"]] = relationship(back_populates="version", cascade="all, delete-orphan")
    prediction_runs: Mapped[list["PredictionRun"]] = relationship(back_populates="version", cascade="all, delete-orphan")
    __table_args__ = (UniqueConstraint("compound_row_id", "version_number", name="uq_compound_version"),)


class PropertyCalculation(Base):
    __tablename__ = "property_calculations"
    id: Mapped[int] = mapped_column(primary_key=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    endpoint: Mapped[str] = mapped_column(String(80), index=True)
    value_type: Mapped[str] = mapped_column(String(20), default="Calculated")
    value: Mapped[str] = mapped_column(Text)
    engine: Mapped[str] = mapped_column(String(60))
    method: Mapped[str] = mapped_column(String(160))
    engine_version: Mapped[str] = mapped_column(String(40))
    calculated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)

    version: Mapped[CompoundVersion] = relationship(back_populates="property_runs")
    __table_args__ = (UniqueConstraint("version_id", "endpoint", "method", name="uq_property_provenance"),)


class StructuralAlert(Base):
    __tablename__ = "structural_alerts"
    id: Mapped[int] = mapped_column(primary_key=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    alert_set: Mapped[str] = mapped_column(String(80))
    alert_name: Mapped[str] = mapped_column(String(200))
    reason: Mapped[str] = mapped_column(Text, default="")
    matched_smiles: Mapped[str] = mapped_column(String(500), default="")
    matched_atoms_json: Mapped[list] = mapped_column(JSON, default=list)

    version: Mapped[CompoundVersion] = relationship(back_populates="structural_alerts")


class PredictionRun(Base):
    """Immutable audit record. Stage 2-5 predictions can reuse this table."""

    __tablename__ = "prediction_runs"
    id: Mapped[int] = mapped_column(primary_key=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    stage: Mapped[str] = mapped_column(String(30), default="stage_1")
    model_name: Mapped[str] = mapped_column(String(120))
    model_version: Mapped[str] = mapped_column(String(60))
    inputs_hash: Mapped[str] = mapped_column(String(64))
    outputs_json: Mapped[dict] = mapped_column(JSON, default=dict)
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)
    confidence: Mapped[str] = mapped_column(String(30), default="High")
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)

    version: Mapped[CompoundVersion] = relationship(back_populates="prediction_runs")
