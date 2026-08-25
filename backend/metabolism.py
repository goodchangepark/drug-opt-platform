"""Persistent Stage 3D metabolic soft spots and metabolite hypotheses."""

from __future__ import annotations

from datetime import datetime

from sqlalchemy import JSON, DateTime, Float, ForeignKey, Integer, String, Text, UniqueConstraint
from sqlalchemy.orm import Mapped, mapped_column, relationship

from .database import Base
from .models import utcnow


class MetabolicPredictionRun(Base):
    __tablename__ = "metabolic_prediction_runs"
    id: Mapped[int] = mapped_column(primary_key=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    inputs_hash: Mapped[str] = mapped_column(String(64), index=True)
    engine_name: Mapped[str] = mapped_column(String(120))
    engine_version: Mapped[str] = mapped_column(String(100), index=True)
    status: Mapped[str] = mapped_column(String(40), default="RUNNING")
    message: Mapped[str] = mapped_column(Text, default="")
    model_status_json: Mapped[dict] = mapped_column(JSON, default=dict)
    liability_summary_json: Mapped[dict] = mapped_column(JSON, default=dict)
    highlighted_svg: Mapped[str] = mapped_column(Text, default="")
    started_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    completed_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    spots = relationship("MetabolicSoftSpot", back_populates="run", cascade="all, delete-orphan")
    metabolites = relationship("PredictedMetabolite", back_populates="run", cascade="all, delete-orphan")


class MetabolicSoftSpot(Base):
    __tablename__ = "metabolic_soft_spots"
    id: Mapped[int] = mapped_column(primary_key=True)
    run_id: Mapped[int] = mapped_column(ForeignKey("metabolic_prediction_runs.id", ondelete="CASCADE"), index=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    rank: Mapped[int] = mapped_column(Integer)
    atom_index: Mapped[int] = mapped_column(Integer)
    atom_environment: Mapped[str] = mapped_column(String(500), default="")
    transformation: Mapped[str] = mapped_column(String(100), index=True)
    phase: Mapped[str] = mapped_column(String(20), index=True)
    cyp_isoform: Mapped[str] = mapped_column(String(80), default="CYP isoform not assigned")
    model_evidence_json: Mapped[dict] = mapped_column(JSON, default=dict)
    rule_evidence_json: Mapped[dict] = mapped_column(JSON, default=dict)
    score: Mapped[float | None] = mapped_column(Float, nullable=True)
    score_type: Mapped[str] = mapped_column(String(180), default="")
    confidence: Mapped[str] = mapped_column(String(30), default="LOW")
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    run = relationship("MetabolicPredictionRun", back_populates="spots")
    metabolites = relationship("PredictedMetabolite", back_populates="soft_spot")
    __table_args__ = (UniqueConstraint("run_id", "rank", name="uq_metabolic_spot_run_rank"),)


class PredictedMetabolite(Base):
    __tablename__ = "predicted_metabolites"
    id: Mapped[int] = mapped_column(primary_key=True)
    run_id: Mapped[int] = mapped_column(ForeignKey("metabolic_prediction_runs.id", ondelete="CASCADE"), index=True)
    soft_spot_id: Mapped[int] = mapped_column(ForeignKey("metabolic_soft_spots.id", ondelete="CASCADE"), index=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    canonical_smiles: Mapped[str] = mapped_column(String(2000), index=True)
    isomeric_smiles: Mapped[str] = mapped_column(String(2000))
    transformation: Mapped[str] = mapped_column(String(100), index=True)
    source_atom: Mapped[int] = mapped_column(Integer)
    phase: Mapped[str] = mapped_column(String(20))
    rank: Mapped[int] = mapped_column(Integer)
    confidence: Mapped[str] = mapped_column(String(30), default="LOW")
    evidence_json: Mapped[dict] = mapped_column(JSON, default=dict)
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    run = relationship("MetabolicPredictionRun", back_populates="metabolites")
    soft_spot = relationship("MetabolicSoftSpot", back_populates="metabolites")
    __table_args__ = (UniqueConstraint("run_id", "canonical_smiles", name="uq_predicted_metabolite_run_structure"),)


class ExperimentalMetabolite(Base):
    __tablename__ = "experimental_metabolites"
    id: Mapped[int] = mapped_column(primary_key=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    canonical_smiles: Mapped[str] = mapped_column(String(2000), default="")
    isomeric_smiles: Mapped[str] = mapped_column(String(2000), default="")
    transformation: Mapped[str] = mapped_column(String(120))
    observed_mass: Mapped[float | None] = mapped_column(Float, nullable=True)
    mass_unit: Mapped[str] = mapped_column(String(30), default="")
    source: Mapped[str] = mapped_column(String(300), default="")
    experiment: Mapped[str] = mapped_column(String(300), default="")
    notes: Mapped[str] = mapped_column(Text, default="")
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


def ensure_metabolism_schema(engine):
    """Idempotent project-local migration for Stage 3D tables."""
    from sqlalchemy import inspect

    if "projects" not in inspect(engine).get_table_names():
        return
    Base.metadata.create_all(bind=engine, tables=[
        MetabolicPredictionRun.__table__, MetabolicSoftSpot.__table__,
        PredictedMetabolite.__table__, ExperimentalMetabolite.__table__,
    ])
