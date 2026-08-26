"""Persistent Stage 4B proposal, candidate, prediction, and ranking records."""

from datetime import datetime

from sqlalchemy import JSON, Boolean, DateTime, Float, ForeignKey, Integer, String, Text
from sqlalchemy.orm import Mapped, mapped_column, relationship

from .database import Base
from .models import utcnow


class OptimizationProposalRun(Base):
    __tablename__ = "optimization_proposal_runs"

    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    optimization_run_id: Mapped[int] = mapped_column(ForeignKey("optimization_runs.id", ondelete="CASCADE"), index=True)
    parent_version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    status: Mapped[str] = mapped_column(String(40), default="PENDING", index=True)
    stage_message: Mapped[str] = mapped_column(Text, default="Queued")
    transformation_library_version: Mapped[str] = mapped_column(String(60), default="4B.1.0")
    model_versions_json: Mapped[dict] = mapped_column(JSON, default=dict)
    endpoint_weights_json: Mapped[dict] = mapped_column(JSON, default=dict)
    hard_constraints_json: Mapped[dict] = mapped_column(JSON, default=dict)
    settings_json: Mapped[dict] = mapped_column(JSON, default=dict)
    random_seed: Mapped[int] = mapped_column(Integer, default=42)
    raw_candidate_count: Mapped[int] = mapped_column(Integer, default=0)
    accepted_count: Mapped[int] = mapped_column(Integer, default=0)
    rejected_count: Mapped[int] = mapped_column(Integer, default=0)
    top_count: Mapped[int] = mapped_column(Integer, default=0)
    summary_json: Mapped[dict] = mapped_column(JSON, default=dict)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    started_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    completed_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)

    candidates: Mapped[list["OptimizationCandidate"]] = relationship(back_populates="proposal_run", cascade="all, delete-orphan")


class OptimizationCandidate(Base):
    __tablename__ = "optimization_candidates"

    id: Mapped[int] = mapped_column(primary_key=True)
    proposal_run_id: Mapped[int] = mapped_column(ForeignKey("optimization_proposal_runs.id", ondelete="CASCADE"), index=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    optimization_run_id: Mapped[int] = mapped_column(ForeignKey("optimization_runs.id", ondelete="CASCADE"), index=True)
    parent_version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    existing_version_id: Mapped[int | None] = mapped_column(ForeignKey("compound_versions.id", ondelete="SET NULL"), nullable=True, index=True)
    candidate_number: Mapped[int] = mapped_column(Integer)
    canonical_smiles: Mapped[str] = mapped_column(String(1200))
    isomeric_smiles: Mapped[str] = mapped_column(String(1200))
    inchikey: Mapped[str] = mapped_column(String(32), index=True)
    generation_priority: Mapped[int] = mapped_column(Integer, default=5)
    generation_source: Mapped[str] = mapped_column(String(100), default="Curated medicinal chemistry transformation")
    hypothesis: Mapped[str] = mapped_column(Text, default="")
    why_generated: Mapped[str] = mapped_column(Text, default="")
    expected_benefit: Mapped[str] = mapped_column(Text, default="")
    generation_timestamp: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    status: Mapped[str] = mapped_column(String(40), default="GENERATED", index=True)
    rejection_stage: Mapped[str] = mapped_column(String(40), default="")
    stage1_json: Mapped[dict] = mapped_column(JSON, default=dict)
    property_delta_json: Mapped[dict] = mapped_column(JSON, default=dict)
    activity_json: Mapped[dict] = mapped_column(JSON, default=dict)
    admet_json: Mapped[dict] = mapped_column(JSON, default=dict)
    soft_spot_json: Mapped[dict] = mapped_column(JSON, default=dict)
    soft_spot_change_json: Mapped[dict] = mapped_column(JSON, default=dict)
    synthetic_feasibility_json: Mapped[dict] = mapped_column(JSON, default=dict)
    parent_similarity: Mapped[float] = mapped_column(Float, default=0.0)
    mcs_coverage: Mapped[float] = mapped_column(Float, default=0.0)
    changed_parent_atoms_json: Mapped[list] = mapped_column(JSON, default=list)
    changed_candidate_atoms_json: Mapped[list] = mapped_column(JSON, default=list)
    structure_svg: Mapped[str] = mapped_column(Text, default="")
    parent_difference_svg: Mapped[str] = mapped_column(Text, default="")
    candidate_difference_svg: Mapped[str] = mapped_column(Text, default="")
    confidence: Mapped[str] = mapped_column(String(30), default="UNKNOWN")
    applicability_domain: Mapped[str] = mapped_column(String(40), default="UNKNOWN")
    objective_vector_json: Mapped[dict] = mapped_column(JSON, default=dict)
    ranking_score: Mapped[float | None] = mapped_column(Float, nullable=True)
    pareto_front: Mapped[int | None] = mapped_column(Integer, nullable=True)
    information_value: Mapped[str] = mapped_column(String(20), default="UNKNOWN")
    main_risk: Mapped[str] = mapped_column(Text, default="")
    selected_top10: Mapped[bool] = mapped_column(Boolean, default=False, index=True)
    user_added: Mapped[bool] = mapped_column(Boolean, default=False)
    user_decision: Mapped[str] = mapped_column(String(30), default="")
    user_decision_reason: Mapped[str] = mapped_column(Text, default="")

    proposal_run: Mapped[OptimizationProposalRun] = relationship(back_populates="candidates")
    transformations: Mapped[list["CandidateTransformation"]] = relationship(back_populates="candidate", cascade="all, delete-orphan")
    predictions: Mapped[list["CandidatePredictionSnapshot"]] = relationship(back_populates="candidate", cascade="all, delete-orphan")
    rankings: Mapped[list["CandidateRanking"]] = relationship(back_populates="candidate", cascade="all, delete-orphan")
    rejection_reasons: Mapped[list["CandidateRejectionReason"]] = relationship(back_populates="candidate", cascade="all, delete-orphan")


class CandidateTransformation(Base):
    __tablename__ = "candidate_transformations"

    id: Mapped[int] = mapped_column(primary_key=True)
    candidate_id: Mapped[int] = mapped_column(ForeignKey("optimization_candidates.id", ondelete="CASCADE"), index=True)
    sequence_number: Mapped[int] = mapped_column(Integer, default=1)
    transformation_id: Mapped[str] = mapped_column(String(100), index=True)
    name: Mapped[str] = mapped_column(String(240))
    reaction_smarts: Mapped[str] = mapped_column(Text, default="")
    transformation_version: Mapped[str] = mapped_column(String(60), default="")
    source: Mapped[str] = mapped_column(Text, default="")
    source_atom_indices_json: Mapped[list] = mapped_column(JSON, default=list)
    changed_parent_atoms_json: Mapped[list] = mapped_column(JSON, default=list)
    execution_status: Mapped[str] = mapped_column(String(40), default="EXECUTED")
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)

    candidate: Mapped[OptimizationCandidate] = relationship(back_populates="transformations")


class CandidatePredictionSnapshot(Base):
    __tablename__ = "candidate_prediction_snapshots"

    id: Mapped[int] = mapped_column(primary_key=True)
    candidate_id: Mapped[int] = mapped_column(ForeignKey("optimization_candidates.id", ondelete="CASCADE"), index=True)
    stage: Mapped[str] = mapped_column(String(40), index=True)
    endpoint: Mapped[str] = mapped_column(String(160), index=True)
    record_type: Mapped[str] = mapped_column(String(40), default="Predicted")
    value_json: Mapped[dict] = mapped_column(JSON, default=dict)
    unit: Mapped[str] = mapped_column(String(60), default="")
    model_name: Mapped[str] = mapped_column(String(180), default="")
    model_version: Mapped[str] = mapped_column(String(100), default="")
    confidence: Mapped[str] = mapped_column(String(30), default="UNKNOWN")
    applicability_domain: Mapped[str] = mapped_column(String(40), default="UNKNOWN")
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)

    candidate: Mapped[OptimizationCandidate] = relationship(back_populates="predictions")


class CandidateRanking(Base):
    __tablename__ = "candidate_rankings"

    id: Mapped[int] = mapped_column(primary_key=True)
    candidate_id: Mapped[int] = mapped_column(ForeignKey("optimization_candidates.id", ondelete="CASCADE"), index=True)
    rank: Mapped[int]
    score: Mapped[float]
    pareto_front: Mapped[int]
    score_formula_version: Mapped[str] = mapped_column(String(60), default="4B.1.0")
    score_breakdown_json: Mapped[dict] = mapped_column(JSON, default=dict)
    diversity_json: Mapped[dict] = mapped_column(JSON, default=dict)
    selected_top10: Mapped[bool] = mapped_column(Boolean, default=False)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)

    candidate: Mapped[OptimizationCandidate] = relationship(back_populates="rankings")


class CandidateRejectionReason(Base):
    __tablename__ = "candidate_rejection_reasons"

    id: Mapped[int] = mapped_column(primary_key=True)
    candidate_id: Mapped[int] = mapped_column(ForeignKey("optimization_candidates.id", ondelete="CASCADE"), index=True)
    code: Mapped[str] = mapped_column(String(80), index=True)
    detail: Mapped[str] = mapped_column(Text)
    stage: Mapped[str] = mapped_column(String(40))
    hard_constraint: Mapped[bool] = mapped_column(Boolean, default=True)
    evidence_type: Mapped[str] = mapped_column(String(60), default="Calculated")
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)

    candidate: Mapped[OptimizationCandidate] = relationship(back_populates="rejection_reasons")


def ensure_proposal_schema(engine):
    """Idempotent project-local Stage 4B migration."""
    from sqlalchemy import inspect

    if "optimization_runs" not in inspect(engine).get_table_names():
        return
    Base.metadata.create_all(bind=engine, tables=[
        OptimizationProposalRun.__table__, OptimizationCandidate.__table__,
        CandidateTransformation.__table__, CandidatePredictionSnapshot.__table__,
        CandidateRanking.__table__, CandidateRejectionReason.__table__,
    ])
