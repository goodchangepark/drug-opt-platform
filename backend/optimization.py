"""Persistent Stage 4A optimization strategy runs (no analog generation)."""

from datetime import datetime

from sqlalchemy import JSON, DateTime, ForeignKey, Integer, String, Text
from sqlalchemy.orm import Mapped, mapped_column

from .database import Base
from .models import utcnow


class OptimizationRun(Base):
    __tablename__ = "optimization_runs"

    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    parent_version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    assay_id: Mapped[int | None] = mapped_column(ForeignKey("assay_definitions.id", ondelete="SET NULL"), nullable=True, index=True)
    objectives_json: Mapped[list] = mapped_column(JSON, default=list)
    custom_objective: Mapped[str] = mapped_column(Text, default="")
    constraints_json: Mapped[dict] = mapped_column(JSON, default=dict)
    endpoint_weights_json: Mapped[dict] = mapped_column(JSON, default=dict)
    manual_overrides_json: Mapped[dict] = mapped_column(JSON, default=dict)
    evidence_json: Mapped[dict] = mapped_column(JSON, default=dict)
    liabilities_json: Mapped[list] = mapped_column(JSON, default=list)
    protected_regions_json: Mapped[list] = mapped_column(JSON, default=list)
    modifiable_regions_json: Mapped[list] = mapped_column(JSON, default=list)
    transformations_json: Mapped[list] = mapped_column(JSON, default=list)
    highlighted_svg: Mapped[str] = mapped_column(Text, default="")
    engine_name: Mapped[str] = mapped_column(String(120), default="Stage 4A deterministic strategy engine")
    engine_version: Mapped[str] = mapped_column(String(60), default="4A.1.0")
    status: Mapped[str] = mapped_column(String(40), default="PENDING")
    message: Mapped[str] = mapped_column(Text, default="")
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    completed_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)


def ensure_optimization_schema(engine):
    """Idempotent Stage 4A project-local migration."""
    from sqlalchemy import inspect

    if "projects" not in inspect(engine).get_table_names():
        return
    Base.metadata.create_all(bind=engine, tables=[OptimizationRun.__table__])
