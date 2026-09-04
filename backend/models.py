from datetime import datetime, timezone
from typing import Optional

from sqlalchemy import Boolean, JSON, DateTime, ForeignKey, Integer, String, Text, UniqueConstraint, inspect, text
from sqlalchemy.orm import Mapped, mapped_column, relationship

from .database import Base


def utcnow():
    return datetime.now(timezone.utc)


class Project(Base):
    __tablename__ = "projects"
    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[str] = mapped_column(String(200), unique=True, index=True)
    target: Mapped[str] = mapped_column(String(300), default="")
    molecule_type: Mapped[str] = mapped_column(String(40), default="Small Molecule")
    indication: Mapped[str] = mapped_column(String(300), default="")
    mechanism_modality: Mapped[str] = mapped_column(String(300), default="")
    description: Mapped[str] = mapped_column(Text, default="")
    is_test_fixture: Mapped[bool] = mapped_column(Boolean, default=False, index=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    updated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow, onupdate=utcnow)

    compounds: Mapped[list["Compound"]] = relationship(back_populates="project", cascade="all, delete-orphan")


class Compound(Base):
    __tablename__ = "compounds"
    id: Mapped[int] = mapped_column(primary_key=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    compound_id: Mapped[str] = mapped_column(String(50), index=True)
    cas_number: Mapped[Optional[str]] = mapped_column(String(12), nullable=True, default=None, index=True)
    name: Mapped[str] = mapped_column(String(200), default="")
    notes: Mapped[str] = mapped_column(Text, default="")
    status: Mapped[str] = mapped_column(String(40), default="CALCULATED", index=True)
    current_version: Mapped[int] = mapped_column(Integer, default=1)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    updated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow, onupdate=utcnow)

    project: Mapped[Project] = relationship(back_populates="compounds")
    versions: Mapped[list["CompoundVersion"]] = relationship(
        back_populates="compound", cascade="all, delete-orphan", order_by="CompoundVersion.version_number"
    )
    __table_args__ = (UniqueConstraint("project_id", "compound_id", name="uq_compound_project_label"),)


class ExternalExperimentalEvidence(Base):
    """Canonical persisted scientific evidence, regardless of external/internal origin.

    The historical table name is retained for migration compatibility.  Its
    lifecycle/origin fields, not its name, define the scientific evidence
    contract used by search, import, and manual capture.
    """
    __tablename__ = "external_experimental_evidence"
    id: Mapped[int] = mapped_column(primary_key=True)
    compound_version_id: Mapped[int] = mapped_column(ForeignKey("compound_versions.id", ondelete="CASCADE"), index=True)
    provenance_key: Mapped[str] = mapped_column(String(64), unique=True, index=True)
    cas_number: Mapped[str] = mapped_column(String(12), default="")
    raw_endpoint_name: Mapped[str] = mapped_column(String(120))
    raw_value: Mapped[str] = mapped_column(Text)
    raw_relation: Mapped[str] = mapped_column(String(12), default="=")
    raw_unit: Mapped[str] = mapped_column(String(80), default="")
    assay_type: Mapped[str] = mapped_column(String(80), default="")
    assay_conditions_json: Mapped[dict] = mapped_column(JSON, default=dict)
    species: Mapped[str] = mapped_column(String(100), default="")
    source_database: Mapped[str] = mapped_column(String(40))
    source_record_id: Mapped[str] = mapped_column(String(160))
    source_assay_id: Mapped[str] = mapped_column(String(160), default="")
    source_document_id: Mapped[str] = mapped_column(String(160), default="")
    reference_text: Mapped[str] = mapped_column(Text, default="")
    source_url: Mapped[str] = mapped_column(Text, default="")
    identity_match_status: Mapped[str] = mapped_column(String(50))
    endpoint_match_status: Mapped[str] = mapped_column(String(50))
    mapping_status: Mapped[str] = mapped_column(String(60), default="EXTERNAL_EVIDENCE_ONLY")
    mapped_assay_id: Mapped[int | None] = mapped_column(ForeignKey("assay_definitions.id", ondelete="SET NULL"), nullable=True, index=True)
    evidence_origin: Mapped[str] = mapped_column(String(40), default="EXPERIMENTAL_EXTERNAL")
    canonical_endpoint_id: Mapped[str] = mapped_column(String(120), default="")
    normalized_value: Mapped[str] = mapped_column(Text, default="")
    normalized_unit: Mapped[str] = mapped_column(String(80), default="")
    normalization_rule: Mapped[str] = mapped_column(String(240), default="")
    normalization_version: Mapped[str] = mapped_column(String(80), default="")
    comparability_status: Mapped[str] = mapped_column(String(60), default="UNSUPPORTED")
    source_quality_class: Mapped[str] = mapped_column(String(4), default="D")
    duplicate_status: Mapped[str] = mapped_column(String(40), default="DISTINCT_MEASUREMENT")
    provenance_fingerprint: Mapped[str] = mapped_column(String(64), default="", index=True)
    # v3.8A: search results are durable candidates before explicit import.
    # ``imported_at`` is retained for compatibility; lifecycle state is the
    # authoritative distinction between searched and accepted evidence.
    evidence_state: Mapped[str] = mapped_column(String(40), default="EXTERNAL_IMPORTED", index=True)
    search_run_id: Mapped[str] = mapped_column(String(80), default="", index=True)
    first_seen_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    last_seen_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    accepted_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    search_version: Mapped[str] = mapped_column(String(80), default="")
    parser_version: Mapped[str] = mapped_column(String(80), default="")
    qualification_version: Mapped[str] = mapped_column(String(80), default="")
    routing_version: Mapped[str] = mapped_column(String(80), default="")
    canonical_endpoint_version: Mapped[str] = mapped_column(String(80), default="")
    unit_normalization_version: Mapped[str] = mapped_column(String(80), default="")
    # v4.1: every unique observation displayed by a completed search is a
    # durable, rehydratable record, including review/related/non-importable
    # observations.  These IDs are display/audit identities, not source data.
    display_evidence_group_id: Mapped[str] = mapped_column(String(120), default="", index=True)
    independent_experiment_group_id: Mapped[str] = mapped_column(String(160), default="", index=True)
    qualification_json: Mapped[dict] = mapped_column(JSON, default=dict)
    qualification_status: Mapped[str] = mapped_column(String(60), default="")
    routing_section: Mapped[str] = mapped_column(String(30), default="")
    routing_reason: Mapped[str] = mapped_column(Text, default="")
    retrieved_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    imported_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    lifecycle_status: Mapped[str] = mapped_column(String(40), default="ACTIVE", index=True)
    revision_number: Mapped[int] = mapped_column(Integer, default=1)
    supersedes_evidence_id: Mapped[int | None] = mapped_column(Integer, nullable=True, index=True)
    updated_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow, onupdate=utcnow)


class EvidenceImportBatch(Base):
    """Auditable explicit acceptance operation for external candidates."""
    __tablename__ = "evidence_import_batches"
    id: Mapped[int] = mapped_column(primary_key=True)
    batch_id: Mapped[str] = mapped_column(String(80), unique=True, index=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    selected_evidence_ids: Mapped[list] = mapped_column(JSON, default=list)
    imported_evidence_ids: Mapped[list] = mapped_column(JSON, default=list)
    skipped_evidence_ids: Mapped[list] = mapped_column(JSON, default=list)
    skip_reasons: Mapped[dict] = mapped_column(JSON, default=dict)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


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


class ExperimentalSearchRun(Base):
    """Immutable audit header for an explicit public evidence search."""
    __tablename__ = "experimental_search_runs"
    id: Mapped[int] = mapped_column(primary_key=True)
    search_run_id: Mapped[str] = mapped_column(String(80), unique=True, index=True)
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id", ondelete="CASCADE"), index=True)
    compound_id: Mapped[int] = mapped_column(ForeignKey("compounds.id", ondelete="CASCADE"), index=True)
    compound_version_id: Mapped[int | None] = mapped_column(ForeignKey("compound_versions.id", ondelete="SET NULL"), nullable=True, index=True)
    query_identity_json: Mapped[dict] = mapped_column(JSON, default=dict)
    identity_graph_version: Mapped[str] = mapped_column(String(80), default="")
    harvester_version: Mapped[str] = mapped_column(String(80), default="")
    parser_version: Mapped[str] = mapped_column(String(80), default="")
    qualification_version: Mapped[str] = mapped_column(String(80), default="")
    routing_version: Mapped[str] = mapped_column(String(80), default="")
    started_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    completed_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    status: Mapped[str] = mapped_column(String(40), default="RUNNING")
    raw_count: Mapped[int] = mapped_column(Integer, default=0)
    unique_count: Mapped[int] = mapped_column(Integer, default=0)
    qualified_count: Mapped[int] = mapped_column(Integer, default=0)
    importable_count: Mapped[int] = mapped_column(Integer, default=0)
    context_qualified_count: Mapped[int] = mapped_column(Integer, default=0)
    persisted_observation_count: Mapped[int] = mapped_column(Integer, default=0)
    display_only_non_persistent_count: Mapped[int] = mapped_column(Integer, default=0)
    source_status_json: Mapped[dict] = mapped_column(JSON, default=dict)
    summary_json: Mapped[dict] = mapped_column(JSON, default=dict)


def ensure_ui_schema(engine):
    """Idempotent, non-destructive migration for the pre-Stage 5 UI workflow."""
    inspector = inspect(engine)
    tables = set(inspector.get_table_names())
    if "projects" not in tables or "compounds" not in tables:
        return
    project_columns = {row["name"] for row in inspector.get_columns("projects")}
    compound_schema = inspector.get_columns("compounds")
    compound_columns = {row["name"] for row in compound_schema}
    cas_column = next((row for row in compound_schema if row["name"] == "cas_number"), None)

    # SQLite cannot alter a column's nullability in place.  The original UI
    # migration created CAS as NOT NULL, while CAS is optional product
    # metadata.  Rebuild only this table, transactionally, preserving every
    # row/ID and the project/version foreign-key contract.  This is an
    # idempotent schema migration, not a database replacement.
    if engine.dialect.name == "sqlite" and cas_column and not cas_column.get("nullable", True):
        raw = engine.raw_connection()
        try:
            raw.execute("PRAGMA foreign_keys=OFF")
            raw.execute("BEGIN")
            raw.execute("""
                CREATE TABLE compounds__cas_nullable (
                    id INTEGER NOT NULL PRIMARY KEY,
                    project_id INTEGER NOT NULL,
                    compound_id VARCHAR(50) NOT NULL,
                    cas_number VARCHAR(12),
                    name VARCHAR(200) NOT NULL,
                    notes TEXT NOT NULL,
                    status VARCHAR(40) NOT NULL,
                    current_version INTEGER NOT NULL,
                    created_at DATETIME NOT NULL,
                    updated_at DATETIME NOT NULL,
                    CONSTRAINT uq_compound_project_label UNIQUE (project_id, compound_id),
                    FOREIGN KEY(project_id) REFERENCES projects (id) ON DELETE CASCADE
                )
            """)
            raw.execute("""
                INSERT INTO compounds__cas_nullable
                (id, project_id, compound_id, cas_number, name, notes, status, current_version, created_at, updated_at)
                SELECT id, project_id, compound_id, NULLIF(trim(cas_number), ''), name, notes, status, current_version, created_at, updated_at
                FROM compounds
            """)
            raw.execute("DROP TABLE compounds")
            raw.execute("ALTER TABLE compounds__cas_nullable RENAME TO compounds")
            raw.execute("CREATE INDEX IF NOT EXISTS ix_compounds_project_id ON compounds (project_id)")
            raw.execute("CREATE INDEX IF NOT EXISTS ix_compounds_status ON compounds (status)")
            raw.execute("CREATE INDEX IF NOT EXISTS ix_compounds_compound_id ON compounds (compound_id)")
            raw.execute("CREATE INDEX IF NOT EXISTS ix_compounds_cas_number ON compounds (cas_number)")
            raw.commit()
        except Exception:
            raw.rollback()
            raise
        finally:
            raw.execute("PRAGMA foreign_keys=ON")
            raw.close()
        inspector = inspect(engine)
        compound_schema = inspector.get_columns("compounds")
        compound_columns = {row["name"] for row in compound_schema}
    with engine.begin() as connection:
        if "is_test_fixture" not in project_columns:
            connection.execute(text("ALTER TABLE projects ADD COLUMN is_test_fixture BOOLEAN NOT NULL DEFAULT 0"))
            connection.execute(text("CREATE INDEX IF NOT EXISTS ix_projects_is_test_fixture ON projects (is_test_fixture)"))
        if "molecule_type" not in project_columns:
            connection.execute(text("ALTER TABLE projects ADD COLUMN molecule_type VARCHAR(40) NOT NULL DEFAULT 'Small Molecule'"))
        if "status" not in compound_columns:
            connection.execute(text("ALTER TABLE compounds ADD COLUMN status VARCHAR(40) NOT NULL DEFAULT 'CALCULATED'"))
        if "cas_number" not in compound_columns:
            connection.execute(text("ALTER TABLE compounds ADD COLUMN cas_number VARCHAR(12) NOT NULL DEFAULT ''"))
        if "external_experimental_evidence" in tables:
            evidence_columns = {row["name"] for row in inspector.get_columns("external_experimental_evidence")}
            if "mapping_status" not in evidence_columns:
                connection.execute(text("ALTER TABLE external_experimental_evidence ADD COLUMN mapping_status VARCHAR(60) NOT NULL DEFAULT 'EXTERNAL_EVIDENCE_ONLY'"))
            if "mapped_assay_id" not in evidence_columns:
                connection.execute(text("ALTER TABLE external_experimental_evidence ADD COLUMN mapped_assay_id INTEGER"))
            for name, definition in {
                "canonical_endpoint_id": "VARCHAR(120) NOT NULL DEFAULT ''",
                "normalized_value": "TEXT NOT NULL DEFAULT ''",
                "normalized_unit": "VARCHAR(80) NOT NULL DEFAULT ''",
                "normalization_rule": "VARCHAR(240) NOT NULL DEFAULT ''",
                "normalization_version": "VARCHAR(80) NOT NULL DEFAULT ''",
                "comparability_status": "VARCHAR(60) NOT NULL DEFAULT 'UNSUPPORTED'",
                "source_quality_class": "VARCHAR(4) NOT NULL DEFAULT 'D'",
                "duplicate_status": "VARCHAR(40) NOT NULL DEFAULT 'DISTINCT_MEASUREMENT'",
                "provenance_fingerprint": "VARCHAR(64) NOT NULL DEFAULT ''",
                "evidence_state": "VARCHAR(40) NOT NULL DEFAULT 'EXTERNAL_IMPORTED'",
                "search_run_id": "VARCHAR(80) NOT NULL DEFAULT ''",
                "first_seen_at": "DATETIME",
                "last_seen_at": "DATETIME",
                "accepted_at": "DATETIME",
                "search_version": "VARCHAR(80) NOT NULL DEFAULT ''",
                "parser_version": "VARCHAR(80) NOT NULL DEFAULT ''",
                "qualification_version": "VARCHAR(80) NOT NULL DEFAULT ''",
                "routing_version": "VARCHAR(80) NOT NULL DEFAULT ''",
                "qualification_status": "VARCHAR(60) NOT NULL DEFAULT ''",
                "routing_section": "VARCHAR(30) NOT NULL DEFAULT ''",
                "routing_reason": "TEXT NOT NULL DEFAULT ''",
                "canonical_endpoint_version": "VARCHAR(80) NOT NULL DEFAULT ''",
                "unit_normalization_version": "VARCHAR(80) NOT NULL DEFAULT ''",
                "display_evidence_group_id": "VARCHAR(120) NOT NULL DEFAULT ''",
                "independent_experiment_group_id": "VARCHAR(160) NOT NULL DEFAULT ''",
                "qualification_json": "JSON NOT NULL DEFAULT '{}'",
                "lifecycle_status": "VARCHAR(40) NOT NULL DEFAULT 'ACTIVE'",
                "revision_number": "INTEGER NOT NULL DEFAULT 1",
                "supersedes_evidence_id": "INTEGER",
                "updated_at": "DATETIME",
            }.items():
                if name not in evidence_columns:
                    connection.execute(text(f"ALTER TABLE external_experimental_evidence ADD COLUMN {name} {definition}"))
            # Existing rows were created only by the explicit import endpoint;
            # preserve that meaning when the lifecycle columns are introduced.
            connection.execute(text("UPDATE external_experimental_evidence SET evidence_state='EXTERNAL_IMPORTED' WHERE evidence_state IS NULL OR trim(evidence_state)=''"))
            connection.execute(text("UPDATE external_experimental_evidence SET lifecycle_status='ACTIVE' WHERE lifecycle_status IS NULL OR trim(lifecycle_status)=''"))
            connection.execute(text("UPDATE external_experimental_evidence SET first_seen_at=COALESCE(first_seen_at, imported_at), last_seen_at=COALESCE(last_seen_at, imported_at)"))
        if "experimental_search_runs" in tables:
            search_columns = {row["name"] for row in inspector.get_columns("experimental_search_runs")}
            for name, definition in {
                "context_qualified_count": "INTEGER NOT NULL DEFAULT 0",
                "persisted_observation_count": "INTEGER NOT NULL DEFAULT 0",
                "display_only_non_persistent_count": "INTEGER NOT NULL DEFAULT 0",
            }.items():
                if name not in search_columns:
                    connection.execute(text(f"ALTER TABLE experimental_search_runs ADD COLUMN {name} {definition}"))
        connection.execute(text("UPDATE projects SET molecule_type='Small Molecule' WHERE molecule_type IS NULL OR trim(molecule_type)=''"))
        connection.execute(text("UPDATE compounds SET status='CALCULATED' WHERE status IS NULL OR trim(status)=''"))


# Register dependent metadata without circular import
import backend.activity_models  # noqa: E402, F401
