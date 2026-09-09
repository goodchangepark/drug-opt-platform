"""Drug-OPT Stable Core v1 persistence and provenance contracts.

Legacy calculation/evidence tables remain intact.  This module adds the
authoritative, append-oriented scientific stores and deterministic migration
adapters that sit above them.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from typing import Any

from sqlalchemy import Boolean, DateTime, Float, Integer, JSON, String, Text, UniqueConstraint, inspect, select, text
from sqlalchemy.orm import Mapped, mapped_column

from .canonical_endpoints import canonicalize_prediction_endpoint
from .database import Base
from .species_registry import normalize_species_code

STABLE_CORE_VERSION = "stable-core-v1"
SCHEMA_VERSION = 5
MIGRATION_ID = "stable-core-v1-005"
UNKNOWN_PROVENANCE = "UNKNOWN_PROVENANCE"
UNKNOWN_PREDICTION_MODE = "UNKNOWN_MODE"
VALID_PREDICTION_MODES = frozenset({"ASSISTED", "HYBRID", "FULL_PREDICTION"})


def utcnow() -> datetime:
    return datetime.now(timezone.utc)


def canonical_context(context: dict[str, Any] | None) -> dict[str, Any]:
    source = context if isinstance(context, dict) else {}
    fields = (
        "matrix", "route", "dose", "dose_unit", "formulation", "fed_state",
        "regimen", "dosing_state", "population", "sex", "analyte", "assay_type",
        "measurement_type", "direction", "ph",
    )
    return {field: source[field] for field in fields if source.get(field) not in (None, "")}


def context_identity(context: dict[str, Any] | None) -> str:
    payload = json.dumps(canonical_context(context), sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(payload.encode()).hexdigest()


class ScientificMutationAudit(Base):
    __tablename__ = "scientific_mutation_audit"
    id: Mapped[int] = mapped_column(primary_key=True)
    occurred_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow, index=True)
    actor: Mapped[str] = mapped_column(String(200), default="UNKNOWN")
    process: Mapped[str] = mapped_column(String(200), default="UNKNOWN")
    operation: Mapped[str] = mapped_column(String(80), index=True)
    object_type: Mapped[str] = mapped_column(String(80), index=True)
    object_id: Mapped[str] = mapped_column(String(120), default="", index=True)
    project_id_snapshot: Mapped[int | None] = mapped_column(Integer, nullable=True, index=True)
    transaction_id: Mapped[str] = mapped_column(String(120), default="", index=True)
    before_identity_json: Mapped[dict] = mapped_column(JSON, default=dict)
    after_identity_json: Mapped[dict] = mapped_column(JSON, default=dict)
    reason: Mapped[str] = mapped_column(Text, default="")
    outcome: Mapped[str] = mapped_column(String(40), default="COMMITTED", index=True)


class HistoricalPrediction(Base):
    """Project-independent immutable prediction provenance."""

    __tablename__ = "historical_predictions"
    id: Mapped[int] = mapped_column(primary_key=True)
    legacy_prediction_run_id: Mapped[int] = mapped_column(Integer, unique=True, index=True)
    compound_version_id_snapshot: Mapped[int] = mapped_column(Integer, index=True)
    compound_row_id_snapshot: Mapped[int | None] = mapped_column(Integer, nullable=True, index=True)
    project_id_snapshot: Mapped[int | None] = mapped_column(Integer, nullable=True, index=True)
    project_name_snapshot: Mapped[str] = mapped_column(String(200), default="")
    compound_label_snapshot: Mapped[str] = mapped_column(String(200), default="")
    structure_revision_snapshot: Mapped[str] = mapped_column(String(120), default="")
    canonical_endpoint: Mapped[str] = mapped_column(String(120), default="MULTI_ENDPOINT", index=True)
    species: Mapped[str] = mapped_column(String(40), default="UNSPECIFIED", index=True)
    context_json: Mapped[dict] = mapped_column(JSON, default=dict)
    context_identity: Mapped[str] = mapped_column(String(64), default="", index=True)
    stage: Mapped[str] = mapped_column(String(40), default="")
    model_id: Mapped[str] = mapped_column(String(160), default=UNKNOWN_PROVENANCE)
    model_version: Mapped[str] = mapped_column(String(100), default=UNKNOWN_PROVENANCE)
    model_artifact_hash: Mapped[str] = mapped_column(String(128), default=UNKNOWN_PROVENANCE)
    engine_version: Mapped[str] = mapped_column(String(160), default=UNKNOWN_PROVENANCE, index=True)
    inputs_hash: Mapped[str] = mapped_column(String(64), default="")
    outputs_json: Mapped[dict] = mapped_column(JSON, default=dict)
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)
    applicability_domain_json: Mapped[dict] = mapped_column(JSON, default=dict)
    uncertainty_json: Mapped[dict] = mapped_column(JSON, default=dict)
    prediction_mode: Mapped[str] = mapped_column(String(40), default=UNKNOWN_PREDICTION_MODE)
    predicted_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), index=True)
    immutable: Mapped[bool] = mapped_column(Boolean, default=True)


class ExperimentalObservation(Base):
    """Canonical append-oriented experimental observation."""

    __tablename__ = "experimental_observations"
    id: Mapped[int] = mapped_column(primary_key=True)
    compound_version_id: Mapped[int] = mapped_column(Integer, index=True)
    canonical_endpoint: Mapped[str] = mapped_column(String(120), index=True)
    species: Mapped[str] = mapped_column(String(40), index=True)
    context_json: Mapped[dict] = mapped_column(JSON, default=dict)
    context_identity: Mapped[str] = mapped_column(String(64), index=True)
    value_text: Mapped[str] = mapped_column(Text, default="")
    numeric_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    unit: Mapped[str] = mapped_column(String(80), default="")
    qualifier: Mapped[str] = mapped_column(String(12), default="=")
    lower_bound: Mapped[float | None] = mapped_column(Float, nullable=True)
    upper_bound: Mapped[float | None] = mapped_column(Float, nullable=True)
    source: Mapped[str] = mapped_column(String(120), default="")
    source_url: Mapped[str] = mapped_column(Text, default="")
    source_record_id: Mapped[str] = mapped_column(String(200), default="")
    reference_id: Mapped[str] = mapped_column(String(200), default="")
    identity_confidence: Mapped[str] = mapped_column(String(60), default="UNKNOWN")
    curation_status: Mapped[str] = mapped_column(String(30), default="CANDIDATE", index=True)
    display_comparable: Mapped[bool] = mapped_column(Boolean, default=False)
    numeric_pairable: Mapped[bool] = mapped_column(Boolean, default=False)
    learning_eligible: Mapped[bool] = mapped_column(Boolean, default=False)
    retrieved_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    accepted_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    provenance_json: Mapped[dict] = mapped_column(JSON, default=dict)
    legacy_record_type: Mapped[str] = mapped_column(String(80), index=True)
    legacy_record_id: Mapped[int] = mapped_column(Integer, index=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    __table_args__ = (
        UniqueConstraint("legacy_record_type", "legacy_record_id", name="uq_experimental_observation_legacy"),
    )


class CurrentPredictionSnapshot(Base):
    """Only authoritative source for a current prediction row."""

    __tablename__ = "current_prediction_snapshots"
    id: Mapped[int] = mapped_column(primary_key=True)
    compound_version_id: Mapped[int] = mapped_column(Integer, index=True)
    canonical_endpoint: Mapped[str] = mapped_column(String(120), index=True)
    species: Mapped[str] = mapped_column(String(40), index=True)
    context_json: Mapped[dict] = mapped_column(JSON, default=dict)
    context_identity: Mapped[str] = mapped_column(String(64), index=True)
    engine_release: Mapped[str] = mapped_column(String(160), index=True)
    value: Mapped[float | None] = mapped_column(Float, nullable=True)
    unit: Mapped[str] = mapped_column(String(80), default="")
    classification: Mapped[str] = mapped_column(String(120), default="")
    model_id: Mapped[str] = mapped_column(String(160), default=UNKNOWN_PROVENANCE)
    model_version: Mapped[str] = mapped_column(String(100), default=UNKNOWN_PROVENANCE)
    model_artifact_hash: Mapped[str] = mapped_column(String(128), default=UNKNOWN_PROVENANCE)
    applicability_domain_json: Mapped[dict] = mapped_column(JSON, default=dict)
    uncertainty_json: Mapped[dict] = mapped_column(JSON, default=dict)
    prediction_mode: Mapped[str] = mapped_column(String(40), default=UNKNOWN_PREDICTION_MODE)
    source_artifact_type: Mapped[str] = mapped_column(String(80), default="")
    source_artifact_id: Mapped[int | None] = mapped_column(Integer, nullable=True)
    structure_revision: Mapped[str] = mapped_column(String(120), default="")
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    replaced_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    is_current: Mapped[bool] = mapped_column(Boolean, default=True, index=True)
    __table_args__ = (
        UniqueConstraint(
            "compound_version_id", "canonical_endpoint", "species", "context_identity", "engine_release",
            name="uq_current_prediction_scientific_key",
        ),
    )


def snapshot_admission_reason(snapshot: Any) -> str | None:
    """Return a fail-closed reason when a snapshot cannot be current.

    This is deliberately provenance-strict: legacy rows are retained for
    audit, but a current engine row is selectable only when its scientific
    identity and generating artifact are explicit and internally consistent.
    """
    if not getattr(snapshot, "is_current", False):
        return "NOT_CURRENT"
    # Resolve the release from the single authoritative registry rather than
    # duplicating a literal in the admission path.
    from .prediction_engine_registry import CURRENT_ENGINE_ID
    if getattr(snapshot, "engine_release", None) != CURRENT_ENGINE_ID:
        return "WRONG_ENGINE"
    endpoint_id = str(getattr(snapshot, "canonical_endpoint", "") or "").strip().upper()
    from .canonical_endpoints import REGISTRY
    definition = REGISTRY.get(endpoint_id)
    if definition is None:
        return "UNKNOWN_ENDPOINT"
    species = str(getattr(snapshot, "species", "") or "").strip().upper()
    if not species or species in {"UNSPECIFIED", "UNKNOWN", "OTHER"}:
        return "UNKNOWN_SPECIES"
    required_species = str(definition.species_requirement or "").strip().upper()
    if required_species and species != required_species:
        return "ENDPOINT_SPECIES_MISMATCH"
    for field in ("model_id", "model_version", "model_artifact_hash"):
        value = str(getattr(snapshot, field, "") or "").strip()
        if not value or value.upper() in {UNKNOWN_PROVENANCE, "UNKNOWN", "NONE", "NULL"}:
            return f"MISSING_{field.upper()}"
    mode = str(getattr(snapshot, "prediction_mode", "") or "").strip().upper()
    if mode not in VALID_PREDICTION_MODES:
        return "INVALID_PREDICTION_MODE"
    value = getattr(snapshot, "value", None)
    classification = str(getattr(snapshot, "classification", "") or "").strip()
    if value is None and not classification:
        return "MISSING_VALUE"
    if not str(getattr(snapshot, "unit", "") or "").strip():
        return "MISSING_UNIT"
    return None


def ensure_stable_core_schema(engine) -> None:
    """Apply the reversible additive Stable Core schema migration."""
    tables = set(inspect(engine).get_table_names())
    if "projects" not in tables:
        return
    Base.metadata.create_all(bind=engine, tables=[
        ScientificMutationAudit.__table__, HistoricalPrediction.__table__,
        ExperimentalObservation.__table__, CurrentPredictionSnapshot.__table__,
    ])
    project_columns = {row["name"] for row in inspect(engine).get_columns("projects")}
    with engine.begin() as connection:
        additions = {
            "lifecycle_status": "VARCHAR(30) NOT NULL DEFAULT 'ACTIVE'",
            "protection_policy": "VARCHAR(40) NOT NULL DEFAULT 'REAL_PROJECT'",
            "archived_at": "DATETIME",
        }
        for name, definition in additions.items():
            if name not in project_columns:
                connection.execute(text(f"ALTER TABLE projects ADD COLUMN {name} {definition}"))
        connection.execute(text("UPDATE projects SET lifecycle_status='ACTIVE' WHERE lifecycle_status IS NULL OR trim(lifecycle_status)=''"))
        connection.execute(text("UPDATE projects SET protection_policy=CASE WHEN is_test_fixture=1 THEN 'SYNTHETIC_TEST' ELSE 'REAL_PROJECT' END WHERE protection_policy IS NULL OR trim(protection_policy)=''"))
        connection.execute(text("UPDATE projects SET protection_policy='PROTECTED_REAL_PROJECT' WHERE id IN (1,3,5,300)"))
        connection.execute(text("CREATE INDEX IF NOT EXISTS ix_projects_lifecycle_status ON projects(lifecycle_status)"))
        connection.execute(text("CREATE INDEX IF NOT EXISTS ix_projects_protection_policy ON projects(protection_policy)"))
        connection.execute(text("CREATE INDEX IF NOT EXISTS ix_experimental_observation_key ON experimental_observations(compound_version_id,canonical_endpoint,species,context_identity,curation_status)"))
        connection.execute(text("CREATE INDEX IF NOT EXISTS ix_current_prediction_key ON current_prediction_snapshots(compound_version_id,canonical_endpoint,species,context_identity,engine_release,is_current)"))
        connection.execute(text("DROP TRIGGER IF EXISTS stable_core_protect_real_project_delete"))
        connection.execute(text("""
            CREATE TRIGGER stable_core_protect_real_project_delete
            BEFORE DELETE ON projects
            WHEN OLD.protection_policy != 'SYNTHETIC_TEST' AND OLD.is_test_fixture != 1
            BEGIN
              SELECT RAISE(ABORT, 'STABLE_CORE: real projects must be archived, not deleted');
            END
        """))
        connection.execute(text("DROP TRIGGER IF EXISTS stable_core_protect_historical_run_delete"))
        connection.execute(text("""
            CREATE TRIGGER stable_core_protect_historical_run_delete
            BEFORE DELETE ON prediction_runs
            WHEN EXISTS (
              SELECT 1 FROM historical_predictions h
              WHERE h.legacy_prediction_run_id=OLD.id AND h.immutable=1
                AND COALESCE((SELECT p.protection_policy FROM projects p WHERE p.id=h.project_id_snapshot),'REAL_PROJECT')!='SYNTHETIC_TEST'
            )
            BEGIN
              SELECT RAISE(ABORT, 'STABLE_CORE: historical predictions are immutable');
            END
        """))
        connection.execute(text("DROP TRIGGER IF EXISTS stable_core_protect_historical_run_update"))
        connection.execute(text("""
            CREATE TRIGGER stable_core_protect_historical_run_update
            BEFORE UPDATE ON prediction_runs
            WHEN EXISTS (
              SELECT 1 FROM historical_predictions h
              WHERE h.legacy_prediction_run_id=OLD.id AND h.immutable=1
                AND COALESCE((SELECT p.protection_policy FROM projects p WHERE p.id=h.project_id_snapshot),'REAL_PROJECT')!='SYNTHETIC_TEST'
            )
            BEGIN
              SELECT RAISE(ABORT, 'STABLE_CORE: historical predictions are immutable');
            END
        """))
        connection.execute(text("DROP TRIGGER IF EXISTS stable_core_protect_real_compound_delete"))
        connection.execute(text("""
            CREATE TRIGGER stable_core_protect_real_compound_delete
            BEFORE DELETE ON compounds
            WHEN EXISTS (
              SELECT 1 FROM projects p
              WHERE p.id=OLD.project_id AND p.protection_policy!='SYNTHETIC_TEST'
            )
            BEGIN
              SELECT RAISE(ABORT, 'STABLE_CORE: compounds in real projects must be archived, not deleted');
            END
        """))
        connection.execute(text("DROP TRIGGER IF EXISTS stable_core_protect_history_update"))
        connection.execute(text("""
            CREATE TRIGGER stable_core_protect_history_update
            BEFORE UPDATE ON historical_predictions
            WHEN OLD.immutable=1 AND COALESCE((SELECT p.protection_policy FROM projects p WHERE p.id=OLD.project_id_snapshot),'REAL_PROJECT')!='SYNTHETIC_TEST'
            BEGIN
              SELECT RAISE(ABORT, 'STABLE_CORE: historical predictions are immutable');
            END
        """))
        connection.execute(text("DROP TRIGGER IF EXISTS stable_core_protect_history_delete"))
        connection.execute(text("""
            CREATE TRIGGER stable_core_protect_history_delete
            BEFORE DELETE ON historical_predictions
            WHEN OLD.immutable=1 AND COALESCE((SELECT p.protection_policy FROM projects p WHERE p.id=OLD.project_id_snapshot),'REAL_PROJECT')!='SYNTHETIC_TEST'
            BEGIN
              SELECT RAISE(ABORT, 'STABLE_CORE: historical predictions are immutable');
            END
        """))
        connection.execute(text("DROP TRIGGER IF EXISTS stable_core_mirror_prediction_run_insert"))
        connection.execute(text("""
            CREATE TRIGGER stable_core_mirror_prediction_run_insert
            AFTER INSERT ON prediction_runs
            WHEN EXISTS (
              SELECT 1 FROM compound_versions cv JOIN compounds c ON c.id=cv.compound_row_id
              JOIN projects p ON p.id=c.project_id
              WHERE cv.id=NEW.version_id AND p.protection_policy!='SYNTHETIC_TEST'
            )
            BEGIN
              INSERT OR IGNORE INTO historical_predictions (
                legacy_prediction_run_id,compound_version_id_snapshot,compound_row_id_snapshot,project_id_snapshot,
                project_name_snapshot,compound_label_snapshot,structure_revision_snapshot,canonical_endpoint,species,
                context_json,context_identity,stage,model_id,model_version,model_artifact_hash,engine_version,
                inputs_hash,outputs_json,provenance_json,applicability_domain_json,uncertainty_json,prediction_mode,
                predicted_at,immutable
              )
              SELECT NEW.id,NEW.version_id,cv.compound_row_id,c.project_id,p.name,COALESCE(NULLIF(c.name,''),c.compound_id),
                cv.inchikey||'@v'||cv.version_number,'MULTI_ENDPOINT','UNSPECIFIED','{}',
                '44136fa355b3678a1146ad16f7e8649e94fb4fc21fe77e8310c060f61caaff8a',NEW.stage,NEW.model_name,
                COALESCE(NULLIF(NEW.model_version,''),'UNKNOWN_PROVENANCE'),
                COALESCE(json_extract(NEW.provenance_json,'$.model_artifact_hash'),json_extract(NEW.provenance_json,'$.artifact_hash'),'UNKNOWN_PROVENANCE'),
                COALESCE(json_extract(NEW.provenance_json,'$.engine_id'),json_extract(NEW.provenance_json,'$.engine_version'),'UNKNOWN_PROVENANCE'),
                NEW.inputs_hash,NEW.outputs_json,NEW.provenance_json,'{}','{}',
                COALESCE(json_extract(NEW.provenance_json,'$.prediction_mode'),'UNKNOWN_MODE'),NEW.created_at,1
              FROM compound_versions cv JOIN compounds c ON c.id=cv.compound_row_id JOIN projects p ON p.id=c.project_id
              WHERE cv.id=NEW.version_id;
              INSERT INTO scientific_mutation_audit (
                occurred_at,actor,process,operation,object_type,object_id,project_id_snapshot,transaction_id,
                before_identity_json,after_identity_json,reason,outcome
              ) SELECT CURRENT_TIMESTAMP,'DATABASE_CONNECTION','SQLITE_TRIGGER','CREATE','HistoricalPrediction',
                CAST(NEW.id AS TEXT),c.project_id,'','{}',json_object('prediction_run_id',NEW.id,'version_id',NEW.version_id),
                'Immutable mirror created from PredictionRun','COMMITTED'
              FROM compound_versions cv JOIN compounds c ON c.id=cv.compound_row_id WHERE cv.id=NEW.version_id;
            END
        """))
        connection.execute(text("DROP TRIGGER IF EXISTS stable_core_audit_project_archive"))
        connection.execute(text("""
            CREATE TRIGGER stable_core_audit_project_archive
            AFTER UPDATE OF lifecycle_status ON projects
            WHEN OLD.lifecycle_status != NEW.lifecycle_status
            BEGIN
              INSERT INTO scientific_mutation_audit (
                occurred_at,actor,process,operation,object_type,object_id,project_id_snapshot,transaction_id,
                before_identity_json,after_identity_json,reason,outcome
              ) VALUES (CURRENT_TIMESTAMP,'DATABASE_CONNECTION','SQLITE_TRIGGER','LIFECYCLE_CHANGE','Project',
                CAST(NEW.id AS TEXT),NEW.id,'',json_object('status',OLD.lifecycle_status),json_object('status',NEW.lifecycle_status),
                'Project lifecycle transition','COMMITTED');
            END
        """))
        connection.execute(text("DROP TRIGGER IF EXISTS stable_core_audit_evidence_state"))
        connection.execute(text("""
            CREATE TRIGGER stable_core_audit_evidence_state
            AFTER UPDATE OF evidence_state,qualification_status,lifecycle_status ON external_experimental_evidence
            WHEN OLD.evidence_state != NEW.evidence_state OR OLD.qualification_status != NEW.qualification_status
                 OR OLD.lifecycle_status != NEW.lifecycle_status
            BEGIN
              INSERT INTO scientific_mutation_audit (
                occurred_at,actor,process,operation,object_type,object_id,project_id_snapshot,transaction_id,
                before_identity_json,after_identity_json,reason,outcome
              ) VALUES (CURRENT_TIMESTAMP,'DATABASE_CONNECTION','SQLITE_TRIGGER','EVIDENCE_STATE_CHANGE',
                'ExternalExperimentalEvidence',CAST(NEW.id AS TEXT),NULL,'',
                json_object('evidence_state',OLD.evidence_state,'qualification_status',OLD.qualification_status,'lifecycle_status',OLD.lifecycle_status),
                json_object('evidence_state',NEW.evidence_state,'qualification_status',NEW.qualification_status,'lifecycle_status',NEW.lifecycle_status),
                'Scientific evidence lifecycle transition','COMMITTED');
            END
        """))
        connection.execute(text("DROP TRIGGER IF EXISTS stable_core_audit_current_prediction_replace"))
        connection.execute(text("""
            CREATE TRIGGER stable_core_audit_current_prediction_replace
            AFTER UPDATE ON current_prediction_snapshots
            WHEN OLD.value IS NOT NEW.value OR OLD.is_current != NEW.is_current OR OLD.engine_release != NEW.engine_release
            BEGIN
              INSERT INTO scientific_mutation_audit (
                occurred_at,actor,process,operation,object_type,object_id,project_id_snapshot,transaction_id,
                before_identity_json,after_identity_json,reason,outcome
              ) VALUES (CURRENT_TIMESTAMP,'DATABASE_CONNECTION','SQLITE_TRIGGER','REPLACE','CurrentPredictionSnapshot',
                CAST(NEW.id AS TEXT),NULL,'',json_object('value',OLD.value,'engine',OLD.engine_release,'current',OLD.is_current),
                json_object('value',NEW.value,'engine',NEW.engine_release,'current',NEW.is_current),
                'Current prediction replacement','COMMITTED');
            END
        """))


def migrate_legacy_scientific_records(connection) -> dict[str, int]:
    """Idempotently materialize authoritative records from legacy sources."""
    migrated_history = 0
    for row in connection.execute(text("""
        SELECT pr.*, cv.compound_row_id, cv.version_number, cv.inchikey,
               c.project_id, c.compound_id, c.name AS compound_name, p.name AS project_name
        FROM prediction_runs pr
        JOIN compound_versions cv ON cv.id=pr.version_id
        JOIN compounds c ON c.id=cv.compound_row_id
        JOIN projects p ON p.id=c.project_id
        WHERE NOT EXISTS (SELECT 1 FROM historical_predictions h WHERE h.legacy_prediction_run_id=pr.id)
        ORDER BY pr.id
    """)).mappings():
        provenance = row["provenance_json"] if isinstance(row["provenance_json"], dict) else json.loads(row["provenance_json"] or "{}")
        engine_version = provenance.get("engine_id") or provenance.get("engine_version") or UNKNOWN_PROVENANCE
        connection.execute(text("""
            INSERT INTO historical_predictions (
              legacy_prediction_run_id,compound_version_id_snapshot,compound_row_id_snapshot,project_id_snapshot,
              project_name_snapshot,compound_label_snapshot,structure_revision_snapshot,canonical_endpoint,species,
              context_json,context_identity,stage,model_id,model_version,model_artifact_hash,engine_version,
              inputs_hash,outputs_json,provenance_json,applicability_domain_json,uncertainty_json,prediction_mode,
              predicted_at,immutable
            ) VALUES (
              :run_id,:version_id,:compound_row_id,:project_id,:project_name,:compound_label,:revision,
              'MULTI_ENDPOINT','UNSPECIFIED','{}',:empty_context,:stage,:model_id,:model_version,:artifact_hash,
              :engine_version,:inputs_hash,:outputs,:provenance,'{}','{}',:prediction_mode,:predicted_at,1
            )
        """), {
            "run_id": row["id"], "version_id": row["version_id"], "compound_row_id": row["compound_row_id"],
            "project_id": row["project_id"], "project_name": row["project_name"],
            "compound_label": row["compound_name"] or row["compound_id"],
            "revision": f"{row['inchikey']}@v{row['version_number']}", "empty_context": context_identity({}),
            "stage": row["stage"], "model_id": row["model_name"] or UNKNOWN_PROVENANCE,
            "model_version": row["model_version"] or UNKNOWN_PROVENANCE,
            "artifact_hash": provenance.get("model_artifact_hash") or provenance.get("artifact_hash") or UNKNOWN_PROVENANCE,
            "engine_version": engine_version, "inputs_hash": row["inputs_hash"],
            "outputs": json.dumps(row["outputs_json"] if isinstance(row["outputs_json"], dict) else json.loads(row["outputs_json"] or "{}")),
            "provenance": json.dumps(provenance), "prediction_mode": provenance.get("prediction_mode", UNKNOWN_PREDICTION_MODE),
            "predicted_at": row["created_at"],
        })
        migrated_history += 1

    migrated_evidence = 0
    for row in connection.execute(text("""
        SELECT * FROM external_experimental_evidence e
        WHERE NOT EXISTS (
          SELECT 1 FROM experimental_observations o
          WHERE o.legacy_record_type='external_experimental_evidence' AND o.legacy_record_id=e.id
        ) ORDER BY e.id
    """)).mappings():
        context = row["assay_conditions_json"] if isinstance(row["assay_conditions_json"], dict) else json.loads(row["assay_conditions_json"] or "{}")
        species = normalize_species_code(row["species"], context)
        try:
            numeric = float(row["normalized_value"] or row["raw_value"])
        except (TypeError, ValueError):
            numeric = None
        accepted = row["accepted_at"] is not None or row["evidence_state"] in {
            "AUTO_QUALIFIED_EXTERNAL", "EXTERNAL_IMPORTED", "INTERNAL_EXPERIMENTAL",
        }
        qualification = row["qualification_json"] if isinstance(row["qualification_json"], dict) else json.loads(row["qualification_json"] or "{}")
        connection.execute(text("""
            INSERT INTO experimental_observations (
              compound_version_id,canonical_endpoint,species,context_json,context_identity,value_text,numeric_value,
              unit,qualifier,lower_bound,upper_bound,source,source_url,source_record_id,reference_id,
              identity_confidence,curation_status,display_comparable,numeric_pairable,learning_eligible,
              retrieved_at,accepted_at,provenance_json,legacy_record_type,legacy_record_id,created_at
            ) VALUES (
              :version_id,:endpoint,:species,:context,:context_id,:value_text,:numeric,:unit,:qualifier,NULL,NULL,
              :source,:source_url,:source_record_id,:reference_id,:identity,:status,:display,:pairable,:learning,
              :retrieved,:accepted_at,:provenance,'external_experimental_evidence',:legacy_id,:created_at
            )
        """), {
            "version_id": row["compound_version_id"], "endpoint": row["canonical_endpoint_id"] or "UNRESOLVED",
            "species": species, "context": json.dumps(canonical_context(context)), "context_id": context_identity(context),
            "value_text": row["normalized_value"] or row["raw_value"], "numeric": numeric,
            "unit": row["normalized_unit"] or row["raw_unit"], "qualifier": row["raw_relation"] or "=",
            "source": row["source_database"], "source_url": row["source_url"], "source_record_id": row["source_record_id"],
            "reference_id": row["source_document_id"] or row["reference_text"], "identity": row["identity_match_status"],
            "status": "ACCEPTED" if accepted else "CANDIDATE", "display": bool(qualification.get("stages", {}).get("ENDPOINT_QUALIFIED")),
            "pairable": bool(qualification.get("stages", {}).get("PREDICTION_PAIRABLE")),
            "learning": bool(qualification.get("adaptation_eligibility", False)), "retrieved": row["retrieved_at"],
            "accepted_at": row["accepted_at"] or (row["imported_at"] if accepted else None),
            "provenance": json.dumps({"legacy_provenance_key": row["provenance_key"], "qualification": qualification}),
            "legacy_id": row["id"], "created_at": row["imported_at"],
        })
        migrated_evidence += 1

    migrated_current = 0
    rows = connection.execute(text("""
        SELECT s.*, cv.inchikey, cv.version_number
        FROM prediction_endpoint_snapshots s
        JOIN compound_versions cv ON cv.id=s.compound_version_id
        WHERE s.base_value IS NOT NULL
        ORDER BY s.compound_version_id,s.endpoint_name,s.created_at DESC,s.id DESC
    """)).mappings()
    seen: set[tuple] = set()
    for row in rows:
        snapshot = row["snapshot_json"] if isinstance(row["snapshot_json"], dict) else json.loads(row["snapshot_json"] or "{}")
        provenance = snapshot.get("provenance") if isinstance(snapshot.get("provenance"), dict) else {}
        species = normalize_species_code(snapshot.get("species"), snapshot)
        mapped = canonicalize_prediction_endpoint(row["endpoint_name"], species=species, route=snapshot.get("route"), context=snapshot)
        context = canonical_context(snapshot)
        engine_release = provenance.get("engine_id") or snapshot.get("engine_id") or snapshot.get("engine_release") or UNKNOWN_PROVENANCE
        key = (row["compound_version_id"], mapped["canonical_endpoint_id"], species, context_identity(context), engine_release)
        if key in seen:
            continue
        seen.add(key)
        exists = connection.execute(text("""
          SELECT 1 FROM current_prediction_snapshots WHERE compound_version_id=:v AND canonical_endpoint=:e
          AND species=:s AND context_identity=:c AND engine_release=:g
        """), {"v": key[0], "e": key[1], "s": key[2], "c": key[3], "g": key[4]}).first()
        if exists:
            continue
        connection.execute(text("""
          INSERT INTO current_prediction_snapshots (
            compound_version_id,canonical_endpoint,species,context_json,context_identity,engine_release,value,unit,
            classification,model_id,model_version,model_artifact_hash,applicability_domain_json,uncertainty_json,
            prediction_mode,source_artifact_type,source_artifact_id,structure_revision,created_at,replaced_at,is_current
          ) VALUES (:v,:e,:s,:j,:c,:g,:value,:unit,'',:model,:model_version,:artifact_hash,:ad,:uncertainty,
                    :mode,'prediction_endpoint_snapshots',:source_id,:revision,:created_at,NULL,1)
        """), {
            "v": row["compound_version_id"], "e": mapped["canonical_endpoint_id"], "s": species,
            "j": json.dumps(context), "c": context_identity(context), "g": engine_release,
            "value": row["project_value"] if row["project_value"] is not None else row["base_value"],
            "unit": row["project_unit"] if row["project_value"] is not None else row["base_unit"],
            "model": provenance.get("engine_name") or snapshot.get("source") or UNKNOWN_PROVENANCE,
            "model_version": provenance.get("engine_version") or UNKNOWN_PROVENANCE,
            "artifact_hash": provenance.get("artifact_hash") or provenance.get("model_artifact_hash") or UNKNOWN_PROVENANCE,
            "ad": json.dumps(snapshot.get("applicability_domain") or {}),
            "uncertainty": json.dumps(snapshot.get("uncertainty") or {}),
            "mode": snapshot.get("prediction_mode") or UNKNOWN_PREDICTION_MODE, "source_id": row["id"],
            "revision": f"{row['inchikey']}@v{row['version_number']}", "created_at": row["created_at"],
        })
        migrated_current += 1
    return {"historical_predictions": migrated_history, "experimental_observations": migrated_evidence, "current_predictions": migrated_current}


def migrate_stable_core_v1_002(connection) -> dict[str, int]:
    """Normalize verified legacy observations without guessing provenance.

    The FDA Orforglipron clearance record explicitly says that route is not
    established.  Three existing user measurements are percent-remaining
    assays and therefore remain scientifically distinct from Clint.
    """
    corrected_orforglipron_cl = 0
    source = connection.execute(text("""
        SELECT * FROM external_experimental_evidence
        WHERE source_record_id='FDA-FOUNDAYO-LABEL-12.3-CL'
          AND canonical_endpoint_id='HUMAN_PK_CL_IV'
    """)).mappings().first()
    if source:
        context = source["assay_conditions_json"] if isinstance(source["assay_conditions_json"], dict) else json.loads(source["assay_conditions_json"] or "{}")
        context["route"] = "UNSPECIFIED"
        context["route_source"] = "UNRESOLVED"
        canonical = canonical_context(context)
        connection.execute(text("""
            UPDATE external_experimental_evidence
            SET canonical_endpoint_id='HUMAN_PK_CL_UNSPECIFIED', assay_conditions_json=:context,
                routing_reason='Stable Core v1: primary source does not establish route',
                updated_at=CURRENT_TIMESTAMP
            WHERE id=:id
        """), {"id": source["id"], "context": json.dumps(context)})
        connection.execute(text("""
            UPDATE experimental_observations
            SET canonical_endpoint='HUMAN_PK_CL_UNSPECIFIED', context_json=:context,
                context_identity=:context_id,
                provenance_json=json_set(COALESCE(provenance_json,'{}'),'$.stable_core_correction','SOURCE_VERIFIED_ROUTE_UNSPECIFIED')
            WHERE legacy_record_type='external_experimental_evidence' AND legacy_record_id=:id
        """), {"id": source["id"], "context": json.dumps(canonical), "context_id": context_identity(context)})
        connection.execute(text("""
            INSERT INTO scientific_mutation_audit (
              occurred_at,actor,process,operation,object_type,object_id,project_id_snapshot,
              transaction_id,before_identity_json,after_identity_json,reason,outcome
            ) VALUES (CURRENT_TIMESTAMP,'stable-core-v1-migration','migrate_stable_core_v1_002',
              'NORMALIZE','ExperimentalObservation',:id,300,'stable-core-v1-002',
              json_object('endpoint','HUMAN_PK_CL_IV','route','IV'),
              json_object('endpoint','HUMAN_PK_CL_UNSPECIFIED','route','UNSPECIFIED'),
              'FDA source describes systemic clearance but does not establish route','COMMITTED')
        """), {"id": str(source["id"])})
        corrected_orforglipron_cl = 1

    migrated_admet = 0
    for row in connection.execute(text("""
        SELECT * FROM admet_measurements m
        WHERE NOT EXISTS (
          SELECT 1 FROM experimental_observations o
          WHERE o.legacy_record_type='admet_measurements' AND o.legacy_record_id=m.id
        ) ORDER BY m.id
    """)).mappings():
        species = normalize_species_code(row["species"], {"matrix": row["matrix"]})
        endpoint = {
            "HUMAN": "HUMAN_MICROSOMAL_STABILITY_PERCENT_REMAINING",
            "RAT": "RAT_MICROSOMAL_STABILITY_PERCENT_REMAINING",
        }.get(species, "UNRESOLVED") if str(row["method"] or "").strip().lower() == "% remaining" else "UNRESOLVED"
        context = canonical_context({
            "matrix": row["matrix"], "assay_type": row["method"],
            "measurement_type": "PERCENT_REMAINING", "species": species,
        })
        provenance = row["provenance_json"] if isinstance(row["provenance_json"], dict) else json.loads(row["provenance_json"] or "{}")
        connection.execute(text("""
            INSERT INTO experimental_observations (
              compound_version_id,canonical_endpoint,species,context_json,context_identity,value_text,
              numeric_value,unit,qualifier,lower_bound,upper_bound,source,source_url,source_record_id,
              reference_id,identity_confidence,curation_status,display_comparable,numeric_pairable,
              learning_eligible,retrieved_at,accepted_at,provenance_json,legacy_record_type,legacy_record_id,created_at
            ) VALUES (:version_id,:endpoint,:species,:context,:context_id,:value_text,:numeric,:unit,
              :qualifier,NULL,NULL,:source,'',:source_record_id,'','EXACT_INTERNAL','ACCEPTED',1,0,0,
              NULL,:accepted_at,:provenance,'admet_measurements',:legacy_id,:created_at)
        """), {
            "version_id": row["version_id"], "endpoint": endpoint, "species": species,
            "context": json.dumps(context), "context_id": context_identity(context),
            "value_text": str(row["value"] if row["value"] is not None else row["qualitative_value"] or ""),
            "numeric": row["value"], "unit": row["unit"] or "", "qualifier": row["qualifier"] or "=",
            "source": row["source"] or "Internal experimental measurement",
            "source_record_id": f"admet_measurements:{row['id']}",
            "accepted_at": row["created_at"],
            "provenance": json.dumps({"legacy_provenance": provenance, "stable_core_semantics": "PERCENT_REMAINING_NOT_CLINT"}),
            "legacy_id": row["id"], "created_at": row["created_at"],
        })
        migrated_admet += 1
    # Legacy AUTO_QUALIFIED status is not accepted here: forensic review found
    # that its PDF parser can mistake nearby times and percentages for PK
    # values.  Stable Core acceptance requires explicit curation or the new
    # deterministic qualification contract below.
    accepted_autoqualified = 0
    return {
        "corrected_orforglipron_clearance": corrected_orforglipron_cl,
        "migrated_admet_measurements": migrated_admet,
        "accepted_deterministically_qualified_observations": accepted_autoqualified,
    }


def migrate_stable_core_v1_003(connection) -> dict[str, int]:
    """Remove provenance claims that were inferred rather than persisted.

    The legacy snapshot and PredictionRun contracts did not require a
    prediction mode.  Stable Core therefore records UNKNOWN_MODE unless the
    originating artifact explicitly stored ASSISTED, HYBRID, or
    FULL_PREDICTION.
    """
    # This is the only controlled migration allowed to amend the newly-built
    # canonical mirror.  The source PredictionRun remains untouched.  The
    # migration driver reinstalls both immutable triggers before returning.
    connection.execute(text("DROP TRIGGER IF EXISTS stable_core_protect_history_update"))
    connection.execute(text("DROP TRIGGER IF EXISTS stable_core_protect_history_delete"))
    current = connection.execute(text("""
        UPDATE current_prediction_snapshots
        SET prediction_mode='UNKNOWN_MODE'
        WHERE source_artifact_type='prediction_endpoint_snapshots'
          AND prediction_mode='FULL_PREDICTION'
          AND NOT EXISTS (
            SELECT 1 FROM prediction_endpoint_snapshots s
            WHERE s.id=current_prediction_snapshots.source_artifact_id
              AND json_extract(s.snapshot_json,'$.prediction_mode') IN ('ASSISTED','HYBRID','FULL_PREDICTION')
          )
    """)).rowcount
    history = connection.execute(text("""
        UPDATE historical_predictions
        SET prediction_mode='UNKNOWN_MODE'
        WHERE prediction_mode='FULL_PREDICTION'
          AND COALESCE(json_extract(provenance_json,'$.prediction_mode'),'') NOT IN ('ASSISTED','HYBRID','FULL_PREDICTION')
    """)).rowcount
    return {
        "current_prediction_modes_corrected": current,
        "historical_prediction_modes_corrected": history,
    }


def migrate_stable_core_v1_004(connection) -> dict[str, int]:
    """Fail closed on legacy automatic evidence qualification.

    No source row is deleted or rewritten. Canonical observations that lack
    explicit acceptance and the Stable Core deterministic acceptance marker
    remain durable CANDIDATE evidence.
    """
    demoted = connection.execute(text("""
      UPDATE experimental_observations
      SET curation_status='CANDIDATE', accepted_at=NULL,
          learning_eligible=0,
          provenance_json=json_set(COALESCE(provenance_json,'{}'),'$.stable_core_acceptance','LEGACY_AUTO_QUALIFICATION_REQUIRES_REVIEW')
      WHERE legacy_record_type='external_experimental_evidence'
        AND curation_status='ACCEPTED'
        AND EXISTS (
          SELECT 1 FROM external_experimental_evidence e
          WHERE e.id=experimental_observations.legacy_record_id
            AND e.evidence_state='AUTO_QUALIFIED_EXTERNAL'
            AND e.accepted_at IS NULL
            AND COALESCE(json_extract(e.qualification_json,'$.stable_core_deterministic_accept'),0) != 1
        )
    """)).rowcount
    return {"legacy_autoqualified_observations_demoted": demoted}


def migrate_stable_core_v1_005(connection) -> dict[str, int]:
    """Verify that curated source semantics survive application restarts.

    Version 005 makes no inferred scientific claim.  It records whether the
    source-backed Orforglipron systemic-clearance correction established by
    migration 002 is present in both the legacy compatibility row and the
    canonical observation.  Application startup no longer runs legacy
    re-indexers, so this invariant is stable after migration.
    """
    source = connection.execute(text("""
      SELECT id,canonical_endpoint_id,assay_conditions_json
      FROM external_experimental_evidence
      WHERE source_record_id='FDA-FOUNDAYO-LABEL-12.3-CL'
    """)).mappings().first()
    if not source:
        return {"orforglipron_clearance_verified": 0, "missing_source_record": 1}
    raw_context = source["assay_conditions_json"]
    context = raw_context if isinstance(raw_context, dict) else json.loads(raw_context or "{}")
    canonical = connection.execute(text("""
      SELECT canonical_endpoint,context_json
      FROM experimental_observations
      WHERE legacy_record_type='external_experimental_evidence'
        AND legacy_record_id=:legacy_id
    """), {"legacy_id": source["id"]}).mappings().first()
    canonical_context_json = canonical["context_json"] if canonical else {}
    canonical_context_value = (
        canonical_context_json
        if isinstance(canonical_context_json, dict)
        else json.loads(canonical_context_json or "{}")
    )
    verified = (
        source["canonical_endpoint_id"] == "HUMAN_PK_CL_UNSPECIFIED"
        and context.get("route") == "UNSPECIFIED"
        and canonical is not None
        and canonical["canonical_endpoint"] == "HUMAN_PK_CL_UNSPECIFIED"
        and canonical_context_value.get("route") == "UNSPECIFIED"
    )
    if not verified:
        raise RuntimeError(
            "Stable Core v1 source-semantics invariant failed for the curated "
            "Orforglipron systemic-clearance observation"
        )
    return {"orforglipron_clearance_verified": 1, "missing_source_record": 0}


def sync_external_observation(connection, evidence) -> None:
    """Publish one legacy evidence row into the canonical observation store."""
    if not inspect(connection).has_table("experimental_observations"):
        return
    raw_context = evidence.assay_conditions_json if isinstance(evidence.assay_conditions_json, dict) else {}
    context = canonical_context(raw_context)
    species = normalize_species_code(evidence.species, raw_context)
    try:
        numeric = float(evidence.normalized_value or evidence.raw_value)
    except (TypeError, ValueError):
        numeric = None
    qualification = evidence.qualification_json if isinstance(evidence.qualification_json, dict) else {}
    stages = qualification.get("stages", {}) if isinstance(qualification.get("stages"), dict) else {}
    deterministic_auto_accept = (
        evidence.evidence_state == "AUTO_QUALIFIED_EXTERNAL"
        and evidence.identity_match_status in {"EXACT_STRUCTURE_MATCH", "EXACT_MATCH"}
        and qualification.get("stable_core_deterministic_accept") is True
        and (
            evidence.qualification_status in {"ENDPOINT_QUALIFIED", "PREDICTION_PAIRABLE"}
            or bool(stages.get("ENDPOINT_QUALIFIED"))
            or bool(stages.get("PREDICTION_PAIRABLE"))
        )
    )
    accepted = bool(evidence.accepted_at) or evidence.evidence_state in {
        "EXTERNAL_IMPORTED", "INTERNAL_EXPERIMENTAL",
    } or deterministic_auto_accept
    rejected = evidence.lifecycle_status in {"INVALIDATED", "DELETED"}
    status = "REJECTED" if rejected else "ACCEPTED" if accepted else "CANDIDATE"
    values = {
        "version_id": evidence.compound_version_id,
        "endpoint": evidence.canonical_endpoint_id or "UNRESOLVED",
        "species": species,
        "context": json.dumps(context),
        "context_id": context_identity(context),
        "value_text": evidence.normalized_value or evidence.raw_value or "",
        "numeric": numeric,
        "unit": evidence.normalized_unit or evidence.raw_unit or "",
        "qualifier": evidence.raw_relation or "=",
        "source": evidence.source_database or "",
        "source_url": evidence.source_url or "",
        "source_record_id": evidence.source_record_id or "",
        "reference_id": evidence.source_document_id or evidence.reference_text or "",
        "identity": evidence.identity_match_status or "UNKNOWN",
        "status": status,
        "display": bool(qualification.get("stages", {}).get("ENDPOINT_QUALIFIED")) or accepted,
        "pairable": bool(qualification.get("stages", {}).get("PREDICTION_PAIRABLE")),
        "learning": bool(qualification.get("adaptation_eligibility", False)) and accepted and not rejected,
        "retrieved": evidence.retrieved_at,
        "accepted_at": evidence.accepted_at or (evidence.imported_at if accepted else None),
        "provenance": json.dumps({
            "legacy_provenance_key": evidence.provenance_key,
            "qualification": qualification,
            "evidence_origin": evidence.evidence_origin,
            "lifecycle_status": evidence.lifecycle_status,
        }),
        "legacy_id": evidence.id,
        "created_at": evidence.imported_at or evidence.first_seen_at or utcnow(),
    }
    existing = connection.execute(text("""
        SELECT id FROM experimental_observations
        WHERE legacy_record_type='external_experimental_evidence' AND legacy_record_id=:legacy_id
    """), {"legacy_id": evidence.id}).first()
    if existing:
        connection.execute(text("""
            UPDATE experimental_observations SET
              compound_version_id=:version_id, canonical_endpoint=:endpoint, species=:species,
              context_json=:context, context_identity=:context_id, value_text=:value_text,
              numeric_value=:numeric, unit=:unit, qualifier=:qualifier, source=:source,
              source_url=:source_url, source_record_id=:source_record_id, reference_id=:reference_id,
              identity_confidence=:identity, curation_status=:status,
              display_comparable=:display, numeric_pairable=:pairable,
              learning_eligible=:learning, retrieved_at=:retrieved, accepted_at=:accepted_at,
              provenance_json=:provenance
            WHERE legacy_record_type='external_experimental_evidence' AND legacy_record_id=:legacy_id
        """), values)
    else:
        connection.execute(text("""
            INSERT INTO experimental_observations (
              compound_version_id,canonical_endpoint,species,context_json,context_identity,value_text,numeric_value,
              unit,qualifier,lower_bound,upper_bound,source,source_url,source_record_id,reference_id,
              identity_confidence,curation_status,display_comparable,numeric_pairable,learning_eligible,
              retrieved_at,accepted_at,provenance_json,legacy_record_type,legacy_record_id,created_at
            ) VALUES (:version_id,:endpoint,:species,:context,:context_id,:value_text,:numeric,:unit,:qualifier,
              NULL,NULL,:source,:source_url,:source_record_id,:reference_id,:identity,:status,:display,:pairable,
              :learning,:retrieved,:accepted_at,:provenance,'external_experimental_evidence',:legacy_id,:created_at)
        """), values)


def sync_admet_measurement_observation(connection, measurement) -> None:
    """Publish an internal ADMET measurement through ExperimentalObservation."""
    if not inspect(connection).has_table("experimental_observations"):
        return
    endpoint_row = connection.execute(text(
        "SELECT name FROM admet_endpoints WHERE id=:id"
    ), {"id": measurement.endpoint_id}).mappings().first()
    if not endpoint_row:
        return
    species = normalize_species_code(measurement.species, {"matrix": measurement.matrix})
    mapped = canonicalize_prediction_endpoint(
        endpoint_row["name"],
        species=species,
        context={"matrix": measurement.matrix, "assay_type": measurement.method},
    )
    endpoint = mapped["canonical_endpoint_id"]
    # Percent-remaining measurements must never be silently recast as Clint.
    if str(measurement.method or "").strip().lower() == "% remaining":
        endpoint = {
            "HUMAN": "HUMAN_MICROSOMAL_STABILITY_PERCENT_REMAINING",
            "RAT": "RAT_MICROSOMAL_STABILITY_PERCENT_REMAINING",
        }.get(species, "UNRESOLVED")
    context = canonical_context({
        "matrix": measurement.matrix, "assay_type": measurement.method,
        "measurement_type": "PERCENT_REMAINING" if str(measurement.method or "").strip().lower() == "% remaining" else endpoint_row["name"],
    })
    values = {
        "version_id": measurement.version_id, "endpoint": endpoint, "species": species,
        "context": json.dumps(context), "context_id": context_identity(context),
        "value_text": str(measurement.value if measurement.value is not None else measurement.qualitative_value or ""),
        "numeric": measurement.value, "unit": measurement.unit or "", "qualifier": measurement.qualifier or "=",
        "source": measurement.source or "Internal experimental measurement",
        "source_record_id": f"admet_measurements:{measurement.id}",
        "provenance": json.dumps({"legacy_provenance": measurement.provenance_json or {}, "stable_core_writer": "ADMETMeasurement"}),
        "legacy_id": measurement.id, "created_at": measurement.created_at or utcnow(),
    }
    connection.execute(text("""
      INSERT INTO experimental_observations (
        compound_version_id,canonical_endpoint,species,context_json,context_identity,value_text,numeric_value,
        unit,qualifier,lower_bound,upper_bound,source,source_url,source_record_id,reference_id,
        identity_confidence,curation_status,display_comparable,numeric_pairable,learning_eligible,
        retrieved_at,accepted_at,provenance_json,legacy_record_type,legacy_record_id,created_at
      ) VALUES (:version_id,:endpoint,:species,:context,:context_id,:value_text,:numeric,:unit,:qualifier,
        NULL,NULL,:source,'',:source_record_id,'','EXACT_INTERNAL','ACCEPTED',1,0,0,NULL,:created_at,
        :provenance,'admet_measurements',:legacy_id,:created_at)
      ON CONFLICT(legacy_record_type,legacy_record_id) DO UPDATE SET
        compound_version_id=excluded.compound_version_id,canonical_endpoint=excluded.canonical_endpoint,
        species=excluded.species,context_json=excluded.context_json,context_identity=excluded.context_identity,
        value_text=excluded.value_text,numeric_value=excluded.numeric_value,unit=excluded.unit,
        qualifier=excluded.qualifier,source=excluded.source,provenance_json=excluded.provenance_json
    """), values)


def sync_activity_measurement_observation(connection, measurement) -> None:
    """Publish target activity without collapsing IC50/Ki/EC50/Kd semantics."""
    if not inspect(connection).has_table("experimental_observations"):
        return
    assay = connection.execute(text("""
        SELECT measurement_type,species,target,cell_line,assay_category FROM assay_definitions WHERE id=:id
    """), {"id": measurement.assay_id}).mappings().first()
    if not assay:
        return
    endpoint = {
        "IC50": "ACTIVITY_IC50", "EC50": "ACTIVITY_EC50",
        "KI": "ACTIVITY_KI", "KD": "ACTIVITY_KD",
    }.get(str(assay["measurement_type"] or "").strip().upper(), "UNRESOLVED")
    species = normalize_species_code(assay["species"], {})
    context = canonical_context({
        "assay_type": assay["assay_category"], "measurement_type": assay["measurement_type"],
        "matrix": assay["cell_line"], "analyte": assay["target"],
    })
    values = {
        "version_id": measurement.version_id, "endpoint": endpoint, "species": species,
        "context": json.dumps(context), "context_id": context_identity(context),
        "value_text": str(measurement.raw_value), "numeric": measurement.normalized_value_nm,
        "unit": "nM", "qualifier": measurement.qualifier or "=", "source": measurement.source or "User experimental",
        "source_record_id": f"activity_measurements:{measurement.id}",
        "provenance": json.dumps({"legacy_provenance": measurement.provenance_json or {}, "original_value": measurement.raw_value, "original_unit": measurement.original_unit}),
        "legacy_id": measurement.id, "created_at": measurement.created_at or utcnow(),
    }
    connection.execute(text("""
      INSERT INTO experimental_observations (
        compound_version_id,canonical_endpoint,species,context_json,context_identity,value_text,numeric_value,
        unit,qualifier,lower_bound,upper_bound,source,source_url,source_record_id,reference_id,
        identity_confidence,curation_status,display_comparable,numeric_pairable,learning_eligible,
        retrieved_at,accepted_at,provenance_json,legacy_record_type,legacy_record_id,created_at
      ) VALUES (:version_id,:endpoint,:species,:context,:context_id,:value_text,:numeric,:unit,:qualifier,
        NULL,NULL,:source,'',:source_record_id,'','EXACT_INTERNAL','ACCEPTED',1,1,0,NULL,:created_at,
        :provenance,'activity_measurements',:legacy_id,:created_at)
      ON CONFLICT(legacy_record_type,legacy_record_id) DO UPDATE SET
        compound_version_id=excluded.compound_version_id,canonical_endpoint=excluded.canonical_endpoint,
        species=excluded.species,context_json=excluded.context_json,context_identity=excluded.context_identity,
        value_text=excluded.value_text,numeric_value=excluded.numeric_value,unit=excluded.unit,
        qualifier=excluded.qualifier,source=excluded.source,provenance_json=excluded.provenance_json
    """), values)
def publish_legacy_endpoint_snapshots(db, version_id: int) -> dict[str, int]:
    """Publish calculation artifacts through the sole current-prediction store."""
    from .admet import PredictionEndpointSnapshot

    artifacts = list(db.scalars(select(PredictionEndpointSnapshot).where(
        PredictionEndpointSnapshot.compound_version_id == version_id,
        PredictionEndpointSnapshot.base_value.is_not(None),
    ).order_by(PredictionEndpointSnapshot.created_at.desc(), PredictionEndpointSnapshot.id.desc())))
    published = replaced = rejected = 0
    seen: set[tuple] = set()
    for artifact in artifacts:
        snapshot = artifact.snapshot_json if isinstance(artifact.snapshot_json, dict) else {}
        provenance = snapshot.get("provenance") if isinstance(snapshot.get("provenance"), dict) else {}
        species = normalize_species_code(snapshot.get("species"), snapshot)
        mapped = canonicalize_prediction_endpoint(artifact.endpoint_name, species=species, route=snapshot.get("route"), context=snapshot)
        # A source species that contradicts an endpoint's canonical species
        # must be rejected, never silently rewritten (e.g. HUMAN_PPB as RAT).
        mapped_species = str(mapped.get("species") or "UNSPECIFIED").upper()
        if species != mapped_species:
            rejected += 1
            continue
        context = canonical_context(snapshot)
        context_id = context_identity(context)
        engine_release = provenance.get("engine_id") or snapshot.get("engine_id") or snapshot.get("engine_release") or UNKNOWN_PROVENANCE
        key = (mapped["canonical_endpoint_id"], mapped_species, context_id, engine_release)
        if key in seen:
            continue
        seen.add(key)
        current = db.scalar(select(CurrentPredictionSnapshot).where(
            CurrentPredictionSnapshot.compound_version_id == version_id,
            CurrentPredictionSnapshot.canonical_endpoint == key[0],
            CurrentPredictionSnapshot.species == key[1],
            CurrentPredictionSnapshot.context_identity == key[2],
            CurrentPredictionSnapshot.engine_release == key[3],
        ))
        values = {
            "value": artifact.project_value if artifact.project_value is not None else artifact.base_value,
            "unit": artifact.project_unit if artifact.project_value is not None else artifact.base_unit,
            "model_id": provenance.get("engine_name") or snapshot.get("source") or UNKNOWN_PROVENANCE,
            "model_version": provenance.get("engine_version") or UNKNOWN_PROVENANCE,
            "model_artifact_hash": provenance.get("artifact_hash") or provenance.get("model_artifact_hash") or UNKNOWN_PROVENANCE,
            "applicability_domain_json": snapshot.get("applicability_domain") or {},
            "uncertainty_json": snapshot.get("uncertainty") or {},
            "prediction_mode": snapshot.get("prediction_mode") or UNKNOWN_PREDICTION_MODE,
            "source_artifact_id": artifact.id,
            "created_at": artifact.created_at,
        }
        if current:
            if current.source_artifact_id != artifact.id:
                for field, value in values.items():
                    setattr(current, field, value)
                replaced += 1
            continue
        db.add(CurrentPredictionSnapshot(
            compound_version_id=version_id, canonical_endpoint=key[0], species=key[1], context_json=context,
            context_identity=key[2], engine_release=key[3], classification="", source_artifact_type="prediction_endpoint_snapshots",
            structure_revision="", is_current=True, **values,
        ))
        published += 1
    if published or replaced:
        db.flush()
    return {"published": published, "replaced": replaced, "rejected": rejected}
