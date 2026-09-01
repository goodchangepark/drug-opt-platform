"""v4.1 persistence invariants for the user-facing scientific workspace."""

from types import SimpleNamespace

from sqlalchemy import create_engine, select
from sqlalchemy.orm import sessionmaker

from backend.database import Base
from backend.main import _persist_harvest_result
from backend.models import Compound, CompoundVersion, ExternalExperimentalEvidence, ExperimentalSearchRun, Project, ensure_ui_schema
from backend.admet import ensure_admet_schema
from backend.metabolism import ensure_metabolism_schema


def _db():
    engine = create_engine("sqlite://")
    Base.metadata.create_all(engine)
    ensure_ui_schema(engine)
    ensure_admet_schema(engine)
    ensure_metabolism_schema(engine)
    return sessionmaker(bind=engine, autoflush=False, autocommit=False)()


def _identity(version):
    return SimpleNamespace(
        inchikey=version.inchikey,
        to_dict=lambda: {"identity_graph_version": "test-identity-v1", "inchikey": version.inchikey},
    )


def _record(group, endpoint, value, *, source_id):
    q = {
        "identity_status": "IDENTITY_QUALIFIED",
        "reference_status": "REFERENCE_QUALIFIED",
        "numeric_status": "NUMERIC_QUALIFIED" if value else "NUMERIC_NOT_QUALIFIED",
        "endpoint_status": "ENDPOINT_QUALIFIED" if endpoint else "ENDPOINT_NOT_QUALIFIED",
        "context_status": "CONTEXT_QUALIFIED" if endpoint else "CONTEXT_NOT_QUALIFIED",
        "stages": {"ENDPOINT_QUALIFIED": bool(endpoint), "CONTEXT_QUALIFIED": bool(endpoint)},
    }
    return {
        "source": "Test source", "source_record_id": source_id, "endpoint": endpoint,
        "value": value, "unit": "% bound" if endpoint else "", "relation": "=",
        "reference": "https://example.test/reference", "reference_status": "REFERENCE_RESOLVED_TEST",
        "identity_match_status": "EXACT_STRUCTURE_MATCH", "display_evidence_group_id": group,
        "independent_experiment_group_id": group + "-independent",
        "display": {"canonical_endpoint_id": "HUMAN_PPB" if endpoint else "UNRESOLVED", "normalized_value": value or None, "normalized_unit": "% bound"},
        "routing": {"section": "ADMET", "canonical_endpoint_id": "HUMAN_PPB", "routing_reason": "test"},
        "qualification": q,
    }


def test_every_displayed_unique_observation_is_persisted_and_repeat_search_is_idempotent():
    db = _db()
    project = Project(name="v41 persistence test")
    db.add(project)
    db.flush()
    compound = Compound(project_id=project.id, compound_id="C1", name="Test compound", cas_number=None)
    db.add(compound)
    db.flush()
    version = CompoundVersion(compound_row_id=compound.id, version_number=1, original_smiles="CCO", canonical_smiles="CCO", isomeric_smiles="CCO", inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
    db.add(version)
    db.flush()
    records = [_record("DISPLAY-a", "PPB", "91.2", source_id="a"), _record("DISPLAY-b", "", "narrative", source_id="b")]
    result = {"records": records, "summary": {"raw_records": 2, "unique_records": 2, "endpoint_qualified": 1, "context_qualified": 1, "importable": 0}, "sources": {"Test source": "COMPLETE"}}
    first = _persist_harvest_result(db, compound, version, _identity(version), result)
    assert first["saved"] is True
    assert first["persisted_candidate_count"] == 2
    assert db.scalar(select(ExternalExperimentalEvidence).where(ExternalExperimentalEvidence.raw_value == "narrative")) is not None
    run = db.scalar(select(ExperimentalSearchRun).order_by(ExperimentalSearchRun.id.desc()))
    assert run.unique_count == 2
    assert run.persisted_observation_count == 2
    assert run.display_only_non_persistent_count == 0

    second_result = {"records": records, "summary": {"raw_records": 2, "unique_records": 2, "endpoint_qualified": 1, "context_qualified": 1, "importable": 0}, "sources": {"Test source": "COMPLETE"}}
    second = _persist_harvest_result(db, compound, version, _identity(version), second_result)
    assert second["persisted_candidate_count"] == 0
    assert second["existing_candidate_count"] == 2
    assert db.query(ExternalExperimentalEvidence).count() == 2
    db.close()


def test_v41_frontend_uses_canonical_db_comparison_on_reload():
    from pathlib import Path

    js = (Path(__file__).parents[1] / "frontend/static/app.js").read_text()
    html = (Path(__file__).parents[1] / "frontend/static/index.html").read_text()
    assert "/scientific-comparison" in js
    assert "persisted_observations" in js
    assert "persistence-v4-1" in html
