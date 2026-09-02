"""v4.4A canonical manual experimental ingestion safeguards."""
from datetime import datetime, timedelta, timezone
from sqlalchemy import create_engine, select
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool
from fastapi.testclient import TestClient

from backend.database import Base, get_db
from backend.endpoint_comparison import build_endpoint_comparison
from backend.evidence_capture import save_internal_evidence
from backend.main import app
from backend.models import Compound, CompoundVersion, ExternalExperimentalEvidence, Project, ensure_ui_schema
from backend.activity_models import AssayDefinition
from backend.admet import ADMETEndpoint, ADMETModelRegistry, ADMETPrediction, ADMETPredictionRun, PredictionExperimentalPairRecord


def _db():
    engine = create_engine("sqlite://", connect_args={"check_same_thread": False}, poolclass=StaticPool)
    Base.metadata.create_all(engine)
    ensure_ui_schema(engine)
    return sessionmaker(bind=engine, autoflush=False, autocommit=False)()


def _client(db):
    def override_get_db():
        yield db
    app.dependency_overrides[get_db] = override_get_db
    return TestClient(app)


def _compound(db):
    project = Project(name="manual evidence test")
    db.add(project); db.flush()
    compound = Compound(project_id=project.id, compound_id="M1", name="Manual", cas_number=None)
    db.add(compound); db.flush()
    version = CompoundVersion(compound_row_id=compound.id, version_number=1, original_smiles="CCO", canonical_smiles="CCO", isomeric_smiles="CCO", inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
    db.add(version); db.commit()
    return project, compound, version


def test_manual_ppb_uses_shared_canonical_evidence_row_and_normalizes_fu():
    db = _db(); project, compound, version = _compound(db)
    evidence = save_internal_evidence(db, project.id, compound.id, {
        "canonical_endpoint_id": "HUMAN_PPB", "raw_value": "0.10", "raw_unit": "fu",
        "species": "Human", "matrix": "plasma", "study_id": "INT-PPB-1"})
    db.commit()
    assert evidence.evidence_state == "INTERNAL_EXPERIMENTAL"
    assert evidence.canonical_endpoint_id == "HUMAN_PPB"
    assert float(evidence.normalized_value) == 90.0
    assert evidence.normalized_unit == "% bound"
    view = build_endpoint_comparison(db, version.id)
    row = next(row for row in view["endpoints"] if row["endpoint_id"] == "HUMAN_PPB")
    assert len(row["experimental_internal"]) == 1
    assert row["experimental_external_candidates"] == []


def test_manual_caco2_without_direction_persists_as_review_not_false_direct_comparison():
    db = _db(); project, compound, version = _compound(db)
    evidence = save_internal_evidence(db, project.id, compound.id, {
        "canonical_endpoint_id": "CACO2_PAPP_AB", "raw_value": "12.4", "raw_unit": "×10^-6 cm/s",
        "matrix": "Caco-2", "study_id": "INT-CACO-1"})
    db.commit()
    assert evidence.evidence_state == "INTERNAL_EXPERIMENTAL"
    assert evidence.comparability_status == "CONDITIONALLY_COMPARABLE"
    assert evidence.qualification_status == "CONTEXT_NOT_QUALIFIED"
    view = build_endpoint_comparison(db, version.id)
    row = next(row for row in view["endpoints"] if row["endpoint_id"] == "CACO2_PAPP_AB")
    assert len(row["needs_review"]) == 1


def test_manual_duplicate_returns_existing_and_revision_preserves_history():
    db = _db(); project, compound, _ = _compound(db)
    payload = {"canonical_endpoint_id": "HUMAN_PPB", "raw_value": "90", "raw_unit": "% bound", "species": "Human", "matrix": "plasma", "study_id": "INT-1"}
    first = save_internal_evidence(db, project.id, compound.id, payload)
    duplicate = save_internal_evidence(db, project.id, compound.id, payload)
    assert duplicate is first
    revised = save_internal_evidence(db, project.id, compound.id, payload | {"raw_value": "91"}, supersedes=first)
    db.commit()
    assert first.lifecycle_status == "SUPERSEDED"
    assert revised.revision_number == 2
    assert revised.supersedes_evidence_id == first.id
    assert db.scalars(select(ExternalExperimentalEvidence)).all()


def test_manual_pk_cmax_preserves_route_dose_and_converts_unit():
    db = _db(); project, compound, version = _compound(db)
    evidence = save_internal_evidence(db, project.id, compound.id, {
        "parameter": "Cmax", "raw_value": "0.42", "raw_unit": "µg/mL",
        "species": "Human", "route": "ORAL", "dose": "200", "dose_unit": "mg",
        "regimen": "Single dose", "analyte": "PARENT", "study_id": "INT-PK-1"})
    db.commit()
    assert evidence.canonical_endpoint_id == "HUMAN_PK_CMAX_ORAL"
    assert float(evidence.normalized_value) == 420.0
    assert evidence.normalized_unit == "ng/mL"
    row = next(row for row in build_endpoint_comparison(db, version.id)["endpoints"] if row["endpoint_id"] == "HUMAN_PK_CMAX_ORAL")
    assert row["experimental_internal"][0]["context"]["dose"] == "200"
    assert row["experimental_internal"][0]["route"] == "ORAL"


def test_manual_cyp_ic50_is_related_without_numeric_comparison():
    db = _db(); project, compound, version = _compound(db)
    evidence = save_internal_evidence(db, project.id, compound.id, {
        "canonical_endpoint_id": "CYP3A4_INHIBITION", "raw_value": "18", "raw_unit": "µM",
        "measurement_type": "IC50", "study_id": "INT-CYP-1"})
    db.commit()
    assert evidence.canonical_endpoint_id == "CYP3A4_INHIBITION"
    assert evidence.comparability_status.startswith("RELATED")
    row = next(row for row in build_endpoint_comparison(db, version.id)["endpoints"] if row["endpoint_id"] == "CYP3A4_INHIBITION")
    assert len(row["related_evidence"]) == 1
    assert row["related_evidence"][0]["comparability"].startswith("RELATED")


def test_manual_activity_requires_matching_project_assay_and_uses_same_row():
    db = _db(); project, compound, version = _compound(db)
    assay = AssayDefinition(project_id=project.id, name="EGFR biochemical", target="EGFR", measurement_type="IC50", unit="nM")
    db.add(assay); db.commit()
    evidence = save_internal_evidence(db, project.id, compound.id, {
        "raw_endpoint": "IC50", "raw_value": "12", "raw_unit": "nM", "assay_id": str(assay.id), "study_id": "INT-ACT-1"})
    db.commit()
    row = next(row for row in build_endpoint_comparison(db, version.id)["endpoints"] if row["endpoint_id"] == f"ACTIVITY_IC50:{assay.id}")
    assert row["experimental_internal"][0]["id"] == evidence.id


def test_manual_api_returns_success_only_after_canonical_row_is_committed():
    db = _db(); project, compound, version = _compound(db)
    client = _client(db)
    try:
        response = client.post(f"/api/projects/{project.id}/compounds/{compound.id}/experimental", json={
            "canonical_endpoint_id": "HUMAN_PPB", "raw_value": "0.096", "raw_unit": "fu",
            "species": "Human", "matrix": "plasma", "study_id": "API-PPB-1"})
        assert response.status_code == 201
        body = response.json()
        assert body["saved"] is True
        evidence_id = body["evidence"]["id"]
        persisted = db.get(ExternalExperimentalEvidence, evidence_id)
        assert persisted and float(persisted.normalized_value) == 90.4
        canonical = client.get(f"/api/compound-versions/{version.id}/scientific-comparison")
        assert canonical.status_code == 200
        row = next(row for row in canonical.json()["scientific_rows"] if row["canonical_endpoint"] == "HUMAN_PPB")
        assert any(item["id"] == evidence_id for item in row["experimental_observations"])
        invalid = client.post(f"/api/projects/{project.id}/compounds/{compound.id}/experimental", json={"raw_value": "1", "raw_unit": "nM"})
        assert invalid.status_code == 400
        assert db.scalars(select(ExternalExperimentalEvidence).where(ExternalExperimentalEvidence.raw_value == "1")).all() == []
    finally:
        app.dependency_overrides.clear()


def test_manual_ppb_creates_eligible_pair_against_immutable_preexperiment_snapshot():
    db = _db(); project, compound, version = _compound(db)
    endpoint = ADMETEndpoint(project_id=project.id, name="Plasma protein binding", preferred_unit="% bound")
    model = ADMETModelRegistry(endpoint_name="Plasma protein binding", model_name="test PPB", implementation_status="READY", output_unit="% bound", is_active=True)
    run = ADMETPredictionRun(version_id=version.id, inputs_hash="frozen-test", status="COMPLETE")
    db.add_all([endpoint, model, run]); db.flush()
    prediction = ADMETPrediction(run_id=run.id, endpoint_id=endpoint.id, version_id=version.id, model_id=model.id,
        predicted_value=91.1, unit="% bound", outputs_json={"prediction_snapshot": {"base_prediction": 91.1}},
        created_at=datetime.now(timezone.utc) - timedelta(minutes=1))
    db.add(prediction); db.commit()
    client = _client(db)
    try:
        response = client.post(f"/api/projects/{project.id}/compounds/{compound.id}/experimental", json={
            "canonical_endpoint_id": "HUMAN_PPB", "raw_value": "90.4", "raw_unit": "% bound",
            "species": "Human", "matrix": "plasma", "study_id": "PAIR-PPB-1"})
        assert response.status_code == 201 and response.json()["pair_created"] is True
        pair = db.scalar(select(PredictionExperimentalPairRecord).where(PredictionExperimentalPairRecord.external_evidence_id == response.json()["evidence"]["id"]))
        assert pair.adaptation_eligibility is True
        assert pair.prediction_record_id == prediction.id
        assert prediction.predicted_value == 91.1
    finally:
        app.dependency_overrides.clear()
