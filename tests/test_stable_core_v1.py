import sqlite3

import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine, delete, select, update
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import sessionmaker

from backend.database import PRODUCTION_DATABASE_PATH, UnsafeDatabaseConfiguration, load_database_settings
from backend.models import Compound, CompoundVersion, ExternalExperimentalEvidence, PredictionRun, Project
from backend.scientific_core_service import compare_scientific_pair
from backend.stable_core import (
    CurrentPredictionSnapshot,
    ExperimentalObservation,
    HistoricalPrediction,
    ScientificMutationAudit,
    context_identity,
    ensure_stable_core_schema,
    admit_current_prediction,
    expected_structure_revision,
    prediction_context_identity,
)
from backend.database import Base
from backend.main import app, bulk_delete_projects, delete_project
from backend.activity_models import ActivityMeasurement, AssayDefinition
from backend.admet import ADMETEndpoint, ADMETMeasurement


def test_test_database_guard_rejects_production_realpath():
    with pytest.raises(UnsafeDatabaseConfiguration, match="TEST DATABASE RESOLVES TO PRODUCTION DATABASE"):
        load_database_settings({
            "DRUGOPT_ENV": "test",
            "DRUGOPT_DATABASE_URL": f"sqlite:///{PRODUCTION_DATABASE_PATH}",
        })


def test_e2e_database_guard_rejects_production_realpath():
    with pytest.raises(UnsafeDatabaseConfiguration, match="TEST DATABASE RESOLVES TO PRODUCTION DATABASE"):
        load_database_settings({
            "DRUGOPT_ENV": "e2e",
            "DRUGOPT_DATABASE_URL": f"sqlite:///{PRODUCTION_DATABASE_PATH}",
        })


def test_pytest_blocks_direct_sqlite_access_to_production_database():
    from backend.database import DATABASE_SETTINGS

    with pytest.raises(RuntimeError, match="TEST DATABASE RESOLVES TO PRODUCTION DATABASE"):
        sqlite3.connect(DATABASE_SETTINGS.production_path)


def test_real_project_and_historical_prediction_are_deletion_proof(tmp_path):
    db_path = tmp_path / "stable-core.sqlite"
    local_engine = create_engine(f"sqlite:///{db_path}")
    Base.metadata.create_all(local_engine)
    ensure_stable_core_schema(local_engine)
    Local = sessionmaker(bind=local_engine)
    with Local() as db:
        project = Project(name="Protected science", protection_policy="PROTECTED_REAL_PROJECT")
        db.add(project); db.flush()
        db.execute(update(Project).where(Project.id == project.id).values(is_test_fixture=False, protection_policy="PROTECTED_REAL_PROJECT"))
        compound = Compound(project_id=project.id, compound_id="P1", name="P1")
        db.add(compound); db.flush()
        version = CompoundVersion(
            compound_row_id=compound.id, version_number=1, original_smiles="CCO", canonical_smiles="CCO",
            isomeric_smiles="CCO", inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        )
        db.add(version); db.flush()
        run = PredictionRun(
            version_id=version.id, stage="test", model_name="fixture", model_version="1",
            inputs_hash="x" * 64, outputs_json={"value": 1}, provenance_json={"engine_id": "historical-engine"},
        )
        db.add(run); db.commit()
        assert db.scalar(select(HistoricalPrediction).where(HistoricalPrediction.legacy_prediction_run_id == run.id))
        with pytest.raises(IntegrityError, match="historical predictions are immutable"):
            db.execute(delete(PredictionRun).where(PredictionRun.id == run.id)); db.commit()
        db.rollback()
        with pytest.raises(IntegrityError, match="historical predictions are immutable"):
            db.execute(update(PredictionRun).where(PredictionRun.id == run.id).values(model_version="rewritten")); db.commit()
        db.rollback()
        with pytest.raises(IntegrityError, match="compounds in real projects must be archived"):
            db.execute(delete(Compound).where(Compound.id == compound.id)); db.commit()
        db.rollback()
        with pytest.raises(IntegrityError, match="real projects must be archived"):
            db.execute(delete(Project).where(Project.id == project.id)); db.commit()


def test_prediction_creation_has_attributable_audit_record(tmp_path):
    db_path = tmp_path / "stable-core-audit.sqlite"
    local_engine = create_engine(f"sqlite:///{db_path}")
    Base.metadata.create_all(local_engine)
    ensure_stable_core_schema(local_engine)
    Local = sessionmaker(bind=local_engine)
    with Local() as db:
        project = Project(name="Audit science", protection_policy="REAL_PROJECT")
        db.add(project); db.flush()
        db.execute(update(Project).where(Project.id == project.id).values(is_test_fixture=False, protection_policy="REAL_PROJECT"))
        compound = Compound(project_id=project.id, compound_id="A1", name="A1")
        db.add(compound); db.flush()
        version = CompoundVersion(compound_row_id=compound.id, version_number=1, original_smiles="CC", canonical_smiles="CC", isomeric_smiles="CC", inchikey="OTMSDBZUPAUEDD-UHFFFAOYSA-N")
        db.add(version); db.flush()
        db.add(PredictionRun(version_id=version.id, stage="test", model_name="m", model_version="1", inputs_hash="y" * 64, outputs_json={}, provenance_json={}))
        db.commit()
        audit = db.scalar(select(ScientificMutationAudit).where(ScientificMutationAudit.object_type == "HistoricalPrediction"))
        assert audit and audit.process == "SQLITE_TRIGGER" and audit.operation == "CREATE"


def test_external_evidence_writer_publishes_candidate_then_accepted(tmp_path):
    local_engine = create_engine(f"sqlite:///{tmp_path / 'stable-core-evidence.sqlite'}")
    Base.metadata.create_all(local_engine)
    ensure_stable_core_schema(local_engine)
    Local = sessionmaker(bind=local_engine)
    with Local() as db:
        project = Project(name="Evidence science")
        db.add(project); db.flush()
        compound = Compound(project_id=project.id, compound_id="E1", name="E1")
        db.add(compound); db.flush()
        version = CompoundVersion(compound_row_id=compound.id, version_number=1, original_smiles="CC", canonical_smiles="CC", isomeric_smiles="CC", inchikey="OTMSDBZUPAUEDD-UHFFFAOYSA-N")
        db.add(version); db.flush()
        evidence = ExternalExperimentalEvidence(
            compound_version_id=version.id, provenance_key="e" * 64, raw_endpoint_name="Cmax",
            raw_value="12", raw_unit="ng/mL", assay_conditions_json={"route": "ORAL", "dose": 10, "dose_unit": "mg"},
            species="Human", source_database="Primary source", source_record_id="source:1",
            identity_match_status="EXACT_STRUCTURE_MATCH", endpoint_match_status="MATCHED",
            canonical_endpoint_id="HUMAN_PK_CMAX_ORAL", normalized_value="12", normalized_unit="ng/mL",
            evidence_state="EXTERNAL_CANDIDATE", evidence_origin="EXPERIMENTAL_EXTERNAL",
        )
        db.add(evidence); db.commit()
        canonical = db.scalar(select(ExperimentalObservation).where(ExperimentalObservation.legacy_record_id == evidence.id))
        assert canonical and canonical.curation_status == "CANDIDATE" and canonical.species == "HUMAN"
        evidence.evidence_state = "EXTERNAL_IMPORTED"
        evidence.accepted_at = evidence.imported_at
        db.commit(); db.refresh(canonical)
        assert canonical.curation_status == "ACCEPTED"


def _pair(endpoint="HUMAN_PK_T_HALF_IV", *, exp_species="HUMAN", pred_species="HUMAN", exp_context=None, pred_context=None, exp_unit="hours", pred_unit="h"):
    experimental = ExperimentalObservation(
        compound_version_id=1, canonical_endpoint=endpoint, species=exp_species,
        context_json=exp_context or {}, context_identity=context_identity(exp_context or {}),
        value_text="10", numeric_value=10, unit=exp_unit, qualifier="=", source="test",
        identity_confidence="EXACT", curation_status="ACCEPTED", display_comparable=True,
        numeric_pairable=True, learning_eligible=False, provenance_json={}, legacy_record_type="test", legacy_record_id=1,
    )
    prediction = CurrentPredictionSnapshot(
        compound_version_id=1, canonical_endpoint=endpoint, species=pred_species,
        context_json=pred_context or {}, context_identity=context_identity(pred_context or {}),
        engine_release="test-engine", value=20, unit=pred_unit, model_id="test", model_version="1",
        model_artifact_hash="hash", applicability_domain_json={}, uncertainty_json={}, prediction_mode="FULL_PREDICTION",
    )
    return experimental, prediction


def test_canonical_comparison_distinguishes_pairability_and_learning_eligibility():
    experimental, prediction = _pair()
    result = compare_scientific_pair(experimental, prediction)
    assert result["display_comparable"] is True
    assert result["numeric_pairable"] is True
    assert result["learning_eligible"] is False
    assert result["fold_error"] == 2
    assert result["comparison_unit"] == "hours"


def test_canonical_comparison_rejects_species_and_context_mismatch():
    experimental, prediction = _pair(pred_species="RAT")
    assert compare_scientific_pair(experimental, prediction)["reason"] == "SPECIES_MISMATCH"
    experimental, prediction = _pair(
        endpoint="HUMAN_PK_CMAX_ORAL",
        exp_context={"route": "ORAL", "dose": 10, "dose_unit": "mg", "formulation": "tablet"},
        pred_context={"route": "ORAL", "dose": 20, "dose_unit": "mg", "formulation": "tablet"},
        exp_unit="ng/mL", pred_unit="ng/mL",
    )
    assert compare_scientific_pair(experimental, prediction)["reason"] == "CONTEXT_MISMATCH:dose"


def test_canonical_comparison_never_pairs_classifier_probability_as_quantity():
    experimental, prediction = _pair(endpoint="CYP3A4_INHIBITOR_CLASS", exp_unit="probability", pred_unit="probability")
    result = compare_scientific_pair(experimental, prediction)
    assert result["numeric_pairable"] is False
    assert result["reason"] == "CLASSIFICATION_OR_NON_CONTINUOUS_ENDPOINT"


def test_internal_admet_and_activity_writers_publish_canonical_observations(tmp_path):
    local_engine = create_engine(f"sqlite:///{tmp_path / 'stable-core-internal.sqlite'}")
    Base.metadata.create_all(local_engine)
    ensure_stable_core_schema(local_engine)
    Local = sessionmaker(bind=local_engine)
    with Local() as db:
        project = Project(name="Internal science")
        db.add(project); db.flush()
        compound = Compound(project_id=project.id, compound_id="I1", name="I1")
        db.add(compound); db.flush()
        version = CompoundVersion(compound_row_id=compound.id, version_number=1, original_smiles="CC", canonical_smiles="CC", isomeric_smiles="CC", inchikey="OTMSDBZUPAUEDD-UHFFFAOYSA-N")
        db.add(version); db.flush()
        endpoint = ADMETEndpoint(project_id=project.id, name="HLM intrinsic clearance", preferred_unit="mL/min/kg")
        assay = AssayDefinition(project_id=project.id, name="Target IC50", measurement_type="IC50", unit="nM", species="Human", target="Target A")
        db.add_all([endpoint, assay]); db.flush()
        db.add(ADMETMeasurement(version_id=version.id, endpoint_id=endpoint.id, species="Human", matrix="microsomes", value=2.5, unit="mL/min/kg", method="intrinsic clearance"))
        db.add(ActivityMeasurement(assay_id=assay.id, version_id=version.id, raw_value=100, original_unit="nM", normalized_value_nm=100))
        db.commit()
        rows = list(db.scalars(select(ExperimentalObservation).where(
            ExperimentalObservation.compound_version_id == version.id,
        ).order_by(ExperimentalObservation.legacy_record_type)))
        assert {(row.legacy_record_type, row.canonical_endpoint, row.curation_status) for row in rows} == {
            ("activity_measurements", "ACTIVITY_IC50", "ACCEPTED"),
            ("admet_measurements", "HLM_CLINT", "ACCEPTED"),
        }


def test_missing_legacy_prediction_mode_is_not_guessed(tmp_path):
    local_engine = create_engine(f"sqlite:///{tmp_path / 'stable-core-mode.sqlite'}")
    Base.metadata.create_all(local_engine)
    ensure_stable_core_schema(local_engine)
    Local = sessionmaker(bind=local_engine)
    with Local() as db:
        project = Project(name="Mode science")
        db.add(project); db.flush()
        db.execute(update(Project).where(Project.id == project.id).values(is_test_fixture=False, protection_policy="REAL_PROJECT"))
        compound = Compound(project_id=project.id, compound_id="M1", name="M1")
        db.add(compound); db.flush()
        version = CompoundVersion(compound_row_id=compound.id, version_number=1, original_smiles="CC", canonical_smiles="CC", isomeric_smiles="CC", inchikey="OTMSDBZUPAUEDD-UHFFFAOYSA-N")
        db.add(version); db.flush()
        run = PredictionRun(version_id=version.id, stage="legacy", model_name="legacy", model_version="1", inputs_hash="z" * 64, outputs_json={}, provenance_json={})
        db.add(run); db.commit()
        history = db.scalar(select(HistoricalPrediction).where(HistoricalPrediction.legacy_prediction_run_id == run.id))
        assert history.prediction_mode == "UNKNOWN_MODE"


def test_recovered_project5_history_and_orforglipron_canonical_pk_contract():
    """Production-shaped isolated fixture protects the repaired scientific identities."""
    client = TestClient(app)
    projects = client.get("/api/projects").json()
    assert {row["id"] for row in projects} >= {1, 3, 5, 300}

    pk = client.get("/api/compound-versions/11/scientific-tabs/pk").json()
    accepted = {
        row["canonical_endpoint"]: row["experimental"]
        for row in pk["rows"] if row.get("experimental")
    }
    assert {"HUMAN_PK_F_ORAL", "HUMAN_PK_VD_IV", "HUMAN_PK_CL_UNSPECIFIED", "HUMAN_PK_CMAX_UNSPECIFIED"} <= set(accepted)
    by_endpoint = {row["canonical_endpoint"]: row for row in pk["rows"] if row.get("experimental")}
    assert by_endpoint["HUMAN_PK_CMAX_UNSPECIFIED"]["context"]["dose"] == 36
    assert by_endpoint["HUMAN_PK_CMAX_UNSPECIFIED"]["context"]["dose_unit"] == "mg"
    assert by_endpoint["HUMAN_PK_CL_UNSPECIFIED"]["context"]["route"] == "UNSPECIFIED"
    assert all(row["species"] == "HUMAN" for row in by_endpoint.values())

    required_runs = {64, 65, 66, 67, 68, 73, 74, 75, 127, 128, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162}
    from backend.database import SessionLocal
    with SessionLocal() as db:
        restored = set(db.scalars(select(HistoricalPrediction.legacy_prediction_run_id).where(
            HistoricalPrediction.project_id_snapshot == 5,
        )))
        assert required_runs <= restored


def test_compound_summary_counts_only_current_engine_snapshots():
    from backend.database import SessionLocal
    from backend.prediction_engine_registry import CURRENT_ENGINE_ID

    client = TestClient(app)
    summary = client.get("/api/compounds/1/summary").json()
    version_ids = [row["version_number"] for row in summary["versions"]]
    assert version_ids
    with SessionLocal() as db:
        compound_version_ids = list(db.scalars(select(CompoundVersion.id).where(
            CompoundVersion.compound_row_id == 1,
        )))
        candidates = list(db.scalars(select(CurrentPredictionSnapshot).where(
            CurrentPredictionSnapshot.compound_version_id.in_(compound_version_ids),
            CurrentPredictionSnapshot.engine_release == CURRENT_ENGINE_ID,
            CurrentPredictionSnapshot.is_current.is_(True),
        )))
        expected = sum(admit_current_prediction(db, row).eligible for row in candidates)
    assert summary["scientific_snapshot"]["current_prediction_snapshot_count"] == expected
    assert summary["prediction_count"] == expected


def test_scientific_row_api_aligns_experiment_and_current_prediction_without_rewriting_history():
    """Exercise the exact UI row contract against the isolated production-shaped DB."""
    from backend.database import SessionLocal
    from backend.model_artifact_authority import artifact_bundle_sha256, model_artifact_registration
    from backend.prediction_engine_registry import CURRENT_ENGINE_ID

    with SessionLocal() as db:
        version = db.get(CompoundVersion, 11)
        context = {}
        observation = ExperimentalObservation(
            compound_version_id=11, canonical_endpoint="MW", species="HUMAN",
            context_json=context, context_identity=prediction_context_identity("MW", context),
            value_text="46.0", numeric_value=46.0, unit="g/mol", qualifier="=",
            source="Stable Core isolated contract fixture", identity_confidence="EXACT",
            curation_status="ACCEPTED", display_comparable=True, numeric_pairable=True,
            learning_eligible=False, provenance_json={"fixture": True},
            legacy_record_type="stable_core_v12_contract", legacy_record_id=1,
        )
        db.add(observation); db.flush()
        observation_id = observation.id
        historical_before = {
            row.legacy_prediction_run_id: (row.engine_version, row.model_version, row.outputs_json)
            for row in db.scalars(select(HistoricalPrediction).where(
                HistoricalPrediction.compound_version_id_snapshot == 11,
            ))
        }
        registration = model_artifact_registration("MW")
        snapshot = CurrentPredictionSnapshot(
            compound_version_id=11,
            canonical_endpoint=observation.canonical_endpoint,
            species=observation.species,
            context_json=observation.context_json,
            context_identity=observation.context_identity,
            engine_release=CURRENT_ENGINE_ID,
            value=47.0,
            unit="g/mol",
            model_id=registration.model_id,
            model_version=registration.model_version,
            model_artifact_hash=artifact_bundle_sha256(registration),
            prediction_mode="FULL_PREDICTION",
            source_artifact_type="PROPERTY_CALCULATION",
            structure_revision=expected_structure_revision(version),
        )
        db.add(snapshot)
        db.commit()
        snapshot_id = snapshot.id

    client = TestClient(app)
    payload = client.get("/api/compound-versions/11/scientific-tabs/properties").json()
    row = next(item for item in payload["rows"] if (item.get("prediction") or {}).get("snapshot_id") == snapshot_id)
    assert row["experimental"]["observation_id"] == observation_id
    assert row["prediction"]["engine_version"] == CURRENT_ENGINE_ID
    assert row["comparison"]["display_comparable"] is True
    assert row["comparison"]["numeric_pairable"] is True
    assert row["comparison"]["fold_error"] == pytest.approx(47.0 / 46.0)

    with SessionLocal() as db:
        current = db.get(CurrentPredictionSnapshot, snapshot_id)
        current.value = 48.0
        db.commit()
        historical_after = {
            item.legacy_prediction_run_id: (item.engine_version, item.model_version, item.outputs_json)
            for item in db.scalars(select(HistoricalPrediction).where(
                HistoricalPrediction.compound_version_id_snapshot == 11,
            ))
        }
        assert historical_after == historical_before


def test_current_engine_and_historical_engine_provenance_are_separate():
    client = TestClient(app)
    health = client.get("/api/health").json()
    assert health["release_status"] in {"PRODUCTION_VALIDATED", "RELEASE_WITH_FINAL_INTEGRITY_BLOCKERS"}
    assert health["rollback_engine"] == "drugopt-prediction-engine-v3@3.3.2"
    assert health["superseded_engine"] == health["rollback_engine"]
    current = client.get("/api/prediction-engine/current").json()["current_production_engine"]
    assert current["engine_id"] == "drugopt-prediction-engine-v3@3.3.3"
    assert current["status"] in {"PRODUCTION_VALIDATED", "RELEASE_WITH_FINAL_INTEGRITY_BLOCKERS"}
    assert current["decision"] in {"STABLE_CORE_V1_2_VALIDATED", "STABLE_CORE_V1_2_VALIDATION_IN_PROGRESS"}
    history = client.get("/api/compound-versions/11/scientific-tabs/history").json()["records"]
    assert history
    assert any(row["engine_version"] != current["engine_id"] for row in history)
    assert all(row["prediction_mode"] in {"ASSISTED", "HYBRID", "FULL_PREDICTION", "UNKNOWN_MODE"} for row in history)
    frontend = (PRODUCTION_DATABASE_PATH.parent / "frontend/static/app.js").read_text(encoding="utf-8")
    assert "sidebar-footer-version'},'v3.3.3'" not in frontend
    assert "currentEngine?.engine_version?'v'+currentEngine.engine_version:'UNKNOWN_CURRENT_ENGINE'" in frontend


def test_real_project_archive_preserves_immutable_history_and_audits_delete_attempt(tmp_path):
    local_engine = create_engine(f"sqlite:///{tmp_path / 'stable-core-archive.sqlite'}")
    Base.metadata.create_all(local_engine)
    ensure_stable_core_schema(local_engine)
    Local = sessionmaker(bind=local_engine)
    with Local() as db:
        project = Project(name="Archive science")
        db.add(project); db.flush()
        db.execute(update(Project).where(Project.id == project.id).values(is_test_fixture=False, protection_policy="REAL_PROJECT"))
        compound = Compound(project_id=project.id, compound_id="ARC", name="ARC")
        db.add(compound); db.flush()
        version = CompoundVersion(compound_row_id=compound.id, version_number=1, original_smiles="CC", canonical_smiles="CC", isomeric_smiles="CC", inchikey="OTMSDBZUPAUEDD-UHFFFAOYSA-N")
        db.add(version); db.flush()
        run = PredictionRun(version_id=version.id, stage="historical", model_name="model", model_version="old", inputs_hash="a" * 64, outputs_json={}, provenance_json={"engine_id": "historical-engine"})
        db.add(run); db.commit()
        archived = delete_project(project.id, {"confirmation_name": project.name}, db)
        assert archived["archived_project_ids"] == [project.id]
        assert db.get(Project, project.id).lifecycle_status == "ARCHIVED"
        assert db.scalar(select(HistoricalPrediction).where(HistoricalPrediction.legacy_prediction_run_id == run.id))
        with pytest.raises(Exception):
            bulk_delete_projects({"projects": [{"id": project.id, "confirmation_name": project.name, "hard_delete": True}]}, db)
        audit = db.scalar(select(ScientificMutationAudit).where(
            ScientificMutationAudit.operation == "DELETE_ATTEMPT",
            ScientificMutationAudit.object_id == str(project.id),
        ))
        assert audit and audit.outcome == "BLOCKED"
