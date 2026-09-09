import re
import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from backend.scientific_core_service import build_scientific_endpoint_rows
from backend.database import Base
from backend.models import Compound, CompoundVersion, Project
from backend.stable_core import CurrentPredictionSnapshot, ensure_stable_core_schema, snapshot_admission_reason


@pytest.fixture
def db_session(tmp_path):
    engine = create_engine(f"sqlite:///{tmp_path / 'stable-core-v11.sqlite'}")
    Base.metadata.create_all(engine)
    ensure_stable_core_schema(engine)
    Session = sessionmaker(bind=engine)
    with Session() as session:
        yield session


@pytest.fixture
def compound_version(db_session):
    project = Project(name="v1.1 fixture", protection_policy="SYNTHETIC_TEST", is_test_fixture=True)
    db_session.add(project)
    db_session.flush()
    compound = Compound(project_id=project.id, compound_id="V11", name="V11")
    db_session.add(compound)
    db_session.flush()
    version = CompoundVersion(
        compound_row_id=compound.id, version_number=1, original_smiles="CCO",
        canonical_smiles="CCO", isomeric_smiles="CCO",
        inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
    )
    db_session.add(version)
    db_session.flush()
    return version


def test_malformed_current_snapshot_is_not_admissible(db_session, compound_version):
    snapshot = CurrentPredictionSnapshot(
        compound_version_id=compound_version.id,
        canonical_endpoint="HUMAN_PPB",
        species="HUMAN",
        context_json={},
        context_identity="x" * 64,
        engine_release="drugopt-prediction-engine-v3@3.3.3",
        value=50.0,
        unit="% bound",
        model_id="UNKNOWN_PROVENANCE",
        model_version="UNKNOWN_PROVENANCE",
        model_artifact_hash="UNKNOWN_PROVENANCE",
        prediction_mode="UNKNOWN_MODE",
        structure_revision="fixture",
    )
    db_session.add(snapshot)
    db_session.commit()

    assert snapshot_admission_reason(snapshot).startswith("MISSING_")
    payload = build_scientific_endpoint_rows(db_session, compound_version.id)
    assert not any((row.get("prediction") or {}).get("snapshot_id") == snapshot.id for row in payload["rows"])


def test_endpoint_species_contradiction_is_not_admissible(db_session, compound_version):
    snapshot = CurrentPredictionSnapshot(
        compound_version_id=compound_version.id,
        canonical_endpoint="HUMAN_PPB",
        species="RAT",
        context_json={},
        context_identity="y" * 64,
        engine_release="drugopt-prediction-engine-v3@3.3.3",
        value=50.0,
        unit="% bound",
        model_id="model",
        model_version="1",
        model_artifact_hash="a" * 64,
        prediction_mode="FULL_PREDICTION",
        structure_revision="fixture",
        is_current=True,
    )
    assert snapshot_admission_reason(snapshot) == "ENDPOINT_SPECIES_MISMATCH"


def test_pk_frontend_has_no_mount_time_simulation_post():
    from pathlib import Path

    js = (Path(__file__).parents[1] / "frontend/static/app.js").read_text(encoding="utf-8")
    # The only simulation POST must be inside the explicit run handler, not
    # the passive loadData effect.
    assert "const autoPayload" not in js
    assert re.search(r"const handleRun\s*=\s*async\(\)\s*=>", js)


def test_legacy_comparison_does_not_publish_current_snapshots():
    from pathlib import Path

    source = (Path(__file__).parents[1] / "backend/endpoint_comparison.py").read_text(encoding="utf-8")
    assert "publication = publish_legacy_endpoint_snapshots" not in source
