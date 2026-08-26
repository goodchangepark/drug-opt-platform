from pathlib import Path

import pytest
from sqlalchemy import create_engine, select
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend import main
from backend.admet import (
    ADMETConsensusPrediction,
    ADMETModelComparison,
    ADMETModelPerformance,
    ADMETModelRegistry,
    ADMETPrediction,
    ensure_admet_schema,
)
from backend.database import Base
from backend.main import (
    create_admet_measurement,
    create_compound,
    create_project,
    list_admet,
    run_admet_predictions,
    run_compound_prediction_workflow,
)
from backend.metabolism import ensure_metabolism_schema
from backend.models import ensure_ui_schema
from backend.schemas import CompoundCreate, ProjectCreate


ROOT = Path(__file__).parents[1]


@pytest.fixture()
def db():
    engine = create_engine("sqlite://", connect_args={"check_same_thread": False}, poolclass=StaticPool)
    Base.metadata.create_all(engine)
    ensure_ui_schema(engine)
    ensure_admet_schema(engine)
    ensure_metabolism_schema(engine)
    session = sessionmaker(bind=engine, autoflush=False, autocommit=False)()
    try:
        yield session
    finally:
        session.close()


def compound(db, project_name="Prediction UX", compound_id="C001"):
    project = create_project(ProjectCreate(name=project_name, target="Test target"), db)
    row = create_compound(project.id, CompoundCreate(
        compound_id=compound_id, name=compound_id, smiles="CCOc1ccccc1", calculate=False,
    ), db)
    return project, row


def test_save_and_predict_workflow_excludes_activity_and_isolates_failures(db, monkeypatch):
    _, row = compound(db)
    called = []
    monkeypatch.setattr(main, "calculate_compound_properties", lambda row_id, session: called.append(("properties", row_id)) or row)
    monkeypatch.setattr(main, "run_admet_predictions", lambda version_id, session: {
        "status": "PARTIAL", "message": "one endpoint unavailable",
        "endpoint_statuses": [{"endpoint": "Solubility", "status": "COMPLETE"},
                              {"endpoint": "BCRP inhibitor", "status": "MODEL_UNAVAILABLE"}],
        "consensus_predictions": [{"endpoint": "Solubility"}],
    })
    monkeypatch.setattr(main, "run_metabolism_predictions", lambda version_id, session: {
        "status": "COMPLETE", "message": "soft spots and metabolites stored",
    })
    result = run_compound_prediction_workflow(row["row_id"], db)
    assert result["status"] == "PARTIAL"
    assert result["activity_excluded"] is True
    assert result["steps"]["activity"]["status"] == "NOT_INCLUDED"
    assert result["steps"]["admet"]["endpoints"][1]["status"] == "MODEL_UNAVAILABLE"
    assert called == [("properties", row["row_id"])]


def test_multi_model_storage_single_endpoint_consensus_and_cache(db, monkeypatch):
    project, row = compound(db)
    for model in db.scalars(select(ADMETModelRegistry)):
        model.is_active = False
    models = [
        ADMETModelRegistry(endpoint_name="Solubility", model_name="Model A", model_version="1",
            implementation_status="READY", output_unit="log10(mol/L)", source="test", training_dataset="set A",
            validation_json={"rmse": 0.7}, license="test", model_priority=10, ensemble_eligible=True,
            species="Not species-specific", output_type="regression", is_active=True),
        ADMETModelRegistry(endpoint_name="Solubility", model_name="Model B", model_version="2",
            implementation_status="READY", output_unit="log10(mol/L)", source="test", training_dataset="set B",
            validation_json={"rmse": 0.8}, license="test", model_priority=20, ensemble_eligible=True,
            species="Not species-specific", output_type="regression", is_active=True),
    ]
    db.add_all(models); db.commit()
    values = iter([-3.0, -2.0])
    monkeypatch.setattr(main, "model_files_available", lambda endpoint: (True, ""))
    monkeypatch.setattr(main, "predict_endpoint", lambda smiles, endpoint: {
        "status": "COMPLETE", "predicted_value": next(values), "unit": "log10(mol/L)",
        "confidence": "MEDIUM", "applicability_domain": {"classification": "IN_DOMAIN"},
        "uncertainty": 0.2, "uncertainty_reason": "test",
    })
    version_id = row["version"]["id"]
    result = run_admet_predictions(version_id, db)
    assert len(result["predictions"]) == 2
    assert len(result["consensus_predictions"]) == 1
    consensus = result["consensus_predictions"][0]
    assert -3.0 < consensus["combined_value"] < -2.0
    assert len(consensus["models"]) == 2
    assert pytest.approx(sum(item["weight"] for item in consensus["models"])) == 1.0
    assert db.query(ADMETPrediction).count() == 2
    assert db.query(ADMETConsensusPrediction).count() == 1
    cached = run_admet_predictions(version_id, db)
    assert cached["status"] == "CACHED" and db.query(ADMETPrediction).count() == 2

    create_admet_measurement(project.id, {"version_id": version_id, "endpoint": "Solubility",
        "value": -2.4, "unit": "log10(mol/L)", "matrix": "aqueous"}, db)
    assert db.query(ADMETModelComparison).count() == 2
    project_metrics = db.scalars(select(ADMETModelPerformance).where(
        ADMETModelPerformance.scope_key == f"PROJECT:{project.id}"
    )).all()
    assert len(project_metrics) == 2 and all(item.sample_size == 1 for item in project_metrics)


def test_single_model_consensus_and_project_compound_isolation(db, monkeypatch):
    first, row = compound(db, "First project", "FIRST")
    second, other = compound(db, "Second project", "SECOND")
    for model in db.scalars(select(ADMETModelRegistry)):
        model.is_active = model.endpoint_name == "Solubility"
    db.commit()
    monkeypatch.setattr(main, "model_files_available", lambda endpoint: (True, ""))
    monkeypatch.setattr(main, "predict_endpoint", lambda smiles, endpoint: {
        "status": "COMPLETE", "predicted_value": -2.5, "unit": "log10(mol/L)",
        "confidence": "LOW", "applicability_domain": {"classification": "BORDERLINE"},
        "uncertainty": 0.3, "uncertainty_reason": "test",
    })
    result = run_admet_predictions(row["version"]["id"], db)
    assert len(result["consensus_predictions"]) == 1
    assert result["consensus_predictions"][0]["models"][0]["weight"] == pytest.approx(1.0)
    assert list_admet(second.id, db)["predictions"] == []
    assert list_admet(second.id, db)["consensus_predictions"] == []
    assert other["version"]["id"] != row["version"]["id"]


def test_prediction_admet_metabolism_and_optimization_ui_contract():
    source = (ROOT / "frontend/static/app.js").read_text()
    ordering = [source.index(text) for text in (
        "1 · EXPERIMENTAL RESULTS", "2 · PREDICTION RESULTS", "3 · EXPERIMENTAL VS PREDICTION",
        "4 · INTEGRATED PROFILE", "5 · MODEL / PROVENANCE DETAILS",
    )]
    assert ordering == sorted(ordering)
    for text in (
        "Save & Predict", "Activity is intentionally excluded", "Individual Models:", "Combined Prediction",
        "Unavailable models (", "Metabolic Stability · Human / Rat / Mouse Liver Microsomes",
        "Supporting ADME Evidence · Permeability and Plasma Protein Binding (PPB)",
        "Optimization Workspace", "Step 1 — Select Project", "Step 2 — Select Compound",
        "Step 3 — Select Optimization Goal", "Analyze Optimization Strategy", "Generate analogs",
    ):
        assert text in source
    assert ("const tabs=['overview','properties','activity','admet','metabolism','history']" in source or
            "const tabs=['overview','properties','activity','admet','metabolism','pk','history']" in source)
