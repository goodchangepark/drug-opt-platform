import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend import admet_predictor
from backend.admet import ADMETPrediction, TRANSPORTER_UNAVAILABLE, ensure_admet_schema
from backend.admet_predictor import (
    MODEL_SPECS, applicability_domain, classification_experimental_evidence,
    comparable_experimental, model_files_available, predict_endpoint,
)
from backend.database import Base
from backend.main import (
    create_admet_measurement, create_compound, create_project, list_admet,
    run_admet_predictions,
)
from backend.schemas import CompoundCreate, ProjectCreate


@pytest.fixture()
def db():
    engine = create_engine("sqlite://", connect_args={"check_same_thread": False}, poolclass=StaticPool)
    Base.metadata.create_all(bind=engine)
    ensure_admet_schema(engine)
    session = sessionmaker(bind=engine, autoflush=False, autocommit=False)()
    try:
        yield session
    finally:
        session.close()


def compound(db, name="Stage 3E", compound_id="T001", smiles="CC(=O)Oc1ccccc1C(=O)O"):
    project = create_project(ProjectCreate(name=name), db)
    created = create_compound(project.id, CompoundCreate(compound_id=compound_id, smiles=smiles), db)
    return project.id, created["version"]["id"]


def measurement(value, unit, endpoint, species="Human", matrix="human P-gp functional assay"):
    return SimpleNamespace(
        id=1, endpoint_id=1, mean_value=None, value=value, qualifier="=", unit=unit,
        species=species, matrix=matrix, method="transporter inhibition assay", notes="", provenance_json={},
    ), endpoint


def test_transporter_registry_activates_only_qualified_human_pgp_inhibitor(db):
    project_id, _ = compound(db)
    models = {row["endpoint"]: row for row in list_admet(project_id, db)["models"]}
    active = models["P-gp inhibitor"]
    assert active["active"] and active["output_unit"] == "probability"
    assert active["details"]["transporter"] == "P-gp / ABCB1"
    assert active["details"]["role"] == "INHIBITOR" and active["details"]["species"] == "Human"
    assert active["details"]["training_n"] == 1275
    for endpoint in TRANSPORTER_UNAVAILABLE:
        row = models[endpoint]
        assert not row["active"] and row["status"] == "MODEL_UNAVAILABLE"
        assert row["details"]["role"] in {"SUBSTRATE", "INHIBITOR"}
        assert row["details"]["species"] == "Human" and row["unavailable_reason"]


def test_pgp_and_bcrp_substrate_inhibitor_endpoints_are_isolated():
    assert "P-gp inhibitor" in MODEL_SPECS and "P-gp substrate" not in MODEL_SPECS
    assert TRANSPORTER_UNAVAILABLE["P-gp substrate"]["role"] == "SUBSTRATE"
    assert TRANSPORTER_UNAVAILABLE["BCRP substrate"]["role"] == "SUBSTRATE"
    assert TRANSPORTER_UNAVAILABLE["BCRP inhibitor"]["role"] == "INHIBITOR"
    result = predict_endpoint("CCO", "P-gp substrate")
    assert result["status"] == "MODEL_UNAVAILABLE" and "predicted_value" not in result


def test_pgp_inhibitor_prediction_is_probability_not_quantitative_potency():
    result = predict_endpoint("CC(=O)Oc1ccccc1C(=O)O", "P-gp inhibitor")
    assert result["status"] == "COMPLETE" and 0 <= result["probability"] <= 1
    assert result["classification"] in {"INHIBITOR", "NON_INHIBITOR"}
    assert result["transporter"] == "P-gp / ABCB1" and result["role"] == "INHIBITOR"
    assert result["species"] == "Human" and result["unit"] == "probability"
    assert "IC50" not in result and "Ki" not in result


def test_experimental_transporter_role_species_and_unit_handling():
    binary, name = measurement(1, "class", "P-gp inhibitor classification")
    assert comparable_experimental("P-gp inhibitor", binary, name)[0] == 1

    substrate, name = measurement(1, "class", "P-gp substrate classification")
    assert comparable_experimental("P-gp inhibitor", substrate, name)[0] is None

    rat, name = measurement(1, "class", "P-gp inhibitor classification", species="Rat", matrix="rat P-gp assay")
    assert comparable_experimental("P-gp inhibitor", rat, name)[0] is None

    ic50, name = measurement(3.0, "µM", "P-gp inhibition IC50")
    assert comparable_experimental("P-gp inhibitor", ic50, name)[0] is None
    evidence = classification_experimental_evidence("P-gp inhibitor", 0.9, [ic50], {1: name})
    assert evidence[0]["evidence_type"] == "QUANTITATIVE"
    assert evidence[0]["comparison"] == "NOT_NUMERICALLY_COMPARABLE"
    assert evidence[0]["absolute_error"] is None and evidence[0]["relative_error_percent"] is None


def test_transporter_applicability_domain_and_confidence_are_endpoint_specific():
    path = admet_predictor.MODEL_ROOT / MODEL_SPECS["P-gp inhibitor"]["model_key"] / "training.csv"
    with path.open() as stream:
        stream.readline()
        exact_smiles = stream.readline().split(",")[0]
    exact = applicability_domain(exact_smiles, "P-gp inhibitor")
    assert exact["nearest_training_similarity"] == 1.0 and exact["classification"] == "IN_DOMAIN"
    out = applicability_domain("C" * 200, "P-gp inhibitor")
    assert out["classification"] == "OUT_OF_DOMAIN" and out["chemical_space_distance"] > 0
    high_probability = predict_endpoint(
        "CC(C)C(CCCN(C)CCC1=CC(=C(C=C1)OC)OC)(C#N)C2=CC(=C(C=C2)OC)OC",
        "P-gp inhibitor",
    )
    assert high_probability["probability"] > 0.99 and high_probability["confidence"] == "LOW"
    assert MODEL_SPECS["P-gp inhibitor"]["independent_validation"]["status"] == "NOT_AVAILABLE"


def test_transporter_experimental_cache_compound_version_and_project_isolation(db):
    project_id, version_id = compound(db, "First")
    other_project, _ = compound(db, "Second")
    create_admet_measurement(project_id, {
        "version_id": version_id, "endpoint": "P-gp inhibitor classification", "species": "Human",
        "matrix": "human P-gp functional assay", "value": 1, "unit": "class",
        "method": "transporter inhibition assay", "source": "Experimental reference",
    }, db)
    first = run_admet_predictions(version_id, db)
    count = db.query(ADMETPrediction).count()
    second = run_admet_predictions(version_id, db)
    assert first["status"] == "COMPLETE" and second["status"] == "CACHED"
    assert db.query(ADMETPrediction).count() == count
    row = next(item for item in list_admet(project_id, db)["predictions"] if item["endpoint"] == "P-gp inhibitor")
    assert row["version_id"] == version_id
    assert row["outputs"]["experimental_evidence"][0]["comparison"] in {"AGREES", "DISAGREES"}
    assert list_admet(other_project, db)["predictions"] == []


def test_transporter_model_assets_are_cpu_compatible_and_sanity_is_honest():
    available, reason = model_files_available("P-gp inhibitor")
    assert available, reason
    artifact = json.loads((Path(__file__).parents[1] / "models/admetica/validation/transporter/pgp_inhibitor_sanity.json").read_text())
    assert artifact["independent_validation"] == "NOT_AVAILABLE" and artifact["metrics"] is None
    assert all(row["direction_correct"] for row in artifact["results"])
    assert any(row["nearest_training_similarity"] == 1.0 for row in artifact["results"])


def test_transporter_ui_has_details_unavailable_and_no_overall_ranking():
    source = (Path(__file__).parents[1] / "frontend/static/app.js").read_text()
    assert "function transporterPredictionTable" in source
    assert "Transporter','Role','Species','Prediction','Probability','Experimental','Domain','Confidence','Model" in source
    assert "unavailableTransporterModels" in source and "MODEL_UNAVAILABLE" in source
    assert "Potential P-gp" not in source  # flag text comes from deterministic backend provenance
    assert "overall candidate score" not in source.lower()
