from types import SimpleNamespace

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend import admet_predictor
from backend.admet import ADMETPrediction, ensure_admet_schema
from backend.admet_predictor import (
    MODEL_VERSION,
    applicability_domain,
    comparable_experimental,
    model_files_available,
    predict_endpoint,
)
from backend.database import Base
from backend.main import (
    create_admet_measurement,
    create_compound,
    create_project,
    list_admet,
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


def _compound(db, smiles="CC(=O)Oc1ccccc1C(=O)O"):
    project = create_project(ProjectCreate(name="Stage 3A"), db)
    compound = create_compound(project.id, CompoundCreate(compound_id="C001", smiles=smiles), db)
    return project.id, compound["version"]["id"]


def test_model_registry_preserves_stage3a_endpoint_specific_models(db):
    project_id, _ = _compound(db)
    registry = list_admet(project_id, db)["models"]
    ready = {row["endpoint"]: row for row in registry if row["active"]}
    assert {"Solubility", "Permeability"}.issubset(ready)
    assert all(ready[name]["model_version"] == MODEL_VERSION for name in ("Solubility", "Permeability"))
    assert ready["Solubility"]["output_unit"] == "log10(mol/L)"
    assert ready["Permeability"]["output_unit"] == "log10(cm/s)"
    assert "direction" in ready["Permeability"]["details"]["limitations"].lower()


def test_prediction_cache_and_compatible_experimental_comparison(db):
    project_id, version_id = _compound(db)
    create_admet_measurement(project_id, {
        "version_id": version_id, "endpoint": "Solubility", "value": 0.01,
        "unit": "mol/L", "matrix": "aqueous", "method": "shake flask",
    }, db)
    create_admet_measurement(project_id, {
        "version_id": version_id, "endpoint": "Caco-2 permeability", "value": 8.4,
        "unit": "10^-6 cm/s", "matrix": "Caco-2", "method": "A→B assay",
    }, db)
    first = run_admet_predictions(version_id, db)
    assert first["status"] == "COMPLETE"
    assert first["cache_hit"] is False
    assert {row["endpoint"] for row in first["predictions"]}.issuperset({"Solubility", "Permeability"})
    first_ids = {row["id"] for row in first["predictions"]}
    assert all(row["confidence"] in {"MEDIUM", "LOW"} for row in first["predictions"])
    assert all(row["applicability_domain"] in {"IN_DOMAIN", "BORDERLINE", "OUT_OF_DOMAIN"} for row in first["predictions"])

    second = run_admet_predictions(version_id, db)
    assert second["status"] == "CACHED" and second["cache_hit"] is True
    assert {row["id"] for row in second["predictions"]} == first_ids
    assert db.query(ADMETPrediction).count() == len(first_ids)

    listing = list_admet(project_id, db)
    assert len(listing["predictions"]) == len(first_ids)
    for prediction in [row for row in listing["predictions"] if row["endpoint"] in {"Solubility", "Permeability"}]:
        comparison = prediction["experimental_comparisons"][0]
        assert comparison["absolute_error"] >= 0
        assert comparison["relative_error_percent_linear_scale"] >= 0


@pytest.mark.parametrize("endpoint,name,matrix,unit,value,expected", [
    ("Solubility", "Solubility", "aqueous", "µM", 10.0, -5.0),
    ("Permeability", "Permeability", "Caco-2", "10^-6 cm/s", 8.4, -5.075721),
    ("Permeability", "Caco-2 permeability", "Caco-2", "µm/s", 1.0, -4.0),
])
def test_unit_handling(endpoint, name, matrix, unit, value, expected):
    measurement = SimpleNamespace(
        mean_value=None, value=value, qualifier="=", unit=unit, matrix=matrix, method="assay",
    )
    converted, _ = comparable_experimental(endpoint, measurement, name)
    assert converted == pytest.approx(expected, abs=1e-6)


def test_incompatible_pH_pampa_and_mdck_are_not_compared():
    base = dict(mean_value=None, value=1.0, qualifier="=", unit="µM", matrix="aqueous", method="assay")
    assert comparable_experimental("Solubility", SimpleNamespace(**base), "Solubility pH 7.4")[0] is None
    for matrix in ("PAMPA", "MDCK"):
        measurement = SimpleNamespace(**{**base, "unit": "10^-6 cm/s", "matrix": matrix})
        assert comparable_experimental("Permeability", measurement, f"{matrix} permeability")[0] is None


def test_applicability_domain_exact_training_and_out_of_domain_compound():
    # First Caco-2 training structure is necessarily an exact nearest neighbour.
    with (admet_predictor.MODEL_ROOT / "caco2" / "training.csv").open() as stream:
        stream.readline()
        exact_smiles = stream.readline().split(",")[1]
    exact = applicability_domain(exact_smiles, "Permeability")
    assert exact["classification"] == "IN_DOMAIN"
    assert exact["nearest_training_similarity"] == 1.0

    out = applicability_domain("C" * 200, "Solubility")
    assert out["classification"] == "OUT_OF_DOMAIN"
    assert out["descriptors_outside_range"]


def test_model_unavailable_is_explicit_and_never_fabricates(tmp_path, monkeypatch):
    monkeypatch.setattr(admet_predictor, "MODEL_ROOT", tmp_path)
    available, reason = model_files_available("Solubility")
    assert available is False and "missing" in reason.lower()
    result = predict_endpoint("CCO", "Solubility")
    assert result["status"] == "MODEL_UNAVAILABLE"
    assert "reason" in result and "predicted_value" not in result


def test_public_reference_compound_directionality_sanity():
    soluble = predict_endpoint("CCO", "Solubility")["predicted_value"]
    hydrophobic = predict_endpoint("Clc1ccc(cc1)c2ccc(Cl)cc2", "Solubility")["predicted_value"]
    assert soluble > hydrophobic + 4.0

    permeable = predict_endpoint("CCO", "Permeability")["predicted_value"]
    low_permeability = predict_endpoint("O=C(O)C(O)C(O)C(O)CO", "Permeability")["predicted_value"]
    assert permeable > low_permeability + 1.0
