import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend import admet_predictor
from backend.admet import ADMETPrediction, ensure_admet_schema
from backend.admet_predictor import (MODEL_SPECS, applicability_domain,
                                     comparable_experimental,
                                     cyp_experimental_evidence,
                                     model_files_available, predict_endpoint)
from backend.database import Base
from backend.main import (create_admet_measurement, create_compound,
                          create_project, list_admet,
                          run_admet_predictions)
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


def compound(db, name="Stage 3C", compound_id="C001", smiles="CC(=O)Oc1ccccc1C(=O)O"):
    project = create_project(ProjectCreate(name=name), db)
    created = create_compound(project.id, CompoundCreate(compound_id=compound_id, smiles=smiles), db)
    return project.id, created["version"]["id"]


def experimental(value, unit, endpoint, matrix="human recombinant CYP"):
    return SimpleNamespace(
        id=1, endpoint_id=1, mean_value=None, value=value, qualifier="=", unit=unit,
        matrix=matrix, method="CYP assay", species="Human", notes="", provenance_json={},
    ), endpoint


def test_cyp_registry_has_five_inhibitors_and_only_qualified_substrates(db):
    project_id, _ = compound(db)
    models = {item["endpoint"]: item for item in list_admet(project_id, db)["models"]}
    for isoform in ("CYP1A2", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A4"):
        row = models[f"{isoform} inhibitor"]
        assert row["active"] and row["output_unit"] == "probability"
        assert row["details"]["role"] == "INHIBITOR"
        assert row["details"]["training_n"] > 12000
    for isoform in ("CYP2C9", "CYP2D6", "CYP3A4"):
        assert models[f"{isoform} substrate"]["active"]
    for isoform in ("CYP1A2", "CYP2C19"):
        row = models[f"{isoform} substrate"]
        assert not row["active"] and row["status"] == "MODEL_UNAVAILABLE"
        assert "validation" in row["details"]["reason"]


@pytest.mark.parametrize("isoform", ["CYP1A2", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A4"])
def test_five_cyp_inhibitor_models_return_probabilities_not_ic50(isoform):
    result = predict_endpoint("CC(=O)Oc1ccccc1C(=O)O", f"{isoform} inhibitor")
    assert result["status"] == "COMPLETE"
    assert 0 <= result["probability"] <= 1
    assert result["predicted_value"] == result["probability"]
    assert result["classification"] in {"INHIBITOR", "NON_INHIBITOR"}
    assert result["role"] == "INHIBITOR" and result["isoform"] == isoform
    assert result["unit"] == "probability" and "IC50" not in result


@pytest.mark.parametrize("isoform", ["CYP2C9", "CYP2D6", "CYP3A4"])
def test_cyp_substrate_models_are_role_isolated(isoform):
    result = predict_endpoint("CC(=O)Oc1ccccc1C(=O)O", f"{isoform} substrate")
    assert 0 <= result["probability"] <= 1
    assert result["classification"] in {"SUBSTRATE", "NON_SUBSTRATE"}
    assert result["role"] == "SUBSTRATE"
    assert MODEL_SPECS[f"{isoform} substrate"]["model_key"] != MODEL_SPECS[f"{isoform} inhibitor"]["model_key"]


def test_classification_and_quantitative_experimental_evidence_are_not_mixed():
    classified, name = experimental(1, "class", "CYP3A4 inhibitor classification")
    converted, _ = comparable_experimental("CYP3A4 inhibitor", classified, name)
    assert converted == 1
    assert comparable_experimental("CYP3A4 substrate", classified, name)[0] is None

    ic50, name = experimental(2.5, "µM", "CYP3A4 inhibition IC50")
    assert comparable_experimental("CYP3A4 inhibitor", ic50, name)[0] is None
    evidence = cyp_experimental_evidence("CYP3A4 inhibitor", 0.8, [ic50], {1: name})
    assert evidence[0]["comparison"] == "NOT_NUMERICALLY_COMPARABLE"
    assert evidence[0]["absolute_error"] is None and evidence[0]["relative_error_percent"] is None


def test_cyp_applicability_domain_exact_training_and_out_of_domain():
    path = admet_predictor.MODEL_ROOT / MODEL_SPECS["CYP1A2 inhibitor"]["model_key"] / "training.csv"
    with path.open() as stream:
        stream.readline(); exact_smiles = stream.readline().split(",")[0]
    assert applicability_domain(exact_smiles, "CYP1A2 inhibitor")["nearest_training_similarity"] == 1.0
    out = applicability_domain("C" * 200, "CYP3A4 substrate")
    assert out["classification"] == "OUT_OF_DOMAIN" and out["chemical_space_distance"] > 0


def test_cyp_confidence_uses_validation_and_domain_not_probability():
    result = predict_endpoint("CCN(CC)CCCC(C)NC1=C2C=CC(Cl)=CC2=NC=C1", "CYP3A4 inhibitor")
    assert 0 <= result["probability"] <= 1
    assert result["confidence"] == "LOW"
    independent = MODEL_SPECS["CYP3A4 inhibitor"]["independent_validation"]
    assert independent["balanced_accuracy"] < 0.70


def test_unavailable_substrate_has_no_fake_model_or_prediction():
    assert "CYP1A2 substrate" not in MODEL_SPECS
    available, reason = model_files_available("CYP1A2 substrate")
    assert not available and "No endpoint" in reason
    result = predict_endpoint("CCO", "CYP1A2 substrate")
    assert result["status"] == "MODEL_UNAVAILABLE" and "predicted_value" not in result


def test_cyp_experimental_input_cache_and_project_isolation(db):
    first_project, first_version = compound(db, "First")
    second_project, _ = compound(db, "Second")
    create_admet_measurement(first_project, {
        "version_id": first_version, "endpoint": "CYP3A4 inhibition IC50", "species": "Human",
        "matrix": "human recombinant CYP3A4", "value": 2.5, "unit": "µM",
    }, db)
    create_admet_measurement(first_project, {
        "version_id": first_version, "endpoint": "CYP3A4 inhibitor classification", "species": "Human",
        "matrix": "human recombinant CYP3A4", "value": 1, "unit": "class",
    }, db)
    first = run_admet_predictions(first_version, db)
    second = run_admet_predictions(first_version, db)
    assert first["status"] == "COMPLETE" and second["status"] == "CACHED"
    assert db.query(ADMETPrediction).count() == len(first["predictions"])
    cyp = next(row for row in list_admet(first_project, db)["predictions"] if row["endpoint"] == "CYP3A4 inhibitor")
    evidence = cyp["outputs"]["experimental_evidence"]
    assert {item["evidence_type"] for item in evidence} == {"CLASSIFICATION", "QUANTITATIVE"}
    quantitative = next(item for item in evidence if item["evidence_type"] == "QUANTITATIVE")
    assert quantitative["comparison"] == "NOT_NUMERICALLY_COMPARABLE"
    assert list_admet(second_project, db)["predictions"] == []


def test_independent_validation_artifact_has_true_probability_metrics_and_overlap_audit():
    path = Path(__file__).parents[1] / "models/admetica/validation/cyp/independent_validation.json"
    result = json.loads(path.read_text())
    for isoform in ("CYP2C9", "CYP2D6", "CYP3A4"):
        metrics = result["inhibitors"][isoform]
        assert metrics["canonical_training_overlap_removed"] == 0
        assert {"AUROC", "AUPRC", "balanced_accuracy", "sensitivity", "specificity", "MCC"}.issubset(metrics)
    sanity = result["substrates"]["CYP3A4"]
    assert sanity["canonical_training_overlap_removed"] == 2 and sanity["n"] == 22
    assert sanity["positive_only"] and "AUROC" not in sanity


def test_all_cyp_model_assets_are_arm64_cpu_runtime_compatible():
    for endpoint in MODEL_SPECS:
        if endpoint.startswith("CYP"):
            available, reason = model_files_available(endpoint)
            assert available, reason
