import csv
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend import admet_predictor
from backend.admet import ADMETPrediction, SAFETY_UNAVAILABLE, ensure_admet_schema
from backend.admet_predictor import (
    MODEL_SPECS, applicability_domain, classification_experimental_evidence,
    comparable_experimental, model_files_available, predict_endpoint,
)
from backend.database import Base
from backend.main import (
    _integrated_admet_profile, create_admet_measurement, create_compound,
    create_project, list_admet, run_admet_predictions,
)
from backend.schemas import CompoundCreate, ProjectCreate

ROOT = Path(__file__).parents[1]
SAFETY = ("hERG liability", "Ames mutagenicity", "DILI clinical liability")


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


def compound(db, name="Stage 3F", compound_id="SAFE-001", smiles="CC(=O)Oc1ccccc1C(=O)O"):
    project = create_project(ProjectCreate(name=name), db)
    created = create_compound(project.id, CompoundCreate(compound_id=compound_id, smiles=smiles), db)
    return project.id, created["version"]["id"]


def measurement(value, unit, endpoint, species="Human", matrix="safety assay"):
    return SimpleNamespace(
        id=1, endpoint_id=1, mean_value=None, value=value, qualifier="=", unit=unit,
        species=species, matrix=matrix, method="experimental safety", notes="", provenance_json={},
    ), endpoint


def test_safety_registry_models_and_explicit_optional_unavailable(db):
    project_id, _ = compound(db)
    models = {row["endpoint"]: row for row in list_admet(project_id, db)["models"]}
    for endpoint in SAFETY:
        row = models[endpoint]
        assert row["active"] and row["output_unit"] == "probability"
        assert row["details"]["endpoint_definition"] and row["details"]["training_dataset"]
        assert row["details"]["validation"] and row["details"]["license"]
    for endpoint, details in SAFETY_UNAVAILABLE.items():
        row = models[endpoint]
        assert not row["active"] and row["status"] == "MODEL_UNAVAILABLE"
        assert row["unavailable_reason"] == details["reason"]


@pytest.mark.parametrize("endpoint,labels", [
    ("hERG liability", {"BLOCKER", "NON_BLOCKER"}),
    ("Ames mutagenicity", {"MUTAGENIC", "NON_MUTAGENIC"}),
    ("DILI clinical liability", {"DILI_CONCERN", "NO_DILI_CONCERN"}),
])
def test_safety_predictions_are_probability_classifications_not_quantitative(endpoint, labels):
    result = predict_endpoint("CC(=O)Oc1ccccc1C(=O)O", endpoint)
    assert result["status"] == "COMPLETE" and 0 <= result["probability"] <= 1
    assert result["classification"] in labels and result["unit"] == "probability"
    assert result["confidence"] == "LOW" and result["safety_endpoint"]
    assert "IC50" not in result and "clinical_probability" not in result
    if endpoint != "hERG liability":
        assert len(result["ensemble_probabilities"]) == 5 and result["uncertainty"] >= 0


def test_herg_endpoint_assay_limit_and_independent_metrics_are_not_hidden():
    spec = MODEL_SPECS["hERG liability"]
    assert "neither a pure binding endpoint nor a pure functional" in spec["endpoint_definition"]
    independent = spec["independent_validation"]
    assert independent["n"] == 728 and independent["both_classes"]
    assert independent["balanced_accuracy"] == pytest.approx(0.5442154170)
    assert independent["specificity"] == pytest.approx(0.1129707113)
    assert predict_endpoint("CN(CCOc1ccc(NS(C)(=O)=O)cc1)CCc1ccc(NS(C)(=O)=O)cc1", "hERG liability")["confidence"] == "LOW"


def test_experimental_safety_classification_and_quantitative_separation():
    binary, name = measurement(1, "class", "hERG blocker classification")
    assert comparable_experimental("hERG liability", binary, name)[0] == 1
    evidence = classification_experimental_evidence("hERG liability", 0.9, [binary], {1: name})
    assert evidence[0]["comparison"] == "AGREES"
    ic50, name = measurement(2.5, "µM", "hERG inhibition IC50")
    assert comparable_experimental("hERG liability", ic50, name)[0] is None
    evidence = classification_experimental_evidence("hERG liability", 0.9, [ic50], {1: name})
    assert evidence[0]["comparison"] == "NOT_NUMERICALLY_COMPARABLE"
    assert evidence[0]["absolute_error"] is None and evidence[0]["relative_error_percent"] is None


def test_safety_applicability_domain_is_endpoint_specific_and_has_ood_case():
    for endpoint in SAFETY:
        available, reason = model_files_available(endpoint)
        assert available, reason
        index_key = MODEL_SPECS[endpoint].get("index_key")
        path = (admet_predictor.ADMET_AI_ROOT / "training" / index_key / "training.csv") if index_key else (admet_predictor.MODEL_ROOT / MODEL_SPECS[endpoint]["model_key"] / "training.csv")
        with path.open() as handle:
            row = next(csv.DictReader(handle))
            exact = row.get("smiles") or row.get("Smiles") or row.get("Drug")
        domain = applicability_domain(exact, endpoint)
        assert domain["nearest_training_similarity"] == 1 and domain["classification"] == "IN_DOMAIN"
    out = applicability_domain("C" * 200, "DILI clinical liability")
    assert out["classification"] == "OUT_OF_DOMAIN" and out["chemical_space_distance"] > 0


def test_safety_cache_experimental_precedence_provenance_and_project_isolation(db):
    project_id, version_id = compound(db, "Safety first")
    other_project, _ = compound(db, "Safety isolated", "SAFE-OTHER")
    create_admet_measurement(project_id, {
        "version_id": version_id, "endpoint": "Ames mutagenicity", "species": "Salmonella typhimurium",
        "matrix": "bacterial reverse mutation", "value": 0, "unit": "class",
        "method": "Ames assay", "source": "Public reference",
    }, db)
    first = run_admet_predictions(version_id, db)
    count = db.query(ADMETPrediction).count()
    second = run_admet_predictions(version_id, db)
    assert first["status"] == "COMPLETE" and second["status"] == "CACHED"
    assert db.query(ADMETPrediction).count() == count
    payload = list_admet(project_id, db)
    row = next(item for item in payload["predictions"] if item["endpoint"] == "Ames mutagenicity")
    assert row["preferred_result"]["source"] == "Experimental" and row["preferred_result"]["prediction_preserved"]
    assert set(row["provenance"]) >= {"record_type", "model_name", "model_version", "endpoint", "unit", "species", "dataset", "license", "validation", "applicability_domain", "confidence", "timestamp", "compound_version_id"}
    profile = payload["integrated_profiles"][str(version_id)]
    assert profile["overall_score"] is None and profile["experimental_precedence"]
    assert profile["provenance_audit"]["status"] == "PASS"
    assert list_admet(other_project, db)["predictions"] == []


def test_integrated_summary_is_deterministic_and_preserves_confidence():
    predictions = [{
        "id": 1, "version_id": 9, "endpoint": "hERG liability", "confidence": "LOW",
        "outputs": {"classification": "BLOCKER", "liability_summary": {"flag": "Potential hERG blocker liability"}},
        "experimental_comparisons": [], "provenance": {key: "x" for key in ("record_type", "model_name", "model_version", "endpoint", "unit", "species", "dataset", "license", "validation", "applicability_domain", "confidence", "timestamp", "compound_version_id")},
    }]
    profile = _integrated_admet_profile(9, predictions, [{"endpoint": "BSEP inhibitor", "active": False, "unavailable_reason": "no qualified model"}])
    assert any("LOW confidence" in item for item in profile["summary"]["concerns"])
    assert profile["summary"]["unknown"] == ["BSEP inhibitor: MODEL_UNAVAILABLE — no qualified model"]
    assert profile["overall_score"] is None


def test_acceptance_dataset_and_validation_artifact_are_honest():
    rows = list(csv.DictReader((ROOT / "validation/stage3_acceptance_dataset.csv").open()))
    assert len(rows) == 6 and {row["endpoint"] for row in rows} == set(SAFETY)
    assert all(row["source_url"].startswith("https://pubchem.ncbi.nlm.nih.gov/") for row in rows)
    assert all("INDEPENDENT" not in row["independent_status"] for row in rows)
    artifact = json.loads((ROOT / "models/admetica/validation/safety/independent_validation.json").read_text())
    assert artifact["overlap_policy"]["exact_canonical_smiles_removed"] == 7249
    assert artifact["metrics"]["specificity"] < 0.2
    assert "LOW confidence" in artifact["interpretation"]


def test_safety_ui_integrated_profile_details_and_no_score():
    source = (ROOT / "frontend/static/app.js").read_text()
    assert "function safetyPredictionTable" in source
    assert "Endpoint','Prediction','Probability','Experimental','Domain','Confidence','Model" in source
    assert "Stage 3 Integrated ADMET Profile" in source and "Experimental values take display precedence" in source
    assert "unavailableSafetyModels" in source and "MODEL_UNAVAILABLE" in source
    assert "No overall ADMET score or candidate ranking is calculated" in source
