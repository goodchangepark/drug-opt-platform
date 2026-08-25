import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend import admet_predictor
from backend.admet import ADMETPrediction, ensure_admet_schema
from backend.admet_predictor import (
    MODEL_SPECS, applicability_domain, comparable_experimental,
    metabolic_stability_assessment, model_files_available, predict_endpoint,
)
from backend.database import Base
from backend.main import (add_measurement, compare, create_admet_measurement, create_assay,
                          create_compound, create_project, list_admet, run_admet_predictions)
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


def compound(db, smiles="CC(=O)Oc1ccccc1C(=O)O"):
    project = create_project(ProjectCreate(name="Stage 3B"), db)
    created = create_compound(project.id, CompoundCreate(compound_id="C001", smiles=smiles), db)
    return project.id, created["version"]["id"]


def measurement(value, unit, endpoint, species="", matrix="", provenance=None):
    return SimpleNamespace(
        mean_value=None, value=value, qualifier="=", unit=unit, matrix=matrix,
        method="microsomal assay", species=species, provenance_json=provenance or {},
    ), endpoint


def test_stage3b_registry_is_species_specific_and_unavailable_species_are_explicit(db):
    project_id, _ = compound(db)
    models = {item["endpoint"]: item for item in list_admet(project_id, db)["models"]}
    assert models["Plasma protein binding"]["output_unit"] == "% bound"
    for code, species in (("HLM", "Human"), ("RLM", "Rat"), ("MLM", "Mouse")):
        row = models[f"{code} intrinsic clearance"]
        assert row["active"] and row["output_unit"] == "log10(mL/min/kg)"
        assert row["details"]["species"] == species
    assert models["Dog liver microsomal intrinsic clearance"]["status"] == "MODEL_UNAVAILABLE"
    assert models["Monkey liver microsomal intrinsic clearance"]["status"] == "MODEL_UNAVAILABLE"
    assert model_files_available("Dog liver microsomal intrinsic clearance")[0] is False


def test_human_ppb_prediction_preserves_bound_and_derives_fu():
    result = predict_endpoint("CC(=O)Oc1ccccc1C(=O)O", "Plasma protein binding")
    assert result["status"] == "COMPLETE" and result["unit"] == "% bound"
    derived = result["derived_outputs"]
    assert derived["fu_fraction"] == pytest.approx(1 - result["predicted_value"] / 100)
    assert derived["fu_percent"] == pytest.approx(100 - result["predicted_value"])
    assert "fu = 1 - fraction bound" in derived["derivation"]


@pytest.mark.parametrize("unit,value,expected", [
    ("% bound", 90, 90), ("fraction bound", 0.9, 90), ("fu", 0.1, 90), ("% unbound", 10, 90),
])
def test_ppb_unit_handling(unit, value, expected):
    row, name = measurement(value, unit, "Human PPB", species="Human", matrix="plasma")
    converted, _ = comparable_experimental("Plasma protein binding", row, name)
    assert converted == pytest.approx(expected)


def test_microsomal_species_isolation_and_units():
    hlm, name = measurement(20, "mL/min/kg", "HLM intrinsic clearance", species="Human", matrix="HLM")
    converted, _ = comparable_experimental("HLM intrinsic clearance", hlm, name)
    assert converted == pytest.approx(1.30103)
    assert comparable_experimental("RLM intrinsic clearance", hlm, name)[0] is None

    raw, name = measurement(10, "µL/min/mg", "HLM intrinsic clearance", species="Human", matrix="HLM")
    assert comparable_experimental("HLM intrinsic clearance", raw, name)[0] is None
    raw.provenance_json = {"microsomal_protein_mg_per_g_liver": 45, "liver_weight_g_per_kg": 26}
    scaled, note = comparable_experimental("HLM intrinsic clearance", raw, name)
    assert scaled == pytest.approx(1.068186, abs=1e-6)
    assert "45" in note and "26" in note


def test_hlm_rlm_mlm_are_distinct_tasks_and_no_half_life_is_fabricated():
    results = {endpoint: predict_endpoint("CCN(CC)CCCC(C)NC1=C2C=CC(Cl)=CC2=NC=C1", endpoint)
               for endpoint in ("HLM intrinsic clearance", "RLM intrinsic clearance", "MLM intrinsic clearance")}
    assert len({round(result["predicted_value"], 5) for result in results.values()}) == 3
    assert all(result["unit"] == "log10(mL/min/kg)" for result in results.values())
    assert all("half" not in json.dumps(result).lower() and "t1/2" not in json.dumps(result).lower() for result in results.values())


def test_experimental_priority_comparison_cache_and_liability_evidence(db):
    project_id, version_id = compound(db)
    create_admet_measurement(project_id, {
        "version_id": version_id, "endpoint": "Human PPB", "species": "Human", "matrix": "plasma",
        "value": 90, "unit": "% bound",
    }, db)
    create_admet_measurement(project_id, {
        "version_id": version_id, "endpoint": "HLM intrinsic clearance", "species": "Human", "matrix": "HLM",
        "value": 100, "unit": "mL/min/kg",
    }, db)
    first = run_admet_predictions(version_id, db)
    assert first["status"] == "COMPLETE" and not first["cache_hit"]
    second = run_admet_predictions(version_id, db)
    assert second["status"] == "CACHED" and second["cache_hit"]
    assert db.query(ADMETPrediction).count() == len(first["predictions"])
    listing = list_admet(project_id, db)["predictions"]
    ppb = next(row for row in listing if row["endpoint"] == "Plasma protein binding")
    hlm = next(row for row in listing if row["endpoint"] == "HLM intrinsic clearance")
    assert ppb["preferred_result"]["source"] == "Experimental"
    assert ppb["predicted_value"] is not None and ppb["experimental_comparisons"][0]["absolute_error"] >= 0
    assert hlm["preferred_result"]["source"] == "Experimental"
    assert hlm["outputs"]["experimental_metabolic_stability_assessment"]["metabolic_liability_flag"] == "METABOLIC STABILITY CONCERN"


def test_endpoint_specific_applicability_domain_and_out_of_domain():
    with (admet_predictor.OPENADMET_ROOT / "X_train.csv").open() as stream:
        stream.readline(); exact_smiles = stream.readline().strip()
    assert applicability_domain(exact_smiles, "HLM intrinsic clearance")["nearest_training_similarity"] == 1.0
    assert applicability_domain("C" * 200, "RLM intrinsic clearance")["classification"] == "OUT_OF_DOMAIN"


def test_metabolic_assessment_thresholds_are_auditable():
    stable = metabolic_stability_assessment("HLM intrinsic clearance", 0.5)
    unstable = metabolic_stability_assessment("HLM intrinsic clearance", 2.0)
    assert stable["category"] == "STABLE" and stable["metabolic_liability_flag"] is None
    assert unstable["category"] == "UNSTABLE"
    assert unstable["metabolic_liability_flag"] == "METABOLIC STABILITY CONCERN"
    assert "25th and 75th percentiles" in unstable["thresholds"]["basis"]


def test_independent_validation_artifact_records_weak_generalization_honestly():
    result = json.loads((Path(__file__).parents[1] / "models/openadmet/microsomal_clearance/independent_validation.json").read_text())
    assert result["clearance"]["HLM"]["n"] > 3000
    assert {"MAE", "RMSE", "R2", "Spearman"}.issubset(result["clearance"]["HLM"])
    assert result["human_ppb_percent_bound"]["n"] == 185
    assert [row["name"] for row in result["reference_compounds"]] == ["Rifampicin", "Isoniazid", "Ethionamide"]
    assert result["reference_directionality_spearman"] == 0.5


def test_compound_comparison_exposes_activity_and_stage3a_3b_without_ranking(db):
    project = create_project(ProjectCreate(name="Comparison"), db)
    first = create_compound(project.id, CompoundCreate(compound_id="C001", smiles="CCO"), db)
    second = create_compound(project.id, CompoundCreate(compound_id="C002", smiles="CCN"), db)
    assay = create_assay(project.id, {"name": "Primary IC50", "measurement_type": "IC50", "unit": "nM"}, db)
    add_measurement(assay["id"], {"version_id": first["version"]["id"], "value": 25, "unit": "nM"}, db)
    run_admet_predictions(first["version"]["id"], db)
    run_admet_predictions(second["version"]["id"], db)
    result = compare(project.id, f"{first['row_id']},{second['row_id']}", db)
    assert [name for name in ("Activity", "HLM", "RLM", "PPB", "Solubility", "Caco-2") if name in result["metrics"]] == [
        "Activity", "HLM", "RLM", "PPB", "Solubility", "Caco-2",
    ]
    assert next(row for row in result["compounds"] if row["compound"] == "C001")["Activity"] == 25
    assert "score" not in result and "ranking" not in result
