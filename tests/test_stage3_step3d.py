import json
from pathlib import Path

import pytest
from rdkit import Chem
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend.admet import ensure_admet_schema
from backend.database import Base
from backend.main import (create_admet_measurement, create_compound,
                          create_experimental_metabolite, create_project,
                          list_metabolism, run_admet_predictions,
                          run_metabolism_predictions, update_compound)
from backend.metabolic_soft_spot import (_sanitized_fragments,
                                         predict_soft_spots)
from backend.metabolism import (ExperimentalMetabolite,
                                MetabolicPredictionRun,
                                PredictedMetabolite)
from backend.schemas import CompoundCreate, CompoundUpdate, ProjectCreate


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


def compound(db, project_name="Stage 3D", compound_id="C001", smiles="COc1ccccc1"):
    project = create_project(ProjectCreate(name=project_name), db)
    created = create_compound(project.id, CompoundCreate(compound_id=compound_id, smiles=smiles), db)
    return project.id, created["row_id"], created["version"]["id"]


@pytest.mark.parametrize(("smiles", "transformation", "expected_metabolite"), [
    ("CC(=O)Nc1ccccc1", "Aromatic hydroxylation", "CC(=O)Nc1ccc(O)cc1"),
    ("Cc1ccccc1", "Benzylic oxidation", "OCc1ccccc1"),
    ("CCN(CC)CC", "N-dealkylation", "CCNCC"),
    ("COc1ccccc1", "O-dealkylation", "Oc1ccccc1"),
    ("CCOC(=O)c1ccccc1", "Ester hydrolysis", "O=C(O)c1ccccc1"),
])
def test_required_transformations_generate_sanitized_atom_mapped_products(smiles, transformation, expected_metabolite):
    result = predict_soft_spots(smiles)
    spot = next(row for row in result["spots"] if row["transformation"] == transformation)
    products = [row for row in result["metabolites"] if row["transformation"] == transformation]
    assert spot["rank"] <= 3 and 0 <= spot["atom_index"] < Chem.MolFromSmiles(smiles).GetNumAtoms()
    assert expected_metabolite in {row["canonical_smiles"] for row in products}
    assert all(Chem.MolFromSmiles(row["isomeric_smiles"]) is not None for row in products)


def test_phase_two_is_separate_and_has_no_fake_atom_probability():
    result = predict_soft_spots("Oc1ccccc1")
    phase_two = {row["transformation"]: row for row in result["spots"] if row["phase"] == "Phase II"}
    assert {"Glucuronidation", "Sulfation"}.issubset(phase_two)
    assert all(row["score_type"].endswith("not an atom probability") for row in phase_two.values())
    assert all(row["confidence"] == "LOW" for row in phase_two.values())


def test_invalid_transformation_product_is_rejected():
    invalid = Chem.MolFromSmiles("[CH5]", sanitize=False)
    invalid.GetAtomWithIdx(0).SetIntProp("react_atom_idx", 0)
    assert _sanitized_fragments(invalid, parent_heavy_atoms=1, source_atom=0) == []


def test_duplicate_metabolites_are_removed_across_symmetric_sites():
    result = predict_soft_spots("CCN(CC)CC")
    canonical = [row["canonical_smiles"] for row in result["metabolites"]]
    assert len(canonical) == len(set(canonical))
    assert canonical.count("CCNCC") == 1


def test_top_three_structure_highlight_and_clickable_ui_contract():
    result = predict_soft_spots("COc1ccccc1")
    svg = result["highlighted_svg"]
    assert "<svg" in svg and "Rank 1" in svg and "Rank 2" in svg and "Rank 3" in svg
    app = (Path(__file__).parents[1] / "frontend/static/app.js").read_text()
    assert "setSelectedSpotId(spot.id)" in app
    assert "PREDICTED METABOLITE HYPOTHESIS" in app


def test_prediction_persistence_cache_and_provenance(db):
    project_id, _, version_id = compound(db)
    first = run_metabolism_predictions(version_id, db)
    second = run_metabolism_predictions(version_id, db)
    assert first["status"] == "COMPLETE" and not first["cache_hit"]
    assert second["status"] == "CACHED" and second["cache_hit"]
    assert db.query(MetabolicPredictionRun).count() == 1
    assert db.query(PredictedMetabolite).count() == len(first["run"]["predicted_metabolites"])
    spot = first["run"]["spots"][0]
    assert spot["provenance"]["compound_version_id"] == version_id
    assert spot["provenance"]["engine"].startswith("SyGMa")
    assert spot["cyp_isoform"] == "CYP isoform not assigned"
    assert spot["model_evidence"]["status"] == "MODEL_UNAVAILABLE"


def test_experimental_metabolite_is_separate_and_validated(db):
    project_id, _, version_id = compound(db)
    row = create_experimental_metabolite(project_id, {
        "version_id": version_id, "smiles": "Oc1ccccc1", "transformation": "O-dealkylation",
        "observed_mass": 94.04, "mass_unit": "Da", "source": "Study A",
        "experiment": "LC-MS/MS", "notes": "Confirmed by standard",
    }, db)
    assert row["type"] == "Experimental" and row["label"] == "EXPERIMENTAL METABOLITE"
    assert row["canonical_smiles"] == "Oc1ccccc1"
    assert db.query(ExperimentalMetabolite).count() == 1
    with pytest.raises(Exception) as error:
        create_experimental_metabolite(project_id, {
            "version_id": version_id, "smiles": "[CH5]", "transformation": "Invalid",
        }, db)
    assert getattr(error.value, "status_code", None) == 400


def test_compound_version_and_project_isolation(db):
    first_project, row_id, first_version = compound(db, "First")
    first = run_metabolism_predictions(first_version, db)
    updated = update_compound(row_id, CompoundUpdate(smiles="Cc1ccccc1", change_note="new version"), db)
    second_version = updated["version"]["id"]
    second = run_metabolism_predictions(second_version, db)
    other_project, _, other_version = compound(db, "Other")
    assert first["run"]["version_id"] != second["run"]["version_id"]
    assert {run["version_id"] for run in list_metabolism(first_project, db)["runs"]} == {first_version, second_version}
    assert list_metabolism(other_project, db)["runs"] == []
    with pytest.raises(Exception) as error:
        create_experimental_metabolite(other_project, {
            "version_id": first_version, "transformation": "O-dealkylation",
        }, db)
    assert getattr(error.value, "status_code", None) == 404 and other_version != first_version


def test_hlm_experimental_priority_and_cyp_linkage_preserve_attribution_limits(db):
    project_id, _, version_id = compound(db, smiles="CCN(CC)Cc1ccccc1")
    create_admet_measurement(project_id, {
        "version_id": version_id, "endpoint": "HLM intrinsic clearance", "species": "Human",
        "matrix": "Human liver microsomes", "value": 85, "unit": "µL/min/mg protein",
    }, db)
    run_admet_predictions(version_id, db)
    result = run_metabolism_predictions(version_id, db)["run"]
    summary = result["liability_summary"]
    assert summary["microsomal_evidence"][0]["source"] == "Experimental"
    assert any(item["endpoint"].endswith("substrate") for item in summary["cyp_evidence"])
    assert "does not assign this atom" in summary["cyp_attribution_limit"]
    assert all(spot["cyp_isoform"] == "CYP isoform not assigned" for spot in result["spots"])


def test_known_drug_sanity_artifact_separates_publisher_and_local_validation():
    path = Path(__file__).parents[1] / "models/sygma/validation/known_drug_sanity.json"
    result = json.loads(path.read_text())
    assert result["metrics"] == {
        "n": 5, "top_1_accuracy": 0.8, "top_2_accuracy": 1.0,
        "top_3_accuracy": 1.0, "atom_level_recall": 1.0,
    }
    assert result["training_overlap_audit"]["status"] == "NOT_ASSESSABLE"
    assert result["publisher_reported_validation"]
    assert all(row["known_metabolite_generated"] for row in result["references"])
