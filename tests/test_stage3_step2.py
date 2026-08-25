import csv
import io

import pytest
from fastapi import HTTPException
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend.database import Base
from backend.admet import ensure_admet_schema
from backend.main import (
    admet_export,
    admet_import,
    admet_import_preview,
    create_admet_measurement,
    create_compound,
    create_project,
    list_admet,
    run_admet_predictions,
)
from backend.schemas import CompoundCreate, ProjectCreate


@pytest.fixture()
def db():
    engine = create_engine(
        "sqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    Base.metadata.create_all(bind=engine)
    ensure_admet_schema(engine)
    database = sessionmaker(bind=engine, autoflush=False, autocommit=False)()
    try:
        yield database
    finally:
        database.close()


def create_project_compound(db, project_name="ADMET Step 2", compound_id="C001", smiles="CCO"):
    project = create_project(ProjectCreate(name=project_name), db)
    compound = create_compound(
        project.id,
        CompoundCreate(compound_id=compound_id, smiles=smiles),
        db,
    )
    return project.id, compound


def test_experimental_admet_save_and_read(db):
    project_id, compound = create_project_compound(db)
    version_id = compound["version"]["id"]
    payload = {
        "version_id": version_id,
        "endpoint": "Solubility",
        "species": "human",
        "matrix": "plasma",
        "value": "12.5",
        "unit": "µM",
        "qualifier": "=",
        "replicate": "R2",
        "mean": "12.0",
        "sd": "0.5",
        "n": "3",
        "method": "shake flask",
        "source": "Study A",
        "date": "2026-08-25",
        "notes": "verified experimental record",
    }

    saved = create_admet_measurement(project_id, payload, db)
    assert {key: saved[key] for key in ("value", "mean", "sd", "n")} == {
        "value": 12.5,
        "mean": 12.0,
        "sd": 0.5,
        "n": 3,
    }
    assert saved["type"] == "Experimental"
    assert saved["version_id"] == version_id

    listing = list_admet(project_id, db)
    assert listing["labels_by_version"][str(version_id)] == ("C001", 1)
    assert len(listing["measurements"]) == 1
    assert listing["measurements"][0]["notes"] == "verified experimental record"
    assert listing["endpoints"][0]["name"] == "Solubility"


def test_experimental_admet_rejects_cross_project_version(db):
    first_project_id, _ = create_project_compound(db, "First", "C001", "CCO")
    _, second_compound = create_project_compound(db, "Second", "C002", "CCN")
    with pytest.raises(HTTPException) as error:
        create_admet_measurement(
            first_project_id,
            {
                "version_id": second_compound["version"]["id"],
                "endpoint": "Solubility",
                "value": 1,
                "unit": "µM",
            },
            db,
        )
    assert error.value.status_code == 404


def test_admet_csv_preview_import_and_export(db):
    project_id, _ = create_project_compound(db)
    csv_text = (
        "compound_id,version_number,endpoint,species,matrix,value,unit,qualifier,replicate,mean,sd,n,method,source,date,notes\n"
        "C001,1,Permeability,human,Caco-2,8.4,10^-6 cm/s,=,R1,,,,Caco-2 assay,Study B,2026-08-24,imported row\n"
    )

    preview = admet_import_preview(project_id, {"csv": csv_text}, db)
    assert preview["valid_count"] == 1
    assert preview["errors"] == []
    assert preview["rows"][0]["row"] == 2

    imported = admet_import(project_id, {"csv": csv_text}, db)
    assert imported["imported"] == 1

    exported = admet_export(project_id, db)
    assert exported.media_type == "text/csv"
    rows = list(csv.DictReader(io.StringIO(exported.body.decode())))
    assert len(rows) == 1
    assert rows[0]["compound_id"] == "C001"
    assert rows[0]["version_number"] == "1"
    assert rows[0]["endpoint"] == "Permeability"
    assert rows[0]["value"] == "8.4"
    assert rows[0]["notes"] == "imported row"


def test_admet_csv_preview_reports_row_errors_without_importing(db):
    project_id, _ = create_project_compound(db)
    csv_text = (
        "compound_id,version_number,endpoint,value,unit\n"
        "C001,1,Solubility,not-a-number,µM\n"
        "UNKNOWN,1,Clearance,4.2,mL/min/kg\n"
    )

    preview = admet_import_preview(project_id, {"csv": csv_text}, db)
    assert preview["valid_count"] == 0
    assert [error["row"] for error in preview["errors"]] == [2, 3]

    with pytest.raises(HTTPException) as error:
        admet_import(project_id, {"csv": csv_text}, db)
    assert error.value.status_code == 400
    assert list_admet(project_id, db)["measurements"] == []


def test_admet_prediction_placeholder_is_auditable_and_project_scoped(db):
    first_project_id, first_compound = create_project_compound(db, "First ADMET", "C001", "CCO")
    first_version_id = first_compound["version"]["id"]
    first_run = run_admet_predictions(first_version_id, db)
    assert first_run["status"] == "NOT_INSTALLED"
    assert first_run["predictions"] == []
    assert first_run["models_available"] == 0
    registry = list_admet(first_project_id, db)["models"]
    assert len(registry) == 4
    assert all(model["status"] == "NOT_INSTALLED" and not model["active"] for model in registry)

    second_project_id, second_compound = create_project_compound(db, "Second ADMET", "C002", "CCN")
    second_version_id = second_compound["version"]["id"]
    run_admet_predictions(second_version_id, db)

    first_listing = list_admet(first_project_id, db)
    second_listing = list_admet(second_project_id, db)
    assert [run["version_id"] for run in first_listing["prediction_runs"]] == [first_version_id]
    assert [run["version_id"] for run in second_listing["prediction_runs"]] == [second_version_id]
