import pytest
from fastapi import HTTPException
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend.main import (
    create_project,
    delete_project,
    get_project,
    update_project,
)
from backend.schemas import ProjectCreate, ProjectUpdate


@pytest.fixture()
def db():
    engine = create_engine("sqlite://", connect_args={"check_same_thread": False}, poolclass=StaticPool)
    from backend.database import Base
    Base.metadata.create_all(bind=engine)
    TestingSession = sessionmaker(bind=engine, autoflush=False, autocommit=False)
    database = TestingSession()
    try:
        yield database
    finally:
        database.close()


def test_project_crud(db):
    project = create_project(ProjectCreate(name="Test Target", target="Kinase X"), db)
    assert project.compound_count == 0
    updated = update_project(project.id, ProjectUpdate(mechanism_modality="inhibitor"), db)
    assert updated.mechanism_modality == "inhibitor"
    deleted = delete_project(project.id, {"confirmation_name": "Test Target"}, db)
    assert deleted["deleted_project_ids"] == [project.id]
    with pytest.raises(HTTPException):
        get_project(project.id, db)


def _add_compound(db, project_id, label, smiles):
    from backend.main import create_compound
    from backend.schemas import CompoundCreate
    return create_compound(project_id, CompoundCreate(compound_id=label, smiles=smiles), db)


def test_compound_validation_versioning_and_duplicates(db):
    project = create_project(ProjectCreate(name="Version Test"), db)
    first = _add_compound(db, project.id, "C001", "CC(=O)Oc1ccccc1C(=O)O")
    assert first["version"]["properties"]["molecular_weight"] > 170
    with pytest.raises(HTTPException) as duplicate:
        _add_compound(db, project.id, "C002", "OC(=O)c1ccccc1OC(C)=O")
    assert duplicate.value.status_code == 409
    with pytest.raises(HTTPException) as invalid:
        _add_compound(db, project.id, "C003", "invalid")
    assert invalid.value.status_code == 400
    from backend.main import update_compound, get_compound
    from backend.schemas import CompoundUpdate
    updated = update_compound(first["row_id"], CompoundUpdate(smiles="CC(C)Cc1ccc(cc1)C(C)C(O)=O"), db)
    assert updated["current_version"] >= 2
    detail = get_compound(first["row_id"], include_versions=True, db=db)
    assert [version["version_number"] for version in detail["versions"]] == [1, 2]
    assert detail["prediction_history"]
    assert detail["version"]["provenance"]["engine"] == "RDKit"


def test_compare_requires_two_project_compounds(db):
    project = create_project(ProjectCreate(name="Compare Test"), db)
    ids=[]
    for label,smiles in [("C001","CCO"),("C002","CC(=O)Oc1ccccc1C(=O)O")]:
        row=_add_compound(db,project.id,label,smiles);ids.append(row["row_id"])
    from backend.main import compare
    from fastapi.params import Query
    with pytest.raises(HTTPException):
        compare(project.id, ids=str(ids[0]), db=db)
    result=compare(project.id, ids=",".join(map(str,ids)), db=db)
    assert len(result["compounds"])==2
    assert set(result["metrics"]) >= {"MW","cLogP","TPSA","HBD","HBA","RotB","Fsp3","QED"}
