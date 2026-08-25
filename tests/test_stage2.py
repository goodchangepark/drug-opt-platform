import math

import numpy as np
import pytest
from fastapi import HTTPException
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend.database import Base
from backend.main import create_assay, create_compound, create_project, list_assays, predict_activity, sar_table, train_assay_model
from backend.main import add_measurement as add_activity
from backend.qsar import normalize_concentration, pactivity, value_from_pactivity
from backend.schemas import CompoundCreate, ProjectCreate


SMILES = [
    "CC(=O)Oc1ccccc1C(=O)O","CC(C)Cc1ccc(C(C)C(=O)O)cc1","Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "CC(=O)Nc1ccc(O)cc1","COc1ccc(C#C)cc1Nc2ncnc3cc(Cl)c(Nc4ccc(C#C)c(OC)c4)cc23",
    "CS(=O)(=O)CCNc1ncnc2cc(Cl)c(Nc3ccc(C#C)c(OC)c3)cc12",
]


@pytest.fixture()
def db():
    engine=create_engine("sqlite://",connect_args={"check_same_thread":False},poolclass=StaticPool)
    Base.metadata.create_all(engine); session=sessionmaker(bind=engine)(); yield session; session.close()


def make_project(db):
    project=create_project(ProjectCreate(name="Stage 2"),db)
    versions=[]
    for i,smi in enumerate(SMILES,1):
        compound=create_compound(project.id,CompoundCreate(compound_id=f"C{i:03d}",smiles=smi),db)
        versions.append(compound["row_id"])
    return project,versions


def test_units_and_pactivity():
    normalized,prov=normalize_concentration(1,"µM")
    assert normalized==1000 and prov["factor_to_nM"]==1000
    assert abs(pactivity(normalized)-6)<1e-10
    assert value_from_pactivity(9)==1


def test_qualifiers_replicates_and_summary(db):
    project,versions=make_project(db)
    assay_id=create_assay(project.id,{"name":"IC50 BaF3","measurement_type":"IC50","unit":"nM"},db)["id"]
    for raw,unit in [(10,"nM"),("0.015","µM")]:
        pass
    for value,unit in [(10,"nM"),(15,"nM"),(11,"nM")]:
        add_activity(assay_id,{"version_id":versions[0],"value":value,"unit":unit},db)
    sar=sar_table(project.id,assay_id,db)
    exp=sar["compounds"][0]["experimental"]
    assert exp["type"]=="Experimental" and exp["n"]==3 and exp["mean_nm"]==12 and exp["sd_nm"]>0


def test_qsar_policy_insufficient_data(db):
    project,versions=make_project(db)
    assay_id=create_assay(project.id,{"name":"Small IC50","measurement_type":"IC50"},db)["id"]
    result=train_assay_model(assay_id,db)
    assert result["policy"]["status"]=="INSUFFICIENT DATA" and result["model"] is None
    with pytest.raises(HTTPException):
            predict_activity(assay_id,versions[4],db)


def test_similarity_prediction_with_five_experiments(db):
    project,versions=make_project(db)
    assay_id=create_assay(project.id,{"name":"Similarity IC50","measurement_type":"IC50"},db)["id"]
    for version,p in zip(versions[:5],[7.5,7,8,6.5,7.2]):
        add_activity(assay_id,{"version_id":version,"value":round(10**(9-p),2),"unit":"nM"},db)
    prediction=predict_activity(assay_id,versions[5],db)
    assert prediction["type"]=="Predicted" and prediction["confidence"] in {"HIGH","MEDIUM","LOW"}
    assert prediction["nearest_neighbors"] and prediction["applicability_domain"]
