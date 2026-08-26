import pytest
from fastapi import HTTPException
from sqlalchemy import create_engine, delete
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

import backend.main as main_module
from backend.activity_models import ActivityMeasurement, ActivityPrediction, AssayDefinition
from backend.admet import (
    ADMETEndpoint, ADMETMeasurement, ADMETModelRegistry, ADMETPrediction, ADMETPredictionRun,
)
from backend.database import Base
from backend.main import (
    add_measurement, bulk_delete_projects, create_admet_measurement, create_assay,
    create_compound, create_project, delete_project,
)
from backend.metabolism import (
    ExperimentalMetabolite, MetabolicPredictionRun, MetabolicSoftSpot, PredictedMetabolite,
)
from backend.models import Compound, CompoundVersion, PredictionRun, Project, PropertyCalculation
from backend.optimization import OptimizationRun
from backend.proposal import (
    CandidatePredictionSnapshot, CandidateRanking, CandidateRejectionReason,
    CandidateTransformation, OptimizationCandidate, OptimizationProposalRun,
)
from backend.schemas import CompoundCreate, ProjectCreate


@pytest.fixture()
def db():
    engine = create_engine("sqlite://", connect_args={"check_same_thread": False}, poolclass=StaticPool)
    Base.metadata.create_all(bind=engine)
    session = sessionmaker(bind=engine, autoflush=False, autocommit=False)()
    try:
        yield session
    finally:
        session.close()


def populated_project(db, name):
    project = create_project(ProjectCreate(name=name, target="Deletion fixture"), db)
    compound = create_compound(
        project.id, CompoundCreate(compound_id=f"{name}-C1", name=f"{name} compound", smiles="CCO"), db,
    )
    version_id = compound["version"]["id"]
    assay = create_assay(project.id, {
        "assay_uid": "AS-" + name.upper().replace(" ", "-"), "name": f"{name} IC50",
        "measurement_type": "IC50", "unit": "nM",
    }, db)
    add_measurement(assay["id"], {"version_id": version_id, "value": 12, "unit": "nM"}, db)
    create_admet_measurement(project.id, {
        "version_id": version_id, "endpoint": "Solubility", "value": 8.5, "unit": "µM",
        "source": f"{name} experimental",
    }, db)
    endpoint = db.query(ADMETEndpoint).filter_by(project_id=project.id, name="Solubility").one()
    model = ADMETModelRegistry(
        endpoint_name=f"{name} delete model", model_name="Delete fixture", model_version="1",
        implementation_status="READY", output_unit="log mol/L", is_active=True,
    )
    db.add(model); db.flush()
    admet_run = ADMETPredictionRun(version_id=version_id, inputs_hash=f"{name}-admet", status="COMPLETE")
    db.add(admet_run); db.flush()
    db.add(ADMETPrediction(
        run_id=admet_run.id, endpoint_id=endpoint.id, version_id=version_id, model_id=model.id,
        predicted_value=-2.1, unit="log mol/L", confidence="LOW", applicability_domain="IN_DOMAIN",
    ))
    metabolic_run = MetabolicPredictionRun(
        version_id=version_id, inputs_hash=f"{name}-metabolism", engine_name="fixture",
        engine_version="1", status="COMPLETE",
    )
    db.add(metabolic_run); db.flush()
    spot = MetabolicSoftSpot(
        run_id=metabolic_run.id, version_id=version_id, rank=1, atom_index=0,
        transformation="aliphatic oxidation", phase="I",
    )
    db.add(spot); db.flush()
    db.add(PredictedMetabolite(
        run_id=metabolic_run.id, soft_spot_id=spot.id, version_id=version_id,
        canonical_smiles="CCO", isomeric_smiles="CCO", transformation="oxidation",
        source_atom=0, phase="I", rank=1,
    ))
    db.add(ExperimentalMetabolite(version_id=version_id, transformation="oxidation", source=name))
    optimization = OptimizationRun(project_id=project.id, parent_version_id=version_id, assay_id=assay["id"], status="COMPLETED")
    db.add(optimization); db.flush()
    proposal = OptimizationProposalRun(
        project_id=project.id, optimization_run_id=optimization.id, parent_version_id=version_id,
        status="COMPLETED",
    )
    db.add(proposal); db.flush()
    candidate = OptimizationCandidate(
        proposal_run_id=proposal.id, project_id=project.id, optimization_run_id=optimization.id,
        parent_version_id=version_id, candidate_number=1, canonical_smiles="CCN",
        isomeric_smiles="CCN", inchikey=f"{name[:8]}-KEY", status="ACCEPTED",
    )
    db.add(candidate); db.flush()
    db.add_all([
        CandidateTransformation(candidate_id=candidate.id, name="fixture", transformation_id="fixture-1", sequence_number=1),
        CandidatePredictionSnapshot(candidate_id=candidate.id, stage="ADMET", endpoint="Solubility"),
        CandidateRanking(candidate_id=candidate.id, rank=1, score=0.8, pareto_front=1),
        CandidateRejectionReason(candidate_id=candidate.id, code="INFO", detail="fixture", stage="FILTERING", hard_constraint=False),
    ])
    db.commit()
    return project, version_id


def assert_project_tree_absent(db, project_id, version_id):
    assert db.get(Project, project_id) is None
    assert db.query(Compound).filter_by(project_id=project_id).count() == 0
    assert db.query(CompoundVersion).filter_by(id=version_id).count() == 0
    assert db.query(AssayDefinition).filter_by(project_id=project_id).count() == 0
    assert db.query(ActivityMeasurement).filter_by(version_id=version_id).count() == 0
    assert db.query(ADMETMeasurement).filter_by(version_id=version_id).count() == 0
    assert db.query(ADMETPrediction).filter_by(version_id=version_id).count() == 0
    assert db.query(ADMETPredictionRun).filter_by(version_id=version_id).count() == 0
    assert db.query(PredictionRun).filter_by(version_id=version_id).count() == 0
    assert db.query(PropertyCalculation).filter_by(version_id=version_id).count() == 0
    assert db.query(MetabolicPredictionRun).filter_by(version_id=version_id).count() == 0
    assert db.query(ExperimentalMetabolite).filter_by(version_id=version_id).count() == 0
    assert db.query(OptimizationRun).filter_by(project_id=project_id).count() == 0
    assert db.query(OptimizationProposalRun).filter_by(project_id=project_id).count() == 0
    assert db.query(OptimizationCandidate).filter_by(project_id=project_id).count() == 0


def test_delete_empty_project_requires_exact_name(db):
    project = create_project(ProjectCreate(name="Empty Delete", target="T"), db)
    result = delete_project(project.id, {"confirmation_name": "Empty Delete"}, db)
    assert result["deleted_project_names"] == ["Empty Delete"]
    assert db.get(Project, project.id) is None


def test_delete_project_with_complete_data_tree(db):
    project, version_id = populated_project(db, "Full Delete")
    result = delete_project(project.id, {"confirmation_name": project.name}, db)
    assert result["deleted_project_ids"] == [project.id]
    assert_project_tree_absent(db, project.id, version_id)


def test_project_delete_isolation_preserves_other_project(db):
    doomed, doomed_version = populated_project(db, "Isolated Delete")
    preserved, preserved_version = populated_project(db, "Preserved Project")
    delete_project(doomed.id, {"confirmation_name": doomed.name}, db)
    assert_project_tree_absent(db, doomed.id, doomed_version)
    assert db.get(Project, preserved.id).name == "Preserved Project"
    assert db.query(ActivityMeasurement).filter_by(version_id=preserved_version).count() == 1
    assert db.query(ADMETMeasurement).filter_by(version_id=preserved_version).count() == 1
    assert db.query(OptimizationCandidate).filter_by(project_id=preserved.id).count() == 1


def test_wrong_confirmation_name_deletes_nothing(db):
    project, version_id = populated_project(db, "Exact Name Required")
    with pytest.raises(HTTPException) as error:
        delete_project(project.id, {"confirmation_name": "exact name required"}, db)
    assert error.value.status_code == 400
    assert db.get(Project, project.id) is not None
    assert db.query(CompoundVersion).filter_by(id=version_id).count() == 1


def test_project_delete_rolls_back_entire_transaction_on_failure(db, monkeypatch):
    project, version_id = populated_project(db, "Rollback Project")

    def injected_failure(session, project_ids):
        session.execute(delete(ActivityMeasurement).where(ActivityMeasurement.version_id == version_id))
        session.execute(delete(ADMETMeasurement).where(ADMETMeasurement.version_id == version_id))
        raise RuntimeError("injected deletion failure")

    monkeypatch.setattr(main_module, "_delete_project_tree_rows", injected_failure)
    with pytest.raises(HTTPException) as error:
        delete_project(project.id, {"confirmation_name": project.name}, db)
    assert error.value.status_code == 500 and "rolled back" in error.value.detail
    assert db.get(Project, project.id) is not None
    assert db.query(ActivityMeasurement).filter_by(version_id=version_id).count() == 1
    assert db.query(ADMETMeasurement).filter_by(version_id=version_id).count() == 1


def test_bulk_delete_requires_each_exact_name_and_preserves_unselected(db):
    first = create_project(ProjectCreate(name="Bulk One"), db)
    second = create_project(ProjectCreate(name="Bulk Two"), db)
    preserved = create_project(ProjectCreate(name="Bulk Preserved"), db)
    with pytest.raises(HTTPException):
        bulk_delete_projects({"projects": [
            {"id": first.id, "confirmation_name": first.name},
            {"id": second.id, "confirmation_name": "wrong"},
        ]}, db)
    assert db.get(Project, first.id) is not None and db.get(Project, second.id) is not None
    result = bulk_delete_projects({"projects": [
        {"id": first.id, "confirmation_name": first.name},
        {"id": second.id, "confirmation_name": second.name},
    ]}, db)
    assert set(result["deleted_project_ids"]) == {first.id, second.id}
    assert db.get(Project, preserved.id) is not None
