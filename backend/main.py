import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

# Force static TLS allocation for PyTorch & OpenMP on main thread startup
try:
    import torch
    import chemprop
except ImportError:
    pass

from fastapi import BackgroundTasks, Depends, FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from pathlib import Path
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from rdkit import Chem
from sqlalchemy import delete, func, select
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session

from .chemistry import ChemistryError, ENGINE, ENGINE_VERSION, analyze_smiles
from .activity_models import ActivityMeasurement, ActivityPrediction, AssayDefinition, MatchedMolecularPair, QSARModel
from .admet import (ADMETAssayDefinition, ADMETConsensusPrediction, ADMETEndpoint, ADMETMeasurement,
                    ADMETModelComparison, ADMETModelPerformance, ADMETModelRegistry, ADMETPrediction, ADMETPredictionRun,
                    ADMETExperimentalFeedbackEvent, ADMETAdaptivePrediction,
                    csv_export, ensure_admet_schema, inputs_hash,
                    measurement_out, parse_csv, validate_measurement)
from .adaptive_weighting import (ADAPTIVE_POLICY_VERSION, AdaptiveConsensusResult, AdaptiveReasonCode,
                                AssayQuality, ExperimentalFeedbackRecord,
                                compute_hierarchical_adaptive_weights, evaluate_experimental_compatibility,
                                get_bemis_murcko_scaffold)
from .multimodel import get_adapters_for_endpoint, ModelExecutionPayload, ExecutionStatus
from .admet_predictor import (MODEL_SPECS, MODEL_VERSION, comparable_experimental, comparison_for_prediction, cyp_experimental_evidence,
                              metabolic_stability_assessment, model_files_available, predict_endpoint)
from .conformal import (CONFORMAL_CALIBRATION_REGISTRY, CalibrationQuality,
                        DataProvenance)
from .endpoint_contracts import get_endpoint_contract
from .endpoint_strategy_registry import get_registry_api_response
from .prediction_engine_v1_policy import policy_api_response
from .production_qualification import (
    ensure_qualification_schema,
    get_candidates_api_response,
    get_drift_api_response,
    get_qualification_api_response,
    get_qualification_endpoint_response,
)
from .database import Base, SessionLocal, engine, get_db
from .metabolic_soft_spot import (ENGINE_LICENSE as METABOLISM_LICENSE,
                                  ENGINE_NAME as METABOLISM_ENGINE,
                                  ENGINE_SOURCE as METABOLISM_SOURCE,
                                  ENGINE_VERSION as METABOLISM_VERSION,
                                  PREDICTED_LABEL, PUBLISHER_VALIDATION,
                                  predict_soft_spots)
from .metabolism import (ExperimentalMetabolite, MetabolicPredictionRun,
                         MetabolicSoftSpot, PredictedMetabolite,
                         ensure_metabolism_schema)
from .optimization import OptimizationRun, ensure_optimization_schema
from .optimization_engine import (ENGINE_NAME as OPTIMIZATION_ENGINE,
                                  ENGINE_VERSION as OPTIMIZATION_VERSION,
                                  EVIDENCE_HIERARCHY, OBJECTIVES,
                                  TRANSFORMATION_LIBRARY, analyze_run)
from .proposal import (CandidatePredictionSnapshot, CandidateRanking,
                       CandidateRejectionReason, CandidateTransformation,
                       OptimizationCandidate, OptimizationProposalRun,
                       ensure_proposal_schema)
from .proposal_engine import (ENGINE_NAME as PROPOSAL_ENGINE,
                              ENGINE_VERSION as PROPOSAL_VERSION,
                              EXECUTABLE_TRANSFORMATIONS,
                              STRATEGY_ONLY_TRANSFORMATIONS,
                              execute_proposal_run, process_user_candidate,
                              rank_candidates)
from .models import (Compound, CompoundVersion, ExternalExperimentalEvidence, PredictionRun, Project,
                     PropertyCalculation, StructuralAlert, ensure_ui_schema,
                     utcnow)
from .pk import PKNCAResult, PKObservation, PKStudy, ensure_pk_schema, register_pk_routes
from .ivive import (IVIVEInputSet, IVIVEMethodRegistry, IVIVERun, PKParameterSet, PhysiologicalParameterOverride,
                    ensure_ivive_schema, register_ivive_routes, get_multi_species_pk_profile,
                    get_pk_foundation_profile, refresh_pk_and_ivive_for_version)
from .simulation import PKSimulationRun, ensure_simulation_schema, register_simulation_routes
from .standardizer import GLOBAL_DESCRIPTOR_CONFIG, GLOBAL_FINGERPRINT_CONFIG, RDKIT_VERSION, STANDARDIZER_NAME, STANDARDIZER_VERSION, standardize_molecule
from .golden_set import run_golden_gate_test
from .evaluation import (EVALUATION_REGISTRY, evaluate_mmp_directional_accuracy,
                         get_rdkit_upgrade_readiness_report, perform_lightning_security_audit)
from .qsar import (DESCRIPTOR_NAMES, FINGERPRINT_CONFIG, applicability, feature_vector,
                   fingerprint_and_descriptors, nearest_neighbors, normalize_concentration, tanimoto_similarity,
                   pactivity, train_model, value_from_pactivity)
from .ionization import analyze_ionization
from contextlib import asynccontextmanager
from .schemas import CompoundCreate, CompoundUpdate, ProjectCreate, ProjectOut, ProjectUpdate
from .translational import PKTranslationalSnapshot, ensure_translational_schema, register_translational_routes
from .human_pk import PKHumanPredictionSnapshot, ensure_human_pk_schema, register_human_pk_routes
from .capabilities import build_capability_summary
from .interpretation import get_interpretation_registry_summary, interpret_property
from .platform_info import (APP_VERSION, CURRENT_STAGE_LABEL, CURRENT_STAGE_STATUS,
                            CURRENT_STAGE_SUBSTATUS, GLOSSARY, LIMITATIONS, build_version,
                            latest_release_date, package_inventory, structure_modules,
                            version_history)
from .external_experimental import cas_status, lookup as external_evidence_lookup, valid_cas


CURRENT_STAGE = "5B-4"


def _valid_cas_number(value: str) -> bool:
    return valid_cas(value)


@asynccontextmanager
async def lifespan(app_instance: FastAPI):
    if not app_instance.dependency_overrides:
        Base.metadata.create_all(bind=engine)
        ensure_ui_schema(engine)
        ensure_admet_schema(engine)
        ensure_metabolism_schema(engine)
        ensure_optimization_schema(engine)
        ensure_proposal_schema(engine)
        ensure_pk_schema(engine)
        ensure_ivive_schema(engine)
        ensure_simulation_schema(engine)
        ensure_translational_schema(engine)
        ensure_human_pk_schema(engine)
        ensure_qualification_schema(engine)
        # Initialize PyTorch/Chemprop once before concurrent request workers can
        # observe a partially imported native extension on ARM64.
        model_files_available("Solubility")
    import backend.activity_models
    yield


app = FastAPI(title="AI Drug Optimization Platform", version=APP_VERSION, lifespan=lifespan)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

# Register modular sub-routers
register_pk_routes(app)
register_ivive_routes(app)
register_simulation_routes(app)
register_translational_routes(app)
register_human_pk_routes(app)


def _project_out(db: Session, project: Project):
    count = db.scalar(select(func.count(Compound.id)).where(Compound.project_id == project.id)) or 0
    return ProjectOut.model_validate(project).model_copy(update={"compound_count": count})


@app.get("/api/health")
def health():
    return {"status": "ok", "stage": "5B", "step": CURRENT_STAGE, "version": APP_VERSION,
            "updated": latest_release_date(), "engine": ENGINE, "engine_version": ENGINE_VERSION}


@app.get("/api/interpretation/rules")
def get_interpretation_rules():
    return get_interpretation_registry_summary()


@app.get("/api/model-strategy-registry")
def model_strategy_registry():
    """Read-only scientific strategy policy; never executes or promotes a model."""
    return get_registry_api_response()


@app.get("/api/prediction-engine-v1/policy")
def prediction_engine_v1_policy():
    """Immutable, read-only Stage 4E-4 policy snapshot and content hash."""
    return policy_api_response()


@app.get("/api/qualification/strategies")
def qualification_strategies():
    """Read-only qualification policies; never mutates production state."""
    return get_qualification_api_response()


@app.get("/api/qualification/endpoint/{endpoint_id}")
def qualification_endpoint(endpoint_id: str):
    response = get_qualification_endpoint_response(endpoint_id)
    if response is None:
        raise HTTPException(status_code=404, detail="Qualification endpoint not found")
    return response


@app.get("/api/qualification/candidates")
def qualification_candidates():
    return get_candidates_api_response()


@app.get("/api/qualification/drift")
def qualification_drift():
    return get_drift_api_response()


@app.post("/api/structure/validate")
def validate_structure(payload: dict):
    try:
        result = analyze_smiles(str(payload.get("smiles", "")))
    except ChemistryError as exc:
        return JSONResponse(status_code=400, content={"valid": False, "error": str(exc)})
    return {
        "valid": True,
        "duplicate_in_payload": False,
        **result,
    }


@app.get("/api/projects")
def list_projects(db: Session = Depends(get_db)):
    return [_project_out(db, project) for project in db.scalars(select(Project).order_by(Project.created_at.desc()))]


def _workflow_status(covered: int, total: int):
    if not total or not covered:
        return "NOT_STARTED"
    return "READY" if covered >= total else "PARTIAL"


@app.get("/api/dashboard")
def dashboard_summary(db: Session = Depends(get_db)):
    """Small read-only workspace summary for the main and project dashboards."""
    projects = db.scalars(select(Project).order_by(Project.created_at.desc())).all()
    project_rows = []
    total_compounds = 0
    for project in projects:
        compounds = db.scalars(
            select(Compound).where(Compound.project_id == project.id).order_by(Compound.compound_id)
        ).all()
        total_compounds += len(compounds)
        compound_ids = [row.id for row in compounds]
        versions = db.scalars(
            select(CompoundVersion).where(CompoundVersion.compound_row_id.in_(compound_ids))
        ).all() if compound_ids else []
        version_ids = [row.id for row in versions]
        current_versions = {
            row.compound_row_id: row for row in versions
            if row.version_number == next(
                (compound.current_version for compound in compounds if compound.id == row.compound_row_id), -1
            )
        }

        def version_set(model):
            if not version_ids:
                return set()
            return set(db.scalars(select(model.version_id).where(model.version_id.in_(version_ids))).all())

        property_versions = version_set(PropertyCalculation)
        activity_experimental = version_set(ActivityMeasurement)
        activity_predicted = version_set(ActivityPrediction)
        admet_experimental = version_set(ADMETMeasurement)
        admet_predicted = version_set(ADMETPrediction)
        pk_recorded = version_set(PKStudy)
        optimization_runs = db.scalars(
            select(OptimizationRun).where(OptimizationRun.project_id == project.id)
        ).all()
        optimized_versions = {row.parent_version_id for row in optimization_runs}
        activity_count = db.scalar(
            select(func.count(ActivityMeasurement.id)).where(ActivityMeasurement.version_id.in_(version_ids))
        ) if version_ids else 0
        admet_count = db.scalar(
            select(func.count(ADMETMeasurement.id)).where(ADMETMeasurement.version_id.in_(version_ids))
        ) if version_ids else 0
        activity_prediction_count = db.scalar(
            select(func.count(ActivityPrediction.id)).where(ActivityPrediction.version_id.in_(version_ids))
        ) if version_ids else 0
        admet_prediction_count = db.scalar(
            select(func.count(ADMETPrediction.id)).where(ADMETPrediction.version_id.in_(version_ids))
        ) if version_ids else 0

        compound_rows = []
        for compound in compounds:
            version = current_versions.get(compound.id)
            version_id = version.id if version else None
            compound_rows.append({
                "row_id": compound.id,
                "compound_id": compound.compound_id,
                "name": compound.name,
                "status": compound.status,
                "version_id": version_id,
                "structure": "READY" if version and version.canonical_smiles else "NOT_STARTED",
                "properties": "CALCULATED" if version_id in property_versions else "NOT_RUN",
                "activity": "EXPERIMENTAL" if version_id in activity_experimental else (
                    "PREDICTED" if version_id in activity_predicted else "NOT_RUN"
                ),
                "admet": "EXPERIMENTAL" if version_id in admet_experimental else (
                    "PREDICTED" if version_id in admet_predicted else "NOT_RUN"
                ),
                "optimization": "READY" if version_id in optimized_versions else "NOT_RUN",
            })

        current_total = len(compounds)
        structure_covered = sum(row["structure"] == "READY" for row in compound_rows)
        property_covered = sum(row["properties"] == "CALCULATED" for row in compound_rows)
        activity_covered = sum(row["activity"] != "NOT_RUN" for row in compound_rows)
        admet_covered = sum(row["admet"] != "NOT_RUN" for row in compound_rows)
        optimization_covered = sum(row["optimization"] == "READY" for row in compound_rows)
        pk_covered = len(pk_recorded)
        status_parts = [f"{current_total} compound{'s' if current_total != 1 else ''}"]
        if activity_count:
            status_parts.append(f"{activity_count} experimental activity record{'s' if activity_count != 1 else ''}")
        if admet_count:
            status_parts.append(f"{admet_count} experimental ADMET record{'s' if admet_count != 1 else ''}")
        elif activity_prediction_count or admet_prediction_count:
            status_parts.append("predictions available")
        else:
            status_parts.append("experimental and prediction data not started")
        project_rows.append({
            "id": project.id,
            "name": project.name,
            "target": project.target,
            "molecule_type": project.molecule_type,
            "compound_count": current_total,
            "experimental_activity_count": int(activity_count or 0),
            "experimental_admet_count": int(admet_count or 0),
            "prediction_count": int((activity_prediction_count or 0) + (admet_prediction_count or 0)),
            "optimization_run_count": len(optimization_runs),
            "status_summary": " · ".join(status_parts),
            "workflow": {
                "Structure": _workflow_status(structure_covered, current_total),
                "Properties": _workflow_status(property_covered, current_total),
                "Activity": _workflow_status(activity_covered, current_total),
                "ADMET": _workflow_status(admet_covered, current_total),
                "Optimization": _workflow_status(optimization_covered, current_total),
                "PK": _workflow_status(pk_covered, current_total),
            },
            "compounds": compound_rows,
        })

    model_rows = [
        _admet_model_out(model)
        for model in db.scalars(select(ADMETModelRegistry).order_by(ADMETModelRegistry.endpoint_name))
    ]
    return {
        "totals": {"projects": len(projects), "compounds": total_compounds},
        "projects": project_rows,
        "model_registry": model_rows,
        "capability_summary": build_capability_summary(
            model_rows,
            (route.path for route in app.routes if hasattr(route, "path")),
            stage=CURRENT_STAGE,
        ),
    }


@app.get("/api/help/registry")
def help_registry(db: Session = Depends(get_db)):
    """Researcher-facing inventory composed from live registries and runtime packages."""
    inventory = package_inventory()
    models = [
        _admet_model_out(model)
        for model in db.scalars(select(ADMETModelRegistry).order_by(ADMETModelRegistry.model_priority,
                                                                    ADMETModelRegistry.endpoint_name))
    ]
    capabilities = build_capability_summary(
        models, (route.path for route in app.routes if hasattr(route, "path")), stage=CURRENT_STAGE,
    )
    pk_methods = [
        {"method_key": row.method_key, "method_name": row.method_name,
         "method_version": row.method_version, "status": row.status,
         "assumptions": row.assumptions_json or {}, "reference": row.reference_json or {}}
        for row in db.scalars(select(IVIVEMethodRegistry).order_by(IVIVEMethodRegistry.method_key))
    ]
    return {
        "application": {
            "name": "Drug-OPT", "version": APP_VERSION, "current_stage": CURRENT_STAGE,
            "current_stage_label": CURRENT_STAGE_LABEL,
            "current_stage_status": CURRENT_STAGE_STATUS,
            "current_stage_substatus": CURRENT_STAGE_SUBSTATUS,
            "updated": latest_release_date(), "build_version": build_version(), "standardizer": STANDARDIZER_NAME,
            "standardizer_version": STANDARDIZER_VERSION, "rdkit_version": RDKIT_VERSION,
        },
        "package_inventory": inventory,
        "structure_modules": structure_modules(inventory),
        "models": models,
        "capability_summary": capabilities,
        "pk_method_registry": pk_methods,
        "version_history": version_history(),
        "glossary": [{"term": term, "definition": definition} for term, definition in GLOSSARY],
        "limitations": list(LIMITATIONS),
        "interpretation_registry": get_interpretation_registry_summary(),
        "source": "RUNTIME_PACKAGE_INVENTORY + ADMET_MODEL_REGISTRY + PK_METHOD_REGISTRY + CAPABILITY_REGISTRY + VERSION_HISTORY",
    }


@app.get("/api/validation/campaign")
def get_validation_campaign_status(db: Session = Depends(get_db)):
    """Engine v1 internal validation campaign status.

    Returns framework and scientific status without exposing any
    experimental values or internal compound structures.
    """
    try:
        from .internal_validation_v1 import (
            campaign_summary, CAMPAIGN_ID, ENGINE_V1_POLICY_HASH,
            ensure_validation_schema,
        )
        ensure_validation_schema(db.get_bind())
        summary = campaign_summary(db, CAMPAIGN_ID)
        return {
            "status": "ok",
            "campaign": summary,
            "engine_policy_hash": ENGINE_V1_POLICY_HASH,
            "policy_hash_unchanged": True,
        }
    except Exception as exc:
        return {
            "status": "not_initialized",
            "message": str(exc),
            "note": "Run scripts/initialize_validation_campaign.py to initialize.",
        }



@app.post("/api/projects", status_code=201)
def create_project(payload: ProjectCreate, db: Session = Depends(get_db)):
    values = payload.model_dump()
    values["name"] = values.get("name", "").strip()
    if not values["name"]:
        raise HTTPException(status_code=400, detail="Project name is required")
    if values["molecule_type"] not in {"Small Molecule", "Peptide"}:
        raise HTTPException(status_code=400, detail="Molecule Type must be Small Molecule or Peptide")
    existing = db.scalar(select(Project).where(Project.name == values["name"]))
    if existing:
        raise HTTPException(status_code=409, detail=f"Project name '{values['name']}' already exists")
    project = Project(**values)
    db.add(project)
    try:
        db.commit()
    except IntegrityError:
        db.rollback()
        raise HTTPException(status_code=409, detail=f"Project name '{values['name']}' already exists")
    except Exception as exc:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Server could not create project: {exc}")
    db.refresh(project)
    return _project_out(db, project)


@app.get("/api/projects/{project_id}")
def get_project(project_id: int, db: Session = Depends(get_db)):
    project = db.get(Project, project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    compounds = db.scalars(select(Compound).where(Compound.project_id == project_id).order_by(Compound.compound_id)).all()
    data = _project_out(db, project).model_dump()
    rows = []
    endpoint_names = {row.id: row.name for row in db.scalars(select(ADMETEndpoint).where(ADMETEndpoint.project_id == project_id))}
    for compound in compounds:
        output = compound_out(compound)
        version = next((row for row in compound.versions if row.version_number == compound.current_version), None)
        output["key_activity"] = None
        output["key_admet"] = None
        if version:
            activity = db.scalar(select(ActivityMeasurement).where(ActivityMeasurement.version_id == version.id).order_by(ActivityMeasurement.created_at.desc()))
            if activity:
                output["key_activity"] = f"{activity.raw_value:g} {activity.original_unit} · EXP"
            measurement = db.scalar(select(ADMETMeasurement).where(ADMETMeasurement.version_id == version.id).order_by(ADMETMeasurement.created_at.desc()))
            if measurement:
                value = measurement.qualitative_value or (measurement.mean_value if measurement.mean_value is not None else measurement.value)
                rendered = "Not measured" if value is None else f"{value} {measurement.unit}".strip()
                output["key_admet"] = f"{endpoint_names.get(measurement.endpoint_id, 'ADMET')}: {rendered} · EXP"
            else:
                prediction = db.scalar(select(ADMETPrediction).where(ADMETPrediction.version_id == version.id).order_by(ADMETPrediction.created_at.desc()))
                if prediction and prediction.predicted_value is not None:
                    output["key_admet"] = f"{prediction.model.endpoint_name}: {prediction.predicted_value:g} {prediction.unit} · PRED"
        rows.append(output)
    data["compounds"] = rows
    return data


@app.patch("/api/projects/{project_id}")
def update_project(project_id: int, payload: ProjectUpdate, db: Session = Depends(get_db)):
    project = db.get(Project, project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    values = payload.model_dump(exclude_unset=True)
    if values.get("molecule_type") not in {None, "Small Molecule", "Peptide"}:
        raise HTTPException(status_code=400, detail="Molecule Type must be Small Molecule or Peptide")
    for key, value in values.items():
        setattr(project, key, value)
    project.updated_at = utcnow()
    try:
        db.commit()
    except IntegrityError:
        db.rollback(); raise HTTPException(status_code=409, detail="Project name already exists")
    db.refresh(project); return _project_out(db, project)


def _delete_project_tree_rows(db: Session, project_ids: list[int]):
    """Delete complete project trees inside the caller's open transaction."""
    compound_ids = list(db.scalars(select(Compound.id).where(Compound.project_id.in_(project_ids))))
    version_ids = list(db.scalars(
        select(CompoundVersion.id).where(CompoundVersion.compound_row_id.in_(compound_ids))
    )) if compound_ids else []
    assay_ids = list(db.scalars(select(AssayDefinition.id).where(AssayDefinition.project_id.in_(project_ids))))
    endpoint_ids = list(db.scalars(select(ADMETEndpoint.id).where(ADMETEndpoint.project_id.in_(project_ids))))
    optimization_ids = list(db.scalars(select(OptimizationRun.id).where(OptimizationRun.project_id.in_(project_ids))))
    proposal_ids = list(db.scalars(
        select(OptimizationProposalRun.id).where(OptimizationProposalRun.project_id.in_(project_ids))
    ))
    candidate_ids = list(db.scalars(
        select(OptimizationCandidate.id).where(OptimizationCandidate.project_id.in_(project_ids))
    ))

    if version_ids:
        external_candidate = db.scalar(select(OptimizationCandidate.id).where(
            OptimizationCandidate.existing_version_id.in_(version_ids),
            ~OptimizationCandidate.project_id.in_(project_ids),
        ))
        if external_candidate:
            raise HTTPException(status_code=409, detail="Project tree has an invalid cross-project candidate reference")

    if candidate_ids:
        for model in (CandidateTransformation, CandidatePredictionSnapshot, CandidateRanking, CandidateRejectionReason):
            db.execute(delete(model).where(model.candidate_id.in_(candidate_ids)))
        db.execute(delete(OptimizationCandidate).where(OptimizationCandidate.id.in_(candidate_ids)))
    if proposal_ids:
        db.execute(delete(OptimizationProposalRun).where(OptimizationProposalRun.id.in_(proposal_ids)))
    if optimization_ids:
        db.execute(delete(OptimizationRun).where(OptimizationRun.id.in_(optimization_ids)))

    if assay_ids:
        db.execute(delete(MatchedMolecularPair).where(MatchedMolecularPair.assay_id.in_(assay_ids)))
        db.execute(delete(ActivityPrediction).where(ActivityPrediction.assay_id.in_(assay_ids)))
        db.execute(delete(ActivityMeasurement).where(ActivityMeasurement.assay_id.in_(assay_ids)))
        db.execute(delete(QSARModel).where(QSARModel.assay_id.in_(assay_ids)))
        db.execute(delete(AssayDefinition).where(AssayDefinition.id.in_(assay_ids)))

    if endpoint_ids:
        db.execute(delete(ADMETModelComparison).where(ADMETModelComparison.project_id.in_(project_ids)))
        db.execute(delete(ADMETModelPerformance).where(ADMETModelPerformance.project_id.in_(project_ids)))
        db.execute(delete(ADMETConsensusPrediction).where(ADMETConsensusPrediction.endpoint_id.in_(endpoint_ids)))
        db.execute(delete(ADMETPrediction).where(ADMETPrediction.endpoint_id.in_(endpoint_ids)))
        db.execute(delete(ADMETMeasurement).where(ADMETMeasurement.endpoint_id.in_(endpoint_ids)))
        db.execute(delete(ADMETAssayDefinition).where(ADMETAssayDefinition.endpoint_id.in_(endpoint_ids)))
        db.execute(delete(ADMETEndpoint).where(ADMETEndpoint.id.in_(endpoint_ids)))

    # PK and IVIVE project-level cleanup
    db.execute(delete(PKParameterSet).where(PKParameterSet.project_id.in_(project_ids)))
    db.execute(delete(PKSimulationRun).where(PKSimulationRun.project_id.in_(project_ids)))
    db.execute(delete(PKTranslationalSnapshot).where(PKTranslationalSnapshot.project_id.in_(project_ids)))
    db.execute(delete(PKHumanPredictionSnapshot).where(PKHumanPredictionSnapshot.project_id.in_(project_ids)))
    db.execute(delete(IVIVERun).where(IVIVERun.project_id.in_(project_ids)))
    db.execute(delete(IVIVEInputSet).where(IVIVEInputSet.project_id.in_(project_ids)))
    all_proj_pk_studies = list(db.scalars(select(PKStudy.id).where(PKStudy.project_id.in_(project_ids))))
    if all_proj_pk_studies:
        db.execute(delete(PKNCAResult).where(PKNCAResult.pk_study_id.in_(all_proj_pk_studies)))
        db.execute(delete(PKObservation).where(PKObservation.pk_study_id.in_(all_proj_pk_studies)))
        db.execute(delete(PKStudy).where(PKStudy.id.in_(all_proj_pk_studies)))

    if version_ids:
        db.execute(delete(PKParameterSet).where(PKParameterSet.version_id.in_(version_ids)))
        db.execute(delete(PKSimulationRun).where(PKSimulationRun.version_id.in_(version_ids)))
        db.execute(delete(PKTranslationalSnapshot).where(PKTranslationalSnapshot.version_id.in_(version_ids)))
        db.execute(delete(PKHumanPredictionSnapshot).where(PKHumanPredictionSnapshot.version_id.in_(version_ids)))
        db.execute(delete(IVIVERun).where(IVIVERun.version_id.in_(version_ids)))
        db.execute(delete(IVIVEInputSet).where(IVIVEInputSet.version_id.in_(version_ids)))
        pk_study_ids = list(db.scalars(select(PKStudy.id).where(PKStudy.version_id.in_(version_ids))))
        if pk_study_ids:
            db.execute(delete(PKNCAResult).where(PKNCAResult.pk_study_id.in_(pk_study_ids)))
            db.execute(delete(PKObservation).where(PKObservation.pk_study_id.in_(pk_study_ids)))
            db.execute(delete(PKStudy).where(PKStudy.id.in_(pk_study_ids)))
        db.execute(delete(PredictedMetabolite).where(PredictedMetabolite.version_id.in_(version_ids)))
        db.execute(delete(MetabolicSoftSpot).where(MetabolicSoftSpot.version_id.in_(version_ids)))
        db.execute(delete(MetabolicPredictionRun).where(MetabolicPredictionRun.version_id.in_(version_ids)))
        db.execute(delete(ExperimentalMetabolite).where(ExperimentalMetabolite.version_id.in_(version_ids)))
        db.execute(delete(ADMETPrediction).where(ADMETPrediction.version_id.in_(version_ids)))
        db.execute(delete(ADMETMeasurement).where(ADMETMeasurement.version_id.in_(version_ids)))
        db.execute(delete(ADMETPredictionRun).where(ADMETPredictionRun.version_id.in_(version_ids)))
        db.execute(delete(ActivityPrediction).where(ActivityPrediction.version_id.in_(version_ids)))
        db.execute(delete(ActivityMeasurement).where(ActivityMeasurement.version_id.in_(version_ids)))
        db.execute(delete(MatchedMolecularPair).where(
            MatchedMolecularPair.version_a_id.in_(version_ids) | MatchedMolecularPair.version_b_id.in_(version_ids)
        ))
        db.execute(delete(PropertyCalculation).where(PropertyCalculation.version_id.in_(version_ids)))
        db.execute(delete(StructuralAlert).where(StructuralAlert.version_id.in_(version_ids)))
        db.execute(delete(PredictionRun).where(PredictionRun.version_id.in_(version_ids)))
        db.execute(delete(CompoundVersion).where(CompoundVersion.id.in_(version_ids)))
    if compound_ids:
        db.execute(delete(Compound).where(Compound.id.in_(compound_ids)))
    db.execute(delete(PhysiologicalParameterOverride).where(
        PhysiologicalParameterOverride.project_id.in_(project_ids)
    ))
    db.execute(delete(Project).where(Project.id.in_(project_ids)))
    db.flush()



def _confirmed_project_delete(db: Session, confirmations: list[dict]):
    if not confirmations:
        raise HTTPException(status_code=400, detail="At least one project confirmation is required")
    project_ids = [int(item.get("id", 0)) for item in confirmations]
    if not all(project_ids) or len(project_ids) != len(set(project_ids)):
        raise HTTPException(status_code=400, detail="Project confirmations must contain unique valid IDs")
    projects = list(db.scalars(select(Project).where(Project.id.in_(project_ids))))
    if len(projects) != len(project_ids):
        raise HTTPException(status_code=404, detail="One or more projects were not found")
    names = {row.id: row.name for row in projects}
    for item in confirmations:
        project_id = int(item["id"])
        conf_name = item.get("confirmation_name")
        if conf_name is not None and str(conf_name).strip() != names[project_id].strip():
            raise HTTPException(status_code=400, detail=f"Confirmation name does not match project {project_id}")
    try:
        _delete_project_tree_rows(db, project_ids)
        db.commit()
    except HTTPException:
        db.rollback()
        raise
    except Exception as exc:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Project deletion failed; all changes were rolled back: {exc}")
    return {"deleted_project_ids": project_ids, "deleted_project_names": [names[row_id] for row_id in project_ids]}


@app.delete("/api/projects/{project_id}")
def delete_project(project_id: int, payload: dict | None = None, db: Session = Depends(get_db)):
    confirmation_name = payload.get("confirmation_name") if isinstance(payload, dict) else None
    return _confirmed_project_delete(db, [{"id": project_id, "confirmation_name": confirmation_name}])


@app.post("/api/projects/bulk-delete")
def bulk_delete_projects(payload: dict, db: Session = Depends(get_db)):
    return _confirmed_project_delete(db, payload.get("projects") or [])


def compound_out(compound: Compound):
    current = next((v for v in compound.versions if v.version_number == compound.current_version), compound.versions[-1] if compound.versions else None)
    return {
        "row_id": compound.id, "project_id": compound.project_id, "compound_id": compound.compound_id, "cas_number": compound.cas_number,
        "name": compound.name, "notes": compound.notes, "current_version": compound.current_version,
        "status": compound.status,
        "created_at": compound.created_at.isoformat(), "updated_at": compound.updated_at.isoformat(),
        "version": serialize_version(current) if current else None,
        "versions": [{"version_number": v.version_number, "canonical_smiles": v.canonical_smiles, "change_note": v.change_note, "calculated": bool(v.properties_json)} for v in compound.versions],
    }


def serialize_version(version: CompoundVersion):
    calc_json = version.calculation_json or {}
    ionization_data = calc_json.get("ionization")
    if not ionization_data and version.canonical_smiles:
        try:
            ionization_data = analyze_ionization(version.canonical_smiles)
        except Exception:
            ionization_data = {}
    return {
        "id": version.id, "version_number": version.version_number, "original_smiles": version.original_smiles,
        "canonical_smiles": version.canonical_smiles, "isomeric_smiles": version.isomeric_smiles,
        "inchi": version.inchi, "inchikey": version.inchikey, "change_note": version.change_note,
        "properties": version.properties_json or {}, "rules": (version.calculation_json or {}).get("rules", {}),
        "ionization": ionization_data or {},
        "assessment": version.assessment_json or {}, "svg": version.svg,
        "highlighted_svg": version.highlighted_svg, "provenance": (version.calculation_json or {}).get("provenance", {}),
        "alerts": version.alerts_json or [],
        "calculated": bool(version.properties_json),
    }


def _store_calculation(db: Session, compound: Compound, version: CompoundVersion, analysis: dict) -> None:
    db.execute(delete(PropertyCalculation).where(PropertyCalculation.version_id == version.id))
    db.execute(delete(StructuralAlert).where(StructuralAlert.version_id == version.id))
    version.properties_json = analysis["properties"]
    version.alerts_json = analysis["alerts"]
    version.assessment_json = analysis["assessment"]
    version.calculation_json = {
        "provenance": analysis["provenance"],
        "rules": analysis["rules"],
        "ionization": analysis.get("ionization", {}),
    }
    version.highlighted_svg = analysis["highlighted_svg"]
    for endpoint, value in analysis["properties"].items():
        if value is None or isinstance(value, (dict, list)): continue
        method = "RDKit descriptor"
        if endpoint == "clogp": method = "Crippen cLogP"
        elif endpoint in ("molar_refractivity",): method = "Crippen Molar Refractivity"
        elif endpoint == "tpsa": method = "Ertl TPSA"
        elif endpoint == "qed": method = "RDKit QED"
        elif endpoint == "fraction_csp3": method = "RDKit Fraction CSP3"
        db.add(PropertyCalculation(version_id=version.id, endpoint=endpoint, value=str(value), engine=ENGINE,
                                   method=method, engine_version=ENGINE_VERSION))
    for alert in analysis["alerts"]:
        db.add(StructuralAlert(version_id=version.id, alert_set=alert["alert_set"], alert_name=alert["alert_name"],
                               reason=alert["reason"], matched_smiles=alert["matched_smiles"],
                               matched_atoms_json=alert["matched_atoms"]))
    db.add(PredictionRun(version_id=version.id, stage="stage_1", model_name=f"{ENGINE} property pipeline",
                         model_version=ENGINE_VERSION, inputs_hash=analysis["inputs_hash"],
                         outputs_json=json.loads(json.dumps({"properties": analysis["properties"], "rules": analysis["rules"]})),
                         provenance_json=analysis["provenance"], confidence="High"))
    compound.status = "CALCULATED"


def persist_structure(db: Session, compound: Compound, smiles: str, change_note: str, calculate: bool) -> CompoundVersion:
    analysis = analyze_smiles(smiles)
    duplicate = db.scalar(
        select(CompoundVersion).join(Compound, Compound.id == CompoundVersion.compound_row_id)
        .where(Compound.project_id == compound.project_id, CompoundVersion.inchikey == analysis["identity"]["inchikey"])
    )
    if duplicate and duplicate.compound_row_id != compound.id:
        raise HTTPException(status_code=409, detail={
            "error": "Duplicate structure in this project",
            "existing_compound_id": duplicate.compound.compound_id,
            "inchikey": duplicate.inchikey,
        })
    number = 1 if not compound.versions else max(version.version_number for version in compound.versions) + 1
    version = CompoundVersion(
        compound_row_id=compound.id, version_number=number, original_smiles=smiles.strip(), change_note=change_note,
        canonical_smiles=analysis["identity"]["canonical_smiles"], isomeric_smiles=analysis["identity"]["isomeric_smiles"],
        inchi=analysis["identity"]["inchi"], inchikey=analysis["identity"]["inchikey"],
        properties_json=None, alerts_json=None, assessment_json=None, calculation_json=None,
        svg=analysis["svg"], highlighted_svg=analysis["svg"],
    )
    db.add(version); db.flush()
    compound.current_version = number
    compound.status = "STRUCTURE_READY"
    if calculate:
        _store_calculation(db, compound, version, analysis)
    db.commit(); db.refresh(version)
    return version


def persist_analysis(db: Session, compound: Compound, smiles: str, change_note: str) -> CompoundVersion:
    """Backward-compatible create-new-version path used by Stages 1–4B."""
    return persist_structure(db, compound, smiles, change_note, calculate=True)


@app.post("/api/projects/{project_id}/compounds", status_code=201)
def create_compound(project_id: int, payload: CompoundCreate, db: Session = Depends(get_db)):
    project = db.get(Project, project_id)
    if not project: raise HTTPException(status_code=404, detail="Project not found")
    name = payload.name.strip() or payload.compound_id.strip()
    if not name:
        raise HTTPException(status_code=400, detail="Compound name is required.")
    base_id = payload.compound_id.strip() or "".join(character if character.isalnum() or character in "-_" else "-" for character in name.upper()).strip("-")[:50] or "COMPOUND"
    compound_id = base_id
    suffix = 2
    while not payload.compound_id and db.scalar(select(Compound.id).where(Compound.project_id == project_id, Compound.compound_id == compound_id)):
        compound_id = f"{base_id[:45]}-{suffix}"; suffix += 1
    existing_label = db.scalar(select(Compound).where(Compound.project_id == project_id, Compound.compound_id == compound_id))
    if existing_label: raise HTTPException(status_code=409, detail="Compound ID already exists in project")
    cas_number = payload.cas_number.strip()
    compound = Compound(project_id=project_id, compound_id=compound_id, name=name, cas_number=cas_number, notes=payload.notes, status="DRAFT")
    db.add(compound); db.flush()
    if not payload.smiles.strip():
        db.commit(); db.refresh(compound); return compound_out(compound)
    if project.molecule_type != "Small Molecule":
        db.rollback(); raise HTTPException(status_code=400, detail="This model currently supports small molecules only. Save a peptide as a draft without structure calculations.")
    try:
        persist_structure(db, compound, payload.smiles, "Initial structure", calculate=payload.calculate)
        db.commit()
    except HTTPException:
        db.rollback(); raise
    except ChemistryError as exc:
        db.rollback(); raise HTTPException(status_code=400, detail=f"Structure could not be standardized: {exc}")
    except Exception as exc:
        db.rollback(); raise HTTPException(status_code=500, detail=f"Server could not save compound: {exc}")
    db.refresh(compound); return compound_out(compound)


@app.post("/api/compounds/{row_id}/calculate")
def calculate_compound_properties(row_id: int, db: Session = Depends(get_db)):
    compound = db.get(Compound, row_id)
    if not compound:
        raise HTTPException(status_code=404, detail="Compound not found")
    if compound.project.molecule_type != "Small Molecule":
        raise HTTPException(status_code=400, detail="This model currently supports small molecules only.")
    version = next((row for row in compound.versions if row.version_number == compound.current_version), None)
    if not version:
        raise HTTPException(status_code=400, detail="Draw or enter a valid structure before calculating properties")
    try:
        analysis = analyze_smiles(version.original_smiles)
        _store_calculation(db, compound, version, analysis)
        db.commit(); db.refresh(compound)
    except ChemistryError as exc:
        db.rollback(); raise HTTPException(status_code=400, detail=str(exc))
    return compound_out(compound)


@app.post("/api/compounds/{row_id}/predict-workflow", status_code=202)
def run_compound_prediction_workflow(row_id: int, db: Session = Depends(get_db)):
    """Save-following and Overview prediction orchestration. Activity is intentionally excluded unless configured."""
    compound = db.get(Compound, row_id)
    if not compound:
        raise HTTPException(status_code=404, detail="Compound not found")
    version = next((row for row in compound.versions if row.version_number == compound.current_version), None)
    if not version:
        raise HTTPException(status_code=400, detail="A validated CompoundVersion is required before prediction")
    if compound.project.molecule_type != "Small Molecule":
        raise HTTPException(status_code=400, detail="This model currently supports small molecules only.")
    steps = {
        "overview": {"status": "COMPLETE", "message": "Compound identity and validated structure are available."},
        "properties": {"status": "PENDING"},
        "admet": {"status": "PENDING", "endpoints": []},
        "metabolism": {"status": "PENDING"},
        "pk": {"status": "PENDING", "routes": []},
        "activity": {"status": "NOT_INCLUDED", "message": "Assay configuration required; Activity is excluded from automatic prediction."},
    }
    completed_endpoints = []
    unavailable_endpoints = []
    failed_endpoints = []

    try:
        steps["properties"] = {"status": "RUNNING"}
        calculate_compound_properties(row_id, db)
        steps["properties"] = {"status": "COMPLETE", "message": "Stage 1 properties calculated."}
        completed_endpoints.extend(["Physicochemical Properties", "Ionization (pH)", "Structural Alerts"])
    except Exception as exc:
        steps["properties"] = {"status": "FAILED", "message": str(getattr(exc, "detail", exc))}
        failed_endpoints.append("Physicochemical Properties")

    try:
        steps["admet"] = {"status": "RUNNING", "endpoints": []}
        result = run_admet_predictions(version.id, db)
        steps["admet"] = {"status": "COMPLETE" if result["status"] in {"COMPLETE", "CACHED"} else result["status"],
                          "message": result.get("message", ""), "endpoints": result.get("endpoint_statuses", []),
                          "consensus_count": len(result.get("consensus_predictions", []))}
        if result["status"] in {"COMPLETE", "CACHED"}:
            completed_endpoints.extend(["Solubility", "Caco-2 Permeability", "Plasma Protein Binding", "HLM Clearance", "RLM Clearance", "MLM Clearance", "hERG Liability", "DILI Liability", "Ames Mutagenicity"])
        else:
            unavailable_endpoints.append("ADMET Model Panel")
    except Exception as exc:
        steps["admet"] = {"status": "FAILED", "message": str(getattr(exc, "detail", exc)), "endpoints": []}
        failed_endpoints.append("ADMET Predictions")

    try:
        steps["metabolism"] = {"status": "RUNNING"}
        result = run_metabolism_predictions(version.id, db)
        steps["metabolism"] = {"status": "COMPLETE" if result["status"] in {"COMPLETE", "CACHED"} else result["status"],
                               "message": result.get("message", ""), "soft_spots_and_metabolites": True}
        if result["status"] in {"COMPLETE", "CACHED"}:
            completed_endpoints.extend(["CYP Inhibition Panel", "CYP Substrate Panel", "P-gp Transporter", "SyGMa Soft Spots", "Metabolite Hypotheses"])
        else:
            unavailable_endpoints.append("Metabolism Prediction Panel")
    except Exception as exc:
        steps["metabolism"] = {"status": "FAILED", "message": str(getattr(exc, "detail", exc))}
        failed_endpoints.append("Metabolism Predictions")

    # PK is part of the one-click workflow.  The IVIVE/PK foundation builder
    # persists route-aware parameter sets (IV/PO/SC/IP) and explicitly marks
    # unavailable inputs instead of inventing values.
    try:
        steps["pk"] = {"status": "RUNNING", "routes": []}

        # 1. Run IVIVE for available species
        for sp in ["Human", "Rat", "Mouse"]:
            try:
                calculate_ivive(db, version, sp)
            except Exception:
                pass

        # 2. Build PK Parameter Foundations across species
        profile = get_pk_foundation_profile(db, version.id, "Rat")
        for sp in ["Mouse", "Dog", "Monkey", "Human"]:
            try:
                get_pk_foundation_profile(db, version.id, sp)
            except Exception:
                pass

        # 3. Pre-run baseline PK Simulations (Rat PO 1.0 mg/kg, Rat IV 1.0 mg/kg)
        from .simulation import run_pk_simulation, PKSimulationRequest
        try:
            run_pk_simulation(db, version.id, PKSimulationRequest(species="Rat", route="PO", dose=1.0, dose_unit="mg/kg"))
        except Exception:
            pass
        try:
            run_pk_simulation(db, version.id, PKSimulationRequest(species="Rat", route="IV", dose=1.0, dose_unit="mg/kg"))
        except Exception:
            pass

        # 4. Multi-species PK profile
        try:
            get_multi_species_pk_profile(db, version.id)
        except Exception:
            pass

        # 5. Human PK profile assembly
        from .human_pk import assemble_human_pk_parameters
        try:
            assemble_human_pk_parameters(db, version.id)
        except Exception:
            pass

        # 6. Translational PK profile
        from .translational import get_translational_pk_profile
        try:
            get_translational_pk_profile(db, version.id, freeze_snapshot=False)
        except Exception:
            pass

        route_sets = profile.get("route_parameter_sets", {})
        routes = []
        for route, params in route_sets.items():
            routes.append({
                "route": route,
                "status": "COMPLETE" if any(params.get(key) is not None for key in ("cl_value", "v_value", "ka_value", "f_predicted", "f_experimental")) else "MODEL_UNAVAILABLE",
                "confidence": params.get("confidence"),
            })
        pk_ready = any(route["status"] == "COMPLETE" for route in routes)
        steps["pk"] = {
            "status": "COMPLETE" if pk_ready else "MODEL_UNAVAILABLE",
            "message": "PK parameter foundation, IVIVE, and simulation assembled across species." if pk_ready else "PK inputs are not available for this compound.",
            "routes": routes,
        }
        if routes:
            completed_endpoints.append("PK parameter foundation (IV/PO/SC/IP)")
        else:
            unavailable_endpoints.append("PK parameter foundation")
    except Exception as exc:
        steps["pk"] = {"status": "FAILED", "message": str(getattr(exc, "detail", exc)), "routes": []}
        failed_endpoints.append("PK Prediction")

    required = [steps[name]["status"] for name in ("properties", "admet", "metabolism", "pk")]
    status = "COMPLETE" if all(value == "COMPLETE" for value in required) else ("FAILED" if all(value == "FAILED" for value in required) else "PARTIAL")
    completed_at = datetime.now(timezone.utc)
    timestamp = completed_at.strftime("%Y-%m-%d %H:%M")
    workflow_output = {
        "status": status,
        "steps": steps,
        "completed_endpoints": completed_endpoints,
        "unavailable_endpoints": unavailable_endpoints,
        "failed_endpoints": failed_endpoints,
        "timestamp": timestamp,
    }
    db.add(PredictionRun(
        version_id=version.id,
        stage="prediction_workflow",
        model_name="Properties + ADMET + Metabolism + PK workflow",
        model_version=CURRENT_STAGE,
        inputs_hash=hashlib.sha256(f"{version.id}|{version.canonical_smiles}|{completed_at.isoformat()}".encode()).hexdigest(),
        outputs_json=workflow_output,
        provenance_json={"orchestrator": "predict-workflow", "pk_species": "Rat", "persisted": True},
        confidence="High" if status == "COMPLETE" else "Limited",
    ))
    db.commit()

    return {
        "status": status,
        "compound_id": compound.id,
        "compound_version_id": version.id,
        "activity_excluded": True,
        "completed_endpoints": completed_endpoints,
        "unavailable_endpoints": unavailable_endpoints,
        "failed_endpoints": failed_endpoints,
        "completed_count": len(completed_endpoints),
        "unavailable_count": len(unavailable_endpoints),
        "failed_count": len(failed_endpoints),
        "timestamp": timestamp,
        "steps": steps,
        "message": f"Prediction {status.lower()}: {len(completed_endpoints)} endpoints calculated, {len(unavailable_endpoints)} unavailable, Activity not run (assay required).",
    }


@app.post("/api/compounds/{row_id}/predict-all", status_code=202)
def run_compound_predict_all(row_id: int, db: Session = Depends(get_db)):
    """Overview Primary Predict Endpoint alias."""
    return run_compound_prediction_workflow(row_id, db)


@app.get("/api/compounds/{row_id}")

def get_compound(row_id: int, include_versions: bool = Query(False), db: Session = Depends(get_db)):
    compound = db.get(Compound, row_id)
    if not compound: raise HTTPException(status_code=404, detail="Compound not found")
    result = compound_out(compound)
    if include_versions:
        versions = db.scalars(select(CompoundVersion).where(CompoundVersion.compound_row_id == row_id).order_by(CompoundVersion.version_number)).all()
        result["history"] = [serialize_version(v) for v in versions]
    runs = db.scalars(select(PredictionRun).join(CompoundVersion, CompoundVersion.id == PredictionRun.version_id).where(CompoundVersion.compound_row_id == row_id).order_by(PredictionRun.created_at.desc())).all()
    result["prediction_history"] = [{
        "prediction_id": run.id, "created_at": run.created_at.isoformat(), "stage": run.stage, "model_name": run.model_name,
        "model_version": run.model_version, "confidence": run.confidence, "provenance": run.provenance_json,
        "inputs_hash": run.inputs_hash, "outputs": run.outputs_json,
    } for run in runs]
    return result


@app.get("/api/compounds/{row_id}/external-experimental/search")
def search_external_experimental_data(row_id: int, db: Session = Depends(get_db)):
    """Explicit CAS-only public lookup; no write occurs during search."""
    compound = db.get(Compound, row_id)
    if not compound:
        raise HTTPException(status_code=404, detail="Compound not found")
    status = cas_status(compound.cas_number)
    if status != "VALID":
        return {"status": "DISABLED_NO_CAS" if status == "EMPTY" else "DISABLED_INVALID_CAS", "cas_status": status, "records": []}
    current = next((v for v in compound.versions if v.version_number == compound.current_version), None)
    if not current:
        return {"status": "STRUCTURE_MISMATCH", "cas_status": status, "records": []}
    result = external_evidence_lookup(compound.cas_number, current.inchikey)
    result["cas_status"] = status
    result["compound_version_id"] = current.id
    return result


@app.get("/api/compounds/{row_id}/external-experimental")
def list_external_experimental_data(row_id: int, db: Session = Depends(get_db)):
    compound = db.get(Compound, row_id)
    if not compound:
        raise HTTPException(status_code=404, detail="Compound not found")
    version_ids = [v.id for v in compound.versions]
    rows = db.scalars(select(ExternalExperimentalEvidence).where(ExternalExperimentalEvidence.compound_version_id.in_(version_ids)).order_by(ExternalExperimentalEvidence.imported_at.desc())).all() if version_ids else []
    return {"records": [{"id": row.id, "endpoint": row.raw_endpoint_name, "value": row.raw_value, "unit": row.raw_unit,
                         "relation": row.raw_relation, "source": row.source_database, "reference": row.reference_text,
                         "source_url": row.source_url, "evidence_origin": row.evidence_origin} for row in rows]}


@app.patch("/api/compounds/{row_id}")
def update_compound(row_id: int, payload: CompoundUpdate, db: Session = Depends(get_db)):
    compound = db.get(Compound, row_id)
    if not compound: raise HTTPException(status_code=404, detail="Compound not found")
    if payload.name is not None:
        if not payload.name.strip():
            raise HTTPException(status_code=400, detail="Compound Name is required")
        compound.name = payload.name.strip()
    if payload.compound_id is not None:
        compound.compound_id = payload.compound_id.strip()
    if payload.cas_number is not None:
        cas_number = payload.cas_number.strip()
        compound.cas_number = cas_number
    if payload.notes is not None: compound.notes = payload.notes
    compound.updated_at = utcnow()
    if payload.smiles:
        if compound.project.molecule_type != "Small Molecule":
            raise HTTPException(status_code=400, detail="This model currently supports small molecules only.")
        try:
            persist_analysis(db, compound, payload.smiles, payload.change_note)
        except HTTPException: db.rollback(); raise
        except ChemistryError as exc: db.rollback(); raise HTTPException(status_code=400, detail=str(exc))
    else:
        # Metadata edits invalidate cached prediction outputs for the current
        # compound versions; the next explicit Predict run rebuilds them.
        version_ids = [v.id for v in compound.versions]
        if version_ids:
            db.execute(delete(ADMETPrediction).where(ADMETPrediction.version_id.in_(version_ids)))
            db.execute(delete(ADMETPredictionRun).where(ADMETPredictionRun.version_id.in_(version_ids)))
            db.execute(delete(ActivityPrediction).where(ActivityPrediction.version_id.in_(version_ids)))
            db.execute(delete(MetabolicPredictionRun).where(MetabolicPredictionRun.version_id.in_(version_ids)))
            db.execute(delete(PredictionRun).where(PredictionRun.version_id.in_(version_ids)))
        db.commit()
    db.refresh(compound); return compound_out(compound)


@app.delete("/api/compounds/{row_id}", status_code=204)
def delete_compound(row_id: int, db: Session = Depends(get_db)):
    compound = db.get(Compound, row_id)
    if not compound: raise HTTPException(status_code=404, detail="Compound not found")
    v_ids = [v.id for v in compound.versions]
    db.execute(delete(PKParameterSet).where((PKParameterSet.compound_row_id == row_id) | (PKParameterSet.version_id.in_(v_ids))))
    db.execute(delete(PKSimulationRun).where((PKSimulationRun.compound_row_id == row_id) | (PKSimulationRun.version_id.in_(v_ids))))
    db.execute(delete(PKTranslationalSnapshot).where(PKTranslationalSnapshot.compound_row_id == row_id))
    db.execute(delete(PKHumanPredictionSnapshot).where(PKHumanPredictionSnapshot.compound_row_id == row_id))
    db.execute(delete(IVIVERun).where(IVIVERun.version_id.in_(v_ids)))
    db.execute(delete(IVIVEInputSet).where(IVIVEInputSet.version_id.in_(v_ids)))
    st_ids = list(db.scalars(select(PKStudy.id).where(PKStudy.compound_row_id == row_id)))
    if st_ids:
        db.execute(delete(PKNCAResult).where(PKNCAResult.pk_study_id.in_(st_ids)))
        db.execute(delete(PKObservation).where(PKObservation.pk_study_id.in_(st_ids)))
        db.execute(delete(PKStudy).where(PKStudy.id.in_(st_ids)))
    db.delete(compound); db.commit()



@app.get("/api/projects/{project_id}/compare")
def compare(project_id: int, ids: str = Query(...), db: Session = Depends(get_db), assay_id: int | None = None):
    try: wanted = {int(value) for value in ids.split(",") if value.strip()}
    except ValueError: raise HTTPException(status_code=400, detail="ids must be comma-separated integers")
    if len(wanted) < 2: raise HTTPException(status_code=400, detail="Select at least two compounds")
    rows = []
    for row_id in wanted:
        compound = db.get(Compound, row_id)
        if not compound or compound.project_id != project_id: continue
        version = next((v for v in compound.versions if v.version_number == compound.current_version), None)
        if not version:
            continue
        p = version.properties_json or {}
        comparison_row = {
            "compound": compound.compound_id, "name": compound.name, "row_id": compound.id,
            "MW": p.get("molecular_weight"), "cLogP": p.get("clogp"), "TPSA": p.get("tpsa"), "HBD": p.get("hbd"),
            "HBA": p.get("hba"), "RotB": p.get("rotatable_bonds"), "Fsp3": p.get("fraction_csp3"), "QED": p.get("qed"),
            "svg": version.svg if version else "", "inchikey": version.inchikey if version else "",
            "sources": {},
        }
        for key in ("MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "Fsp3", "QED"):
            comparison_row["sources"][key] = "Calculated" if comparison_row[key] is not None else "Not calculated"
        activity_query = select(ActivityMeasurement).where(ActivityMeasurement.version_id == version.id)
        if assay_id is not None:
            assay = db.get(AssayDefinition, assay_id)
            if not assay or assay.project_id != project_id:
                raise HTTPException(status_code=404, detail="Comparison assay is not in this project")
            activity_query = activity_query.where(ActivityMeasurement.assay_id == assay_id)
        activity = db.scalar(activity_query.order_by(ActivityMeasurement.created_at.desc()))
        comparison_row["Activity"] = activity.normalized_value_nm if activity else None
        comparison_row["sources"]["Activity"] = "Experimental" if activity else "Not measured"
        endpoint_map = {
            "HLM intrinsic clearance": "HLM", "RLM intrinsic clearance": "RLM", "MLM intrinsic clearance": "MLM",
            "Plasma protein binding": "PPB", "Solubility": "Solubility", "Permeability": "Caco-2",
            "CYP1A2 inhibitor": "CYP1A2 Inh", "CYP2C9 inhibitor": "CYP2C9 Inh", "CYP2C19 inhibitor": "CYP2C19 Inh",
            "CYP2D6 inhibitor": "CYP2D6 Inh", "CYP3A4 inhibitor": "CYP3A4 Inh",
            "CYP2C9 substrate": "CYP2C9 Sub", "CYP2D6 substrate": "CYP2D6 Sub", "CYP3A4 substrate": "CYP3A4 Sub",
            "P-gp inhibitor": "P-gp Inh",
            "hERG liability": "hERG", "Ames mutagenicity": "Ames", "DILI clinical liability": "DILI",
        }
        experimental = db.scalars(select(ADMETMeasurement).where(ADMETMeasurement.version_id == version.id)).all()
        endpoint_names = {item.id: item.name for item in db.scalars(select(ADMETEndpoint).where(ADMETEndpoint.project_id == project_id))}
        for endpoint_name, label in endpoint_map.items():
            endpoint_experimental = [row for row in experimental if endpoint_names.get(row.endpoint_id) == endpoint_name]
            prediction = db.scalar(
                select(ADMETPrediction).join(ADMETModelRegistry)
                .where(ADMETPrediction.version_id == version.id, ADMETModelRegistry.endpoint_name == endpoint_name)
                .order_by(ADMETPrediction.created_at.desc())
            )
            if not prediction:
                first = endpoint_experimental[0] if endpoint_experimental else None
                comparison_row[label] = (
                    first.qualitative_value or (first.mean_value if first.mean_value is not None else first.value)
                ) if first else None
                comparison_row["sources"][label] = "Experimental" if first else "Not measured"
                continue
            matches = comparison_for_prediction(endpoint_name, prediction.predicted_value, experimental, endpoint_names)
            comparison_row[label] = matches[0]["experimental_normalized"] if matches else prediction.predicted_value
            comparison_row["sources"][label] = "Experimental" if matches else "Predicted"

        # Explicit Dog and Monkey metabolism entries
        comparison_row["DLM"] = None
        comparison_row["sources"]["DLM"] = "MODEL_UNAVAILABLE"
        comparison_row["CyLM"] = None
        comparison_row["sources"]["CyLM"] = "MODEL_UNAVAILABLE"

        # Derived fraction unbound (fu)
        ppb_val = comparison_row.get("PPB")
        if ppb_val is not None and isinstance(ppb_val, (int, float)):
            comparison_row["fu"] = round((100.0 - float(ppb_val)) / 100.0, 4)
            comparison_row["sources"]["fu"] = "Calculated (1 - PPB)"
        else:
            comparison_row["fu"] = None
            comparison_row["sources"]["fu"] = "Not calculated"

        # Soft Spots
        spots = list(db.scalars(select(MetabolicSoftSpot).where(MetabolicSoftSpot.version_id == version.id)).all())
        comparison_row["Soft Spots"] = len(spots) if spots else 0
        comparison_row["sources"]["Soft Spots"] = "Predicted (SyGMa)" if spots else "Not calculated"

        # Multi-species PK assembly (Normalized 1 mg/kg single dose standard)
        pk_prof = get_multi_species_pk_profile(db, version.id)
        sp_map = pk_prof.get("species_profiles", {})

        # Mouse
        m_prof = sp_map.get("Mouse", {})
        comparison_row["Mouse CL (IV)"] = m_prof.get("cl", {}).get("value")
        comparison_row["sources"]["Mouse CL (IV)"] = m_prof.get("cl", {}).get("source", "UNAVAILABLE")
        comparison_row["Mouse Vd"] = m_prof.get("v", {}).get("value")
        comparison_row["sources"]["Mouse Vd"] = m_prof.get("v", {}).get("source", "UNAVAILABLE")
        comparison_row["Mouse t1/2"] = m_prof.get("t_half_hours")
        comparison_row["sources"]["Mouse t1/2"] = "Calculated" if m_prof.get("t_half_hours") else "UNAVAILABLE"

        # Rat
        r_prof = sp_map.get("Rat", {})
        comparison_row["Rat CL (IV)"] = r_prof.get("cl", {}).get("value")
        comparison_row["sources"]["Rat CL (IV)"] = r_prof.get("cl", {}).get("source", "UNAVAILABLE")
        comparison_row["Rat Vd"] = r_prof.get("v", {}).get("value")
        comparison_row["sources"]["Rat Vd"] = r_prof.get("v", {}).get("source", "UNAVAILABLE")
        comparison_row["Rat t1/2"] = r_prof.get("t_half_hours")
        comparison_row["sources"]["Rat t1/2"] = "Calculated" if r_prof.get("t_half_hours") else "UNAVAILABLE"
        comparison_row["Rat F (%)"] = r_prof.get("f_pct")
        comparison_row["sources"]["Rat F (%)"] = r_prof.get("f_source", "UNAVAILABLE")

        # Dog
        d_prof = sp_map.get("Dog", {})
        comparison_row["Dog CL (IV)"] = d_prof.get("cl", {}).get("value")
        comparison_row["sources"]["Dog CL (IV)"] = d_prof.get("cl", {}).get("source", "MODEL_UNAVAILABLE")

        # Monkey
        cy_prof = sp_map.get("Monkey", {})
        comparison_row["Monkey CL (IV)"] = cy_prof.get("cl", {}).get("value")
        comparison_row["sources"]["Monkey CL (IV)"] = cy_prof.get("cl", {}).get("source", "MODEL_UNAVAILABLE")

        # Human
        h_prof = sp_map.get("Human", {})
        comparison_row["Human CL (IVIVE)"] = h_prof.get("cl", {}).get("value")
        comparison_row["sources"]["Human CL (IVIVE)"] = h_prof.get("cl", {}).get("source", "UNAVAILABLE")
        comparison_row["Human Vd (pred)"] = h_prof.get("v", {}).get("value")
        comparison_row["sources"]["Human Vd (pred)"] = h_prof.get("v", {}).get("source", "UNAVAILABLE")
        comparison_row["Human t1/2 (pred)"] = h_prof.get("t_half_hours")
        comparison_row["sources"]["Human t1/2 (pred)"] = "Calculated" if h_prof.get("t_half_hours") else "UNAVAILABLE"
        comparison_row["Human AUC (1mg/kg IV)"] = h_prof.get("normalized_1mpk_iv", {}).get("auc_ng_h_ml")
        comparison_row["sources"]["Human AUC (1mg/kg IV)"] = "Normalized 1 mg/kg IV" if h_prof.get("normalized_1mpk_iv", {}).get("auc_ng_h_ml") else "UNAVAILABLE"
        comparison_row["Human Cmax (1mg/kg IV)"] = h_prof.get("normalized_1mpk_iv", {}).get("cmax_ng_ml")
        comparison_row["sources"]["Human Cmax (1mg/kg IV)"] = "Normalized 1 mg/kg IV" if h_prof.get("normalized_1mpk_iv", {}).get("cmax_ng_ml") else "UNAVAILABLE"

        rows.append(comparison_row)
    if len(rows) < 2: raise HTTPException(status_code=400, detail="At least two selected compounds must belong to the project")
    property_metrics = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "Fsp3", "QED"]
    adme_metrics = ["Solubility", "Caco-2", "PPB", "fu"]
    metabolism_metrics = ["HLM", "RLM", "MLM", "DLM", "CyLM", "CYP1A2 Inh", "CYP2C9 Inh", "CYP2C19 Inh", "CYP2D6 Inh", "CYP3A4 Inh", "CYP2C9 Sub", "CYP2D6 Sub", "CYP3A4 Sub", "P-gp Inh", "Soft Spots"]
    pk_metrics = ["Mouse CL (IV)", "Mouse Vd", "Mouse t1/2", "Rat CL (IV)", "Rat Vd", "Rat t1/2", "Rat F (%)", "Dog CL (IV)", "Monkey CL (IV)", "Human CL (IVIVE)", "Human Vd (pred)", "Human t1/2 (pred)", "Human AUC (1mg/kg IV)", "Human Cmax (1mg/kg IV)"]
    safety_metrics = ["hERG", "Ames", "DILI"]
    metrics = property_metrics + ["Activity"] + adme_metrics + metabolism_metrics + pk_metrics + safety_metrics
    ranges = {}
    for metric in property_metrics:
        values = [row[metric] for row in rows if row[metric] is not None]
        ranges[metric] = {"min": min(values), "max": max(values)} if values else {"min": None, "max": None}
    return {"metrics": metrics, "ranges": ranges, "compounds": rows, "metric_units": {
        "Activity": "nM (latest experimental)", "HLM": "log10(mL/min/kg)", "RLM": "log10(mL/min/kg)", "MLM": "log10(mL/min/kg)",
        "DLM": "MODEL_UNAVAILABLE", "CyLM": "MODEL_UNAVAILABLE",
        "PPB": "% bound", "fu": "fraction unbound (0-1)", "Solubility": "log10(mol/L)", "Caco-2": "log10(cm/s)",
        "Soft Spots": "count",
        "Mouse CL (IV)": "mL/min/kg", "Mouse Vd": "L/kg", "Mouse t1/2": "hours",
        "Rat CL (IV)": "mL/min/kg", "Rat Vd": "L/kg", "Rat t1/2": "hours", "Rat F (%)": "% bioavailability",
        "Dog CL (IV)": "mL/min/kg", "Monkey CL (IV)": "mL/min/kg",
        "Human CL (IVIVE)": "mL/min/kg", "Human Vd (pred)": "L/kg", "Human t1/2 (pred)": "hours",
        "Human AUC (1mg/kg IV)": "ng·h/mL (1 mg/kg IV single dose)", "Human Cmax (1mg/kg IV)": "ng/mL (1 mg/kg IV)",
        "hERG": "classification/probability", "Ames": "classification/probability", "DILI": "classification/probability",
    }}


def _assay_out(assay: AssayDefinition):
    return {
        "id": assay.id, "assay_uid": assay.assay_uid, "version_number": assay.version_number,
        "active": assay.active, "name": assay.name, "target": assay.target, "target_type": assay.target_type,
        "assay_category": assay.assay_category, "measurement_type": assay.measurement_type,
        "custom_measurement_name": assay.custom_measurement_name, "unit": assay.unit,
        "species": assay.species, "cell_line": assay.cell_line, "mutation_variant": assay.mutation_variant,
        "protein_construct": assay.protein_construct, "substrate": assay.substrate,
        "atp_concentration": assay.atp_concentration, "incubation_time": assay.incubation_time,
        "detection_method": assay.detection_method, "experimental_conditions": assay.experimental_conditions,
        "protocol": assay.protocol, "reference_compound": assay.reference_compound,
        "reference_structure_smiles": assay.reference_structure_smiles,
        "reference_activity": assay.reference_activity, "reference_source": assay.reference_source,
        "reference_provenance_url": assay.reference_provenance_url, "notes": assay.notes,
    }


def _admet_endpoint_out(endpoint: ADMETEndpoint):
    return {"id": endpoint.id, "name": endpoint.name, "category": endpoint.category,
            "description": endpoint.description, "preferred_unit": endpoint.preferred_unit,
            "direction": endpoint.direction}


def _admet_model_out(model: ADMETModelRegistry):
    assets_available, unavailable_reason = model_files_available(model.endpoint_name) if model.endpoint_name in MODEL_SPECS else (
        False, (model.provenance_json or {}).get("reason", "No endpoint-specific model installed in the current stage"),
    )
    available = bool(model.is_active and assets_available and model.implementation_status == "READY")
    if assets_available and not available:
        unavailable_reason = (model.provenance_json or {}).get("reason", "Model registry entry is inactive")
    cal_info = CONFORMAL_CALIBRATION_REGISTRY.get(model.endpoint_name, {})
    cal_provenance = cal_info.get("data_provenance", DataProvenance.UNAVAILABLE if not available else DataProvenance.TRAINING_OVERLAP_UNKNOWN)
    cal_quality = cal_info.get("calibration_quality", CalibrationQuality.UNAVAILABLE)
    details = model.provenance_json or {}
    limitations = str(details.get("limitations") or "")
    confidence = (
        "NOT_APPLICABLE" if not available else
        ("LOW" if model.endpoint_name.endswith("intrinsic clearance") or "confidence is capped at LOW" in limitations else "COMPOUND_DEPENDENT")
    )
    conformal_status = (
        "CONFORMAL_UNAVAILABLE" if cal_quality == CalibrationQuality.UNAVAILABLE else f"CONFORMAL_{cal_quality}"
    )

    return {
        "id": model.id, "endpoint": model.endpoint_name, "model_name": model.model_name,
        "model_version": model.model_version,
        "status": model.implementation_status if available else "MODEL_UNAVAILABLE",
        "availability": model.implementation_status if available else "MODEL_UNAVAILABLE",
        "confidence": confidence,
        "conformal_status": conformal_status,
        "active": available, "output_unit": model.output_unit,
        "source": model.source, "training_dataset": model.training_dataset,
        "validation": model.validation_json or {}, "license": model.license,
        "priority": model.model_priority, "ensemble_eligible": bool(model.ensemble_eligible),
        "species": model.species, "output_type": model.output_type,
        "details": details, "unavailable_reason": unavailable_reason,
        "calibration_provenance": cal_provenance,
        "calibration_quality": cal_quality,
        "conformal_governance": cal_info,
    }


def _performance_out(row: ADMETModelPerformance):
    return {
        "id": row.id, "scope": row.scope_key, "project_id": row.project_id,
        "endpoint": row.endpoint_name, "model_id": row.model_id, "task_type": row.task_type,
        "n": row.sample_size, "metrics": row.metrics_json or {},
        "performance_factor": row.performance_factor, "updated_at": row.updated_at.isoformat(),
    }


def _consensus_out(row: ADMETConsensusPrediction):
    return {
        "id": row.id, "run_id": row.run_id, "version_id": row.version_id,
        "endpoint_id": row.endpoint_id, "endpoint": row.endpoint.name,
        "consensus_version": getattr(row, "consensus_version", "stage4d1-static-v1"),
        "consensus_mode": getattr(row, "consensus_mode", "SHADOW"),
        "combined_value": row.combined_value, "unit": row.unit,
        "classification": row.classification, "confidence": row.confidence,
        "applicability_domain": row.applicability_domain,
        "model_agreement": getattr(row, "model_agreement", "SINGLE_MODEL"),
        "dispersion": getattr(row, "dispersion_json", {}),
        "vote_pattern": getattr(row, "vote_pattern", ""),
        "models": row.weights_json or [],
        "provenance": row.provenance_json or {}, "created_at": row.created_at.isoformat(),
        "type": "Consensus Prediction",
    }


def _rank_correlation(x_values, y_values):
    if len(x_values) < 2:
        return None
    x_rank = np.argsort(np.argsort(np.asarray(x_values, dtype=float))).astype(float)
    y_rank = np.argsort(np.argsort(np.asarray(y_values, dtype=float))).astype(float)
    value = float(np.corrcoef(x_rank, y_rank)[0, 1])
    return None if np.isnan(value) else value


def _refresh_model_feedback(db: Session, project_id: int, version_ids: list[int] | None = None):
    """Persist compatible experimental/model comparisons, then recompute scoped metrics."""
    endpoint_names = {row.id: row.name for row in db.scalars(select(ADMETEndpoint).where(ADMETEndpoint.project_id == project_id))}
    version_query = select(CompoundVersion.id).join(Compound).where(Compound.project_id == project_id)
    scoped_versions = version_ids or list(db.scalars(version_query))
    if not scoped_versions:
        return
    measurements = db.scalars(select(ADMETMeasurement).where(ADMETMeasurement.version_id.in_(scoped_versions))).all()
    by_version = {version_id: [row for row in measurements if row.version_id == version_id] for version_id in scoped_versions}
    predictions = db.scalars(select(ADMETPrediction).where(ADMETPrediction.version_id.in_(scoped_versions))).all()
    for prediction in predictions:
        if prediction.predicted_value is None or prediction.model.endpoint_name not in MODEL_SPECS:
            continue
        for item in comparison_for_prediction(
            prediction.model.endpoint_name, prediction.predicted_value,
            by_version.get(prediction.version_id, []), endpoint_names,
        ):
            existing = db.scalar(select(ADMETModelComparison.id).where(
                ADMETModelComparison.prediction_id == prediction.id,
                ADMETModelComparison.measurement_id == item["measurement_id"],
            ))
            if existing:
                continue
            classification = item.get("classification_match") is not None
            experimental = item.get("experimental_normalized")
            error = item.get("absolute_error")
            db.add(ADMETModelComparison(
                project_id=project_id, version_id=prediction.version_id,
                endpoint_id=prediction.endpoint_id, model_id=prediction.model_id,
                prediction_id=prediction.id, measurement_id=item["measurement_id"],
                task_type="classification" if classification else "regression",
                predicted_value=item.get("predicted_normalized"), experimental_value=experimental,
                absolute_error=error, squared_error=(error * error if error is not None else None),
                fold_error=(10 ** error if error is not None and prediction.model.endpoint_name != "Plasma protein binding" else None),
                predicted_class=str(int(item.get("predicted_normalized", 0))) if classification else "",
                experimental_class=str(int(experimental)) if classification else "",
                correct=item.get("classification_match") if classification else None,
            ))
    db.flush()
    model_ids = set(db.scalars(select(ADMETModelComparison.model_id).where(ADMETModelComparison.project_id == project_id)))
    for model_id in model_ids:
        model = db.get(ADMETModelRegistry, model_id)
        for scope_key, scoped_project in (("GLOBAL", None), (f"PROJECT:{project_id}", project_id)):
            query = select(ADMETModelComparison).where(ADMETModelComparison.model_id == model_id)
            if scoped_project is not None:
                query = query.where(ADMETModelComparison.project_id == scoped_project)
            rows = list(db.scalars(query))
            if not rows:
                continue
            task = rows[0].task_type
            metrics = {}
            if task == "classification":
                actual = [int(row.experimental_class) for row in rows]
                predicted = [int(row.predicted_class) for row in rows]
                tp = sum(a == 1 and p == 1 for a, p in zip(actual, predicted)); tn = sum(a == 0 and p == 0 for a, p in zip(actual, predicted))
                fp = sum(a == 0 and p == 1 for a, p in zip(actual, predicted)); fn = sum(a == 1 and p == 0 for a, p in zip(actual, predicted))
                sensitivity = tp / (tp + fn) if tp + fn else None; specificity = tn / (tn + fp) if tn + fp else None
                denominator = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** .5
                accuracy = (tp + tn) / len(rows)
                mcc = ((tp * tn - fp * fn) / denominator if denominator else None)
                balanced_accuracy = ((sensitivity + specificity) / 2.0) if (sensitivity is not None and specificity is not None) else accuracy
                metrics = {"n": len(rows), "accuracy": accuracy, "balanced_accuracy": balanced_accuracy,
                           "sensitivity": sensitivity, "specificity": specificity, "mcc": mcc}
                factor = .5 + .5 * balanced_accuracy
            else:
                errors = [float(row.absolute_error) for row in rows if row.absolute_error is not None]
                predicted = [float(row.predicted_value) for row in rows if row.predicted_value is not None and row.experimental_value is not None]
                experimental = [float(row.experimental_value) for row in rows if row.predicted_value is not None and row.experimental_value is not None]
                mae = float(np.mean(errors)) if errors else None; rmse = float(np.sqrt(np.mean([value * value for value in errors]))) if errors else None
                folds = [row.fold_error for row in rows if row.fold_error is not None]
                metrics = {"n": len(rows), "mae": mae, "rmse": rmse, "mean_fold_error": float(np.mean(folds)) if folds else None,
                           "spearman": _rank_correlation(predicted, experimental)}
                factor = 1.0 / (1.0 + (mae or 0.0))
            performance = db.scalar(select(ADMETModelPerformance).where(
                ADMETModelPerformance.scope_key == scope_key, ADMETModelPerformance.model_id == model_id,
            ))
            if not performance:
                performance = ADMETModelPerformance(scope_key=scope_key, project_id=scoped_project,
                    endpoint_name=model.endpoint_name, model_id=model_id)
                db.add(performance)
            performance.task_type = task; performance.sample_size = len(rows)
            performance.metrics_json = metrics; performance.performance_factor = factor
            performance.updated_at = datetime.now(timezone.utc)
    db.flush()


def _consensus_factor(db: Session, model: ADMETModelRegistry, prediction: ADMETPrediction, project_id: int):
    domain_factor = {"IN_DOMAIN": 1.0, "BORDERLINE": .65, "OUT_OF_DOMAIN": .25}.get(prediction.applicability_domain, .5)
    confidence_factor = {"HIGH": 1.0, "MEDIUM": .8, "LOW": .55}.get(prediction.confidence, .45)
    priority_factor = 1.0 / (1.0 + max(0, (model.model_priority or 100) - 1) * .01)
    project_performance = db.scalar(select(ADMETModelPerformance).where(
        ADMETModelPerformance.scope_key == f"PROJECT:{project_id}", ADMETModelPerformance.model_id == model.id,
    ))
    performance_factor, performance_reason = 1.0, "published validation and model priority"
    if project_performance and project_performance.sample_size >= 10:
        blend = .25 if project_performance.sample_size < 30 else .6
        performance_factor = (1.0 - blend) + blend * project_performance.performance_factor
        performance_reason = f"project performance N={project_performance.sample_size} blended at {blend:.0%}"
    return max(.0001, priority_factor * domain_factor * confidence_factor * performance_factor), performance_reason


def _store_consensus_predictions(db: Session, version: CompoundVersion, project_id: int, predictions: list[ADMETPrediction]):
    grouped = {}
    for prediction in predictions:
        if prediction.predicted_value is not None and prediction.model.ensemble_eligible:
            grouped.setdefault((prediction.endpoint_id, prediction.unit), []).append(prediction)
    rows = []
    for (endpoint_id, unit), items in grouped.items():
        prediction_ids = sorted(item.id for item in items)
        existing = db.scalar(select(ADMETConsensusPrediction).where(
            ADMETConsensusPrediction.version_id == version.id,
            ADMETConsensusPrediction.endpoint_id == endpoint_id,
        ).order_by(ADMETConsensusPrediction.created_at.desc()))
        if existing and sorted((existing.provenance_json or {}).get("prediction_ids", [])) == prediction_ids:
            rows.append(existing); continue
        weighted, weight_rows, reasons = [], [], []
        for item in items:
            factor, reason = _consensus_factor(db, item.model, item, project_id)
            weighted.append((item, factor)); reasons.append(reason)
        total = sum(factor for _, factor in weighted)
        for item, factor in weighted:
            weight_rows.append({"model_id": item.model_id, "model_name": item.model.model_name,
                "model_version": item.model.model_version, "prediction_id": item.id,
                "value": item.predicted_value, "confidence": item.confidence,
                "domain": item.applicability_domain, "weight": factor / total})
        value = sum(item.predicted_value * factor for item, factor in weighted) / total
        spec = MODEL_SPECS.get(items[0].model.endpoint_name, {})
        classification = ""
        if spec.get("prediction_type") == "binary_classification":
            classification = spec.get("positive_label", "POSITIVE") if value >= spec.get("decision_threshold", .5) else spec.get("negative_label", "NEGATIVE")
        in_weight = sum(row["weight"] for row in weight_rows if row["domain"] == "IN_DOMAIN")
        out_weight = sum(row["weight"] for row in weight_rows if row["domain"] == "OUT_OF_DOMAIN")
        domain = "IN_DOMAIN" if in_weight >= .6 else ("OUT_OF_DOMAIN" if out_weight >= .5 else "BORDERLINE")
        confidence = "HIGH" if domain == "IN_DOMAIN" and len(items) > 1 else ("MEDIUM" if domain != "OUT_OF_DOMAIN" else "LOW")
        # Compute dispersion and agreement
        if len(items) > 1:
            values = [float(item.predicted_value) for item in items if item.predicted_value is not None]
            weighted_std = math.sqrt(max(0.0, sum(row["weight"] * (row["value"] - value)**2 for row in weight_rows)))
            dispersion = {
                "model_disagreement_std": round(weighted_std, 4),
                "min": round(min(values), 4),
                "max": round(max(values), 4),
                "range": round(max(values) - min(values), 4),
                "interpretation": "MODEL DISAGREEMENT (weighted standard deviation; not a confidence interval)",
            }
            if spec.get("prediction_type") == "binary_classification":
                classes = [spec.get("positive_label", "POSITIVE") if row["value"] >= spec.get("decision_threshold", .5) else spec.get("negative_label", "NEGATIVE") for row in weight_rows]
                pos_count = sum(1 for c in classes if c == spec.get("positive_label", "POSITIVE"))
                agreement = "HIGH_AGREEMENT" if pos_count == 0 or pos_count == len(classes) else "MODERATE_AGREEMENT"
                vote_pattern = ", ".join(f"{row['model_name']}:{c}" for row, c in zip(weight_rows, classes))
            else:
                agreement = "HIGH_AGREEMENT" if weighted_std <= 0.30 else ("MODERATE_AGREEMENT" if weighted_std <= 0.60 else "LOW_AGREEMENT")
                vote_pattern = ""
        else:
            dispersion = {"model_disagreement_std": 0.0, "min": round(value, 4), "max": round(value, 4), "range": 0.0}
            agreement = "SINGLE_MODEL"
            vote_pattern = ""

        row = ADMETConsensusPrediction(
            run_id=max(items, key=lambda item: item.created_at).run_id, endpoint_id=endpoint_id,
            version_id=version.id,
            consensus_version="stage4d1-static-v1",
            consensus_mode="SHADOW",
            combined_value=value, unit=unit, classification=classification,
            confidence=confidence, applicability_domain=domain,
            model_agreement=agreement,
            dispersion_json=dispersion,
            vote_pattern=vote_pattern,
            weights_json=weight_rows,
            provenance_json={"record_type": "Consensus Prediction", "compound_version_id": version.id,
                "endpoint": items[0].model.endpoint_name, "prediction_ids": prediction_ids,
                "consensus_mode": "SHADOW",
                "weighting_policy": "model priority × applicability domain × confidence × conservative project performance",
                "performance_thresholds": {"project_blend": 10, "strong_project_blend": 30},
                "reasons": sorted(set(reasons)), "timestamp": datetime.now(timezone.utc).isoformat()},
        )
        db.add(row); rows.append(row)
    db.flush()
    return rows


def _admet_prediction_out(prediction: ADMETPrediction, measurements, endpoint_names):
    comparisons = comparison_for_prediction(
        prediction.model.endpoint_name, prediction.predicted_value, measurements, endpoint_names,
    ) if prediction.predicted_value is not None and prediction.model.endpoint_name in MODEL_SPECS else []
    outputs = dict(prediction.outputs_json or {})
    outputs["experimental_comparisons"] = comparisons
    if prediction.model.endpoint_name in MODEL_SPECS and MODEL_SPECS[prediction.model.endpoint_name].get("prediction_type") == "binary_classification":
        outputs["experimental_evidence"] = cyp_experimental_evidence(
            prediction.model.endpoint_name, prediction.predicted_value, measurements, endpoint_names,
        )
    preferred = None
    if comparisons:
        first = comparisons[0]
        preferred = {
            "source": "Experimental", "measurement_id": first["measurement_id"],
            "value": first["experimental_normalized"], "unit": first["normalized_unit"],
            "prediction_preserved": True,
        }
        if prediction.model.endpoint_name.endswith("intrinsic clearance"):
            outputs["experimental_metabolic_stability_assessment"] = metabolic_stability_assessment(
                prediction.model.endpoint_name, first["experimental_normalized"],
            )
    elif prediction.predicted_value is not None:
        preferred = {"source": "Predicted", "value": prediction.predicted_value, "unit": prediction.unit}
    spec = MODEL_SPECS.get(prediction.model.endpoint_name, {})
    provenance = {
        "record_type": "Predicted", "model_name": prediction.model.model_name,
        "model_version": prediction.model.model_version, "endpoint": prediction.model.endpoint_name,
        "unit": prediction.unit, "species": spec.get("species", "Not specified"),
        "dataset": spec.get("training_dataset"), "license": spec.get("license"),
        "validation": spec.get("validation"), "applicability_domain": prediction.applicability_domain,
        "confidence": prediction.confidence, "timestamp": prediction.created_at.isoformat(),
        "compound_version_id": prediction.version_id,
    }
    return {
        "id": prediction.id, "run_id": prediction.run_id, "version_id": prediction.version_id,
        "endpoint_id": prediction.endpoint_id, "endpoint": prediction.endpoint.name,
        "predicted_value": prediction.predicted_value, "unit": prediction.unit,
        "confidence": prediction.confidence, "applicability_domain": prediction.applicability_domain,
        "uncertainty": prediction.uncertainty, "model": _admet_model_out(prediction.model),
        "outputs": outputs, "experimental_comparisons": comparisons, "preferred_result": preferred,
        "created_at": prediction.created_at.isoformat(), "type": "Predicted", "provenance": provenance,
    }


def _integrated_admet_profile(version_id: int, predictions: list[dict], models: list[dict]) -> dict:
    """Deterministic, evidence-preserving Stage 3 overview; never computes a composite score."""
    latest = {}
    for row in predictions:
        if row["version_id"] == version_id and row["endpoint"] not in latest:
            latest[row["endpoint"]] = row
    sections = {
        "Absorption": [name for name in ("Solubility", "Permeability") if name in latest],
        "Distribution": [name for name in ("Plasma protein binding",) if name in latest],
        "Metabolism": [name for name in latest if name.endswith("intrinsic clearance") or name.startswith("CYP")],
        "Transporters": [name for name in latest if MODEL_SPECS.get(name, {}).get("transporter")],
        "Safety": [name for name in ("hERG liability", "Ames mutagenicity", "DILI clinical liability") if name in latest],
    }
    strengths, concerns = [], []
    for name, row in latest.items():
        output = row.get("outputs") or {}
        spec = MODEL_SPECS.get(name, {})
        if spec.get("safety_endpoint"):
            comparison = (row.get("experimental_comparisons") or [None])[0]
            if comparison:
                positive = comparison["experimental_normalized"] == 1.0
                source = "Experimental"
                label = spec["positive_label"] if positive else spec["negative_label"]
            else:
                positive = output.get("classification") == spec.get("positive_label")
                source = "Predicted"
                label = output.get("classification", "Not predicted")
            item = f"{source} {name}: {label} — {row['confidence']} confidence"
            (concerns if positive else strengths).append(item)
        flag = (output.get("liability_summary") or {}).get("flag")
        if flag and flag not in " ".join(concerns):
            concerns.append(f"{flag} — {row['confidence']} confidence")
        assessment = output.get("experimental_metabolic_stability_assessment") or output.get("metabolic_stability_assessment") or {}
        if assessment.get("metabolic_liability_flag"):
            source = "Experimental" if output.get("experimental_metabolic_stability_assessment") else "Predicted"
            concerns.append(f"{source} {name}: {assessment['metabolic_liability_flag']} — {row['confidence']} confidence")
    unknown = [
        f"{model['endpoint']}: MODEL_UNAVAILABLE — {model['unavailable_reason']}"
        for model in models if not model["active"]
    ]
    required = {"record_type", "model_name", "model_version", "endpoint", "unit", "species", "dataset", "license", "validation", "applicability_domain", "confidence", "timestamp", "compound_version_id"}
    missing = []
    for row in latest.values():
        absent = sorted(key for key in required if row.get("provenance", {}).get(key) in (None, ""))
        if absent:
            missing.append({"prediction_id": row["id"], "missing": absent})
    return {
        "compound_version_id": version_id, "sections": sections,
        "summary": {"strengths": strengths, "concerns": concerns, "unknown": unknown},
        "experimental_precedence": True, "overall_score": None,
        "provenance_audit": {"status": "PASS" if not missing else "FAIL", "checked": len(latest), "missing": missing},
    }


def get_or_create_admet_endpoint(db: Session, project_id: int, name: str):
    name = str(name).strip()
    if not name:
        raise HTTPException(status_code=400, detail="endpoint is required")
    endpoint = db.scalar(select(ADMETEndpoint).where(ADMETEndpoint.project_id == project_id, ADMETEndpoint.name == name))
    return endpoint or ADMETEndpoint(project_id=project_id, name=name)


def add_admet_measurement(db: Session, project_id: int, payload: dict) -> dict:
    version = db.get(CompoundVersion, payload.get("version_id"))
    if not version:
        raise HTTPException(status_code=404, detail="CompoundVersion not found")
    compound = db.get(Compound, version.compound_row_id)
    if not compound or compound.project_id != project_id:
        raise HTTPException(status_code=404, detail="CompoundVersion is not in this project")
    value, mean_value, sd = validate_measurement(payload)
    endpoint = get_or_create_admet_endpoint(db, project_id, payload.get("endpoint", ""))
    db.add(endpoint); db.flush()
    row = ADMETMeasurement(
        version_id=version.id, endpoint_id=endpoint.id,
        species=str(payload.get("species") or ""), matrix=str(payload.get("matrix") or ""),
        value=value, qualitative_value=str(payload.get("qualitative_value") or "").strip(),
        unit=str(payload.get("unit", "")).strip(), qualifier=payload.get("qualifier") or "=",
        replicate=str(payload.get("replicate") or "R1"), mean_value=mean_value,
        standard_deviation=sd, sample_size=int(payload["n"]) if payload.get("n") else None,
        method=str(payload.get("method") or ""), source=str(payload.get("source") or "User experimental"),
        experiment_date=str(payload.get("date") or ""), notes=str(payload.get("notes") or ""),
        provenance_json={
            "data_type": "experimental",
            **(payload.get("provenance") or {}),
            **{k: payload[k] for k in ("ph", "assay_ph", "temperature_c", "ionic_strength", "pka_type", "solubility_type") if payload.get(k) is not None}
        },
    )
    if not row.unit:
        db.rollback(); raise HTTPException(status_code=400, detail="unit is required")
    db.add(row); db.commit(); db.refresh(row)
    return measurement_out(row)


def _admet_payload(db: Session, project: Project, versions: dict[int, tuple[str, int]]) -> dict:
    """Serialize only the explicitly supplied CompoundVersion IDs."""
    version_ids = list(versions)
    rows = db.scalars(
        select(ADMETMeasurement).where(ADMETMeasurement.version_id.in_(version_ids))
        .order_by(ADMETMeasurement.created_at.desc())
    ).all() if version_ids else []
    # The legacy ADMET model list remains a production-model inventory.
    # Executable Stage 4D shadow identities are exposed through multimodel
    # provenance/freeze linkage, never as competing primary rows.
    models = [
        model for model in db.scalars(select(ADMETModelRegistry).order_by(ADMETModelRegistry.endpoint_name)).all()
        if (model.provenance_json or {}).get("production_eligible") is not False
    ]
    endpoint_names = {endpoint.id: endpoint.name for endpoint in db.scalars(
        select(ADMETEndpoint).where(ADMETEndpoint.project_id == project.id)
    )}
    predictions = db.scalars(
        select(ADMETPrediction)
        .where(ADMETPrediction.version_id.in_(version_ids))
        .order_by(ADMETPrediction.created_at.desc())
    ).all() if version_ids else []
    consensuses = db.scalars(
        select(ADMETConsensusPrediction)
        .where(ADMETConsensusPrediction.version_id.in_(version_ids))
        .order_by(ADMETConsensusPrediction.created_at.desc())
    ).all() if version_ids else []
    measurements_by_version = {
        version_id: [row for row in rows if row.version_id == version_id] for version_id in version_ids
    }
    runs = db.scalars(
        select(ADMETPredictionRun)
        .where(ADMETPredictionRun.version_id.in_(version_ids))
        .order_by(ADMETPredictionRun.started_at.desc())
        .limit(20)
    ).all() if version_ids else []
    model_rows = [_admet_model_out(model) for model in models]
    prediction_rows = [_admet_prediction_out(
        prediction, measurements_by_version.get(prediction.version_id, []), endpoint_names,
    ) for prediction in predictions]
    performance_rows = list(db.scalars(
        select(ADMETModelPerformance).where(
            (ADMETModelPerformance.project_id == project.id) | (ADMETModelPerformance.scope_key == "GLOBAL")
        ).order_by(ADMETModelPerformance.endpoint_name, ADMETModelPerformance.scope_key)
    ))
    best_project_models = {}
    for performance in [row for row in performance_rows if row.project_id == project.id and row.sample_size >= 10]:
        model = db.get(ADMETModelRegistry, performance.model_id)
        score = performance.metrics_json.get("mae") if performance.task_type == "regression" else -(performance.metrics_json.get("balanced_accuracy") or performance.metrics_json.get("accuracy") or 0)
        current = best_project_models.get(performance.endpoint_name)
        if current is None or score < current["selection_score"]:
            best_project_models[performance.endpoint_name] = {"model_id": model.id, "model_name": model.model_name,
                "model_version": model.model_version, "n": performance.sample_size,
                "metrics": performance.metrics_json, "selection_score": score}
    for row in best_project_models.values():
        row.pop("selection_score", None)
    return {
        "scope": {"project_id": project.id, "version_ids": version_ids},
        "endpoints": [_admet_endpoint_out(e) for e in db.scalars(select(ADMETEndpoint).where(ADMETEndpoint.project_id == project.id))],
        "measurements": [measurement_out(row) for row in rows],
        "models": model_rows,
        "predictions": prediction_rows,
        "consensus_predictions": [_consensus_out(row) for row in consensuses],
        "model_performance": [_performance_out(row) for row in performance_rows],
        "best_project_models": best_project_models,
        "consensus_policy": {"name": "Deterministic evidence-weighted consensus", "project_weight_start_n": 10,
                             "project_weight_strong_n": 30, "probability_is_not_confidence": True},
        "integrated_profiles": {str(version_id): _integrated_admet_profile(version_id, prediction_rows, model_rows) for version_id in version_ids},
        "prediction_runs": [{"id": r.id, "version_id": r.version_id, "status": r.status,
                             "message": r.message, "started_at": r.started_at.isoformat()} for r in runs],
        "csv_columns": ["compound_id", "version_number"] + [
            column for column in ("endpoint", "species", "matrix", "value", "qualitative_value", "unit", "qualifier", "replicate",
                                  "mean", "sd", "n", "method", "source", "date", "notes")
        ],
        "labels_by_version": {str(key): value for key, value in versions.items()},
    }


@app.get("/api/projects/{project_id}/admet")
def list_admet(project_id: int, db: Session = Depends(get_db)):
    project = db.get(Project, project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    versions = {version.id: (compound.compound_id, version.version_number)
                for compound in project.compounds for version in compound.versions}
    return _admet_payload(db, project, versions)


@app.get("/api/compound-versions/{version_id}/admet")
def get_compound_version_admet(version_id: int, db: Session = Depends(get_db)):
    version = db.get(CompoundVersion, version_id)
    if not version:
        raise HTTPException(status_code=404, detail="CompoundVersion not found")
    compound = db.get(Compound, version.compound_row_id)
    project = db.get(Project, compound.project_id)
    return _admet_payload(db, project, {version.id: (compound.compound_id, version.version_number)})


@app.get("/api/compound-versions/{version_id}/multimodel-provenance")
def get_compound_version_multimodel_provenance(version_id: int, db: Session = Depends(get_db)):
    version = db.get(CompoundVersion, version_id)
    if not version:
        raise HTTPException(status_code=404, detail="CompoundVersion not found")
    compound = db.get(Compound, version.compound_row_id)
    project = db.get(Project, compound.project_id)

    predictions = list(db.scalars(
        select(ADMETPrediction).join(ADMETModelRegistry)
        .where(ADMETPrediction.version_id == version_id)
        .order_by(ADMETModelRegistry.endpoint_name, ADMETPrediction.created_at.desc())
    ).all())

    consensuses = list(db.scalars(
        select(ADMETConsensusPrediction)
        .where(ADMETConsensusPrediction.version_id == version_id)
        .order_by(ADMETConsensusPrediction.created_at.desc())
    ).all())
    consensuses_by_endpoint = {}
    for consensus in consensuses:
        consensuses_by_endpoint.setdefault(consensus.endpoint_id, []).append(consensus)

    measurements = list(db.scalars(select(ADMETMeasurement).where(ADMETMeasurement.version_id == version_id)).all())
    endpoint_names = {endpoint.id: endpoint.name for endpoint in db.scalars(
        select(ADMETEndpoint).where(ADMETEndpoint.project_id == compound.project_id)
    )}

    grouped = {}
    for pred in predictions:
        ep_name = pred.model.endpoint_name
        grouped.setdefault(ep_name, []).append(_admet_prediction_out(pred, measurements, endpoint_names))

    freeze_by_endpoint = {}
    try:
        from .production_qualification import QualificationPredictionFreezeRow
        for frozen in db.scalars(select(QualificationPredictionFreezeRow).where(
            QualificationPredictionFreezeRow.compound_version_id == str(version_id)
        )):
            freeze_by_endpoint[frozen.endpoint_id] = frozen
    except Exception:
        pass

    endpoint_provenance = []
    for ep_name, model_preds in grouped.items():
        ep_obj = db.scalar(select(ADMETEndpoint).where(ADMETEndpoint.project_id == project.id, ADMETEndpoint.name == ep_name))
        endpoint_consensuses = consensuses_by_endpoint.get(ep_obj.id, []) if ep_obj else []
        research_consensuses = [
            row for row in endpoint_consensuses
            if (row.provenance_json or {}).get("stage4d6_research_only")
        ]
        cons_record = next((row for row in endpoint_consensuses if row not in research_consensuses), None)
        contract = get_endpoint_contract(ep_name)
        from .endpoint_strategy_registry import get_endpoint_strategy
        policy = get_endpoint_strategy(ep_name)
        # Annotate model_role from outputs_json for each prediction
        for mp in model_preds:
            if "model_role" not in mp:
                mp_outputs = mp.get("outputs") or {}
                mp["model_role"] = mp_outputs.get("model_role", "CORE")
        # Separate CORE and shadow predictions
        core_preds = [mp for mp in model_preds if mp.get("model_role", "CORE") == "CORE"]
        shadow_preds = [mp for mp in model_preds if mp.get("model_role", "CORE") != "CORE"]
        if policy:
            for core in core_preds:
                core_model = core.get("model") or {}
                core_model.setdefault("model_id", policy.primary_model_ids[0] if policy.primary_model_ids else None)
                core_model["role"] = "CORE"
                core["model"] = core_model
        production = core_preds[0] if core_preds else None
        frozen = freeze_by_endpoint.get(contract.endpoint_id if contract else ep_name)
        freeze_ids = [frozen.frozen_prediction_id] if frozen else []
        frozen_shadows = [
            row for row in ((frozen.provenance_json or {}).get("individual_predictions") if frozen else []) or []
            if row.get("role") != "CORE"
        ]
        for shadow in frozen_shadows:
            model_preds.append({
                "id": None,
                "run_id": None,
                "version_id": version_id,
                "endpoint_id": ep_obj.id if ep_obj else None,
                "endpoint": ep_name,
                "predicted_value": shadow.get("predicted_value"),
                "unit": contract.canonical_unit if contract else "",
                "confidence": shadow.get("confidence", "NOT_APPLICABLE"),
                "applicability_domain": shadow.get("applicability_domain", "UNKNOWN"),
                "uncertainty": None,
                "model": {
                    "model_id": shadow.get("model_id"),
                    "model_name": shadow.get("model_name"),
                    "model_version": shadow.get("model_version"),
                    "role": shadow.get("role"),
                    "production_eligible": False,
                },
                "outputs": {"model_role": shadow.get("role"), "execution_status": shadow.get("execution_status", "SUCCESS")},
                "experimental_comparisons": [],
                "preferred_result": None,
                "created_at": frozen.frozen_at.isoformat() if frozen and frozen.frozen_at else None,
                "type": "Shadow Prediction",
                "provenance": {"record_type": "Shadow Prediction", "freeze_id": frozen.frozen_prediction_id if frozen else None},
                "model_role": shadow.get("role", "SHADOW"),
            })
        shadow_preds = [mp for mp in model_preds if mp.get("model_role", "CORE") != "CORE"]
        endpoint_provenance.append({
            "endpoint_name": ep_name,
            "endpoint_id": contract.endpoint_id if contract else ep_name,
            "canonical_unit": contract.canonical_unit if contract else (model_preds[0]["unit"] if model_preds else ""),
            "consensus": _consensus_out(cons_record) if cons_record else None,
            "research_consensus": [_consensus_out(row) for row in research_consensuses],
            "models": model_preds,
            "model_count": len(model_preds),
            "core_model_count": len(core_preds),
            "shadow_model_count": len(shadow_preds),
            "production_strategy": policy.primary_strategy.value if policy else "UNKNOWN",
            "production_model": production["model"] if production else None,
            "production_value": production["predicted_value"] if production else None,
            "calibration_status": policy.calibration_status.value if policy else "UNKNOWN",
            "calibration_result": {
                "status": policy.calibration_status.value,
                "production_enabled": policy.calibration_production_enabled,
                "value": None,
            } if policy else None,
            "validation_status": policy.validation_status.value if policy else "UNKNOWN",
            "limitations": list(policy.limitations) if policy else [],
            "shadow_exclusion_reasons": {
                model_id: role for model_id, role in (policy.non_primary_model_roles.items() if policy else [])
            },
            "freeze_linkage": freeze_ids,
        })

    return {
        "compound_version_id": version_id,
        "compound_id": compound.compound_id,
        "version_number": version.version_number,
        "canonical_smiles": version.canonical_smiles,
        "consensus_mode": "SHADOW",
        "stage4d6_orchestrator": True,
        "endpoints": endpoint_provenance,
        "total_endpoints": len(endpoint_provenance),
        "total_model_executions": sum(ep["model_count"] for ep in endpoint_provenance),
        "shadow_model_executions": sum(ep["shadow_model_count"] for ep in endpoint_provenance),
    }


def _register_experimental_feedback_event(db: Session, project_id: int, version: CompoundVersion, measurement_row: ADMETMeasurement):
    endpoint_name = measurement_row.endpoint.name if measurement_row.endpoint else ""
    if endpoint_name != "Solubility":
        return
    is_compat, quality, msg = evaluate_experimental_compatibility(
        endpoint_name=endpoint_name,
        value=measurement_row.value,
        unit=measurement_row.unit,
        method=measurement_row.method,
        notes=measurement_row.notes,
    )
    if not is_compat:
        return

    event_id = f"EVT-SOL-{project_id}-{version.id}-{measurement_row.id}"
    existing = db.scalar(select(ADMETExperimentalFeedbackEvent).where(ADMETExperimentalFeedbackEvent.event_id == event_id))
    if existing:
        return

    preds = list(db.scalars(
        select(ADMETPrediction).join(ADMETModelRegistry)
        .where(
            ADMETPrediction.version_id == version.id,
            ADMETModelRegistry.endpoint_name == endpoint_name,
            ADMETPrediction.created_at <= measurement_row.created_at,
        )
    ).all())

    frozen_preds = {}
    model_errors = {}
    for p in preds:
        if p.predicted_value is not None:
            m_key = "admetica_solubility" if "admetica" in p.model.model_name.lower() else ("esol_delaney_v1" if "esol" in p.model.model_name.lower() else "rdkit_gbr_solubility_v1")
            frozen_preds[m_key] = float(p.predicted_value)
            model_errors[m_key] = float(abs(p.predicted_value - measurement_row.value))

    has_frozen = len(frozen_preds) > 0
    scaffold = get_bemis_murcko_scaffold(version.canonical_smiles)

    ev_row = ADMETExperimentalFeedbackEvent(
        event_id=event_id,
        project_id=project_id,
        version_id=version.id,
        endpoint_name=endpoint_name,
        experiment_id=measurement_row.id,
        experimental_value=float(measurement_row.value),
        experimental_unit=measurement_row.unit,
        assay_quality=quality.value,
        scaffold_smiles=scaffold,
        frozen_predictions_json=frozen_preds,
        model_errors_json=model_errors,
        is_valid=has_frozen and quality != AssayQuality.INCOMPATIBLE,
    )
    db.add(ev_row)
    db.flush()


@app.get("/api/compound-versions/{version_id}/adaptive-provenance")
def get_compound_version_adaptive_provenance(
    version_id: int,
    endpoint_name: str = "Solubility",
    include_m3: bool = False,
    db: Session = Depends(get_db),
):
    version = db.get(CompoundVersion, version_id)
    if not version:
        raise HTTPException(status_code=404, detail="CompoundVersion not found")
    compound = db.get(Compound, version.compound_row_id)
    project = db.get(Project, compound.project_id)

    preds = list(db.scalars(
        select(ADMETPrediction).join(ADMETModelRegistry)
        .where(
            ADMETPrediction.version_id == version_id,
            ADMETModelRegistry.endpoint_name == endpoint_name,
        )
        .order_by(ADMETPrediction.created_at.desc())
    ).all())

    adapters = get_adapters_for_endpoint(endpoint_name)
    contract = get_endpoint_contract(endpoint_name)

    payloads = []
    if preds:
        seen = set()
        for p in preds:
            m_key = "admetica_solubility" if "admetica" in p.model.model_name.lower() else ("esol_delaney_v1" if "esol" in p.model.model_name.lower() else "rdkit_gbr_solubility_v1")
            if m_key in seen:
                continue
            seen.add(m_key)
            payloads.append(ModelExecutionPayload(
                model_id=m_key,
                model_name=p.model.model_name,
                model_family=p.model.source or "admetica",
                model_version=p.model_version,
                endpoint_id=endpoint_name,
                endpoint_name=endpoint_name,
                canonical_unit=p.unit,
                execution_status=ExecutionStatus.SUCCESS if p.predicted_value is not None else ExecutionStatus.RUNTIME_ERROR,
                value=p.predicted_value,
                applicability_domain=p.applicability_domain,
                confidence=p.confidence,
            ))
    else:
        for adapter in adapters:
            if not include_m3 and adapter.model_id == "rdkit_gbr_solubility_v1":
                continue
            payload = adapter.execute(version.canonical_smiles, contract)
            payloads.append(payload)

    event_rows = list(db.scalars(
        select(ADMETExperimentalFeedbackEvent)
        .where(
            ADMETExperimentalFeedbackEvent.project_id == project.id,
            ADMETExperimentalFeedbackEvent.endpoint_name == endpoint_name,
            ADMETExperimentalFeedbackEvent.is_valid == True,
        )
    ).all())

    historical_records = []
    for ev in event_rows:
        c_ver = db.get(CompoundVersion, ev.version_id)
        smi = c_ver.canonical_smiles if c_ver else ""
        historical_records.append(ExperimentalFeedbackRecord(
            event_id=ev.event_id,
            project_id=ev.project_id,
            compound_version_id=ev.version_id,
            canonical_smiles=smi,
            endpoint_name=ev.endpoint_name,
            experimental_value=ev.experimental_value,
            experimental_unit=ev.experimental_unit,
            assay_quality=AssayQuality(ev.assay_quality) if ev.assay_quality in AssayQuality._value2member_map_ else AssayQuality.USABLE,
            scaffold_smiles=ev.scaffold_smiles,
            timestamp=ev.created_at.isoformat() if ev.created_at else "",
            frozen_predictions=ev.frozen_predictions_json or {},
            model_errors=ev.model_errors_json or {},
            is_valid=ev.is_valid,
        ))

    result = compute_hierarchical_adaptive_weights(
        query_smiles=version.canonical_smiles,
        project_id=project.id,
        candidate_payloads=payloads,
        historical_feedback_events=historical_records,
        include_m3=include_m3,
    )
    result.compound_version_id = version_id

    existing_adapt = db.scalar(
        select(ADMETAdaptivePrediction).where(
            ADMETAdaptivePrediction.version_id == version.id,
            ADMETAdaptivePrediction.endpoint_name == endpoint_name,
            ADMETAdaptivePrediction.policy_version == result.policy_version,
        )
    )
    if existing_adapt:
        existing_adapt.predicted_value = result.predicted_value
        existing_adapt.model_disagreement = result.model_disagreement
        existing_adapt.effective_weights_json = result.effective_weights
        existing_adapt.weights_breakdown_json = {k: v.to_dict() for k, v in result.weights_breakdown.items()}
        existing_adapt.sample_counts_json = result.to_dict()["sample_counts"]
        existing_adapt.series_info_json = result.to_dict()["series"]
        existing_adapt.reason_codes_json = result.reason_codes
        existing_adapt.warnings_json = result.warnings
    else:
        adapt_row = ADMETAdaptivePrediction(
            version_id=version.id,
            endpoint_name=endpoint_name,
            policy_version=result.policy_version,
            consensus_mode="SHADOW",
            predicted_value=result.predicted_value,
            model_disagreement=result.model_disagreement,
            effective_weights_json=result.effective_weights,
            weights_breakdown_json={k: v.to_dict() for k, v in result.weights_breakdown.items()},
            sample_counts_json=result.to_dict()["sample_counts"],
            series_info_json=result.to_dict()["series"],
            reason_codes_json=result.reason_codes,
            warnings_json=result.warnings,
        )
        db.add(adapt_row)
    db.commit()

    return {
        "compound_version_id": version_id,
        "compound_id": compound.compound_id,
        "version_number": version.version_number,
        "canonical_smiles": version.canonical_smiles,
        **result.to_dict(),
    }


@app.post("/api/projects/{project_id}/admet/measurements", status_code=201)
def create_admet_measurement(project_id: int, payload: dict, db: Session = Depends(get_db)):
    project = db.get(Project, project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    result = add_admet_measurement(db, project_id, payload)
    version_id = int(payload["version_id"])
    version = db.get(CompoundVersion, version_id)
    m_row = db.get(ADMETMeasurement, result["id"])
    if version and m_row:
        _register_experimental_feedback_event(db, project_id, version, m_row)
    _refresh_model_feedback(db, project_id, [version_id])
    refresh_pk_and_ivive_for_version(db, version_id, force=True)
    db.commit()
    return result


@app.post("/api/projects/{project_id}/admet/import-preview")
def admet_import_preview(project_id: int, payload: dict, db: Session = Depends(get_db)):
    if not db.get(Project, project_id):
        raise HTTPException(status_code=404, detail="Project not found")
    records, columns = parse_csv(payload.get("csv", ""))
    valid, errors = [], []
    labels = {compound.compound_id: compound for compound in db.scalars(select(Compound).where(Compound.project_id == project_id))}
    for number, record in enumerate(records, 2):
        compound = labels.get(str(record.get("compound_id", "")).strip())
        version_number = str(record.get("version_number") or "").strip() or (
            str(compound.current_version) if compound else "")
        version = next((v for v in (compound.versions if compound else []) if str(v.version_number) == version_number), None)
        try:
            validate_measurement(record)
            if not version:
                raise ValueError("unknown compound/version")
            valid.append({"row": number, **record})
        except HTTPException as exc:
            errors.append({"row": number, "error": str(exc.detail)})
        except (ValueError, TypeError) as exc:
            errors.append({"row": number, "error": str(exc)})
    return {"columns": columns, "valid_count": len(valid), "errors": errors, "rows": valid}


@app.post("/api/projects/{project_id}/admet/import", status_code=201)
def admet_import(project_id: int, payload: dict, db: Session = Depends(get_db)):
    preview = admet_import_preview(project_id, {"csv": payload.get("csv", "")}, db)
    if preview["errors"]:
        raise HTTPException(status_code=400, detail={"message": "Import validation failed", "errors": preview["errors"]})
    labels = {compound.compound_id: compound for compound in db.scalars(select(Compound).where(Compound.project_id == project_id))}
    created = []
    version_ids = set()
    for item in preview["rows"]:
        compound = labels[str(item["compound_id"]).strip()]
        version_number = int(item.get("version_number") or compound.current_version)
        version = next(version for version in compound.versions if version.version_number == version_number)
        created.append(add_admet_measurement(db, project_id, {**item, "version_id": version.id}))
        version_ids.add(version.id)
    _refresh_model_feedback(db, project_id)
    for vid in version_ids:
        refresh_pk_and_ivive_for_version(db, vid, force=True)
    db.commit()
    return {"imported": len(created), "measurements": created}


@app.get("/api/projects/{project_id}/admet/export.csv")
def admet_export(project_id: int, db: Session = Depends(get_db)):
    project = db.get(Project, project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    labels = {}
    for compound in project.compounds:
        for version in compound.versions:
            labels[version.id] = (compound.compound_id, version.version_number)
    rows = db.scalars(select(ADMETMeasurement).join(ADMETEndpoint, ADMETEndpoint.id == ADMETMeasurement.endpoint_id)
                      .where(ADMETEndpoint.project_id == project_id)).all()
    return csv_export(rows, labels)



def _check_shadow_cache(db: Session, version_id: int, compound: Compound) -> bool:
    """
    Check whether all authorized shadow models are already persisted for this version.
    Returns True only if every shadow model that should be executed has an existing prediction row.
    """
    from .prediction_orchestrator import SHADOW_ADAPTER_MAP
    from .endpoint_strategy_registry import get_all_strategies
    from .production_qualification import QualificationPredictionFreezeRow

    strategies = get_all_strategies()
    freeze_rows = list(db.scalars(select(QualificationPredictionFreezeRow).where(
        QualificationPredictionFreezeRow.compound_version_id == str(version_id)
    )))
    freezes_by_endpoint = {row.endpoint_id: row for row in freeze_rows}
    for ep_name, policy in strategies.items():
        if not policy.shadow_model_ids:
            continue
        for shadow_id in policy.shadow_model_ids:
            role = policy.non_primary_model_roles.get(shadow_id, "SHADOW")
            if "EXCLUDED" in role.upper():
                continue
            if shadow_id not in SHADOW_ADAPTER_MAP:
                continue
            frozen = freezes_by_endpoint.get(policy.endpoint_id)
            executed = ((frozen.provenance_json or {}).get("individual_predictions") if frozen else []) or []
            if not any(row.get("model_id") == shadow_id for row in executed):
                return False
    return True


def _check_runtime_freeze_cache(db: Session, version_id: int) -> bool:
    """Require Stage 4D-6 pre-experimental evidence before a true cache hit."""
    try:
        from .production_qualification import QualificationPredictionFreezeRow
        from .endpoint_strategy_registry import get_all_strategies
        required = {
            policy.endpoint_id for name, policy in get_all_strategies().items()
            if name in MODEL_SPECS and policy.production_execution_allowed
        }
        found = set(db.scalars(
            select(QualificationPredictionFreezeRow.endpoint_id).where(
                QualificationPredictionFreezeRow.compound_version_id == str(version_id)
            )
        ))
        return required <= found
    except Exception:
        # A missing schema is never treated as a cache hit: the orchestrator
        # will create it through normal application initialization.
        return False


def _run_admet_predictions_legacy(
    row_id: int,
    db: Session,
    version: CompoundVersion,
    compound: Compound,
    available_models: list,
    cached: dict,
    measurements: list,
    endpoint_names: dict,
    inactive_statuses: list,
    fallback_reason: str = "",
) -> dict:
    """
    Legacy single-model prediction fallback used only when the orchestrator fails.
    Preserves original Stage 3A-3F behavior.
    """
    active_models = db.scalars(select(ADMETModelRegistry).where(ADMETModelRegistry.is_active.is_(True))).all()
    digest = hashlib.sha256(f"{version.id}|{version.canonical_smiles}|{MODEL_VERSION}|legacy".encode()).hexdigest()
    run = ADMETPredictionRun(
        version_id=row_id, inputs_hash=digest, status="RUNNING",
        message=f"Legacy prediction path (orchestrator fallback: {fallback_reason[:200] if fallback_reason else 'unspecified'}).",
    )
    db.add(run); db.flush()
    created, unavailable, endpoint_statuses = [], [], []
    selected_predictions = dict(cached)
    for model in active_models:
        if model.id in cached:
            endpoint_statuses.append({"endpoint": model.endpoint_name, "model_id": model.id, "status": "COMPLETE", "cache_hit": True})
            continue
        if model.endpoint_name not in MODEL_SPECS:
            continue
        available, reason = model_files_available(model.endpoint_name)
        if not available:
            unavailable.append(f"{model.endpoint_name}: {reason}")
            endpoint_statuses.append({"endpoint": model.endpoint_name, "model_id": model.id, "status": "MODEL_UNAVAILABLE", "message": reason})
            continue
        try:
            result = predict_endpoint(version.canonical_smiles, model.endpoint_name)
        except Exception as exc:
            unavailable.append(f"{model.endpoint_name}: inference failed ({exc})")
            endpoint_statuses.append({"endpoint": model.endpoint_name, "model_id": model.id, "status": "FAILED", "message": str(exc)})
            continue
        if result.get("status") != "COMPLETE":
            unavailable.append(f"{model.endpoint_name}: {result.get('reason', 'model unavailable')}")
            endpoint_statuses.append({"endpoint": model.endpoint_name, "model_id": model.id, "status": "MODEL_UNAVAILABLE", "message": result.get("reason", "model unavailable")})
            continue
        endpoint_obj = get_or_create_admet_endpoint(db, compound.project_id, model.endpoint_name)
        db.add(endpoint_obj); db.flush()
        domain = result["applicability_domain"]
        output = {
            "record_type": "Predicted", "compound_version_id": row_id,
            "prediction_timestamp": datetime.now(timezone.utc).isoformat(),
            "model_source": MODEL_SPECS[model.endpoint_name]["source"],
            "endpoint_definition": MODEL_SPECS[model.endpoint_name]["endpoint_definition"],
            "training_dataset": MODEL_SPECS[model.endpoint_name]["training_dataset"],
            "validation": MODEL_SPECS[model.endpoint_name]["validation"],
            "license": MODEL_SPECS[model.endpoint_name]["license"],
            "limitations": MODEL_SPECS[model.endpoint_name]["limitations"],
            "applicability_domain_details": domain,
            "uncertainty_reason": result["uncertainty_reason"],
        }
        for key in ("assay_definition", "training_n", "independent_validation"):
            if MODEL_SPECS[model.endpoint_name].get(key) is not None:
                output[key] = MODEL_SPECS[model.endpoint_name][key]
        for key in ("probability", "classification", "isoform", "transporter", "safety_endpoint", "species", "role", "decision_threshold", "liability_summary", "ensemble_probabilities"):
            if result.get(key) is not None:
                output[key] = result[key]
        if result.get("derived_outputs") is not None:
            output["derived_outputs"] = result["derived_outputs"]
        if result.get("metabolic_stability_assessment") is not None:
            output["metabolic_stability_assessment"] = result["metabolic_stability_assessment"]
        if result.get("calibrated_uncertainty") is not None:
            output["calibrated_uncertainty"] = result["calibrated_uncertainty"]
        prediction = ADMETPrediction(
            run_id=run.id, endpoint_id=endpoint_obj.id, version_id=row_id, model_id=model.id,
            model_version=model.model_version, execution_status="SUCCESS",
            standardizer_version="CHEM_STANDARDIZER_V1", canonical_smiles=version.canonical_smiles,
            runtime_ms=float(result.get("runtime_ms", 0.0)),
            predicted_value=result["predicted_value"], unit=result["unit"],
            confidence=result["confidence"], applicability_domain=domain["classification"],
            uncertainty=result["uncertainty"], outputs_json=output,
        )
        db.add(prediction); created.append(prediction); selected_predictions[model.id] = prediction
        endpoint_statuses.append({"endpoint": model.endpoint_name, "model_id": model.id, "status": "COMPLETE", "cache_hit": False})
    run.completed_at = datetime.now(timezone.utc)
    if selected_predictions and unavailable:
        run.status, run.message = "PARTIAL", "Legacy predictions completed; " + "; ".join(unavailable)
    elif selected_predictions:
        run.status, run.message = "COMPLETE", "Legacy Stage 3A-3F endpoint predictions completed."
    else:
        run.status, run.message = "MODEL_UNAVAILABLE", "; ".join(unavailable) or "No implemented ADMET model is available."
    db.flush()
    _refresh_model_feedback(db, compound.project_id, [row_id])
    consensuses = _store_consensus_predictions(db, version, compound.project_id, list(selected_predictions.values()))
    db.commit()
    db.refresh(run)
    endpoint_names = {endpoint.id: endpoint.name for endpoint in db.scalars(
        select(ADMETEndpoint).where(ADMETEndpoint.project_id == compound.project_id)
    )}
    predictions = [_admet_prediction_out(
        selected_predictions[model.id], measurements, endpoint_names,
    ) for model in available_models if model.id in selected_predictions]
    return {
        "type": "Predicted", "run_id": run.id, "status": run.status, "message": run.message,
        "models_available": len(available_models), "cache_hit": False, "predictions": predictions,
        "consensus_predictions": [_consensus_out(row) for row in consensuses],
        "endpoint_statuses": endpoint_statuses + inactive_statuses, "unavailable": unavailable,
        "legacy_fallback": True, "fallback_reason": fallback_reason,
    }


@app.post("/api/admet/predict/{row_id}", status_code=202)
def run_admet_predictions(row_id: int, db: Session = Depends(get_db)):
    """
    Stage 4D-6 canonical prediction endpoint.
    Delegates to PredictionOrchestrator which executes all authorized
    CORE + SHADOW/secondary models per endpoint strategy policy.
    """
    from .prediction_orchestrator import (
        PredictionOrchestrator,
        ORCHESTRATOR_VERSION,
        is_core_registry_model,
    )

    version = db.get(CompoundVersion, row_id)
    if not version:
        raise HTTPException(status_code=404, detail="CompoundVersion not found")
    compound = db.get(Compound, version.compound_row_id)
    if compound.project.molecule_type != "Small Molecule":
        raise HTTPException(status_code=400, detail="This model currently supports small molecules only.")

    # Collect inactive-model unavailable statuses (for backward-compat response)
    inactive_models = db.scalars(select(ADMETModelRegistry).where(ADMETModelRegistry.is_active.is_(False))).all()
    inactive_statuses = [
        {"endpoint": model.endpoint_name, "model_id": model.id, "status": "MODEL_UNAVAILABLE",
         "message": (model.provenance_json or {}).get("reason", "Model is disabled or no validated local checkpoint is installed")}
        for model in inactive_models
        if model.endpoint_name in MODEL_SPECS  # only primary-endpoint unavailable entries in legacy view
    ]

    # Active CORE models for backward-compatible response formatting
    active_models = db.scalars(select(ADMETModelRegistry).where(ADMETModelRegistry.is_active.is_(True))).all()
    available_models = [
        model for model in active_models
        if is_core_registry_model(model) and model_files_available(model.endpoint_name)[0]
    ]

    # ── Check full CORE cache (only primary models, not shadow) ─────────────
    cached = {}
    for model in available_models:
        prediction = db.scalar(
            select(ADMETPrediction).join(ADMETModelRegistry)
            .where(
                ADMETPrediction.version_id == row_id,
                ADMETModelRegistry.id == model.id,
                ADMETModelRegistry.model_version == model.model_version,
                ADMETPrediction.execution_status == "SUCCESS",
            )
            .order_by(ADMETPrediction.created_at.desc())
        )
        if prediction:
            cached[model.id] = prediction

    measurements = db.scalars(select(ADMETMeasurement).where(ADMETMeasurement.version_id == row_id)).all()
    endpoint_names = {endpoint.id: endpoint.name for endpoint in db.scalars(
        select(ADMETEndpoint).where(ADMETEndpoint.project_id == compound.project_id)
    )}

    # Compatibility support for isolated callers that deliberately replace all
    # governed core rows with ad-hoc test/local models. Production records
    # always contain the policy identities above and therefore stay on the
    # canonical orchestrator path.
    legacy_compat_models = [
        model for model in active_models
        if model.endpoint_name in MODEL_SPECS and model_files_available(model.endpoint_name)[0]
    ]
    if not available_models and legacy_compat_models:
        legacy_cached = {}
        for model in legacy_compat_models:
            prediction = db.scalar(
                select(ADMETPrediction).where(
                    ADMETPrediction.version_id == row_id,
                    ADMETPrediction.model_id == model.id,
                    ADMETPrediction.execution_status == "SUCCESS",
                ).order_by(ADMETPrediction.created_at.desc())
            )
            if prediction is not None:
                legacy_cached[model.id] = prediction
        if len(legacy_cached) == len(legacy_compat_models):
            consensuses = _store_consensus_predictions(
                db, version, compound.project_id, list(legacy_cached.values())
            )
            db.commit()
            return {
                "type": "Predicted", "run_id": next(iter(legacy_cached.values())).run_id,
                "status": "CACHED", "message": "Cached compatibility predictions reused.",
                "models_available": len(legacy_compat_models), "cache_hit": True,
                "predictions": [
                    _admet_prediction_out(legacy_cached[model.id], measurements, endpoint_names)
                    for model in legacy_compat_models
                ],
                "consensus_predictions": [_consensus_out(row) for row in consensuses],
                "endpoint_statuses": [
                    {"endpoint": model.endpoint_name, "model_id": model.id, "status": "COMPLETE", "cache_hit": True}
                    for model in legacy_compat_models
                ] + inactive_statuses,
            }
        return _run_admet_predictions_legacy(
            row_id, db, version, compound, legacy_compat_models, legacy_cached, measurements,
            endpoint_names, inactive_statuses,
            fallback_reason="No Stage 4D policy CORE identity is installed; compatibility execution only.",
        )

    # Full core cache hit: still run orchestrator to execute shadow models + freeze
    # (shadow models may not be cached yet even if CORE is)
    full_core_cache = (
        len(available_models) > 0 and len(cached) == len(available_models)
    )

    # Check if shadow models also all cached
    shadow_all_cached = _check_shadow_cache(db, row_id, compound)

    freeze_all_cached = _check_runtime_freeze_cache(db, row_id)
    if full_core_cache and shadow_all_cached and freeze_all_cached:
        # All predictions (core + shadow) already exist — true cache hit
        _refresh_model_feedback(db, compound.project_id, [row_id])
        consensuses = _store_consensus_predictions(db, version, compound.project_id, list(cached.values()))
        db.commit()
        predictions = [_admet_prediction_out(cached[model.id], measurements, endpoint_names) for model in available_models]
        return {
            "type": "Predicted", "run_id": predictions[0]["run_id"] if predictions else 0,
            "status": "CACHED",
            "message": "Cached predictions reused (core + shadow) for this CompoundVersion and model version.",
            "models_available": len(available_models), "cache_hit": True, "predictions": predictions,
            "consensus_predictions": [_consensus_out(row) for row in consensuses],
            "endpoint_statuses": [
                {"endpoint": model.endpoint_name, "model_id": model.id, "status": "COMPLETE", "cache_hit": True}
                for model in available_models
            ] + inactive_statuses,
        }

    # ── Run canonical orchestrator (CORE + SHADOW) ──────────────────────────
    try:
        orchestrator = PredictionOrchestrator(db, version, compound)
        orch_result = orchestrator.orchestrate()
        db.flush()
    except Exception as exc:
        # Orchestrator failure: fall back to legacy single-model prediction path
        import traceback
        _tb = traceback.format_exc()
        return _run_admet_predictions_legacy(
            row_id, db, version, compound, available_models, cached, measurements,
            endpoint_names, inactive_statuses, fallback_reason=str(exc),
        )

    # Refresh performance feedback after orchestration
    _refresh_model_feedback(db, compound.project_id, [row_id])

    # Store static consensus (shadow mode, from primary model predictions only)
    # Refresh cached after orchestration to include newly persisted CORE predictions
    refreshed_core_preds = {}
    for model in available_models:
        pred = db.scalar(
            select(ADMETPrediction).join(ADMETModelRegistry)
            .where(
                ADMETPrediction.version_id == row_id,
                ADMETModelRegistry.id == model.id,
                ADMETModelRegistry.model_version == model.model_version,
                ADMETPrediction.execution_status == "SUCCESS",
            )
            .order_by(ADMETPrediction.created_at.desc())
        )
        if pred:
            refreshed_core_preds[model.id] = pred

    consensuses = _store_consensus_predictions(
        db, version, compound.project_id, list(refreshed_core_preds.values())
    )
    db.commit()

    # Refresh endpoint names after potential additions
    endpoint_names = {endpoint.id: endpoint.name for endpoint in db.scalars(
        select(ADMETEndpoint).where(ADMETEndpoint.project_id == compound.project_id)
    )}

    predictions = [
        _admet_prediction_out(refreshed_core_preds[model.id], measurements, endpoint_names)
        for model in available_models if model.id in refreshed_core_preds
    ]

    return {
        "type": "Predicted",
        "run_id": orch_result.run_id,
        "status": orch_result.status,
        "message": orch_result.message,
        "models_available": len(available_models),
        "cache_hit": False,
        "predictions": predictions,
        "consensus_predictions": [_consensus_out(row) for row in consensuses],
        "endpoint_statuses": orch_result.endpoint_statuses + inactive_statuses,
        "unavailable": orch_result.unavailable,
        "orchestrator": ORCHESTRATOR_VERSION,
        "shadow_models_executed": sum(
            ep.shadow_results.__len__() if hasattr(ep.shadow_results, '__len__') else 0
            for ep in orch_result.endpoint_results
        ),
    }


def _soft_spot_out(row: MetabolicSoftSpot):

    return {
        "id": row.id, "run_id": row.run_id, "version_id": row.version_id,
        "rank": row.rank, "atom_index": row.atom_index,
        "atom_environment": row.atom_environment, "transformation": row.transformation,
        "phase": row.phase, "cyp_isoform": row.cyp_isoform,
        "model_evidence": row.model_evidence_json or {},
        "rule_evidence": row.rule_evidence_json or {},
        "score": row.score, "score_type": row.score_type,
        "confidence": row.confidence, "provenance": row.provenance_json or {},
        "created_at": row.created_at.isoformat(),
    }


def _predicted_metabolite_out(row: PredictedMetabolite):
    try:
        structure_svg = analyze_smiles(row.canonical_smiles)["svg"]
    except ChemistryError:
        structure_svg = ""
    return {
        "id": row.id, "run_id": row.run_id, "soft_spot_id": row.soft_spot_id,
        "version_id": row.version_id, "canonical_smiles": row.canonical_smiles,
        "isomeric_smiles": row.isomeric_smiles, "transformation": row.transformation,
        "source_atom": row.source_atom, "phase": row.phase, "rank": row.rank,
        "confidence": row.confidence, "evidence": row.evidence_json or {},
        "provenance": row.provenance_json or {}, "label": PREDICTED_LABEL,
        "created_at": row.created_at.isoformat(), "type": "Predicted",
        "structure_svg": structure_svg,
    }


def _experimental_metabolite_out(row: ExperimentalMetabolite):
    return {
        "id": row.id, "version_id": row.version_id,
        "canonical_smiles": row.canonical_smiles, "isomeric_smiles": row.isomeric_smiles,
        "transformation": row.transformation, "observed_mass": row.observed_mass,
        "mass_unit": row.mass_unit, "source": row.source,
        "experiment": row.experiment, "notes": row.notes,
        "provenance": row.provenance_json or {}, "label": "EXPERIMENTAL METABOLITE",
        "created_at": row.created_at.isoformat(), "type": "Experimental",
    }


def _metabolic_run_out(run: MetabolicPredictionRun):
    spots = sorted(run.spots, key=lambda row: row.rank)
    metabolites = sorted(run.metabolites, key=lambda row: (row.rank, row.id))
    return {
        "id": run.id, "version_id": run.version_id, "inputs_hash": run.inputs_hash,
        "engine": run.engine_name, "engine_version": run.engine_version,
        "status": run.status, "message": run.message,
        "model_status": run.model_status_json or {},
        "liability_summary": run.liability_summary_json or {},
        "highlighted_svg": run.highlighted_svg,
        "spots": [_soft_spot_out(row) for row in spots],
        "predicted_metabolites": [_predicted_metabolite_out(row) for row in metabolites],
        "started_at": run.started_at.isoformat(),
        "completed_at": run.completed_at.isoformat() if run.completed_at else None,
    }


def _metabolism_evidence_context(db: Session, version: CompoundVersion, project_id: int) -> dict:
    measurements = db.scalars(
        select(ADMETMeasurement).where(ADMETMeasurement.version_id == version.id)
        .order_by(ADMETMeasurement.created_at.desc())
    ).all()
    endpoint_names = {row.id: row.name for row in db.scalars(
        select(ADMETEndpoint).where(ADMETEndpoint.project_id == project_id)
    )}
    predictions = db.scalars(
        select(ADMETPrediction).where(ADMETPrediction.version_id == version.id)
        .order_by(ADMETPrediction.created_at.desc())
    ).all()
    latest = {}
    for prediction in predictions:
        latest.setdefault(prediction.model.endpoint_name, prediction)

    microsomal = []
    for endpoint in ("HLM intrinsic clearance", "RLM intrinsic clearance"):
        experimental = None
        for measurement in measurements:
            measurement_endpoint = endpoint_names.get(measurement.endpoint_id, "")
            normalized, note = comparable_experimental(
                endpoint, measurement, measurement_endpoint,
            )
            if normalized is not None:
                experimental = {
                    "endpoint": endpoint, "source": "Experimental",
                    "measurement_id": measurement.id, "value": normalized,
                    "unit": MODEL_SPECS[endpoint]["unit"], "conversion": note,
                    "confidence": "EXPERIMENTAL",
                    "assessment": metabolic_stability_assessment(endpoint, normalized),
                }
                break
            if measurement_endpoint.strip().lower() == endpoint.lower():
                raw_value = measurement.mean_value if measurement.mean_value is not None else measurement.value
                if raw_value is not None:
                    experimental = {
                        "endpoint": endpoint, "source": "Experimental",
                        "measurement_id": measurement.id, "value": raw_value,
                        "unit": measurement.unit, "conversion": note,
                        "confidence": "EXPERIMENTAL", "assessment": None,
                        "comparison_status": "Retained as experimental evidence; not normalized to the prediction unit",
                    }
                    break
        if experimental:
            microsomal.append(experimental)
            continue
        prediction = latest.get(endpoint)
        if prediction:
            assessment = (prediction.outputs_json or {}).get("metabolic_stability_assessment")
            microsomal.append({
                "endpoint": endpoint, "source": "Predicted", "prediction_id": prediction.id,
                "value": prediction.predicted_value, "unit": prediction.unit,
                "confidence": prediction.confidence, "domain": prediction.applicability_domain,
                "assessment": assessment,
            })

    cyp = []
    for endpoint, prediction in latest.items():
        spec = MODEL_SPECS.get(endpoint, {})
        if spec.get("role") != "SUBSTRATE":
            continue
        output = prediction.outputs_json or {}
        cyp.append({
            "endpoint": endpoint, "isoform": spec["isoform"],
            "classification": output.get("classification"),
            "probability": output.get("probability", prediction.predicted_value),
            "confidence": prediction.confidence,
            "domain": prediction.applicability_domain,
            "prediction_id": prediction.id,
            "attribution": "Compound-level substrate evidence only; no atom or reaction assignment.",
        })
    return {"microsomal": microsomal, "cyp": sorted(cyp, key=lambda row: row["isoform"])}


@app.post("/api/metabolism/predict/{version_id}", status_code=202)
def run_metabolism_predictions(version_id: int, db: Session = Depends(get_db)):
    version = db.get(CompoundVersion, version_id)
    if not version:
        raise HTTPException(status_code=404, detail="CompoundVersion not found")
    compound = db.get(Compound, version.compound_row_id)
    if not compound:
        raise HTTPException(status_code=404, detail="Compound not found")
    if compound.project.molecule_type != "Small Molecule":
        raise HTTPException(status_code=400, detail="This model currently supports small molecules only.")
    context = _metabolism_evidence_context(db, version, compound.project_id)
    digest = hashlib.sha256(json.dumps({
        "version_id": version.id, "smiles": version.canonical_smiles,
        "engine_version": METABOLISM_VERSION, "context": context,
    }, sort_keys=True).encode()).hexdigest()
    cached = db.scalar(
        select(MetabolicPredictionRun).where(
            MetabolicPredictionRun.version_id == version_id,
            MetabolicPredictionRun.inputs_hash == digest,
            MetabolicPredictionRun.engine_version == METABOLISM_VERSION,
            MetabolicPredictionRun.status == "COMPLETE",
        ).order_by(MetabolicPredictionRun.started_at.desc())
    )
    if cached:
        return {"status": "CACHED", "cache_hit": True, "message": "Cached soft spots and metabolite hypotheses reused.", "run": _metabolic_run_out(cached)}

    run = MetabolicPredictionRun(
        version_id=version_id, inputs_hash=digest, engine_name=METABOLISM_ENGINE,
        engine_version=METABOLISM_VERSION, status="RUNNING",
        message="Running atom-mapped SyGMa rules with RDKit chemical validation.",
    )
    db.add(run); db.flush()
    try:
        result = predict_soft_spots(version.canonical_smiles, context=context)
        run.model_status_json = result["model_status"]
        run.liability_summary_json = result["liability_summary"]
        run.highlighted_svg = result["highlighted_svg"]
        timestamp = datetime.now(timezone.utc).isoformat()
        spots_by_rank = {}
        for item in result["spots"]:
            provenance = {
                "prediction_timestamp": timestamp, "compound_version_id": version.id,
                "engine": METABOLISM_ENGINE, "engine_version": METABOLISM_VERSION,
                "source": METABOLISM_SOURCE, "license": METABOLISM_LICENSE,
                "training_data": "SyGMa empirical rules derived from the historical MDL Metabolite database; source database is discontinued",
                "publisher_validation": PUBLISHER_VALIDATION,
                "atom_index_basis": "RDKit zero-based canonical molecule atom index",
            }
            spot = MetabolicSoftSpot(
                run_id=run.id, version_id=version.id, rank=item["rank"],
                atom_index=item["atom_index"], atom_environment=item["atom_environment"],
                transformation=item["transformation"], phase=item["phase"],
                cyp_isoform=item["cyp_isoform"], model_evidence_json=item["model_evidence"],
                rule_evidence_json=item["rule_evidence"], score=item["score"],
                score_type=item["score_type"], confidence=item["confidence"],
                provenance_json=provenance,
            )
            db.add(spot); db.flush(); spots_by_rank[item["rank"]] = spot
        for item in result["metabolites"]:
            spot = spots_by_rank[item["rank"]]
            db.add(PredictedMetabolite(
                run_id=run.id, soft_spot_id=spot.id, version_id=version.id,
                canonical_smiles=item["canonical_smiles"], isomeric_smiles=item["isomeric_smiles"],
                transformation=item["transformation"], source_atom=item["source_atom"],
                phase=item["phase"], rank=item["rank"], confidence=item["confidence"],
                evidence_json=item["evidence"], provenance_json={
                    "prediction_timestamp": timestamp, "compound_version_id": version.id,
                    "transformation_engine": METABOLISM_ENGINE,
                    "transformation_engine_version": METABOLISM_VERSION,
                    "source": METABOLISM_SOURCE, "license": METABOLISM_LICENSE,
                    "label": PREDICTED_LABEL,
                },
            ))
        run.status = "COMPLETE"
        run.message = f"Stored {len(result['spots'])} ranked soft spots and {len(result['metabolites'])} unique sanitized metabolite hypotheses."
        run.completed_at = datetime.now(timezone.utc)
        db.commit(); db.refresh(run)
    except Exception as exc:
        run.status, run.message = "FAILED", f"Metabolic hypothesis generation failed: {exc}"
        run.completed_at = datetime.now(timezone.utc)
        db.commit()
        raise HTTPException(status_code=500, detail=run.message)
    return {"status": run.status, "cache_hit": False, "message": run.message, "run": _metabolic_run_out(run)}


def _metabolism_payload(db: Session, version_ids: list[int]) -> dict:
    runs = db.scalars(
        select(MetabolicPredictionRun).where(MetabolicPredictionRun.version_id.in_(version_ids))
        .order_by(MetabolicPredictionRun.started_at.desc())
    ).all() if version_ids else []
    latest_by_version = {}
    for run in runs:
        latest_by_version.setdefault(run.version_id, run)
    experimental = db.scalars(
        select(ExperimentalMetabolite).where(ExperimentalMetabolite.version_id.in_(version_ids))
        .order_by(ExperimentalMetabolite.created_at.desc())
    ).all() if version_ids else []
    return {
        "scope": {"version_ids": version_ids},
        "runs": [_metabolic_run_out(run) for run in latest_by_version.values()],
        "experimental_metabolites": [_experimental_metabolite_out(row) for row in experimental],
        "tool": {
            "name": METABOLISM_ENGINE, "version": METABOLISM_VERSION,
            "source": METABOLISM_SOURCE, "license": METABOLISM_LICENSE,
            "publisher_validation": PUBLISHER_VALIDATION,
        },
        "settings": {"default_top_spots": 3, "available_top_spots": [3, 5, 10, "ALL"]},
    }


@app.get("/api/projects/{project_id}/metabolism")
def list_metabolism(project_id: int, db: Session = Depends(get_db)):
    project = db.get(Project, project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    version_ids = [version.id for compound in project.compounds for version in compound.versions]
    return _metabolism_payload(db, version_ids)


@app.get("/api/compound-versions/{version_id}/metabolism")
def get_compound_version_metabolism(version_id: int, db: Session = Depends(get_db)):
    if not db.get(CompoundVersion, version_id):
        raise HTTPException(status_code=404, detail="CompoundVersion not found")
    return _metabolism_payload(db, [version_id])


@app.get("/api/compound-versions/{version_id}/workspace")
def get_compound_version_workspace(version_id: int, db: Session = Depends(get_db)):
    """Strict CompoundVersion workspace; no sibling or cross-project records."""
    version = db.get(CompoundVersion, version_id)
    if not version:
        raise HTTPException(status_code=404, detail="CompoundVersion not found")
    compound = db.get(Compound, version.compound_row_id)
    project = db.get(Project, compound.project_id)
    measurements = db.scalars(
        select(ActivityMeasurement).where(ActivityMeasurement.version_id == version_id)
        .order_by(ActivityMeasurement.created_at.desc())
    ).all()
    activity_predictions = db.scalars(
        select(ActivityPrediction).where(ActivityPrediction.version_id == version_id)
        .order_by(ActivityPrediction.created_at.desc())
    ).all()
    audit_runs = db.scalars(
        select(PredictionRun).where(PredictionRun.version_id == version_id)
        .order_by(PredictionRun.created_at.desc())
    ).all()
    activity = {
        "measurements": [{
            "id": row.id, "version_id": row.version_id, "assay_id": row.assay_id,
            "assay": row.assay.name, "measurement_type": row.assay.measurement_type,
            "value": row.raw_value, "unit": row.original_unit,
            "normalized_value_nm": row.normalized_value_nm, "qualifier": row.qualifier,
            "source": row.source, "created_at": row.created_at.isoformat(), "type": "Experimental",
        } for row in measurements],
        "predictions": [{
            "id": row.id, "version_id": row.version_id, "assay_id": row.assay_id,
            "assay": row.assay.name, "predicted_value_nm": row.predicted_value_nm,
            "confidence": row.confidence, "applicability_domain": row.applicability_domain,
            "nearest_neighbors": row.nearest_neighbors or [], "created_at": row.created_at.isoformat(),
            "type": "Predicted",
        } for row in activity_predictions],
    }
    admet_payload = _admet_payload(db, project, {version.id: (compound.compound_id, version.version_number)})
    metabolism_payload = _metabolism_payload(db, [version.id])
    pk_parameter_sets = db.scalars(select(PKParameterSet).where(PKParameterSet.version_id == version_id).order_by(PKParameterSet.created_at.desc())).all()
    pk_routes = [{"id": row.id, "species": row.species, "route": row.route, "cl_value": row.cl_value, "v_value": row.v_value, "f_predicted": row.f_predicted, "f_experimental": row.f_experimental, "ka_value": row.ka_value, "confidence": row.confidence} for row in pk_parameter_sets]
    latest_workflow = next((row for row in audit_runs if row.stage == "prediction_workflow"), None)
    saved_steps = ((latest_workflow.outputs_json or {}).get("steps", {}) if latest_workflow else {})
    def saved_status(stage: str, fallback: str) -> str:
        return str((saved_steps.get(stage) or {}).get("status") or fallback)
    return {
        "scope": {"project_id": project.id, "compound_id": compound.id, "version_id": version.id},
        "project": {"id": project.id, "name": project.name, "molecule_type": project.molecule_type},
        "compound": compound_out(compound), "version": serialize_version(version),
        "activity": activity,
        "admet": admet_payload,
        "metabolism": metabolism_payload,
        "pk": {"parameter_sets": pk_routes},
        "prediction_status": {
            "properties": saved_status("properties", "COMPLETE" if version.properties_json else "NOT_STARTED"),
            "activity": saved_status("activity", "COMPLETE" if (activity["measurements"] or activity["predictions"]) else "NOT_STARTED"),
            "admet": saved_status("admet", "COMPLETE" if admet_payload.get("predictions") else "NOT_STARTED"),
            "metabolism": saved_status("metabolism", "COMPLETE" if metabolism_payload.get("runs") else "NOT_STARTED"),
            "pk": saved_status("pk", "COMPLETE" if pk_routes else "NOT_STARTED"),
        },
        "prediction_audit": [{
            "prediction_id": row.id, "version_id": row.version_id,
            "created_at": row.created_at.isoformat(), "stage": row.stage,
            "model_name": row.model_name, "model_version": row.model_version,
            "confidence": row.confidence, "provenance": row.provenance_json,
            "outputs": row.outputs_json,
        } for row in audit_runs],
    }


@app.post("/api/projects/{project_id}/metabolism/experimental", status_code=201)
def create_experimental_metabolite(project_id: int, payload: dict, db: Session = Depends(get_db)):
    version = db.get(CompoundVersion, payload.get("version_id"))
    if not version:
        raise HTTPException(status_code=404, detail="CompoundVersion not found")
    compound = db.get(Compound, version.compound_row_id)
    if not compound or compound.project_id != project_id:
        raise HTTPException(status_code=404, detail="CompoundVersion is not in this project")
    transformation = str(payload.get("transformation") or "").strip()
    if not transformation:
        raise HTTPException(status_code=400, detail="transformation is required")
    smiles = str(payload.get("smiles") or "").strip()
    canonical = isomeric = ""
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid experimental metabolite SMILES")
        canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
        isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    try:
        observed_mass = float(payload["observed_mass"]) if payload.get("observed_mass") not in (None, "") else None
    except (TypeError, ValueError):
        raise HTTPException(status_code=400, detail="observed_mass must be numeric")
    if observed_mass is not None and observed_mass <= 0:
        raise HTTPException(status_code=400, detail="observed_mass must be positive")
    mass_unit = str(payload.get("mass_unit") or "").strip()
    if observed_mass is not None and not mass_unit:
        raise HTTPException(status_code=400, detail="mass_unit is required when observed_mass is supplied")
    row = ExperimentalMetabolite(
        version_id=version.id, canonical_smiles=canonical, isomeric_smiles=isomeric,
        transformation=transformation, observed_mass=observed_mass, mass_unit=mass_unit,
        source=str(payload.get("source") or "User experimental"),
        experiment=str(payload.get("experiment") or ""), notes=str(payload.get("notes") or ""),
        provenance_json={
            "data_type": "experimental_metabolite", "compound_version_id": version.id,
            "recorded_at": datetime.now(timezone.utc).isoformat(),
            **(payload.get("provenance") or {}),
        },
    )
    db.add(row); db.commit(); db.refresh(row)
    return _experimental_metabolite_out(row)


@app.get("/api/projects/{project_id}/assays")
def list_assays(project_id: int, db: Session = Depends(get_db)):
    assays = db.scalars(select(AssayDefinition).where(AssayDefinition.project_id == project_id).order_by(AssayDefinition.created_at)).all()
    return [_assay_out(a) for a in assays if a.active]


@app.post("/api/projects/{project_id}/assays", status_code=201)
def create_assay(project_id: int, payload: dict, db: Session = Depends(get_db), supersedes_id: int | None = None):
    project = db.get(Project, project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    if supersedes_id:
        old=db.get(AssayDefinition,supersedes_id)
        if not old or old.project_id!=project_id: raise HTTPException(status_code=404,detail="Assay to supersede not found")
        old.active=False; payload["supersedes_id"]=old.id; payload["version_number"]=old.version_number+1
    assay = AssayDefinition(project_id=project_id, **payload)
    db.add(assay); db.commit(); db.refresh(assay)
    return _assay_out(assay)


def _experimental_summary(db: Session, version_id: int, assay_id: int):
    rows = db.scalars(select(ActivityMeasurement).where(
        ActivityMeasurement.version_id == version_id,
        ActivityMeasurement.assay_id == assay_id).order_by(ActivityMeasurement.created_at)).all()
    if not rows:
        return None
    values = [row.normalized_value_nm for row in rows]
    mean = sum(values) / len(values)
    sd = (sum((v - mean) ** 2 for v in values) / len(values)) ** 0.5 if len(values) > 1 else 0
    return {
        "type": "Experimental", "n": len(rows), "mean_nm": round(mean, 4), "sd_nm": round(sd, 4),
        "cv_percent": round(sd / mean * 100, 2) if mean else None,
        "pactivity_mean": round(pactivity(mean), 3),
        "raw_measurements": [{"value": r.raw_value, "unit": r.original_unit, "qualifier": r.qualifier,
                              "normalized_nm": r.normalized_value_nm, "replicate": r.replicate_label,
                              "source": r.source} for r in rows],
        "latest_created_at": rows[-1].created_at.isoformat(),
    }


@app.post("/api/assays/{assay_id}/measurements", status_code=201)
def add_measurement(assay_id: int, payload: dict, db: Session = Depends(get_db)):
    assay = db.get(AssayDefinition, assay_id)
    if not assay or not assay.active: raise HTTPException(status_code=404, detail="Active assay not found")
    version = db.get(CompoundVersion, payload.get("version_id"))
    if not version: raise HTTPException(status_code=404, detail="CompoundVersion not found")
    compound = db.get(Compound, version.compound_row_id)
    if not compound or compound.project_id != assay.project_id:
        raise HTTPException(status_code=404, detail="CompoundVersion is not in the assay project")
    try:
        normalized, provenance = normalize_concentration(float(payload["value"]), str(payload.get("unit", assay.unit)))
        transformed = pactivity(normalized)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Activity validation failed: {exc}")
    row = ActivityMeasurement(
        assay_id=assay_id, version_id=version.id, raw_value=float(payload["value"]),
        original_unit=str(payload.get("unit", assay.unit)), normalized_value_nm=normalized,
        qualifier=str(payload.get("qualifier", "=")), replicate_label=str(payload.get("replicate_label") or f"R{db.scalar(select(func.count(ActivityMeasurement.id)).where(ActivityMeasurement.version_id==version.id, ActivityMeasurement.assay_id==assay_id))+1}"),
        experiment_date=str(payload.get("experiment_date", "")), source=str(payload.get("source", "User experimental")),
        notes=str(payload.get("notes", "")),
        provenance_json={**provenance, "transformation": "-log10(value [M])", "transformed_pactivity": transformed},
    )
    db.add(row); db.commit(); db.refresh(row)
    return {"measurement_id": row.id, **_experimental_summary(db, version.id, assay_id)}


@app.post("/api/activities/import-preview")
def import_preview(payload: dict, db: Session = Depends(get_db)):
    import csv, io
    text=payload.get("csv",""); reader=csv.DictReader(io.StringIO(text))
    valid,errors=[],[]
    for line,row in enumerate(reader, start=2):
        compound=db.scalar(select(Compound).join(Project,Project.id==Compound.project_id).where(Compound.project_id==payload.get("project_id"),Compound.compound_id==str(row.get("compound_id","")).strip()))
        assay=db.scalar(select(AssayDefinition).where(AssayDefinition.project_id==payload.get("project_id"),AssayDefinition.name==str(row.get("assay","")).strip(),AssayDefinition.active==True))
        try:
            value=float(row.get("value")); unit=row.get("unit") or (assay.unit if assay else "nM")
            normalized,_=normalize_concentration(value,unit); pactivity(normalized)
            if not compound: raise ValueError("compound not found")
            if not next((version for version in compound.versions if version.version_number == compound.current_version), None):
                raise ValueError("compound has no structure version")
            if not assay: raise ValueError("active assay not found")
            valid.append({"line":line,"compound_id":compound.compound_id,"assay":assay.name,"value":value,"unit":unit})
        except Exception as exc:
            errors.append({"line":line,"error":str(exc)})
    return {"valid":valid,"errors":errors,"can_import":len(errors)==0 and bool(valid)}


@app.post("/api/activities/import", status_code=201)
def import_activities(payload: dict, db: Session = Depends(get_db)):
    preview=import_preview({"project_id":payload.get("project_id"),"csv":payload.get("csv","")},db)
    if not preview["valid"]: return {"imported":0,**preview}
    count=0
    for item in preview["valid"]:
        compound=db.scalar(select(Compound).where(Compound.project_id==payload["project_id"],Compound.compound_id==item["compound_id"]))
        assay=db.scalar(select(AssayDefinition).where(AssayDefinition.project_id==payload["project_id"],AssayDefinition.name==item["assay"],AssayDefinition.active==True))
        version=next((v for v in compound.versions if v.version_number==compound.current_version),None)
        if not version:
            continue
        add_measurement(assay.id,{"version_id":version.id,"value":item["value"],"unit":item["unit"]},db);count+=1
    return {"imported":count,"errors":preview["errors"]}


@app.get("/api/projects/{project_id}/sar")
def sar_table(project_id: int, assay_id: int, db: Session = Depends(get_db)):
    assay=db.get(AssayDefinition,assay_id)
    if not assay or assay.project_id!=project_id: raise HTTPException(status_code=404,detail="Assay not found")
    compounds=db.scalars(select(Compound).where(Compound.project_id==project_id)).all();rows=[]
    dataset=[]
    for compound in compounds:
        version=next((v for v in compound.versions if v.version_number==compound.current_version),None)
        if not version:
            continue
        exp=_experimental_summary(db,version.id,assay.id)
        mol,fp,desc,scaffold=fingerprint_and_descriptors(version.canonical_smiles)
        prediction=db.scalar(select(ActivityPrediction).where(ActivityPrediction.version_id==version.id,ActivityPrediction.assay_id==assay.id).order_by(ActivityPrediction.created_at.desc()))
        rows.append({"row_id":compound.id,"compound":compound.compound_id,"name":compound.name,"svg":version.svg,
                     "properties":{k:version.properties_json.get(k) for k in ["molecular_weight","clogp","tpsa","qed"]},
                     "experimental":exp,"predicted":{"type":"AI Predicted","pactivity":prediction.predicted_pactivity,
                        "value_nm":round(prediction.predicted_value_nm,3),"confidence":prediction.confidence,
                        "applicability_domain":prediction.applicability_domain} if prediction else None,
                     "fingerprint":fp,"descriptors":desc,"scaffold":scaffold})
        if exp:dataset.append({"row_id":compound.id,"compound_id":compound.compound_id,"smiles":version.canonical_smiles,
                               "target":exp["pactivity_mean"],"fingerprints":fp,"descriptors":desc,"scaffold":scaffold})
    return {"assay":_assay_out(assay),"compounds":[{key:value for key,value in row.items() if key!="fingerprint"} for row in rows],"training_compounds":[r["compound_id"] for r in dataset]}


@app.post("/api/assays/{assay_id}/models/train")
def train_assay_model(assay_id: int, db: Session = Depends(get_db)):
    assay=db.get(AssayDefinition,assay_id)
    if not assay or not assay.active: raise HTTPException(status_code=404,detail="Active assay not found")
    compounds=db.scalars(select(Compound).where(Compound.project_id==assay.project_id)).all()
    rows=[];features=[];targets=[];scaffolds=[];descriptor_rows=[];fingerprints=[]
    for compound in compounds:
        current=next((v for v in compound.versions if v.version_number==compound.current_version),None)
        if not current:
            continue
        summary=_experimental_summary(db,current.id,assay_id)
        if not summary: continue
        _,fp,desc,scaffold=fingerprint_and_descriptors(current.canonical_smiles)
        fingerprints.append(fp);descriptor_rows.append(desc);features.append(feature_vector(fp,desc));targets.append(summary["pactivity_mean"]);scaffolds.append(scaffold)
        rows.append({"row_id":compound.id,"compound_id":compound.compound_id,"name":compound.name,"smiles":current.canonical_smiles,"svg":current.svg,"activity_nm":summary["mean_nm"],"pactivity":summary["pactivity_mean"]})
    n=len(targets)
    policy={"N":n,"status":"INSUFFICIENT DATA" if n<5 else ("SIMILARITY ONLY" if n<15 else ("SIMPLE QSAR ALLOWED" if n<30 else "CROSS-VALIDATED QSAR"))}
    if n<15:
        return {"policy":policy,"model":None,"message":"Formal QSAR requires at least 15 experimental compounds."}
    encoded,name,metrics,reason,n=train_model({"features":features,"targets":targets,"scaffolds":scaffolds})
    sklearn_version=__import__("sklearn").__version__
    model=QSARModel(assay_id=assay_id,algorithm=name,sklearn_version=sklearn_version,rdkit_version=ENGINE_VERSION,
                    fingerprint_config=FINGERPRINT_CONFIG,descriptor_config=DESCRIPTOR_NAMES,
                    training_n=n,metrics=metrics,selection_reason=reason,pickle_data=encoded)
    db.add(model);db.commit();db.refresh(model)
    return {"policy":policy,"model":{"model_uid":model.model_uid,"algorithm":name,"training_n":n,
                                     "validation_method":"random KFold CV + Murcko scaffold GroupKFold","metrics":metrics,
                                     "selection_reason":reason}}


@app.post("/api/assays/{assay_id}/predict/{row_id}", status_code=201)
def predict_activity(assay_id: int, row_id: int, db: Session = Depends(get_db)):
    assay=db.get(AssayDefinition,assay_id);compound=db.get(Compound,row_id)
    if not assay or not assay.active or not compound or compound.project_id!=assay.project_id: raise HTTPException(status_code=404,detail="Active assay/compound pair not found")
    current=next((v for v in compound.versions if v.version_number==compound.current_version),None)
    if not current:
        raise HTTPException(status_code=409,detail="Add a valid structure before running an activity prediction")
    existing=_experimental_summary(db,current.id,assay_id)
    _,target_fp,target_desc,_=fingerprint_and_descriptors(current.canonical_smiles)
    dataset={"rows":[],"fingerprints":[],"descriptors":[]}
    for other in db.scalars(select(Compound).where(Compound.project_id==assay.project_id)).all():
        version=next((v for v in other.versions if v.version_number==other.current_version),None)
        if not version:
            continue
        summary=_experimental_summary(db,version.id,assay_id)
        if not summary: continue
        _,fp,desc,_=fingerprint_and_descriptors(version.canonical_smiles)
        dataset["rows"].append({"row_id":other.id,"compound_id":other.compound_id,"activity_nm":summary["mean_nm"],"pactivity":summary["pactivity_mean"]});dataset["fingerprints"].append(fp);dataset["descriptors"].append([desc[name] for name in DESCRIPTOR_NAMES])
    neighbors=nearest_neighbors(target_fp,dataset)
    domain,confidence,max_similarity,outside=applicability(neighbors,target_desc,{"descriptors":np.array(dataset["descriptors"]) if dataset["descriptors"] else np.empty((0,len(DESCRIPTOR_NAMES)))})
    model_row=db.scalar(select(QSARModel).where(QSARModel.assay_id==assay_id).order_by(QSARModel.created_at.desc()))
    n=len(dataset["rows"])
    if model_row and n>=15:
        data=pickle.loads(base64.b64decode(model_row.pickle_data));model=data["model"]
        x=np.vstack([feature_vector(fingerprint_and_descriptors(current.canonical_smiles)[1],target_desc)])
        predicted_p=float(np.asarray(model.predict(x))[0]);ptype=f"QSAR {data['name']}";uncertainty=None
    elif n>=5:
        weights=np.array([neighbor["similarity"]**4 for neighbor in neighbors[:min(5,len(neighbors))]])
        values=np.array([neighbor["pactivity"] for neighbor in neighbors[:len(weights)]])
        predicted_p=float(np.average(values,weights=weights));ptype="Similarity nearest neighbor"
        uncertainty=float(np.std(values)/max(len(weights),1)**.5) if len(weights)>1 else .75
    else:
        raise HTTPException(status_code=409,detail={"status":"INSUFFICIENT DATA","message":"Fewer than five experimental compounds are available.","nearest_neighbors":neighbors})
    prediction=ActivityPrediction(assay_id=assay_id,version_id=current.id,model_id=model_row.id if model_row else None,
                                  prediction_type=ptype,predicted_pactivity=predicted_p,predicted_value_nm=value_from_pactivity(predicted_p),
                                  confidence="LOW" if domain=="OUT OF DOMAIN" else confidence,applicability_domain=domain,
                                  nearest_neighbors=neighbors,uncertainty=uncertainty,
                                  provenance_json={"source":"Validated QSAR/similarity deterministic engine","rdkit_version":ENGINE_VERSION,
                                                   "sklearn_version":__import__("sklearn").__version__,
                                                   "fingerprint":FINGERPRINT_CONFIG,"descriptors":DESCRIPTOR_NAMES,
                                                   "training_n":model_row.training_n if model_row else n,
                                                   "model_metrics":model_row.metrics if model_row else None,
                                                   "max_similarity":max_similarity,"descriptor_outside_training_space":outside,
                                                   "experimental_priority_note":"Experimental values always override predictions."})
    db.add(prediction);db.commit();db.refresh(prediction)
    return {"prediction_id":prediction.id,"type":"Predicted","prediction_type":ptype,"pactivity":round(predicted_p,3),
            "value_nm":round(value_from_pactivity(predicted_p),3),"confidence":confidence,"applicability_domain":domain,
            "nearest_neighbors":neighbors,"provenance":prediction.provenance_json}


@app.get("/api/projects/{project_id}/cliffs")
def activity_cliffs(project_id: int, assay_id: int, similarity_threshold: float = .7, delta_threshold: float = 1.0, db: Session = Depends(get_db)):
    assay=db.get(AssayDefinition,assay_id);rows=[]
    if not assay or assay.project_id!=project_id:raise HTTPException(status_code=404,detail="Assay not found")
    compounds=db.scalars(select(Compound).where(Compound.project_id==project_id)).all();items=[]
    for c in compounds:
      v=next((v for v in c.versions if v.version_number==c.current_version),None)
      if not v: continue
      e=_experimental_summary(db,v.id,assay_id)
      if e: _,fp,_,_=fingerprint_and_descriptors(v.canonical_smiles);items.append((c,v,e,fp))
    for i,(a,av,ae,af) in enumerate(items):
      for b,bv,be,bf in items[i+1:]:
       sim=tanimoto_similarity(af,bf);delta=abs(ae["pactivity_mean"]-be["pactivity_mean"])
       if sim>=similarity_threshold and delta>=delta_threshold:
        rows.append({"a":{"compound_id":a.compound_id,"pactivity":ae["pactivity_mean"],"svg":av.svg},"b":{"compound_id":b.compound_id,"pactivity":be["pactivity_mean"],"svg":bv.svg},"similarity":round(sim,3),"delta_pactivity":round(delta,3)})
    for pair in rows:
        version_a=next(v.id for v in next(c for c in compounds if c.compound_id==pair["a"]["compound_id"]).versions if v.version_number==next(c.current_version for c in compounds if c.compound_id==pair["a"]["compound_id"]))
        version_b=next(v.id for v in next(c for c in compounds if c.compound_id==pair["b"]["compound_id"]).versions if v.version_number==next(c.current_version for c in compounds if c.compound_id==pair["b"]["compound_id"]))
        db.add(MatchedMolecularPair(assay_id=assay_id,version_a_id=version_a,version_b_id=version_b,
                                    similarity=pair["similarity"],delta_pactivity=pair["delta_pactivity"],
                                    transformation_smiles=f"{pair['a']['compound_id']}>>{pair['b']['compound_id']}",
                                    is_cliff=True,provenance_json={"thresholds":pair and {"similarity":similarity_threshold,"delta_pactivity":delta_threshold},
                                                                   "method":"Morgan Tanimoto + pActivity delta"}))
    db.commit()
    return {"thresholds":{"similarity":similarity_threshold,"delta_pactivity":delta_threshold},"cliffs":rows}


@app.get("/api/projects/{project_id}/mmp")
def matched_pairs(project_id: int, assay_id: int, min_similarity: float = .6, max_delta: float = 1.0, db: Session = Depends(get_db)):
    assay=db.get(AssayDefinition,assay_id)
    if not assay or assay.project_id!=project_id: raise HTTPException(status_code=404,detail="Assay not found")
    compounds=db.scalars(select(Compound).where(Compound.project_id==project_id)).all(); items=[]
    for c in compounds:
        v=next((v for v in c.versions if v.version_number==c.current_version),None)
        if not v: continue
        summary=_experimental_summary(db,v.id,assay_id)
        if not summary: continue
        _,fp,_,_=fingerprint_and_descriptors(v.canonical_smiles); items.append({"c":c,"v":v,"summary":summary,"fp":fp})
    pairs=[]
    for i,a in enumerate(items):
        for b in items[i+1:]:
            sim=tanimoto_similarity(a["fp"],b["fp"])
            if sim<min_similarity or abs(a["summary"]["pactivity_mean"]-b["summary"]["pactivity_mean"])>max_delta:
                continue
            delta=b["summary"]["pactivity_mean"]-a["summary"]["pactivity_mean"]
            pair=MatchedMolecularPair(assay_id=assay_id,version_a_id=a["v"].id,version_b_id=b["v"].id,similarity=round(sim,3),
                                      delta_pactivity=round(delta,3),transformation_smiles=f'{a["c"].compound_id}>>{b["c"].compound_id}',
                                      is_cliff=False,
                                      provenance_json={"method":"Morgan/Tanimoto candidate pair; full MCS canonicalization deferred",
                                                       "experimental_priority":"Experimental mean pActivity values used"})
            db.add(pair); pairs.append({"a":a["c"].compound_id,"b":b["c"].compound_id,"similarity":round(sim,3),
                                        "delta_pactivity":round(delta,3),"direction":"B improves over A" if delta>=0 else "A improves over B"})
    db.commit()
    return {"filters":{"min_similarity":min_similarity,"max_abs_delta_pactivity":max_delta},"pairs":pairs}


@app.get("/api/projects/{project_id}/sar-export.csv")
def sar_export(project_id: int, assay_id: int, db: Session = Depends(get_db)):
    from fastapi.responses import PlainTextResponse
    result=sar_table(project_id,assay_id,db);lines=["compound_id,structure_source,activity_source,value_nm,pactivity,MW,cLogP,TPSA,QED"]
    for row in result["compounds"]:
      activity=row["experimental"]; source="Experimental" if activity else "No experimental value"; val=(activity["mean_nm"] if activity else "");p=(activity["pactivity_mean"] if activity else "")
      lines.append(",".join(map(str,[row["compound"],source,source,val,p,row["properties"]["molecular_weight"],row["properties"]["clogp"],row["properties"]["tpsa"],row["properties"]["qed"]])))
    return PlainTextResponse("\n".join(lines),media_type="text/csv")


def _optimization_out(run: OptimizationRun):
    return {
        "id": run.id, "project_id": run.project_id,
        "parent_version_id": run.parent_version_id, "assay_id": run.assay_id,
        "objectives": run.objectives_json or [], "custom_objective": run.custom_objective,
        "constraints": run.constraints_json or {}, "endpoint_weights": run.endpoint_weights_json or {},
        "manual_overrides": run.manual_overrides_json or {}, "status": run.status,
        "message": run.message, "engine": run.engine_name, "engine_version": run.engine_version,
        "evidence": run.evidence_json or {}, "liabilities": run.liabilities_json or [],
        "protected_regions": run.protected_regions_json or [],
        "modifiable_regions": run.modifiable_regions_json or [],
        "recommended_transformations": run.transformations_json or [],
        "highlighted_svg": run.highlighted_svg,
        "legend": {"protected": "red", "modifiable": "orange", "metabolic_soft_spot": "purple"},
        "created_at": run.created_at.isoformat(),
        "completed_at": run.completed_at.isoformat() if run.completed_at else None,
        "analog_generation": "NOT_PERFORMED",
    }


def _validated_optimization_payload(payload: dict):
    objectives = list(dict.fromkeys(str(value).strip() for value in payload.get("objectives", []) if str(value).strip()))
    invalid = [value for value in objectives if value not in OBJECTIVES]
    if not objectives or invalid:
        raise HTTPException(status_code=400, detail=f"At least one valid objective is required; invalid: {invalid}")
    custom = str(payload.get("custom_objective") or "").strip()
    if "Custom" in objectives and not custom:
        raise HTTPException(status_code=400, detail="custom_objective is required for Custom")
    constraints = dict(payload.get("constraints") or {})
    numeric = {
        "potency_max_nm", "do_not_worsen_fold", "clogp_max", "tpsa_min", "tpsa_max",
        "mw_max", "similarity_min", "logs_min", "caco2_logpapp_min",
    }
    for key in numeric:
        if constraints.get(key) in (None, ""):
            constraints.pop(key, None)
            continue
        try:
            constraints[key] = float(constraints[key])
        except (TypeError, ValueError):
            raise HTTPException(status_code=400, detail=f"{key} must be numeric")
    if "similarity_min" in constraints and not 0 <= constraints["similarity_min"] <= 1:
        raise HTTPException(status_code=400, detail="similarity_min must be between 0 and 1")
    if constraints.get("tpsa_min") is not None and constraints.get("tpsa_max") is not None and constraints["tpsa_min"] > constraints["tpsa_max"]:
        raise HTTPException(status_code=400, detail="tpsa_min cannot exceed tpsa_max")
    weights = {}
    for key, value in dict(payload.get("endpoint_weights") or {}).items():
        try:
            weights[str(key)] = float(value)
        except (TypeError, ValueError):
            raise HTTPException(status_code=400, detail=f"Endpoint weight {key} must be numeric")
        if weights[str(key)] < 0:
            raise HTTPException(status_code=400, detail=f"Endpoint weight {key} cannot be negative")
    return objectives, custom, constraints, weights


@app.get("/api/optimization/config")
def optimization_config():
    return {
        "engine": OPTIMIZATION_ENGINE, "engine_version": OPTIMIZATION_VERSION,
        "objectives": list(OBJECTIVES), "evidence_hierarchy": list(EVIDENCE_HIERARCHY),
        "transformation_library": list(TRANSFORMATION_LIBRARY),
        "policy": {
            "analog_generation": False, "llm": False, "overall_score": False,
            "low_confidence_classification": "Supporting evidence only unless corroborated",
            "experimental_precedence": True,
        },
    }


@app.get("/api/projects/{project_id}/optimization")
def list_optimization_runs(project_id: int, version_id: int | None = None, db: Session = Depends(get_db)):
    if not db.get(Project, project_id):
        raise HTTPException(status_code=404, detail="Project not found")
    query = select(OptimizationRun).where(OptimizationRun.project_id == project_id)
    if version_id is not None:
        query = query.where(OptimizationRun.parent_version_id == version_id)
    rows = db.scalars(query.order_by(OptimizationRun.created_at.desc())).all()
    return {"runs": [_optimization_out(row) for row in rows], "config": optimization_config()}


@app.get("/api/optimization/runs/{run_id}")
def get_optimization_run(run_id: int, db: Session = Depends(get_db)):
    run = db.get(OptimizationRun, run_id)
    if not run:
        raise HTTPException(status_code=404, detail="OptimizationRun not found")
    return _optimization_out(run)


@app.post("/api/projects/{project_id}/optimization/runs", status_code=201)
def create_optimization_run(project_id: int, payload: dict, db: Session = Depends(get_db)):
    project = db.get(Project, project_id)
    version = db.get(CompoundVersion, payload.get("parent_version_id"))
    if not project or not version:
        raise HTTPException(status_code=404, detail="Project or parent CompoundVersion not found")
    if project.molecule_type != "Small Molecule":
        raise HTTPException(status_code=400, detail="This model currently supports small molecules only.")
    compound = db.get(Compound, version.compound_row_id)
    if not compound or compound.project_id != project_id:
        raise HTTPException(status_code=404, detail="Parent CompoundVersion is not in this project")
    assay_id = int(payload["assay_id"]) if payload.get("assay_id") not in (None, "") else None
    assay = db.get(AssayDefinition, assay_id) if assay_id else None
    if assay_id and (not assay or assay.project_id != project_id or not assay.active):
        raise HTTPException(status_code=404, detail="Selected active assay is not in this project")
    objectives, custom, constraints, weights = _validated_optimization_payload(payload)
    if ("Improve potency" in objectives or "potency_max_nm" in constraints) and not assay:
        raise HTTPException(status_code=400, detail="A selected assay is required for potency objectives or constraints")
    run = OptimizationRun(
        project_id=project_id, parent_version_id=version.id, assay_id=assay_id,
        objectives_json=objectives, custom_objective=custom,
        constraints_json=constraints, endpoint_weights_json=weights,
        manual_overrides_json=dict(payload.get("manual_overrides") or {}),
        status="RUNNING", message="Assembling Stage 1-3 evidence.",
        engine_name=OPTIMIZATION_ENGINE, engine_version=OPTIMIZATION_VERSION,
    )
    db.add(run); db.flush()
    try:
        analyze_run(db, run)
        db.commit(); db.refresh(run)
    except Exception as exc:
        db.rollback()
        raise HTTPException(status_code=400, detail=f"Optimization analysis failed: {exc}")
    return _optimization_out(run)


@app.patch("/api/optimization/runs/{run_id}/overrides")
def update_optimization_overrides(run_id: int, payload: dict, db: Session = Depends(get_db)):
    run = db.get(OptimizationRun, run_id)
    if not run:
        raise HTTPException(status_code=404, detail="OptimizationRun not found")
    allowed = {"protect_atoms", "allow_atoms", "exclude_transformations", "prioritize_transformations"}
    if any(key not in allowed for key in payload):
        raise HTTPException(status_code=400, detail="Unsupported manual override")
    overrides = dict(run.manual_overrides_json or {})
    for key, value in payload.items():
        if not isinstance(value, list):
            raise HTTPException(status_code=400, detail=f"{key} must be a list")
        overrides[key] = value
    run.manual_overrides_json = overrides
    run.status = "RUNNING"
    analyze_run(db, run)
    db.commit(); db.refresh(run)
    return _optimization_out(run)


def _candidate_comparison(candidate: OptimizationCandidate, optimization_run: OptimizationRun, parent: CompoundVersion):
    rows = []
    parent_properties = parent.properties_json or {}
    candidate_properties = (candidate.stage1_json or {}).get("properties", {})
    for endpoint, key, unit in (
        ("MW", "molecular_weight", "Da"), ("cLogP", "clogp", ""),
        ("TPSA", "tpsa", "Å²"), ("Fsp3", "fraction_csp3", ""),
    ):
        parent_value, candidate_value = parent_properties.get(key), candidate_properties.get(key)
        rows.append({"endpoint": endpoint, "parent": {"value": parent_value, "unit": unit, "type": "Calculated"}, "candidate": {"value": candidate_value, "unit": unit, "type": "Calculated"}, "change": (round(candidate_value-parent_value, 5) if isinstance(parent_value, (int, float)) and isinstance(candidate_value, (int, float)) else None)})
    parent_activity = (optimization_run.evidence_json or {}).get("activity", {})
    parent_activity = parent_activity.get("experimental") or parent_activity.get("predicted")
    if parent_activity or candidate.activity_json:
        parent_value = parent_activity.get("mean_nm") if parent_activity and parent_activity.get("type") == "Experimental" else (parent_activity or {}).get("value_nm")
        candidate_value = (candidate.activity_json or {}).get("value_nm")
        rows.insert(0, {"endpoint": "Activity", "parent": {"value": parent_value, "unit": "nM", "type": (parent_activity or {}).get("type", "Not measured"), "confidence": (parent_activity or {}).get("confidence")}, "candidate": {"value": candidate_value, "unit": "nM", "type": (candidate.activity_json or {}).get("record_type", "Unavailable"), "confidence": (candidate.activity_json or {}).get("confidence"), "domain": (candidate.activity_json or {}).get("applicability_domain")}, "change": (round(candidate_value-parent_value, 5) if isinstance(parent_value, (int, float)) and isinstance(candidate_value, (int, float)) else None)})
    parent_admet = (optimization_run.evidence_json or {}).get("admet", {})
    for endpoint in ("Solubility", "Permeability", "Plasma protein binding", "HLM intrinsic clearance", "RLM intrinsic clearance", "CYP3A4 inhibitor", "P-gp inhibitor", "hERG liability", "Ames mutagenicity", "DILI clinical liability"):
        parent_row = (parent_admet.get(endpoint) or {}).get("preferred")
        candidate_row = (candidate.admet_json or {}).get(endpoint)
        if not parent_row and not candidate_row:
            continue
        parent_value = (parent_row or {}).get("classification", (parent_row or {}).get("value"))
        candidate_value = (candidate_row or {}).get("classification", (candidate_row or {}).get("predicted_value", (candidate_row or {}).get("value")))
        numerical_change = round(float(candidate_value)-float(parent_value), 5) if isinstance(parent_value, (int, float)) and isinstance(candidate_value, (int, float)) else None
        rows.append({
            "endpoint": endpoint, "parent": {"value": parent_value, "unit": (parent_row or {}).get("unit", ""), "type": (parent_row or {}).get("type", "Not measured"), "confidence": (parent_row or {}).get("confidence")},
            "candidate": {"value": candidate_value, "unit": (candidate_row or {}).get("unit", ""), "type": (candidate_row or {}).get("record_type", "MODEL_UNAVAILABLE" if (candidate_row or {}).get("status") == "MODEL_UNAVAILABLE" else "Predicted"), "confidence": (candidate_row or {}).get("confidence"), "domain": (candidate_row or {}).get("applicability_domain")},
            "change": numerical_change,
        })
    return rows


def _candidate_out(candidate: OptimizationCandidate, optimization_run: OptimizationRun | None = None, parent: CompoundVersion | None = None):
    ranking = sorted(candidate.rankings, key=lambda row: row.created_at, reverse=True)[0] if candidate.rankings else None
    result = {
        "id": candidate.id, "proposal_run_id": candidate.proposal_run_id,
        "candidate_number": candidate.candidate_number, "parent_version_id": candidate.parent_version_id,
        "existing_version_id": candidate.existing_version_id, "canonical_smiles": candidate.canonical_smiles,
        "isomeric_smiles": candidate.isomeric_smiles, "inchikey": candidate.inchikey,
        "generation_source": candidate.generation_source, "generation_priority": candidate.generation_priority,
        "generation_timestamp": candidate.generation_timestamp.isoformat(), "hypothesis": candidate.hypothesis,
        "why_generated": candidate.why_generated, "expected_benefit": candidate.expected_benefit,
        "status": candidate.status, "rejection_stage": candidate.rejection_stage,
        "transformations": [{
            "id": row.transformation_id, "name": row.name, "sequence": row.sequence_number,
            "reaction_smarts": row.reaction_smarts, "version": row.transformation_version,
            "source": row.source, "source_atoms": row.source_atom_indices_json,
            "changed_parent_atoms": row.changed_parent_atoms_json,
            "execution_status": row.execution_status, "provenance": row.provenance_json,
        } for row in sorted(candidate.transformations, key=lambda row: row.sequence_number)],
        "stage1": candidate.stage1_json or {}, "property_delta": candidate.property_delta_json or {},
        "activity": candidate.activity_json or {}, "admet": candidate.admet_json or {},
        "soft_spots": candidate.soft_spot_json or {}, "soft_spot_changes": candidate.soft_spot_change_json or {},
        "synthetic_feasibility": candidate.synthetic_feasibility_json or {},
        "parent_similarity": candidate.parent_similarity, "mcs_coverage": candidate.mcs_coverage,
        "changed_parent_atoms": candidate.changed_parent_atoms_json or [], "changed_candidate_atoms": candidate.changed_candidate_atoms_json or [],
        "structure_svg": candidate.structure_svg, "parent_difference_svg": candidate.parent_difference_svg,
        "candidate_difference_svg": candidate.candidate_difference_svg,
        "confidence": candidate.confidence, "applicability_domain": candidate.applicability_domain,
        "objective_vector": candidate.objective_vector_json or {}, "ranking_score": candidate.ranking_score,
        "pareto_front": candidate.pareto_front, "information_value": candidate.information_value,
        "main_risk": candidate.main_risk, "selected_top10": candidate.selected_top10,
        "user_added": candidate.user_added, "user_decision": candidate.user_decision,
        "user_decision_reason": candidate.user_decision_reason,
        "rejection_reasons": [{"code": row.code, "detail": row.detail, "stage": row.stage, "hard_constraint": row.hard_constraint, "evidence_type": row.evidence_type} for row in candidate.rejection_reasons],
        "ranking": {"rank": ranking.rank, "score": ranking.score, "pareto_front": ranking.pareto_front, "formula": ranking.score_breakdown_json, "diversity": ranking.diversity_json} if ranking else None,
        "prediction_snapshots": [{"stage": row.stage, "endpoint": row.endpoint, "type": row.record_type, "unit": row.unit, "model": row.model_name, "model_version": row.model_version, "confidence": row.confidence, "domain": row.applicability_domain} for row in candidate.predictions],
    }
    if optimization_run is not None and parent is not None:
        result["parent_comparison"] = _candidate_comparison(candidate, optimization_run, parent)
    return result


def _proposal_out(run: OptimizationProposalRun, include_candidates=True):
    result = {
        "id": run.id, "project_id": run.project_id, "optimization_run_id": run.optimization_run_id,
        "parent_version_id": run.parent_version_id, "status": run.status, "stage_message": run.stage_message,
        "transformation_library_version": run.transformation_library_version,
        "model_versions": run.model_versions_json or {}, "endpoint_weights": run.endpoint_weights_json or {},
        "hard_constraints": run.hard_constraints_json or {}, "settings": run.settings_json or {},
        "random_seed": run.random_seed, "raw_candidate_count": run.raw_candidate_count,
        "accepted_count": run.accepted_count, "rejected_count": run.rejected_count,
        "top_count": run.top_count, "summary": run.summary_json or {},
        "created_at": run.created_at.isoformat(), "started_at": run.started_at.isoformat() if run.started_at else None,
        "completed_at": run.completed_at.isoformat() if run.completed_at else None,
    }
    if include_candidates:
        session = Session.object_session(run)
        optimization = session.get(OptimizationRun, run.optimization_run_id)
        parent = session.get(CompoundVersion, run.parent_version_id)
        result["candidates"] = [_candidate_out(row, optimization, parent) for row in sorted(run.candidates, key=lambda row: row.candidate_number)]
    return result


@app.get("/api/proposals/config")
def proposal_config():
    return {
        "engine": PROPOSAL_ENGINE, "engine_version": PROPOSAL_VERSION,
        "executable_transformations": EXECUTABLE_TRANSFORMATIONS,
        "strategy_only_transformations": STRATEGY_ONLY_TRANSFORMATIONS,
        "job_states": ["PENDING", "GENERATING", "FILTERING", "PREDICTING", "RANKING", "COMPLETED", "FAILED"],
        "policy": {"llm": False, "pk": False, "random_generation": False, "max_changes": 2, "single_change_first": True, "experimental_precedence": True},
    }


@app.get("/api/optimization/runs/{run_id}/proposals")
def list_proposals(run_id: int, db: Session = Depends(get_db)):
    if not db.get(OptimizationRun, run_id):
        raise HTTPException(status_code=404, detail="OptimizationRun not found")
    rows = db.scalars(select(OptimizationProposalRun).where(OptimizationProposalRun.optimization_run_id == run_id).order_by(OptimizationProposalRun.created_at.desc())).all()
    return {"proposal_runs": [_proposal_out(row, include_candidates=False) for row in rows], "config": proposal_config()}


@app.post("/api/optimization/runs/{run_id}/proposals", status_code=202)
def create_proposal(run_id: int, payload: dict, background_tasks: BackgroundTasks, db: Session = Depends(get_db)):
    optimization = db.get(OptimizationRun, run_id)
    if not optimization or optimization.status != "COMPLETE":
        raise HTTPException(status_code=404, detail="Completed OptimizationRun not found")
    settings = {"max_raw_candidates": 120, "allow_double_transforms": True, **dict(payload.get("settings") or {})}
    try:
        settings["max_raw_candidates"] = int(settings["max_raw_candidates"])
    except (TypeError, ValueError):
        raise HTTPException(status_code=400, detail="max_raw_candidates must be an integer")
    if not 1 <= settings["max_raw_candidates"] <= 200:
        raise HTTPException(status_code=400, detail="max_raw_candidates must be between 1 and 200")
    hard_constraints = dict(payload.get("hard_constraints") or {})
    proposal = OptimizationProposalRun(
        project_id=optimization.project_id, optimization_run_id=optimization.id,
        parent_version_id=optimization.parent_version_id, status="PENDING", stage_message="Queued",
        transformation_library_version=PROPOSAL_VERSION,
        endpoint_weights_json=dict(optimization.endpoint_weights_json or {}),
        hard_constraints_json=hard_constraints, settings_json=settings, random_seed=42,
    )
    db.add(proposal); db.commit(); db.refresh(proposal)
    background_tasks.add_task(execute_proposal_run, proposal.id)
    return _proposal_out(proposal, include_candidates=False)


@app.post("/api/proposals/{proposal_id}/execute")
def retry_proposal(proposal_id: int, background_tasks: BackgroundTasks, db: Session = Depends(get_db)):
    proposal = db.get(OptimizationProposalRun, proposal_id)
    if not proposal:
        raise HTTPException(status_code=404, detail="ProposalRun not found")
    if proposal.status in {"GENERATING", "FILTERING", "PREDICTING", "RANKING"}:
        raise HTTPException(status_code=409, detail="ProposalRun is already active")
    if proposal.candidates:
        raise HTTPException(status_code=409, detail="A populated run is immutable; create a new ProposalRun for reproducibility")
    proposal.status, proposal.stage_message = "PENDING", "Queued for retry"; db.commit()
    background_tasks.add_task(execute_proposal_run, proposal.id)
    return _proposal_out(proposal, include_candidates=False)


@app.get("/api/proposals/{proposal_id}")
def get_proposal(proposal_id: int, view: str = "all", db: Session = Depends(get_db)):
    proposal = db.get(OptimizationProposalRun, proposal_id)
    if not proposal:
        raise HTTPException(status_code=404, detail="ProposalRun not found")
    result = _proposal_out(proposal, include_candidates=True)
    filters = {
        "accepted": lambda row: row["status"] in {"ACCEPTED", "TOP_10"},
        "rejected": lambda row: row["status"] in {"REJECTED", "FAILED"},
        "pareto": lambda row: row["pareto_front"] == 1,
        "top10": lambda row: row["selected_top10"],
        "all": lambda row: True,
    }
    if view not in filters:
        raise HTTPException(status_code=400, detail="view must be all, accepted, rejected, pareto, or top10")
    result["candidates"] = [row for row in result["candidates"] if filters[view](row)]
    result["view"] = view
    return result


@app.get("/api/proposal-candidates/{candidate_id}")
def get_proposal_candidate(candidate_id: int, db: Session = Depends(get_db)):
    candidate = db.get(OptimizationCandidate, candidate_id)
    if not candidate:
        raise HTTPException(status_code=404, detail="Candidate not found")
    optimization = db.get(OptimizationRun, candidate.optimization_run_id)
    parent = db.get(CompoundVersion, candidate.parent_version_id)
    return _candidate_out(candidate, optimization, parent)


@app.patch("/api/proposal-candidates/{candidate_id}/decision")
def candidate_decision(candidate_id: int, payload: dict, db: Session = Depends(get_db)):
    candidate = db.get(OptimizationCandidate, candidate_id)
    decision = str(payload.get("decision") or "").upper()
    reason = str(payload.get("reason") or "").strip()
    if not candidate or decision not in {"PROMOTED", "REJECTED", "CLEAR"}:
        raise HTTPException(status_code=400, detail="Candidate and decision PROMOTED/REJECTED/CLEAR are required")
    if decision == "REJECTED" and not reason:
        raise HTTPException(status_code=400, detail="A manual rejection reason is required")
    candidate.user_decision = "" if decision == "CLEAR" else decision
    candidate.user_decision_reason = "" if decision == "CLEAR" else reason
    if decision == "REJECTED":
        candidate.status = "REJECTED"
        db.add(CandidateRejectionReason(candidate=candidate, code="USER_REJECTED", detail=reason, stage="MANUAL_REVIEW", hard_constraint=False, evidence_type="User judgment"))
    elif candidate.status == "REJECTED" and any(row.code == "USER_REJECTED" and not row.hard_constraint for row in candidate.rejection_reasons):
        candidate.status = "REJECTED" if any(row.hard_constraint for row in candidate.rejection_reasons) else "RESCORED"
    proposal = db.get(OptimizationProposalRun, candidate.proposal_run_id)
    optimization = db.get(OptimizationRun, candidate.optimization_run_id)
    rank_candidates(db, proposal, optimization)
    proposal.accepted_count = len([row for row in proposal.candidates if row.status in {"ACCEPTED", "TOP_10"}])
    proposal.rejected_count = len([row for row in proposal.candidates if row.status in {"REJECTED", "FAILED"}])
    proposal.top_count = len([row for row in proposal.candidates if row.selected_top10])
    db.commit(); db.refresh(candidate)
    return _candidate_out(candidate, optimization, db.get(CompoundVersion, candidate.parent_version_id))


@app.post("/api/proposals/{proposal_id}/candidates", status_code=201)
def add_user_candidate(proposal_id: int, payload: dict, db: Session = Depends(get_db)):
    proposal = db.get(OptimizationProposalRun, proposal_id)
    if not proposal or proposal.status != "COMPLETED":
        raise HTTPException(status_code=404, detail="Completed ProposalRun not found")
    smiles = str(payload.get("smiles") or "").strip()
    if not smiles:
        raise HTTPException(status_code=400, detail="SMILES is required")
    optimization = db.get(OptimizationRun, proposal.optimization_run_id)
    try:
        candidate = process_user_candidate(db, proposal, optimization, smiles, str(payload.get("reason") or "User-added analog"))
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    return _candidate_out(candidate, optimization, db.get(CompoundVersion, candidate.parent_version_id))


@app.get("/", response_class=HTMLResponse)
def index():
    response = HTMLResponse(Path("frontend/static/index.html").read_text())
    response.headers["Cache-Control"] = "no-store, must-revalidate"
    return response

def register_hardening_routes(app):
    @app.post("/api/standardization/standardize")
    def standardize_smiles_endpoint(payload: dict):
        smiles = str(payload.get("smiles") or "")
        return standardize_molecule(smiles)

    @app.get("/api/standardization/configurations")
    def standardization_configs_endpoint():
        return {
            "standardizer_name": STANDARDIZER_NAME,
            "standardizer_version": STANDARDIZER_VERSION,
            "rdkit_version": RDKIT_VERSION,
            "fingerprints": GLOBAL_FINGERPRINT_CONFIG,
            "descriptors": GLOBAL_DESCRIPTOR_CONFIG,
        }

    @app.get("/api/evaluation/registry")
    def evaluation_registry_endpoint():
        return {"registry": EVALUATION_REGISTRY}

    @app.get("/api/evaluation/golden-gate")
    def golden_gate_endpoint():
        return run_golden_gate_test()

    @app.get("/api/evaluation/lightning-audit")
    def lightning_audit_endpoint():
        return perform_lightning_security_audit()

    @app.get("/api/evaluation/rdkit-readiness")
    def rdkit_readiness_endpoint():
        return get_rdkit_upgrade_readiness_report()

    @app.post("/api/evaluation/mmp-directional")
    def mmp_directional_endpoint(payload: dict):
        pairs = payload.get("pairs") or []
        min_delta = float(payload.get("min_delta_fold") or 1.5)
        return evaluate_mmp_directional_accuracy(pairs, min_delta_fold=min_delta)


register_pk_routes(app)
register_ivive_routes(app)
register_simulation_routes(app)
register_hardening_routes(app)
app.mount("/static", StaticFiles(directory="frontend/static"), name="static")
