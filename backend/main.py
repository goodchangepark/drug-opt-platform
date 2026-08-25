import hashlib
import json

import numpy as np

from fastapi import Depends, FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from rdkit import Chem
from sqlalchemy import func, select
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session

from .chemistry import ChemistryError, ENGINE, ENGINE_VERSION, analyze_smiles
from .activity_models import ActivityMeasurement, ActivityPrediction, AssayDefinition, MatchedMolecularPair, QSARModel
from .admet import (ADMETEndpoint, ADMETMeasurement, ADMETModelRegistry, ADMETPrediction, ADMETPredictionRun,
                    csv_export, ensure_admet_schema, inputs_hash,
                    measurement_out, parse_csv, validate_measurement)
from .admet_predictor import (MODEL_SPECS, MODEL_VERSION, comparison_for_prediction, cyp_experimental_evidence,
                              metabolic_stability_assessment, model_files_available, predict_endpoint)
from .database import Base, SessionLocal, engine, get_db
from .models import Compound, CompoundVersion, PredictionRun, Project, PropertyCalculation, StructuralAlert, utcnow
from .qsar import (DESCRIPTOR_NAMES, FINGERPRINT_CONFIG, applicability, feature_vector,
                   fingerprint_and_descriptors, nearest_neighbors, normalize_concentration, tanimoto_similarity,
                   pactivity, train_model, value_from_pactivity)
from .schemas import CompoundCreate, CompoundUpdate, ProjectCreate, ProjectOut, ProjectUpdate

app = FastAPI(title="AI Drug Optimization Platform", version="0.3.0-stage3c")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])


@app.on_event("startup")
def startup():
    if not app.dependency_overrides:
        Base.metadata.create_all(bind=engine)
        ensure_admet_schema(engine)
    import backend.activity_models


def _project_out(db: Session, project: Project):
    count = db.scalar(select(func.count(Compound.id)).where(Compound.project_id == project.id)) or 0
    return ProjectOut.model_validate(project).model_copy(update={"compound_count": count})


@app.get("/api/health")
def health():
    return {"status": "ok", "stage": 3, "engine": ENGINE, "engine_version": ENGINE_VERSION}


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


@app.post("/api/projects", status_code=201)
def create_project(payload: ProjectCreate, db: Session = Depends(get_db)):
    project = Project(**payload.model_dump())
    db.add(project)
    try:
        db.commit()
    except IntegrityError:
        db.rollback()
        raise HTTPException(status_code=409, detail="Project name already exists")
    db.refresh(project)
    return _project_out(db, project)


@app.get("/api/projects/{project_id}")
def get_project(project_id: int, db: Session = Depends(get_db)):
    project = db.get(Project, project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    compounds = db.scalars(select(Compound).where(Compound.project_id == project_id).order_by(Compound.compound_id)).all()
    data = _project_out(db, project).model_dump()
    data["compounds"] = [compound_out(compound) for compound in compounds]
    return data


@app.patch("/api/projects/{project_id}")
def update_project(project_id: int, payload: ProjectUpdate, db: Session = Depends(get_db)):
    project = db.get(Project, project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    for key, value in payload.model_dump(exclude_unset=True).items():
        setattr(project, key, value)
    project.updated_at = utcnow()
    try:
        db.commit()
    except IntegrityError:
        db.rollback(); raise HTTPException(status_code=409, detail="Project name already exists")
    db.refresh(project); return _project_out(db, project)


@app.delete("/api/projects/{project_id}", status_code=204)
def delete_project(project_id: int, db: Session = Depends(get_db)):
    project = db.get(Project, project_id)
    if not project: raise HTTPException(status_code=404, detail="Project not found")
    db.delete(project); db.commit()


def compound_out(compound: Compound):
    current = next((v for v in compound.versions if v.version_number == compound.current_version), compound.versions[-1])
    return {
        "row_id": compound.id, "project_id": compound.project_id, "compound_id": compound.compound_id,
        "name": compound.name, "notes": compound.notes, "current_version": compound.current_version,
        "created_at": compound.created_at.isoformat(), "updated_at": compound.updated_at.isoformat(),
        "version": serialize_version(current), "versions": [{"version_number": v.version_number, "canonical_smiles": v.canonical_smiles, "change_note": v.change_note} for v in compound.versions],
    }


def serialize_version(version: CompoundVersion):
    return {
        "id": version.id, "version_number": version.version_number, "original_smiles": version.original_smiles,
        "canonical_smiles": version.canonical_smiles, "isomeric_smiles": version.isomeric_smiles,
        "inchi": version.inchi, "inchikey": version.inchikey, "change_note": version.change_note,
        "properties": version.properties_json or {}, "rules": (version.calculation_json or {}).get("rules", {}),
        "assessment": version.assessment_json or {}, "svg": version.svg,
        "highlighted_svg": version.highlighted_svg, "provenance": (version.calculation_json or {}).get("provenance", {}),
        "alerts": version.alerts_json or [],
    }


def persist_analysis(db: Session, compound: Compound, smiles: str, change_note: str) -> CompoundVersion:
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
        properties_json=analysis["properties"], alerts_json=analysis["alerts"], assessment_json=analysis["assessment"],
        calculation_json={"provenance": analysis["provenance"], "rules": analysis["rules"]},
        svg=analysis["svg"], highlighted_svg=analysis["highlighted_svg"],
    )
    db.add(version); db.flush()
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
    compound.current_version = number; db.commit(); db.refresh(version)
    return version


@app.post("/api/projects/{project_id}/compounds", status_code=201)
def create_compound(project_id: int, payload: CompoundCreate, db: Session = Depends(get_db)):
    project = db.get(Project, project_id)
    if not project: raise HTTPException(status_code=404, detail="Project not found")
    existing_label = db.scalar(select(Compound).where(Compound.project_id == project_id, Compound.compound_id == payload.compound_id))
    if existing_label: raise HTTPException(status_code=409, detail="Compound ID already exists in project")
    compound = Compound(project_id=project_id, compound_id=payload.compound_id, name=payload.name, notes=payload.notes)
    db.add(compound); db.flush()
    try:
        version = persist_analysis(db, compound, payload.smiles, "Initial structure")
    except HTTPException:
        db.rollback(); raise
    except ChemistryError as exc:
        db.rollback(); raise HTTPException(status_code=400, detail=str(exc))
    db.refresh(compound); return compound_out(compound)


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


@app.patch("/api/compounds/{row_id}")
def update_compound(row_id: int, payload: CompoundUpdate, db: Session = Depends(get_db)):
    compound = db.get(Compound, row_id)
    if not compound: raise HTTPException(status_code=404, detail="Compound not found")
    if payload.name is not None: compound.name = payload.name
    if payload.notes is not None: compound.notes = payload.notes
    compound.updated_at = utcnow()
    if payload.smiles:
        try:
            persist_analysis(db, compound, payload.smiles, payload.change_note)
        except HTTPException: db.rollback(); raise
        except ChemistryError as exc: db.rollback(); raise HTTPException(status_code=400, detail=str(exc))
    else:
        db.commit()
    db.refresh(compound); return compound_out(compound)


@app.delete("/api/compounds/{row_id}", status_code=204)
def delete_compound(row_id: int, db: Session = Depends(get_db)):
    compound = db.get(Compound, row_id)
    if not compound: raise HTTPException(status_code=404, detail="Compound not found")
    db.delete(compound); db.commit()


@app.get("/api/projects/{project_id}/compare")
def compare(project_id: int, ids: str = Query(...), db: Session = Depends(get_db)):
    try: wanted = {int(value) for value in ids.split(",") if value.strip()}
    except ValueError: raise HTTPException(status_code=400, detail="ids must be comma-separated integers")
    if len(wanted) < 2: raise HTTPException(status_code=400, detail="Select at least two compounds")
    rows = []
    for row_id in wanted:
        compound = db.get(Compound, row_id)
        if not compound or compound.project_id != project_id: continue
        version = next((v for v in compound.versions if v.version_number == compound.current_version), None)
        p = version.properties_json or {}
        comparison_row = {
            "compound": compound.compound_id, "name": compound.name, "row_id": compound.id,
            "MW": p.get("molecular_weight"), "cLogP": p.get("clogp"), "TPSA": p.get("tpsa"), "HBD": p.get("hbd"),
            "HBA": p.get("hba"), "RotB": p.get("rotatable_bonds"), "Fsp3": p.get("fraction_csp3"), "QED": p.get("qed"),
            "svg": version.svg if version else "", "inchikey": version.inchikey if version else "",
        }
        activity = db.scalar(select(ActivityMeasurement).where(ActivityMeasurement.version_id == version.id).order_by(ActivityMeasurement.created_at.desc()))
        comparison_row["Activity"] = activity.normalized_value_nm if activity else None
        endpoint_map = {
            "HLM intrinsic clearance": "HLM", "RLM intrinsic clearance": "RLM",
            "Plasma protein binding": "PPB", "Solubility": "Solubility", "Permeability": "Caco-2",
        }
        experimental = db.scalars(select(ADMETMeasurement).where(ADMETMeasurement.version_id == version.id)).all()
        endpoint_names = {item.id: item.name for item in db.scalars(select(ADMETEndpoint).where(ADMETEndpoint.project_id == project_id))}
        for endpoint_name, label in endpoint_map.items():
            prediction = db.scalar(
                select(ADMETPrediction).join(ADMETModelRegistry)
                .where(ADMETPrediction.version_id == version.id, ADMETModelRegistry.endpoint_name == endpoint_name)
                .order_by(ADMETPrediction.created_at.desc())
            )
            if not prediction:
                comparison_row[label] = None
                continue
            matches = comparison_for_prediction(endpoint_name, prediction.predicted_value, experimental, endpoint_names)
            comparison_row[label] = matches[0]["experimental_normalized"] if matches else prediction.predicted_value
        rows.append(comparison_row)
    if len(rows) < 2: raise HTTPException(status_code=400, detail="At least two selected compounds must belong to the project")
    property_metrics = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "Fsp3", "QED"]
    metrics = property_metrics + ["Activity", "HLM", "RLM", "PPB", "Solubility", "Caco-2"]
    ranges = {metric: {"min": min(r[metric] for r in rows if r[metric] is not None),
                       "max": max(r[metric] for r in rows if r[metric] is not None)} for metric in property_metrics}
    return {"metrics": metrics, "ranges": ranges, "compounds": rows, "metric_units": {
        "Activity": "nM (latest experimental)", "HLM": "log10(mL/min/kg)", "RLM": "log10(mL/min/kg)",
        "PPB": "% bound", "Solubility": "log10(mol/L)", "Caco-2": "log10(cm/s)",
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
    available, unavailable_reason = model_files_available(model.endpoint_name) if model.endpoint_name in MODEL_SPECS else (
        False, (model.provenance_json or {}).get("reason", "No endpoint-specific model installed in the current stage"),
    )
    return {
        "id": model.id, "endpoint": model.endpoint_name, "model_name": model.model_name,
        "model_version": model.model_version,
        "status": model.implementation_status if available else "MODEL_UNAVAILABLE",
        "active": bool(model.is_active and available), "output_unit": model.output_unit,
        "details": model.provenance_json or {}, "unavailable_reason": unavailable_reason,
    }


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
    return {
        "id": prediction.id, "run_id": prediction.run_id, "version_id": prediction.version_id,
        "endpoint_id": prediction.endpoint_id, "endpoint": prediction.endpoint.name,
        "predicted_value": prediction.predicted_value, "unit": prediction.unit,
        "confidence": prediction.confidence, "applicability_domain": prediction.applicability_domain,
        "uncertainty": prediction.uncertainty, "model": _admet_model_out(prediction.model),
        "outputs": outputs, "experimental_comparisons": comparisons, "preferred_result": preferred,
        "created_at": prediction.created_at.isoformat(), "type": "Predicted",
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
        value=value, unit=str(payload.get("unit", "")).strip(), qualifier=payload.get("qualifier") or "=",
        replicate=str(payload.get("replicate") or "R1"), mean_value=mean_value,
        standard_deviation=sd, sample_size=int(payload["n"]) if payload.get("n") else None,
        method=str(payload.get("method") or ""), source=str(payload.get("source") or "User experimental"),
        experiment_date=str(payload.get("date") or ""), notes=str(payload.get("notes") or ""),
        provenance_json={"data_type": "experimental", **(payload.get("provenance") or {})},
    )
    if not row.unit:
        db.rollback(); raise HTTPException(status_code=400, detail="unit is required")
    db.add(row); db.commit(); db.refresh(row)
    return measurement_out(row)


@app.get("/api/projects/{project_id}/admet")
def list_admet(project_id: int, db: Session = Depends(get_db)):
    project = db.get(Project, project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    versions = {version.id: (compound.compound_id, version.version_number)
                for compound in project.compounds for version in compound.versions}
    rows = db.scalars(select(ADMETMeasurement).join(
        ADMETEndpoint, ADMETEndpoint.id == ADMETMeasurement.endpoint_id
    ).where(ADMETEndpoint.project_id == project_id).order_by(ADMETMeasurement.created_at.desc())).all()
    models = db.scalars(select(ADMETModelRegistry).order_by(ADMETModelRegistry.endpoint_name)).all()
    version_ids = list(versions)
    endpoint_names = {endpoint.id: endpoint.name for endpoint in db.scalars(
        select(ADMETEndpoint).where(ADMETEndpoint.project_id == project_id)
    )}
    predictions = db.scalars(
        select(ADMETPrediction)
        .where(ADMETPrediction.version_id.in_(version_ids))
        .order_by(ADMETPrediction.created_at.desc())
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
    return {
        "endpoints": [_admet_endpoint_out(e) for e in db.scalars(select(ADMETEndpoint).where(ADMETEndpoint.project_id == project_id))],
        "measurements": [measurement_out(row) for row in rows],
        "models": [_admet_model_out(model) for model in models],
        "predictions": [_admet_prediction_out(
            prediction, measurements_by_version.get(prediction.version_id, []), endpoint_names,
        ) for prediction in predictions],
        "prediction_runs": [{"id": r.id, "version_id": r.version_id, "status": r.status,
                             "message": r.message, "started_at": r.started_at.isoformat()} for r in runs],
        "csv_columns": ["compound_id", "version_number"] + [
            column for column in ("endpoint", "species", "matrix", "value", "unit", "qualifier", "replicate",
                                  "mean", "sd", "n", "method", "source", "date", "notes")
        ],
        "labels_by_version": {str(key): value for key, value in versions.items()},
    }


@app.post("/api/projects/{project_id}/admet/measurements", status_code=201)
def create_admet_measurement(project_id: int, payload: dict, db: Session = Depends(get_db)):
    project = db.get(Project, project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    return add_admet_measurement(db, project_id, payload)


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
    for item in preview["rows"]:
        compound = labels[str(item["compound_id"]).strip()]
        version_number = int(item.get("version_number") or compound.current_version)
        version = next(version for version in compound.versions if version.version_number == version_number)
        created.append(add_admet_measurement(db, project_id, {**item, "version_id": version.id}))
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


@app.post("/api/admet/predict/{row_id}", status_code=202)
def run_admet_predictions(row_id: int, db: Session = Depends(get_db)):
    version = db.get(CompoundVersion, row_id)
    if not version:
        raise HTTPException(status_code=404, detail="CompoundVersion not found")
    active_models = db.scalars(select(ADMETModelRegistry).where(ADMETModelRegistry.is_active.is_(True))).all()
    implemented_models = [model for model in active_models if model.endpoint_name in MODEL_SPECS]
    available_models = [model for model in active_models if model.endpoint_name in MODEL_SPECS and model_files_available(model.endpoint_name)[0]]
    cached = {}
    for model in available_models:
        prediction = db.scalar(
            select(ADMETPrediction).join(ADMETModelRegistry)
            .where(ADMETPrediction.version_id == row_id,
                   ADMETModelRegistry.id == model.id,
                   ADMETModelRegistry.model_version == model.model_version)
            .order_by(ADMETPrediction.created_at.desc())
        )
        if prediction:
            cached[model.id] = prediction
    compound = db.get(Compound, version.compound_row_id)
    measurements = db.scalars(select(ADMETMeasurement).where(ADMETMeasurement.version_id == row_id)).all()
    endpoint_names = {endpoint.id: endpoint.name for endpoint in db.scalars(
        select(ADMETEndpoint).where(ADMETEndpoint.project_id == compound.project_id)
    )}
    if implemented_models and len(available_models) == len(implemented_models) and len(cached) == len(available_models):
        predictions = [_admet_prediction_out(cached[model.id], measurements, endpoint_names) for model in available_models]
        return {
            "type": "Predicted", "run_id": predictions[0]["run_id"], "status": "CACHED",
            "message": "Cached predictions reused for this CompoundVersion and model version.",
            "models_available": len(available_models), "cache_hit": True, "predictions": predictions,
        }

    digest = hashlib.sha256(f"{version.id}|{version.canonical_smiles}|{MODEL_VERSION}".encode()).hexdigest()
    run = ADMETPredictionRun(
        version_id=row_id, inputs_hash=digest, status="RUNNING",
        message="Running endpoint-specific ADMET predictions through Stage 3C.",
    )
    db.add(run); db.flush()
    created, unavailable = [], []
    selected_predictions = dict(cached)
    for model in active_models:
        if model.id in cached:
            continue
        if model.endpoint_name not in MODEL_SPECS:
            continue
        available, reason = model_files_available(model.endpoint_name)
        if not available:
            unavailable.append(f"{model.endpoint_name}: {reason}")
            continue
        try:
            result = predict_endpoint(version.canonical_smiles, model.endpoint_name)
        except Exception as exc:
            unavailable.append(f"{model.endpoint_name}: inference failed ({exc})")
            continue
        if result.get("status") != "COMPLETE":
            unavailable.append(f"{model.endpoint_name}: {result.get('reason', 'model unavailable')}")
            continue
        endpoint = get_or_create_admet_endpoint(db, compound.project_id, model.endpoint_name)
        db.add(endpoint); db.flush()
        domain = result["applicability_domain"]
        output = {
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
        for key in ("probability", "classification", "isoform", "role", "decision_threshold", "liability_summary"):
            if result.get(key) is not None:
                output[key] = result[key]
        if result.get("derived_outputs") is not None:
            output["derived_outputs"] = result["derived_outputs"]
        if result.get("metabolic_stability_assessment") is not None:
            output["metabolic_stability_assessment"] = result["metabolic_stability_assessment"]
        prediction = ADMETPrediction(
            run_id=run.id, endpoint_id=endpoint.id, version_id=row_id, model_id=model.id,
            predicted_value=result["predicted_value"], unit=result["unit"],
            confidence=result["confidence"], applicability_domain=domain["classification"],
            uncertainty=result["uncertainty"], outputs_json=output,
        )
        db.add(prediction); created.append(prediction); selected_predictions[model.id] = prediction
    from datetime import datetime, timezone
    run.completed_at = datetime.now(timezone.utc)
    if selected_predictions and unavailable:
        run.status, run.message = "PARTIAL", "Predictions completed; " + "; ".join(unavailable)
    elif selected_predictions:
        run.status, run.message = "COMPLETE", "Stage 3A/3B/3C endpoint predictions completed." + (
            f" Reused {len(cached)} cached endpoint." if cached else ""
        )
    else:
        run.status, run.message = "MODEL_UNAVAILABLE", "; ".join(unavailable) or "No implemented ADMET model is available."
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
        "unavailable": unavailable,
    }


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
        version=next(v for v in compound.versions if v.version_number==compound.current_version)
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
    current=next(v for v in compound.versions if v.version_number==compound.current_version)
    existing=_experimental_summary(db,current.id,assay_id)
    _,target_fp,target_desc,_=fingerprint_and_descriptors(current.canonical_smiles)
    dataset={"rows":[],"fingerprints":[],"descriptors":[]}
    for other in db.scalars(select(Compound).where(Compound.project_id==assay.project_id)).all():
        version=next(v for v in other.versions if v.version_number==other.current_version)
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
      v=next(v for v in c.versions if v.version_number==c.current_version);e=_experimental_summary(db,v.id,assay_id)
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
        v=next((v for v in c.versions if v.version_number==c.current_version),None); summary=_experimental_summary(db,v.id,assay_id)
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


@app.get("/", response_class=HTMLResponse)
def index():
    with open("frontend/static/index.html") as handle: return handle.read()

app.mount("/static", StaticFiles(directory="frontend/static"), name="static")
