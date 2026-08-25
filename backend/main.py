import json

from fastapi import Depends, FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from sqlalchemy import func, select
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session

from .chemistry import ChemistryError, ENGINE, ENGINE_VERSION, analyze_smiles
from .database import Base, SessionLocal, engine, get_db
from .models import Compound, CompoundVersion, PredictionRun, Project, PropertyCalculation, StructuralAlert, utcnow
from .schemas import CompoundCreate, CompoundUpdate, ProjectCreate, ProjectOut, ProjectUpdate

app = FastAPI(title="AI Drug Optimization Platform", version="0.1.0-stage1")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])


@app.on_event("startup")
def startup():
    if not app.dependency_overrides:
        Base.metadata.create_all(bind=engine)


def _project_out(db: Session, project: Project):
    count = db.scalar(select(func.count(Compound.id)).where(Compound.project_id == project.id)) or 0
    return ProjectOut.model_validate(project).model_copy(update={"compound_count": count})


@app.get("/api/health")
def health():
    return {"status": "ok", "stage": 1, "engine": ENGINE, "engine_version": ENGINE_VERSION}


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
        rows.append({
            "compound": compound.compound_id, "name": compound.name, "row_id": compound.id,
            "MW": p.get("molecular_weight"), "cLogP": p.get("clogp"), "TPSA": p.get("tpsa"), "HBD": p.get("hbd"),
            "HBA": p.get("hba"), "RotB": p.get("rotatable_bonds"), "Fsp3": p.get("fraction_csp3"), "QED": p.get("qed"),
            "svg": version.svg if version else "", "inchikey": version.inchikey if version else "",
        })
    if len(rows) < 2: raise HTTPException(status_code=400, detail="At least two selected compounds must belong to the project")
    metrics = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "Fsp3", "QED"]
    ranges = {metric: {"min": min(r[metric] for r in rows if r[metric] is not None),
                       "max": max(r[metric] for r in rows if r[metric] is not None)} for metric in metrics}
    return {"metrics": metrics, "ranges": ranges, "compounds": rows}


@app.get("/", response_class=HTMLResponse)
def index():
    with open("frontend/static/index.html") as handle: return handle.read()

app.mount("/static", StaticFiles(directory="frontend/static"), name="static")
