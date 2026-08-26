from pathlib import Path

import pytest
from fastapi import HTTPException
from sqlalchemy import create_engine, inspect, text
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend.admet import ADMETEndpoint, ADMETModelRegistry, ADMETPrediction, ADMETPredictionRun, ensure_admet_schema
from backend.database import Base
from backend.main import (
    add_measurement as add_activity_measurement,
    calculate_compound_properties,
    compare,
    create_admet_measurement,
    create_assay,
    create_compound,
    create_experimental_metabolite,
    create_project,
    dashboard_summary,
    get_compound_version_admet,
    get_compound_version_workspace,
    get_project,
    matched_pairs,
    sar_table,
    train_assay_model,
)
from backend.metabolism import MetabolicPredictionRun, ensure_metabolism_schema
from backend.models import ensure_ui_schema
from backend.schemas import CompoundCreate, ProjectCreate


ROOT = Path(__file__).parents[1]


@pytest.fixture()
def db():
    engine = create_engine("sqlite://", connect_args={"check_same_thread": False}, poolclass=StaticPool)
    Base.metadata.create_all(engine)
    ensure_ui_schema(engine)
    ensure_admet_schema(engine)
    ensure_metabolism_schema(engine)
    session = sessionmaker(bind=engine, autoflush=False, autocommit=False)()
    try:
        yield session
    finally:
        session.close()


def add(db, project_id, name, smiles="", calculate=True, compound_id=""):
    return create_compound(
        project_id,
        CompoundCreate(compound_id=compound_id, name=name, smiles=smiles, calculate=calculate),
        db,
    )


def test_non_destructive_ui_schema_migration_backfills_existing_rows():
    engine = create_engine("sqlite://")
    with engine.begin() as connection:
        connection.execute(text("CREATE TABLE projects (id INTEGER PRIMARY KEY, name VARCHAR(200))"))
        connection.execute(text("CREATE TABLE compounds (id INTEGER PRIMARY KEY, project_id INTEGER, name VARCHAR(200))"))
        connection.execute(text("INSERT INTO projects (id, name) VALUES (1, 'Legacy')"))
        connection.execute(text("INSERT INTO compounds (id, project_id, name) VALUES (1, 1, 'Legacy compound')"))
    ensure_ui_schema(engine)
    assert {row["name"] for row in inspect(engine).get_columns("projects")} >= {"molecule_type"}
    assert {row["name"] for row in inspect(engine).get_columns("compounds")} >= {"status"}
    with engine.connect() as connection:
        assert connection.execute(text("SELECT name, molecule_type FROM projects WHERE id=1")).one() == ("Legacy", "Small Molecule")
        assert connection.execute(text("SELECT name, status FROM compounds WHERE id=1")).one() == ("Legacy compound", "CALCULATED")


def test_project_molecule_type_and_draft_name_only_save(db):
    project = create_project(ProjectCreate(name="Simple Project", target="EGFR", molecule_type="Small Molecule"), db)
    draft = add(db, project.id, "HIT-001")
    assert project.target == "EGFR" and project.molecule_type == "Small Molecule"
    assert draft["name"] == "HIT-001" and draft["compound_id"] == "HIT-001"
    assert draft["status"] == "DRAFT" and draft["version"] is None
    assert get_project(project.id, db)["compounds"][0]["status"] == "DRAFT"


def test_save_structure_without_calculation_then_calculate_same_version(db):
    project = create_project(ProjectCreate(name="Save Calculate", target="EGFR"), db)
    saved = add(db, project.id, "Lead-023", "CCO", calculate=False)
    assert saved["status"] == "STRUCTURE_READY"
    assert saved["version"]["canonical_smiles"] == "CCO"
    assert saved["version"]["properties"] == {}
    assert saved["version"]["calculated"] is False
    version_id = saved["version"]["id"]
    calculated = calculate_compound_properties(saved["row_id"], db)
    assert calculated["status"] == "CALCULATED"
    assert calculated["version"]["id"] == version_id
    assert calculated["version"]["calculated"] is True
    assert calculated["version"]["properties"]["molecular_weight"] > 40


def test_invalid_structure_not_saved_as_calculated_and_peptide_isolated(db):
    project = create_project(ProjectCreate(name="Validation", target="T"), db)
    with pytest.raises(HTTPException) as invalid:
        add(db, project.id, "Bad", "invalid", calculate=False)
    assert invalid.value.status_code == 400
    peptide = create_project(ProjectCreate(name="Peptide", target="Target", molecule_type="Peptide"), db)
    draft = add(db, peptide.id, "PEP-001")
    assert draft["status"] == "DRAFT"
    with pytest.raises(HTTPException, match="small molecules only"):
        add(db, peptide.id, "PEP-002", "NCC(=O)O", calculate=False)


def test_compound_workspace_isolates_project_compound_and_version_data(db):
    project_a = create_project(ProjectCreate(name="Project A", target="A"), db)
    project_b = create_project(ProjectCreate(name="Project B", target="B"), db)
    a1 = add(db, project_a.id, "A1", "CCO")
    a2 = add(db, project_a.id, "A2", "CCN")
    b1 = add(db, project_b.id, "B1", "CCC")
    for compound, project_id, value in ((a1, project_a.id, 11), (a2, project_a.id, 22), (b1, project_b.id, 33)):
        create_admet_measurement(project_id, {
            "version_id": compound["version"]["id"], "endpoint": "Solubility",
            "value": value, "unit": "µM", "source": compound["name"],
        }, db)
        create_experimental_metabolite(project_id, {
            "version_id": compound["version"]["id"], "transformation": f"marker-{compound['name']}",
            "source": compound["name"],
        }, db)
    endpoint = db.query(ADMETEndpoint).filter_by(project_id=project_a.id, name="Solubility").one()
    model = db.query(ADMETModelRegistry).filter_by(endpoint_name="Solubility").one()
    a2_run = ADMETPredictionRun(version_id=a2["version"]["id"], inputs_hash="a2-only", status="COMPLETE", message="A2 marker")
    db.add(a2_run); db.flush()
    db.add(ADMETPrediction(
        run_id=a2_run.id, endpoint_id=endpoint.id, version_id=a2["version"]["id"], model_id=model.id,
        predicted_value=-2.2, unit=model.output_unit, confidence="LOW", applicability_domain="BORDERLINE",
        outputs_json={"marker": "A2 prediction"},
    ))
    db.add(MetabolicPredictionRun(
        version_id=a2["version"]["id"], inputs_hash="a2-metabolism", engine_name="fixture",
        engine_version="1", status="COMPLETE", message="A2 metabolism marker",
    ))
    db.commit()
    assay = create_assay(project_a.id, {"name": "A IC50", "measurement_type": "IC50", "unit": "nM"}, db)
    add_activity_measurement(assay["id"], {"version_id": a1["version"]["id"], "value": 7, "unit": "nM", "source": "A1-only"}, db)

    workspace = get_compound_version_workspace(a1["version"]["id"], db)
    assert workspace["scope"] == {"project_id": project_a.id, "compound_id": a1["row_id"], "version_id": a1["version"]["id"]}
    assert {row["version_id"] for row in workspace["admet"]["measurements"]} == {a1["version"]["id"]}
    assert [row["value"] for row in workspace["admet"]["measurements"]] == [11.0]
    assert {row["version_id"] for row in workspace["activity"]["measurements"]} == {a1["version"]["id"]}
    assert {row["version_id"] for row in workspace["metabolism"]["experimental_metabolites"]} == {a1["version"]["id"]}
    assert workspace["admet"]["predictions"] == []
    assert workspace["metabolism"]["runs"] == []
    assert all(row["version_id"] == a1["version"]["id"] for row in workspace["prediction_audit"])
    exact_admet = get_compound_version_admet(a2["version"]["id"], db)
    assert [row["value"] for row in exact_admet["measurements"]] == [22.0]
    assert {row["version_id"] for row in exact_admet["predictions"]} == {a2["version"]["id"]}


def test_activity_measurement_rejects_cross_project_version(db):
    first = create_project(ProjectCreate(name="Assay owner", target="A"), db)
    second = create_project(ProjectCreate(name="Other project", target="B"), db)
    assay = create_assay(first.id, {"name": "IC50", "measurement_type": "IC50", "unit": "nM"}, db)
    foreign = add(db, second.id, "FOREIGN", "CCO")
    with pytest.raises(HTTPException) as error:
        add_activity_measurement(assay["id"], {"version_id": foreign["version"]["id"], "value": 5, "unit": "nM"}, db)
    assert error.value.status_code == 404


def test_qualitative_experimental_measurements_keep_endpoint_roles_separate(db):
    project = create_project(ProjectCreate(name="Qualitative", target="CYP"), db)
    compound = add(db, project.id, "Q-1", "CCO")
    version_id = compound["version"]["id"]
    inhibitor = create_admet_measurement(project.id, {
        "version_id": version_id, "endpoint": "CYP3A4 inhibitor", "qualitative_value": "POSITIVE",
        "unit": "classification", "method": "classification · Inhibition",
    }, db)
    substrate = create_admet_measurement(project.id, {
        "version_id": version_id, "endpoint": "CYP3A4 substrate", "qualitative_value": "NEGATIVE",
        "unit": "classification", "method": "classification · Substrate",
    }, db)
    assert inhibitor["qualitative_value"] == "POSITIVE"
    assert substrate["qualitative_value"] == "NEGATIVE"
    listing = get_compound_version_admet(version_id, db)
    names = {row["id"]: row["name"] for row in listing["endpoints"]}
    assert {names[row["endpoint_id"]] for row in listing["measurements"]} == {"CYP3A4 inhibitor", "CYP3A4 substrate"}


def test_compare_contains_only_selected_same_project_compounds_and_handles_uncalculated(db):
    project = create_project(ProjectCreate(name="Compare isolation", target="T"), db)
    first = add(db, project.id, "C1", "CCO", calculate=False)
    second = add(db, project.id, "C2", "CCN", calculate=False)
    unselected = add(db, project.id, "C3", "CCC")
    result = compare(project.id, ids=f"{first['row_id']},{second['row_id']}", db=db)
    assert {row["name"] for row in result["compounds"]} == {"C1", "C2"}
    assert unselected["row_id"] not in {row["row_id"] for row in result["compounds"]}
    assert result["ranges"]["MW"] == {"min": None, "max": None}


def test_dashboard_summary_is_project_isolated_and_uses_current_version_status(db):
    project_a = create_project(ProjectCreate(name="Dashboard A", target="EGFR"), db)
    project_b = create_project(ProjectCreate(name="Dashboard B", target="BRAF"), db)
    compound_a = add(db, project_a.id, "A-1", "CCO")
    add(db, project_b.id, "B-1", "CCN", calculate=False)
    assay = create_assay(project_a.id, {"name": "A IC50", "measurement_type": "IC50", "unit": "nM"}, db)
    add_activity_measurement(assay["id"], {"version_id": compound_a["version"]["id"], "value": 9, "unit": "nM"}, db)
    create_admet_measurement(project_a.id, {
        "version_id": compound_a["version"]["id"], "endpoint": "Solubility", "value": 15,
        "unit": "µM", "source": "Dashboard A only",
    }, db)

    result = dashboard_summary(db)
    rows = {row["id"]: row for row in result["projects"]}
    assert result["totals"] == {"projects": 2, "compounds": 2}
    assert rows[project_a.id]["experimental_activity_count"] == 1
    assert rows[project_a.id]["experimental_admet_count"] == 1
    assert rows[project_a.id]["workflow"]["Properties"] == "READY"
    assert rows[project_a.id]["compounds"][0]["activity"] == "EXPERIMENTAL"
    assert rows[project_b.id]["experimental_activity_count"] == 0
    assert rows[project_b.id]["experimental_admet_count"] == 0
    assert rows[project_b.id]["compounds"][0]["properties"] == "NOT_RUN"
    assert "Dashboard A only" not in str(rows[project_b.id])


def test_drafts_do_not_break_stage2_project_analysis(db):
    project = create_project(ProjectCreate(name="Draft safe Stage 2", target="T"), db)
    add(db, project.id, "Draft only")
    structured = add(db, project.id, "Structured", "CCO")
    assay = create_assay(project.id, {"name": "IC50", "measurement_type": "IC50", "unit": "nM"}, db)
    add_activity_measurement(assay["id"], {"version_id": structured["version"]["id"], "value": 10, "unit": "nM"}, db)
    assert [row["compound"] for row in sar_table(project.id, assay["id"], db)["compounds"]] == ["STRUCTURED"]
    assert train_assay_model(assay["id"], db)["policy"]["N"] == 1
    assert matched_pairs(project.id, assay["id"], db=db)["pairs"] == []


def test_ui_contract_has_compound_editor_selector_comparison_and_no_pk_tab():
    source = (ROOT / "frontend/static/app.js").read_text()
    for phrase in (
        "Draw Chemical Structure", "Or Enter SMILES", "Save Compound", "Save & Calculate",
        "Calculate Properties", "Add Experimental Data", "What experimental data do you want to add?",
        "Microsomal Stability", "Caco-2 Permeability", "Plasma Protein Binding (PPB)",
        "Compare Selected", "Strict scope", "MODEL UNAVAILABLE",
    ):
        assert phrase in source
    assert "const tabs=['overview','properties','activity','admet','metabolism','optimization','history']" in source
    assert "const tabs=['overview','properties','activity','admet','metabolism','optimization','history','pk']" not in source.lower()
    assert "PK / DMPK" in source and "PLANNED" in source
    assert (ROOT / "frontend/static/ketcher/standalone/index.html").exists()
    assert "Ketcher 3.5.0" in (ROOT / "frontend/static/ketcher/NOTICE.md").read_text()


def test_main_dashboard_contract_and_responsive_layout():
    source = (ROOT / "frontend/static/app.js").read_text()
    styles = (ROOT / "frontend/static/app.css").read_text()
    for phrase in (
        "PLATFORM OVERVIEW", "Drug Optimization Platform", "Available Scientific Modules",
        "Structure & Chemistry", "Activity & SAR", "CYP & Transporters", "Safety / Toxicology",
        "PK / DMPK", "PLANNED", "Create New Project", "Typical Workflow", "Default Workspace Settings",
        "Current Project Status", "Compound Status", "WORKFLOW STATUS", "Continue Current Project", "Open Project",
    ):
        assert phrase in source
    assert "useState('dashboard')" in source
    assert "function MainDashboard()" in source
    assert "function ProjectWorkspace()" in source
    assert "projectTab==='dashboard'?MainDashboard():ProjectWorkspace()" in source
    assert ".module-grid" in styles and ".dashboard-project-grid" in styles and ".global-nav" in styles
    assert ".sidebar.open .sidebar-body" in styles and ".project-overview-grid" in styles
    assert "@media(max-width:900px)" in styles and "@media(max-width:620px)" in styles


def test_global_sidebar_contains_only_six_top_level_workflow_items():
    source = (ROOT / "frontend/static/app.js").read_text()
    sidebar = source[source.index("const sidebarItems="):source.index("const sidebar=e('aside'")]
    for label in ("Dashboard", "New Project", "Projects", "Optimization", "Settings", "Help"):
        assert f"['{label}'" in sidebar
    for scientific_item in (
        "Structure", "Properties", "Activity / SAR", "Absorption", "Distribution", "Metabolism",
        "CYP", "Transporters", "Toxicology", "PK / DMPK",
    ):
        assert f"['{scientific_item}'" not in sidebar
    assert "sidebarGroups" not in source
    assert "aria-label':'Primary navigation'" in source
    assert "Where scientific functions live" in source


def test_project_delete_ui_requires_typed_confirmation_and_manual_bulk_selection():
    source = (ROOT / "frontend/static/app.js").read_text()
    styles = (ROOT / "frontend/static/app.css").read_text()
    for phrase in (
        "Delete Selected", "Delete Project…", "This action permanently deletes all project-linked data.",
        "Type ", " to confirm", "Delete Project Permanently", "Delete Selected Projects Permanently",
        "confirmation_name", "projectSelection", "experimental_activity_count", "prediction_count",
    ):
        assert phrase in source
    assert "useState([]),[deleteProjects" in source
    assert "disabled:deleteBusy||!deleteNamesMatch" in source
    assert ".project-delete-modal" in styles and ".delete-count-grid" in styles
