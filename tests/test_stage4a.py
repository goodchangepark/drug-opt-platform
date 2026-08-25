from pathlib import Path

import pytest
from fastapi import HTTPException
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend.activity_models import MatchedMolecularPair
from backend.admet import ensure_admet_schema
from backend.database import Base
from backend.main import (
    add_measurement, create_admet_measurement, create_assay, create_compound,
    create_optimization_run, create_project, list_optimization_runs,
    update_optimization_overrides,
)
from backend.metabolism import MetabolicPredictionRun, MetabolicSoftSpot, ensure_metabolism_schema
from backend.optimization import OptimizationRun, ensure_optimization_schema
from backend.optimization_engine import EVIDENCE_HIERARCHY, TRANSFORMATION_LIBRARY
from backend.schemas import CompoundCreate, ProjectCreate


ROOT = Path(__file__).parents[1]


@pytest.fixture()
def db():
    engine = create_engine("sqlite://", connect_args={"check_same_thread": False}, poolclass=StaticPool)
    Base.metadata.create_all(bind=engine)
    ensure_admet_schema(engine)
    ensure_metabolism_schema(engine)
    ensure_optimization_schema(engine)
    session = sessionmaker(bind=engine, autoflush=False, autocommit=False)()
    try:
        yield session
    finally:
        session.close()


def setup_parent(db, project_name="Stage 4A", compound_id="OPT-001", smiles="CCN(CC)Cc1ccc(OC)cc1"):
    project = create_project(ProjectCreate(name=project_name), db)
    compound = create_compound(project.id, CompoundCreate(compound_id=compound_id, smiles=smiles), db)
    return project, compound


def analyze(db, project, compound, **payload):
    data = {
        "parent_version_id": compound["version"]["id"],
        "objectives": ["Balanced optimization"],
        "constraints": {},
        **payload,
    }
    return create_optimization_run(project.id, data, db)


def add_hlm(db, project_id, version_id, value):
    return create_admet_measurement(project_id, {
        "version_id": version_id, "endpoint": "HLM intrinsic clearance", "species": "Human",
        "matrix": "HLM", "value": value, "unit": "log10(mL/min/kg)",
        "method": "human liver microsome", "source": "Public experimental reference",
    }, db)


def add_project_pair(db, project, parent, other_smiles="CCN(CC)Cc1ccncc1", parent_nm=10, other_nm=12, cliff=False):
    other = create_compound(project.id, CompoundCreate(compound_id="OPT-002", smiles=other_smiles), db)
    assay = create_assay(project.id, {"name": "Target IC50", "measurement_type": "IC50", "unit": "nM"}, db)
    add_measurement(assay["id"], {"version_id": parent["version"]["id"], "value": parent_nm, "unit": "nM"}, db)
    add_measurement(assay["id"], {"version_id": other["version"]["id"], "value": other_nm, "unit": "nM"}, db)
    pair = MatchedMolecularPair(
        assay_id=assay["id"], version_a_id=parent["version"]["id"],
        version_b_id=other["version"]["id"], similarity=0.8,
        delta_pactivity=-2.0 if cliff else -0.08,
        transformation_smiles=parent["version"]["canonical_smiles"] + ">>" + other["version"]["canonical_smiles"],
        is_cliff=cliff, provenance_json={"source": "Project experimental activity"},
    )
    db.add(pair); db.commit()
    return other, assay, pair


def test_architecture_objectives_constraints_and_no_analog_generation(db):
    project, parent = setup_parent(db)
    result = analyze(db, project, parent, objectives=["Improve solubility", "Reduce hERG liability"], constraints={
        "clogp_max": 4, "tpsa_min": 40, "tpsa_max": 100, "mw_max": 550,
        "similarity_min": 0.6, "herg_do_not_increase": True,
    })
    assert result["status"] == "COMPLETE" and result["analog_generation"] == "NOT_PERFORMED"
    assert result["objectives"] == ["Improve solubility", "Reduce hERG liability"]
    assert result["constraints"]["similarity_min"] == 0.6
    assert db.get(OptimizationRun, result["id"]).parent_version_id == parent["version"]["id"]


def test_explicit_evidence_hierarchy_and_experimental_hlm_precedence(db):
    project, parent = setup_parent(db)
    add_hlm(db, project.id, parent["version"]["id"], 2.2)
    result = analyze(db, project, parent, objectives=["Improve metabolic stability"])
    assert [row["type"] for row in EVIDENCE_HIERARCHY] == [
        "Experimental", "Project-specific validated model/SAR", "External validated quantitative model",
        "External classification model", "Rule-based hypothesis",
    ]
    liability = next(row for row in result["liabilities"] if row["liability_type"] == "metabolic_stability")
    assert liability["evidence_type"] == "Experimental"
    assert result["evidence"]["admet"]["HLM intrinsic clearance"]["predicted"] is None
    assert any(row["purpose"] == "Metabolism" for row in result["recommended_transformations"])


def test_potency_constraint_uses_selected_assay(db):
    project, parent = setup_parent(db)
    assay = create_assay(project.id, {"name": "Potency", "measurement_type": "IC50", "unit": "nM"}, db)
    add_measurement(assay["id"], {"version_id": parent["version"]["id"], "value": 100, "unit": "nM"}, db)
    result = analyze(db, project, parent, assay_id=assay["id"], objectives=["Improve potency"], constraints={"potency_max_nm": 30, "do_not_worsen_fold": 2})
    potency = next(row for row in result["liabilities"] if row["id"] == "LIAB_POTENCY")
    assert potency["evidence_type"] == "Experimental" and "30" in potency["rationale"]
    with pytest.raises(HTTPException):
        analyze(db, project, parent, objectives=["Improve potency"], constraints={"potency_max_nm": 30})


def test_project_mmp_has_priority_and_modifiable_region(db):
    project, parent = setup_parent(db)
    _, assay, pair = add_project_pair(db, project, parent)
    result = analyze(db, project, parent, assay_id=assay["id"], objectives=["Balanced optimization"])
    first = result["recommended_transformations"][0]
    assert first["id"] == f"MMP_PROJECT_OBSERVED_{pair.id}" and first["confidence"] == "HIGH"
    assert any(row["id"] == f"MMP_{pair.id}" and row["risk"] == "LOW" for row in result["modifiable_regions"])


def test_project_mmp_records_experimental_hlm_improvement(db):
    project, parent = setup_parent(db)
    other, assay, pair = add_project_pair(db, project, parent)
    add_hlm(db, project.id, parent["version"]["id"], 2.2)
    add_hlm(db, project.id, other["version"]["id"], 1.4)
    result = analyze(db, project, parent, assay_id=assay["id"], objectives=["Improve metabolic stability"])
    evidence_pair = next(row for row in result["evidence"]["activity"]["mmp"] if row["pair_id"] == pair.id)
    effect = next(row for row in evidence_pair["endpoint_effects"] if row["endpoint"] == "HLM intrinsic clearance")
    assert effect["direction"] == "IMPROVED" and effect["delta"] == pytest.approx(-0.8)
    strategy = next(row for row in result["recommended_transformations"] if row["id"] == f"MMP_PROJECT_OBSERVED_{pair.id}")
    assert strategy["score"] == 100 and any("HLM intrinsic clearance IMPROVED" in item for item in strategy["evidence"])


def test_activity_cliff_protects_parent_region(db):
    project, parent = setup_parent(db)
    _, assay, pair = add_project_pair(db, project, parent, other_smiles="CCN(CC)Cc1ccc(F)cc1", parent_nm=5, other_nm=500, cliff=True)
    result = analyze(db, project, parent, assay_id=assay["id"], objectives=["Improve potency"])
    region = next(row for row in result["protected_regions"] if row["id"] == f"CLIFF_{pair.id}")
    assert region["status"] == "HIGH-RISK TO MODIFY" and region["atom_indices"]
    assert region["source"] == "Project experimental SAR/activity cliff"


def test_soft_spot_becomes_modifiable_only_as_hypothesis(db):
    project, parent = setup_parent(db, smiles="CCc1ccccc1")
    run = MetabolicPredictionRun(version_id=parent["version"]["id"], inputs_hash="fixture", status="COMPLETE", message="fixture", engine_name="Rules", engine_version="test")
    db.add(run); db.flush()
    db.add(MetabolicSoftSpot(
        run_id=run.id, version_id=parent["version"]["id"], rank=1, atom_index=1,
        atom_environment="benzylic carbon", transformation="benzylic oxidation", phase="Phase I",
        cyp_isoform="CYP isoform not assigned", model_evidence_json={}, rule_evidence_json={},
        score=None, confidence="MEDIUM", provenance_json={"source": "rule fixture"},
    ))
    db.commit()
    add_hlm(db, project.id, parent["version"]["id"], 2.3)
    result = analyze(db, project, parent, objectives=["Improve metabolic stability"])
    assert any(row["id"].startswith("SOFT_SPOT_") and row["atom_indices"] == [1] for row in result["modifiable_regions"])
    assert any(row["liability_type"] == "metabolic_soft_spot" for row in result["liabilities"])


def test_low_confidence_classification_alone_is_supporting_only(db, monkeypatch):
    project, parent = setup_parent(db, smiles="CCO")
    evidence = {
        "activity": {}, "properties": {"clogp": {"value": 0.0}}, "structural_alerts": [],
        "metabolism": {"soft_spots": []}, "admet": {"hERG liability": {"preferred": {
            "type": "External classification model", "classification": "BLOCKER", "confidence": "LOW",
        }}},
    }
    from backend.optimization_engine import identify_liabilities
    rows = identify_liabilities(evidence, ["Reduce hERG liability"], {}, {})
    assert rows[0]["actionability"] == "SUPPORTING_ONLY" and rows[0]["score"] <= 42
    result = analyze(db, project, parent, objectives=["Balanced optimization"])
    assert result["status"] == "COMPLETE"


def test_manual_overrides_are_saved_and_rerank_without_analog(db):
    project, parent = setup_parent(db, smiles="CCCCc1ccccc1")
    result = analyze(db, project, parent, objectives=["Improve solubility"], constraints={"clogp_max": 1})
    available = result["recommended_transformations"]
    assert available
    selected_id = available[-1]["id"]
    updated = update_optimization_overrides(result["id"], {
        "protect_atoms": [[0]], "allow_atoms": [[1]],
        "prioritize_transformations": [selected_id], "exclude_transformations": [available[0]["id"]],
    }, db)
    assert updated["manual_overrides"]["protect_atoms"] == [[0]]
    assert any(row["id"] == selected_id and row["manual_priority"] for row in updated["recommended_transformations"])
    assert all(row["id"] != available[0]["id"] for row in updated["recommended_transformations"])
    assert updated["analog_generation"] == "NOT_PERFORMED"


def test_project_and_compound_version_isolation(db):
    first_project, first = setup_parent(db, "First", "FIRST")
    second_project, second = setup_parent(db, "Second", "SECOND", "CCO")
    result = analyze(db, first_project, first)
    assert len(list_optimization_runs(first_project.id, first["version"]["id"], db)["runs"]) == 1
    assert list_optimization_runs(second_project.id, second["version"]["id"], db)["runs"] == []
    with pytest.raises(HTTPException):
        analyze(db, first_project, second)
    assert result["evidence"]["parent"]["version_id"] == first["version"]["id"]


def test_transformation_library_provenance_and_required_families():
    purposes = {row["purpose"] for row in TRANSFORMATION_LIBRARY}
    ids = {row["id"] for row in TRANSFORMATION_LIBRARY}
    assert {"Metabolism", "Lipophilicity", "Solubility", "Safety"} <= purposes
    assert {"MET_F_FLUORINATION", "MET_N_DEALK_BLOCK", "MET_O_DEALK_BLOCK", "POT_BIOISOSTERE", "POT_LINKER_REPLACE", "SAFE_ALERT_REMOVAL"} <= ids
    assert all(row["reaction_smarts"] and row["source"] and row["version"] and row["possible_risk"] for row in TRANSFORMATION_LIBRARY)


def test_ui_contract_has_optimization_workflow_and_no_analog_claim():
    source = (ROOT / "frontend/static/app.js").read_text()
    for text in ("Optimization", "Current profile", "Main liabilities", "Protected regions", "Modifiable regions", "Recommended transformations"):
        assert text in source
    assert "no analog" in source.lower()


def test_public_acceptance_examples_include_expected_strategy_direction():
    from scripts.validate_stage4a_engine import validate
    result = validate()
    assert result["passed"] == result["total"] == 3
    assert all(row["observed_rank"] is not None and row["analog_generation"] == "NOT_PERFORMED" for row in result["results"])
    assert "not an independent" in result["scope"]
