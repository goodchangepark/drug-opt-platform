from pathlib import Path

import pytest
from rdkit import Chem
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend.admet import ensure_admet_schema
from backend.database import Base
from backend.main import (
    candidate_decision, create_admet_measurement, create_compound,
    create_optimization_run, create_project,
)
from backend.metabolism import ensure_metabolism_schema
from backend.optimization import OptimizationRun, ensure_optimization_schema
from backend.optimization_engine import TRANSFORMATION_LIBRARY
from backend.proposal import (
    CandidatePredictionSnapshot, CandidateRanking, CandidateRejectionReason,
    CandidateTransformation, OptimizationCandidate, OptimizationProposalRun,
    ensure_proposal_schema,
)
from backend import proposal_engine
from backend.proposal_engine import (
    EXECUTABLE_TRANSFORMATIONS, STRATEGY_ONLY_TRANSFORMATIONS,
    _cheap_constraints, chemical_validation, execute_proposal_run,
    execute_strategy, pareto_fronts, process_user_candidate, synthetic_feasibility,
)
from backend.schemas import CompoundCreate, ProjectCreate


ROOT = Path(__file__).parents[1]


@pytest.fixture()
def db():
    engine = create_engine("sqlite://", connect_args={"check_same_thread": False}, poolclass=StaticPool)
    Base.metadata.create_all(bind=engine)
    ensure_admet_schema(engine); ensure_metabolism_schema(engine); ensure_optimization_schema(engine); ensure_proposal_schema(engine)
    session = sessionmaker(bind=engine, autoflush=False, autocommit=False)()
    try:
        yield session
    finally:
        session.close()


def parent_and_strategy(db, smiles="CCN(CC)C(=O)c1c(C)cccc1C", constraints=None):
    project = create_project(ProjectCreate(name="Stage 4B"), db)
    parent = create_compound(project.id, CompoundCreate(compound_id="PARENT", smiles=smiles), db)
    create_admet_measurement(project.id, {
        "version_id": parent["version"]["id"], "endpoint": "HLM intrinsic clearance",
        "species": "Human", "matrix": "HLM", "value": 2.2, "unit": "log10(mL/min/kg)",
        "method": "public fixture", "source": "Directional test fixture",
    }, db)
    optimization = create_optimization_run(project.id, {
        "parent_version_id": parent["version"]["id"], "objectives": ["Improve metabolic stability"],
        "constraints": {"similarity_min": 0.45, **(constraints or {})},
    }, db)
    return project, parent, db.get(OptimizationRun, optimization["id"])


def proposal(db, project, parent, optimization, max_raw=12, doubles=False, hard=None):
    row = OptimizationProposalRun(
        project_id=project.id, optimization_run_id=optimization.id,
        parent_version_id=parent["version"]["id"], status="PENDING",
        hard_constraints_json=hard or {}, settings_json={"max_raw_candidates": max_raw, "allow_double_transforms": doubles},
    )
    db.add(row); db.commit(); db.refresh(row)
    return row


def fake_admet(db, candidate, project_id):
    result = {
        "Solubility": {"status": "COMPLETE", "scored": True, "predicted_value": -3.0, "unit": "log10(mol/L)", "confidence": "MEDIUM", "applicability_domain": "IN_DOMAIN", "record_type": "Predicted"},
        "Permeability": {"status": "COMPLETE", "scored": True, "predicted_value": -5.0, "unit": "log10(cm/s)", "confidence": "MEDIUM", "applicability_domain": "IN_DOMAIN", "record_type": "Predicted"},
        "HLM intrinsic clearance": {"status": "COMPLETE", "scored": True, "predicted_value": 1.2, "unit": "log10(mL/min/kg)", "confidence": "LOW", "applicability_domain": "BORDERLINE", "record_type": "Predicted"},
        "hERG liability": {"status": "COMPLETE", "scored": True, "predicted_value": 0.3, "classification": "NON_BLOCKER", "unit": "probability", "confidence": "LOW", "applicability_domain": "IN_DOMAIN", "record_type": "Predicted"},
        "BSEP inhibitor": {"status": "MODEL_UNAVAILABLE", "reason": "No public checkpoint", "scored": False},
    }
    return result


def fake_soft(smiles, context=None, max_spots=12):
    return {"spots": [{"rank": 1, "atom_index": 0, "transformation": "Aliphatic oxidation", "confidence": "LOW"}], "metabolites": [], "engine": "fixture"}


def run_fast(db, monkeypatch, max_raw=12, doubles=False, hard=None):
    project, parent, optimization = parent_and_strategy(db)
    row = proposal(db, project, parent, optimization, max_raw=max_raw, doubles=doubles, hard=hard)
    monkeypatch.setattr(proposal_engine, "rescore_admet", fake_admet)
    monkeypatch.setattr(proposal_engine, "predict_soft_spots", fake_soft)
    execute_proposal_run(row.id, session=db); db.refresh(row)
    return project, parent, optimization, row


def test_proposal_schema_and_required_entities(db):
    tables = set(inspect(db.bind).get_table_names())
    assert {"optimization_proposal_runs", "optimization_candidates", "candidate_transformations", "candidate_prediction_snapshots", "candidate_rankings", "candidate_rejection_reasons"} <= tables
    assert all(entity.__tablename__ in tables for entity in (OptimizationProposalRun, OptimizationCandidate, CandidateTransformation, CandidatePredictionSnapshot, CandidateRanking, CandidateRejectionReason))


def test_executable_and_strategy_only_are_explicit_and_complete():
    library_ids = {row["id"] for row in TRANSFORMATION_LIBRARY}
    assert set(EXECUTABLE_TRANSFORMATIONS).isdisjoint(STRATEGY_ONLY_TRANSFORMATIONS)
    assert set(EXECUTABLE_TRANSFORMATIONS) | set(STRATEGY_ONLY_TRANSFORMATIONS) == library_ids
    assert {"MET_F_FLUORINATION", "MET_METHYL_REMOVAL", "LIPO_ALKYL_REDUCTION", "LIPO_PHENYL_HETEROARYL", "POT_LINKER_REPLACE"} <= set(EXECUTABLE_TRANSFORMATIONS)


def test_targeted_fluorination_generates_sanitized_single_change():
    rule = next(row for row in TRANSFORMATION_LIBRARY if row["id"] == "MET_F_FLUORINATION")
    products = execute_strategy("CCc1ccccc1", rule, allowed_atoms=list(range(8)))
    assert products and all(product is not None for product, _ in products)
    assert all("F" in Chem.MolToSmiles(product) for product, _ in products)


def test_protected_region_modification_rejected():
    initial = chemical_validation("CCc1ccccc1", "CC(F)c1ccccc1", set())
    assert initial["valid"] and initial["changed_parent_atoms"]
    blocked = chemical_validation("CCc1ccccc1", "CC(F)c1ccccc1", {initial["changed_parent_atoms"][0]})
    assert not blocked["valid"] and blocked["code"] == "PROTECTED_REGION_MODIFIED"


def test_invalid_valence_and_fragmentation_rejected():
    invalid = chemical_validation("CCO", "C[CH5]O", set())
    fragmented = chemical_validation("CCO", "CCO.Cl", set())
    assert not invalid["valid"] and invalid["code"] == "INVALID_VALENCE_OR_SANITIZATION"
    assert not fragmented["valid"] and fragmented["code"] == "FRAGMENTED_STRUCTURE"


def test_stereochemistry_loss_rejected():
    result = chemical_validation("C[C@H](O)F", "CC(O)F", set())
    assert not result["valid"] and result["code"] == "STEREOCHEMISTRY_LOSS"
    inverted = chemical_validation("C[C@H](O)F", "C[C@@H](O)F", set())
    assert not inverted["valid"] and inverted["code"] == "STEREOCHEMISTRY_LOSS"


def test_synthetic_feasibility_is_surrogate_not_probability():
    result = synthetic_feasibility("CCO", "CCOc1ccccc1")
    assert result["classification"] in {"LOW SYNTHETIC COMPLEXITY", "MODERATE SYNTHETIC COMPLEXITY", "HIGH SYNTHETIC COMPLEXITY"}
    assert result["not_synthesis_success_probability"] and 1 <= result["sa_score"] <= 10


def test_duplicate_analog_gets_persisted_rejection(db):
    project, parent, optimization = parent_and_strategy(db, smiles="CCc1ccccc1")
    row = proposal(db, project, parent, optimization)
    strategy = {"id": "TEST", "name": "test", "reaction_smarts": "", "version": "1", "source": "test", "purpose": "test", "expected_effect": "test", "evidence": ["test"]}
    first = proposal_engine._new_candidate(db, row, optimization, 1, Chem.MolFromSmiles("CC(F)c1ccccc1"), strategy, [1])
    second = proposal_engine._new_candidate(db, row, optimization, 2, Chem.MolFromSmiles("CC(F)c1ccccc1"), strategy, [1])
    parent_analysis = __import__("backend.chemistry", fromlist=["analyze_smiles"]).analyze_smiles(parent["version"]["canonical_smiles"])
    seen = set()
    assert _cheap_constraints(db, first, parent_analysis, {"similarity_min": 0}, set(), seen)
    assert not _cheap_constraints(db, second, parent_analysis, {"similarity_min": 0}, set(), seen)
    assert second.rejection_reasons[0].code == "DUPLICATE_ANALOG"


@pytest.mark.parametrize("hard,code", [
    ({"mw_max": 50}, "EXCESSIVE_MOLECULAR_WEIGHT"),
    ({"similarity_min": 0.99}, "LOW_PARENT_SIMILARITY"),
])
def test_mw_and_similarity_hard_gates(db, hard, code):
    project, parent, optimization = parent_and_strategy(db, smiles="CCc1ccccc1")
    row = proposal(db, project, parent, optimization)
    strategy = {"id": "TEST", "name": "test", "reaction_smarts": "", "version": "1", "source": "test", "purpose": "test", "expected_effect": "test", "evidence": ["test"]}
    candidate = proposal_engine._new_candidate(db, row, optimization, 1, Chem.MolFromSmiles("CC(F)c1ccccc1"), strategy, [1])
    parent_analysis = __import__("backend.chemistry", fromlist=["analyze_smiles"]).analyze_smiles(parent["version"]["canonical_smiles"])
    assert not _cheap_constraints(db, candidate, parent_analysis, hard, set(), set())
    assert candidate.rejection_reasons[0].code == code


def test_new_structural_alert_hard_gate(db):
    project, parent, optimization = parent_and_strategy(db, smiles="CCc1ccccc1")
    row = proposal(db, project, parent, optimization)
    strategy = {"id": "TEST", "name": "test", "reaction_smarts": "", "version": "1", "source": "test", "purpose": "test", "expected_effect": "test", "evidence": ["test"]}
    candidate = proposal_engine._new_candidate(db, row, optimization, 1, Chem.MolFromSmiles("CCc1ccc([N+](=O)[O-])cc1"), strategy, [4])
    parent_analysis = __import__("backend.chemistry", fromlist=["analyze_smiles"]).analyze_smiles(parent["version"]["canonical_smiles"])
    assert not _cheap_constraints(db, candidate, parent_analysis, {"similarity_min": 0, "no_new_structural_alert": True}, set(), set())
    assert candidate.rejection_reasons[0].code == "NEW_STRUCTURAL_ALERT"


def test_staged_pipeline_rescores_only_survivors_and_records_unavailable(db, monkeypatch):
    _, _, _, row = run_fast(db, monkeypatch, max_raw=8)
    assert row.status == "COMPLETED" and row.raw_candidate_count >= 1
    accepted = [item for item in row.candidates if item.status in {"ACCEPTED", "TOP_10"}]
    assert accepted and all(item.stage1_json and item.admet_json and item.soft_spot_json for item in accepted)
    assert all(not item.admet_json["BSEP inhibitor"]["scored"] for item in accepted)
    assert row.summary_json["llm_used"] is False and row.summary_json["pk_run"] is False


def test_out_of_domain_activity_is_penalized_not_automatic_rejection(db, monkeypatch):
    monkeypatch.setattr(proposal_engine, "rescore_admet", fake_admet)
    monkeypatch.setattr(proposal_engine, "predict_soft_spots", fake_soft)
    monkeypatch.setattr(proposal_engine, "predict_candidate_activity", lambda *args, **kwargs: {"status": "COMPLETE", "record_type": "Predicted", "value_nm": 20, "pactivity": 7.7, "unit": "nM", "confidence": "LOW", "applicability_domain": "OUT OF DOMAIN", "nearest_neighbors": []})
    project, parent, optimization = parent_and_strategy(db)
    optimization.assay_id = None; db.commit()
    row = proposal(db, project, parent, optimization, max_raw=5)
    execute_proposal_run(row.id, session=db); db.refresh(row)
    accepted = [item for item in row.candidates if item.status in {"ACCEPTED", "TOP_10"}]
    assert accepted and all(item.applicability_domain == "OUT_OF_DOMAIN" for item in accepted)
    assert all(item.ranking_score is not None for item in accepted)


def test_low_confidence_herg_does_not_trigger_hard_rejection(db):
    project, parent, optimization = parent_and_strategy(db)
    row = proposal(db, project, parent, optimization)
    candidate = OptimizationCandidate(
        proposal_run_id=row.id, project_id=project.id, optimization_run_id=optimization.id,
        parent_version_id=parent["version"]["id"], candidate_number=1,
        canonical_smiles="CCN(CC)C(=O)c1c(F)cccc1C", isomeric_smiles="CCN(CC)C(=O)c1c(F)cccc1C",
        inchikey="fixture", admet_json={"hERG liability": {"status": "COMPLETE", "predicted_value": 0.99, "confidence": "LOW", "applicability_domain": "IN_DOMAIN", "record_type": "Predicted"}},
    )
    db.add(candidate); db.flush()
    assert proposal_engine._post_prediction_constraints(db, candidate, optimization, {"herg_do_not_increase": True})
    assert candidate.status != "REJECTED"


def test_conflicting_objectives_remain_as_pareto_tradeoff():
    from types import SimpleNamespace
    potency = SimpleNamespace(id=1, objective_vector_json={"values": {"Activity": 0.95, "Solubility": 0.25}}, pareto_front=None)
    soluble = SimpleNamespace(id=2, objective_vector_json={"values": {"Activity": 0.45, "Solubility": 0.9}}, pareto_front=None)
    pareto_fronts([potency, soluble])
    assert potency.pareto_front == soluble.pareto_front == 1


def test_pareto_diversity_information_value_and_transparent_score(db, monkeypatch):
    _, _, _, row = run_fast(db, monkeypatch, max_raw=12, doubles=True)
    accepted = [item for item in row.candidates if item.status in {"ACCEPTED", "TOP_10"}]
    assert accepted and all(item.pareto_front and item.information_value in {"HIGH", "MEDIUM", "LOW"} for item in accepted)
    assert sum(item.selected_top10 for item in accepted) <= 10
    ranking = next(item.rankings[0] for item in accepted if item.rankings)
    assert "100 × max" in ranking.score_breakdown_json["formula"]


def test_manual_promote_and_reject_are_persisted(db, monkeypatch):
    _, _, optimization, row = run_fast(db, monkeypatch, max_raw=8)
    candidates = [item for item in row.candidates if item.status in {"ACCEPTED", "TOP_10"}]
    target = candidates[-1]
    promoted = candidate_decision(target.id, {"decision": "PROMOTED", "reason": "Discriminates hypotheses"}, db)
    assert promoted["user_decision"] == "PROMOTED" and promoted["selected_top10"]
    rejected = candidate_decision(target.id, {"decision": "REJECTED", "reason": "Synthetic route unavailable"}, db)
    assert rejected["status"] == "REJECTED" and rejected["rejection_reasons"][-1]["code"] == "USER_REJECTED"
    db.refresh(row)
    assert row.rejected_count == len([item for item in row.candidates if item.status in {"REJECTED", "FAILED"}])
    assert row.accepted_count == len([item for item in row.candidates if item.status in {"ACCEPTED", "TOP_10"}])


def test_user_added_analog_runs_same_pipeline(db, monkeypatch):
    project, parent, optimization, row = run_fast(db, monkeypatch, max_raw=5)
    candidate = process_user_candidate(db, row, optimization, "CCN(CC)C(=O)c1c(F)cccc1C", "ChemDraw proposal")
    assert candidate.user_added and candidate.stage1_json
    assert candidate.activity_json and candidate.admet_json and candidate.soft_spot_json
    assert candidate.transformations[0].execution_status == "USER_DEFINED"


def test_no_valid_analog_completed_honestly(db, monkeypatch):
    project, parent, optimization = parent_and_strategy(db)
    atom_count = Chem.MolFromSmiles(parent["version"]["canonical_smiles"]).GetNumAtoms()
    optimization.protected_regions_json = [{"id": "ALL", "atom_indices": list(range(atom_count)), "status": "DO NOT MODIFY"}]
    db.commit()
    row = proposal(db, project, parent, optimization, max_raw=10)
    monkeypatch.setattr(proposal_engine, "rescore_admet", fake_admet); monkeypatch.setattr(proposal_engine, "predict_soft_spots", fake_soft)
    execute_proposal_run(row.id, session=db); db.refresh(row)
    assert row.status == "COMPLETED" and row.accepted_count == 0 and row.summary_json["no_valid_analog"]


def test_project_and_parent_version_isolation(db):
    first_project, first_parent, first_optimization = parent_and_strategy(db)
    second_project = create_project(ProjectCreate(name="Other"), db)
    second_parent = create_compound(second_project.id, CompoundCreate(compound_id="OTHER", smiles="CCO"), db)
    row = proposal(db, first_project, first_parent, first_optimization)
    assert row.project_id == first_project.id and row.parent_version_id == first_parent["version"]["id"]
    assert row.parent_version_id != second_parent["version"]["id"]


def test_stage4b_ui_contract_strings():
    source = (ROOT / "frontend/static/app.js").read_text()
    for text in ("Generate analogs", "Show rejected", "Pareto", "Top 10", "Parent vs Candidate", "Promote"):
        assert text in source
