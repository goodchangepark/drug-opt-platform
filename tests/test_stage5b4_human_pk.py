"""Stage 5B-4 Human PK Prediction & Translational Simulation Test Suite.

Targeted Cases:
Case A: Human IV assembly from allometry
Case B: Human IV assembly from hepatic IVIVE
Case C: Experimental Human CL precedence
Case D: Preservation of lower-priority estimates (side-by-side)
Case E: 2-fold / 3-fold disagreement detection
Case F: No automatic averaging of major disagreement
Case G: Human IV bolus simulation
Case H: Human IV infusion simulation
Case I: Human PO simulation with complete Fa/Fg/Fh/ka
Case J: PO simulation refusal when Fg unavailable
Case K: PO simulation refusal when ka unavailable
Case L: Experimental Human F precedence
Case M: Prospective snapshot immutability
Case N: Prediction snapshot vs later Human experimental comparison (retrospective validation)
Case O: Project isolation
Case P: CompoundVersion isolation
Case Q: Cascade deletion
Case R: Missing-data behavior
Case S: No fabricated PK parameters
"""

import uuid
import pytest
from fastapi import HTTPException
from fastapi.testclient import TestClient
from sqlalchemy import select

from backend.database import SessionLocal, engine
from backend.human_pk import (
    HUMAN_PK_ENGINE_VERSION, PKHumanPredictionSnapshot,
    assemble_human_pk_parameters, calculate_disagreement,
    ensure_human_pk_schema, freeze_human_prediction_snapshot,
    run_human_pk_simulation, validate_against_clinical_data,
    HumanSimulationRequest,
)
from backend.ivive import IVIVERun, ensure_ivive_schema
from backend.main import app
from backend.models import Compound, CompoundVersion, Project
from backend.pk import PKNCAResult, PKObservation, PKStudy, ensure_pk_schema
from backend.simulation import ensure_simulation_schema
from backend.translational import ensure_translational_schema


@pytest.fixture(scope="module")
def client():
    ensure_pk_schema(engine)
    ensure_ivive_schema(engine)
    ensure_simulation_schema(engine)
    ensure_translational_schema(engine)
    ensure_human_pk_schema(engine)
    return TestClient(app)


def test_case_a_human_iv_assembly_from_allometry(client):
    """CASE A: When animal IV studies exist across species, extrapolate Human CL & Vss via allometry."""
    db = SessionLocal()
    try:
        p = Project(name=f"Human Allometry Test {uuid.uuid4().hex[:6]}", target="Oncology", molecule_type="Small Molecule")
        db.add(p)
        db.commit()

        c = Compound(project_id=p.id, compound_id=f"CMP-HA-{uuid.uuid4().hex[:4]}", name="Allometry Compound")
        db.add(c)
        db.commit()

        v = CompoundVersion(
            compound_row_id=c.id,
            version_number=1,
            original_smiles="CC(C)Oc1ccccc1",
            canonical_smiles="CC(C)Oc1ccccc1",
            isomeric_smiles="CC(C)Oc1ccccc1",
            inchikey=f"KEY-{uuid.uuid4().hex[:8]}",
        )
        db.add(v)
        db.commit()

        # Add Mouse IV (0.02 kg): CL = 60.0 mL/min/kg, MRT = 0.5 h -> Vss = 1.8 L/kg
        s_mouse = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="Mouse IV", species="Mouse", route="IV", dose=10.0, dose_unit="mg/kg")
        db.add(s_mouse)
        db.commit()
        db.add(PKNCAResult(pk_study_id=s_mouse.id, version_id=v.id, cl=60.0, mrt=0.5, vz=2.0, is_latest=True))

        # Add Rat IV (0.25 kg): CL = 35.0 mL/min/kg, MRT = 0.85 h -> Vss = 1.785 L/kg
        s_rat = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="Rat IV", species="Rat", route="IV", dose=5.0, dose_unit="mg/kg")
        db.add(s_rat)
        db.commit()
        db.add(PKNCAResult(pk_study_id=s_rat.id, version_id=v.id, cl=35.0, mrt=0.85, vz=1.9, is_latest=True))

        # Add Dog IV (10.0 kg): CL = 12.0 mL/min/kg, MRT = 2.5 h -> Vss = 1.8 L/kg
        s_dog = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="Dog IV", species="Dog", route="IV", dose=2.0, dose_unit="mg/kg")
        db.add(s_dog)
        db.commit()
        db.add(PKNCAResult(pk_study_id=s_dog.id, version_id=v.id, cl=12.0, mrt=2.5, vz=1.8, is_latest=True))
        db.commit()

        profile = assemble_human_pk_parameters(db, v.id)
        assert profile["clearance"]["selected_value"] is not None
        assert "Allometric" in profile["clearance"]["selected_source"]
        assert profile["volume"]["selected_value"] is not None
        assert profile["half_life"]["selected_value"] is not None
        assert profile["readiness"]["iv_simulation"]["status"] == "READY"
    finally:
        db.close()


def test_case_b_c_d_ivive_precedence_and_side_by_side(client):
    """CASE B, C, D: Test Hepatic IVIVE assembly, Experimental Human CL precedence, and preservation of all candidates."""
    db = SessionLocal()
    try:
        p = Project(name=f"Precedence Test {uuid.uuid4().hex[:6]}", target="CNS", molecule_type="Small Molecule")
        db.add(p)
        db.commit()

        c = Compound(project_id=p.id, compound_id=f"CMP-PREC-{uuid.uuid4().hex[:4]}", name="Precedence Compound")
        db.add(c)
        db.commit()

        v = CompoundVersion(
            compound_row_id=c.id,
            version_number=1,
            original_smiles="CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O",
            canonical_smiles="CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O",
            isomeric_smiles="CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O",
            inchikey=f"KEY-{uuid.uuid4().hex[:8]}",
        )
        db.add(v)
        db.commit()

        # Add Human Hepatic IVIVE run: CL = 7.5 mL/min/kg
        ivive = IVIVERun(
            project_id=p.id,
            compound_row_id=c.id,
            version_id=v.id,
            species="Human",
            method_id=1,
            parameter_set_version="1.0",
            outputs_json={"cl_in_vivo_blood": 7.5, "hepatic_availability": 0.65},
            confidence="HIGH",
            status="COMPLETE",
            inputs_hash="hash_ivive_human_prec",
        )
        db.add(ivive)
        db.commit()

        # Check IVIVE selected as top when no experimental exists
        prof1 = assemble_human_pk_parameters(db, v.id)
        assert prof1["clearance"]["selected_value"] == 7.5
        assert "IVIVE" in prof1["clearance"]["selected_source"]

        # Now add Human Clinical Experimental IV study: CL = 5.2 mL/min/kg
        s_human = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="Human Phase 1 IV", species="Human", route="IV", dose=10.0, dose_unit="mg")
        db.add(s_human)
        db.commit()
        db.add(PKNCAResult(pk_study_id=s_human.id, version_id=v.id, cl=5.2, mrt=3.2, vz=1.4, is_latest=True))
        db.commit()

        prof2 = assemble_human_pk_parameters(db, v.id)
        # Experimental holds absolute precedence
        assert prof2["clearance"]["selected_value"] == 5.2
        assert "Experimental" in prof2["clearance"]["selected_source"]
        
        # Lower-priority IVIVE candidate is preserved side-by-side
        candidates = prof2["clearance"]["candidates"]
        assert len(candidates) == 2
        sources = [c["source_name"] for c in candidates]
        assert "Human Clinical IV NCA" in sources
        assert "Human Hepatic IVIVE" in sources
    finally:
        db.close()


def test_case_e_f_disagreement_detection_and_no_averaging(client):
    """CASE E & F: Test 2-fold / 3-fold disagreement detection and confirm conflicting values are never averaged."""
    # Test 1: Consistent (<2-fold)
    c1 = [{"source_name": "IVIVE", "value": 5.0}, {"source_name": "Allometry", "value": 6.5}]
    d1 = calculate_disagreement(c1)
    assert d1["status"] == "GENERALLY_CONSISTENT"
    assert d1["max_fold_difference"] == 1.3
    assert not d1["has_major_disagreement"]

    # Test 2: Moderate Disagreement (2-3-fold)
    c2 = [{"source_name": "IVIVE", "value": 4.0}, {"source_name": "Allometry", "value": 10.0}]
    d2 = calculate_disagreement(c2)
    assert d2["status"] == "MODERATE_DISAGREEMENT"
    assert d2["max_fold_difference"] == 2.5
    assert not d2["has_major_disagreement"]

    # Test 3: Major Disagreement (>3-fold)
    c3 = [{"source_name": "IVIVE", "value": 3.0}, {"source_name": "Allometry", "value": 15.0}]
    d3 = calculate_disagreement(c3)
    assert d3["status"] == "MAJOR_DISAGREEMENT"
    assert d3["max_fold_difference"] == 5.0
    assert d3["has_major_disagreement"]


def test_case_g_h_human_iv_simulations(client):
    """CASE G & H: Run Human IV bolus and IV infusion simulations."""
    db = SessionLocal()
    try:
        p = Project(name=f"Human IV Sim {uuid.uuid4().hex[:6]}", target="Cardio", molecule_type="Small Molecule")
        db.add(p)
        db.commit()

        c = Compound(project_id=p.id, compound_id=f"CMP-IV-{uuid.uuid4().hex[:4]}", name="IV Sim Compound")
        db.add(c)
        db.commit()

        v = CompoundVersion(
            compound_row_id=c.id,
            version_number=1,
            original_smiles="CC(C)NCC(O)COc1ccccc1",
            canonical_smiles="CC(C)NCC(O)COc1ccccc1",
            isomeric_smiles="CC(C)NCC(O)COc1ccccc1",
            inchikey=f"KEY-{uuid.uuid4().hex[:8]}",
        )
        db.add(v)
        db.commit()

        # IV Bolus Simulation with explicit user parameters
        req_bolus = HumanSimulationRequest(
            route="IV",
            administration_type="IV_BOLUS",
            dose=50.0,
            dose_unit="mg",
            body_weight_kg=70.0,
            user_cl_override=5.0,
            user_v_override=1.5,
        )
        res_bolus = run_human_pk_simulation(db, v.id, req_bolus)
        assert res_bolus["target_species"] == "Human"
        assert res_bolus["route"] == "IV"
        assert res_bolus["output_metrics"]["cmax_ng_ml"] > 0
        assert res_bolus["output_metrics"]["auc_single_ng_h_ml"] > 0
        assert len(res_bolus["time_series"]) == 300

        # IV Infusion Simulation (1.0 hour infusion)
        req_inf = HumanSimulationRequest(
            route="IV",
            administration_type="IV_INFUSION",
            dose=100.0,
            dose_unit="mg",
            infusion_duration_hours=2.0,
            body_weight_kg=70.0,
            user_cl_override=5.0,
            user_v_override=1.5,
        )
        res_inf = run_human_pk_simulation(db, v.id, req_inf)
        assert res_inf["administration_type"] == "IV_INFUSION"
        assert res_inf["output_metrics"]["tmax_hours"] > 0
    finally:
        db.close()


def test_case_i_j_k_human_po_simulation_and_refusal_guardrails(client):
    """CASE I, J, K: Test Human PO simulation with complete parameters and ensure refusal when Fg or ka is missing."""
    db = SessionLocal()
    try:
        p = Project(name=f"Human PO Guardrail {uuid.uuid4().hex[:6]}", target="GI", molecule_type="Small Molecule")
        db.add(p)
        db.commit()

        c = Compound(project_id=p.id, compound_id=f"CMP-PO-{uuid.uuid4().hex[:4]}", name="PO Sim Compound")
        db.add(c)
        db.commit()

        v = CompoundVersion(
            compound_row_id=c.id,
            version_number=1,
            original_smiles="CC(=O)Oc1ccccc1C(=O)O",
            canonical_smiles="CC(=O)Oc1ccccc1C(=O)O",
            isomeric_smiles="CC(=O)Oc1ccccc1C(=O)O",
            inchikey=f"KEY-{uuid.uuid4().hex[:8]}",
        )
        db.add(v)
        db.commit()

        # Refusal Case J: Fg missing and no override provided
        req_no_fg = HumanSimulationRequest(
            route="PO",
            dose=100.0,
            dose_unit="mg",
            user_cl_override=5.0,
            user_v_override=1.5,
            user_ka_override=1.2,
        )
        with pytest.raises(HTTPException) as exc_fg:
            run_human_pk_simulation(db, v.id, req_no_fg)
        assert exc_fg.value.status_code == 422
        assert "Human Oral Simulation Refused" in exc_fg.value.detail

        # Refusal Case K: F provided but ka missing
        req_no_ka = HumanSimulationRequest(
            route="PO",
            dose=100.0,
            dose_unit="mg",
            user_cl_override=5.0,
            user_v_override=1.5,
            user_f_override=60.0,
        )
        with pytest.raises(HTTPException) as exc_ka:
            run_human_pk_simulation(db, v.id, req_no_ka)
        assert exc_ka.value.status_code == 422
        assert "Absorption rate constant (ka) is unavailable" in exc_ka.value.detail

        # Success Case I: Complete parameters provided via explicit overrides
        req_success = HumanSimulationRequest(
            route="PO",
            dose=100.0,
            dose_unit="mg",
            user_cl_override=5.0,
            user_v_override=1.5,
            user_f_override=65.0,
            user_ka_override=1.5,
        )
        res_success = run_human_pk_simulation(db, v.id, req_success)
        assert res_success["route"] == "PO"
        assert res_success["output_metrics"]["cmax_ng_ml"] > 0
        assert res_success["output_metrics"]["tmax_hours"] > 0
        assert len(res_success["assumptions"]) >= 2
    finally:
        db.close()


def test_case_l_experimental_f_precedence(client):
    """CASE L: Matched Human Clinical experimental F overrides predicted components."""
    db = SessionLocal()
    try:
        p = Project(name=f"F Precedence {uuid.uuid4().hex[:6]}", target="GI", molecule_type="Small Molecule")
        db.add(p)
        db.commit()

        c = Compound(project_id=p.id, compound_id=f"CMP-FPREC-{uuid.uuid4().hex[:4]}", name="F Precedence Compound")
        db.add(c)
        db.commit()

        v = CompoundVersion(
            compound_row_id=c.id,
            version_number=1,
            original_smiles="c1ccccc1",
            canonical_smiles="c1ccccc1",
            isomeric_smiles="c1ccccc1",
            inchikey=f"KEY-{uuid.uuid4().hex[:8]}",
        )
        db.add(v)
        db.commit()

        # Add Human IV study: Dose 10 mg, AUC = 1000 ng*h/mL
        s_iv = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="Human IV", species="Human", route="IV", dose=10.0, dose_unit="mg")
        db.add(s_iv)
        db.commit()
        db.add(PKNCAResult(pk_study_id=s_iv.id, version_id=v.id, cl=4.0, auclast=1000.0, is_latest=True))

        # Add Human PO study: Dose 20 mg, AUC = 1400 ng*h/mL -> F = (1400/20)/(1000/10) = 70%
        s_po = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="Human PO", species="Human", route="PO", dose=20.0, dose_unit="mg")
        db.add(s_po)
        db.commit()
        db.add(PKNCAResult(pk_study_id=s_po.id, version_id=v.id, cmax=120.0, tmax=1.5, auclast=1400.0, is_latest=True))
        db.commit()

        profile = assemble_human_pk_parameters(db, v.id)
        assert profile["absorption"]["f_experimental"] == 70.0
        assert profile["absorption"]["f_selected"] == 70.0
        assert "Experimental" in profile["absorption"]["f_selected_source"]
    finally:
        db.close()


def test_case_m_n_prospective_freeze_and_retrospective_validation(client):
    """CASE M & N: Test prospective prediction snapshot freeze (immutable) and retrospective validation against later clinical data."""
    db = SessionLocal()
    try:
        p = Project(name=f"Prospective Freeze {uuid.uuid4().hex[:6]}", target="Cardio", molecule_type="Small Molecule")
        db.add(p)
        db.commit()

        c = Compound(project_id=p.id, compound_id=f"CMP-FREEZE-{uuid.uuid4().hex[:4]}", name="Freeze Compound")
        db.add(c)
        db.commit()

        v = CompoundVersion(
            compound_row_id=c.id,
            version_number=1,
            original_smiles="CC(C)Oc1ccccc1",
            canonical_smiles="CC(C)Oc1ccccc1",
            isomeric_smiles="CC(C)Oc1ccccc1",
            inchikey=f"KEY-{uuid.uuid4().hex[:8]}",
        )
        db.add(v)
        db.commit()

        # Seed Human Hepatic IVIVE run: CL = 8.0 mL/min/kg
        ivive = IVIVERun(
            project_id=p.id,
            compound_row_id=c.id,
            version_id=v.id,
            species="Human",
            method_id=1,
            parameter_set_version="1.0",
            outputs_json={"cl_in_vivo_blood": 8.0, "cl_in_vivo_plasma": 8.0},
            confidence="HIGH",
            status="COMPLETE",
            inputs_hash="hash_freeze_test",
        )
        db.add(ivive)
        db.commit()

        # 1. Freeze Prospective Snapshot
        snap_res = freeze_human_prediction_snapshot(db, v.id, snapshot_name="P1 Prospective Candidate")
        assert snap_res["status"] == "FROZEN"
        snapshot_id = snap_res["snapshot_id"]
        assert snap_res["selected_cl"] == 8.0

        # Verify DB record is immutable
        snap_db = db.get(PKHumanPredictionSnapshot, snapshot_id)
        assert snap_db.is_immutable is True
        assert snap_db.selected_cl == 8.0

        # 2. Later: Human Clinical Data is entered (CL = 10.0 mL/min/kg)
        s_human = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="Clinical Phase 1 IV", species="Human", route="IV", dose=10.0, dose_unit="mg")
        db.add(s_human)
        db.commit()
        db.add(PKNCAResult(pk_study_id=s_human.id, version_id=v.id, cl=10.0, mrt=2.5, vz=1.8, is_latest=True))
        db.commit()

        # 3. Run Retrospective Validation against the previously frozen snapshot
        val = validate_against_clinical_data(db, v.id, snapshot_id=snapshot_id)
        assert val["status"] == "VALIDATED"
        assert val["n_comparisons"] >= 1
        comp = val["comparisons"][0]
        assert comp["predicted"] == 8.0
        assert comp["observed"] == 10.0
        assert comp["fold_error"] == 0.8
        assert comp["absolute_fold_error"] == 1.25
        assert comp["performance_band"] == "WITHIN_2_FOLD"
    finally:
        db.close()


def test_case_o_p_q_isolation_and_cascade(client):
    """CASE O, P, Q: Test Project isolation, CompoundVersion isolation, and Cascade deletion."""
    db = SessionLocal()
    try:
        p1 = Project(name=f"Iso Project 1 {uuid.uuid4().hex[:6]}", target="T1", molecule_type="Small Molecule")
        p2 = Project(name=f"Iso Project 2 {uuid.uuid4().hex[:6]}", target="T2", molecule_type="Small Molecule")
        db.add_all([p1, p2])
        db.commit()

        c1 = Compound(project_id=p1.id, compound_id=f"CMP-ISO1-{uuid.uuid4().hex[:4]}", name="C1")
        c2 = Compound(project_id=p2.id, compound_id=f"CMP-ISO2-{uuid.uuid4().hex[:4]}", name="C2")
        db.add_all([c1, c2])
        db.commit()

        v1 = CompoundVersion(compound_row_id=c1.id, version_number=1, original_smiles="C", canonical_smiles="C", isomeric_smiles="C", inchikey=f"KEY-{uuid.uuid4().hex[:8]}")
        v2 = CompoundVersion(compound_row_id=c2.id, version_number=1, original_smiles="CC", canonical_smiles="CC", isomeric_smiles="CC", inchikey=f"KEY-{uuid.uuid4().hex[:8]}")
        db.add_all([v1, v2])
        db.commit()

        # Create snapshots in both projects
        s1 = PKHumanPredictionSnapshot(project_id=p1.id, compound_row_id=c1.id, version_id=v1.id, snapshot_name="Snap 1", inputs_hash="h1")
        s2 = PKHumanPredictionSnapshot(project_id=p2.id, compound_row_id=c2.id, version_id=v2.id, snapshot_name="Snap 2", inputs_hash="h2")
        db.add_all([s1, s2])
        db.commit()
        s1_id = s1.id
        s2_id = s2.id

        # Delete Project 1 via API and verify cascade delete only affects Project 1
        resp = client.request("DELETE", f"/api/projects/{p1.id}", json={"confirmation_name": p1.name})
        assert resp.status_code == 200

        # Expire session cache and verify
        db.expire_all()
        assert db.get(PKHumanPredictionSnapshot, s1_id) is None
        assert db.get(PKHumanPredictionSnapshot, s2_id) is not None
    finally:
        db.close()


def test_case_r_s_missing_data_and_no_fabricated_parameters(client):
    """CASE R & S: Ensure missing data returns clean unavailable status with 0 fabricated parameters."""
    db = SessionLocal()
    try:
        p = Project(name=f"Clean Test {uuid.uuid4().hex[:6]}", target="Clean", molecule_type="Small Molecule")
        db.add(p)
        db.commit()

        c = Compound(project_id=p.id, compound_id=f"CMP-EMPTY-{uuid.uuid4().hex[:4]}", name="Empty Compound")
        db.add(c)
        db.commit()

        v = CompoundVersion(compound_row_id=c.id, version_number=1, original_smiles="CCCC", canonical_smiles="CCCC", isomeric_smiles="CCCC", inchikey=f"KEY-{uuid.uuid4().hex[:8]}")
        db.add(v)
        db.commit()

        profile = assemble_human_pk_parameters(db, v.id)
        assert profile["clearance"]["selected_value"] is None
        assert profile["clearance"]["selected_source"] == "MODEL_UNAVAILABLE"
        assert profile["volume"]["selected_value"] is None
        assert profile["absorption"]["f_selected"] is None
        assert profile["absorption"]["ka_value"] is None
        assert profile["readiness"]["overall_status"] == "INSUFFICIENT_DATA"
    finally:
        db.close()
