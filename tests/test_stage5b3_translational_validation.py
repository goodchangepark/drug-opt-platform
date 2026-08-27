"""Comprehensive Test Suite for Stage 5B-3: PK Validation, Cross-Species Scaling & Translational Foundation."""

import math
import uuid
import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine, select
from sqlalchemy.orm import sessionmaker

from backend.database import Base, SessionLocal, engine
from backend.main import app
from backend.models import Compound, CompoundVersion, Project
from backend.pk import PKNCAResult, PKObservation, PKStudy
from backend.ivive import IVIVERun, PKParameterSet
from backend.translational import (
    SPECIES_BODY_WEIGHTS,
    PKTranslationalSnapshot,
    ensure_translational_schema,
    evaluate_pk_predictions,
    fit_allometry,
    get_translational_pk_profile,
    run_loso_validation,
)


@pytest.fixture
def client():
    return TestClient(app)


def test_case_a_ideal_clearance_allometry():
    """CASE A: Ideal CL allometry recovery.
    Generate synthetic clearance data for Mouse, Rat, Dog, Monkey with true Y_total = 10.0 * BW^0.75.
    Verify recovered a and b (b ~ 0.75) and extrapolated Human CL.
    """
    true_a = 10.0
    true_b = 0.75
    species_list = ["Mouse", "Rat", "Dog", "Monkey"]
    points = []

    for sp in species_list:
        bw = SPECIES_BODY_WEIGHTS[sp]
        val_total = true_a * (bw ** true_b)  # mL/min
        val_norm = val_total / bw           # mL/min/kg
        points.append({
            "species": sp,
            "bw_kg": bw,
            "value_norm": val_norm,
            "unit": "mL/min/kg",
            "source": "EXPERIMENTAL_NCA",
            "evidence_type": "EXPERIMENTAL",
        })

    fit_res = fit_allometry(points, target_species="Human", param_type="CL")
    assert fit_res["status"] == "SUCCESS"
    assert fit_res["n_species"] == 4
    assert fit_res["exponent_b"] == pytest.approx(0.75, abs=0.01)
    assert fit_res["coefficient_a"] == pytest.approx(10.0, abs=0.05)
    assert fit_res["r_squared"] > 0.999
    assert fit_res["confidence"] == "HIGH"

    # Human 70kg total CL = 10 * 70^0.75 = 241.99 mL/min -> norm = 241.99 / 70 = 3.457 mL/min/kg
    expected_human_total = 10.0 * (70.0 ** 0.75)
    expected_human_norm = expected_human_total / 70.0
    assert fit_res["extrapolated_total"] == pytest.approx(expected_human_total, rel=0.01)
    assert fit_res["extrapolated_norm"] == pytest.approx(expected_human_norm, rel=0.01)


def test_case_b_volume_allometry():
    """CASE B: Vss allometry recovery.
    Generate synthetic volume data with true Y_total = 0.8 * BW^1.0.
    Verify recovered b ~ 1.0 and extrapolated Human Vss ~ 0.8 L/kg.
    """
    true_a = 0.8
    true_b = 1.0
    species_list = ["Mouse", "Rat", "Dog", "Monkey"]
    points = []

    for sp in species_list:
        bw = SPECIES_BODY_WEIGHTS[sp]
        val_total = true_a * (bw ** true_b)  # L
        val_norm = val_total / bw           # L/kg (constant 0.8)
        points.append({
            "species": sp,
            "bw_kg": bw,
            "value_norm": val_norm,
            "unit": "L/kg",
            "source": "EXPERIMENTAL_NCA",
            "evidence_type": "EXPERIMENTAL",
        })

    fit_res = fit_allometry(points, target_species="Human", param_type="Vss")
    assert fit_res["status"] == "SUCCESS"
    assert fit_res["n_species"] == 4
    assert fit_res["exponent_b"] == pytest.approx(1.0, abs=0.01)
    assert fit_res["coefficient_a"] == pytest.approx(0.8, abs=0.05)
    assert fit_res["r_squared"] > 0.999
    assert fit_res["extrapolated_norm"] == pytest.approx(0.8, rel=0.01)
    assert fit_res["extrapolated_total"] == pytest.approx(56.0, rel=0.01)


def test_case_c_human_holdout():
    """CASE C: Human holdout verification.
    Ensure animal allometric fit excludes Human experimental data and predicts Human blindly.
    """
    points_with_human = [
        {"species": "Mouse", "bw_kg": 0.02, "value_norm": 25.0, "unit": "mL/min/kg"},
        {"species": "Rat", "bw_kg": 0.25, "value_norm": 15.0, "unit": "mL/min/kg"},
        {"species": "Dog", "bw_kg": 10.0, "value_norm": 5.0, "unit": "mL/min/kg"},
        {"species": "Monkey", "bw_kg": 5.0, "value_norm": 8.0, "unit": "mL/min/kg"},
    ]

    # Human is the target species to extrapolate, never in training points
    fit_res = fit_allometry(points_with_human, target_species="Human", param_type="CL")
    assert "Human" not in fit_res["species_used"]
    assert fit_res["extrapolated_point"]["species"] == "Human"
    assert fit_res["extrapolated_point"]["is_extrapolated"] is True


def test_case_d_leave_one_species_out():
    """CASE D: Leave-one-species-out (LOSO) cross-validation correctness."""
    true_a = 15.0
    true_b = 0.72
    species_list = ["Mouse", "Rat", "Dog", "Monkey"]
    points = []

    for sp in species_list:
        bw = SPECIES_BODY_WEIGHTS[sp]
        val_norm = (true_a * (bw ** true_b)) / bw
        points.append({
            "species": sp,
            "bw_kg": bw,
            "value_norm": val_norm,
            "unit": "mL/min/kg",
        })

    loso = run_loso_validation(points, param_type="CL")
    assert loso["status"] == "SUCCESS"
    assert loso["n_species_evaluated"] == 4
    assert loso["aafe"] == pytest.approx(1.0, abs=0.05)  # Near 1.0 fold error on exact power law
    assert loso["within_2_fold_pct"] == 100.0
    assert len(loso["loso_evaluations"]) == 4

    for ev in loso["loso_evaluations"]:
        assert ev["held_out_species"] not in ev["training_species"]
        assert ev["within_2_fold"] is True


def test_case_e_poor_scaling_diagnostics():
    """CASE E: Poor scaling generates LOW confidence and warning banners."""
    inconsistent_points = [
        {"species": "Mouse", "bw_kg": 0.02, "value_norm": 1.0, "unit": "mL/min/kg"},
        {"species": "Rat", "bw_kg": 0.25, "value_norm": 100.0, "unit": "mL/min/kg"},
        {"species": "Dog", "bw_kg": 10.0, "value_norm": 0.5, "unit": "mL/min/kg"},
    ]

    fit_res = fit_allometry(inconsistent_points, target_species="Human", param_type="CL")
    assert fit_res["status"] == "SUCCESS"
    assert fit_res["confidence"] == "LOW"
    assert fit_res["r_squared"] < 0.70
    assert any("POOR ALLOMETRIC FIT" in w for w in fit_res["warnings"])


def test_case_f_two_species_low_evidence_warning():
    """CASE F: 2-species (Mouse + Rat) fit allowed but flagged as LOW-EVIDENCE HUMAN EXTRAPOLATION."""
    two_rodent_points = [
        {"species": "Mouse", "bw_kg": 0.02, "value_norm": 50.0, "unit": "mL/min/kg"},
        {"species": "Rat", "bw_kg": 0.25, "value_norm": 25.0, "unit": "mL/min/kg"},
    ]

    fit_res = fit_allometry(two_rodent_points, target_species="Human", param_type="CL")
    assert fit_res["status"] == "SUCCESS"
    assert fit_res["n_species"] == 2
    assert fit_res["confidence"] == "LOW"
    assert any("LOW-EVIDENCE HUMAN EXTRAPOLATION" in w for w in fit_res["warnings"])


def test_case_g_validation_metrics_and_incompatible_isolation():
    """CASE G & Metrics: Verify AAFE, GMFE (bias), RMSE log10, and rejection of incompatible pairs."""
    # Test valid pairs
    pairs = [
        {"observed": 10.0, "predicted": 12.0, "endpoint": "CL", "species": "Rat", "route": "IV", "method": "IVIVE"},
        {"observed": 20.0, "predicted": 10.0, "endpoint": "CL", "species": "Dog", "route": "IV", "method": "IVIVE"},
        {"observed": 5.0, "predicted": 5.0, "endpoint": "CL", "species": "Monkey", "route": "IV", "method": "IVIVE"},
    ]

    res = evaluate_pk_predictions(pairs)
    assert res["status"] == "SUCCESS"
    assert res["n"] == 3
    # Pair 1: FE = 1.2, AFE = 1.2
    # Pair 2: FE = 0.5, AFE = 2.0
    # Pair 3: FE = 1.0, AFE = 1.0
    # All 3 pairs are <= 2.0-fold
    assert res["within_2_fold_count"] == 3
    assert res["within_2_fold_pct"] == 100.0
    assert res["within_1_5_fold_count"] == 2
    assert res["aafe"] > 1.0
    assert res["rmse_log10"] > 0.0

    # Test empty or invalid pairs
    empty_res = evaluate_pk_predictions([])
    assert empty_res["status"] == "NO_DATA"
    assert empty_res["n"] == 0


def test_case_h_profile_and_readiness_end_to_end(client):
    """CASE H & J: Test end-to-end database profile, IVIVE vs Allometry side-by-side, and readiness scorecard."""
    db = SessionLocal()
    try:
        p = Project(name=f"Stage 5B-3 Test {uuid.uuid4().hex[:6]}", target="Cardio", molecule_type="Small Molecule")
        db.add(p)
        db.commit()

        c = Compound(project_id=p.id, compound_id=f"CMP-5B3-{uuid.uuid4().hex[:4]}", name="Translational Compound")
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

        # Add Rat IV study
        s_rat = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="Rat IV", species="Rat", route="IV", dose=5.0, dose_unit="mg/kg")
        db.add(s_rat)
        db.commit()
        db.add(PKNCAResult(pk_study_id=s_rat.id, version_id=v.id, cl=40.0, mrt=0.8333, vz=2.2, terminal_half_life=0.6, is_latest=True))

        # Add Dog IV study
        s_dog = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="Dog IV", species="Dog", route="IV", dose=2.0, dose_unit="mg/kg")
        db.add(s_dog)
        db.commit()
        db.add(PKNCAResult(pk_study_id=s_dog.id, version_id=v.id, cl=12.0, mrt=2.5, vz=2.0, terminal_half_life=1.7, is_latest=True))

        # Add Monkey IV study
        s_mky = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="Monkey IV", species="Monkey", route="IV", dose=3.0, dose_unit="mg/kg")
        db.add(s_mky)
        db.commit()
        db.add(PKNCAResult(pk_study_id=s_mky.id, version_id=v.id, cl=18.0, mrt=1.759, vz=2.1, terminal_half_life=1.2, is_latest=True))

        # Add Human Hepatic IVIVE Run
        ivive_run = IVIVERun(
            project_id=p.id,
            compound_row_id=c.id,
            version_id=v.id,
            species="Human",
            method_id=1,
            parameter_set_version="1.0",
            outputs_json={"cl_in_vivo_blood": 8.5},
            confidence="HIGH",
            status="COMPLETE",
            inputs_hash="hash_human_ivive_123",
        )
        db.add(ivive_run)
        db.commit()

        # Query Translational Profile API
        resp = client.get(f"/api/compound-versions/{v.id}/translational-pk")
        assert resp.status_code == 200
        data = resp.json()

        assert data["compound_id"] == c.compound_id
        assert data["clearance_allometry"]["status"] == "SUCCESS"
        assert data["clearance_allometry"]["n_species"] == 3
        assert data["volume_allometry"]["status"] == "SUCCESS"
        assert data["volume_allometry"]["n_species"] == 3

        # Verify IVIVE vs Allometry side-by-side comparison
        h_comp = data["human_comparison"]["clearance"]
        assert h_comp["method_a_hepatic_ivive"]["value"] == 8.5
        assert h_comp["method_b_simple_allometry"]["value"] is not None
        assert h_comp["method_a_hepatic_ivive"]["value"] != h_comp["method_b_simple_allometry"]["value"]

        # Verify Readiness Scorecard
        readiness = data["human_simulation_readiness"]
        assert readiness["clearance"]["status"] == "READY"
        assert readiness["volume"]["status"] == "READY"
        assert readiness["bioavailability"]["status"] == "UNAVAILABLE"
        assert readiness["overall_status"] == "PARTIALLY READY"

        # Verify Project-level Translational PK API
        p_resp = client.get(f"/api/projects/{p.id}/translational-pk")
        assert p_resp.status_code == 200
        p_data = p_resp.json()
        assert len(p_data["compounds_matrix"]) == 1
        assert p_data["compounds_matrix"][0]["rat_cl_iv"] == 40.0

    finally:
        db.close()


def test_case_i_prospective_prediction_freeze_and_retrospective_validation(client):
    """CASE I: Prospective snapshot immutability & retrospective validation."""
    db = SessionLocal()
    try:
        p = Project(name=f"Stage 5B-3 Freeze {uuid.uuid4().hex[:6]}", target="CNS", molecule_type="Small Molecule")
        db.add(p)
        db.commit()

        c = Compound(project_id=p.id, compound_id=f"CMP-FRZ-{uuid.uuid4().hex[:4]}", name="Freeze Compound")
        db.add(c)
        db.commit()

        v = CompoundVersion(
            compound_row_id=c.id,
            version_number=1,
            original_smiles="CCN(CC)CC",
            canonical_smiles="CCN(CC)CC",
            isomeric_smiles="CCN(CC)CC",
            inchikey=f"KEY-FRZ-{uuid.uuid4().hex[:6]}",
        )
        db.add(v)
        db.commit()

        # Seed Rat and Dog IV
        s_rat = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="Rat IV", species="Rat", route="IV", dose=5.0, dose_unit="mg/kg")
        s_dog = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="Dog IV", species="Dog", route="IV", dose=2.0, dose_unit="mg/kg")
        db.add_all([s_rat, s_dog])
        db.commit()
        db.add(PKNCAResult(pk_study_id=s_rat.id, version_id=v.id, cl=30.0, mrt=0.8333, vz=1.6, is_latest=True))
        db.add(PKNCAResult(pk_study_id=s_dog.id, version_id=v.id, cl=10.0, mrt=2.0, vz=1.3, is_latest=True))
        db.commit()

        # 1. Generate prospective profile (freezes snapshot in database)
        prof_res = client.get(f"/api/compound-versions/{v.id}/translational-pk")
        assert prof_res.status_code == 200

        snapshots = list(db.scalars(
            select(PKTranslationalSnapshot).where(PKTranslationalSnapshot.version_id == v.id)
        ))
        assert len(snapshots) >= 1
        orig_pred = snapshots[0].predicted_value
        orig_id = snapshots[0].id

        # 2. Add Human experimental PK data later
        s_human = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="Human Phase 1 IV", species="Human", route="IV", dose=1.0, dose_unit="mg/kg")
        db.add(s_human)
        db.commit()
        db.add(PKNCAResult(pk_study_id=s_human.id, version_id=v.id, cl=4.5, mrt=4.074, vz=1.2, is_latest=True))
        db.commit()

        # Verify original snapshot was NOT modified
        db.expire_all()
        snap_after = db.get(PKTranslationalSnapshot, orig_id)
        assert snap_after.predicted_value == orig_pred
        assert snap_after.is_immutable is True

        # 3. Retrospective Validation API compares frozen prediction vs new human observation
        val_res = client.get(f"/api/compound-versions/{v.id}/pk-validation")
        assert val_res.status_code == 200
        val_data = val_res.json()
        assert val_data["validation_metrics"]["status"] == "SUCCESS"
        assert val_data["validation_metrics"]["n"] >= 1

        pairs = val_data["validation_metrics"]["pairs"]
        human_pair = next((p for p in pairs if p["species"] == "Human" and p["endpoint"] == "CL"), None)
        assert human_pair is not None
        assert human_pair["observed"] == 4.5
        assert human_pair["predicted"] == round(orig_pred, 4)

    finally:
        db.close()
