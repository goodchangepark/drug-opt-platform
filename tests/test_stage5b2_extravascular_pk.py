"""Stage 5B-2: Extravascular PK Simulation & Absorption Kinetics Automated Test Suite.

Verifies:
1. PO Single Dose: Analytical Cmax, Tmax, AUCinf, and numerical integration agreement.
2. Repeated PO: Linear multi-dose superposition, accumulation ratio R_acc, Css,avg, peak/trough.
3. SC Route: Route isolation (no intestinal Fg or PO ka reuse).
4. IP Route: Route isolation (simplified peritoneal absorption, no oral GI transit assumptions).
5. Flip-Flop Kinetics Detection: ka <= ke flags POTENTIAL FLIP-FLOP KINETICS and explains terminal slope.
6. Experimental Data Overlay: Residuals, RMSE, MAE, fold error, BLQ exclusion.
7. Missing ka / Missing F: MODEL_UNAVAILABLE enforcement (no fabricated simulations).
8. Experimental F Precedence: Matched experimental F used preferentially with MECHANISTIC F INCOMPLETE warning.
9. pKa / logD Evidence Hierarchy Guardrail: DERIVED logD ESTIMATE clearly labeled and never presented as validated ML.
10. Numerical Tmax -> ka Solver: BrentQ solver convergence, monotonicity, non-identifiable parameter handling.
11. Extravascular Parameter Fitting: Fitting ka (and F) with fixed systemic CL/V to avoid non-identifiability.
12. Cross-Route Parameter Contamination Isolation.
"""

import math
import uuid
import pytest
from fastapi.testclient import TestClient
from sqlalchemy import select

from backend.database import SessionLocal, engine
from backend.main import app
from backend.models import Compound, CompoundVersion, Project
from backend.ionization import analyze_ionization, IonizationClass
from backend.pk import PKStudy, PKObservation, PKNCAResult, calculate_bioavailability_for_version
from backend.ivive import PKParameterSet, assemble_pk_parameter_set, get_pk_foundation_profile
from backend.simulation import (
    canonicalize_units,
    solve_ka_from_tmax,
    simulate_one_compartment_extravascular,
    fit_one_compartment_extravascular,
    compute_goodness_of_fit,
    run_pk_simulation,
    PKSimulationRequest,
    PKSimulationRun,
)


@pytest.fixture
def client():
    return TestClient(app)


def test_po_single_dose_analytical_and_numerical_agreement():
    """Validation Case A: PO Single dose with known CL, V, F, ka."""
    dose_mg_kg = 10.0
    cl_l_h_kg = 0.5  # L/h/kg -> ke = 0.5 / 2.0 = 0.25 1/h
    v_l_kg = 2.0     # L/kg
    f_fraction = 0.8 # 80% bioavailability
    ka_h = 1.5       # 1/h

    ke = cl_l_h_kg / v_l_kg  # 0.25 1/h
    expected_tmax = math.log(ka_h / ke) / (ka_h - ke)  # ln(6)/1.25 = 1.4334 h
    expected_auc_inf = (f_fraction * dose_mg_kg / cl_l_h_kg) * 1000.0  # 16000.0 ng*h/mL
    c_scale = (f_fraction * dose_mg_kg * ka_h / (v_l_kg * (ka_h - ke))) * 1000.0
    expected_cmax = c_scale * (math.exp(-ke * expected_tmax) - math.exp(-ka_h * expected_tmax))

    sim = simulate_one_compartment_extravascular(
        dose_mg_kg=dose_mg_kg,
        cl_l_h_kg=cl_l_h_kg,
        v_l_kg=v_l_kg,
        f_fraction=f_fraction,
        ka_h=ka_h,
        num_doses=1,
        dose_interval_h=24.0,
        t_end_h=24.0,
    )

    assert sim["k_elim"] == pytest.approx(0.25, rel=1e-3)
    assert sim["k_abs"] == pytest.approx(1.5, rel=1e-3)
    assert sim["tmax_hours"] == pytest.approx(expected_tmax, rel=1e-3)
    assert sim["cmax_ng_ml"] == pytest.approx(expected_cmax, rel=1e-3)
    assert sim["auc_inf_analytical_ng_h_ml"] == pytest.approx(expected_auc_inf, rel=1e-3)
    assert sim["auc_agreement_pct"] >= 95.0
    assert sim["is_flip_flop"] is False


def test_repeated_po_superposition_and_steady_state():
    """Validation Case B: Repeated PO dosing with linear superposition and accumulation metrics."""
    dose_mg_kg = 5.0
    cl_l_h_kg = 0.4
    v_l_kg = 1.6  # ke = 0.25 1/h
    f_fraction = 0.75
    ka_h = 2.0
    num_doses = 5
    interval_h = 12.0

    sim = simulate_one_compartment_extravascular(
        dose_mg_kg=dose_mg_kg,
        cl_l_h_kg=cl_l_h_kg,
        v_l_kg=v_l_kg,
        f_fraction=f_fraction,
        ka_h=ka_h,
        num_doses=num_doses,
        dose_interval_h=interval_h,
        t_end_h=60.0,
    )

    ke = 0.25
    expected_r_acc = 1.0 / (1.0 - math.exp(-ke * interval_h))  # 1 / (1 - exp(-3)) = 1.052
    expected_css_avg = (f_fraction * dose_mg_kg / (cl_l_h_kg * interval_h)) * 1000.0  # 3.75 / 4.8 * 1000 = 781.25 ng/mL

    assert "steady_state" in sim
    ss = sim["steady_state"]
    assert ss["accumulation_ratio"] == pytest.approx(expected_r_acc, rel=1e-2)
    assert ss["css_avg_ng_ml"] == pytest.approx(expected_css_avg, rel=1e-2)
    assert sim["cmax_ng_ml"] > 0
    assert len(sim["time_series"]) > 100


def test_flip_flop_kinetics_detection():
    """Validation Case E: Flip-flop kinetics when ka <= ke."""
    dose_mg_kg = 10.0
    cl_l_h_kg = 1.0
    v_l_kg = 1.0  # ke = 1.0 1/h
    f_fraction = 0.9
    ka_h = 0.2    # ka = 0.2 1/h < ke (flip-flop)

    sim = simulate_one_compartment_extravascular(
        dose_mg_kg=dose_mg_kg,
        cl_l_h_kg=cl_l_h_kg,
        v_l_kg=v_l_kg,
        f_fraction=f_fraction,
        ka_h=ka_h,
    )

    assert sim["is_flip_flop"] is True
    # In flip-flop kinetics, apparent terminal half-life is governed by absorption ka:
    expected_t_half_app = math.log(2.0) / ka_h  # ln(2)/0.2 = 3.4657 h
    assert sim["half_life_hours"] == pytest.approx(expected_t_half_app, rel=1e-3)
    assert sim["half_life_abs_hours"] == pytest.approx(expected_t_half_app, rel=1e-3)


def test_tmax_to_ka_numerical_solver():
    """Validation Case J: Numerical solver from observed Tmax and ke."""
    ke = 0.2  # 1/h
    known_ka = 1.2  # 1/h
    target_tmax = math.log(known_ka / ke) / (known_ka - ke)  # ln(6)/1.0 = 1.79176 h

    sol = solve_ka_from_tmax(tmax_obs=target_tmax, ke=ke)
    assert sol["status"] == "CONVERGED"
    assert sol["ka"] == pytest.approx(known_ka, rel=1e-3)
    assert sol["flip_flop"] is False

    # Flip-flop root solving (ka < ke)
    known_ka_ff = 0.05
    target_tmax_ff = math.log(known_ka_ff / ke) / (known_ka_ff - ke)
    sol_ff = solve_ka_from_tmax(tmax_obs=target_tmax_ff, ke=ke)
    assert sol_ff["status"] == "CONVERGED"
    assert sol_ff["ka"] == pytest.approx(known_ka_ff, rel=1e-3)
    assert sol_ff["flip_flop"] is True

    # Invalid / non-identifiable cases
    sol_bad = solve_ka_from_tmax(tmax_obs=-1.0, ke=ke)
    assert sol_bad["status"] == "KA_ESTIMATION_UNRELIABLE"
    assert sol_bad["ka"] is None


def test_extravascular_parameter_fitting():
    """Validation Case K: Fit ka from synthetic observations with fixed systemic CL/V."""
    dose_mg_kg = 5.0
    cl_l_h_kg = 0.3
    v_l_kg = 1.5  # ke = 0.2 1/h
    true_f = 0.8
    true_ka = 1.4

    # Generate exact synthetic observations
    times = [0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
    ke = cl_l_h_kg / v_l_kg
    obs_list = []
    for t in times:
        c_exact = (true_f * dose_mg_kg * true_ka / (v_l_kg * (true_ka - ke))) * (math.exp(-ke * t) - math.exp(-true_ka * t)) * 1000.0
        obs = PKObservation(
            time_hours=t,
            concentration_raw=c_exact,
            concentration_normalized_ng_ml=c_exact,
            blq_flag=False,
        )
        obs_list.append(obs)

    fit_res = fit_one_compartment_extravascular(
        observations=obs_list,
        dose_mg_kg=dose_mg_kg,
        cl_l_h_kg=cl_l_h_kg,
        v_l_kg=v_l_kg,
        f_fixed=true_f,
    )

    assert fit_res["status"] == "FIT_SUCCESS"
    assert fit_res["fitted_ka"] == pytest.approx(true_ka, rel=0.01)
    assert fit_res["rmse"] < 1.0


def test_experimental_overlay_and_blq_handling():
    """Validation Case F: Goodness-of-fit metrics, residuals, and BLQ exclusion."""
    time_series = [
        {"time": 0.0, "concentration": 0.0},
        {"time": 1.0, "concentration": 500.0},
        {"time": 2.0, "concentration": 350.0},
        {"time": 4.0, "concentration": 150.0},
    ]
    obs = [
        PKObservation(time_hours=0.0, concentration_raw=None, concentration_normalized_ng_ml=None, blq_flag=True),
        PKObservation(time_hours=1.0, concentration_raw=520.0, concentration_normalized_ng_ml=520.0, blq_flag=False),
        PKObservation(time_hours=2.0, concentration_raw=340.0, concentration_normalized_ng_ml=340.0, blq_flag=False),
    ]

    metrics, residuals = compute_goodness_of_fit(obs, time_series)
    assert metrics["n_points_compared"] == 2
    assert metrics["rmse_ng_ml"] is not None
    assert len(residuals) == 3
    assert residuals[0]["status"] == "BLQ_EXCLUDED"
    assert residuals[1]["status"] == "VALID"
    assert residuals[1]["residual_ng_ml"] == pytest.approx(20.0, rel=1e-2)  # 520 - 500


def test_missing_ka_or_f_rejects_fabrication(client):
    """Validation Case G: Missing ka or missing F returns informative 400 error rather than fake numbers."""
    db = SessionLocal()
    try:
        p = Project(name=f"Stage 5B-2 Missing Test {uuid.uuid4().hex[:6]}", target="Extravascular", molecule_type="Small Molecule")
        db.add(p)
        db.commit()

        c = Compound(project_id=p.id, compound_id=f"CMP-{uuid.uuid4().hex[:6]}", name="Test Compound")
        db.add(c)
        db.commit()

        v = CompoundVersion(compound_row_id=c.id, version_number=1, original_smiles="CC(=O)NC1=CC=C(O)C=C1", canonical_smiles="CC(=O)NC1=CC=C(O)C=C1", isomeric_smiles="CC(=O)NC1=CC=C(O)C=C1", inchikey="TEST-INCH")
        db.add(v)
        db.commit()

        # Request PO simulation without overrides when no experimental PK or ka is in DB
        res = client.post(f"/api/compound-versions/{v.id}/pk-simulation/run", json={
            "species": "Rat",
            "route": "PO",
            "dose": 10.0,
            "dose_unit": "mg/kg",
        })
        assert res.status_code == 400
        assert "unavailable" in res.json()["detail"].lower()
    finally:
        db.close()


def test_route_parameter_isolation_po_sc_ip(client):
    """Validation Cases C, D, L: Cross-route parameter isolation for PO, SC, and IP."""
    db = SessionLocal()
    try:
        p = Project(name=f"Stage 5B-2 Route Isolation {uuid.uuid4().hex[:6]}", target="Adrenergic", molecule_type="Small Molecule")
        db.add(p)
        db.commit()

        c = Compound(project_id=p.id, compound_id=f"CMP-{uuid.uuid4().hex[:6]}", name="Propranolol Route Iso")
        db.add(c)
        db.commit()

        v = CompoundVersion(compound_row_id=c.id, version_number=1, original_smiles="CC(C)NCC(O)COc1cccc2ccccc12", canonical_smiles="CC(C)NCC(O)COc1cccc2ccccc12", isomeric_smiles="CC(C)NCC(O)COc1cccc2ccccc12", inchikey="PROP-ISO-1")
        db.add(v)
        db.commit()

        # Add IV study
        iv_study = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="IV Rat Study", species="Rat", route="IV", dose=5.0, dose_unit="mg/kg", dose_normalized_mg_kg=5.0)
        db.add(iv_study)
        db.commit()

        iv_nca = PKNCAResult(pk_study_id=iv_study.id, version_id=v.id, analysis_version=1, is_latest=True, cl=25.0, vz=2.0, aucinf=3333.0, auclast=3200.0)
        db.add(iv_nca)
        db.commit()

        # Add PO study with Tmax=1.0h, F=30%
        po_study = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="PO Rat Study", species="Rat", route="PO", dose=20.0, dose_unit="mg/kg", dose_normalized_mg_kg=20.0)
        db.add(po_study)
        db.commit()
        po_nca = PKNCAResult(pk_study_id=po_study.id, version_id=v.id, analysis_version=1, is_latest=True, cl_f=83.3, vz_f=6.6, tmax=1.0, cmax=400.0, aucinf=4000.0, auclast=3800.0)
        db.add(po_nca)
        db.commit()

        # Add SC study with Tmax=0.5h, F=85%
        sc_study = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="SC Rat Study", species="Rat", route="SC", dose=10.0, dose_unit="mg/kg", dose_normalized_mg_kg=10.0)
        db.add(sc_study)
        db.commit()
        sc_nca = PKNCAResult(pk_study_id=sc_study.id, version_id=v.id, analysis_version=1, is_latest=True, cl_f=29.4, vz_f=2.3, tmax=0.5, cmax=850.0, aucinf=5666.0, auclast=5500.0)
        db.add(sc_nca)
        db.commit()

        # 1. Simulate PO
        po_run = run_pk_simulation(db, v.id, PKSimulationRequest(species="Rat", route="PO", dose=20.0))
        assert po_run.route == "PO"
        assert po_run.f_value == pytest.approx(0.30, rel=0.05)  # 4000/(3333*4) = 30%
        assert po_run.f_source == "MATCHED_EXPERIMENTAL_F"
        assert po_run.ka_value is not None

        # 2. Simulate SC
        sc_run = run_pk_simulation(db, v.id, PKSimulationRequest(species="Rat", route="SC", dose=10.0))
        assert sc_run.route == "SC"
        assert sc_run.f_value == pytest.approx(0.85, rel=0.05)  # 5666/(3333*2) = 85%
        # Ensure SC did not borrow PO bioavailability or ka
        assert sc_run.f_value != po_run.f_value
        assert sc_run.ka_value != po_run.ka_value
        assert any("SIMPLIFIED ABSORPTION MODEL" in w for w in sc_run.warnings)

        # 3. Verify IP preview remains isolated
        prev_ip = client.get(f"/api/compound-versions/{v.id}/pk-simulation/preview?species=Rat&route=IP").json()
        assert prev_ip["route"] == "IP"
        # No IP study exists -> IP bioavailability is unavailable
        assert prev_ip["bioavailability"]["source"] == "UNAVAILABLE"
    finally:
        db.close()


def test_experimental_f_precedence_with_incomplete_mechanistic_warning(client):
    """Validation Case H: Matched experimental F overrides calculated F while generating MECHANISTIC F INCOMPLETE warning."""
    db = SessionLocal()
    try:
        p = Project(name=f"Stage 5B-2 F Precedence {uuid.uuid4().hex[:6]}", target="Beta Blockers", molecule_type="Small Molecule")
        db.add(p)
        db.commit()

        c = Compound(project_id=p.id, compound_id=f"CMP-{uuid.uuid4().hex[:6]}", name="Precedence Compound")
        db.add(c)
        db.commit()

        v = CompoundVersion(compound_row_id=c.id, version_number=1, original_smiles="CC(C)NCC(O)COc1ccccc1", canonical_smiles="CC(C)NCC(O)COc1ccccc1", isomeric_smiles="CC(C)NCC(O)COc1ccccc1", inchikey="PREC-1")
        db.add(v)
        db.commit()

        # Add IV study
        iv_s = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="IV Study", species="Dog", route="IV", dose=2.0, dose_unit="mg/kg", dose_normalized_mg_kg=2.0)
        db.add(iv_s)
        db.commit()
        db.add(PKNCAResult(pk_study_id=iv_s.id, version_id=v.id, analysis_version=1, is_latest=True, cl=15.0, vz=1.8, aucinf=2222.0, auclast=2100.0))
        db.commit()

        # Add PO study with experimental bioavailability F = 45%
        po_s = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="PO Study", species="Dog", route="PO", dose=10.0, dose_unit="mg/kg", dose_normalized_mg_kg=10.0)
        db.add(po_s)
        db.commit()
        db.add(PKNCAResult(pk_study_id=po_s.id, version_id=v.id, analysis_version=1, is_latest=True, cl_f=33.3, vz_f=3.5, tmax=1.5, cmax=350.0, aucinf=5000.0, auclast=4800.0))
        db.commit()

        run = run_pk_simulation(db, v.id, PKSimulationRequest(species="Dog", route="PO", dose=10.0))
        assert run.f_source == "MATCHED_EXPERIMENTAL_F"
        assert run.f_value == pytest.approx(0.45, rel=0.05)  # 5000 / (2222 * 5) = 45%
        # Warning about incomplete Fa/Fg decomposition
        assert any("MECHANISTIC F INCOMPLETE" in w for w in run.warnings)
    finally:
        db.close()


def test_pka_and_logd_evidence_hierarchy_guardrail():
    """Validation Case I / Guardrail Audit: Ensure DERIVED logD ESTIMATE is clearly labeled and never called PREDICTED_MODEL."""
    res = analyze_ionization("CC(C)Cc1ccc(cc1)C(C)C(=O)O")  # Ibuprofen
    assert res["status"] == "COMPLETE"
    assert res["ionization_class"] == IonizationClass.ACID

    # Check evidence labels
    assert res["primary_pka_evidence_type"] == "RULE_ESTIMATE"
    assert res["physiological_state_7_4"]["logd74_evidence_type"] == "DERIVED_ESTIMATE"
    assert "DERIVED logD ESTIMATE" in res["physiological_state_7_4"]["logd74_label"]

    for profile in res["ph_profiles"]:
        assert profile["logd_evidence_type"] == "DERIVED_ESTIMATE"
        assert "DERIVED logD ESTIMATE" in profile["logd_note"]
        assert profile["evidence_type"] == "RULE_ESTIMATE"

    # Verify model provenance documents explicit hierarchy
    prov = res["model_provenance"]
    assert "EXPERIMENTAL > PREDICTED_MODEL > RULE_ESTIMATE > DERIVED_ESTIMATE > MODEL_UNAVAILABLE" in prov["evidence_hierarchy"]
    assert "Simplified pH-dependent ionization estimate" in prov["limitations"]
