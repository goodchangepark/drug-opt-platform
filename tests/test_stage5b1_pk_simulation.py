"""Targeted unit & integration tests for Stage 5B-1 IV PK Simulation Engine."""

import math
import pytest
from fastapi.testclient import TestClient

from backend.database import Base, SessionLocal, engine
from backend.main import app
from backend.models import Compound, CompoundVersion, Project
from backend.pk import PKNCAResult, PKObservation, PKStudy
from backend.ivive import PKParameterSet
from backend.simulation import (
    PKSimulationRun,
    canonicalize_units,
    compute_goodness_of_fit,
    ensure_simulation_schema,
    fit_two_compartment_experimental,
    simulate_one_compartment_iv_bolus,
    simulate_one_compartment_iv_infusion,
    simulate_two_compartment_iv_bolus,
)


@pytest.fixture(autouse=True)
def setup_db():
    Base.metadata.create_all(bind=engine)
    ensure_simulation_schema(engine)
    yield


@pytest.fixture
def client():
    return TestClient(app)


def test_canonical_unit_conversions():
    # 5 mg/kg, 16.67 mL/min/kg (= 1.0 L/h/kg), 2 L/kg
    c = canonicalize_units(5.0, "mg/kg", 16.66666667, "mL/min/kg", 2.0, "L/kg")
    assert pytest.approx(c["dose_mg_kg"], rel=1e-3) == 5.0
    assert pytest.approx(c["cl_l_h_kg"], rel=1e-3) == 1.0
    assert pytest.approx(c["v_l_kg"], rel=1e-3) == 2.0

    # Test µg/kg conversion
    c2 = canonicalize_units(5000.0, "µg/kg", 1.0, "L/h/kg", 2000.0, "mL/kg")
    assert pytest.approx(c2["dose_mg_kg"], rel=1e-3) == 5.0
    assert pytest.approx(c2["v_l_kg"], rel=1e-3) == 2.0

    # Dimensional checks on invalid inputs
    with pytest.raises(ValueError, match="strictly positive"):
        canonicalize_units(-1.0, "mg/kg", 1.0, "L/h/kg", 2.0, "L/kg")
    with pytest.raises(ValueError, match="strictly positive"):
        canonicalize_units(5.0, "mg/kg", 0.0, "L/h/kg", 2.0, "L/kg")
    with pytest.raises(ValueError, match="Unsupported dose unit"):
        canonicalize_units(5.0, "invalid_unit", 1.0, "L/h/kg", 2.0, "L/kg")


def test_one_compartment_iv_bolus_synthetic_case():
    """Synthetic test case:

    Dose = 10 mg/kg
    CL = 1.0 L/h/kg
    V = 2.0 L/kg
    k = CL / V = 0.5 h^-1
    t1/2 = ln(2)/0.5 = 1.3863 h
    C0 = 10 / 2 = 5 mg/L = 5000 ng/mL
    AUCinf = Dose / CL = 10 mg*h/L = 10000 ng*h/mL
    """
    res = simulate_one_compartment_iv_bolus(
        dose_mg_kg=10.0,
        cl_l_h_kg=1.0,
        v_l_kg=2.0,
        t_end_h=24.0,
        num_points=500,
    )
    assert res["k_elim"] == 0.5
    assert pytest.approx(res["half_life_hours"], rel=1e-3) == 1.3863
    assert res["c0_ng_ml"] == 5000.0
    assert res["auc_inf_analytical_ng_h_ml"] == 10000.0
    # Numerical AUC agreement check (< 1.0%)
    assert res["auc_agreement_pct"] >= 99.0


def test_one_compartment_iv_infusion_case():
    """Infusion test case:

    Dose = 10 mg/kg
    Tinf = 2.0 h -> R0 = 5.0 mg/kg/h
    CL = 1.0 L/h/kg
    V = 2.0 L/kg -> k = 0.5 h^-1
    Css_inf = 5 / 1 * 1000 = 5000 ng/mL
    C(Tinf) = 5000 * (1 - exp(-0.5 * 2)) = 5000 * (1 - e^-1) = 3160.6 ng/mL
    Analytical AUCinf = 10000 ng*h/mL
    """
    res = simulate_one_compartment_iv_infusion(
        dose_mg_kg=10.0,
        infusion_duration_h=2.0,
        cl_l_h_kg=1.0,
        v_l_kg=2.0,
        t_end_h=24.0,
        num_points=500,
    )
    assert res["r0_mg_kg_h"] == 5.0
    assert pytest.approx(res["c_tinf_ng_ml"], rel=1e-3) == 3160.6028
    assert res["tmax_hours"] == 2.0
    assert pytest.approx(res["cmax_ng_ml"], rel=1e-3) == 3160.6028
    assert res["auc_inf_analytical_ng_h_ml"] == 10000.0
    assert res["auc_agreement_pct"] >= 99.0


def test_two_compartment_iv_bolus_synthetic_case():
    """2-Compartment forward simulation case."""
    res = simulate_two_compartment_iv_bolus(
        dose_mg_kg=10.0,
        cl_l_h_kg=1.0,
        vc_l_kg=1.0,
        q_l_h_kg=2.0,
        vp_l_kg=2.0,
        t_end_h=24.0,
        num_points=500,
    )
    assert res["c0_ng_ml"] == 10000.0
    assert res["auc_inf_analytical_ng_h_ml"] == 10000.0
    assert res["alpha"] > res["beta"]
    assert len(res["time_series"]) == 501


def test_experimental_two_compartment_fitting():
    class DummyObs:
        def __init__(self, t, c, blq=False):
            self.time_hours = t
            self.concentration_normalized_ng_ml = c
            self.blq_flag = blq

    # Generate synthetic 2-comp points: C(t) = 4000*exp(-2*t) + 1000*exp(-0.2*t)
    obs = [
        DummyObs(0.1, 4000 * math.exp(-0.2) + 1000 * math.exp(-0.02)),
        DummyObs(0.5, 4000 * math.exp(-1.0) + 1000 * math.exp(-0.1)),
        DummyObs(1.0, 4000 * math.exp(-2.0) + 1000 * math.exp(-0.2)),
        DummyObs(2.0, 4000 * math.exp(-4.0) + 1000 * math.exp(-0.4)),
        DummyObs(4.0, 4000 * math.exp(-8.0) + 1000 * math.exp(-0.8)),
        DummyObs(8.0, 4000 * math.exp(-16.0) + 1000 * math.exp(-1.6)),
    ]

    fit = fit_two_compartment_experimental(obs, dose_mg_kg=10.0)
    assert fit["status"] == "FIT_SUCCESS"
    assert fit["rmse"] < 100.0
    assert fit["cl_l_h_kg"] > 0
    assert fit["vc_l_kg"] > 0


import uuid

def test_api_simulation_preview_and_run(client):
    db = SessionLocal()
    try:
        p = Project(name=f"Stage 5B-1 Test Project {uuid.uuid4().hex[:6]}", target="IV Simulation", molecule_type="Small Molecule")
        db.add(p)
        db.commit()
        db.refresh(p)

        c = Compound(project_id=p.id, compound_id="SIM-001", name="Sim Compound")
        db.add(c)
        db.commit()
        db.refresh(c)

        v = CompoundVersion(compound_row_id=c.id, version_number=1, original_smiles="CCO", canonical_smiles="CCO", isomeric_smiles="CCO", inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
        db.add(v)
        db.commit()
        db.refresh(v)

        # Seed experimental IV study + NCA
        study = PKStudy(project_id=p.id, compound_row_id=c.id, version_id=v.id, study_name="Rat IV PK", species="Rat", route="IV", dose=5.0, dose_unit="mg/kg")
        db.add(study)
        db.commit()
        db.refresh(study)

        obs = [
            PKObservation(pk_study_id=study.id, version_id=v.id, time_raw=0.083, time_hours=0.083, concentration_raw=100.0, concentration_normalized_ng_ml=100.0),
            PKObservation(pk_study_id=study.id, version_id=v.id, time_raw=1.0, time_hours=1.0, concentration_raw=50.0, concentration_normalized_ng_ml=50.0),
            PKObservation(pk_study_id=study.id, version_id=v.id, time_raw=4.0, time_hours=4.0, concentration_raw=12.5, concentration_normalized_ng_ml=12.5),
        ]
        db.add_all(obs)

        nca = PKNCAResult(pk_study_id=study.id, version_id=v.id, cl=15.0, vz=1.5, is_latest=True)
        db.add(nca)
        db.commit()

        # Seed IV parameter set
        pset = PKParameterSet(project_id=p.id, compound_row_id=c.id, version_id=v.id, species="Rat", route="IV", cl_value=15.0, cl_unit="mL/min/kg", cl_source_type="EXPERIMENTAL_NCA", v_value=1.5, v_unit="L/kg", v_type="Vss", v_source_type="EXPERIMENTAL_NCA", confidence="HIGH")
        db.add(pset)
        db.commit()

        # 1. Preview API
        prev_res = client.get(f"/api/compound-versions/{v.id}/pk-simulation/preview?species=Rat")
        assert prev_res.status_code == 200
        prev_data = prev_res.json()
        assert prev_data["clearance"]["value"] == 15.0
        assert prev_data["volume"]["value"] == 1.5

        # 2. Run IV Bolus Simulation API
        run_payload = {
            "species": "Rat",
            "administration_type": "IV_BOLUS",
            "dose": 5.0,
            "dose_unit": "mg/kg",
            "model_type": "ONE_COMPARTMENT",
        }
        run_res = client.post(f"/api/compound-versions/{v.id}/pk-simulation/run", json=run_payload)
        assert run_res.status_code == 200
        run_data = run_res.json()
        assert run_data["administration_type"] == "IV_BOLUS"
        assert run_data["output_metrics"]["c0_ng_ml"] == 3333.3333
        assert len(run_data["time_series"]) > 0
        assert len(run_data["residuals"]) == 3
        run_id = run_data["id"]

        # 3. Get History & Get Run API
        hist_res = client.get(f"/api/compound-versions/{v.id}/pk-simulation/history?species=Rat")
        assert hist_res.status_code == 200
        assert len(hist_res.json()) >= 1

        get_res = client.get(f"/api/pk-simulation-runs/{run_id}")
        assert get_res.status_code == 200
        assert get_res.json()["id"] == run_id

        # 4. Delete Run API
        del_res = client.delete(f"/api/pk-simulation-runs/{run_id}")
        assert del_res.status_code == 200

        # Clean up
        db.delete(p)
        db.commit()
    finally:
        db.close()
