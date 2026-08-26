import math
from pathlib import Path

import pytest
from sqlalchemy import create_engine, select
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend.database import Base
from backend.models import Compound, CompoundVersion, Project, ensure_ui_schema
from backend.pk import (
    PKNCAResult,
    PKObservation,
    PKStudy,
    calculate_bioavailability_for_version,
    ensure_pk_schema,
    normalize_conc_to_ng_ml,
    normalize_dose_to_mg_kg,
    normalize_time_to_hours,
    parse_pk_csv,
    run_nca_calculation,
    serialize_nca,
    serialize_observation,
    serialize_study,
)
from backend.schemas import CompoundCreate, ProjectCreate
from backend.main import create_compound, create_project

ROOT = Path(__file__).parents[1]


@pytest.fixture()
def db():
    engine = create_engine("sqlite://", connect_args={"check_same_thread": False}, poolclass=StaticPool)
    Base.metadata.create_all(engine)
    ensure_ui_schema(engine)
    ensure_pk_schema(engine)
    session = sessionmaker(bind=engine, autoflush=False, autocommit=False)()
    try:
        yield session
    finally:
        session.close()


def setup_compound(db, project_name="PK Test Project", compound_id="PK-001"):
    project = create_project(ProjectCreate(name=project_name, target="Test Target"), db)
    row = create_compound(
        project.id,
        CompoundCreate(
            compound_id=compound_id,
            name=compound_id,
            smiles="CCOc1ccccc1",  # MW = 122.16
            calculate=True,
        ),
        db,
    )
    version = db.get(CompoundVersion, row["version"]["id"])
    return project, row, version


def test_unit_conversion_engine():
    # Time
    assert normalize_time_to_hours(60, "min") == pytest.approx(1.0)
    assert normalize_time_to_hours(3600, "sec") == pytest.approx(1.0)
    assert normalize_time_to_hours(2, "h") == pytest.approx(2.0)
    assert normalize_time_to_hours(1, "day") == pytest.approx(24.0)

    # Concentration mass
    assert normalize_conc_to_ng_ml(1000, "pg/mL") == pytest.approx(1.0)
    assert normalize_conc_to_ng_ml(5.0, "ng/mL") == pytest.approx(5.0)
    assert normalize_conc_to_ng_ml(2.5, "µg/mL") == pytest.approx(2500.0)

    # Concentration molar (MW = 200 g/mol)
    assert normalize_conc_to_ng_ml(1.0, "µM", mw=200.0) == pytest.approx(200.0)
    assert normalize_conc_to_ng_ml(500.0, "nM", mw=200.0) == pytest.approx(100.0)

    # Dose
    assert normalize_dose_to_mg_kg(10.0, "mg/kg") == pytest.approx(10.0)
    assert normalize_dose_to_mg_kg(5000.0, "µg/kg") == pytest.approx(5.0)


def test_nca_analytical_synthetic_iv_bolus():
    """
    Synthetic 1-compartment IV bolus: C(t) = 100 * exp(-0.1 * t).
    Dose = 10 mg/kg.
    Exact theoretical values:
    Cmax = 100 ng/mL, Tmax = 0 h, lambda_z = 0.1 h^-1, t1/2 = 6.9315 h.
    AUCinf = 1000 ng*h/mL.
    CL = 10 mg/kg / 1000 ng*h/mL = 166.667 mL/min/kg.
    Vz = 10 L/h/kg / 0.1 h^-1 = 100 L/kg.
    """
    times = [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
    concs = [100.0 * math.exp(-0.1 * t) for t in times]
    obs = [
        PKObservation(id=i + 1, time_hours=t, concentration_normalized_ng_ml=c)
        for i, (t, c) in enumerate(zip(times, concs))
    ]

    res = run_nca_calculation(obs, route="IV", dose_mg_kg=10.0, dose_unit="mg/kg")

    assert res["status"] == "COMPLETE"
    assert res["cmax"] == pytest.approx(100.0, rel=1e-3)
    assert res["tmax"] == pytest.approx(0.0)
    assert res["lambda_z"] == pytest.approx(0.1, rel=1e-3)
    assert res["terminal_half_life"] == pytest.approx(math.log(2.0) / 0.1, rel=1e-3)
    assert res["aucinf"] == pytest.approx(1000.0, rel=1e-2)
    assert res["cl"] == pytest.approx(166.6667, rel=1e-2)
    assert res["vz"] == pytest.approx(100.0, rel=1e-2)
    assert res["cl_f"] is None  # IV route must not populate CL/F or Vz/F
    assert res["vz_f"] is None


def test_nca_analytical_synthetic_oral():
    """
    Synthetic 1-compartment Oral (PO): C(t) = 150 * (exp(-0.1 * t) - exp(-1.0 * t)).
    Dose = 20 mg/kg.
    Observed peak point in discrete sample is at t = 2.5584 h (or nearest sampled point t=2.0 h).
    AUCinf = 1350 ng*h/mL.
    CL/F = 246.914 mL/min/kg, Vz/F = 148.148 L/kg.
    """
    times = [0.0, 0.5, 1.0, 2.0, 2.5584, 4.0, 8.0, 12.0, 24.0, 36.0]
    concs = [max(0.0, 150.0 * (math.exp(-0.1 * t) - math.exp(-1.0 * t))) for t in times]
    obs = [
        PKObservation(id=i + 1, time_raw=t, time_hours=t, concentration_raw=c, concentration_normalized_ng_ml=c)
        for i, (t, c) in enumerate(zip(times, concs))
    ]

    res = run_nca_calculation(obs, route="PO", dose_mg_kg=20.0, dose_unit="mg/kg")

    assert res["status"] == "COMPLETE"
    assert res["tmax"] == pytest.approx(2.5584, abs=0.1)
    assert res["lambda_z"] == pytest.approx(0.1, rel=5e-2)
    assert res["terminal_half_life"] == pytest.approx(math.log(2.0) / 0.1, rel=5e-2)
    assert res["aucinf"] == pytest.approx(1350.0, rel=5e-2)
    assert res["cl_f"] == pytest.approx(246.914, rel=5e-2)
    assert res["vz_f"] == pytest.approx(148.148, rel=5e-2)
    assert res["cl"] is None  # PO route must not populate absolute CL or Vz
    assert res["vz"] is None


def test_bioavailability_matching_and_isolation(db):
    project, row, version = setup_compound(db)

    # 1. Create IV Study in Rat
    iv_study = PKStudy(
        project_id=project.id,
        compound_row_id=row["row_id"],
        version_id=version.id,
        study_name="Rat IV 5 mg/kg",
        species="Rat",
        route="IV",
        dose=5.0,
        dose_unit="mg/kg",
        dose_normalized_mg_kg=5.0,
    )
    db.add(iv_study)
    db.commit()

    # Add IV observations & run NCA
    times_iv = [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
    concs_iv = [50.0 * math.exp(-0.1 * t) for t in times_iv]
    for t, c in zip(times_iv, concs_iv):
        db.add(PKObservation(pk_study_id=iv_study.id, version_id=version.id, time_raw=t, time_hours=t, concentration_raw=c, concentration_normalized_ng_ml=c))
    db.commit()

    def clean_kwargs(res):
        return {
            "selection_mode": res.get("selection_mode", "AUTO"),
            "cmax": res.get("cmax"),
            "cmax_unit": res.get("cmax_unit", "ng/mL"),
            "tmax": res.get("tmax"),
            "tmax_unit": res.get("tmax_unit", "h"),
            "auclast": res.get("auclast"),
            "auclast_unit": res.get("auclast_unit", "ng*h/mL"),
            "aucinf": res.get("aucinf"),
            "aucinf_unit": res.get("aucinf_unit", "ng*h/mL"),
            "lambda_z": res.get("lambda_z"),
            "terminal_half_life": res.get("terminal_half_life"),
            "mrt": res.get("mrt"),
            "cl": res.get("cl"),
            "cl_unit": res.get("cl_unit", "mL/min/kg"),
            "cl_f": res.get("cl_f"),
            "cl_f_unit": res.get("cl_f_unit", "mL/min/kg"),
            "vz": res.get("vz"),
            "vz_unit": res.get("vz_unit", "L/kg"),
            "vz_f": res.get("vz_f"),
            "vz_f_unit": res.get("vz_f_unit", "L/kg"),
            "aumclast": res.get("aumclast"),
            "aumcinf": res.get("aumcinf"),
            "auc_extrapolated_pct": res.get("auc_extrapolated_pct"),
            "terminal_point_count": res.get("terminal_point_count", 0),
            "terminal_points_json": res.get("terminal_points", []),
            "r_squared": res.get("r_squared"),
            "adjusted_r2": res.get("adjusted_r2"),
            "warnings_json": res.get("warnings", []),
            "blq_policy_json": res.get("blq_policy", {}),
            "nca_engine": res.get("nca_engine"),
            "nca_engine_version": res.get("nca_engine_version"),
            "calculation_method": res.get("calculation_method"),
        }

    iv_obs = db.scalars(select(PKObservation).where(PKObservation.pk_study_id == iv_study.id)).all()
    iv_nca_res = run_nca_calculation(iv_obs, route="IV", dose_mg_kg=5.0, dose_unit="mg/kg")
    db.add(PKNCAResult(pk_study_id=iv_study.id, version_id=version.id, **clean_kwargs(iv_nca_res)))
    db.commit()

    # 2. Check Bioavailability before PO study exists
    ba_before = calculate_bioavailability_for_version(version.id, db)
    assert ba_before["bioavailability"] == []

    # 3. Create PO Study in Rat (same species & version)
    po_study = PKStudy(
        project_id=project.id,
        compound_row_id=row["row_id"],
        version_id=version.id,
        study_name="Rat PO 10 mg/kg",
        species="Rat",
        route="PO",
        dose=10.0,
        dose_unit="mg/kg",
        dose_normalized_mg_kg=10.0,
    )
    db.add(po_study)
    db.commit()

    # Add PO observations (AUCinf = 500 ng*h/mL -> F = (500/10) / (500/5) * 100 = 50%)
    times_po = [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
    concs_po = [55.55 * (math.exp(-0.1 * t) - math.exp(-1.0 * t)) for t in times_po]
    for t, c in zip(times_po, concs_po):
        db.add(PKObservation(pk_study_id=po_study.id, version_id=version.id, time_raw=t, time_hours=t, concentration_raw=c, concentration_normalized_ng_ml=c))
    db.commit()

    po_obs = db.scalars(select(PKObservation).where(PKObservation.pk_study_id == po_study.id)).all()
    po_nca_res = run_nca_calculation(po_obs, route="PO", dose_mg_kg=10.0, dose_unit="mg/kg")
    db.add(PKNCAResult(pk_study_id=po_study.id, version_id=version.id, **clean_kwargs(po_nca_res)))
    db.commit()

    # 4. Check Bioavailability with matched PO + IV
    ba_after = calculate_bioavailability_for_version(version.id, db)
    assert len(ba_after["bioavailability"]) == 1
    match = ba_after["bioavailability"][0]
    assert match["status"] == "MATCHED"
    assert match["label"] == "F_PO"
    assert match["species"] == "Rat"
    assert match["bioavailability_pct"] == pytest.approx(50.0, rel=5e-2)

    # 5. CompoundVersion & Project Isolation
    other_proj, other_row, other_ver = setup_compound(db, "Other Project", "PK-002")
    assert calculate_bioavailability_for_version(other_ver.id, db)["bioavailability"] == []


def test_csv_import_with_custom_mapping_and_blq():
    csv_content = """TIME_HR,CONC_NG_ML,SUBJECT,BLQ
0.0,0.0,Subj_1,0
0.5,45.2,Subj_1,0
1.0,88.4,Subj_1,0
2.0,65.1,Subj_1,0
4.0,22.0,Subj_1,0
8.0,BLQ,Subj_1,1
12.0,BLQ,Subj_1,1
"""
    mapping = {"TIME_HR": "time", "CONC_NG_ML": "concentration", "SUBJECT": "subject", "BLQ": "blq"}
    valid_rows, errors, fieldnames = parse_pk_csv(csv_content, mapping)

    assert not errors
    assert len(valid_rows) == 7
    assert valid_rows[0]["time_raw"] == 0.0
    assert valid_rows[1]["concentration_raw"] == 45.2
    assert valid_rows[5]["blq_flag"] is True


def test_manual_terminal_override_and_warnings():
    times = [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0, 24.0]
    concs = [0.0, 100.0, 80.0, 50.0, 30.0, 20.0, 15.0, 10.0]
    obs = [
        PKObservation(id=i + 1, time_hours=t, concentration_normalized_ng_ml=c)
        for i, (t, c) in enumerate(zip(times, concs))
    ]

    # Auto selection
    auto_res = run_nca_calculation(obs, route="PO", dose_mg_kg=10.0, dose_unit="mg/kg")
    assert auto_res["selection_mode"] == "AUTO"

    # Manual selection override specifying observation IDs 6, 7, 8 (times 12.0, 16.0, 24.0)
    manual_res = run_nca_calculation(
        obs, route="PO", dose_mg_kg=10.0, dose_unit="mg/kg", manual_terminal_indices=[6, 7, 8]
    )
    assert manual_res["selection_mode"] == "MANUAL_OVERRIDE"
    assert len(manual_res["terminal_points"]) == 3
    assert manual_res["terminal_points"] == [6, 7, 8]
