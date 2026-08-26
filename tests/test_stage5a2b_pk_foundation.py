"""Stage 5A-2B — Vd, Absorption Foundation & Route-Aware PK Assembly Tests."""

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

import backend.main as main_module
from backend.database import Base
from backend.ivive import (
    PKParameterSet, assemble_pk_parameter_set, ensure_ivive_schema,
    estimate_absorption_components, estimate_volume_of_distribution,
    get_pk_foundation_profile,
)
from backend.models import Compound, CompoundVersion, Project
from backend.pk import PKNCAResult, PKStudy, ensure_pk_schema


@pytest.fixture()
def db_engine():
    engine = create_engine("sqlite://", connect_args={"check_same_thread": False}, poolclass=StaticPool)
    Base.metadata.create_all(engine)
    ensure_pk_schema(engine)
    ensure_ivive_schema(engine)
    session = sessionmaker(bind=engine, autoflush=False, autocommit=False)()
    try:
        yield session, engine
    finally:
        session.close()


def make_version(db, project_name="Stage 5A-2B Project", label="PKF-001", version_number=1, clogp=2.5):
    full_project_name = f"{project_name} {label}"
    project = Project(name=full_project_name, target="PK Foundation validation")
    db.add(project); db.flush()
    compound = Compound(project_id=project.id, compound_id=label, name=label, current_version=version_number)
    db.add(compound); db.flush()
    version = CompoundVersion(
        compound_row_id=compound.id, version_number=version_number, original_smiles="CCOc1ccccc1",
        canonical_smiles="CCOc1ccccc1", isomeric_smiles="CCOc1ccccc1", inchikey=f"{label}-{version_number}",
        properties_json={"clogp": clogp, "mw": 138.18},
    )
    db.add(version); db.commit()
    return project, compound, version


def add_pk_study(db, version, species="Rat", route="IV", cl=50.0, vz=5.0, cl_f=150.0, vz_f=15.0, mrt=2.0, auc=1000.0):
    study = PKStudy(
        project_id=version.compound.project_id,
        compound_row_id=version.compound_row_id,
        version_id=version.id,
        study_name=f"{species} {route} Study",
        species=species,
        route=route,
        dose=10.0,
        dose_unit="mg/kg",
    )
    db.add(study); db.flush()

    nca = PKNCAResult(
        pk_study_id=study.id,
        version_id=version.id,
        cmax=100.0,
        tmax=0.5 if route != "IV" else 0.0,
        auclast=auc,
        aucinf=auc,
        cl=cl if route == "IV" else None,
        vz=vz if route == "IV" else None,
        cl_f=cl_f if route != "IV" else None,
        vz_f=vz_f if route != "IV" else None,
        mrt=mrt if route == "IV" else None,
        r_squared=0.99,
        auc_extrapolated_pct=5.0,
        warnings_json=[],
    )
    db.add(nca); db.commit()
    return study, nca


def test_vd_architecture_separation(db_engine):
    """Verify Vz != Vz/F, Vss != Vz, and PO V/F is not treated as absolute V."""
    db, _ = db_engine
    project, compound, version = make_version(db)

    # 1. Add PO study only (Vz/F = 15.0 L/kg)
    add_pk_study(db, version, species="Rat", route="PO", cl_f=120.0, vz_f=15.0)
    vd_po = estimate_volume_of_distribution(db, project.id, version.id, "Rat")

    assert vd_po["v_value"] is None or vd_po["v_source_type"] != "EXPERIMENTAL_VZ"
    assert vd_po["apparent_vzf"] is not None
    assert vd_po["apparent_vzf"]["v_value"] == 15.0
    assert vd_po["apparent_vzf"]["v_type"] == "Vz_F"
    assert "Not absolute Vd or Vss" in vd_po["apparent_vzf"]["message"]

    # 2. Add IV study (Vz = 4.0 L/kg, CL = 40 mL/min/kg, MRT = 3.0 h -> Vss = 7.2 L/kg)
    add_pk_study(db, version, species="Rat", route="IV", cl=40.0, vz=4.0, mrt=3.0)
    vd_iv = estimate_volume_of_distribution(db, project.id, version.id, "Rat")

    assert vd_iv["v_type"] == "Vss"
    assert vd_iv["v_source_type"] == "EXPERIMENTAL_VSS"
    assert vd_iv["v_value"] == pytest.approx(40.0 * 60.0 / 1000.0 * 3.0)  # 7.2 L/kg
    assert vd_iv["v_value"] != 4.0  # Vss != Vz
    assert vd_iv["v_value"] != 15.0  # Vss != Vz/F


def test_experimental_v_priority_over_predicted(db_engine):
    """Verify experimental IV V is used over empirical prediction."""
    db, _ = db_engine
    project, compound, version = make_version(db, clogp=3.0)

    # Without IV study -> predicts Lombardo empirical Vd
    res_pred = estimate_volume_of_distribution(db, project.id, version.id, "Rat")
    assert res_pred["v_source_type"] == "MODEL_UNAVAILABLE" or res_pred["v_source_type"] == "PREDICTED_VD"

    # Add IV study
    add_pk_study(db, version, species="Rat", route="IV", cl=30.0, vz=2.5)
    res_exp = estimate_volume_of_distribution(db, project.id, version.id, "Rat")
    assert res_exp["v_source_type"] in {"EXPERIMENTAL_VSS", "EXPERIMENTAL_VZ"}
    assert res_exp["confidence"] == "HIGH"


def test_absorption_components_fa_fg_fh(db_engine):
    """Verify Fa/Fg/Fh separation and predicted F calculation policy."""
    db, _ = db_engine
    project, compound, version = make_version(db)

    abs_info = estimate_absorption_components(db, project.id, version.id, "Rat")
    assert abs_info["fh_value"] is None or isinstance(abs_info["fh_value"], float)
    assert abs_info["fg_status"] == "MODEL_UNAVAILABLE"
    assert abs_info["fg_value"] is None
    # Predicted F unavailable because Fa/Fg not quantitatively supported
    assert abs_info["f_predicted"] is None
    assert "unavailable" in abs_info["f_predicted_message"]


def test_po_cl_f_not_treated_as_iv_cl(db_engine):
    """Verify PO route assembly keeps CL/F and V/F distinct from IV CL and V."""
    db, _ = db_engine
    project, compound, version = make_version(db)

    add_pk_study(db, version, species="Rat", route="PO", cl_f=180.0, vz_f=22.0)

    pset_po = assemble_pk_parameter_set(db, project.id, version.id, "Rat", "PO")
    assert pset_po.cl_value == 180.0
    assert pset_po.cl_source_type == "EXPERIMENTAL_NCA"
    assert pset_po.v_value == 22.0
    assert pset_po.v_type == "Vz_F"
    assert pset_po.v_source_type == "EXPERIMENTAL_VZ_F"

    pset_iv = assemble_pk_parameter_set(db, project.id, version.id, "Rat", "IV")
    assert pset_iv.cl_value is None or pset_iv.cl_source_type != "EXPERIMENTAL_NCA"
    assert pset_iv.v_type != "Vz_F"


def test_species_route_version_isolation(db_engine):
    """Verify complete scope isolation across species, routes, and versions."""
    db, _ = db_engine
    p1, c1, v1 = make_version(db, label="ISO-001")
    p2, c2, v2 = make_version(db, label="ISO-002")

    add_pk_study(db, v1, species="Rat", route="IV", cl=25.0, vz=3.0)
    add_pk_study(db, v2, species="Mouse", route="PO", cl_f=100.0, vz_f=12.0)

    prof1 = get_pk_foundation_profile(db, v1.id, "Rat")
    prof2 = get_pk_foundation_profile(db, v2.id, "Mouse")

    assert prof1["scope"]["compound_id"] == c1.id
    assert prof2["scope"]["compound_id"] == c2.id
    assert prof1["route_parameter_sets"]["IV"]["cl_value"] == 25.0
    assert prof2["route_parameter_sets"]["IV"]["cl_value"] is None
    assert prof2["route_parameter_sets"]["PO"]["cl_value"] == 100.0


def test_confidence_ceiling_and_missing_inputs(db_engine):
    """Verify confidence ceiling rule and missing data handling without dummy values."""
    db, _ = db_engine
    project, compound, version = make_version(db)

    pset_iv = assemble_pk_parameter_set(db, project.id, version.id, "Dog", "IV")
    assert pset_iv.confidence in {"LOW", "MODEL_UNAVAILABLE"}
    assert pset_iv.ka_source_type == "MODEL_UNAVAILABLE"
    assert pset_iv.ka_value is None
    assert pset_iv.fg_status == "NOT_REQUIRED"

    pset_po = assemble_pk_parameter_set(db, project.id, version.id, "Dog", "PO")
    assert pset_po.fg_status == "MODEL_UNAVAILABLE"
    assert pset_po.fg_value is None
