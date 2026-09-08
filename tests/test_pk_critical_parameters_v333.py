"""
Unit & Scientific Regression Test Suite: PK Critical Parameters & Candidate v3.3.3.
===================================================================================

Verifies:
1. Canonical fu calculation and physical lower bound (fu >= 0.0001).
2. Well-stirred hepatic clearance calculation and extraction ratio bounds.
3. 1-compartment IV bolus and oral disposition simulations.
4. Full PKParameterSet dataclass assembly and Monte Carlo uncertainty intervals.
5. Scientific maturity integrity: pKa & logD7.4 remain Level 1 (★☆☆☆☆), VDss remains Level 3 (★★★☆☆).
6. Total 50-endpoint maturity distribution (11 L4, 1 L3, 15 L2, 23 L1).
7. DrugBank 250 reference catalog and database integrity (100% collision-free, 100% CAS hydrated).
8. Candidate v3.3.3 policy hash and baseline preservation (v3.3.2 active, v3.3.1 preserved).
9. Clinical PK gate evaluation (PASS_WITH_LIMITATIONS).
"""
import json
import math
import sqlite3
from pathlib import Path
import pytest
from starlette.testclient import TestClient

from backend.main import app
from backend.pk_parameter_set import (
    HUMAN_DEFAULT_BW_KG,
    HUMAN_QH_ML_MIN_KG,
    HUMAN_LIVER_WT_G_KG,
    HUMAN_MPPGL_MG_G,
    canonical_fu_from_ppb,
    calculate_well_stirred_clearance,
    simulate_one_compartment_disposition,
    build_pk_parameter_set,
)
from backend.prediction_maturity import (
    ENDPOINT_MATURITY_REGISTRY,
    get_maturity_statistics,
)
from backend.prediction_engine_registry import (
    CURRENT_ENGINE_ID,
    CURRENT_POLICY_HASH,
    CANDIDATE_ENGINE_ID,
    CANDIDATE_POLICY_HASH,
    get_current_production_engine_info,
    get_candidate_prediction_engine_info,
    get_prediction_model_history,
)
from backend.prediction_engine_v3_3_3_policy import (
    get_v3_3_3_policy_hash,
    ENGINE_V3_3_NAME,
)

client = TestClient(app)


def test_canonical_fu_bounds():
    """Directive: fu must have physical lower bound 0.0001 and upper bound 1.0."""
    # Normal case: 90% bound -> fu = 0.10
    fu = canonical_fu_from_ppb(90.0)
    assert math.isclose(fu, 0.10, rel_tol=1e-5)

    # Moderate bound: 50% bound -> fu = 0.50
    assert math.isclose(canonical_fu_from_ppb(50.0), 0.50, rel_tol=1e-5)

    # Unbound: 0% bound -> fu = 1.0
    assert math.isclose(canonical_fu_from_ppb(0.0), 1.0, rel_tol=1e-5)

    # Extreme bound: 99.999% bound -> clamped to 0.0001 (fail-safe lower bound)
    assert canonical_fu_from_ppb(99.999) == 0.0001
    assert canonical_fu_from_ppb(100.0) == 0.0001

    # Over 100% (experimental noise) -> clamped to 0.0001
    assert canonical_fu_from_ppb(105.0) == 0.0001

    # Negative bound (experimental artifact) -> clamped to 1.0
    assert canonical_fu_from_ppb(-5.0) == 1.0


def test_well_stirred_clearance_calculation():
    """Directive: Hepatic clearance must respect well-stirred model and blood flow limit."""
    cl_int_hlm = 35.0  # mL/min/kg
    fu = 0.05          # 95% bound

    res = calculate_well_stirred_clearance(
        cl_int_input=cl_int_hlm,
        fu=fu,
    )

    # CLh must be strictly less than Qh (20.714 mL/min/kg)
    assert res["cl_h_plasma_ml_min_kg"] < HUMAN_QH_ML_MIN_KG
    assert 0.0 <= res["extraction_ratio_eh"] < 1.0
    assert res["cl_plasma_l_hr_70kg"] > 0.0
    # Whole-body clearance in L/hr = cl_h_plasma_ml_min_kg * 70 * 60 / 1000
    expected_l_hr = res["cl_h_plasma_ml_min_kg"] * 70.0 * 60.0 / 1000.0
    assert math.isclose(res["cl_plasma_l_hr_70kg"], expected_l_hr, rel_tol=1e-3)


def test_one_compartment_pk_simulation():
    """Directive: 1-compartment IV and oral simulations must be physically consistent."""
    # IV bolus test
    iv_sim = simulate_one_compartment_disposition(
        dose_mg=100.0,
        route="IV",
        cl_plasma_ml_min_kg=2.38,
        vdss_l_kg=0.714,
        time_points_hr=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
    )
    assert iv_sim["route"] == "IV"
    assert iv_sim["c0_ng_ml"] > 0.0
    assert iv_sim["auc_inf_ng_hr_ml"] > 0.0
    assert iv_sim["half_life_hr"] > 0.0

    # Oral test
    oral_sim = simulate_one_compartment_disposition(
        dose_mg=100.0,
        route="ORAL",
        cl_plasma_ml_min_kg=2.38,
        vdss_l_kg=0.714,
        f_oral=0.8,
        ka_hr_inv=1.0,
        time_points_hr=[0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
    )
    assert oral_sim["route"] == "ORAL"
    # Oral at t=0 must be 0
    assert oral_sim["c0_ng_ml"] == 0.0
    # Cmax must be positive
    assert oral_sim["cmax_ng_ml"] > 0.0
    # Tmax must be positive
    assert oral_sim["tmax_hr"] > 0.0
    assert oral_sim["auc_inf_ng_hr_ml"] > 0.0


def test_pk_parameter_set_assembly():
    """Directive: Full PKParameterSet assembly with Monte Carlo intervals."""
    pk_set = build_pk_parameter_set(
        smiles="CC(=O)Nc1ccc(O)cc1",
        compound_name="TestCompound",
        dose_mg=500.0,
        route="ORAL",
        ppb_percent=20.0,
        cl_int_hlm_ml_min_kg=15.0,
        vdss_l_kg=0.80,
        caco2_log10_cm_s=-4.8,
    )
    assert pk_set.compound_name == "TestCompound"
    assert math.isclose(pk_set.fu_plasma, 0.80, rel_tol=1e-5)
    assert pk_set.cl_plasma_l_hr_70kg > 0.0
    assert pk_set.vdss_l_70kg == 0.80 * 70.0
    assert pk_set.auc_inf_ng_hr_ml > 0.0
    assert pk_set.cmax_ng_ml > 0.0
    assert pk_set.half_life_hr > 0.0
    assert "cl_plasma_ml_min_kg" in pk_set.confidence_interval_90
    assert "auc_inf_ng_hr_ml" in pk_set.confidence_interval_90


def test_endpoint_maturity_scientific_gates():
    """Directive: Strictly prevent artificial star elevations and preserve honest maturity."""
    # pKa must remain Level 1
    pka = ENDPOINT_MATURITY_REGISTRY["PKA"]
    assert pka["maturity_level"] == 1
    assert pka["stars"] == "★☆☆☆☆"
    assert pka["status"] == "DETERMINISTIC_PROPERTY"

    # logD7.4 must remain Level 1
    logd = ENDPOINT_MATURITY_REGISTRY["LOGD_7_4"]
    assert logd["maturity_level"] == 1
    assert logd["stars"] == "★☆☆☆☆"
    assert logd["status"] == "DETERMINISTIC_PROPERTY"

    # VDss must remain Level 3
    vdss = ENDPOINT_MATURITY_REGISTRY["VDSS"]
    assert vdss["maturity_level"] == 3
    assert vdss["stars"] == "★★★☆☆"

    # Verify overall 50-endpoint distribution
    stats = get_maturity_statistics()
    assert stats["total_endpoints"] == 50
    assert stats["level_breakdown"]["level_4_production_validated"] == 11
    assert stats["level_breakdown"]["level_3_validated_multi_model"] == 1
    assert stats["level_breakdown"]["level_2_validated_base"] == 15
    assert stats["level_breakdown"]["level_1_base"] == 23


def test_drugbank_250_catalog_and_database_integrity():
    """Directive: DrugBank catalog has exactly 250 approved reference drugs with zero collision."""
    catalog_path = Path(__file__).resolve().parent.parent / "backend" / "reference_drugs_250.json"
    assert catalog_path.exists(), "reference_drugs_250.json must exist"

    with open(catalog_path, "r", encoding="utf-8") as f:
        drugs = json.load(f)

    assert len(drugs) == 250, f"Expected exactly 250 drugs, got {len(drugs)}"

    # Check collision-free
    names = [d["name"] for d in drugs]
    db_ids = [d["drugbank_id"] for d in drugs]
    cas_nos = [d["cas_number"] for d in drugs]
    smiles_list = [d["smiles"] for d in drugs]

    assert len(set(names)) == 250, "Drug names must be unique"
    assert len(set(db_ids)) == 250, "DrugBank IDs must be unique"
    assert len(set(cas_nos)) == 250, "CAS numbers must be unique"
    assert len(set(smiles_list)) == 250, "SMILES must be unique"

    # Verify SQLite database Project 300
    db_path = Path(__file__).resolve().parent.parent / "drug_opt.db"
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    cur.execute("SELECT COUNT(*) FROM compounds WHERE project_id = 300;")
    c_count = cur.fetchone()[0]
    assert c_count == 1000, f"Expected 1000 reference compounds in Project 300, got {c_count}"

    # Verify CAS hydration
    cur.execute("""
        SELECT COUNT(DISTINCT compound_id)
        FROM compound_identifiers
        WHERE identifier_type = 'CAS' AND compound_id IN (
            SELECT id FROM compounds WHERE project_id = 300
        );
    """)
    cas_hydrated = cur.fetchone()[0]
    assert cas_hydrated == 250, f"Expected 250 DrugBank CAS identifiers hydrated, got {cas_hydrated}"

    # Verify InChIKey uniqueness in compound_versions for Project 300
    cur.execute("""
        SELECT COUNT(DISTINCT cv.inchikey)
        FROM compound_versions cv
        JOIN compounds c ON c.id = cv.compound_row_id
        WHERE c.project_id = 300;
    """)
    unique_inchikeys = cur.fetchone()[0]
    assert unique_inchikeys == 1000, f"Expected 1000 unique reference InChIKeys, got {unique_inchikeys}"

    # Verify zero orphans
    cur.execute("SELECT COUNT(*) FROM compound_versions WHERE compound_row_id NOT IN (SELECT id FROM compounds);")
    cv_orphans = cur.fetchone()[0]
    assert cv_orphans == 0, "No orphaned compound_versions"

    conn.close()


def test_engine_v3_3_3_candidate_policy_hash():
    """Directive: Candidate v3.3.3 policy hash and baseline preservation."""
    candidate_hash = get_v3_3_3_policy_hash()
    assert candidate_hash == "2ba75ad8813cafd84173369dfbda8abd4190789c16f52f90a905750e620e43d2"
    assert CANDIDATE_POLICY_HASH == candidate_hash

    # v3.3.3 is the current release; its endpoint claims remain independently governed.
    assert CURRENT_ENGINE_ID == "drugopt-prediction-engine-v3@3.3.3"
    assert CURRENT_POLICY_HASH == candidate_hash

    # Historical v3.3.1 strictly preserved
    history = get_prediction_model_history()
    v331 = next((h for h in history if h["version"] == "v3.3.1"), None)
    assert v331 is not None
    assert v331["policy_hash"] == "4647810a58bdbdbc700e4f5c26c5a187032e5cebc80bee6b0d64738f640954a9"

    # Candidate v3.3.3 registered
    v333 = next((h for h in history if h["version"] == "v3.3.3"), None)
    assert v333 is not None
    assert v333["production_status"] == "PRODUCTION_DEFAULT"
    assert v333["reference_compound_N"] == 250


def test_clinical_pk_gate_validation_status():
    """Directive: PK_CRITICAL_GATE formal assessment must be recorded."""
    report_path = Path(__file__).resolve().parent.parent / "validation" / "pk_critical_parameter_validation_report.json"
    assert report_path.exists(), "pk_critical_parameter_validation_report.json must exist"

    with open(report_path, "r", encoding="utf-8") as f:
        report = json.load(f)

    gate_verdict = report.get("report_metadata", {}).get("gate_verdict")
    assert gate_verdict == "PASS_WITH_LIMITATIONS"
    assert report.get("gate_checks", {}).get("canonical_fu_bounds_enforced") is True
    assert len(report.get("clinical_validation", {}).get("drug_reports", [])) >= 6
    assert report.get("internal_projects_validation", {}).get("total_internal_compounds") == 15
    assert report.get("internal_projects_validation", {}).get("successful_pk_simulations") == 15
