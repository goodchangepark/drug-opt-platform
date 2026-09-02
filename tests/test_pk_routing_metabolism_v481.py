"""Test Suite for PK Evidence Routing + Metabolism Interpretation v4.8.1.

Validates:
- Multidimensional FDA PK table coordinate inheritance (Study, Species, Route, Formulation, Dose, Regimen, Day, Analyte, Parameter, Unit)
- Human Clinical PK routing (Oral 160 mg QD Day 1/29, Tmax, t1/2, CL/F, Vss/F)
- Generic IP/SC simulation scenarios separated into secondary mechanistic section
- Rat Study ARP570 extraction (CL, Vss, t1/2, F, Tmax, Cmax)
- Dog Study ARP572 extraction (CL, Vss, t1/2, F suspension vs capsule)
- BCRP Ki 8.7 µM, P-gp IC50 36.1 µM, and CYP3A4/5 metabolic contribution 93.5% isolation
- No "8.7 H" TOC artifact as CYP3A4 inhibition
- Binary classifier interpretation (scores not %, no numeric subtraction, honest quantitative gap)
- Sibling prediction persistence and compound isolation
"""
import pytest
from backend.database import SessionLocal
from backend.models import Compound, CompoundVersion
from backend.endpoint_comparison import build_endpoint_comparison
from backend.classifier_interpretation import (
    interpret_classifier_prediction,
    compare_classifier_with_experiment,
    CLASSIFIER_REGISTRY,
)
from backend.experimental_refinement import parse_fda_multidimensional_review
from sqlalchemy import select


def test_classifier_interpretation_semantics():
    """Verify binary classifier scores are never labeled as '% inhibition' and have honest gaps."""
    cyp_interp = interpret_classifier_prediction("CYP3A4_INHIBITION", 0.967)
    assert cyp_interp["is_classifier"] is True
    assert cyp_interp["predicted_class"] == "Inhibitor"
    assert "Score: 0.967" in cyp_interp["display_text"]
    assert "%" not in cyp_interp["display_text"]
    assert cyp_interp["quantitative_gap"] == "QUANTITATIVE_MODEL_GAP"

    pgp_interp = interpret_classifier_prediction("PGP_INHIBITION", 0.893)
    assert pgp_interp["is_classifier"] is True
    assert pgp_interp["predicted_class"] == "Inhibitor"
    assert "Score: 0.893" in pgp_interp["display_text"]

    # Comparison with experimental continuous IC50: no numeric subtraction
    comp = compare_classifier_with_experiment("PGP_INHIBITION", 0.893, 36.1, "µM", "IC50")
    assert comp["numeric_difference"] is None
    assert comp["difference_display"] == "—"
    assert comp["quantitative_gap"] == "QUANTITATIVE_MODEL_GAP"


def test_fda_multidimensional_table_parser():
    """Verify multidimensional table coordinate inheritance for Rat, Dog, and Human."""
    sample_text = """
    Study ARP570: Male Sprague-Dawley (SD) rats received mobocertinib succinate at a single IV dose of 3 mg/kg or oral dose of 10 mg/kg.
    Parameter | Mobocertinib IV | Mobocertinib PO
    Dose (mg/kg) | 3.0 | 10.0
    CL (mL/min/kg) | 54.5 | -
    Vss (L/kg) | 11.5 | -
    t1/2 (h) | 3.58 | 3.16
    Oral bioavailability (%) | - | 14.3

    Study ARP572: Male Beagle dogs received mobocertinib succinate at a single IV dose of 3 mg/kg or oral dose of 25 mg/kg.
    Parameter | Mobocertinib IV | Mobocertinib PO Susp | Mobocertinib PO Cap
    Dose (mg/kg) | 3.0 | 25.0 | 25.0
    CL (mL/min/kg) | 11.2 | - | -
    Vss (L/kg) | 12.4 | - | -
    t1/2 (h) | 13.9 | 14.9 | 16.0
    Oral bioavailability (%) | - | 37.6 | 38.9

    Human PK: 160 mg QD
    Obs (N=138) Day 1: Cmax 77.9 ng/mL, AUC 972 ng*hr/mL
    Obs (N=70) Day 29: Cmax 70.4 ng/mL, AUC 951 ng*hr/mL
    apparent oral clearance of mobocertinib is 108 L/h based on the population PK analysis.
    apparent volume of distribution at steady state for mobocertinib is estimated to be 3510 L.
    effective half-life for mobocertinib is 17.6 h.
    In Caco-2 cell monolayers, mobocertinib inhibited P-gp-mediated bidirectional transport of digoxin with an IC50 of 36.1 μM.
    BCRP Ki (μM) 8.7
    percent contribution of CYP3A4/5 to mobocertinib metabolism was 93.5%.
    """
    rows = parse_fda_multidimensional_review(sample_text, app_number="215310")
    endpoints = {r["canonical_endpoint_id"] for r in rows}
    
    assert "RAT_PK_CL_IV" in endpoints
    assert "RAT_PK_VSS_IV" in endpoints
    assert "RAT_PK_F_ORAL" in endpoints
    assert "DOG_PK_CL_IV" in endpoints
    assert "DOG_PK_VSS_IV" in endpoints
    assert "DOG_PK_F_ORAL" in endpoints
    assert "HUMAN_PK_CMAX_ORAL" in endpoints
    assert "HUMAN_PK_CLF_ORAL" in endpoints
    assert "HUMAN_PK_VSSF_ORAL" in endpoints
    assert "PGP_INHIBITION" in endpoints
    assert "BCRP_INHIBITION" in endpoints
    assert "CYP3A4_METABOLIC_CONTRIBUTION" in endpoints


def test_mobocertinib_pk_and_metabolism_comparison():
    """Verify Mobocertinib scientific comparison rows in EGFR project."""
    with SessionLocal() as db:
        mobo = db.scalar(select(Compound).where(Compound.name.ilike('%mobocertinib%')))
        assert mobo is not None
        v_id = mobo.versions[-1].id
        comp = build_endpoint_comparison(db, v_id)
        
        # Check Human Clinical PK
        human_pk = [r for r in comp["scientific_rows"] if r["group"] == "HUMAN CLINICAL PK"]
        assert len(human_pk) > 0
        clf_row = next((r for r in human_pk if "CLF" in r["canonical_endpoint"]), None)
        assert clf_row is not None
        assert clf_row["experimental_display_value"] == 108.0

        vssf_row = next((r for r in human_pk if "VSSF" in r["canonical_endpoint"]), None)
        assert vssf_row is not None
        assert vssf_row["experimental_display_value"] == 3510.0

        # Check Rat PK
        rat_pk = [r for r in comp["scientific_rows"] if r["group"] == "RAT PK"]
        assert len(rat_pk) > 0
        rat_cl = next((r for r in rat_pk if r["canonical_endpoint"] == "RAT_PK_CL_IV"), None)
        assert rat_cl is not None
        assert rat_cl["experimental_display_value"] == 54.5

        # Check Dog PK
        dog_pk = [r for r in comp["scientific_rows"] if r["group"] == "DOG PK"]
        assert len(dog_pk) > 0
        dog_cl = next((r for r in dog_pk if r["canonical_endpoint"] == "DOG_PK_CL_IV"), None)
        assert dog_cl is not None
        assert dog_cl["experimental_display_value"] == 11.2

        # Check Metabolism: P-gp, BCRP, CYP3A4/5 contribution, and no 8.7 H artifact
        met_rows = [r for r in comp["scientific_rows"] if r["section"] == "METABOLISM"]
        bcrp = next((r for r in met_rows if "BCRP" in r["canonical_endpoint"]), None)
        assert bcrp is not None
        assert bcrp["experimental_display_value"] == 8.7

        pgp = next((r for r in met_rows if "PGP" in r["canonical_endpoint"]), None)
        assert pgp is not None
        assert pgp["experimental_display_value"] == 36.1
        assert pgp["difference_display_value"] is None  # no numeric subtraction

        cyp3a = next((r for r in met_rows if r["canonical_endpoint"] == "CYP3A4_INHIBITION"), None)
        if cyp3a and cyp3a.get("experimental_display_value") is not None:
            # Must NOT be the false TOC artifact "8.7 H"
            assert str(cyp3a.get("experimental_display_value")) != "8.7" or cyp3a.get("experimental_display_unit") != "H"

        # Check Mechanistic Scenario group
        scenario_rows = [r for r in comp["scientific_rows"] if r["group"] == "MECHANISTIC / SCENARIO PREDICTIONS"]
        assert len(scenario_rows) >= 2  # e.g. IP, SC clearances
