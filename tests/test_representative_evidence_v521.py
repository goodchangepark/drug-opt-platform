"""Unit and integration tests for Representative Evidence Correction v5.2.1."""
import pytest
from backend.evidence_qualification_v51 import (
    qualify_evidence_record_v51,
    STATE_AUTO_QUALIFIED,
    STATE_RELATED,
)
from backend.representative_experimental import (
    representative_rank,
    select_representative,
    REPRESENTATIVE_EXPERIMENTAL_VERSION,
)
from backend.canonical_endpoints import normalize_experimental_observation
from backend.database import SessionLocal
from backend.models import Compound
from backend.endpoint_comparison import build_endpoint_comparison


def test_orforglipron_ratio_ec50_classified_as_ratio():
    """Ratio EC50 fold shift is classified as Ratio/Selectivity_Ratio, not molar EC50."""
    rec = {
        "raw_endpoint_name": "Ratio EC50",
        "raw_value": "1681.0",
        "raw_unit": "",
        "species": "HUMAN",
        "source_database": "ChEMBL",
        "assay_conditions_json": {
            "conditions": "Positive allosteric modulator activity at human GLP-1R expressed in PSC-HEK293 cells assessed as potentiation of GLP1(9-36)NH2-induced cAMP accumulation by measuring shift in EC50 of endogenous GLP1(9-36)NH2 at 10 uM incubated for 30 mins by HTRF cAMP assay relative to control",
            "target": "CHEMBL1784"
        }
    }
    decision = qualify_evidence_record_v51(rec)
    assert decision.measurement_type == "Selectivity_Ratio"
    assert decision.canonical_endpoint_id == "ACTIVITY_SELECTIVITY_RATIO" or decision.canonical_endpoint_id == "ACTIVITY_RATIO"
    assert decision.normalized_unit == "ratio"
    assert decision.evidence_state == STATE_RELATED
    assert decision.comparability_status == "RELATED_NOT_SAME_ENDPOINT"


def test_orforglipron_direct_agonist_ec50_precedence():
    """Direct human GLP-1R agonist cAMP EC50 (1.2 nM) takes precedence over allosteric modulation (600 nM)."""
    # 1. Direct agonist cAMP assay
    exp_direct = {
        "id": 144,
        "raw_endpoint": "EC50",
        "normalized_value": 1.2,
        "normalized_unit": "nM",
        "origin": "AUTO_QUALIFIED_EXTERNAL",
        "comparability": "DIRECTLY_COMPARABLE",
        "display": {"value": 1.2, "unit": "nM"},
        "context": {
            "species": "HUMAN",
            "conditions": "Agonist activity at human GLP-1 expressed in HEK293 cells assessed as reduction in intracellular cAMP accumulation incubated for 60 mins by HTRF assay",
            "target": "GLP-1R"
        }
    }
    # 2. Allosteric modulation assay in presence of peptide
    exp_pam = {
        "id": 141,
        "raw_endpoint": "EC50",
        "normalized_value": 600.0,
        "normalized_unit": "nM",
        "origin": "AUTO_QUALIFIED_EXTERNAL",
        "comparability": "DIRECTLY_COMPARABLE",
        "display": {"value": 600.0, "unit": "nM"},
        "context": {
            "species": "HUMAN",
            "conditions": "Positive allosteric modulator activity at human GLP-1R expressed in PSC-HEK293 cells in presence of EC20 level of GLP1(9-36)NH2 incubated for 30 mins by HTRF cAMP assay",
            "target": "GLP-1R"
        }
    }

    rank_direct = representative_rank(exp_direct)
    rank_pam = representative_rank(exp_pam)
    assert rank_direct < rank_pam, f"Direct agonist rank ({rank_direct}) must be better than PAM rank ({rank_pam})"

    selected, reason = select_representative([exp_pam, exp_direct])
    assert selected["id"] == 144
    assert selected["normalized_value"] == 1.2


def test_mobocertinib_day1_vs_day29_separation():
    """Mobocertinib 160 mg QD Day 1 and Day 29 are separated into distinct PK rows."""
    rec_day1 = {
        "raw_endpoint_name": "Cmax",
        "raw_value": "77.9",
        "raw_unit": "ng/mL",
        "species": "HUMAN",
        "source_database": "Drugs@FDA",
        "reference_text": "Drugs@FDA NDA215310 · Table 2.b Clinical Pharmacology (Day 1, 160 mg QD)",
        "assay_conditions_json": {
            "species": "HUMAN", "route": "ORAL", "dose": 160.0, "dose_unit": "mg", "regimen": "QD", "day": 1
        }
    }
    dec_day1 = qualify_evidence_record_v51(rec_day1)
    assert dec_day1.canonical_endpoint_id == "HUMAN_PK_CMAX_ORAL_DAY1"
    assert "Day 1" in dec_day1.display_name
    assert dec_day1.normalized_value == 77.9

    rec_day29 = {
        "raw_endpoint_name": "Cmax",
        "raw_value": "70.4",
        "raw_unit": "ng/mL",
        "species": "HUMAN",
        "source_database": "Drugs@FDA",
        "reference_text": "Drugs@FDA NDA215310 · Table 2.b Clinical Pharmacology (Day 29, 160 mg QD)",
        "assay_conditions_json": {
            "species": "HUMAN", "route": "ORAL", "dose": 160.0, "dose_unit": "mg", "regimen": "QD", "day": 29
        }
    }
    dec_day29 = qualify_evidence_record_v51(rec_day29)
    assert dec_day29.canonical_endpoint_id == "HUMAN_PK_CMAX_ORAL_DAY29"
    assert "Day 29" in dec_day29.display_name
    assert dec_day29.normalized_value == 70.4


def test_zero_prediction_proximity_influence():
    """Representative selection never accepts or uses predicted numerical values."""
    item1 = {
        "id": 1,
        "origin": "AUTO_QUALIFIED_EXTERNAL",
        "comparability": "DIRECTLY_COMPARABLE",
        "display": {"value": 10.0, "unit": "nM"},
        "context": {"species": "HUMAN", "target": "EGFR Exon20ins"}
    }
    item2 = {
        "id": 2,
        "origin": "AUTO_QUALIFIED_EXTERNAL",
        "comparability": "DIRECTLY_COMPARABLE",
        "display": {"value": 1000.0, "unit": "nM"},
        "context": {"species": "MOUSE", "target": "EGFR Exon20ins"}
    }
    # Rank depends purely on species / context, not prediction
    rank1 = representative_rank(item1)
    rank2 = representative_rank(item2)
    assert rank1 < rank2


def test_database_endpoint_comparison_orforglipron_mobocertinib():
    """Integration check on actual database records for Orforglipron and Mobocertinib."""
    db = SessionLocal()
    try:
        orf = db.query(Compound).filter(Compound.name.ilike("%Orforglipron%")).first()
        assert orf is not None
        res_orf = build_endpoint_comparison(db, orf.versions[-1].id)
        ec50_row = next((r for r in res_orf["scientific_rows"] if r["section"] == "ACTIVITY" and r.get("display_name") == "EC50"), None)
        assert ec50_row is not None
        assert ec50_row["primary_experimental_display"]["value"] == 1.2
        assert ec50_row["primary_experimental_display"]["unit"] == "nM"

        mob = db.query(Compound).filter(Compound.name.ilike("%Mobocertinib%")).first()
        assert mob is not None
        res_mob = build_endpoint_comparison(db, mob.versions[-1].id)
        cmax_day1 = next((r for r in res_mob["scientific_rows"] if "Day 1" in r.get("display_name", "") and "Cmax" in r.get("display_name", "")), None)
        assert cmax_day1 is not None
        assert cmax_day1["primary_experimental_display"]["value"] == 77.9

        cmax_day29 = next((r for r in res_mob["scientific_rows"] if "Day 29" in r.get("display_name", "") and "Cmax" in r.get("display_name", "")), None)
        assert cmax_day29 is not None
        assert cmax_day29["primary_experimental_display"]["value"] == 70.4
    finally:
        db.close()
