"""
Drug-OPT — Global Experimental Evidence Qualification Test Suite v5.1
Policy Version: drugopt-evidence-qualification-v5.1
"""

import pytest
from backend.evidence_qualification_v51 import (
    qualify_evidence_record_v51,
    resolve_species_v51,
    resolve_route_v51,
    resolve_analyte_v51,
    resolve_target_context_v51,
    parse_numeric_strict,
    STATE_AUTO_QUALIFIED,
    STATE_RELATED,
    STATE_REVIEW_REQUIRED,
    STATE_UNUSABLE,
    REASON_LITERATURE_CITATION_ONLY,
    REASON_NON_NUMERIC_OBSERVATION,
    REASON_TOC_OR_FOOTNOTE_ARTIFACT,
    REASON_RELATIVE_RATIO_NOT_ABSOLUTE_EXPOSURE,
)


def test_evidence_funnel_auto_qualified_flow():
    """Verify 7-stage evidence funnel progression for directly auto-qualified measurement."""
    record = {
        "raw_endpoint_name": "IC50",
        "raw_value": "0.4",
        "raw_unit": "nM",
        "species": "Homo sapiens",
        "source_database": "ChEMBL",
        "assay_conditions_json": {"target": "EGFR Exon20ins", "assay": "Biochemical kinase"},
    }
    decision = qualify_evidence_record_v51(record)
    
    assert decision.funnel["source_found"] is True
    assert decision.funnel["raw_evidence"] is True
    assert decision.funnel["observation_extracted"] is True
    assert decision.funnel["endpoint_classified"] is True
    assert decision.funnel["unit_normalized"] is True
    assert decision.funnel["qualification_state"] == STATE_AUTO_QUALIFIED
    assert decision.funnel["displayed"] is True
    assert decision.funnel["drop_stage"] is None
    assert decision.funnel["drop_reason"] is None
    
    assert decision.canonical_endpoint_id == "ACTIVITY_IC50"
    assert decision.normalized_value == 0.4
    assert decision.normalized_unit == "nM"
    assert decision.evidence_state == STATE_AUTO_QUALIFIED
    assert decision.comparability_status == "DIRECTLY_COMPARABLE"


def test_evidence_funnel_literature_citation_drop():
    """Verify literature candidate without numeric value is captured in stage 3 with explicit reason."""
    record = {
        "raw_endpoint_name": "Literature candidate",
        "raw_value": "Title: EGFR inhibitors review",
        "raw_unit": "",
        "species": "Human",
        "source_database": "Europe PMC",
        "assay_conditions_json": {},
    }
    decision = qualify_evidence_record_v51(record)
    
    assert decision.funnel["source_found"] is True
    assert decision.funnel["raw_evidence"] is True
    assert decision.funnel["observation_extracted"] is False
    assert decision.funnel["drop_stage"] == "OBSERVATION_EXTRACTED"
    assert decision.funnel["drop_reason"] == REASON_LITERATURE_CITATION_ONLY
    assert decision.evidence_state == STATE_REVIEW_REQUIRED
    assert decision.unresolved_reason == REASON_LITERATURE_CITATION_ONLY


def test_evidence_funnel_toc_footnote_artifact_drop():
    """Verify TOC / Footnote index text is dropped as UNUSABLE."""
    record = {
        "raw_endpoint_name": "Metabolite",
        "raw_value": "6",
        "raw_unit": "H",
        "species": "Human",
        "source_database": "FDA / Regulatory",
        "assay_conditions_json": {"conditions": "16&/& ZLWK (*)5 ([RQ LQV LQ 6HFRQG-line Settings ............................. 26 Table 2"},
    }
    decision = qualify_evidence_record_v51(record)
    
    assert decision.funnel["drop_stage"] == "OBSERVATION_EXTRACTED"
    assert decision.funnel["drop_reason"] == REASON_TOC_OR_FOOTNOTE_ARTIFACT
    assert decision.evidence_state == STATE_UNUSABLE


def test_pk_absolute_vs_ddi_relative_distinction():
    """Verify that absolute PK exposure and relative DDI % changes are correctly classified."""
    # 1. Absolute PK Cmax
    abs_pk = {
        "raw_endpoint_name": "Cmax",
        "raw_value": "60",
        "raw_unit": "ng/mL",
        "species": "Human",
        "source_database": "FDA / Regulatory",
        "assay_conditions_json": {"conditions": "Mean steady state Cmax following 160 mg once daily oral administration"},
    }
    d_abs = qualify_evidence_record_v51(abs_pk)
    assert d_abs.section == "PK"
    assert d_abs.canonical_endpoint_id == "HUMAN_PK_CMAX_ORAL"
    assert d_abs.normalized_value == 60.0
    assert d_abs.normalized_unit == "ng/mL"
    assert d_abs.evidence_state == STATE_AUTO_QUALIFIED

    # 2. DDI relative % change
    ddi_pk = {
        "raw_endpoint_name": "Cmax",
        "raw_value": "32",
        "raw_unit": "%",
        "species": "Human",
        "source_database": "FDA / Regulatory",
        "assay_conditions_json": {"conditions": "Itraconazole increased Cmax by 32% following concomitant administration"},
    }
    d_ddi = qualify_evidence_record_v51(ddi_pk)
    assert d_ddi.section == "PK"
    assert d_ddi.canonical_endpoint_id == "HUMAN_PK_DDI_RELATIVE_RATIO"
    assert d_ddi.evidence_state == STATE_RELATED
    assert d_ddi.unresolved_reason == REASON_RELATIVE_RATIO_NOT_ABSOLUTE_EXPOSURE


def test_unit_normalization_conversions():
    """Verify scientific unit normalization across potencies, PPB, Caco-2, and PK."""
    # µM to nM
    d1 = qualify_evidence_record_v51({
        "raw_endpoint_name": "Ki",
        "raw_value": "0.015",
        "raw_unit": "µM",
        "source_database": "ChEMBL",
        "assay_conditions_json": {},
    })
    assert d1.normalized_value == 15.0
    assert d1.normalized_unit == "nM"

    # fu to % bound
    d2 = qualify_evidence_record_v51({
        "raw_endpoint_name": "protein binding",
        "raw_value": "0.004",
        "raw_unit": "fu",
        "species": "Human",
        "source_database": "FDA / Regulatory",
        "assay_conditions_json": {},
    })
    assert d2.normalized_value == 99.6
    assert d2.normalized_unit == "% bound"

    # Caco-2 10^-6 cm/s
    d3 = qualify_evidence_record_v51({
        "raw_endpoint_name": "Papp A->B",
        "raw_value": "12.5",
        "raw_unit": "10^-6 cm/s",
        "source_database": "ChEMBL",
        "assay_conditions_json": {"conditions": "Caco-2 apical to basolateral"},
    })
    assert d3.normalized_unit == "log10(cm/s)"
    assert round(d3.normalized_value, 2) == -4.90


def test_metabolism_cyp_and_excretion_classification():
    """Verify CYP isoforms and fecal/urinary excretion recovery are classified into METABOLISM."""
    # Fecal excretion %
    d_feces = qualify_evidence_record_v51({
        "raw_endpoint_name": "feces",
        "raw_value": "79",
        "raw_unit": "%",
        "species": "Human",
        "source_database": "FDA / Regulatory",
        "assay_conditions_json": {"conditions": "79% recovered in feces"},
    })
    assert d_feces.section == "METABOLISM"
    assert d_feces.canonical_endpoint_id == "EXCRETION_FECAL"
    assert d_feces.normalized_value == 79.0
    assert d_feces.normalized_unit == "% dose"

    # CYP3A4 metabolic contribution
    d_cyp = qualify_evidence_record_v51({
        "raw_endpoint_name": "CYP3A4",
        "raw_value": "54",
        "raw_unit": "%",
        "species": "Human",
        "source_database": "FDA / Regulatory",
        "assay_conditions_json": {"conditions": "primarily metabolized by CYP3A4"},
    })
    assert d_cyp.section == "METABOLISM"
    assert d_cyp.canonical_endpoint_id == "CYP3A4_METABOLIC_CONTRIBUTION"
    assert d_cyp.evidence_state == STATE_AUTO_QUALIFIED
