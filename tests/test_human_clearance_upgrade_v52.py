from backend.clearance_architecture import (
    CL_HEPATIC, CL_ORAL_APPARENT, CL_RENAL, CL_TOTAL_IV, CL_UNRESOLVED,
    TOTAL_CL_INCOMPLETE, ClearanceComponent, assemble_total_clearance,
    classify_clearance_semantics, conversion_cl_l_h_kg_to_ml_min_kg,
    renal_readiness,
)


def test_clearance_semantics_are_not_interchanged():
    assert classify_clearance_semantics("CL", "IV") == CL_TOTAL_IV
    assert classify_clearance_semantics("HLM Clint", "IV") == CL_HEPATIC
    assert classify_clearance_semantics("CL/F", "ORAL") == CL_ORAL_APPARENT
    assert classify_clearance_semantics("renal clearance", "IV") == CL_RENAL
    assert classify_clearance_semantics("clearance", "") == CL_UNRESOLVED


def test_total_clearance_fails_closed_without_renal_and_other_components():
    hepatic = ClearanceComponent("CL_H", 2.0, "mL/min/kg", "IVIVE", {"source": "MODEL"})
    result = assemble_total_clearance(hepatic=hepatic)
    assert result["status"] == TOTAL_CL_INCOMPLETE
    assert result["value"] is None
    assert "CL_R" in result["missing_components"]
    assert result["assumptions"] == []


def test_total_clearance_sums_only_explicit_components():
    c = lambda n, v: ClearanceComponent(n, v, "mL/min/kg", "EXPERIMENTAL", {})
    result = assemble_total_clearance(hepatic=c("CL_H", 1), renal=c("CL_R", 2), other=c("CL_OTHER", 3))
    assert result["status"] == "TOTAL_CL_COMPLETE"
    assert result["value"] == 6


def test_renal_unknown_never_becomes_zero():
    result = renal_readiness()
    assert result["status"] == "RENAL_UNRESOLVED"
    assert result["renal_clearance"] is None
    assert result["silent_zero_assumption"] is False


def test_clearance_unit_conversion_is_dimensional():
    assert conversion_cl_l_h_kg_to_ml_min_kg(1.2) == 20.0
