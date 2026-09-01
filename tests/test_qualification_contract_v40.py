from backend.qualification_contract import (
    ADAPTATION_ELIGIBLE, CONTEXT_QUALIFIED, DIRECTLY_COMPARABLE,
    ENDPOINT_QUALIFIED, IMPORTABLE, PREDICTION_PAIRABLE,
    QUALIFICATION_VERSION, RELATED_SAME_GROUP, aggregate_qualification,
    qualify_record,
)


def _record(endpoint="PPB", value="91.2", *, comparison="DIRECTLY_COMPARABLE", source="FDA / Regulatory", imported=False):
    return {
        "source": source, "endpoint": endpoint, "value": value,
        "raw_value": value, "unit": "% bound", "raw_unit": "% bound",
        "identity_match_status": "EXACT_STRUCTURE_MATCH",
        "reference_status": "REFERENCE_RESOLVED_REGULATORY",
        "canonical_endpoint_id": "HUMAN_PPB" if endpoint == "PPB" else "HUMAN_PK_CMAX_UNSPECIFIED", "endpoint_qualified": True,
        "comparability_status": comparison,
        "importable": not imported, "state": "EXTERNAL_IMPORTED" if imported else "EXTERNAL_CANDIDATE",
    }


def test_contract_distinguishes_endpoint_qualification_from_pairability():
    q = qualify_record(_record(), prediction_endpoints=())
    assert q["qualification_version"] == QUALIFICATION_VERSION
    assert q["stages"][ENDPOINT_QUALIFIED]
    assert q["stages"][CONTEXT_QUALIFIED]
    assert not q["stages"][PREDICTION_PAIRABLE]
    assert not q["stages"][DIRECTLY_COMPARABLE]
    assert q["primary_gap_reason"] == "NO_CURRENT_PREDICTION_ENDPOINT"


def test_pairability_and_direct_comparison_require_actual_prediction():
    q = qualify_record(_record(), prediction_endpoints={"HUMAN_PPB"})
    assert q["stages"][PREDICTION_PAIRABLE]
    assert q["stages"][DIRECTLY_COMPARABLE]
    assert q["comparability_status"] == DIRECTLY_COMPARABLE
    assert q["stages"][IMPORTABLE]
    imported = qualify_record(_record(imported=True), prediction_endpoints={"HUMAN_PPB"})
    assert imported["stages"][ADAPTATION_ELIGIBLE] is False


def test_related_same_group_is_not_direct():
    q = qualify_record(_record(comparison="RELATED_NOT_SAME_ENDPOINT"), prediction_endpoints={"HUMAN_PPB"})
    assert q["stages"][PREDICTION_PAIRABLE]
    assert q["stages"][RELATED_SAME_GROUP]
    assert not q["stages"][DIRECTLY_COMPARABLE]
    assert q["primary_gap_reason"] == "RELATED_MEASUREMENT_SEMANTICS"


def test_source_aggregation_is_stage_explicit_and_reconciles():
    records = [_record(), _record(endpoint="Cmax", value="412", source="FDA / Regulatory", comparison="UNSUPPORTED")]
    result = aggregate_qualification(records, prediction_endpoints={"HUMAN_PPB"}, raw_source_counts={"FDA / Regulatory": 3})
    source = result["sources"]["FDA / Regulatory"]
    assert source["found"] == 3
    assert source["unique"] == 2
    assert source["endpoint_qualified"] == 2
    assert source["prediction_pairable"] == 1
    assert source["direct"] == 1
    assert "qualified" not in source


def test_needs_review_has_precise_reason():
    row = _record(endpoint="Caco-2 Papp", comparison="CONDITIONALLY_COMPARABLE")
    row["display"] = {"reason": "Caco-2 direction is not recorded", "canonical_endpoint_id": "CACO2_PAPP_AB"}
    q = qualify_record(row, prediction_endpoints={"CACO2_PAPP_AB"})
    assert not q["stages"][CONTEXT_QUALIFIED]
    assert q["primary_gap_reason"] == "DIRECTION_MISSING"
