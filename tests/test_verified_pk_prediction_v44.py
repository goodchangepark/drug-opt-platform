from backend.pk_engine_v1 import INSUFFICIENT_INPUT, estimate_one_compartment, request_fingerprint
from pathlib import Path


def test_pk_overlay_refuses_silent_defaults_and_incomplete_oral_f():
    result=estimate_one_compartment(species="HUMAN", route="ORAL", inputs={"dose_mg_per_kg":1,"cl_l_per_h_per_kg":1,"v_l_per_kg":1,"sources":{}})
    assert result["status"]==INSUFFICIENT_INPUT
    assert "f_fraction" in result["missing_inputs"]


def test_complete_experiment_informed_oral_context_is_mechanistic_not_model():
    values={"dose_mg_per_kg":1,"f_fraction":0.5,"ka_per_h":1,"cl_l_per_h_per_kg":1,"v_l_per_kg":1,
            "sources":{"dose_mg_per_kg":"USER_SUPPLIED","f_fraction":"EXPERIMENTAL_IMPORTED","ka_per_h":"USER_SUPPLIED","cl_l_per_h_per_kg":"EXPERIMENTAL_INTERNAL","v_l_per_kg":"EXPERIMENTAL_INTERNAL"}}
    result=estimate_one_compartment(species="RAT",route="ORAL",inputs=values)
    assert result["status"]=="COMPLETE"
    assert result["outputs"]["t_half"]["prediction_type"]=="DERIVED_ESTIMATE"
    assert result["outputs"]["cmax"]["prediction_type"]=="MECHANISTIC_ESTIMATE"


def test_pk_request_fingerprint_is_deterministic_and_context_sensitive():
    a={"species":"RAT","route":"ORAL","dose":1}
    assert request_fingerprint(a)==request_fingerprint(dict(a))
    assert request_fingerprint(a)!=request_fingerprint(a|{"dose":2})


def test_review_evidence_is_separated_from_primary_scientific_tables():
    js=(Path(__file__).resolve().parents[1] / "frontend" / "static" / "app.js").read_text()
    assert "function scientificReviewQueue(rows)" in js
    assert "Evidence Requiring Review" in js
    assert "const primaryRows=scientificRows.filter(row=>!scientificReviewRow(row));" in js
