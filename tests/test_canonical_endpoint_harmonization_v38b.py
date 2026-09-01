from math import isclose

from backend.canonical_endpoints import (
    COMPARISON_UNIT_VERSION,
    CANONICAL_ENDPOINT_VERSION,
    CONVERTED,
    DIRECT,
    RELATED,
    canonicalize_prediction_endpoint,
    normalize_experimental_observation,
    registry_report,
)


def test_semantic_aliases_join_to_same_canonical_endpoint():
    assert normalize_experimental_observation(
        "Plasma protein binding", 0.083, "fu", species="Human"
    )["canonical_endpoint_id"] == canonicalize_prediction_endpoint("PPB", species="Human")["canonical_endpoint_id"] == "HUMAN_PPB"
    assert normalize_experimental_observation(
        "Caco-2 Papp A-B", 9.3e-6, "cm/s", context={"direction": "A->B"}
    )["canonical_endpoint_id"] == "CACO2_PAPP_AB"
    assert normalize_experimental_observation(
        "human microsomal intrinsic clearance", 12, "mL/min/kg", species="Human"
    )["canonical_endpoint_id"] == "HLM_CLINT"
    assert normalize_experimental_observation(
        "terminal half-life", 120, "min", species="Rat", context={"route": "PO"}
    )["canonical_endpoint_id"] == "RAT_PK_T_HALF_ORAL"
    assert normalize_experimental_observation(
        "Apparent oral clearance", 1, "L/h/kg", species="Rat", context={"route": "oral"}
    )["canonical_endpoint_id"] == "RAT_PK_CLF_ORAL_ORAL"


def test_ppb_and_caco2_units_are_deterministic():
    ppb = normalize_experimental_observation("PPB", 0.083, "fu", species="Human")
    assert isclose(ppb["normalized_value"], 91.7)
    assert ppb["comparability_status"] == CONVERTED
    caco = normalize_experimental_observation("Caco-2 Papp A-B", 9.3e-6, "cm/s")
    assert isclose(caco["normalized_value"], -5.031517051446065)
    assert caco["normalized_unit"] == "log10(cm/s)"
    assert caco["comparability_status"] == CONVERTED


def test_pk_units_and_context_are_not_merged_unsafely():
    cl = normalize_experimental_observation("CL", 10, "mL/min/kg", species="Rat", context={"route": "IV"})
    clf = normalize_experimental_observation("CL/F", 10, "mL/min/kg", species="Rat", context={"route": "PO"})
    assert cl["canonical_endpoint_id"] == "RAT_PK_CL_IV"
    assert clf["canonical_endpoint_id"] == "RAT_PK_CLF_ORAL_ORAL"
    assert cl["canonical_endpoint_id"] != clf["canonical_endpoint_id"]
    assert normalize_experimental_observation("CL", 10, "mL/min/kg", species="Human", context={"route": "IV"})["canonical_endpoint_id"] != cl["canonical_endpoint_id"]
    assert isclose(normalize_experimental_observation("CL", 1, "L/h/kg", species="Rat", context={"route": "IV"})["normalized_value"], 16.6666666667)
    assert isclose(normalize_experimental_observation("F", 0.4, "fraction", species="Rat", context={"route": "PO"})["normalized_value"], 40.0)


def test_pk_display_unit_contract_and_analyte_identity():
    cmax = normalize_experimental_observation("Cmax", 1, "µg/L", species="Rat", context={"route": "IV"})
    auc = normalize_experimental_observation("AUC0-inf", 1, "µg*h/L", species="Rat", context={"route": "IV"})
    half_life = normalize_experimental_observation("t1/2", 30, "min", species="Rat", context={"route": "IV"})
    metabolite = normalize_experimental_observation("Cmax metabolite", 1, "ng/mL", species="Rat", context={"route": "IV"})
    assert cmax["normalized_value"] == 1 and cmax["normalized_unit"] == "ng/mL"
    assert auc["normalized_value"] == 1 and auc["normalized_unit"] == "ng*h/mL"
    assert half_life["normalized_value"] == 0.5 and half_life["normalized_unit"] == "hours"
    assert metabolite["analyte"] == "METABOLITE"
    assert metabolite["comparison_key"] != cmax["comparison_key"]


def test_hepatocyte_and_cyp_semantics_remain_distinct():
    hlm = normalize_experimental_observation("HLM Clint", 12, "mL/min/kg", species="Human")
    hepatocyte = normalize_experimental_observation("CLH", 12, "µL/min/10^6 cells", species="Human", context={"matrix": "hepatocytes"})
    assert hlm["canonical_endpoint_id"] == "HLM_CLINT"
    assert hepatocyte["canonical_endpoint_id"] == "HEPATOCYTE_CLINT"
    cyp_exp = normalize_experimental_observation("CYP3A4 IC50", 18.2, "µM", species="Human")
    cyp_pred = canonicalize_prediction_endpoint("CYP3A4 inhibitor", species="Human")
    assert cyp_exp["canonical_endpoint_id"] == cyp_pred["canonical_endpoint_id"] == "CYP3A4_INHIBITION"
    assert cyp_exp["comparability_status"] == RELATED
    assert canonicalize_prediction_endpoint("CYP1A2 inhibitor")["canonical_endpoint_id"] == "CYP1A2_INHIBITION"


def test_registry_is_versioned():
    report = registry_report()
    assert report["canonical_endpoint_version"] == CANONICAL_ENDPOINT_VERSION
    assert report["comparison_unit_version"] == COMPARISON_UNIT_VERSION
    assert any(row["canonical_endpoint_id"] == "HUMAN_PPB" for row in report["endpoints"])
