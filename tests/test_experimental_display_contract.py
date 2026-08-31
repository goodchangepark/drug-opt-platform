from backend.experimental_display import (
    COMPARABLE_AFTER_DETERMINISTIC_CONVERSION, CONDITIONALLY_COMPARABLE,
    DIRECTLY_COMPARABLE, NOT_COMPARABLE, RELATED_NOT_SAME_ENDPOINT,
    NORMALIZATION_VERSION, normalize_experimental,
)


def test_activity_types_are_preserved_and_classifier_potency_is_related_only():
    row = normalize_experimental("IC50", 10, "nM", target="CYP3A4")
    assert row["comparability_status"] == RELATED_NOT_SAME_ENDPOINT
    assert row["raw_value"] == "10" and row["raw_unit"] == "nM"


def test_solubility_molar_and_mass_rules_preserve_context():
    molar = normalize_experimental("aqueous solubility", 10, "µmol/L", conditions="pH 7.4")
    assert molar["comparability_status"] == CONDITIONALLY_COMPARABLE
    assert molar["normalized_value"] == -5
    mass = normalize_experimental("solubility", 1, "mg/mL", mw=180)
    assert mass["normalized_value"] is not None
    assert mass["normalization_version"] == NORMALIZATION_VERSION
    assert normalize_experimental("solubility", 1, "mg/mL")["comparability_status"] == NOT_COMPARABLE


def test_caco2_direction_and_log_conversion():
    ab = normalize_experimental("Caco-2 Papp A->B", 1e-5, "cm/s")
    assert ab["comparability_status"] == COMPARABLE_AFTER_DETERMINISTIC_CONVERSION
    assert ab["normalized_value"] == -5
    assert normalize_experimental("Caco-2 Papp B->A", 1e-5, "cm/s")["comparability_status"] == NOT_COMPARABLE
    assert normalize_experimental("PAMPA", 1e-5, "cm/s")["comparability_status"] == NOT_COMPARABLE


def test_ppb_species_and_fu_conversion():
    assert normalize_experimental("PPB", 95, "% bound", species="Human")["comparability_status"] == DIRECTLY_COMPARABLE
    fu = normalize_experimental("PPB fu", .05, "fraction", species="Human")
    assert fu["normalized_value"] == 95
    assert normalize_experimental("PPB", 95, "% bound", species="Rat")["comparability_status"] == NOT_COMPARABLE


def test_classifier_and_logd_guardrails():
    assert normalize_experimental("hERG IC50", 2.3, "µM")["comparability_status"] == RELATED_NOT_SAME_ENDPOINT
    assert normalize_experimental("logD", 2, "logD", conditions="pH 7.4")["comparability_status"] == DIRECTLY_COMPARABLE
    assert normalize_experimental("logP", 2, "logP")["comparability_status"] == RELATED_NOT_SAME_ENDPOINT


def test_clearance_isolation_and_categorical_ames_support():
    assert normalize_experimental("HLM intrinsic clearance", 20, "µL/min/mg protein", species="Human")["comparability_status"] != DIRECTLY_COMPARABLE
    assert normalize_experimental("RLM intrinsic clearance", 20, "µL/min/mg protein", species="Rat")["comparability_status"] != DIRECTLY_COMPARABLE
    assert normalize_experimental("Ames mutagenicity", "Positive")["comparability_status"] == DIRECTLY_COMPARABLE


def test_pka_microstate_is_condition_dependent_and_relations_are_not_changed():
    pka = normalize_experimental("micro-pKa", 6.2, "pKa")
    assert pka["comparability_status"] == CONDITIONALLY_COMPARABLE
    assert normalize_experimental("IC50", 2, "µM")["raw_value"] == "2"
