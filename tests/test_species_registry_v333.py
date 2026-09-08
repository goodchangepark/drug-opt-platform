from backend.canonical_endpoints import normalize_species
from backend.main import get_species_registry
from backend.pk_parameter_set import calculate_well_stirred_clearance


def test_v333_species_registry_has_required_pk_species():
    payload = get_species_registry()
    assert payload["registry_version"] == "DRUGOPT_SPECIES_V1"
    assert payload["canonical_species"] == ["HUMAN", "RAT", "MOUSE", "DOG", "MONKEY"]
    assert all(row["physiology_available"] for row in payload["species"])


def test_species_aliases_normalize_without_cross_species_substitution():
    assert normalize_species("Homo sapiens") == "HUMAN"
    assert normalize_species("Rattus norvegicus") == "RAT"
    assert normalize_species("Mus musculus") == "MOUSE"
    assert normalize_species("beagle dog") == "DOG"
    assert normalize_species("cynomolgus monkey") == "MONKEY"
    assert normalize_species("") == "UNSPECIFIED"
    assert normalize_species("Candida albicans") == "OTHER"


def test_well_stirred_physiology_is_species_specific_and_never_human_fallback():
    dog = calculate_well_stirred_clearance(20.0, 0.1, species="DOG")
    monkey = calculate_well_stirred_clearance(20.0, 0.1, species="MONKEY")
    human = calculate_well_stirred_clearance(20.0, 0.1, species="HUMAN")
    assert dog["species"] == "DOG"
    assert monkey["species"] == "MONKEY"
    assert dog["cl_h_plasma_ml_min_kg"] != human["cl_h_plasma_ml_min_kg"]
    assert monkey["cl_h_plasma_ml_min_kg"] != human["cl_h_plasma_ml_min_kg"]


def test_unknown_species_fails_closed_instead_of_using_human_physiology():
    import pytest

    with pytest.raises(ValueError, match="Unsupported species physiology"):
        calculate_well_stirred_clearance(20.0, 0.1, species="FERRET")
