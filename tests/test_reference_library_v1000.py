import json


def test_reference_library_v1_has_identity_safe_checkpoints():
    with open("validation/reference_library_v1_1000.json") as handle:
        artifact = json.load(handle)
    assert artifact["target_n"] == 1000
    assert artifact["identity_qc"]["qualified_compounds"] == 1000
    assert artifact["identity_qc"]["duplicate_inchikey"] == 0
    assert [artifact["checkpoints"][str(n)]["n"] for n in (500, 750, 1000)] == [500, 750, 1000]


def test_external_missing_identifiers_are_not_fabricated():
    with open("validation/reference_library_v1_1000.json") as handle:
        artifact = json.load(handle)
    external = [x for x in artifact["compounds"] if x["role"] == "DIVERSITY_SELECTED_EXTERNAL"]
    assert len(external) == 750
    assert all(x["inchikey"] and x["smiles"] for x in external)
    assert all(x["drugbank_id"] is None for x in external)
