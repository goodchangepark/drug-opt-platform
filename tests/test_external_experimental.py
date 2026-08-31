from backend import external_experimental as ex


def test_cas_checksum_states():
    assert ex.cas_status("50-78-2") == "VALID"
    assert ex.cas_status("50-78-3") == "INVALID_CHECKSUM"
    assert ex.cas_status("") == "EMPTY"


def test_chembl_enriches_document_and_preserves_measurements(monkeypatch):
    def fake(url):
        if "molecule.json" in url:
            return {"molecules": [{"molecule_chembl_id": "CHEMBL1"}]}
        if "activity.json" in url:
            return {"activities": [
                {"activity_id": 10, "assay_chembl_id": "CHEMBL_ASSAY1", "document_chembl_id": "CHEMBL_DOC1",
                 "target_chembl_id": "CHEMBL_T1", "standard_type": "IC50", "standard_value": "12", "standard_units": "nM",
                 "standard_relation": "<", "assay_type": "B"},
                {"activity_id": 11, "standard_type": "Ki", "standard_value": "4", "standard_units": "nM"},
            ]}
        if "assay/" in url:
            return {"assay_chembl_id": "CHEMBL_ASSAY1", "description": "binding assay", "assay_type": "B", "cell_line_name": "HEK293"}
        if "target/" in url:
            return {"target_chembl_id": "CHEMBL_T1", "pref_name": "Target", "organism": "Homo sapiens"}
        if "document/" in url:
            return {"document_chembl_id": "CHEMBL_DOC1", "title": "Paper", "journal": "Journal", "year": 2024, "doi": "10.1000/test", "pubmed_id": 123}
        return {}
    monkeypatch.setattr(ex, "_get_json", fake)
    rows = ex._chembl_activities("AAAA-BBBBBB-C", "50-78-2")
    assert rows[0]["doi"] == "10.1000/test"
    assert rows[0]["pmid"] == "123"
    assert rows[0]["reference_status"] == "REFERENCE_RESOLVED_DOI"
    assert rows[0]["cell_line"] == "HEK293"
    assert rows[0]["relation"] == "<"
    assert rows[1]["endpoint"] == "Ki"


def test_private_smiles_not_used_in_lookup(monkeypatch):
    seen = []
    monkeypatch.setattr(ex, "_get_json", lambda url: (seen.append(url) or {"_not_found": True}))
    ex.lookup("50-78-2", "PUBLIC-INCHIKEY")
    assert seen and "PUBLIC-INCHIKEY" not in seen[0]
    assert all("C=" not in url for url in seen)


def test_pubchem_annotations_require_reference(monkeypatch):
    monkeypatch.setattr(ex, "_get_json", lambda url: {"Record": {"Section": [{"TOCHeading": "Experimental Properties", "Information": [
        {"Name": "Solubility", "Value": {"StringWithMarkup": [{"String": "1 mg/mL"}]}, "Reference": [{"SourceName": "Paper"}]},
        {"Name": "Computed X", "Value": {"StringWithMarkup": [{"String": "2"}]}}
    ]}]}})
    rows = ex._pubchem_experimental_annotations(1, "50-78-2")
    assert rows[0]["reference_status"] == "REFERENCE_RESOLVED"
    assert len(rows) == 2
