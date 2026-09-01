from backend import experimental_harvester as h
import json
from pathlib import Path


def test_public_identity_does_not_require_or_send_local_structure(monkeypatch):
    urls = []
    monkeypatch.setattr(h, "_get_json", lambda url: (urls.append(url) or {"PropertyTable": {"Properties": [{"CID": 2244, "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N", "CanonicalSMILES": "CC(=O)Oc1ccccc1C(=O)O"}]}}))
    identity = h.resolve_public_identity(cas="50-78-2", local_inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
    assert identity.pubchem_cid == "2244" and identity.identity_status == "EXACT_STRUCTURE_MATCH"
    assert all("CC(=O)" not in url for url in urls)


def test_chembl_keeps_assay_type_and_extended_measurements(monkeypatch):
    def fake(url):
        if "molecule.json" in url: return {"molecules": [{"molecule_chembl_id": "CHEMBL1"}]}
        return {"activities": [{"activity_id": 1, "assay_type": "A", "standard_type": "Solubility", "standard_value": "12", "standard_units": "uM"}, {"activity_id": 2, "assay_type": "X", "standard_type": "IC50", "standard_value": "1"}]}
    monkeypatch.setattr(h, "_get_json", fake)
    rows = h.ChEMBLAdapter().harvest(h.PublicIdentity(inchikey="AAAA"))
    assert len(rows) == 1 and rows[0]["assay_type"] == "A" and rows[0]["endpoint"] == "Solubility"


def test_comptox_is_safe_without_api_key(monkeypatch):
    monkeypatch.delenv("COMPTOX_API_KEY", raising=False)
    assert h.CompToxAdapter().status() == "NOT_CONFIGURED"
    assert h.CompToxAdapter().harvest(h.PublicIdentity(cas="50-78-2")) == []


def test_deduplication_marks_without_silently_dropping():
    identity = h.PublicIdentity(inchikey="X")
    rows = h.deduplicate(identity, [h._record("A", "1", "IC50", 2, "nM", doi="x"), h._record("B", "2", "IC50", 2, "nM", doi="x")])
    assert rows[0]["duplicate_status"] == "DISTINCT_MEASUREMENT"
    assert rows[1]["duplicate_status"] == "SAME_MEASUREMENT"


def test_public_identity_expands_verified_pubchem_aliases(monkeypatch):
    def fake(url):
        if "synonyms" in url:
            return {"InformationList": {"Information": [{"Synonym": ["DZD9008", "DZD-9008", "ZEGFROVY", "CHEMBL5314564", "DTXSID701376536", "L1Q2K5JYO8"]}]}}
        return {"PropertyTable": {"Properties": [{"CID": 1, "InChIKey": "KEY", "CanonicalSMILES": "CC"}]}}
    monkeypatch.setattr(h, "_get_json", fake)
    identity = h.resolve_public_identity(name="public drug")
    assert {"DZD9008", "ZEGFROVY"}.issubset(identity.synonyms)
    assert identity.chembl_id == "CHEMBL5314564" and identity.dtxsid == "DTXSID701376536" and identity.unii == "L1Q2K5JYO8"


def test_regulatory_documents_and_candidates_remain_separate(monkeypatch):
    adapter = h.RegulatoryAdapter()
    monkeypatch.setattr(h, "_get_document_text", lambda _url: "Cmax 412 ng/mL at steady state.")
    monkeypatch.setattr(h, "_get_json", lambda _url: {"results": [{"application_number": "NDA123456", "submissions": [{"submission_type": "ORIG", "submission_number": "1", "submission_status": "AP", "submission_status_date": "20250101", "application_docs": [{"id": "1", "url": "https://example.test/label.pdf", "type": "Label", "date": "20250101"}]}]}]})
    rows = adapter.harvest(h.PublicIdentity(name="PUBLIC"))
    assert any(row["record_status"] == "DOCUMENT_DISCOVERED" for row in rows)
    assert any(row["record_status"] == "REGULATORY_CANDIDATE" and row["value"] for row in rows)


def test_regulatory_table_context_prefers_ppb_over_incidental_cmax():
    endpoint, value, unit = h.RegulatoryAdapter.endpoint_value_from_context(
        "91.62%, respectively, bound to human plasma protein at 1 µM; mean Cmax was 412 ng/mL"
    )
    assert (endpoint, value, unit) == ("protein binding", "91.62", "%")


def test_supplement_csv_parser_is_in_memory_and_reports_parse_state(monkeypatch):
    monkeypatch.setattr(h, "_get_document_bytes", lambda _url, max_bytes=0: (b"endpoint,value\nCmax,12 ng/mL\n", "text/csv"))
    text, file_type, state = h._supplement_text("https://example.test/table.csv")
    assert file_type == "CSV" and state == "SUPPLEMENT_PARSED" and "Cmax" in text


def test_supplement_download_failure_is_explicit(monkeypatch):
    monkeypatch.setattr(h, "_get_document_bytes", lambda _url, max_bytes=0: (b"", ""))
    _text, _file_type, state = h._supplement_text("https://example.test/table.pdf")
    assert state == "SUPPLEMENT_DOWNLOAD_FAILED"


def test_nmpa_official_approval_is_not_inferred_from_a_generic_hit(monkeypatch):
    monkeypatch.setattr(h, "_get_json", lambda _url: {"data": [{"id": 1, "title": "Public drug approved for marketing", "abstractdesc": "Official approval notice", "pubUrl": "https://nmpa.test/1", "pubTime": "2023-08-23"}]})
    identity = h.PublicIdentity(name="Public drug")
    rows = h.NMPAAdapter().harvest(identity)
    assert identity.approval["NMPA"]["status"] == "APPROVED"
    assert rows[0]["record_status"] == "NMPA_APPROVAL_CONFIRMED_DOCUMENT_NOT_PUBLICLY_ACCESSIBLE"


def test_approved_drug_coverage_artifact_has_five_public_controls():
    artifact = json.loads((Path(__file__).parents[1] / "validation" / "approved_drug_harvester_v3_coverage.json").read_text())
    names = {row["identity"]["name"] for row in artifact["drugs"]}
    assert len(artifact["drugs"]) >= 5
    assert {"sunvozertinib", "osimertinib", "midazolam", "warfarin", "metformin"}.issubset(names)


def test_evidence_router_assigns_one_primary_section_and_preserves_gap_reason():
    from backend.experimental_evidence_router import route_evidence
    activity = route_evidence({"endpoint": "IC50", "target": "EGFR", "reference_status": "REFERENCE_RESOLVED_SOURCE_RECORD"}, {"comparability_status": "RELATED_NOT_SAME_ENDPOINT", "reason": "Needs assay mapping"})
    pk = route_evidence({"endpoint": "Cmax", "reference_status": "REFERENCE_RESOLVED_REGULATORY"}, {"comparability_status": "UNSUPPORTED", "reason": "PK context is not a model endpoint"})
    assert activity["section"] == "ACTIVITY" and pk["section"] == "PK"
    assert activity["qualification_status"] == "QUALIFIED_RELATED"
    assert pk["qualification_status"] == "UNSUPPORTED" and pk["qualification_label"] == "PK context is not a model endpoint"
    assert activity["section"] in {"ACTIVITY", "ADMET", "METABOLISM", "PK", "TOXICITY", "UNCLASSIFIED"}
