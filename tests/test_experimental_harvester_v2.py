from backend import experimental_harvester as h


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
