import sys
import json
import urllib.request
from sqlalchemy import create_engine, select, func
from backend.models import Project, Compound, CompoundVersion, ExternalExperimentalEvidence, PredictionRun, CompoundIdentifier

BASE_URL = "http://127.0.0.1:8765"

def get_json(url):
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read().decode("utf-8"))

def main():
    engine = create_engine("sqlite:///drug_opt.db")
    with engine.connect() as conn:
        projects = conn.execute(select(Project.id, Project.name)).fetchall()
        proj_map = {p[0]: p[1] for p in projects}
        print(f"Projects in DB: {proj_map}")
        assert set(proj_map.keys()) == {1, 3, 5, 300}, f"Unexpected projects: {proj_map.keys()}"

        # Compound counts
        for pid, expected in [(1, 4), (3, 7), (5, 4), (300, 100)]:
            count = conn.execute(select(func.count(Compound.id)).where(Compound.project_id == pid)).scalar()
            print(f"Project {pid} ({proj_map[pid]}): {count} compounds (expected {expected})")
            assert count == expected, f"Project {pid} count mismatch: {count} != {expected}"

        # DrugBank audit
        db_comps = conn.execute(select(Compound.id, Compound.name, Compound.cas_number).where(Compound.project_id == 300)).fetchall()
        assert len(db_comps) == 100
        for row in db_comps:
            cid, name, cas = row
            assert cas and cas != "UNKNOWN", f"Missing CAS for {name} ({cid})"

            # Check identifier rows
            id_rows = conn.execute(select(CompoundIdentifier.identifier_type, CompoundIdentifier.identifier_value, CompoundIdentifier.verification_status).where(CompoundIdentifier.compound_id == cid)).fetchall()
            id_dict = {r[0]: r[1] for r in id_rows}
            status_dict = {r[0]: r[2] for r in id_rows}
            assert len(id_rows) == 7, f"Compound {name} ({cid}) has {len(id_rows)} identifiers, expected 7"
            assert "CAS" in id_dict and id_dict["CAS"] == cas
            assert "CANONICAL_SMILES" in id_dict
            assert "INCHIKEY" in id_dict
            assert "PUBCHEM_CID" in id_dict
            assert "CHEMBL_ID" in id_dict
            assert "DRUGBANK_ID" in id_dict
            assert "UNII" in id_dict
            for itype, st in status_dict.items():
                assert st == "VERIFIED", f"Identifier {itype} not VERIFIED for {name}: {st}"

        # Evidence count in DrugBank
        db_ev_count = conn.execute(select(func.count(ExternalExperimentalEvidence.id)).join(CompoundVersion, CompoundVersion.id == ExternalExperimentalEvidence.compound_version_id).join(Compound, Compound.id == CompoundVersion.compound_row_id).where(Compound.project_id == 300)).scalar()
        print(f"DrugBank qualified evidence rows: {db_ev_count}")
        assert db_ev_count >= 650, f"Evidence count too low: {db_ev_count}"

    # Test API endpoints for all projects
    print("Verifying API parity for all projects...")
    for pid in [1, 3, 5, 300]:
        pdata = get_json(f"{BASE_URL}/api/projects/{pid}")
        comps = pdata.get("compounds", [])
        print(f"API Project {pid} has {len(comps)} compounds")
        for c in comps:
            row_id = c["row_id"]
            c_detail = get_json(f"{BASE_URL}/api/compounds/{row_id}")
            assert c_detail["name"] == c["name"]
            assert c_detail["cas_number"] == c.get("cas_number")
            version = c_detail.get("version")
            if version:
                assert version.get("svg") and "<svg" in version["svg"], f"Missing SVG in compound {row_id}"
                assert version.get("canonical_smiles"), f"Missing SMILES in compound {row_id}"
                assert version.get("inchikey"), f"Missing InChIKey in compound {row_id}"
            
            if pid == 300:
                assert len(c_detail.get("identifiers", [])) == 7, f"API identifiers count mismatch in compound {row_id}"
                assert c_detail.get("verification_status") == "VERIFIED"
                assert c_detail.get("evidence_count", 0) > 0, f"Expected evidence for DrugBank compound {row_id}"
                assert c.get("cas_number"), f"Missing CAS in list view for {c['name']}"
                assert c.get("drugbank_id"), f"Missing DrugBank ID in list view for {c['name']}"
                assert c.get("chembl_id"), f"Missing ChEMBL ID in list view for {c['name']}"
                assert c.get("pubchem_cid"), f"Missing PubChem CID in list view for {c['name']}"
                assert c.get("unii"), f"Missing UNII in list view for {c['name']}"

    print("ALL AUDIT AND VERIFICATION CHECKS PASSED PERFECTLY!")

if __name__ == "__main__":
    main()
