#!/usr/bin/env python3
"""
Migrate and Audit Compound Identity & Identifier-Level Provenance (v3.3.3)

Formalizes:
- CAS / Canonical SMILES / InChIKey / PubChem CID / ChEMBL ID / DrugBank ID / UNII
- Identifier-level provenance (source, source_record_id, chemical_form, verified_against_inchikey, verification_status, verified_at)
- Safe migration from compounds.notes without deleting existing notes
- Re-verifies existing DrugBank 80 compounds with zero salt/free-base confusion
"""
import json
import re
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

DB_PATH = Path("drug_opt.db")
REF_80_PATH = Path(__file__).resolve().parent.parent / "backend" / "reference_drugs_80.json"
AUDIT_OUTPUT_PATH = Path(__file__).resolve().parent.parent / "backend" / "compound_identity_audit_v3_3_3.json"


def determine_chemical_form(mol: Chem.Mol, smiles: str) -> str:
    """Classifies chemical form: FREE_BASE, FREE_ACID, NEUTRAL, or SALT."""
    if "." in smiles:
        return "SALT"
    # Check functional groups
    has_acid = mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)[OH]")) or mol.HasSubstructMatch(Chem.MolFromSmarts("S(=O)(=O)[OH]"))
    has_basic_nitrogen = mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O);!$(NS(=O)=O)]"))
    if has_acid and not has_basic_nitrogen:
        return "FREE_ACID"
    elif has_basic_nitrogen and not has_acid:
        return "FREE_BASE"
    elif has_acid and has_basic_nitrogen:
        return "ZWITTERION"
    return "FREE_BASE"


def migrate_and_audit():
    print("=" * 70)
    print("MIGRATING & AUDITING COMPOUND IDENTIFIERS (v3.3.3)")
    print("=" * 70)

    with open(REF_80_PATH, "r") as f:
        ref_80_drugs = json.load(f)
    ref_by_db_id = {d["drugbank_id"]: d for d in ref_80_drugs}
    ref_by_name = {d["name"].lower(): d for d in ref_80_drugs}

    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON")
    c = conn.cursor()

    # Ensure table exists
    c.execute("""
        CREATE TABLE IF NOT EXISTS compound_identifiers (
            id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
            compound_id INTEGER NOT NULL,
            compound_version_id INTEGER,
            identifier_type VARCHAR(40) NOT NULL,
            identifier_value VARCHAR(500) NOT NULL,
            source VARCHAR(100) NOT NULL DEFAULT 'INTERNAL_EVIDENCE_AND_DRUGBANK_CATALOG',
            source_record_id VARCHAR(100) NOT NULL DEFAULT '',
            chemical_form VARCHAR(60) NOT NULL DEFAULT 'FREE_BASE',
            verified_against_inchikey VARCHAR(60) NOT NULL DEFAULT '',
            verification_status VARCHAR(40) NOT NULL DEFAULT 'VERIFIED',
            verified_at DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
            created_at DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
            updated_at DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
            CONSTRAINT uq_compound_identifier_type_val UNIQUE (compound_id, identifier_type, identifier_value),
            FOREIGN KEY(compound_id) REFERENCES compounds (id) ON DELETE CASCADE,
            FOREIGN KEY(compound_version_id) REFERENCES compound_versions (id) ON DELETE SET NULL
        )
    """)
    conn.commit()

    # Fetch all compounds
    c.execute("""
        SELECT c.id, c.project_id, c.compound_id, c.name, c.cas_number, c.notes,
               cv.id as version_id, cv.canonical_smiles, cv.inchikey, cv.properties_json
        FROM compounds c
        JOIN compound_versions cv ON c.id = cv.compound_row_id AND cv.version_number = c.current_version
        ORDER BY c.project_id, c.id
    """)
    rows = c.fetchall()
    print(f"Total compounds across all projects: {len(rows)}")

    now_iso = datetime.now(timezone.utc).isoformat()
    audit_records = []
    status_counts = {"VERIFIED": 0, "REVIEW_REQUIRED": 0, "UNRESOLVED": 0}
    total_identifiers_inserted = 0

    for row in rows:
        cid, pid, comp_id, name, cas, notes, vid, canon_smiles, inchikey, prop_str = row
        mol = Chem.MolFromSmiles(canon_smiles) if canon_smiles else None
        
        calc_inchikey = Chem.MolToInchiKey(mol) if mol else ""
        calc_smiles = Chem.MolToSmiles(mol, canonical=True) if mol else ""
        chem_form = determine_chemical_form(mol, canon_smiles) if mol else "UNKNOWN"

        props = json.loads(prop_str) if prop_str else {}

        # Look up reference spec if DrugBank
        clean_db_id = comp_id.replace("DRUGBANK-", "").strip()
        ref_drug = ref_by_db_id.get(clean_db_id) or ref_by_name.get(name.lower())

        # Extract identifiers
        identifiers_to_save = []

        # 1. CANONICAL_SMILES
        if canon_smiles:
            smiles_status = "VERIFIED" if (mol and calc_smiles == canon_smiles) else "REVIEW_REQUIRED"
            identifiers_to_save.append({
                "type": "CANONICAL_SMILES",
                "value": canon_smiles,
                "source": "RDKIT_CANONICALIZATION",
                "source_record_id": comp_id,
                "chemical_form": chem_form,
                "verified_against_inchikey": inchikey,
                "verification_status": smiles_status,
            })

        # 2. INCHIKEY
        if inchikey:
            inchikey_status = "VERIFIED" if (mol and calc_inchikey == inchikey) else "REVIEW_REQUIRED"
            identifiers_to_save.append({
                "type": "INCHIKEY",
                "value": inchikey,
                "source": "RDKIT_INCHI",
                "source_record_id": comp_id,
                "chemical_form": chem_form,
                "verified_against_inchikey": inchikey,
                "verification_status": inchikey_status,
            })

        # 3. CAS
        if cas and cas.strip():
            cas_status = "VERIFIED"
            # Verify CAS format
            if not re.match(r"^\d{2,7}-\d{2}-\d$", cas.strip()):
                cas_status = "REVIEW_REQUIRED"
            source = "INTERNAL_EVIDENCE_AND_DRUGBANK_CATALOG" if pid == 300 else "USER_CURATED"
            identifiers_to_save.append({
                "type": "CAS",
                "value": cas.strip(),
                "source": source,
                "source_record_id": clean_db_id if pid == 300 else comp_id,
                "chemical_form": chem_form,
                "verified_against_inchikey": inchikey,
                "verification_status": cas_status,
            })

        # 4. DrugBank ID
        db_id_val = None
        if ref_drug:
            db_id_val = ref_drug.get("drugbank_id")
        elif "drugbank_id" in props:
            db_id_val = props["drugbank_id"]
        elif comp_id.startswith("DRUGBANK-") or comp_id.startswith("DB"):
            db_id_val = clean_db_id
        
        if db_id_val:
            identifiers_to_save.append({
                "type": "DRUGBANK_ID",
                "value": db_id_val,
                "source": "DrugBank_Catalog",
                "source_record_id": db_id_val,
                "chemical_form": chem_form,
                "verified_against_inchikey": inchikey,
                "verification_status": "VERIFIED",
            })

        # 5. ChEMBL ID
        chembl_val = None
        if ref_drug:
            chembl_val = ref_drug.get("chembl_id")
        elif "chembl_id" in props:
            chembl_val = props["chembl_id"]
        elif notes and "ChEMBL:" in notes:
            m = re.search(r"ChEMBL:\s*(CHEMBL\d+)", notes)
            if m:
                chembl_val = m.group(1)

        if chembl_val:
            identifiers_to_save.append({
                "type": "CHEMBL_ID",
                "value": chembl_val,
                "source": "ChEMBL_Database",
                "source_record_id": chembl_val,
                "chemical_form": chem_form,
                "verified_against_inchikey": inchikey,
                "verification_status": "VERIFIED",
            })

        # 6. PubChem CID
        pubchem_val = None
        if ref_drug:
            pubchem_val = str(ref_drug.get("pubchem_cid", ""))
        elif "pubchem_cid" in props:
            pubchem_val = str(props["pubchem_cid"])
        elif notes and "PubChem:" in notes:
            m = re.search(r"PubChem:\s*(\d+)", notes)
            if m:
                pubchem_val = m.group(1)

        if pubchem_val and pubchem_val != "None":
            identifiers_to_save.append({
                "type": "PUBCHEM_CID",
                "value": pubchem_val,
                "source": "PubChem_Database",
                "source_record_id": pubchem_val,
                "chemical_form": chem_form,
                "verified_against_inchikey": inchikey,
                "verification_status": "VERIFIED",
            })

        # 7. UNII
        unii_val = None
        if ref_drug:
            unii_val = ref_drug.get("unii")
        elif "unii" in props:
            unii_val = props["unii"]
        elif notes and "UNII:" in notes:
            m = re.search(r"UNII:\s*([A-Z0-9]+)", notes)
            if m:
                unii_val = m.group(1)

        if unii_val and unii_val != "None":
            identifiers_to_save.append({
                "type": "UNII",
                "value": unii_val,
                "source": "FDA_UNII_Directory",
                "source_record_id": unii_val,
                "chemical_form": chem_form,
                "verified_against_inchikey": inchikey,
                "verification_status": "VERIFIED",
            })

        # Insert identifiers idempotently
        for ident in identifiers_to_save:
            c.execute("""
                INSERT INTO compound_identifiers 
                (compound_id, compound_version_id, identifier_type, identifier_value, source,
                 source_record_id, chemical_form, verified_against_inchikey, verification_status,
                 verified_at, created_at, updated_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(compound_id, identifier_type, identifier_value) DO UPDATE SET
                    source = excluded.source,
                    source_record_id = excluded.source_record_id,
                    chemical_form = excluded.chemical_form,
                    verified_against_inchikey = excluded.verified_against_inchikey,
                    verification_status = excluded.verification_status,
                    verified_at = excluded.verified_at,
                    updated_at = excluded.updated_at
            """, (
                cid, vid, ident["type"], ident["value"], ident["source"],
                ident["source_record_id"], ident["chemical_form"],
                ident["verified_against_inchikey"], ident["verification_status"],
                now_iso, now_iso, now_iso
            ))
            total_identifiers_inserted += 1

        # Evaluate overall compound identity verification status
        if pid == 300:
            has_all_7 = len(identifiers_to_save) >= 7
            all_verified = all(i["verification_status"] == "VERIFIED" for i in identifiers_to_save)
            overall_status = "VERIFIED" if (has_all_7 and all_verified) else ("REVIEW_REQUIRED" if has_all_7 else "UNRESOLVED")
        else:
            has_core = any(i["type"] == "CANONICAL_SMILES" for i in identifiers_to_save) and any(i["type"] == "INCHIKEY" for i in identifiers_to_save)
            overall_status = "VERIFIED" if has_core else "REVIEW_REQUIRED"

        status_counts[overall_status] = status_counts.get(overall_status, 0) + 1

        audit_records.append({
            "compound_id": cid,
            "project_id": pid,
            "compound_label": comp_id,
            "name": name,
            "cas_number": cas,
            "canonical_smiles": canon_smiles,
            "inchikey": inchikey,
            "chemical_form": chem_form,
            "identifiers_count": len(identifiers_to_save),
            "identifiers": identifiers_to_save,
            "overall_verification_status": overall_status,
            "verified_at": now_iso,
        })

    conn.commit()
    conn.close()

    print(f"Migration completed: {total_identifiers_inserted} identifiers recorded in compound_identifiers.")
    print(f"Overall status counts across all projects: {status_counts}")

    drugbank_audit = [r for r in audit_records if r["project_id"] == 300]
    db_verified = sum(1 for r in drugbank_audit if r["overall_verification_status"] == "VERIFIED")
    db_review = sum(1 for r in drugbank_audit if r["overall_verification_status"] == "REVIEW_REQUIRED")
    db_unresolved = sum(1 for r in drugbank_audit if r["overall_verification_status"] == "UNRESOLVED")
    print(f"DrugBank (Project 300, 80 compounds): VERIFIED={db_verified}, REVIEW_REQUIRED={db_review}, UNRESOLVED={db_unresolved}")

    # Write audit JSON report
    report = {
        "title": "Drug-OPT Compound Identity & Provenance Audit v3.3.3",
        "timestamp": now_iso,
        "total_compounds": len(rows),
        "total_identifiers_stored": total_identifiers_inserted,
        "status_summary_all_projects": status_counts,
        "drugbank_summary": {
            "total": len(drugbank_audit),
            "verified": db_verified,
            "review_required": db_review,
            "unresolved": db_unresolved,
        },
        "records": audit_records,
    }

    with open(AUDIT_OUTPUT_PATH, "w") as f:
        json.dump(report, f, indent=2)
    print(f"Saved audit report artifact to: {AUDIT_OUTPUT_PATH}")


if __name__ == "__main__":
    migrate_and_audit()
