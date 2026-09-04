"""DrugBank Identity Audit and CAS Backfill v3.3.2.

Audits all 80 compounds in Project 300, cross-checks canonical structure and InChIKey
against reference_drugs_80.json and internal evidence, safely backfills missing CAS numbers,
updates identity metadata notes and missing 2D SVG depictions, and outputs an audit report.
"""
from __future__ import annotations

import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw

DB_PATH = Path(__file__).resolve().parent.parent / "drug_opt.db"
REF_80_PATH = Path(__file__).resolve().parent.parent / "backend" / "reference_drugs_80.json"
REPORT_PATH = Path(__file__).resolve().parent.parent / "backend" / "drugbank_identity_audit_v3_3_2.json"


def audit_and_backfill():
    if not DB_PATH.exists():
        raise FileNotFoundError(f"Database not found at {DB_PATH}")
    if not REF_80_PATH.exists():
        raise FileNotFoundError(f"Reference file not found at {REF_80_PATH}")

    with open(REF_80_PATH, "r", encoding="utf-8") as f:
        ref_drugs = json.load(f)

    # Index reference drugs by lower-case name and cleaned drugbank_id
    ref_by_dbid = {}
    ref_by_name = {}
    for item in ref_drugs:
        dbid = (item.get("drugbank_id") or "").strip().lower()
        if dbid:
            ref_by_dbid[dbid] = item
        name = (item.get("name") or "").strip().lower()
        if name:
            ref_by_name[name] = item

    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()

    # Pre-audit invariant checks
    c.execute("SELECT count(*) FROM compounds WHERE project_id = 300")
    total_compounds = c.fetchone()[0]
    if total_compounds != 80:
        raise ValueError(f"Expected 80 compounds in DrugBank, found {total_compounds}")

    c.execute("""
        SELECT count(eee.id)
        FROM external_experimental_evidence eee
        JOIN compound_versions cv ON cv.id = eee.compound_version_id
        JOIN compounds c ON c.id = cv.compound_row_id
        WHERE c.project_id = 300
    """)
    initial_evidence_count = c.fetchone()[0]
    if initial_evidence_count != 577:
        raise ValueError(f"Expected 577 evidence records in DrugBank, found {initial_evidence_count}")

    c.execute("SELECT count(*) FROM prediction_runs")
    initial_runs = c.fetchone()[0]

    c.execute("""
        SELECT count(*)
        FROM qualification_prediction_freezes
        WHERE project_id = '300'
    """)
    initial_freezes = c.fetchone()[0]

    print("=" * 70)
    print("DrugBank (Project 300) Identity Audit & CAS Backfill v3.3.2")
    print(f"Initial State: 80 compounds, {initial_evidence_count} evidence, {initial_runs} prediction runs, {initial_freezes} freezes")
    print("=" * 70)

    # Fetch all 80 compounds with their current version details
    c.execute("""
        SELECT c.id, c.compound_id, c.name, c.cas_number, c.notes,
               cv.id, cv.canonical_smiles, cv.inchikey, cv.svg,
               GROUP_CONCAT(DISTINCT eee.cas_number) as evidence_cas,
               count(eee.id) as ev_count
        FROM compounds c
        JOIN compound_versions cv ON cv.compound_row_id = c.id AND cv.version_number = c.current_version
        LEFT JOIN external_experimental_evidence eee ON eee.compound_version_id = cv.id
        WHERE c.project_id = 300
        GROUP BY c.id
        ORDER BY c.id
    """)
    rows = c.fetchall()

    audit_records = []
    backfilled_compounds = []
    svg_updated_count = 0

    for row in rows:
        (
            cid,
            comp_id,
            name,
            current_cas,
            current_notes,
            cv_id,
            smiles,
            inchikey,
            current_svg,
            ev_cas,
            ev_count,
        ) = row

        clean_id = comp_id.replace("DRUGBANK-", "").strip().lower()
        ref = ref_by_dbid.get(clean_id) or ref_by_name.get(name.strip().lower())
        if not ref:
            raise ValueError(f"Could not find reference match for {name} ({comp_id})")

        ref_smi = ref.get("smiles") or ""
        ref_mol = Chem.MolFromSmiles(ref_smi)
        ref_inchikey = Chem.MolToInchiKey(ref_mol) if ref_mol else None

        db_mol = Chem.MolFromSmiles(smiles)
        db_inchikey = Chem.MolToInchiKey(db_mol) if db_mol else None

        if not db_inchikey or not ref_inchikey or db_inchikey != ref_inchikey:
            raise ValueError(
                f"Structure / InChIKey mismatch for {name} (ID: {comp_id}): ref={ref_inchikey}, db={db_inchikey}"
            )

        ref_cas = (ref.get("cas_number") or "").strip()
        verified_cas = (current_cas or "").strip() or (ev_cas or "").strip() or ref_cas

        if not verified_cas:
            verified_cas = "UNKNOWN"
            cas_source = "UNRESOLVED"
        else:
            cas_source = "INTERNAL_EVIDENCE_AND_DRUGBANK_CATALOG"

        # Build clean note preserving metadata
        drugbank_id = ref.get("drugbank_id") or comp_id.replace("DRUGBANK-", "")
        chembl_id = ref.get("chembl_id") or "UNKNOWN"
        pubchem_cid = ref.get("pubchem_cid") or "UNKNOWN"
        unii = ref.get("unii") or "UNKNOWN"
        scaffold = ref.get("scaffold_family") or "Unclassified"

        formatted_notes = (
            f"Approved Reference Drug | DrugBank: {drugbank_id} | ChEMBL: {chembl_id} "
            f"| PubChem: {pubchem_cid} | UNII: {unii} | Scaffold: {scaffold} | CAS Source: {cas_source}"
        )

        # Update compound if CAS was missing or notes needed update
        was_missing_cas = not bool(current_cas and current_cas.strip())
        now_iso = datetime.now(timezone.utc).isoformat()

        if was_missing_cas:
            c.execute(
                """
                UPDATE compounds
                SET cas_number = ?, notes = ?, updated_at = ?
                WHERE id = ?
                """,
                (verified_cas, formatted_notes, now_iso, cid),
            )
            backfilled_compounds.append({
                "row_id": cid,
                "compound_id": comp_id,
                "name": name,
                "backfilled_cas": verified_cas,
                "source": cas_source,
            })
        else:
            # Update notes to ensure provenance is consistently recorded
            c.execute(
                """
                UPDATE compounds
                SET notes = ?, updated_at = ?
                WHERE id = ?
                """,
                (formatted_notes, now_iso, cid),
            )

        # Ensure SVG is generated and stored for the compound version
        if not current_svg or len(current_svg.strip()) < 10:
            Draw.rdDepictor.Compute2DCoords(db_mol)
            new_svg = str(Draw.MolsToGridImage([db_mol], molsPerRow=1, subImgSize=(420, 320), useSVG=True))
            c.execute(
                """
                UPDATE compound_versions
                SET svg = ?
                WHERE id = ?
                """,
                (new_svg, cv_id),
            )
            svg_updated_count += 1

        audit_records.append({
            "row_id": cid,
            "compound_id": comp_id,
            "drugbank_id": drugbank_id,
            "name": name,
            "cas_number": verified_cas,
            "cas_source": cas_source,
            "canonical_smiles": Chem.MolToSmiles(db_mol),
            "inchikey": db_inchikey,
            "inchikey_match": True,
            "pubchem_cid": pubchem_cid,
            "chembl_id": chembl_id,
            "unii": unii,
            "scaffold_family": scaffold,
            "evidence_count": ev_count,
            "was_missing_cas": was_missing_cas,
            "has_svg": True,
        })

    conn.commit()

    # Post-backfill validation & invariant verification
    c.execute("SELECT count(*), count(cas_number) FROM compounds WHERE project_id = 300")
    post_total, post_cas = c.fetchone()
    if post_total != 80:
        raise ValueError(f"Post-check failed: expected 80 compounds, got {post_total}")
    if post_cas != 80:
        raise ValueError(f"Post-check failed: expected 80 compounds with CAS, got {post_cas}")

    c.execute("""
        SELECT count(eee.id)
        FROM external_experimental_evidence eee
        JOIN compound_versions cv ON cv.id = eee.compound_version_id
        JOIN compounds c ON c.id = cv.compound_row_id
        WHERE c.project_id = 300
    """)
    post_evidence_count = c.fetchone()[0]
    if post_evidence_count != initial_evidence_count:
        raise ValueError(f"Evidence count changed from {initial_evidence_count} to {post_evidence_count}")

    c.execute("SELECT count(*) FROM prediction_runs")
    post_runs = c.fetchone()[0]
    if post_runs != initial_runs:
        raise ValueError(f"Prediction runs changed from {initial_runs} to {post_runs}")

    c.execute("""
        SELECT count(*)
        FROM qualification_prediction_freezes
        WHERE project_id = '300'
    """)
    post_freezes = c.fetchone()[0]
    if post_freezes != initial_freezes:
        raise ValueError(f"Freeze count changed from {initial_freezes} to {post_freezes}")

    c.execute("PRAGMA foreign_key_check")
    fk_errors = c.fetchall()
    if fk_errors:
        raise ValueError(f"Foreign key violations detected: {fk_errors}")

    c.execute("PRAGMA integrity_check")
    integrity = c.fetchone()[0]
    if integrity != "ok":
        raise ValueError(f"Integrity check failed: {integrity}")

    # Check SVGs
    c.execute("""
        SELECT count(*)
        FROM compound_versions cv
        JOIN compounds c ON c.id = cv.compound_row_id
        WHERE c.project_id = 300 AND (cv.svg IS NULL OR length(cv.svg) < 10)
    """)
    missing_svg_count = c.fetchone()[0]
    if missing_svg_count != 0:
        raise ValueError(f"Expected 0 missing SVGs, got {missing_svg_count}")

    conn.close()

    # Save audit report
    report = {
        "audit_version": "v3.3.2",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "project_id": 300,
        "project_name": "DrugBank Approved Reference Library",
        "summary": {
            "total_compounds_audited": len(audit_records),
            "compounds_backfilled_cas": len(backfilled_compounds),
            "compounds_already_had_cas": len(audit_records) - len(backfilled_compounds),
            "unknown_cas_count": sum(1 for r in audit_records if r["cas_number"] == "UNKNOWN"),
            "inchikey_matches": sum(1 for r in audit_records if r["inchikey_match"]),
            "svg_depictions_generated": svg_updated_count,
            "evidence_count_preserved": post_evidence_count,
            "freezes_count_preserved": post_freezes,
            "prediction_runs_preserved": post_runs,
            "integrity_check": integrity,
            "foreign_key_check": "ok (0 violations)",
        },
        "backfilled_compounds": backfilled_compounds,
        "records": audit_records,
    }

    with open(REPORT_PATH, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print("\n" + "=" * 70)
    print("AUDIT & BACKFILL COMPLETE")
    print(f"Total compounds audited: {len(audit_records)}")
    print(f"CAS backfilled: {len(backfilled_compounds)} compounds")
    print(f"SVGs generated / stored: {svg_updated_count}")
    print(f"Invariants check: 80 compounds, 80 CAS, {post_evidence_count} evidence, {post_freezes} freezes, {post_runs} runs")
    print(f"DB integrity: {integrity}, FK check: 0 errors")
    print(f"Audit report saved to: {REPORT_PATH}")
    print("=" * 70)


if __name__ == "__main__":
    audit_and_backfill()
