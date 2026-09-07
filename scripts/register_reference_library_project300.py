"""Safely register the file reference library in runtime Project 300.

Only identity/structure rows are inserted.  No prediction calculation is
run, so historical PredictionRuns and experimental evidence are untouched.
Run with --apply for the explicitly requested additive registration.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from backend.database import SessionLocal  # noqa: E402
from backend.models import Compound, CompoundIdentifier, CompoundVersion, Project  # noqa: E402


def identity(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("invalid SMILES")
    return {
        "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
        "isomeric_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
        "inchi": Chem.MolToInchi(mol),
        "inchikey": Chem.MolToInchiKey(mol),
    }


def manifest(session, project):
    rows = []
    for compound in sorted(project.compounds, key=lambda x: x.id):
        version = sorted(compound.versions, key=lambda x: x.version_number)[-1] if compound.versions else None
        rows.append({
            "id": compound.id, "compound_id": compound.compound_id, "name": compound.name,
            "cas_number": compound.cas_number, "current_version": compound.current_version,
            "inchikey": version.inchikey if version else None,
            "canonical_smiles": version.canonical_smiles if version else None,
        })
    return {"artifact": "PROJECT300_IDENTITY_MANIFEST_BEFORE_V70", "created_at": datetime.now(timezone.utc).isoformat(), "project_id": project.id, "compound_n": len(rows), "compounds": rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    library = json.loads((ROOT / "validation/reference_library_v1_1000.json").read_text())["compounds"]
    session = SessionLocal()
    try:
        project = session.get(Project, 300)
        if project is None:
            raise RuntimeError("Project 300 not found")
        before = manifest(session, project)
        backup = ROOT / "validation/project300_identity_manifest_before_v70.json"
        if not backup.exists():
            backup.write_text(json.dumps(before, indent=2, sort_keys=True) + "\n")
        existing = {v.inchikey for c in project.compounds for v in c.versions if v.inchikey}
        plan = [x for x in library if x["inchikey"] not in existing]
        if not args.apply:
            print(json.dumps({"mode": "DRY_RUN", "before": len(project.compounds), "to_add": len(plan), "after": len(project.compounds) + len(plan)}, indent=2))
            return
        added = 0
        for item in plan:
            ident = identity(item["smiles"])
            if ident["inchikey"] in existing:
                continue
            label = str(item.get("library_id") or ("REFERENCE:" + ident["inchikey"]))[:50]
            # Compound labels are source-backed library keys, not invented
            # DrugBank IDs.  The source family is preserved in notes/identifiers.
            compound = Compound(project_id=300, compound_id=label,
                                name=item.get("name") or ("Reference compound " + ident["inchikey"]),
                                cas_number=item.get("cas"),
                                notes=json.dumps({"source_name": item.get("source_family"), "source_record_id": item.get("source_record"), "library_id": item.get("library_id"), "role": item.get("role")}, sort_keys=True),
                                status="STRUCTURE_READY", current_version=1)
            session.add(compound)
            session.flush()
            version = CompoundVersion(compound_row_id=compound.id, version_number=1,
                original_smiles=item["smiles"], canonical_smiles=ident["canonical_smiles"],
                isomeric_smiles=ident["isomeric_smiles"], inchi=ident["inchi"], inchikey=ident["inchikey"],
                change_note="Reference Library V1 additive identity registration; no prediction run",
                properties_json=item.get("descriptor") or {}, alerts_json=None,
                assessment_json=None, calculation_json=None, svg="", highlighted_svg="")
            session.add(version)
            session.flush()
            for kind, value in (("INCHIKEY", ident["inchikey"]), ("SOURCE_RECORD", str(item.get("source_record") or item.get("library_id") or ident["inchikey"]))):
                session.add(CompoundIdentifier(compound_id=compound.id, compound_version_id=version.id,
                    identifier_type=kind, identifier_value=value[:500], source=str(item.get("source_family") or "REFERENCE_LIBRARY")[:100],
                    source_record_id=str(item.get("source_record") or item.get("library_id") or "")[:100],
                    chemical_form="UNRESOLVED_SOURCE_FORM" if item.get("role") == "DIVERSITY_SELECTED_EXTERNAL" else "REFERENCE_FORM",
                    verified_against_inchikey=ident["inchikey"], verification_status="VERIFIED"))
            existing.add(ident["inchikey"])
            added += 1
        session.commit()
        final_count = session.query(Compound).filter(Compound.project_id == 300).count()
        final_keys = [v.inchikey for c in session.query(Compound).filter(Compound.project_id == 300).all() for v in c.versions if v.inchikey]
        if final_count != 1000 or len(final_keys) != len(set(final_keys)):
            raise RuntimeError(f"registration invariant failed: compounds={final_count} keys={len(final_keys)} unique={len(set(final_keys))}")
        print(json.dumps({"mode": "APPLIED", "before": len(before["compounds"]), "added": added, "after": final_count, "unique_inchikeys": len(set(final_keys))}, indent=2))
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()


if __name__ == "__main__":
    main()
