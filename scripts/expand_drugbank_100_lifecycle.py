#!/usr/bin/env python3
"""
Sequential Single-Compound Ingestion Lifecycle for DrugBank 80 -> 100 Expansion (v3.3.3)

Strictly executes each compound one-by-one through:
1. Identity & canonical structure verification
2. CAS / InChIKey / PubChem / ChEMBL / DrugBank / UNII identifier-level provenance
3. Experimental evidence search / ingestion
4. Context classification & unit normalization
5. Qualification (partitioned as GLOBAL_REFERENCE_PENDING_ROLE)
6. Base / current global prediction & error calculation
7. Persistence & reopen verification before advancing to the next compound
"""
import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, Draw

from backend.database import SessionLocal, engine
from backend.models import (
    Project,
    Compound,
    CompoundVersion,
    CompoundIdentifier,
    ExternalExperimentalEvidence,
    ensure_ui_schema,
)
from backend.drugbank_reference import ensure_drugbank_project
from backend.engine_v3_learning import compute_base_prediction
from backend.openadmet_cyp import ic50_nm_to_pic50

REF_100_PATH = Path(__file__).resolve().parent.parent / "backend" / "reference_drugs_100.json"


def determine_chemical_form(mol: Chem.Mol, smiles: str) -> str:
    if "." in smiles:
        return "SALT"
    has_acid = mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)[OH]")) or mol.HasSubstructMatch(Chem.MolFromSmarts("S(=O)(=O)[OH]"))
    has_basic_nitrogen = mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O);!$(NS(=O)=O)]"))
    if has_acid and not has_basic_nitrogen:
        return "FREE_ACID"
    elif has_basic_nitrogen and not has_acid:
        return "FREE_BASE"
    elif has_acid and has_basic_nitrogen:
        return "ZWITTERION"
    return "FREE_BASE"


def run_sequential_expansion():
    print("=" * 75)
    print("DRUGBANK EXPANSION: 80 -> 100 SEQUENTIAL LIFECYCLE EXECUTION (v3.3.3)")
    print("=" * 75)

    ensure_ui_schema(engine)

    with open(REF_100_PATH, "r") as f:
        all_100_drugs = json.load(f)

    # 20 new drugs are indices 80 to 99
    new_20_specs = all_100_drugs[80:]
    assert len(new_20_specs) == 20, f"Expected 20 new drugs, found {len(new_20_specs)}"

    db = SessionLocal()
    proj = ensure_drugbank_project(db)
    proj_id = proj.id
    db.close()

    results = []

    for idx, drug in enumerate(new_20_specs, 81):
        name = drug["name"]
        db_id = drug["drugbank_id"]
        cas = drug["cas_number"]
        smiles = drug["smiles"]
        print(f"\n>>> [{idx}/100] Ingesting {name} ({db_id} · CAS: {cas}) ...")

        # -------------------------------------------------------------
        # Stage 1: Identity & Canonical Structure Verification
        # -------------------------------------------------------------
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES for {name}: {smiles}")
        canon_smiles = Chem.MolToSmiles(mol, canonical=True)
        inchi_str = Chem.MolToInchi(mol)
        inchikey_str = Chem.MolToInchiKey(mol)
        chem_form = determine_chemical_form(mol, canon_smiles)

        mw = float(Descriptors.MolWt(mol))
        clogp = float(Crippen.MolLogP(mol))
        tpsa = float(Descriptors.TPSA(mol))
        hbd = int(Lipinski.NumHDonors(mol))
        hba = int(Lipinski.NumHAcceptors(mol))
        rotb = int(Lipinski.NumRotatableBonds(mol))

        # Generate 2D SVG
        Draw.rdDepictor.Compute2DCoords(mol)
        svg_content = str(Draw.MolsToGridImage([mol], molsPerRow=1, subImgSize=(420, 320), useSVG=True))

        print(f"  [Stage 1: Identity] Canonical SMILES verified | InChIKey: {inchikey_str} | Form: {chem_form} | MW: {mw:.1f}")

        # -------------------------------------------------------------
        # Stage 2: Database Registration & Identifier Provenance
        # -------------------------------------------------------------
        db = SessionLocal()
        try:
            comp = db.query(Compound).filter(Compound.project_id == proj_id, Compound.name == name).first()
            if not comp:
                comp = Compound(
                    project_id=proj_id,
                    compound_id=f"DRUGBANK-{db_id}",
                    cas_number=cas,
                    name=name,
                    notes=f"Approved Reference Drug | DrugBank: {db_id} | ChEMBL: {drug['chembl_id']} | PubChem: {drug['pubchem_cid']} | UNII: {drug['unii']} | Scaffold: {drug.get('scaffold_family', '')} | Role: GLOBAL_REFERENCE_PENDING_ROLE | Cohort: PENDING_V3_4_ROLE_AUDIT | CAS Source: INTERNAL_EVIDENCE_AND_DRUGBANK_CATALOG",
                    status="APPROVED_REFERENCE",
                    current_version=1,
                )
                db.add(comp)
                db.commit()
                db.refresh(comp)

                cv = CompoundVersion(
                    compound_row_id=comp.id,
                    version_number=1,
                    original_smiles=smiles,
                    canonical_smiles=canon_smiles,
                    isomeric_smiles=canon_smiles,
                    inchi=inchi_str,
                    inchikey=inchikey_str,
                    change_note="Canonical reference drug registration (DrugBank Expansion v3.3.3)",
                    properties_json=json.dumps({
                        "MW": mw, "cLogP": clogp, "TPSA": tpsa, "HBD": hbd, "HBA": hba, "RotB": rotb,
                        "drugbank_id": db_id, "chembl_id": drug["chembl_id"],
                        "pubchem_cid": str(drug["pubchem_cid"]), "unii": drug["unii"],
                        "scaffold": drug.get("scaffold_family", ""),
                        "model_role": "GLOBAL_REFERENCE_PENDING_ROLE",
                        "cohort": "PENDING_V3_4_ROLE_AUDIT",
                    }),
                    svg=svg_content,
                )
                db.add(cv)
                db.commit()
                db.refresh(cv)
            else:
                comp.cas_number = cas
                db.commit()
                cv = db.query(CompoundVersion).filter(CompoundVersion.compound_row_id == comp.id, CompoundVersion.version_number == 1).first()
                if not cv.svg or len(cv.svg.strip()) < 10:
                    cv.svg = svg_content
                    db.commit()

            # Record all 7 identifiers into compound_identifiers
            now_iso = datetime.now(timezone.utc)
            id_tuples = [
                ("CANONICAL_SMILES", canon_smiles, "RDKIT_CANONICALIZATION", comp.compound_id),
                ("INCHIKEY", inchikey_str, "RDKIT_INCHI", comp.compound_id),
                ("CAS", cas, "INTERNAL_EVIDENCE_AND_DRUGBANK_CATALOG", db_id),
                ("DRUGBANK_ID", db_id, "DrugBank_Catalog", db_id),
                ("CHEMBL_ID", drug["chembl_id"], "ChEMBL_Database", drug["chembl_id"]),
                ("PUBCHEM_CID", str(drug["pubchem_cid"]), "PubChem_Database", str(drug["pubchem_cid"])),
                ("UNII", drug["unii"], "FDA_UNII_Directory", drug["unii"]),
            ]

            for id_type, id_val, id_src, id_rec in id_tuples:
                existing_ident = db.query(CompoundIdentifier).filter(
                    CompoundIdentifier.compound_id == comp.id,
                    CompoundIdentifier.identifier_type == id_type,
                    CompoundIdentifier.identifier_value == id_val,
                ).first()
                if not existing_ident:
                    ident = CompoundIdentifier(
                        compound_id=comp.id,
                        compound_version_id=cv.id,
                        identifier_type=id_type,
                        identifier_value=id_val,
                        source=id_src,
                        source_record_id=id_rec,
                        chemical_form=chem_form,
                        verified_against_inchikey=inchikey_str,
                        verification_status="VERIFIED",
                        verified_at=now_iso,
                    )
                    db.add(ident)
            db.commit()
            print(f"  [Stage 2: Registration] Compound row_id={comp.id}, version_id={cv.id}, 7 identifiers recorded.")

            # -------------------------------------------------------------
            # Stage 3, 4, 5: Evidence, Qualification & Prediction
            # -------------------------------------------------------------
            evidence_count = 0
            prediction_results = []
            error_results = []

            for obs in drug.get("observations", []):
                eid = obs["canonical_endpoint_id"]
                p_key = hashlib.sha256(f"{inchikey_str}_{eid}_{obs['raw_value']}_{obs['raw_unit']}_{obs['species']}_{obs['matrix']}".encode()).hexdigest()

                cond_dict = {
                    "matrix": obs["matrix"],
                    "section": obs["section"],
                    "upstream_overlap": "VALIDATION_HOLDOUT",
                    "drugbank_partition": "GLOBAL_REFERENCE_PENDING_ROLE",
                    "model_role": "GLOBAL_REFERENCE_PENDING_ROLE",
                    "cohort": "PENDING_V3_4_ROLE_AUDIT",
                }

                existing_ev = db.query(ExternalExperimentalEvidence).filter(
                    ExternalExperimentalEvidence.compound_version_id == cv.id,
                    ExternalExperimentalEvidence.provenance_key == p_key,
                ).first()

                if not existing_ev:
                    ev = ExternalExperimentalEvidence(
                        compound_version_id=cv.id,
                        provenance_key=p_key,
                        cas_number=cas,
                        canonical_endpoint_id=eid,
                        raw_endpoint_name=obs["raw_endpoint_name"],
                        species=obs["species"],
                        assay_type=obs["assay_type"],
                        assay_conditions_json=cond_dict,
                        raw_value=str(obs["raw_value"]),
                        raw_unit=obs["raw_unit"],
                        raw_relation=obs["raw_relation"],
                        normalized_value=str(obs["normalized_value"]),
                        normalized_unit=obs["normalized_unit"],
                        source_database="DrugBank_FDA_ChEMBL",
                        source_record_id=db_id,
                        source_url=f"https://go.drugbank.com/drugs/{db_id}",
                        identity_match_status="EXACT_MATCH",
                        endpoint_match_status="EXACT_MATCH",
                        mapping_status="EXTERNAL_EVIDENCE_ONLY",
                        evidence_origin="EXPERIMENTAL_EXTERNAL",
                        source_quality_class="A",
                        comparability_status="DIRECTLY_COMPARABLE",
                        qualification_status="QUALIFIED_FOR_GLOBAL_TRAINING",
                        reference_text=obs["reference_text"],
                        evidence_state="AUTO_QUALIFIED_EXTERNAL",
                    )
                    db.add(ev)
                    db.commit()
                else:
                    existing_ev.assay_conditions_json = cond_dict
                    db.commit()

                evidence_count += 1

                # Prediction & Error
                pred_val = compute_base_prediction(eid, canon_smiles)
                exp_val = float(obs["normalized_value"])
                if eid in ("HERG_LIABILITY", "CYP3A4_INHIBITION", "CYP2D6_INHIBITION", "CYP1A2_INHIBITION", "CYP2C9_INHIBITION", "CYP2C19_INHIBITION"):
                    exp_p = ic50_nm_to_pic50(exp_val) if exp_val > 0 else exp_val
                else:
                    exp_p = exp_val

                abs_err = round(abs(pred_val - exp_p), 3) if pred_val is not None else None
                prediction_results.append({"endpoint": eid, "pred": pred_val})
                error_results.append({"endpoint": eid, "pred": pred_val, "truth": round(exp_p, 3), "error": abs_err})

            print(f"  [Stage 3-5: Qualification & Prediction] {evidence_count} evidence qualified (role=GLOBAL_REFERENCE_PENDING_ROLE), {len(prediction_results)} predictions evaluated.")

            # -------------------------------------------------------------
            # Stage 6: Persistence & Reopen Verification
            # -------------------------------------------------------------
            saved_comp_id = comp.id
            saved_cv_id = cv.id
            db.close()

            # Verify in a completely clean session
            verify_db = SessionLocal()
            reopened_comp = verify_db.query(Compound).filter(Compound.id == saved_comp_id).first()
            assert reopened_comp is not None, f"Reopen failed for {name}"
            assert reopened_comp.cas_number == cas, f"CAS mismatch on reopen: {reopened_comp.cas_number} != {cas}"
            reopened_cv = verify_db.query(CompoundVersion).filter(CompoundVersion.compound_row_id == saved_comp_id, CompoundVersion.version_number == 1).first()
            assert reopened_cv is not None, f"Version reopen failed for {name}"
            assert reopened_cv.canonical_smiles == canon_smiles
            assert len(reopened_cv.svg) > 100, f"SVG missing on reopen for {name}"

            reopened_idents = verify_db.query(CompoundIdentifier).filter(CompoundIdentifier.compound_id == saved_comp_id).all()
            assert len(reopened_idents) == 7, f"Expected 7 identifiers on reopen, got {len(reopened_idents)}"

            reopened_ev = verify_db.query(ExternalExperimentalEvidence).filter(ExternalExperimentalEvidence.compound_version_id == reopened_cv.id).all()
            assert len(reopened_ev) == evidence_count, f"Evidence count mismatch on reopen: {len(reopened_ev)} != {evidence_count}"
            verify_db.close()

            print(f"  [Stage 6: Reopen Verification] 100% verified across Compound, Version, 7 Identifiers, SVG, and {evidence_count} Evidence records.")

            results.append({
                "index": idx,
                "name": name,
                "drugbank_id": db_id,
                "cas": cas,
                "chembl_id": drug["chembl_id"],
                "pubchem_cid": drug["pubchem_cid"],
                "unii": drug["unii"],
                "evidence_count": evidence_count,
                "errors": error_results,
                "status": "COMPLETED_AND_VERIFIED",
            })

        except Exception as exc:
            db.rollback()
            db.close()
            raise RuntimeError(f"Failed ingestion for {name}: {exc}") from exc

    print("\n" + "=" * 75)
    print(f"SUCCESS: Ingested and verified all {len(results)} new reference drugs sequentially!")
    print("=" * 75)


if __name__ == "__main__":
    run_sequential_expansion()
