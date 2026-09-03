"""
DrugBank Reference Drug Library & Incremental Learning (Drug-OPT Stage 6 / v3.0.0 Global Completion).

Provides:
- Canonical 'DrugBank' project management (GLOBAL_MODEL_DEVELOPMENT mode)
- Step-by-step sequential reference drug ingestion across 40 distinct approved reference drugs:
    1-10: Gefitinib, Imatinib, Propranolol, Atorvastatin, Midazolam, Verapamil, Fluoxetine, Ketoconazole, Sildenafil, Quinidine
    11-15: Dextromethorphan, Amiodarone, Clarithromycin, Duloxetine, Haloperidol
    16-20: Paroxetine, Metoprolol, Terbinafine, Ritonavir, Cimetidine
    21-30 (Dev Training): Bupropion, Carvedilol, Clopidogrel, Diltiazem, Erythromycin, Flecainide, Lansoprazole, Nifedipine, Omeprazole, Simvastatin
    31-35 (Model Selection Validation): Celecoxib, Diazepam, Diclofenac, Indomethacin, Warfarin
    36-40 (Locked Final Test Cohort 2): Atenolol, Caffeine, Ibuprofen, Lorcaserin, Rosuvastatin
- 4-Tier Partitioning per endpoint:
    * DEVELOPMENT_TRAINING: Used for fitting candidate calibration models
    * MODEL_SELECTION_VALIDATION: Used for candidate selection & hyperparameter tuning
    * FINAL_TEST_COHORT_1_CONSUMED: Cimetidine (already evaluated; frozen from tuning)
    * LOCKED_FINAL_TEST_COHORT_2: Pristine locked external holdout; evaluated strictly once
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski
from sqlalchemy import select
from sqlalchemy.orm import Session

from backend.database import SessionLocal
from backend.models import (
    Project,
    Compound,
    CompoundVersion,
    ExternalExperimentalEvidence,
)
import backend.activity_models  # Ensure AssayDefinition is in metadata
from backend.openadmet_cyp import predict_chemeleon_cyp_pic50, ic50_nm_to_pic50, compute_fold_error
from backend.quantitative_safety_transporters import (
    predict_quantitative_herg_pic50,
    predict_quantitative_pgp_pic50,
    evaluate_safety_applicability_domain,
)

DRUGBANK_PROJECT_NAME = "DrugBank"
DRUGBANK_PROJECT_INDICATION = "Global Reference Drug Library (GLOBAL_MODEL_DEVELOPMENT)"
DRUGBANK_PROJECT_DESC = "Canonical reference drug library curated for Drug-OPT Global Prediction Engine v3.0 training and multi-tiered benchmarking."

ROLE_DEVELOPMENT_TRAINING = "DEVELOPMENT_TRAINING"
ROLE_MODEL_SELECTION_VALIDATION = "MODEL_SELECTION_VALIDATION"
ROLE_FINAL_TEST_COHORT_1_CONSUMED = "FINAL_TEST_COHORT_1_CONSUMED"
ROLE_LOCKED_FINAL_TEST_COHORT_2 = "LOCKED_FINAL_TEST_COHORT_2"

# Load full 40 reference drugs catalog
CATALOG_PATH = Path(__file__).parent / "reference_drugs_40.json"
if CATALOG_PATH.exists():
    with open(CATALOG_PATH, "r") as f:
        REFERENCE_DRUGS_CATALOG = json.load(f)
else:
    REFERENCE_DRUGS_CATALOG = []


def ensure_drugbank_project(db: Session) -> Project:
    """Ensures the canonical DrugBank project exists in the database."""
    proj = db.scalar(select(Project).where(Project.name == DRUGBANK_PROJECT_NAME))
    if not proj:
        proj = Project(
            name=DRUGBANK_PROJECT_NAME,
            target="PAN_TARGET_REFERENCE",
            molecule_type="Small Molecule",
            indication=DRUGBANK_PROJECT_INDICATION,
            mechanism_modality="SMALL_MOLECULE_DRUG",
            description=DRUGBANK_PROJECT_DESC,
        )
        db.add(proj)
        db.commit()
        db.refresh(proj)
    return proj


def ingest_reference_drug_by_spec(db: Session, drug_spec: Dict[str, Any]) -> Dict[str, Any]:
    """
    Ingests a single approved reference drug into the DrugBank project following the strict
    Identity -> Search -> Qualification -> Base Prediction -> Error pipeline.
    """
    proj = ensure_drugbank_project(db)
    mol = Chem.MolFromSmiles(drug_spec["smiles"])
    if mol is None:
        raise ValueError(f"Invalid SMILES for {drug_spec['name']}")

    canon_smiles = Chem.MolToSmiles(mol, canonical=True)
    inchi_str = Chem.MolToInchi(mol)
    inchikey_str = Chem.MolToInchiKey(mol)

    comp = db.scalar(select(Compound).where(Compound.project_id == proj.id, Compound.name == drug_spec["name"]))
    if not comp:
        comp = Compound(
            project_id=proj.id,
            compound_id=f"DRUGBANK-{drug_spec['drugbank_id']}",
            cas_number=drug_spec["cas_number"],
            name=drug_spec["name"],
            notes=f"Approved Reference Drug | DrugBank: {drug_spec['drugbank_id']} | ChEMBL: {drug_spec['chembl_id']} | PubChem: {drug_spec['pubchem_cid']} | UNII: {drug_spec['unii']} | Scaffold: {drug_spec.get('scaffold_family', '')} | Role: {drug_spec.get('model_role', ROLE_MODEL_SELECTION_VALIDATION)} | Cohort: {drug_spec.get('cohort', 'VALIDATION_COHORT_1')}",
            status="APPROVED_REFERENCE",
            current_version=1,
        )
        db.add(comp)
        db.commit()
        db.refresh(comp)

        mw = float(Descriptors.MolWt(mol))
        clogp = float(Crippen.MolLogP(mol))
        tpsa = float(Descriptors.TPSA(mol))
        hbd = int(Lipinski.NumHDonors(mol))
        hba = int(Lipinski.NumHAcceptors(mol))
        rotb = int(Lipinski.NumRotatableBonds(mol))

        cv = CompoundVersion(
            compound_row_id=comp.id,
            version_number=1,
            original_smiles=drug_spec["smiles"],
            canonical_smiles=canon_smiles,
            isomeric_smiles=canon_smiles,
            inchi=inchi_str,
            inchikey=inchikey_str,
            change_note="Canonical reference drug registration",
            properties_json=json.dumps({
                "MW": mw, "cLogP": clogp, "TPSA": tpsa, "HBD": hbd, "HBA": hba, "RotB": rotb,
                "drugbank_id": drug_spec["drugbank_id"], "chembl_id": drug_spec["chembl_id"],
                "pubchem_cid": drug_spec["pubchem_cid"], "unii": drug_spec["unii"],
                "scaffold": drug_spec.get("scaffold_family", ""),
                "model_role": drug_spec.get("model_role", ROLE_MODEL_SELECTION_VALIDATION),
                "cohort": drug_spec.get("cohort", "VALIDATION_COHORT_1"),
            }),
        )
        db.add(cv)
        db.commit()
        db.refresh(cv)
    else:
        cv = db.scalar(select(CompoundVersion).where(CompoundVersion.compound_row_id == comp.id, CompoundVersion.version_number == 1))

    persisted_records = []
    ad_status, nearest_sim, violations, metrics, ad_reason = evaluate_safety_applicability_domain(mol)
    upstream_overlap = drug_spec.get("upstream_overlap", {})
    model_role = drug_spec.get("model_role", ROLE_MODEL_SELECTION_VALIDATION)
    cohort = drug_spec.get("cohort", "VALIDATION_COHORT_1")

    for obs in drug_spec["observations"]:
        eid = obs["canonical_endpoint_id"]
        overlap_status = upstream_overlap.get(eid, "VALIDATION_HOLDOUT" if obs["training_eligible"] else "NOT_ELIGIBLE")

        # Partitioning
        if not obs["training_eligible"]:
            partition = "NOT_ELIGIBLE"
        elif overlap_status == "EXACT_STRUCTURE_OVERLAP":
            partition = "TRAINING_ELIGIBLE"
        elif model_role == ROLE_DEVELOPMENT_TRAINING:
            partition = "DEVELOPMENT_TRAINING"
        elif model_role == ROLE_FINAL_TEST_COHORT_1_CONSUMED:
            partition = "FINAL_TEST_COHORT_1_CONSUMED"
        elif model_role == ROLE_LOCKED_FINAL_TEST_COHORT_2:
            partition = "LOCKED_FINAL_TEST_COHORT_2"
        else:
            partition = "MODEL_SELECTION_VALIDATION"

        p_key = hashlib.sha256(f"{cv.inchikey}_{eid}_{obs['raw_value']}_{obs['raw_unit']}_{obs['species']}_{obs['matrix']}".encode()).hexdigest()
        existing_ev = db.scalar(select(ExternalExperimentalEvidence).where(
            ExternalExperimentalEvidence.compound_version_id == cv.id,
            ExternalExperimentalEvidence.provenance_key == p_key
        ))

        cond_dict = {
            "matrix": obs["matrix"],
            "section": obs["section"],
            "upstream_overlap": overlap_status,
            "drugbank_partition": partition,
            "model_role": model_role,
            "cohort": cohort,
        }

        if not existing_ev:
            ev = ExternalExperimentalEvidence(
                compound_version_id=cv.id,
                provenance_key=p_key,
                cas_number=drug_spec["cas_number"],
                canonical_endpoint_id=eid,
                raw_endpoint_name=obs["raw_endpoint_name"],
                species=obs["species"],
                assay_type=obs["assay_type"],
                assay_conditions_json=cond_dict,
                raw_value=obs["raw_value"],
                raw_unit=obs["raw_unit"],
                raw_relation=obs["raw_relation"],
                normalized_value=obs["normalized_value"],
                normalized_unit=obs["normalized_unit"],
                source_database="DrugBank_FDA_ChEMBL",
                source_record_id=drug_spec["drugbank_id"],
                source_url=f"https://go.drugbank.com/drugs/{drug_spec['drugbank_id']}",
                identity_match_status="EXACT_MATCH",
                endpoint_match_status="EXACT_MATCH",
                mapping_status="EXTERNAL_EVIDENCE_ONLY",
                evidence_origin="EXPERIMENTAL_EXTERNAL",
                source_quality_class="A",
                comparability_status="DIRECTLY_COMPARABLE",
                qualification_status="QUALIFIED_FOR_GLOBAL_TRAINING" if obs["training_eligible"] else "CLINICAL_PK_COMPOSITE",
                reference_text=obs["reference_text"],
                evidence_state="AUTO_QUALIFIED_EXTERNAL",
            )
            db.add(ev)
            db.commit()
            db.refresh(ev)
            persisted_records.append(ev)
        else:
            existing_ev.assay_conditions_json = cond_dict
            db.commit()
            persisted_records.append(existing_ev)

    return {
        "status": "SUCCESS",
        "compound_name": drug_spec["name"],
        "drugbank_id": drug_spec["drugbank_id"],
        "records_ingested_n": len(persisted_records),
    }


def ingest_gefitinib_reference_drug(db: Session) -> Dict[str, Any]:
    """Ingests Gefitinib (Drug 1)."""
    return ingest_reference_drug_by_spec(db, REFERENCE_DRUGS_CATALOG[0])


def ingest_all_drugbank_reference_drugs(db: Session) -> List[Dict[str, Any]]:
    """
    Ingests all 40 reference drugs sequentially.
    """
    results = []
    for spec in REFERENCE_DRUGS_CATALOG:
        res = ingest_reference_drug_by_spec(db, spec)
        results.append(res)
    return results
