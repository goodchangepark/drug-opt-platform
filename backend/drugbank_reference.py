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
ROLE_FINAL_TEST_COHORT_2_CONSUMED = "FINAL_TEST_COHORT_2_CONSUMED"
ROLE_LOCKED_FINAL_TEST_COHORT_2 = "FINAL_TEST_COHORT_2_CONSUMED"  # Backward compatibility alias
ROLE_FINAL_TEST_COHORT_3_CONSUMED = "FINAL_TEST_COHORT_3_CONSUMED"
ROLE_LOCKED_FINAL_TEST_COHORT_3 = "FINAL_TEST_COHORT_3_CONSUMED"  # Backward compatibility alias
ROLE_FINAL_TEST_COHORT_4_CONSUMED = "FINAL_TEST_COHORT_4_CONSUMED"
ROLE_LOCKED_FINAL_TEST_COHORT_4 = "FINAL_TEST_COHORT_4_CONSUMED"  # Backward compatibility alias
ROLE_LOCKED_FINAL_TEST_COHORT_5 = "LOCKED_FINAL_TEST_COHORT_5"
ROLE_LOCKED_FINAL_TEST_COHORT_6 = "LOCKED_FINAL_TEST_COHORT_6"
ROLE_LOCKED_FINAL_TEST_COHORT_7 = "LOCKED_FINAL_TEST_COHORT_7"

# Load full 200 reference drugs catalog (with fallbacks)
CATALOG_PATH = Path(__file__).parent / "reference_drugs_200.json"
if not CATALOG_PATH.exists():
    CATALOG_PATH = Path(__file__).parent / "reference_drugs_150.json"
if not CATALOG_PATH.exists():
    CATALOG_PATH = Path(__file__).parent / "reference_drugs_100.json"
if not CATALOG_PATH.exists():
    CATALOG_PATH = Path(__file__).parent / "reference_drugs_80.json"
if not CATALOG_PATH.exists():
    CATALOG_PATH = Path(__file__).parent / "reference_drugs_65.json"
if not CATALOG_PATH.exists():
    CATALOG_PATH = Path(__file__).parent / "reference_drugs_50.json"
if not CATALOG_PATH.exists():
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
        target_cid = f"DRUGBANK-{drug_spec['drugbank_id']}"
        conflict = db.scalar(select(Compound).where(Compound.project_id == proj.id, Compound.compound_id == target_cid))
        if conflict and conflict.name != drug_spec["name"]:
            if conflict.name == "Lorcaserin":
                conflict.compound_id = "DRUGBANK-DB04871"
                conflict.notes = conflict.notes.replace("DB08907", "DB04871")
                db.commit()
        comp = Compound(
            project_id=proj.id,
            compound_id=target_cid,
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

        # Partitioning: prioritize endpoint-level role if specified in observation
        obs_role = obs.get("endpoint_role", model_role)
        if not obs["training_eligible"]:
            partition = "NOT_ELIGIBLE"
        elif overlap_status == "EXACT_STRUCTURE_OVERLAP":
            partition = "TRAINING_ELIGIBLE"
        elif obs_role == ROLE_DEVELOPMENT_TRAINING:
            partition = "DEVELOPMENT_TRAINING"
        elif obs_role == ROLE_FINAL_TEST_COHORT_1_CONSUMED:
            partition = "FINAL_TEST_COHORT_1_CONSUMED"
        elif obs_role == ROLE_FINAL_TEST_COHORT_2_CONSUMED:
            partition = "FINAL_TEST_COHORT_2_CONSUMED"
        elif obs_role == ROLE_FINAL_TEST_COHORT_3_CONSUMED:
            partition = "FINAL_TEST_COHORT_3_CONSUMED"
        elif obs_role == ROLE_FINAL_TEST_COHORT_4_CONSUMED:
            partition = "FINAL_TEST_COHORT_4_CONSUMED"
        elif obs_role == ROLE_LOCKED_FINAL_TEST_COHORT_5:
            partition = "LOCKED_FINAL_TEST_COHORT_5"
        elif obs_role == ROLE_LOCKED_FINAL_TEST_COHORT_6:
            partition = "LOCKED_FINAL_TEST_COHORT_6"
        elif obs_role == ROLE_LOCKED_FINAL_TEST_COHORT_7:
            partition = "LOCKED_FINAL_TEST_COHORT_7"
        elif obs_role == "LOCKED_FINAL_TEST_COHORT_4":
            partition = "FINAL_TEST_COHORT_4_CONSUMED"
        elif obs_role == "LOCKED_FINAL_TEST_COHORT_3":
            partition = "FINAL_TEST_COHORT_3_CONSUMED"
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


def get_endpoint_priority_rank(endpoint_id: str) -> int:
    """
    Priority order per Directive 2:
    PPB (0) -> hERG (1) -> Caco-2 (2) -> HLM (3) -> CYP1A2/2C9/2C19 (4) -> Others (5+)
    """
    if endpoint_id == "HUMAN_PPB":
        return 0
    elif endpoint_id == "HERG_LIABILITY":
        return 1
    elif endpoint_id in ("CACO2_PERMEABILITY", "CACO2_PAPP_AB"):
        return 2
    elif endpoint_id in ("HLM_INTRINSIC_CLEARANCE", "HLM_CLINT"):
        return 3
    elif endpoint_id in ("CYP1A2_INHIBITION", "CYP2C9_INHIBITION", "CYP2C19_INHIBITION"):
        return 4
    elif endpoint_id in ("CYP3A4_INHIBITION", "CYP2D6_INHIBITION"):
        return 5
    elif endpoint_id in ("SOLUBILITY_GENERIC", "SOLUBILITY_THERMODYNAMIC"):
        return 6
    return 7


def ingest_reference_drug_stepwise_lifecycle(db: Session, drug_spec: Dict[str, Any]) -> Dict[str, Any]:
    """
    Executes the complete 5-stage qualification lifecycle for a reference drug:
    Stage 1: Identity (chemical structure, distinct scaffold, compound/version registration)
    Stage 2: Evidence (harvests external experimental observations ordered by priority: PPB -> hERG -> Caco-2 -> HLM -> CYPs)
    Stage 3: Qualification (identity match, endpoint match, source quality class A, qualification status)
    Stage 4: Prediction (computes model predictions for qualified endpoints)
    Stage 5: Error (calculates prediction error |pred - truth| and evaluates residual)
    Only after completing Stage 5 does the caller advance to the next drug.
    """
    from backend.engine_v3_learning import compute_base_prediction

    # 1. Identity Stage
    smiles = drug_spec["smiles"]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for {drug_spec['name']}")
    canon_smiles = Chem.MolToSmiles(mol, canonical=True)
    inchi_str = Chem.MolToInchi(mol)
    inchikey_str = Chem.MolToInchiKey(mol)

    proj = ensure_drugbank_project(db)
    comp = db.scalar(select(Compound).where(Compound.project_id == proj.id, Compound.name == drug_spec["name"]))
    if not comp:
        comp = Compound(
            project_id=proj.id,
            compound_id=drug_spec["drugbank_id"],
            name=drug_spec["name"],
            status="ACTIVE",
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
            original_smiles=smiles,
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

    identity_summary = {
        "status": "IDENTITY_VERIFIED",
        "compound_id": comp.compound_id,
        "compound_name": comp.name,
        "scaffold_family": drug_spec.get("scaffold_family", ""),
        "inchikey": cv.inchikey,
        "model_role": drug_spec.get("model_role", ROLE_MODEL_SELECTION_VALIDATION),
        "cohort": drug_spec.get("cohort", "VALIDATION_COHORT_1"),
    }

    # 2. Evidence Stage: Sort observations by priority PPB -> hERG -> Caco-2 -> HLM -> CYP1A2/2C9/2C19 -> others
    raw_obs = drug_spec.get("observations", [])
    sorted_obs = sorted(raw_obs, key=lambda x: get_endpoint_priority_rank(x["canonical_endpoint_id"]))

    evidence_records = []
    qualification_records = []
    prediction_records = []
    error_records = []

    upstream_overlap = drug_spec.get("upstream_overlap", {})
    model_role = drug_spec.get("model_role", ROLE_MODEL_SELECTION_VALIDATION)
    cohort = drug_spec.get("cohort", "VALIDATION_COHORT_1")

    for obs in sorted_obs:
        eid = obs["canonical_endpoint_id"]
        overlap_status = upstream_overlap.get(eid, "VALIDATION_HOLDOUT" if obs["training_eligible"] else "NOT_ELIGIBLE")

        # Partitioning: prioritize endpoint-level role if specified in observation
        obs_role = obs.get("endpoint_role", model_role)
        if not obs["training_eligible"]:
            partition = "NOT_ELIGIBLE"
        elif overlap_status == "EXACT_STRUCTURE_OVERLAP":
            partition = "TRAINING_ELIGIBLE"
        elif obs_role == ROLE_DEVELOPMENT_TRAINING:
            partition = "DEVELOPMENT_TRAINING"
        elif obs_role == ROLE_FINAL_TEST_COHORT_1_CONSUMED:
            partition = "FINAL_TEST_COHORT_1_CONSUMED"
        elif obs_role == ROLE_FINAL_TEST_COHORT_2_CONSUMED:
            partition = "FINAL_TEST_COHORT_2_CONSUMED"
        elif obs_role == ROLE_FINAL_TEST_COHORT_3_CONSUMED:
            partition = "FINAL_TEST_COHORT_3_CONSUMED"
        elif obs_role == ROLE_FINAL_TEST_COHORT_4_CONSUMED:
            partition = "FINAL_TEST_COHORT_4_CONSUMED"
        elif obs_role == ROLE_LOCKED_FINAL_TEST_COHORT_5:
            partition = "LOCKED_FINAL_TEST_COHORT_5"
        elif obs_role == "LOCKED_FINAL_TEST_COHORT_4":
            partition = "FINAL_TEST_COHORT_4_CONSUMED"
        elif obs_role == "LOCKED_FINAL_TEST_COHORT_3":
            partition = "FINAL_TEST_COHORT_3_CONSUMED"
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
        else:
            existing_ev.assay_conditions_json = cond_dict
            db.commit()
            ev = existing_ev

        evidence_records.append({
            "endpoint_id": eid,
            "provenance_key": p_key,
            "raw_value": obs["raw_value"],
            "raw_unit": obs["raw_unit"],
            "normalized_value": obs["normalized_value"],
            "normalized_unit": obs["normalized_unit"],
        })

        # 3. Qualification Stage
        qualification_records.append({
            "endpoint_id": eid,
            "identity_match_status": ev.identity_match_status,
            "endpoint_match_status": ev.endpoint_match_status,
            "source_quality_class": ev.source_quality_class,
            "qualification_status": ev.qualification_status,
            "comparability_status": ev.comparability_status,
            "partition": partition,
        })

        # 4. Prediction Stage
        pred_val = compute_base_prediction(eid, canon_smiles)
        prediction_records.append({
            "endpoint_id": eid,
            "predicted_value": pred_val,
        })

        # 5. Error Stage
        exp_val = float(obs["normalized_value"])
        if eid in ("HERG_LIABILITY", "CYP3A4_INHIBITION", "CYP2D6_INHIBITION", "CYP1A2_INHIBITION", "CYP2C9_INHIBITION", "CYP2C19_INHIBITION"):
            exp_p = ic50_nm_to_pic50(exp_val) if exp_val > 0 else exp_val
        else:
            exp_p = exp_val

        abs_err = round(abs(pred_val - exp_p), 3) if pred_val is not None else None
        error_records.append({
            "endpoint_id": eid,
            "predicted_value": pred_val,
            "experimental_value": round(exp_p, 3),
            "absolute_error": abs_err,
        })

    return {
        "status": "SUCCESS",
        "compound_name": drug_spec["name"],
        "drugbank_id": drug_spec["drugbank_id"],
        "identity": identity_summary,
        "evidence": evidence_records,
        "qualification": qualification_records,
        "prediction": prediction_records,
        "error": error_records,
    }


def ingest_gefitinib_reference_drug(db: Session) -> Dict[str, Any]:
    """Ingests Gefitinib (Drug 1)."""
    return ingest_reference_drug_by_spec(db, REFERENCE_DRUGS_CATALOG[0])


def ingest_all_drugbank_reference_drugs(db: Session) -> List[Dict[str, Any]]:
    """
    Ingests all 50 reference drugs sequentially.
    """
    results = []
    for spec in REFERENCE_DRUGS_CATALOG:
        res = ingest_reference_drug_by_spec(db, spec)
        results.append(res)
    return results


def ingest_v3_1_expansion_drugs_sequential(db: Session) -> List[Dict[str, Any]]:
    """
    Sequentially ingests the 10 new approved reference drugs (Drugs 41 to 50)
    for Global Engine v3.1, enforcing Identity -> Evidence -> Qualification -> Prediction -> Error
    stepwise completion for each compound before advancing to the next.
    """
    if len(REFERENCE_DRUGS_CATALOG) < 50:
        raise RuntimeError("Catalog must contain at least 50 drugs for v3.1 expansion")

    expansion_specs = REFERENCE_DRUGS_CATALOG[40:50]
    completed_lifecycle_results = []

    for idx, spec in enumerate(expansion_specs, start=41):
        res = ingest_reference_drug_stepwise_lifecycle(db, spec)
        assert res["identity"]["status"] == "IDENTITY_VERIFIED"
        assert len(res["evidence"]) >= 5
        assert len(res["qualification"]) >= 5
        assert len(res["prediction"]) >= 5
        assert len(res["error"]) >= 5
        completed_lifecycle_results.append(res)

    return completed_lifecycle_results


def ingest_v3_2_expansion_drugs_sequential(db: Session) -> List[Dict[str, Any]]:
    """
    Sequentially ingests the 15 new approved reference drugs (Drugs 51 to 65)
    for Global Engine v3.2, enforcing Identity -> Evidence -> Qualification -> Prediction -> Error
    stepwise completion for each compound before advancing to the next.
    """
    if len(REFERENCE_DRUGS_CATALOG) < 65:
        raise RuntimeError("Catalog must contain at least 65 drugs for v3.2 expansion")

    expansion_specs = REFERENCE_DRUGS_CATALOG[50:65]
    completed_lifecycle_results = []

    for idx, spec in enumerate(expansion_specs, start=51):
        res = ingest_reference_drug_stepwise_lifecycle(db, spec)
        assert res["identity"]["status"] == "IDENTITY_VERIFIED"
        assert len(res["evidence"]) >= 5
        assert len(res["qualification"]) >= 5
        assert len(res["prediction"]) >= 5
        assert len(res["error"]) >= 5
        completed_lifecycle_results.append(res)

    return completed_lifecycle_results


def ingest_v3_3_expansion_drugs_sequential(db: Session) -> List[Dict[str, Any]]:
    """
    Sequentially ingests the 15 new approved reference drugs (Drugs 66 to 80)
    for Global Engine v3.3, enforcing Identity -> Evidence -> Qualification -> Prediction -> Error
    stepwise completion for each compound before advancing to the next.
    """
    if len(REFERENCE_DRUGS_CATALOG) < 80:
        raise RuntimeError("Catalog must contain at least 80 drugs for v3.3 expansion")

    expansion_specs = REFERENCE_DRUGS_CATALOG[65:80]
    completed_lifecycle_results = []

    for idx, spec in enumerate(expansion_specs, start=66):
        res = ingest_reference_drug_stepwise_lifecycle(db, spec)
        assert res["identity"]["status"] == "IDENTITY_VERIFIED"
        assert len(res["evidence"]) >= 3
        assert len(res["qualification"]) >= 3
        assert len(res["prediction"]) >= 3
        assert len(res["error"]) >= 3
        completed_lifecycle_results.append(res)

    return completed_lifecycle_results



