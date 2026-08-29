"""Initialize the Engine v1 internal validation campaign.

This script:
1. Creates validation DB tables
2. Creates the campaign record
3. Enrolls GLP-1 internal compounds (ORFORGLIPRON, ALENIGLIPRON, ELECOGLIPRON)
4. Registers their existing qualification_prediction_freezes as campaign freezes
5. Writes the campaign artifact JSON

Run once.  Idempotent — safe to re-run.

Scientific rules enforced:
- No experimental values are read or written here
- Policy hash verified before enrollment
- Blinding check active during freeze registration
- No model fitting, threshold change, or Engine v1 modification
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.database import SessionLocal, engine
from backend.prediction_engine_v1_policy import policy_hash as live_policy_hash
from backend.internal_validation_v1 import (
    ENGINE_V1_POLICY_HASH,
    ENGINE_V1_POLICY_VERSION,
    CAMPAIGN_ID,
    CAMPAIGN_PROTOCOL_ID,
    ensure_validation_schema,
    get_or_create_campaign,
    register_cohort_entry,
    register_prediction_freeze,
    InternalValidationCampaignRow,
)

EXPECTED_HASH = "12757ab197b5a70d8ea1754678d9a342ab0b6ea0d82f2896bebb767d686bbdeb"


def _sha256(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def verify_engine_v1() -> None:
    h = live_policy_hash()
    if h != EXPECTED_HASH:
        raise RuntimeError(
            f"ENGINE V1 POLICY HASH MISMATCH.\n"
            f"  Expected: {EXPECTED_HASH}\n"
            f"  Got:      {h}\n"
            "Engine v1 has been modified. Halt validation."
        )
    print(f"[OK] Engine v1 policy hash verified: {h}")


def run() -> None:
    print("=" * 60)
    print("ENGINE V1 INTERNAL VALIDATION CAMPAIGN INITIALIZATION")
    print("=" * 60)

    # 1. Verify policy hash
    verify_engine_v1()

    # 2. Create validation schema
    print("[INFO] Ensuring validation schema ...")
    ensure_validation_schema(engine)
    print("[OK] Validation schema ready")

    session = SessionLocal()
    try:
        # 3. Create campaign
        campaign = get_or_create_campaign(session)
        print(f"[OK] Campaign: {campaign.campaign_id}")
        print(f"     Framework status: {campaign.framework_status}")
        print(f"     Scientific status: {campaign.scientific_status}")

        # 4. Enroll GLP-1 compounds
        import sqlalchemy as sa
        from backend.production_qualification import QualificationPredictionFreezeRow

        # Load compound version data from DB
        result = session.execute(
            sa.text("""
                SELECT c.id, c.compound_id, c.name, cv.id as cv_id,
                       cv.canonical_smiles, cv.inchikey, cv.created_at,
                       c.project_id
                FROM compounds c
                JOIN compound_versions cv
                  ON cv.compound_row_id=c.id AND cv.version_number=c.current_version
                WHERE c.project_id=1
                ORDER BY c.id
            """)
        ).fetchall()

        compound_entries = {}
        for row in result:
            cid, compound_id, name, cv_id, smiles, inchikey, created_at, project_id = row
            # Structure hash is SHA-256 of canonical SMILES
            struct_hash = _sha256(smiles) if smiles else ""
            # Murcko scaffold hash (use RDKit if available, else empty)
            try:
                from rdkit import Chem
                from rdkit.Chem.Scaffolds import MurckoScaffold
                mol = Chem.MolFromSmiles(smiles)
                scaffold = MurckoScaffold.GetScaffoldForMol(mol) if mol else None
                scaffold_smi = Chem.MolToSmiles(scaffold) if scaffold else ""
                scaffold_hash = _sha256(scaffold_smi) if scaffold_smi else ""
            except Exception:
                scaffold_hash = ""

            entry = register_cohort_entry(
                session=session,
                campaign_id=CAMPAIGN_ID,
                compound_version_id=str(cv_id),
                compound_label=compound_id,  # Public compound ID
                compound_identifier=compound_id,
                inchikey=inchikey,
                structure_hash=struct_hash,
                project_label="GLP1-SM",  # De-identified project label
                chemical_series_label="GLP1-SM-PYRIDINONE",
                murcko_scaffold_hash=scaffold_hash,
            )
            compound_entries[str(cv_id)] = {
                "entry": entry,
                "compound_id": compound_id,
                "name": name,
                "inchikey": inchikey,
                "struct_hash": struct_hash,
            }
            print(f"[OK] Enrolled: {compound_id} (cv_id={cv_id}, entry_id={entry.entry_id})")

        # 5. Register prediction freezes for each compound
        # Get all existing qualification_prediction_freezes for GLP-1 project (project_id='1')
        freezes = session.execute(
            sa.text("""
                SELECT frozen_prediction_id, compound_version_id, endpoint_id,
                       strategy, prediction_value, probability, unit,
                       applicability_domain, policy_version, standardizer_version,
                       frozen_at, models_json, provenance_json
                FROM qualification_prediction_freezes
                WHERE project_id = '1'
                ORDER BY compound_version_id, endpoint_id
            """)
        ).fetchall()

        print(f"\n[INFO] Registering {len(freezes)} prediction freezes ...")

        # Load endpoint strategy info for reliability
        from backend.prediction_engine_v1_policy import policy_rows
        endpoint_info = {r["endpoint_id"]: r for r in policy_rows()}

        freeze_count = 0
        for fr in freezes:
            (fpid, cv_id, endpoint_id, strategy, pred_val, prob,
             unit, ad, policy_ver, std_ver, frozen_at_str, models_json_str, prov_str) = fr

            cv_key = str(cv_id)
            if cv_key not in compound_entries:
                continue

            entry_data = compound_entries[cv_key]
            entry = entry_data["entry"]
            ep_info = endpoint_info.get(endpoint_id, {})
            reliability = ep_info.get("reliability", "")
            limitations = ep_info.get("limitations", [])

            # Parse models for core model info
            models = json.loads(models_json_str) if models_json_str else []
            core_model_id = ""
            core_model_version = ""
            for m in models:
                if m.get("role") == "CORE":
                    core_model_id = m.get("model_id", "")
                    core_model_version = m.get("model_version", "")
                    break
            if not core_model_id and models:
                core_model_id = models[0].get("model_id", "")
                core_model_version = models[0].get("model_version", "")

            # Parse frozen_at
            if isinstance(frozen_at_str, datetime):
                frozen_at = frozen_at_str
            else:
                frozen_at_str_clean = str(frozen_at_str).replace(" ", "T")
                try:
                    frozen_at = datetime.fromisoformat(frozen_at_str_clean)
                except Exception:
                    frozen_at = datetime.now(timezone.utc)
            if frozen_at.tzinfo is None:
                frozen_at = frozen_at.replace(tzinfo=timezone.utc)

            vfreeze = register_prediction_freeze(
                session=session,
                campaign_id=CAMPAIGN_ID,
                entry_id=entry.entry_id,
                upstream_frozen_prediction_id=fpid,
                compound_version_id=cv_key,
                inchikey=entry_data["inchikey"],
                structure_hash=entry_data["struct_hash"],
                endpoint_id=endpoint_id,
                strategy=strategy,
                evidence_class=ep_info.get("evidence_class", "MODEL_PREDICTION"),
                prediction_value=pred_val,
                probability=prob,
                unit=unit,
                applicability_domain=ad,
                reliability=reliability,
                freeze_timestamp=frozen_at,
                core_model_id=core_model_id,
                core_model_version=core_model_version,
                limitations=limitations,
            )
            freeze_count += 1

        print(f"[OK] Registered {freeze_count} campaign prediction freezes")

        # 6. Update campaign counts
        campaign.compound_count = len(compound_entries)
        campaign.prediction_freeze_complete = True
        session.commit()
        print(f"[OK] Campaign prediction_freeze_complete = True")

        # 7. Write campaign artifact
        _write_campaign_artifact(session, campaign, compound_entries)
        print(f"[OK] Campaign artifact written")

    finally:
        session.close()

    print()
    print("=" * 60)
    print("INITIALIZATION COMPLETE")
    print(f"Campaign: {CAMPAIGN_ID}")
    print("Next step: collect internal experimental data and import via")
    print("  scripts/import_validation_experiments.py")
    print("=" * 60)


def _write_campaign_artifact(session, campaign, compound_entries):
    """Write validation/internal_validation_v1_campaign.json"""
    import sqlalchemy as sa

    entries_out = []
    for cv_id, data in compound_entries.items():
        entry = data["entry"]
        entries_out.append({
            "entry_id": entry.entry_id,
            "campaign_id": entry.campaign_id,
            "compound_version_id": entry.compound_version_id,
            "compound_label": entry.compound_label,
            "inchikey": entry.inchikey,
            "project_label": entry.project_label,
            "chemical_series_label": entry.chemical_series_label,
            "eligibility_status": entry.eligibility_status,
            "enrolled_at": entry.enrolled_at.isoformat() if entry.enrolled_at else None,
        })

    artifact = {
        "campaign_id": campaign.campaign_id,
        "name": campaign.name,
        "protocol_id": campaign.protocol_id,
        "engine_policy_id": campaign.engine_policy_id,
        "engine_policy_version": campaign.engine_policy_version,
        "engine_policy_hash": campaign.engine_policy_hash,
        "standardizer_version": campaign.standardizer_version,
        "framework_status": campaign.framework_status,
        "scientific_status": campaign.scientific_status,
        "status": campaign.status,
        "bootstrap_seed": campaign.bootstrap_seed,
        "prediction_freeze_complete": campaign.prediction_freeze_complete,
        "experiment_import_complete": campaign.experiment_import_complete,
        "analysis_complete": campaign.analysis_complete,
        "compound_count": campaign.compound_count,
        "created_at": campaign.created_at.isoformat() if campaign.created_at else None,
        "cohort": entries_out,
        "_note": (
            "Internal compound structures are NOT stored in this artifact. "
            "InChIKeys are stable, non-reversible identifiers sufficient for matching."
        ),
    }

    out_path = ROOT / "validation" / "internal_validation_v1_campaign.json"
    out_path.write_text(json.dumps(artifact, indent=2))
    print(f"[OK] Written: {out_path}")


if __name__ == "__main__":
    run()
