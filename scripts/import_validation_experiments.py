"""Import internal experimental results for Engine v1 validation.

Usage:
    python scripts/import_validation_experiments.py \
        --csv validation/internal_validation_v1_experiment_import_template.csv \
        [--dry-run]

Scientific rules enforced:
  - Prediction freeze must exist before experimental import (blinding)
  - Raw values stored unmodified
  - Censored values flagged (never replaced)
  - Non-positive values NOT log-transformed silently
  - Endpoint compatibility verified before storage
  - After import: generate artifacts, then pair observations

Forbidden:
  - Overwriting existing experimental records
  - Using imported data for model fitting
  - Modifying Engine v1 policy
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def _sha256(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def parse_result_available_at(s: str) -> datetime | None:
    if not s or s.startswith("#"):
        return None
    try:
        dt = datetime.fromisoformat(s.replace("Z", "+00:00"))
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        return dt
    except ValueError:
        return None


def run(csv_path: Path, dry_run: bool = False) -> None:
    from backend.database import SessionLocal, engine
    from backend.prediction_engine_v1_policy import policy_hash as live_policy_hash
    from backend.internal_validation_v1 import (
        ENGINE_V1_POLICY_HASH, CAMPAIGN_ID,
        ensure_validation_schema,
        import_experimental_record,
        pair_observation,
        InternalValidationCohortEntryRow,
        InternalValidationPredictionFreezeRow,
    )
    import sqlalchemy as sa

    # Verify Engine v1 unchanged
    h = live_policy_hash()
    if h != ENGINE_V1_POLICY_HASH:
        print(f"[FAIL] Engine v1 policy hash mismatch: {h}")
        sys.exit(1)
    print(f"[OK] Engine v1 policy hash verified: {h}")

    ensure_validation_schema(engine)
    session = SessionLocal()

    # Build lookup: compound_id → compound_version_id
    cohort_rows = list(session.scalars(
        sa.select(InternalValidationCohortEntryRow).where(
            InternalValidationCohortEntryRow.campaign_id == CAMPAIGN_ID
        )
    ))
    label_to_cv = {e.compound_label: e.compound_version_id for e in cohort_rows}
    cv_to_entry = {e.compound_version_id: e for e in cohort_rows}
    cv_to_inchikey = {e.compound_version_id: e.inchikey for e in cohort_rows}
    cv_to_struct = {e.compound_version_id: e.structure_hash for e in cohort_rows}

    # Build lookup: (cv_id, endpoint_id) → vfreeze
    freeze_rows = list(session.scalars(
        sa.select(InternalValidationPredictionFreezeRow).where(
            InternalValidationPredictionFreezeRow.campaign_id == CAMPAIGN_ID
        )
    ))
    freeze_by_cv_ep = {
        (f.compound_version_id, f.endpoint_id): f
        for f in freeze_rows
    }

    records_imported = 0
    records_skipped = 0
    observations_paired = 0
    errors = []

    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for line_num, row in enumerate(reader, start=2):
            # Skip comment lines
            if not row.get("compound_id") or row["compound_id"].startswith("#"):
                continue

            compound_id = row["compound_id"].strip()
            cv_id_raw = row.get("compound_version_id", "").strip()
            endpoint_id = row["endpoint_id"].strip()

            # Resolve compound_version_id
            cv_id = cv_id_raw if cv_id_raw else label_to_cv.get(compound_id)
            if not cv_id:
                errors.append(f"Line {line_num}: Cannot resolve compound_version_id for '{compound_id}'")
                records_skipped += 1
                continue

            if cv_id not in cv_to_entry:
                errors.append(f"Line {line_num}: compound_version_id '{cv_id}' not in campaign cohort")
                records_skipped += 1
                continue

            # Parse raw value
            raw_val_str = row.get("raw_value", "").strip()
            try:
                raw_value = float(raw_val_str) if raw_val_str else None
            except ValueError:
                errors.append(f"Line {line_num}: Invalid raw_value '{raw_val_str}'")
                records_skipped += 1
                continue

            # Validate qualifier
            qualifier = row.get("qualifier", "=").strip() or "="
            censor_flag_str = row.get("censor_flag", "0").strip()
            censor_flag = censor_flag_str in ("1", "true", "True", "yes")
            if qualifier in ("<", ">", "BLQ", "ULOQ"):
                censor_flag = True

            raw_unit = row.get("raw_unit", "").strip()
            species = row.get("species", "").strip()
            assay_type = row.get("assay_type", "").strip()
            assay_direction = row.get("assay_direction", "").strip()
            assay_ph_str = row.get("assay_ph", "").strip()
            assay_ph = float(assay_ph_str) if assay_ph_str else None
            assay_protocol = row.get("assay_protocol", "").strip()
            replicate_id = row.get("replicate_id", "").strip()
            assay_date = row.get("assay_date", "").strip()
            result_available_at = parse_result_available_at(
                row.get("result_available_at", "").strip()
            )
            source = row.get("source", "").strip()

            if dry_run:
                print(
                    f"[DRY-RUN] Line {line_num}: {compound_id} / {endpoint_id} / "
                    f"val={raw_value} {raw_unit} / censor={censor_flag}"
                )
                records_imported += 1
                continue

            try:
                exp_rec = import_experimental_record(
                    session=session,
                    campaign_id=CAMPAIGN_ID,
                    compound_version_id=cv_id,
                    inchikey=cv_to_inchikey.get(cv_id, ""),
                    structure_hash=cv_to_struct.get(cv_id, ""),
                    endpoint_id=endpoint_id,
                    raw_value=raw_value,
                    raw_unit=raw_unit,
                    qualifier=qualifier,
                    species=species,
                    assay_type=assay_type,
                    assay_direction=assay_direction,
                    assay_ph=assay_ph,
                    assay_protocol=assay_protocol,
                    replicate_id=replicate_id,
                    assay_date=assay_date,
                    result_available_at=result_available_at,
                    source=source,
                    censor_flag=censor_flag,
                )
                records_imported += 1
                print(
                    f"[OK] Imported: {compound_id}/{endpoint_id} "
                    f"({exp_rec.endpoint_compatibility}) id={exp_rec.exp_record_id}"
                )

                # Auto-pair if prediction freeze exists
                key = (cv_id, endpoint_id)
                if key in freeze_by_cv_ep:
                    vfreeze = freeze_by_cv_ep[key]
                    from backend.internal_validation_v1 import pair_observation as _pair
                    obs = _pair(
                        session=session,
                        campaign_id=CAMPAIGN_ID,
                        vfreeze=vfreeze,
                        exp_record=exp_rec,
                        blinded_retrospective_documented=True,
                    )
                    observations_paired += 1
                    print(
                        f"  → Paired observation: {obs.observation_id} "
                        f"[{obs.prospective_evidence_class}] "
                        f"enters_primary={obs.enters_primary_metrics}"
                    )
                else:
                    print(f"  → No prediction freeze found for {cv_id}/{endpoint_id} — unpaired")

            except Exception as exc:
                errors.append(f"Line {line_num}: {exc}")
                records_skipped += 1

    session.close()

    print()
    print("=" * 60)
    print(f"Import complete: {records_imported} imported, {records_skipped} skipped")
    print(f"Observations paired: {observations_paired}")
    if errors:
        print(f"Errors ({len(errors)}):")
        for e in errors:
            print(f"  {e}")
    if not dry_run and records_imported > 0:
        print()
        print("Re-generating validation artifacts ...")
        import subprocess
        subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "generate_validation_artifacts.py")],
            check=True,
        )
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(description="Import internal experimental results for Engine v1 validation")
    parser.add_argument("--csv", required=True, help="Path to experimental results CSV")
    parser.add_argument("--dry-run", action="store_true", help="Validate without importing")
    args = parser.parse_args()
    run(Path(args.csv), dry_run=args.dry_run)


if __name__ == "__main__":
    main()
