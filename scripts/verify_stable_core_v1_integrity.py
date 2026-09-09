#!/usr/bin/env python3
"""Read-only physical and business-integrity verification for Stable Core v1."""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

PROTECTED_PROJECTS = {1: "GLP-1 (small molecule)", 3: "EGFR", 5: "AMYR (small molecules)", 300: "DrugBank"}
PROJECT5_COMPOUNDS = (11, 12, 13, 14)
PROJECT5_VERSIONS = (14, 15, 16, 17)
RECOVERED_RUNS = (64, 65, 66, 67, 68, 73, 74, 75, 127, 128, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162)
LEGITIMATE_POST_BASELINE_RUNS = (167, 168)
PRE_STABLE_CORE_UNINTENDED_RUNS = (169, 170)
CURRENT_ENGINE_ID = "drugopt-prediction-engine-v3@3.3.3"
VALID_MODES = {"ASSISTED", "HYBRID", "FULL_PREDICTION"}


def rows(connection, sql, parameters=()):
    return [list(row) for row in connection.execute(sql, parameters)]


def verify(database: Path) -> dict:
    connection = sqlite3.connect(f"file:{database.resolve()}?mode=ro", uri=True)
    try:
        projects = dict(connection.execute("SELECT id,name FROM projects ORDER BY id"))
        p5_compounds = tuple(row[0] for row in connection.execute("SELECT id FROM compounds WHERE project_id=5 ORDER BY id"))
        p5_versions = tuple(row[0] for row in connection.execute(
            "SELECT id FROM compound_versions WHERE compound_row_id IN (11,12,13,14) ORDER BY id"
        ))
        prediction_ids = tuple(row[0] for row in connection.execute("SELECT id FROM prediction_runs ORDER BY id"))
        historical_ids = tuple(row[0] for row in connection.execute(
            "SELECT legacy_prediction_run_id FROM historical_predictions ORDER BY legacy_prediction_run_id"
        ))
        project300_count = connection.execute("SELECT COUNT(*) FROM compounds WHERE project_id=300").fetchone()[0]
        project300_unique = connection.execute("""
          SELECT COUNT(DISTINCT cv.inchikey)
          FROM compounds c JOIN compound_versions cv
            ON cv.compound_row_id=c.id AND cv.version_number=c.current_version
          WHERE c.project_id=300 AND cv.inchikey IS NOT NULL AND trim(cv.inchikey)!=''
        """).fetchone()[0]
        orfor = rows(connection, """
          SELECT o.canonical_endpoint,o.species,o.value_text,o.unit,json_extract(o.context_json,'$.route'),
                 o.curation_status,o.source_record_id
          FROM experimental_observations o
          WHERE o.compound_version_id=11 AND o.curation_status='ACCEPTED'
          ORDER BY o.canonical_endpoint
        """)
        representative_samples = rows(connection, """
          WITH selected AS (
            SELECT c.project_id, MIN(c.id) AS compound_id
            FROM compounds c
            WHERE c.project_id IN (1,3,5,300)
            GROUP BY c.project_id
          )
          SELECT c.project_id,c.id,c.compound_id,c.name,cv.id,cv.inchikey,
                 (SELECT COUNT(*) FROM experimental_observations o
                   WHERE o.compound_version_id=cv.id AND o.curation_status='ACCEPTED') AS accepted_observations,
                 (SELECT COUNT(*) FROM current_prediction_snapshots s
                   WHERE s.compound_version_id=cv.id AND s.is_current=1
                     AND s.engine_release='drugopt-prediction-engine-v3@3.3.3') AS current_predictions,
                 (SELECT COUNT(*) FROM historical_predictions h
                   WHERE h.compound_version_id_snapshot=cv.id) AS historical_predictions
          FROM selected x
          JOIN compounds c ON c.id=x.compound_id
          JOIN compound_versions cv ON cv.compound_row_id=c.id AND cv.version_number=c.current_version
          ORDER BY c.project_id
        """)
        current_snapshot_total = connection.execute(
            "SELECT COUNT(*) FROM current_prediction_snapshots"
        ).fetchone()[0]
        current_v333_snapshot_total = connection.execute(
            """SELECT COUNT(*) FROM current_prediction_snapshots
               WHERE is_current=1 AND engine_release=?""",
            (CURRENT_ENGINE_ID,),
        ).fetchone()[0]
        duplicate_current_snapshot_keys = rows(connection, """
          SELECT compound_version_id,canonical_endpoint,species,context_identity,engine_release,COUNT(*)
          FROM current_prediction_snapshots
          GROUP BY compound_version_id,canonical_endpoint,species,context_identity,engine_release
          HAVING COUNT(*) > 1
        """)
        invalid_v333_current = rows(connection, """
          SELECT id,canonical_endpoint,species,model_id,model_version,model_artifact_hash,prediction_mode
          FROM current_prediction_snapshots
          WHERE is_current=1 AND engine_release=?
            AND (
              model_id IS NULL OR trim(model_id)='' OR upper(model_id)='UNKNOWN_PROVENANCE' OR
              model_version IS NULL OR trim(model_version)='' OR upper(model_version)='UNKNOWN_PROVENANCE' OR
              model_artifact_hash IS NULL OR trim(model_artifact_hash)='' OR upper(model_artifact_hash)='UNKNOWN_PROVENANCE' OR
              upper(coalesce(prediction_mode,'')) NOT IN ('ASSISTED','HYBRID','FULL_PREDICTION')
            )
        """, (CURRENT_ENGINE_ID,))
        from backend.canonical_endpoints import REGISTRY
        species_contradictions = []
        for row in connection.execute("SELECT id,canonical_endpoint,species FROM current_prediction_snapshots WHERE is_current=1 AND engine_release=?", (CURRENT_ENGINE_ID,)):
            definition = REGISTRY.get(str(row[1]).upper())
            if definition and definition.species_requirement and str(row[2]).upper() != definition.species_requirement:
                species_contradictions.append(list(row))
        unknown_provenance_current_total = connection.execute(
            """SELECT COUNT(*) FROM current_prediction_snapshots
               WHERE is_current=1 AND engine_release='UNKNOWN_PROVENANCE'"""
        ).fetchone()[0]
        triggers = {row[0] for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type='trigger' AND name LIKE 'stable_core_%'"
        )}
        required_triggers = {
            "stable_core_protect_real_project_delete",
            "stable_core_protect_real_compound_delete",
            "stable_core_protect_historical_run_delete",
            "stable_core_protect_historical_run_update",
            "stable_core_protect_history_delete",
            "stable_core_protect_history_update",
            "stable_core_mirror_prediction_run_insert",
        }
        checks = {
            "required_protected_projects_present": all(pid in projects for pid in PROTECTED_PROJECTS),
            "project5_compounds_exact": p5_compounds == PROJECT5_COMPOUNDS,
            "project5_versions_exact": p5_versions == PROJECT5_VERSIONS,
            "project5_runs_restored": set(RECOVERED_RUNS).issubset(prediction_ids),
            "later_runs_167_168_preserved": set(LEGITIMATE_POST_BASELINE_RUNS).issubset(prediction_ids),
            "pre_stable_core_unintended_runs_preserved_for_provenance": set(PRE_STABLE_CORE_UNINTENDED_RUNS).issubset(prediction_ids),
            "history_identity_parity": prediction_ids == historical_ids,
            "recovered_identity_baseline_present": len(prediction_ids) >= 142,
            "project300_compounds_1000": project300_count == 1000,
            "project300_unique_current_inchikeys_1000": project300_unique == 1000,
            "orforglipron_four_accepted_pk": {
                row[0] for row in orfor
            } >= {"HUMAN_PK_F_ORAL", "HUMAN_PK_VD_IV", "HUMAN_PK_CL_UNSPECIFIED", "HUMAN_PK_CMAX_UNSPECIFIED"},
            "orforglipron_clearance_route_unspecified": any(
                row[0] == "HUMAN_PK_CL_UNSPECIFIED" and row[4] == "UNSPECIFIED" for row in orfor
            ),
            "required_safety_triggers": required_triggers.issubset(triggers),
            "representative_all_protected_projects": {row[0] for row in representative_samples} == {1, 3, 5, 300},
            "current_prediction_scientific_key_unique": not duplicate_current_snapshot_keys,
            "v333_invalid_current_snapshot_zero": not invalid_v333_current,
            "v333_current_species_consistent": not species_contradictions,
        }
        result = {
            "contract": "StableCoreBusinessIntegrity/v1",
            "database": str(database.resolve()),
            "sqlite_user_version": connection.execute("PRAGMA user_version").fetchone()[0],
            "integrity_check": connection.execute("PRAGMA integrity_check").fetchone()[0],
            "foreign_key_check": rows(connection, "PRAGMA foreign_key_check"),
            "projects": projects,
            "project5_compound_ids": p5_compounds,
            "project5_version_ids": p5_versions,
            "prediction_run_total": len(prediction_ids),
            "historical_prediction_total": len(historical_ids),
            "prediction_run_identity_policy": {
                "required_recovered_ids": list(RECOVERED_RUNS),
                "legitimate_post_baseline_ids": list(LEGITIMATE_POST_BASELINE_RUNS),
                "pre_stable_core_unintended_ids_preserved_for_provenance": list(PRE_STABLE_CORE_UNINTENDED_RUNS),
                "note": "Stable Core validates immutable identities rather than forcing a historical count. Runs 169/170 predate the final isolated-test gate and are retained pending provenance review.",
            },
            "project300": {"compound_count": project300_count, "unique_current_inchikey_count": project300_unique},
            "orforglipron_accepted_observations": orfor,
            "representative_parity_samples": representative_samples,
            "current_prediction_snapshots": {
                "total": current_snapshot_total,
                "current_v333": current_v333_snapshot_total,
                "duplicate_scientific_keys": duplicate_current_snapshot_keys,
                "invalid_v333_current": invalid_v333_current,
                "species_contradictions": species_contradictions,
                "legacy_unknown_provenance_current": unknown_provenance_current_total,
                "selection_policy": "Only exact current engine_release rows are eligible for Current Prediction; UNKNOWN_PROVENANCE is excluded.",
            },
            "safety_triggers": sorted(triggers),
            "checks": checks,
        }
        result["status"] = "PASS" if all(checks.values()) and result["integrity_check"] == "ok" and not result["foreign_key_check"] else "FAIL"
        if result["status"] != "PASS":
            raise RuntimeError(json.dumps(result, indent=2))
        return result
    finally:
        connection.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", default=str(Path(__file__).resolve().parents[1] / "drug_opt.db"))
    parser.add_argument("--output")
    args = parser.parse_args()
    result = verify(Path(args.database))
    rendered = json.dumps(result, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        Path(args.output).write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
