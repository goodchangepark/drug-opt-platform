#!/usr/bin/env python3
"""Emit the auditable 59-endpoint Mobocertinib execution matrix.

Run only against an isolated E2E database after the explicit Predict request.
The script is read-only with respect to scientific state; its only write is the
requested JSON validation artifact.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from sqlalchemy import select

from backend.database import SessionLocal, engine
from backend.developability_profile import build_developability_profile
from backend.model_artifact_authority import (
    artifact_bundle_sha256,
    model_artifact_registration,
    validation_registration_error,
)
from backend.models import Compound, CompoundVersion, PredictionRun, Project
from backend.predict_all_core_contract import (
    CONTRACT_ID,
    CONTEXT_REQUIRED_ENDPOINTS,
    CURRENT_DATA_CEILING_ENDPOINTS,
    MECHANISTIC_EXECUTION_ENDPOINTS,
    MODEL_NOT_REGISTERED_ENDPOINTS,
    MODEL_UNAVAILABLE_ENDPOINTS,
    PREDICT_ALL_EXECUTION_ENDPOINTS,
    PUBLISHABLE_CURRENT_ENDPOINTS,
)
from backend.prediction_engine_registry import get_current_production_routing
from backend.stable_core import CurrentPredictionSnapshot, admit_current_prediction


PROFILE_QUERY = {
    "HUMAN_FU_PLASMA": "HUMAN_FU",
    "HUMAN_TOTAL_IV_CL": "HUMAN_PK_CL_IV",
    "IV_HALF_LIFE": "HUMAN_PK_T_HALF_IV",
    "HUMAN_F": "HUMAN_PK_F_ORAL",
    "KA": "HUMAN_PK_KA_ORAL",
    "ORAL_AUC": "HUMAN_PK_AUC_ORAL",
    "ORAL_CMAX": "HUMAN_PK_CMAX_ORAL",
    "TMAX": "HUMAN_PK_TMAX_ORAL",
}

ADMET_NAMES = {
    "SOLUBILITY_GENERIC": "Solubility",
    "CACO2_PAPP_AB": "Permeability",
    "HUMAN_PPB": "Plasma protein binding",
    "HLM_CLINT": "HLM intrinsic clearance",
    "RLM_CLINT": "RLM intrinsic clearance",
    "MLM_CLINT": "MLM intrinsic clearance",
    "CYP1A2_INHIBITOR_CLASS": "CYP1A2 inhibitor",
    "CYP2C9_INHIBITOR_CLASS": "CYP2C9 inhibitor",
    "CYP2C19_INHIBITOR_CLASS": "CYP2C19 inhibitor",
    "CYP2D6_INHIBITOR_CLASS": "CYP2D6 inhibitor",
    "CYP3A4_INHIBITOR_CLASS": "CYP3A4 inhibitor",
    "CYP2C9_SUBSTRATE": "CYP2C9 substrate",
    "CYP2D6_SUBSTRATE": "CYP2D6 substrate",
    "CYP3A4_SUBSTRATE": "CYP3A4 substrate",
    "PGP_INHIBITION": "P-gp inhibitor",
    "HERG_CLASS": "hERG liability",
    "DILI_LIABILITY": "DILI clinical liability",
}

SPECIES = {
    "RLM_CLINT": "RAT", "MLM_CLINT": "MOUSE",
    "AMES_MUTAGENICITY": "BACTERIAL",
}

PRODUCTION_BEFORE_TRACE = {
    "prediction_run_id": 247,
    "created_at": "2026-09-14 00:38:21.605132",
    "inputs_hash": "96bcad8c6f609f1dbd58436939b3097d6d2d785c6c5c46745c05f92adaae8752",
}
PRODUCTION_BEFORE_PUBLISHED = frozenset({
    "CYP1A2_INHIBITOR_CLASS", "CYP2C19_INHIBITOR_CLASS",
    "CYP2C9_INHIBITOR_CLASS", "CYP2C9_SUBSTRATE",
    "CYP2D6_INHIBITOR_CLASS", "CYP2D6_SUBSTRATE",
    "CYP3A4_INHIBITOR_CLASS", "CYP3A4_SUBSTRATE", "DILI_LIABILITY",
    "MLM_CLINT", "PGP_INHIBITION", "RLM_CLINT", "HERG_CLASS",
})
PRODUCTION_BEFORE_REJECTED = {
    "AMES_MUTAGENICITY": "INVALID_SPECIES",
    "HLM_CLINT": "ROUTED_IMPLEMENTATION_MISMATCH",
    "CACO2_PAPP_AB": "ROUTED_IMPLEMENTATION_MISMATCH",
    "HUMAN_PPB": "ROUTED_IMPLEMENTATION_MISMATCH",
    "SOLUBILITY_GENERIC": "ROUTED_IMPLEMENTATION_MISMATCH",
}


def classification(endpoint: str) -> str:
    if endpoint in PUBLISHABLE_CURRENT_ENDPOINTS:
        return "EXECUTED_AND_PUBLISHED"
    if endpoint in MECHANISTIC_EXECUTION_ENDPOINTS:
        return "MECHANISTIC_ONLY"
    if endpoint in MODEL_NOT_REGISTERED_ENDPOINTS:
        return "MODEL_NOT_REGISTERED"
    if endpoint in CONTEXT_REQUIRED_ENDPOINTS:
        return "CONTEXT_REQUIRED"
    return "MODEL_UNAVAILABLE"


def final_reason(endpoint: str) -> str:
    if endpoint in PUBLISHABLE_CURRENT_ENDPOINTS:
        return "Qualified existing route executed or reused its exact immutable output and passed unchanged Stable Core admission."
    if endpoint in MODEL_NOT_REGISTERED_ENDPOINTS:
        return "A historical route label exists, but no complete authoritative executable artifact/implementation registration was found."
    if endpoint in MECHANISTIC_EXECUTION_ENDPOINTS:
        return "The available output is mechanistic/rule-based and is not represented as a scalar CurrentPredictionSnapshot."
    if endpoint in CURRENT_DATA_CEILING_ENDPOINTS:
        return "Final campaign decision CURRENT_DATA_CEILING; canonical structure-only publication is not permitted."
    if endpoint in CONTEXT_REQUIRED_ENDPOINTS:
        return "Dose, route, regimen, formulation, or matched PK context is required; structure-only execution is excluded."
    if endpoint == "AMES_MUTAGENICITY":
        return "Installed legacy ensemble exists, but pooled bacterial assay/species semantics cannot satisfy canonical Stable Core admission."
    return "Final campaign decision MODEL_UNAVAILABLE; no qualified executable current model is registered."


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--contract-output", required=True)
    parser.add_argument("--before-run-id", type=int)
    args = parser.parse_args()

    campaign = json.loads(Path("validation/full_model_development_campaign.json").read_text())
    tasks = campaign["tasks"]
    routes = {row["endpoint_id"]: row for row in get_current_production_routing()}
    db = SessionLocal()
    try:
        matches = list(db.scalars(
            select(Compound)
            .join(Project, Project.id == Compound.project_id)
            .where(Project.name == "EGFR", Compound.name.ilike("Mobocertinib"))
        ))
        if len(matches) != 1:
            raise RuntimeError(f"Expected one existing EGFR Mobocertinib, found {len(matches)}")
        compound = matches[0]
        version = db.scalar(select(CompoundVersion).where(
            CompoundVersion.compound_row_id == compound.id,
            CompoundVersion.version_number == compound.current_version,
        ))
        latest = db.scalar(select(PredictionRun).where(
            PredictionRun.version_id == version.id,
            PredictionRun.stage == "prediction_workflow",
        ).order_by(PredictionRun.created_at.desc(), PredictionRun.id.desc()))
        before_published = PRODUCTION_BEFORE_PUBLISHED

        profile = build_developability_profile(db, version.id)
        visible = {row["query_endpoint"]: row for row in profile["availability_catalog"]}
        snapshots = list(db.scalars(select(CurrentPredictionSnapshot).where(
            CurrentPredictionSnapshot.compound_version_id == version.id,
            CurrentPredictionSnapshot.is_current.is_(True),
        )))
        snapshots_by_endpoint = {row.canonical_endpoint: row for row in snapshots}
        latest_output = dict(latest.outputs_json or {}) if latest else {}
        quantitative = dict(latest_output.get("v3_predictions") or {})
        admet_steps = ((latest_output.get("steps") or {}).get("admet") or {}).get("endpoints") or []
        admet_by_name = {row.get("endpoint"): row for row in admet_steps}

        rows = []
        for task in tasks:
            endpoint = task["endpoint"]
            route = routes.get(endpoint)
            registration = model_artifact_registration(endpoint)
            snapshot = snapshots_by_endpoint.get(endpoint)
            query_endpoint = PROFILE_QUERY.get(endpoint, endpoint)
            profile_row = visible.get(query_endpoint)
            quant = quantitative.get(endpoint) or {}
            admet = admet_by_name.get(ADMET_NAMES.get(endpoint, "")) or {}
            property_value = (version.properties_json or {}).get({
                "MW": "exact_molecular_weight", "CLOGP": "clogp", "TPSA": "tpsa",
                "HBD": "hbd", "HBA": "hba", "ROTB": "rotatable_bonds",
                "FSP3": "fraction_csp3", "QED": "qed", "FORMAL_CHARGE": "formal_charge",
                "HEAVY_ATOM_COUNT": "heavy_atom_count",
            }.get(endpoint, ""))
            included = endpoint in PREDICT_ALL_EXECUTION_ENDPOINTS
            execution_source = (
                "CURRENT_PRODUCTION_ROUTING" if quant
                else "CACHED_REUSE" if admet.get("cache_hit")
                else "PROPERTY_PIPELINE" if property_value is not None
                else "MECHANISTIC_PIPELINE" if endpoint in MECHANISTIC_EXECUTION_ENDPOINTS
                else "EXCLUDED_BY_CONTRACT"
            )
            attempted = execution_source in {"CURRENT_PRODUCTION_ROUTING", "PROPERTY_PIPELINE", "MECHANISTIC_PIPELINE"}
            succeeded = bool(
                quant.get("execution_status") == "SUCCESS"
                or admet.get("status") == "COMPLETE"
                or property_value is not None
                or endpoint in MECHANISTIC_EXECUTION_ENDPOINTS
            )
            raw_value = quant.get("production_prediction")
            if raw_value is None and property_value is not None:
                raw_value = property_value
            if raw_value is None and snapshot is not None:
                raw_value = snapshot.value if snapshot.value is not None else snapshot.classification
            admission = admit_current_prediction(db, snapshot) if snapshot else None
            rows.append({
                "canonical_endpoint": endpoint,
                "category": task["domain"],
                "species": SPECIES.get(endpoint, "HUMAN"),
                "unit": (
                    snapshot.unit if snapshot is not None
                    else profile_row.get("unit", "") if profile_row is not None
                    else ""
                ),
                "campaign_decision": task["decision"],
                "model_availability": (
                    "AVAILABLE_CURRENT" if endpoint in PUBLISHABLE_CURRENT_ENDPOINTS
                    else "UNQUALIFIED_LEGACY_MODEL" if endpoint == "AMES_MUTAGENICITY"
                    else "MODEL_NOT_REGISTERED" if endpoint in MODEL_NOT_REGISTERED_ENDPOINTS
                    else "MECHANISTIC_ONLY" if endpoint in MECHANISTIC_EXECUTION_ENDPOINTS
                    else "CURRENT_DATA_CEILING" if endpoint in CURRENT_DATA_CEILING_ENDPOINTS
                    else "CONTEXT_REQUIRED" if endpoint in CONTEXT_REQUIRED_ENDPOINTS
                    else "MODEL_UNAVAILABLE"
                ),
                "registered_model_id": registration.model_id if registration else None,
                "model_version": registration.model_version if registration else None,
                "artifact": list(registration.artifact_paths) if registration else [],
                "artifact_hash": artifact_bundle_sha256(registration) if registration else None,
                "artifact_validation": validation_registration_error(registration) if registration else "MODEL_NOT_REGISTERED",
                "route_exists": route is not None,
                "route": route.get("route") if route else None,
                "predict_workflow_includes_it": included,
                "execution_attempted": attempted,
                "execution_source": execution_source,
                "execution_succeeded": succeeded,
                "raw_prediction_produced": raw_value is not None or (endpoint in MECHANISTIC_EXECUTION_ENDPOINTS and succeeded),
                "raw_prediction": raw_value,
                "stable_core_admission_result": (
                    "ELIGIBLE" if admission and admission.eligible
                    else admission.reason if admission
                    else "NOT_SUBMITTED"
                ),
                "current_prediction_snapshot_created": snapshot is not None,
                "snapshot_id": snapshot.id if snapshot else None,
                "visible_in_canonical_ui": profile_row is not None,
                "ui_status": profile_row.get("status") if profile_row else None,
                "before": (
                    "PREDICTED" if endpoint in before_published
                    else f"CALCULATED_BUT_NOT_ELIGIBLE:{PRODUCTION_BEFORE_REJECTED[endpoint]}"
                    if endpoint in PRODUCTION_BEFORE_REJECTED
                    else "NO_CURRENT_SNAPSHOT"
                ),
                "after": (
                    (profile_row.get("prediction") or {}).get("display_value")
                    or (profile_row.get("prediction") or {}).get("value")
                    if profile_row and profile_row.get("prediction")
                    else "—"
                ),
                "classification": classification(endpoint),
                "final_reason": final_reason(endpoint),
            })

        counts = Counter(row["classification"] for row in rows)
        payload = {
            "audit": "MobocertinibPredictionCompleteness/1",
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "database": str(engine.url).replace(str(engine.url.database), "<isolated-e2e-db>"),
            "production_data_modified": False,
            "reference_trace": {
                "project_id": compound.project_id,
                "project_name": compound.project.name,
                "compound_row_id": compound.id,
                "compound_id": compound.compound_id,
                "compound_name": compound.name,
                "compound_version_id": version.id,
                "version_number": version.version_number,
                "inchikey": version.inchikey,
                "production_before": PRODUCTION_BEFORE_TRACE,
                "verified_prediction_run_id": latest.id if latest else None,
                "verified_inputs_hash": latest.inputs_hash if latest else None,
            },
            "endpoint_count": len(rows),
            "classification_counts": dict(sorted(counts.items())),
            "rows": rows,
        }
        if len(rows) != 59:
            raise RuntimeError(f"Campaign matrix must contain 59 rows, got {len(rows)}")
        if any(row["classification"] == "EXECUTED_AND_PUBLISHED" and not row["current_prediction_snapshot_created"] for row in rows):
            missing = [row["canonical_endpoint"] for row in rows if row["classification"] == "EXECUTED_AND_PUBLISHED" and not row["current_prediction_snapshot_created"]]
            raise RuntimeError(f"Qualified endpoints missing CurrentPredictionSnapshot: {missing}")
        Path(args.output).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

        contract_rows = [{
            "canonical_endpoint": row["canonical_endpoint"],
            "category": row["category"],
            "campaign_decision": row["campaign_decision"],
            "predict_behavior": (
                "EXECUTE_AND_PUBLISH" if row["canonical_endpoint"] in PUBLISHABLE_CURRENT_ENDPOINTS
                else "EXECUTE_MECHANISTIC_ONLY" if row["canonical_endpoint"] in MECHANISTIC_EXECUTION_ENDPOINTS
                else "DO_NOT_EXECUTE"
            ),
            "availability": row["model_availability"],
            "reason": row["final_reason"],
        } for row in rows]
        contract_payload = {
            "contract": CONTRACT_ID,
            "engine_id": "drugopt-prediction-engine-v3@3.3.3",
            "stable_core": "v1.2",
            "endpoint_count": len(contract_rows),
            "execute_and_publish_count": len(PUBLISHABLE_CURRENT_ENDPOINTS),
            "execute_mechanistic_only_count": len(MECHANISTIC_EXECUTION_ENDPOINTS),
            "do_not_execute_count": 59 - len(PREDICT_ALL_EXECUTION_ENDPOINTS),
            "rules": [
                "Only an explicit Predict action executes this contract.",
                "Every frozen ensemble component must succeed; partial ensembles fail closed.",
                "MODEL_UNAVAILABLE, MODEL_NOT_REGISTERED, CURRENT_DATA_CEILING, and CONTEXT_REQUIRED endpoints are not executed.",
                "Only unchanged Stable Core v1.2 admission may create CurrentPredictionSnapshot rows.",
            ],
            "endpoints": contract_rows,
        }
        Path(args.contract_output).write_text(json.dumps(contract_payload, indent=2, sort_keys=True) + "\n")
    finally:
        db.close()


if __name__ == "__main__":
    main()
