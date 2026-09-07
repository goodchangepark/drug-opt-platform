"""Create the v5.3 model-health inventory without changing production models."""
from __future__ import annotations

import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.prediction_engine_registry import (  # noqa: E402
    CURRENT_ENGINE_ID, CURRENT_ENGINE_VERSION, CURRENT_POLICY_HASH,
    CANDIDATE_ENGINE_ID, CANDIDATE_ENGINE_VERSION, CANDIDATE_POLICY_HASH,
    get_current_production_routing,
)
from backend.prediction_maturity import get_endpoint_maturity_registry  # noqa: E402

OUT_BASELINE = ROOT / "validation/model_stabilization_baseline_v5_3.json"
OUT_INVENTORY = ROOT / "validation/model_performance_inventory_v5_3.json"
OUT_GAPS = ROOT / "validation/information_gap_map_v5_3.json"

REPORT = ROOT / "validation/expanded_pk_validation_report_v1.json"


def _read(name, default=None):
    path = ROOT / "validation" / name
    try:
        return json.loads(path.read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return default


def _triage(endpoint_id: str, route: str, maturity: int, unavailable: bool) -> tuple[str, str]:
    if endpoint_id in {"HUMAN_PK_CLF_ORAL", "HUMAN_PK_VDF_ORAL"}:
        return "CONTEXT_LIMITED", "Keep route/dose/regimen-specific PK and fail closed when F or context is absent."
    if endpoint_id in {"PK_AUC", "PK_CMAX", "PK_TMAX", "PK_HALF_LIFE", "HUMAN_HEPATIC_CL", "HUMAN_TOTAL_CL"}:
        return "MODEL_LIMITED", "Retain current mechanistic route; require independent endpoint validation before promotion."
    if endpoint_id in {"VDSS"}:
        return "IMPROVEMENT_CANDIDATE", "Evaluate only on a new leakage-safe independent cohort; do not use observed VDss as a scored prediction."
    if endpoint_id in {"PKA", "LOGD_7_4"}:
        return "DATA_LIMITED", "Confirm Level-1 rule/derived status; stop repeated model hunting until a qualified cohort/checkpoint exists."
    if endpoint_id in {"CYP2C19_INHIBITION", "PGP_INHIBITION_QUANT", "BCRP_INHIBITOR_QUANT", "OATP1B1_INHIBITOR", "OATP1B3_INHIBITOR", "OCT1_INHIBITOR", "OCT2_INHIBITOR"}:
        return "MODEL_UNAVAILABLE", "Remain fail closed; no executable, independently validated quantitative production model."
    if unavailable:
        return "MODEL_UNAVAILABLE", "Remain explicitly unavailable."
    if maturity >= 4:
        return "STABLE_PRODUCTION", "Regression monitoring only; no retraining without a stronger independent holdout."
    if maturity == 3:
        return "IMPROVEMENT_CANDIDATE", "Review only with new independent evidence and locked-cohort protection."
    return "DATA_LIMITED", "No promotion without endpoint-specific independent validation."


def main() -> None:
    now = datetime.now(timezone.utc).isoformat()
    maturity_rows = get_endpoint_maturity_registry()
    routing = get_current_production_routing()
    route_map = {row.get("endpoint_id"): row for row in routing}
    report = _read("expanded_pk_validation_report_v1.json", {})
    pk_metrics = report.get("pk_parameter_performance", {})
    db = sqlite3.connect(ROOT / "drug_opt.db")
    try:
        protected = db.execute("SELECT id,name FROM projects WHERE id IN (1,3,5,300) ORDER BY id").fetchall()
        drugbank_n = db.execute("SELECT COUNT(*) FROM compounds WHERE project_id=300").fetchone()[0]
        prediction_runs = db.execute("SELECT COUNT(*) FROM prediction_runs").fetchone()[0]
    finally:
        db.close()

    inventory = []
    for row in maturity_rows:
        endpoint_id = row["endpoint_id"]
        triage, action = _triage(endpoint_id, row.get("model_route", ""), row.get("maturity_level", 1), row.get("is_unavailable", False))
        inventory.append({
            "endpoint": endpoint_id,
            "display_name": row.get("endpoint_name"),
            "section": "PK" if endpoint_id.startswith(("HUMAN_PK", "VDSS")) else "ADMET",
            "current_model": row.get("model_route"),
            "model_version": row.get("engine_version"),
            "maturity": row.get("maturity_label"),
            "maturity_level": row.get("maturity_level"),
            "experimental_n": row.get("validation_n", 0),
            "development_n": None,
            "validation_n": row.get("validation_n", 0),
            "locked_n": row.get("locked_test_n", 0),
            "metric": row.get("maturity_reason", "Registry metric not specified"),
            "error": None,
            "bias": None,
            "ad_status": row.get("ad_status", "UNKNOWN"),
            "ood_performance": "NOT_REPORTED_SEPARATELY",
            "triage": triage,
            "recommended_action": action,
            "status": row.get("status"),
        })

    custom = [
        ("HUMAN_HEPATIC_CL", "Human hepatic clearance", "HEPATIC_IVIVE_APPARENT", "CL AAFE 3.295; N=30", "MODEL_LIMITED"),
        ("HUMAN_TOTAL_CL", "Human total systemic clearance", "NO_PRODUCTION_MODEL", "Independent predicted total CL = 0/30", "SEMANTICS_LIMITED"),
        ("HUMAN_RENAL_CL", "Human renal clearance readiness", "ROUTING_ONLY", "15/30 renal CL observations; no universal renal model", "DATA_LIMITED"),
        ("PK_AUC", "Oral/clinical AUC", "MECHANISTIC_1COMPARTMENT", "AAFE 3.282; N=30", "MODEL_LIMITED"),
        ("PK_CMAX", "Oral/clinical Cmax", "MECHANISTIC_1COMPARTMENT", "AAFE 2.540; N=30", "MODEL_LIMITED"),
        ("PK_TMAX", "Oral/clinical Tmax", "MECHANISTIC_1COMPARTMENT", "Cohort endpoint; no promotion claim", "MODEL_LIMITED"),
        ("PK_HALF_LIFE", "Clinical half-life", "MECHANISTIC_1COMPARTMENT", "AAFE 2.619; N=30", "MODEL_LIMITED"),
    ]
    for endpoint, name, model, metric, triage in custom:
        inventory.append({
            "endpoint": endpoint, "display_name": name, "section": "PK", "current_model": model,
            "model_version": "drugopt-pk-engine-v1", "maturity": "Base / Mechanistic Estimate", "maturity_level": 1,
            "experimental_n": 30 if endpoint != "HUMAN_RENAL_CL" else 15, "development_n": None,
            "validation_n": 30 if endpoint != "HUMAN_RENAL_CL" else 15, "locked_n": 30,
            "metric": metric, "error": metric, "bias": None, "ad_status": "NOT_CALIBRATED",
            "ood_performance": "NOT_REPORTED_SEPARATELY", "triage": triage,
            "recommended_action": "Retain fail-closed semantics; improve only with independent development/validation data.",
            "status": "LIMITED_BENCHMARK",
        })

    baseline = {
        "artifact": "model_stabilization_baseline_v5_3",
        "created_at": now,
        "production_engine": {"id": CURRENT_ENGINE_ID, "version": CURRENT_ENGINE_VERSION, "policy_hash": CURRENT_POLICY_HASH},
        "candidate_engine": {"id": CANDIDATE_ENGINE_ID, "version": CANDIDATE_ENGINE_VERSION, "policy_hash": CANDIDATE_POLICY_HASH},
        "frozen_prediction_engine_v1": {"id": "drugopt-prediction-engine-v1@1.0.0", "hash": "12757ab197b5a70d8ea1754678d9a342ab0b6ea0d82f2896bebb767d686bbdeb"},
        "pk_layer": {"id": "drugopt-pk-engine-v1", "clearance_architecture": "drugopt-clearance-architecture-v1"},
        "production_routes": routing,
        "endpoint_maturity_registry_sha256": hashlib.sha256(json.dumps(maturity_rows, sort_keys=True).encode()).hexdigest(),
        "validation_cohorts": {"clinical_pk": {"n": 30, "artifact": "validation/expanded_pk_validation_report_v1.json", "verdict": report.get("final_scientific_verdict")}},
        "metrics": pk_metrics,
        "protected_projects": protected,
        "drugbank_compounds": drugbank_n,
        "prediction_runs": prediction_runs,
        "production_unchanged": True,
    }
    gaps = {
        "artifact": "information_gap_map_v5_3", "created_at": now,
        "decision": "DO_NOT_EXPAND_YET",
        "reason": "Current weakness is primarily semantic/component completeness and lack of qualified human clearance checkpoints, not an unqualified need for 50 generic compounds.",
        "gaps": [
            {"endpoint": "HUMAN_TOTAL_CL", "missing": ["independent renal/other component measurements", "route-matched IV total CL", "validated executable human total-CL model"], "useful_additional_compounds": "20-30 with IV CL and component provenance", "chemical_space": "high/low extraction, high PPB, renal-dominant and mixed elimination"},
            {"endpoint": "VDSS", "missing": ["independent predicted-vs-observed evaluation without observed VDss input"], "useful_additional_compounds": "20+ diverse IV volume studies", "chemical_space": "high VDss bases, high PPB acids, ampholytes"},
            {"endpoint": "PK oral AUC/Cmax/Tmax", "missing": ["dose/formulation/food/regimen-complete cohorts", "validated F and ka provenance"], "useful_additional_compounds": "30 context-complete oral PK records", "chemical_space": "absorption-limited and permeability/solubility-diverse drugs"},
            {"endpoint": "pKa/logD7.4", "missing": ["site-resolved continuous labels and deployable checkpoint"], "useful_additional_compounds": "not actionable until a qualified model/data source exists", "chemical_space": "acid/base/macro/micro ionization diversity"},
            {"endpoint": "CYP2C19/P-gp/BCRP quantitative", "missing": ["validated quantitative checkpoint and exact assay semantics"], "useful_additional_compounds": "endpoint-specific, not generic DrugBank expansion", "chemical_space": "substrate/inhibitor mechanism and concentration-complete assays"},
        ],
        "next_expansion_if_revisited": "Targeted 50 compounds stratified by clearance route, extraction, PPB, VDss, ionization, logD and clinical PK richness; no generic count expansion.",
    }
    OUT_BASELINE.write_text(json.dumps(baseline, indent=2) + "\n")
    OUT_INVENTORY.write_text(json.dumps({"artifact": "model_performance_inventory_v5_3", "created_at": now, "rows": inventory}, indent=2) + "\n")
    OUT_GAPS.write_text(json.dumps(gaps, indent=2) + "\n")
    print(json.dumps({"inventory_rows": len(inventory), "protected_projects": protected, "drugbank": drugbank_n, "prediction_runs": prediction_runs, "outputs": [str(OUT_BASELINE), str(OUT_INVENTORY), str(OUT_GAPS)]}, indent=2))


if __name__ == "__main__":
    main()
