"""Build a read-only project-learning readiness artifact from live projects."""
from __future__ import annotations

import hashlib
import json
import sqlite3
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from backend.project_learning_contract import adapter_capability, GLOBAL_ENGINE_HASH, GLOBAL_ENGINE_POLICY, PROJECT_LEARNING_VERSION


def main():
    con = sqlite3.connect(ROOT / "drug_opt.db")
    projects = list(con.execute("select id,name from projects where id in (1,3,5) order by id"))
    rows = []
    supported_names = {"Solubility", "Permeability", "Plasma protein binding", "HLM intrinsic clearance", "RLM intrinsic clearance", "MLM intrinsic clearance"}
    for project_id, name in projects:
        compound_n = con.execute("select count(distinct id) from compounds where project_id=?", (project_id,)).fetchone()[0]
        endpoints = []
        for endpoint, n in con.execute("""select e.name,count(distinct cv.compound_row_id)
            from compounds c join compound_versions cv on cv.compound_row_id=c.id
            join admet_measurements m on m.version_id=cv.id join admet_endpoints e on e.id=m.endpoint_id
            where c.project_id=? group by e.name""", (project_id,)):
            endpoints.append({"endpoint": endpoint, "experimental_compounds": n, "global_model": "existing core route" if endpoint in supported_names else "endpoint-specific/unknown", "adapter": adapter_capability(endpoint, global_status="PRODUCTION_STABLE"), "status": "INSUFFICIENT_REAL_N"})
        for endpoint, n in con.execute("""select e.canonical_endpoint_id,count(distinct cv.compound_row_id)
            from compounds c join compound_versions cv on cv.compound_row_id=c.id
            join external_experimental_evidence e on e.compound_version_id=cv.id
            where c.project_id=? group by e.canonical_endpoint_id""", (project_id,)):
            if n:
                endpoints.append({"endpoint": endpoint, "experimental_compounds": n, "global_model": "canonical route if supported", "adapter": adapter_capability(endpoint, global_status="PRODUCTION_STABLE"), "status": "INSUFFICIENT_REAL_N"})
        rows.append({"project_id": project_id, "project": name, "compound_n": compound_n, "endpoint_readiness": endpoints, "trained_real_adapter_n": 0, "reason": "No project/endpoint has sufficient validated LOCO improvement evidence; no real adapter was trained."})
    registry_n = con.execute("select count(*) from project_adapter_versions").fetchone()[0]
    baseline = json.loads((ROOT / "validation/global_developability_v1_baseline.json").read_text())
    artifact = {
        "artifact": "project_learning_foundation_v7_0",
        "created_at": str(date.today()), "policy_version": PROJECT_LEARNING_VERSION,
        "global_engine": {"policy": GLOBAL_ENGINE_POLICY, "hash": GLOBAL_ENGINE_HASH, "status": "UNCHANGED"},
        "projects": rows,
        "project_adapter_registry": {"table": "project_adapter_versions", "existing_rows": registry_n, "historical_immutability": True},
        "real_pilot": {"trained": False, "status": "NO_SCIENTIFIC_ADAPTER_TRAINED", "reason": "Real projects do not currently meet leakage-safe validation gate."},
        "synthetic_fixture_policy": "TEST_ONLY; never inserted into project measurements, evidence warehouse, or runtime DB",
        "global_baseline_sha256": hashlib.sha256(json.dumps(baseline, sort_keys=True).encode()).hexdigest(),
        "adapter_contract": ["global_prediction", "global_uncertainty", "GLOBAL_AD", "PROJECT_AD", "model_version", "model_status", "recommended_source", "residual"],
        "recommendation": "ACCUMULATE_INTERNAL_EXPERIMENTAL_DATA_AND_ACTIVATE_PROJECT_ADAPTERS",
    }
    out = ROOT / "validation/project_learning_foundation_v7_0.json"
    out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"projects": [(x["project_id"], x["compound_n"]) for x in rows], "registry_rows": registry_n, "real_adapter_trained": False}, indent=2))


if __name__ == "__main__":
    main()
