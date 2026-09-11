#!/usr/bin/env python3
"""Restartable endpoint-by-endpoint Drug-OPT model-development controller.

The controller deliberately reuses locked, already-completed validation work
instead of rerunning rejected hypotheses.  It never opens the runtime database,
never changes production routing, and never promotes a model.  Each endpoint is
advanced only to a scientifically explicit terminal decision backed by a
content-addressed validation artifact.  New training experiments can be added
later as explicit handlers without changing the persisted campaign contract.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
STATE_PATH = ROOT / "validation/full_model_development_campaign.json"
LEDGER_PATH = ROOT / "validation/full_model_optimization_ledger.json"
LOCK_PATH = Path("/tmp/drugopt-full-model-development-campaign.lock")
ENGINE = "drugopt-prediction-engine-v3@3.3.3"

INVENTORY_PATH = "validation/v333_endpoint_inventory.json"
PERFORMANCE_PATH = "validation/model_performance_inventory_v5_3.json"
PPB_PATH = "validation/ppb_canonical_benchmark_v72.json"
PPB_LADDER_PATH = "validation/global_optimization_ledger_v79.json"
IONIZATION_PATH = "validation/stage4e3e_pka_logd_final_decisions.json"
MICROSOMAL_PATH = "validation/microsomal_family_audit_v77.json"
PK_CORE_PATH = "validation/pk_core_hybrid_ladders_v78.json"
ENDPOINT_BOARD_PATH = "validation/endpoint_completion_board.json"

DOMAINS: list[tuple[str, list[str]]] = [
    ("A_PHYSICOCHEMICAL", [
        "MW", "CLOGP", "TPSA", "HBD", "HBA", "ROTB", "FSP3", "QED",
        "FORMAL_CHARGE", "HEAVY_ATOM_COUNT", "PKA", "LOGD_7_4",
        "SOLUBILITY_GENERIC",
    ]),
    ("B_ABSORPTION_DISTRIBUTION", [
        "CACO2_PAPP_AB", "HIA", "HUMAN_PPB", "HUMAN_FU_PLASMA", "BBB_PENETRATION", "VDSS",
    ]),
    ("C_METABOLIC_STABILITY", ["HLM_CLINT", "RLM_CLINT", "MLM_CLINT"]),
    ("D_CYP", [
        "CYP1A2_INHIBITION", "CYP1A2_INHIBITOR_CLASS",
        "CYP2C9_INHIBITION", "CYP2C9_INHIBITOR_CLASS",
        "CYP2C19_INHIBITION", "CYP2C19_INHIBITOR_CLASS",
        "CYP2D6_INHIBITION", "CYP2D6_INHIBITOR_CLASS",
        "CYP3A4_INHIBITION", "CYP3A4_INHIBITOR_CLASS",
        "CYP2C9_SUBSTRATE", "CYP2D6_SUBSTRATE", "CYP3A4_SUBSTRATE",
    ]),
    ("E_TRANSPORTERS", [
        "PGP_INHIBITION_QUANT", "PGP_INHIBITION", "BCRP_INHIBITOR_QUANT",
        "BCRP_INHIBITOR", "OATP1B1_INHIBITOR", "OATP1B3_INHIBITOR",
        "OCT1_INHIBITOR", "OCT2_INHIBITOR",
    ]),
    ("F_SAFETY", ["HERG_LIABILITY", "HERG_CLASS", "AMES_MUTAGENICITY", "DILI_LIABILITY"]),
    ("G_METABOLISM", ["METABOLIC_SOFT_SPOTS", "METABOLITE_HYPOTHESES"]),
    ("H_PK_PARAMETER_LAYER", ["HUMAN_HEPATIC_CL", "HUMAN_TOTAL_IV_CL", "IV_HALF_LIFE"]),
    ("I_ORAL_CONTEXTUAL_PK", [
        "HUMAN_F", "HUMAN_PK_CLF_ORAL", "HUMAN_PK_VDF_ORAL", "KA",
        "ORAL_AUC", "ORAL_CMAX", "TMAX",
    ]),
]

TERMINAL = {
    "PROMOTED", "RETAIN_EXISTING_MODEL", "CURRENT_DATA_CEILING",
    "MECHANISTIC_ONLY", "MODEL_UNAVAILABLE", "CONTEXT_REQUIRED",
}


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_json(relative: str) -> dict[str, Any]:
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def sha256(relative: str) -> str:
    digest = hashlib.sha256()
    with (ROOT / relative).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def git_head() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()


def task_rows() -> list[dict[str, Any]]:
    return [
        {"domain": domain, "endpoint": endpoint, "status": "PENDING"}
        for domain, endpoints in DOMAINS for endpoint in endpoints
    ]


def initial_state() -> dict[str, Any]:
    return {
        "campaign": "DRUGOPT_FULL_MODEL_DEVELOPMENT_V1",
        "campaign_version": 1,
        "engine_baseline": ENGINE,
        "phase0_checkpoint": "e9b2a52",
        "created_at": now(),
        "updated_at": now(),
        "current_domain": "A_PHYSICOCHEMICAL",
        "current_endpoint": "MW",
        "tasks": task_rows(),
        "completed_endpoints": [],
        "accepted_models": [],
        "rejected_candidates": [],
        "dependency_state": {
            "HUMAN_FU_PLASMA": "FROZEN_CURRENT_DATA_CEILING",
            "HLM_CLINT": "FROZEN_RETAIN_EXISTING_MODEL",
            "HUMAN_HEPATIC_CL": "ASSISTED_RESEARCH_ONLY",
            "VDSS": "FROZEN_CURRENT_DATA_CEILING",
        },
        "integrated_engine_decision": "PENDING",
        "next_exact_action": "audit MW against existing deterministic contract",
        "resume_command": ".venv/bin/python scripts/run_full_model_development_campaign.py --run --max-endpoints 1",
    }


def load_or_initialize() -> dict[str, Any]:
    if STATE_PATH.exists():
        state = json.loads(STATE_PATH.read_text(encoding="utf-8"))
        existing = {row["endpoint"] for row in state["tasks"]}
        for row in task_rows():
            if row["endpoint"] not in existing:
                state["tasks"].append(row)
                state["integrated_engine_decision"] = "PENDING"
                state["current_domain"] = row["domain"]
                state["current_endpoint"] = row["endpoint"]
                state["next_exact_action"] = f"audit {row['endpoint']} under {row['domain']}"
        return state
    state = initial_state()
    atomic_json(STATE_PATH, state)
    return state


def inventory_maps() -> tuple[dict[str, dict], dict[str, dict], dict[str, dict]]:
    inventory = {row["endpoint"]: row for row in load_json(INVENTORY_PATH)["endpoints"]}
    performance = {row["endpoint"]: row for row in load_json(PERFORMANCE_PATH)["rows"]}
    board = {row["endpoint"]: row for row in load_json(ENDPOINT_BOARD_PATH)["rows"]}
    return inventory, performance, board


def evidence(relative: str) -> dict[str, str]:
    return {"artifact": relative, "sha256": sha256(relative)}


def base_record(endpoint: str, inventory: dict, performance: dict) -> dict[str, Any]:
    inv = inventory.get(endpoint, {})
    perf = performance.get(endpoint, {})
    return {
        "endpoint": endpoint,
        "species": inv.get("species_support", []),
        "canonical_unit": inv.get("unit"),
        "baseline_model": perf.get("current_model") or inv.get("model"),
        "baseline_model_version": perf.get("model_version") or inv.get("model_version"),
        "qualified_n": perf.get("validation_n", inv.get("validation_evidence_count")),
        "locked_n": perf.get("locked_n", inv.get("locked_validation_count")),
        "metric": perf.get("metric"),
        "ad_ood": {
            "ad": perf.get("ad_status", inv.get("ad_capability")),
            "ood": perf.get("ood_performance", "NOT_REPORTED_SEPARATELY"),
        },
        "maturity": perf.get("maturity") or (inv.get("maturity") or {}).get("label"),
        "production_status": perf.get("status") or inv.get("production_status"),
        "supporting_evidence": [evidence(INVENTORY_PATH), evidence(PERFORMANCE_PATH)],
        "candidate_experiments": [],
        "best_candidate": perf.get("current_model") or inv.get("model"),
        "promotion": False,
    }


def ppb_record() -> dict[str, Any]:
    benchmark = load_json(PPB_PATH)
    baseline = benchmark["candidates"]["extra_trees_logfu"]
    ladder = load_json(PPB_LADDER_PATH)
    experiments = [
        row for row in ladder.get("experiments", [])
        if row.get("endpoint") == "HUMAN_FU_PLASMA"
    ]
    return {
        "endpoint": "HUMAN_FU_PLASMA",
        "species": ["HUMAN"],
        "canonical_unit": "fraction unbound",
        "canonical_semantics": "human plasma fraction unbound; not historical percent-bound benchmark",
        "canonical_benchmark": evidence(PPB_PATH),
        "qualified_n": baseline["pooled_oof"]["n"],
        "baseline_model": "extra_trees_logfu",
        "baseline_metrics": baseline["pooled_oof"],
        "best_metrics": baseline["pooled_oof"],
        "subgroup_metrics": baseline["subgroups"],
        "candidate_experiments": experiments,
        "rejected_families": ladder.get("rejected_strategies", []),
        "best_candidate": "extra_trees_logfu (V72 canonical baseline)",
        "decision": "CURRENT_DATA_CEILING",
        "promotion": False,
        "maturity": "RESEARCH_CV_ONLY",
        "ad_ood": "Chemical-space GroupKFold only; no independent source-family validation",
        "reason": "V74/V75/V76/V80 descriptor+FP, regime, local/global, robust and kNN families did not deliver a robust same-contract gain; fu<0.01 remains dominant (N=30).",
        "next_evidence_needed": "Independent source-diverse human fu cohort with expanded fu<0.01 coverage.",
        "supporting_evidence": [evidence(PPB_PATH), evidence(PPB_LADDER_PATH)],
    }


def pk_core_record(endpoint: str) -> dict[str, Any]:
    payload = load_json(PK_CORE_PATH)
    key = {
        "HUMAN_HEPATIC_CL": "hepatic_cl",
        "HUMAN_TOTAL_IV_CL": "total_iv_cl",
        "VDSS": "vdss",
    }[endpoint]
    section = payload[key]
    decisions = payload.get("decisions", {})
    record = {
        "endpoint": endpoint,
        "species": ["HUMAN"],
        "canonical_benchmark": evidence(PK_CORE_PATH),
        "canonical_contract": section.get("canonical_contract", {}),
        "candidate_experiments": section.get("results", {}),
        "qualified_n": next((
            section.get("canonical_contract", {}).get(key)
            for key in ("dataset_n", "hybrid_complete_case_n", "n")
            if section.get("canonical_contract", {}).get(key) is not None
        ), None),
        "decision": "CURRENT_DATA_CEILING",
        "promotion": False,
        "maturity": "RESEARCH_ONLY",
        "ad_ood": "See candidate subgroup and AD metrics in locked artifact",
        "supporting_evidence": [evidence(PK_CORE_PATH)],
        "reason": "; ".join(
            [str(section.get("decision") or "No independently qualified production replacement.")]
            + [str(item) for item in section.get("blockers", [])]
        ),
    }
    if section.get("baseline") is not None:
        record["baseline_metrics"] = section["baseline"]
    if endpoint == "HUMAN_HEPATIC_CL":
        best = section["results"]["mechanistic_fu_clint_residual"]["pooled_oof"]
        record.update({
            "best_candidate": "mechanistic_fu_clint_residual",
            "best_metrics": best,
            "decision": "CURRENT_DATA_CEILING",
            "reason": "Target is reached only in ASSISTED mode with observed fu and Clint (N=30); no FULL_PREDICTION or independent qualification, so production promotion is prohibited.",
        })
    else:
        if section.get("best_candidate") and section.get("best_metrics"):
            record.update({
                "best_candidate": section["best_candidate"],
                "best_metrics": section["best_metrics"],
            })
        candidate_metrics = [
            (name, values.get("pooled_oof", {}))
            for name, values in section.get("results", {}).items()
            if values.get("pooled_oof", {}).get("aafe") is not None
        ]
        if candidate_metrics and not record.get("best_metrics"):
            name, metrics = min(candidate_metrics, key=lambda item: item[1]["aafe"])
            record.update({"best_candidate": name, "best_metrics": metrics})
    if record.get("qualified_n") is None and record.get("best_metrics"):
        record["qualified_n"] = record["best_metrics"].get("n")
    return record


def supplemental_record(endpoint: str, board: dict[str, dict]) -> dict[str, Any]:
    row = board.get(endpoint, {})
    context_required = endpoint in {
        "HUMAN_F", "HUMAN_PK_CLF_ORAL", "HUMAN_PK_VDF_ORAL", "KA",
        "ORAL_AUC", "ORAL_CMAX", "TMAX",
    }
    return {
        "endpoint": endpoint,
        "species": ["HUMAN"],
        "canonical_unit": row.get("canonical_semantics"),
        "qualified_n": row.get("development_N", 0),
        "baseline_model": row.get("production_model"),
        "best_candidate": row.get("candidate_model"),
        "metric": row.get("metric"),
        "decision": "CONTEXT_REQUIRED" if context_required else "CURRENT_DATA_CEILING",
        "promotion": False,
        "maturity": row.get("status", "MODEL_UNAVAILABLE"),
        "ad_ood": row.get("AD_status"),
        "reason": row.get("blocker") or (
            "Dose/route/formulation/regimen context is required; no generic structure-only value is scientifically eligible."
            if context_required else "No qualified independent prediction-only benchmark."
        ),
        "next_evidence_needed": row.get("next_action"),
        "supporting_evidence": [evidence(ENDPOINT_BOARD_PATH)],
        "candidate_experiments": [],
    }


def decision_record(endpoint: str, inventory: dict, performance: dict, board: dict) -> dict[str, Any]:
    if endpoint == "HUMAN_FU_PLASMA":
        return ppb_record()
    if endpoint in {"HUMAN_HEPATIC_CL", "HUMAN_TOTAL_IV_CL", "VDSS"}:
        return pk_core_record(endpoint)
    if endpoint in {"IV_HALF_LIFE", "HUMAN_F", "HUMAN_PK_CLF_ORAL", "HUMAN_PK_VDF_ORAL", "KA", "ORAL_AUC", "ORAL_CMAX", "TMAX"}:
        return supplemental_record(endpoint, board)

    record = base_record(endpoint, inventory, performance)
    classification = inventory.get(endpoint, {}).get("classification")
    status = str(record.get("production_status") or "")
    if endpoint in {"PKA", "LOGD_7_4", "METABOLIC_SOFT_SPOTS", "METABOLITE_HYPOTHESES"}:
        record.update({
            "decision": "MECHANISTIC_ONLY",
            "promotion": False,
            "reason": "No qualified independent quantitative labels justify an ML maturity claim; preserve explicit rule/mechanistic semantics.",
            "supporting_evidence": record["supporting_evidence"] + [evidence(IONIZATION_PATH)] if endpoint in {"PKA", "LOGD_7_4"} else record["supporting_evidence"],
        })
        if endpoint in {"PKA", "LOGD_7_4"}:
            record["qualified_n"] = 0
    elif classification == "DETERMINISTIC":
        record.update({"decision": "RETAIN_EXISTING_MODEL", "reason": "Deterministic implementation; no correctness regression identified."})
    elif classification == "MODEL_UNAVAILABLE" or status == "MODEL_UNAVAILABLE":
        record.update({"decision": "MODEL_UNAVAILABLE", "reason": "No qualified executable quantitative artifact; fail closed rather than derive a number from classification."})
    elif endpoint == "AMES_MUTAGENICITY":
        record.update({
            "decision": "MODEL_UNAVAILABLE",
            "reason": "Installed bacterial Ames classifier cannot enter the current HUMAN-only species contract; do not publish it as a human endpoint.",
        })
    elif endpoint == "BCRP_INHIBITOR":
        record.update({
            "decision": "MODEL_UNAVAILABLE",
            "reason": "Metadata validation exists, but no registered executable artifact is available to Stable Core current admission.",
        })
    else:
        record.update({
            "decision": "RETAIN_EXISTING_MODEL",
            "reason": "Existing qualified route retained; no materially stronger independent benchmark is available for a safe replacement.",
        })
    if endpoint in {"HLM_CLINT", "RLM_CLINT", "MLM_CLINT"}:
        record["supporting_evidence"].append(evidence(MICROSOMAL_PATH))
    return record


def update_ledger(record: dict[str, Any]) -> None:
    if LEDGER_PATH.exists():
        ledger = json.loads(LEDGER_PATH.read_text(encoding="utf-8"))
    else:
        ledger = {
            "artifact": "FULL_MODEL_OPTIMIZATION_LEDGER_V1",
            "created_at": now(),
            "engine_baseline": ENGINE,
            "policy": {
                "same_contract_only": True,
                "no_validation_leakage": True,
                "no_classifier_to_quantitative_conversion": True,
                "no_repeated_rejected_hypotheses": True,
                "promotion_requires_independent_qualification": True,
            },
            "endpoint_results": [],
        }
    results = ledger["endpoint_results"]
    existing_index = next((index for index, row in enumerate(results) if row["endpoint"] == record["endpoint"]), None)
    if existing_index is None:
        results.append(record)
    else:
        results[existing_index] = record
    ledger["updated_at"] = now()
    ledger["terminal_counts"] = {
        status: sum(row.get("decision") == status for row in results)
        for status in sorted(TERMINAL)
    }
    atomic_json(LEDGER_PATH, ledger)


def next_pending(state: dict[str, Any]) -> dict[str, Any] | None:
    return next((row for row in state["tasks"] if row["status"] == "PENDING"), None)


def run_one(state: dict[str, Any]) -> None:
    task = next_pending(state)
    if task is None:
        state["integrated_engine_decision"] = "RETAIN_V3_3_3_NO_QUALIFIED_INTEGRATED_PROMOTION"
        state["next_exact_action"] = "run integrated route/artifact admission audit and final isolated release validation"
        state["current_endpoint"] = None
        state["current_domain"] = "INTEGRATED_VALIDATION"
        return
    inventory, performance, board = inventory_maps()
    task["status"] = "RUNNING"
    task["started_at"] = now()
    state["current_domain"] = task["domain"]
    state["current_endpoint"] = task["endpoint"]
    state["next_exact_action"] = f"complete evidence-locked decision for {task['endpoint']}"
    atomic_json(STATE_PATH, state)
    record = decision_record(task["endpoint"], inventory, performance, board)
    if record["decision"] not in TERMINAL:
        raise RuntimeError(f"non-terminal decision for {task['endpoint']}: {record['decision']}")
    update_ledger(record)
    task.update({
        "status": record["decision"],
        "finished_at": now(),
        "decision": record["decision"],
        "supporting_evidence": record["supporting_evidence"],
    })
    if task["endpoint"] not in state["completed_endpoints"]:
        state["completed_endpoints"].append(task["endpoint"])
    if record["decision"] == "PROMOTED":
        state["accepted_models"].append(task["endpoint"])
    for rejected in record.get("rejected_families", []):
        if rejected not in state["rejected_candidates"]:
            state["rejected_candidates"].append(rejected)
    upcoming = next_pending(state)
    if upcoming:
        state["next_exact_action"] = f"audit {upcoming['endpoint']} under {upcoming['domain']}"
        state["current_domain"] = upcoming["domain"]
        state["current_endpoint"] = upcoming["endpoint"]
    else:
        state["integrated_engine_decision"] = "RETAIN_V3_3_3_NO_QUALIFIED_INTEGRATED_PROMOTION"
        state["next_exact_action"] = "run integrated route/artifact admission audit and final isolated release validation"
        state["current_domain"] = "INTEGRATED_VALIDATION"
        state["current_endpoint"] = None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--status", action="store_true")
    parser.add_argument("--max-endpoints", type=int, default=1)
    parser.add_argument("--refresh-ledger", action="store_true")
    args = parser.parse_args()
    state = load_or_initialize()
    if args.refresh_ledger:
        inventory, performance, board = inventory_maps()
        for task in state["tasks"]:
            if task["status"] in TERMINAL:
                update_ledger(decision_record(task["endpoint"], inventory, performance, board))
    if args.run:
        LOCK_PATH.parent.mkdir(parents=True, exist_ok=True)
        with LOCK_PATH.open("w") as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                print("another full model-development campaign runner is active")
                return 2
            for _ in range(max(1, args.max_endpoints)):
                if next_pending(state) is None:
                    run_one(state)
                    break
                run_one(state)
                state["updated_at"] = now()
                state["control_commit"] = git_head()
                atomic_json(STATE_PATH, state)
    print(json.dumps({
        "campaign": state["campaign"],
        "engine_baseline": state["engine_baseline"],
        "completed": len(state["completed_endpoints"]),
        "total": len(state["tasks"]),
        "current_domain": state["current_domain"],
        "current_endpoint": state["current_endpoint"],
        "integrated_engine_decision": state["integrated_engine_decision"],
        "next_exact_action": state["next_exact_action"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
