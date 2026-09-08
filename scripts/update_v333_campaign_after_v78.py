#!/usr/bin/env python3
"""Persist the v3.3.3 campaign control plane after the V75-V78 batch."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def main() -> None:
    board_path = ROOT / "validation/endpoint_completion_board.json"
    board = json.loads(board_path.read_text())
    board["production_engine_version"] = "drugopt-prediction-engine-v3@3.3.3"
    board["control_plane_version"] = "V78"
    board["last_updated"] = now()
    for row in board["rows"]:
        endpoint = row["endpoint"]
        if endpoint == "HUMAN_FU_PLASMA":
            row.update({
                "status": "CURRENT_DATA_CEILING",
                "current_value": 4.760840445074518,
                "candidate_model": "extra_trees_logfu (V72 canonical baseline)",
                "last_experiment": "validation/ppb_weighted_ladder_v76.json",
                "next_action": "Acquire independent/source-diverse human fu validation; revisit local/regime strategies only after new evidence.",
                "blocker": "Canonical N=476 remains dominated by fu<0.01 errors; no independent source-family validation.",
            })
        elif endpoint == "HUMAN_HEPATIC_CL":
            row.update({
                "status": "RESEARCH_CANDIDATE",
                "current_value": 2.097108825628647,
                "candidate_model": "mechanistic_fu_clint_residual (V78 assisted)",
                "last_experiment": "validation/pk_core_hybrid_ladders_v78.json",
                "development_N": 30,
                "next_action": "Obtain independent prediction-only hepatic-clearance cohort; do not route assisted candidate.",
                "blocker": "V78 is assisted (observed fu + HLM Clint), N=30; independent validation absent.",
            })
        elif endpoint == "HUMAN_VDSS":
            row.update({
                "status": "CURRENT_DATA_CEILING",
                "current_value": 3.7057366603955706,
                "candidate_model": "direct_structure (V78 research)",
                "last_experiment": "validation/pk_core_hybrid_ladders_v78.json",
                "next_action": "Acquire independent record-level true Vss/Vdss context; preserve mechanistic production route.",
                "blocker": "Strict N=50 single-source cohort; no independent true-Vss validation.",
            })
        elif endpoint == "HUMAN_TOTAL_IV_CL":
            row.update({
                "status": "CURRENT_DATA_CEILING",
                "current_value": 7.341110193181601,
                "candidate_model": "random_forest (V78 within-source research)",
                "last_experiment": "validation/pk_core_hybrid_ladders_v78.json",
                "next_action": "Acquire record-level independent human IV total-CL cohort; do not tune consumed stress cohort.",
                "blocker": "Compact source markers lack record-level IV context and independent source-family validation.",
            })
        elif endpoint == "HUMAN_HEPATOCYTE_CLINT":
            row["next_action"] = "Acquire qualified hepatocyte-Clint labels; multitask transfer remains blocked without endpoint-specific gold."
        elif endpoint in {"IV_HALF_LIFE", "ORAL_AUC", "ORAL_CMAX"}:
            row["next_action"] = "Requeue after species-matched CL/Vdss/F/ka snapshots; context is required and no generic value is emitted."
    board_path.write_text(json.dumps(board, indent=2, sort_keys=True) + "\n")

    ledger = {
        "artifact": "GLOBAL_OPTIMIZATION_LEDGER_V79",
        "created_at": now(),
        "production_engine": "drugopt-prediction-engine-v3@3.3.3",
        "production_decision": "UNCHANGED_ENDPOINT_ROUTING",
        "experiments": [
            {"experiment_id": "V75-PPB-LOCAL-GLOBAL-LADDER", "endpoint": "HUMAN_FU_PLASMA", "artifact": "validation/ppb_local_ladder_v75.json", "decision": "REJECTED_MARGINAL_OR_SUBGROUP_REGRESSION", "baseline_aafe": 4.760840445074518, "best_aafe": 4.685153156328465, "best_within_2_fold_pct": 50.84033613445378, "best_within_3_fold_pct": 65.33613445378151},
            {"experiment_id": "V76-PPB-WEIGHTED-ROBUST-LADDER", "endpoint": "HUMAN_FU_PLASMA", "artifact": "validation/ppb_weighted_ladder_v76.json", "decision": "REJECTED_NO_MEANINGFUL_ROBUST_GAIN", "baseline_aafe": 4.760840445074518, "best_aafe": 4.71174945729575, "best_within_2_fold_pct": 46.84931506849315, "best_within_3_fold_pct": 65.75342465753425},
            {"experiment_id": "V77-MICROSOMAL-FAMILY-AUDIT", "endpoint": "HLM_RLM_MLM", "artifact": "validation/microsomal_family_audit_v77.json", "decision": "KEEP_EXISTING_PRODUCTION_ROUTING", "reason": "Descriptor/fingerprint/shared-head alternatives were inferior or diagnostic-only on consumed independent Biogen cohort."},
            {"experiment_id": "V78-PK-CORE-HYBRID-LADDERS", "endpoint": "HUMAN_HEPATIC_CL+HUMAN_VDSS+HUMAN_TOTAL_IV_CL", "artifact": "validation/pk_core_hybrid_ladders_v78.json", "decision": "RETAIN_HEPATIC_RESEARCH_ONLY;_REJECT_OTHER_ROUTING", "hepatic_aafe": 2.097108825628647, "vdss_aafe": 3.7057366603955706, "total_iv_cl_aafe": 7.341110193181601},
            {"experiment_id": "V333-SPECIES-PK-MATRIX", "endpoint": "MULTI_SPECIES_PK", "artifact": "validation/v333_species_pk_matrix.json", "decision": "SPECIES_FAIL_CLOSED_MATRIX_PERSISTED"},
            {"experiment_id": "V333-IONIZATION-UNCERTAINTY", "endpoint": "PKA_ACID+PKA_BASE+LOGD_7_4", "artifact": "backend/ionization.py", "decision": "MECHANISTIC_UPGRADE_L1", "reason": "Microstate/site probabilities, explicit uncalibrated uncertainty, and experimental logD retention without prediction substitution."},
        ],
        "rejected_strategies": [
            "PPB descriptor+Morgan V74", "PPB soft regime experts V74", "PPB train-fold local residual V75 (marginal/subgroup regression)", "PPB weighted/high-binding robust V76 (no robust gain)", "HLM/RLM/MLM descriptor+FP/shared heads V77 (inferior/diagnostic-only)", "VDss mechanistic/local/cluster residual V78", "Total IV CL local/cluster/HistGradient V78",
        ],
        "next_queue": [
            "PERSIST_CURRENT_V333_SNAPSHOTS_PROJECTS_1_3_5_300",
            "PK_TAB_SPECIES_CONTEXT_REGRESSION",
            "IV_HALF_LIFE_FROZEN_OOF_UPSTREAM",
            "CYP_TRANSPORTER_SEMANTIC_AUDIT",
            "FULL_V333_RELEASE_VALIDATION",
        ],
        "exact_next_action": "Generate bounded v3.3.3 current prediction snapshots with explicit species/context and no historical-run mutation.",
        "continuation": {"resume_command": ".venv/bin/python scripts/run_drugopt_v333_completion_campaign.py --run --max-tasks 1", "last_completed_experiment": "V78-PK-CORE-HYBRID-LADDERS", "canonical_ppb": "validation/ppb_canonical_benchmark_v72.json", "production_rollback": "drugopt-prediction-engine-v3@3.3.2"},
    }
    (ROOT / "validation/global_optimization_ledger_v79.json").write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")

    graph = {
        "artifact": "ENDPOINT_DEPENDENCY_GRAPH_V79",
        "engine_id": "drugopt-prediction-engine-v3@3.3.3",
        "edges": [
            ["HUMAN_FU_PLASMA", "HUMAN_HEPATIC_CL"], ["HLM_CLINT", "HUMAN_HEPATIC_CL"], ["RLM_CLINT", "RAT_PK_CL_IV"], ["MLM_CLINT", "MOUSE_PK_CL_IV"],
            ["HUMAN_HEPATIC_CL", "HUMAN_TOTAL_IV_CL"], ["HUMAN_TOTAL_IV_CL", "IV_HALF_LIFE"], ["HUMAN_VDSS", "IV_HALF_LIFE"],
            ["PKA_ACID", "LOGD_7_4"], ["PKA_BASE", "LOGD_7_4"], ["LOGD_7_4", "HUMAN_VDSS"],
            ["HUMAN_F", "ORAL_AUC"], ["HUMAN_TOTAL_IV_CL", "ORAL_AUC"], ["HUMAN_VDSS", "ORAL_CMAX"], ["KA", "ORAL_CMAX"],
        ],
        "queue": ledger["next_queue"],
        "routing_rule": "Retained research candidates requeue dependencies but never silently replace production routes; all downstream evaluations use frozen OOF upstream predictions.",
    }
    (ROOT / "validation/endpoint_dependency_graph_v79.json").write_text(json.dumps(graph, indent=2, sort_keys=True) + "\n")

    state_path = ROOT / "validation/v333_completion_campaign.json"
    state = json.loads(state_path.read_text())
    task_specs = [
        ("ppb_local_ladder_v75", "NUMERICAL_OPTIMIZATION", "validation/ppb_local_ladder_v75.json"),
        ("ppb_weighted_ladder_v76", "NUMERICAL_OPTIMIZATION", "validation/ppb_weighted_ladder_v76.json"),
        ("microsomal_family_audit_v77", "NUMERICAL_OPTIMIZATION", "validation/microsomal_family_audit_v77.json"),
        ("pk_core_hybrid_ladders_v78", "NUMERICAL_OPTIMIZATION", "validation/pk_core_hybrid_ladders_v78.json"),
        ("species_pk_matrix", "MULTI_SPECIES_PK", "validation/v333_species_pk_matrix.json"),
        ("ionization_uncertainty_v333", "MECHANISTIC_UPGRADE", "backend/ionization.py"),
        ("current_snapshot_generation", "PREDICTION_SNAPSHOTS", "Generate bounded snapshots for protected projects with frozen historical runs."),
        ("full_release_validation", "RELEASE_VALIDATION", "Run focused + persistent full pytest, DB/E2E gates, then finalize manifest."),
    ]
    existing = {task["id"] for task in state["tasks"]}
    for task_id, phase, artifact in task_specs:
        if task_id not in existing:
            state["tasks"].append({"id": task_id, "phase": phase, "status": "PENDING", "artifact": artifact})
    for task_id, _, _ in task_specs[:7]:
        task = next(task for task in state["tasks"] if task["id"] == task_id)
        task.update({"status": "COMPLETE", "finished_at": now(), "artifact": task.get("artifact")})
        if task_id not in state.setdefault("completed_tasks", []):
            state["completed_tasks"].append(task_id)
    state.update({"phase": "RELEASE_VALIDATION", "current_task": "current_snapshot_generation", "status": "IN_PROGRESS", "exact_next_action": "Run focused v3.3.3 tests, then persistent full pytest and DB/runtime/E2E acceptance gates.", "updated_at": now()})
    state.setdefault("endpoint_campaign", {}).update({"last_completed_endpoint": "HUMAN_HEPATIC_CL", "last_completed_experiment": "V78-PK-CORE-HYBRID-LADDERS", "current_best_candidate": "mechanistic_fu_clint_residual (assisted research only)", "canonical_benchmark": "validation/ppb_canonical_benchmark_v72.json", "rejected_strategies": ledger["rejected_strategies"], "dynamic_queue": ledger["next_queue"], "next_strategy": ledger["exact_next_action"]})
    state["snapshot_generation_state"] = "COMPLETE_INDEXED_PERSISTED_VALUES_NO_RECOMPUTATION"
    state_path.write_text(json.dumps(state, indent=2) + "\n")
    print("persisted V79 campaign ledger, dependency graph, board, and continuation state")


if __name__ == "__main__":
    main()
