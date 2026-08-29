#!/usr/bin/env python3
"""Generate the fail-closed Stage 4E-2R blocker-resolution record."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation"


def dump(name: str, body: dict) -> None:
    (OUT / name).write_text(json.dumps(body, indent=2, sort_keys=True) + "\n")


def main() -> None:
    blockers = [
        {"candidate": "MODEL_CARDIOGENAI_HERG", "original_blocker": "LEGAL_REVIEW_REQUIRED", "actions_taken": ["Official repository/README/release path inspected"], "sources_checked": ["https://github.com/gregory-kyro/cardiogenai"], "resolved": False, "resolution": "External unversioned parameters and exact hERG head/license remain unresolved.", "final_decision": "BLOCKER_UNRESOLVED_CARDIOGENAI"},
        {"candidate": "MODEL_METABOGNN_CLEARANCE", "original_blocker": "LEGAL_REVIEW_REQUIRED", "actions_taken": ["Official repository, checkpoint path, species/units inspected"], "sources_checked": ["https://github.com/qwon135/MetaboGNN"], "resolved": False, "resolution": "License and direct species/unit contract are unresolved.", "final_decision": "BLOCKER_UNRESOLVED_METABOGNN"},
        {"candidate": "MODEL_PKASOLVER_LITE", "original_blocker": "ARM64_WORKAROUND_REQUIRED", "actions_taken": ["Official source snapshot pinned", "PyG isolated environment installed on aarch64", "CPU model-load attempted"], "sources_checked": ["https://github.com/mayrf/pkasolver"], "resolved": False, "resolution": "Source-contained lite checkpoint is incompatible with current released PyG architecture; strict state-dict load fails.", "final_decision": "BLOCKER_UNRESOLVED_REPRODUCIBILITY"},
        {"candidate": "MODEL_PKALEARN_GNN", "original_blocker": "LEGAL_REVIEW_REQUIRED", "actions_taken": ["Official repository/weight layout/environment inspected"], "sources_checked": ["https://github.com/MoitessierLab/pKaLearn"], "resolved": False, "resolution": "Weight/data lineage remains legally unresolved; CUDA/Windows-oriented environment is not ARM64 evidence.", "final_decision": "BLOCKER_UNRESOLVED_PKALEARN"},
        {"candidate": "DATA_BIOGEN_PROSPECTIVE", "original_blocker": "NO_GO_DATA_ACCESS", "actions_taken": ["Local lineage and public source chain rechecked"], "sources_checked": ["validation/stage4d0_dataset_lineage.json"], "resolved": False, "resolution": "No immutable public raw structure/value release or reuse grant located.", "final_decision": "NO_GO_DATA_ACCESS"},
        {"candidate": "DATA_EXPANSIONRX", "original_blocker": "LICENSE_REVIEW_REQUIRED", "actions_taken": ["Official OpenADMET/Hugging Face release inspected", "Pinned raw CSV acquired", "SHA256 and standardization intake completed"], "sources_checked": ["https://openadmet.org/datasetsmodels/", "https://huggingface.co/datasets/openadmet/openadmet-expansionrx-challenge-data"], "resolved": True, "resolution": "CC-BY-4.0 official release with raw SMILES, Caco-2 A→B values, HLM/RLM/MLM columns, DOI, pinned revision, and SHA256.", "final_decision": "PASS_WITH_EXCLUSIONS"},
        {"candidate": "DATA_LOGD74_1130", "original_blocker": "LICENSE_REVIEW_REQUIRED", "actions_taken": ["Official repository/source-license path rechecked"], "sources_checked": ["https://github.com/nanxstats/logd74"], "resolved": False, "resolution": "Reusable raw-data license remains unresolved.", "final_decision": "BLOCKER_UNRESOLVED_LICENSE"},
    ]
    dump("stage4e2r_blocker_resolution.json", {"stage": "4E-2R", "production_changed": False, "records": blockers})
    models = [
        {"candidate_id": "MODEL_ADMET_AI_V2_HERG", "endpoint": "safety_herg_blocker_prob", "source": "https://pypi.org/project/admet-ai/", "checkpoint": "admet-ai==2.0.1 wheel, SHA256 fef3527f637abb00d272cf824e8eef0136fe31ebde6c56881f1a8c02c0417806; bundled model assets", "license": "MIT", "architecture": "Chemprop v2 message-passing GNN", "training_lineage": "TDC hERG (648); high overlap risk with TDC/Wang-derived audits", "endpoint_semantics": "Binary hERG blocker classification", "arm64": "Package dependency resolution selected CUDA-heavy Torch 2.13 stack; no CPU smoke completed", "external_benchmark_path": "Requires a non-TDC hERG external dataset", "decision": "ARM64_WORKAROUND_REQUIRED"},
    ]
    dump("stage4e2r_replacement_candidates.json", {"replacement_search_scope": "hERG only after CardioGenAI blocker", "candidates": models})
    datasets = [
        {"dataset_id": "DATA_OPENADMET_EXPANSIONRX_CACO2_PAPP_AB", "endpoint": "permeability_caco2_logpapp", "source": "https://huggingface.co/datasets/openadmet/openadmet-expansionrx-challenge-data", "revision": "6b898ccc43d10d25b230fb09e22a6e30c30022b5", "doi": "10.57967/hf/9687", "n": 7618, "structures": "SMILES; 7,618 valid after current standardizer", "values": "Caco-2 Permeability Papp A>B; 3,773 numeric and 33 censored", "unit": "Papp A→B, source numeric unit recorded as source unit; log transform is deferred to Stage 4E-3 contract reconciliation", "license": "CC-BY-4.0", "lineage": "Expansion Therapeutics post-challenge public release", "overlap_risk": "UNKNOWN_FOR_ADMETICA_CORE; exact canonical-SMILES exclusions required", "acquisition": {"filename": "expansion_data_raw.csv", "sha256": "f674ec74cca1146bc386f832a32d4b8d921d3c312f92cb436cc005901c724a3c", "outside_git": True}, "decision": "PASS_WITH_EXCLUSIONS"}
    ]
    dump("stage4e2r_dataset_candidates.json", {"datasets": datasets})
    dump("stage4e2r_pkasolver_arm64.json", {"candidate_id": "MODEL_PKASOLVER_LITE", "architecture": "aarch64", "environment": {"python": "3.11", "torch": "2.8.0+cpu (read-only production site packages)", "torch_geometric": "2.8.0.post1 isolated"}, "build_method": "isolated environment with official ARM64 wheels for torch_geometric and its Python dependencies", "load_success": False, "test_panel": "not run because QueryModel construction failed before inference", "cold_latency_seconds": None, "warm_latency_seconds": None, "memory_mb": None, "determinism": "not testable", "errors": ["Checkpoint state-dict keys expect legacy PyG GIN layout; current released PyG model layout differs."], "decision": "NO_GO_REPRODUCIBILITY"})
    gate = [
        {"asset": "MODEL_ADMET_AI_V2_HERG", "license_pass": True, "checkpoint_or_data_pass": True, "endpoint_pass": True, "runtime_pass": False, "overlap_plan_pass": False, "external_dataset_available": False, "stage4e3_eligible": False, "final": "FAIL_ARM64_RUNTIME"},
        {"asset": "DATA_OPENADMET_EXPANSIONRX_CACO2_PAPP_AB", "license_pass": True, "checkpoint_or_data_pass": True, "endpoint_pass": True, "runtime_pass": True, "overlap_plan_pass": True, "external_dataset_available": True, "stage4e3_eligible": True, "final": "PASS_AS_STAGE4E3_EXTERNAL_WITH_EXCLUSIONS"},
    ]
    dump("stage4e2r_stage4e3_entry_gate.json", {"production_changed": False, "assets": gate})
    dump("stage4e2r_stage4e3_plan.json", {"ready_state": "PARTIAL_READY_FOR_STAGE_4E3", "eligible_model_ids": [], "eligible_dataset_ids": ["DATA_OPENADMET_EXPANSIONRX_CACO2_PAPP_AB"], "pilot_size": {"models": 0, "datasets": 1}, "plans": [{"asset": "DATA_OPENADMET_EXPANSIONRX_CACO2_PAPP_AB", "current_core_comparator": "admetica_caco2", "purpose": "External Caco-2 A→B benchmark only; no model activation", "exact_overlap_exclusion": "Current standardizer canonical SMILES; remove duplicates and any available core-training/Stage4D cohort matches", "scaffold_strategy": "Murcko scaffold split after exclusion; chemically meaningful grouping for acyclic compounds", "metrics": ["MAE", "RMSE", "bias", "Spearman", "within-2-fold", "within-3-fold"], "bootstrap": "paired 1,000 replicate comparison against current core", "ad_analysis": "report in-domain, borderline, and OOD separately", "complementarity": "only if a future independently qualified secondary Caco-2 model exists"}]})


if __name__ == "__main__":
    main()
