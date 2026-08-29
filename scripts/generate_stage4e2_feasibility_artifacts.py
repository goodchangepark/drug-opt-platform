#!/usr/bin/env python3
"""Emit the Stage 4E-2 technical/legal gate record.

This deliberately has no imports from the prediction runtime.  It records a
fail-closed acquisition review, not a model registry or a scientific result.
"""
from __future__ import annotations

import hashlib
import json
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation"
ACCESS_DATE = "2026-08-29"


def dump(name: str, body: dict) -> None:
    (OUT / name).write_text(json.dumps(body, indent=2, sort_keys=True) + "\n")


def record_sha256(path: str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main() -> None:
    # Repository snapshots were intentionally inspected outside the repository
    # and then deleted.  The two hashes below identify the small, selected
    # files encountered in those immutable source snapshots; no model binary is
    # retained or committed by Drug-OPT.
    source_rows = [
        {"source_id": "SRC_CARDIOGENAI_E2", "candidate_id": "MODEL_CARDIOGENAI_HERG", "source_type": "OFFICIAL_REPOSITORY", "title": "CardioGenAI", "url": "https://github.com/gregory-kyro/cardiogenai", "commit": "2f403b84bb6fb4a44762ef8d118377a2b1529416", "accessed_at": ACCESS_DATE, "license_evidence": "Repository MIT; no separate checkpoint grant located.", "notes": "Required parameter files are partly external Google Drive assets without a release manifest/checksum."},
        {"source_id": "SRC_METABOGNN_E2", "candidate_id": "MODEL_METABOGNN_CLEARANCE", "source_type": "OFFICIAL_REPOSITORY", "title": "MetaboGNN", "url": "https://github.com/qwon135/MetaboGNN", "commit": "c442f715356c5554902c1fc1a53521f28fa39172", "accessed_at": ACCESS_DATE, "license_evidence": "No unambiguous code, checkpoint, or data license evidence found in the release path.", "notes": "Google Drive checkpoint archive is not a versioned release manifest; HLM/MLM units are unresolved and RLM is not claimed."},
        {"source_id": "SRC_PKASOLVER_E2", "candidate_id": "MODEL_PKASOLVER_LITE", "source_type": "OFFICIAL_REPOSITORY", "title": "pkasolver", "url": "https://github.com/mayrf/pkasolver", "commit": "a6ec86e0337474efc9a6797b5a6374e6961a2748", "accessed_at": ACCESS_DATE, "license_evidence": "LICENSE.md is MIT; README explicitly excludes Epik-transfer weights but distributes pkasolver-lite.", "notes": "Repository snapshot contained versioned lite ensemble weights and a microstate/site-oriented API."},
        {"source_id": "SRC_PKALEARN_E2", "candidate_id": "MODEL_PKALEARN_GNN", "source_type": "OFFICIAL_REPOSITORY", "title": "pKaLearn", "url": "https://github.com/MoitessierLab/pKaLearn", "commit": "ba07b6887a32f027e70510b9565bf61d48cada66", "accessed_at": ACCESS_DATE, "license_evidence": "Repository LICENSE is MIT; repository-contained weight provenance and training-data rights still need legal/lineage review.", "notes": "README documents ionizable-center/pH inference and requires torch_geometric; supplied environment is CUDA/Windows oriented."},
        {"source_id": "SRC_BIOGEN_E2", "candidate_id": "DATA_BIOGEN_PROSPECTIVE", "source_type": "LOCAL_LINEAGE_METADATA", "title": "Biogen Public ADME Prospective Benchmark", "url": "validation/stage4d0_dataset_lineage.json", "accessed_at": ACCESS_DATE, "license_evidence": "Local lineage metadata does not grant raw-data acquisition or reuse rights.", "notes": "Reported N=3,521; no immutable public raw-data release was available to acquire in this gate."},
        {"source_id": "SRC_EXPANSIONRX_E2", "candidate_id": "DATA_EXPANSIONRX", "source_type": "OFFICIAL_DOCUMENTATION", "title": "OpenADMET ExpansionRx challenge", "url": "https://openadmet.org/", "accessed_at": ACCESS_DATE, "license_evidence": "Endpoint listing found; raw-data terms and reusable structured release not established.", "notes": "Caco-2 A→B, HLM, and MLM must remain separate contracts."},
        {"source_id": "SRC_LOGD74_E2", "candidate_id": "DATA_LOGD74_1130", "source_type": "OFFICIAL_REPOSITORY", "title": "logD7.4 curated dataset", "url": "https://github.com/nanxstats/logd74", "accessed_at": ACCESS_DATE, "license_evidence": "Repository source did not establish a reusable data license.", "notes": "pH 7.4 claim alone is not a redistribution/internal-use grant."},
    ]
    dump("stage4e2_source_manifest.json", {"stage": "4E-2", "accessed_at": ACCESS_DATE, "sources": source_rows})

    licenses = [
        {"candidate_id": "MODEL_CARDIOGENAI_HERG", "code": "PASS_INTERNAL_RESEARCH", "weights": "LICENSE_UNCLEAR", "data": "LICENSE_UNCLEAR", "commercial_internal_research": "LICENSE_UNCLEAR", "redistribution": "LICENSE_UNCLEAR", "source_ids": ["SRC_CARDIOGENAI_E2"], "decision": "LEGAL_REVIEW_REQUIRED"},
        {"candidate_id": "MODEL_METABOGNN_CLEARANCE", "code": "LICENSE_UNCLEAR", "weights": "LICENSE_UNCLEAR", "data": "LICENSE_UNCLEAR", "commercial_internal_research": "LICENSE_UNCLEAR", "redistribution": "LICENSE_UNCLEAR", "source_ids": ["SRC_METABOGNN_E2"], "decision": "LEGAL_REVIEW_REQUIRED"},
        {"candidate_id": "MODEL_PKASOLVER_LITE", "code": "PASS_INTERNAL_RESEARCH", "weights": "PASS_INTERNAL_RESEARCH", "data": "PASS_WITH_LEGAL_REVIEW", "commercial_internal_research": "PASS_INTERNAL_RESEARCH", "redistribution": "MIT_NOTICE_REQUIRED", "source_ids": ["SRC_PKASOLVER_E2"], "decision": "PASS_INTERNAL_RESEARCH"},
        {"candidate_id": "MODEL_PKALEARN_GNN", "code": "PASS_INTERNAL_RESEARCH", "weights": "PASS_WITH_LEGAL_REVIEW", "data": "PASS_WITH_LEGAL_REVIEW", "commercial_internal_research": "PASS_WITH_LEGAL_REVIEW", "redistribution": "MIT_NOTICE_REQUIRED", "source_ids": ["SRC_PKALEARN_E2"], "decision": "PASS_WITH_LEGAL_REVIEW"},
        {"candidate_id": "DATA_BIOGEN_PROSPECTIVE", "code": "NOT_APPLICABLE", "weights": "NOT_APPLICABLE", "data": "LICENSE_UNCLEAR", "commercial_internal_research": "LICENSE_UNCLEAR", "redistribution": "LICENSE_UNCLEAR", "source_ids": ["SRC_BIOGEN_E2"], "decision": "LICENSE_REVIEW_REQUIRED"},
        {"candidate_id": "DATA_EXPANSIONRX", "code": "NOT_APPLICABLE", "weights": "NOT_APPLICABLE", "data": "LICENSE_UNCLEAR", "commercial_internal_research": "LICENSE_UNCLEAR", "redistribution": "LICENSE_UNCLEAR", "source_ids": ["SRC_EXPANSIONRX_E2"], "decision": "LICENSE_REVIEW_REQUIRED"},
        {"candidate_id": "DATA_LOGD74_1130", "code": "NOT_APPLICABLE", "weights": "NOT_APPLICABLE", "data": "LICENSE_UNCLEAR", "commercial_internal_research": "LICENSE_UNCLEAR", "redistribution": "LICENSE_UNCLEAR", "source_ids": ["SRC_LOGD74_E2"], "decision": "LICENSE_REVIEW_REQUIRED"},
    ]
    dump("stage4e2_license_matrix.json", {"policy": "FAIL_CLOSED", "entries": licenses})

    contracts = [
        {"candidate_id": "MODEL_CARDIOGENAI_HERG", "drug_opt_endpoint": "safety_herg_blocker_prob", "candidate_endpoint": "cardiac-ion-channel output; exact hERG head label/threshold unresolved", "species": "not release-qualified", "assay": "unresolved", "unit": "probability/risk semantics unresolved", "output_type": "unresolved classification head", "threshold": "unresolved", "transformation_required": "none may be assumed", "compatibility_status": "ENDPOINT_TRANSFORM_REQUIRES_VALIDATION", "notes": "A continuous score or alternate assay cutoff cannot be silently mapped to the blocker probability contract."},
        {"candidate_id": "MODEL_METABOGNN_CLEARANCE", "drug_opt_endpoint": "hlm_intrinsic_clearance_scaled_log10; mlm_intrinsic_clearance_scaled_log10", "candidate_endpoint": "HLM/MLM only; RLM absent", "species": "human, mouse claimed; rat not claimed", "assay": "intrinsic clearance context unresolved", "unit": "unresolved", "output_type": "regression", "threshold": "not applicable", "transformation_required": "ASSAY_CONTEXT_REQUIRED", "compatibility_status": "ASSAY_CONTEXT_REQUIRED", "notes": "No deterministic conversion is valid until source units, protein/scaling, and log transform are released."},
        {"candidate_id": "MODEL_PKASOLVER_LITE", "drug_opt_endpoint": "physchem_pka_quantitative_ml (currently unavailable); existing ionization_pka_estimated is rule based", "candidate_endpoint": "microstate/site pKa ensemble", "species": "chemical", "assay": "pKa", "unit": "pKa", "output_type": "site-specific regression with uncertainty", "threshold": "not applicable", "transformation_required": "macro-pKa/endpoint policy validation", "compatibility_status": "ENDPOINT_COMPATIBLE_WITH_NEW_CONTRACT", "notes": "Must not replace the current rule estimate before ionization-site, polyprotic, and zwitterion qualification."},
        {"candidate_id": "MODEL_PKALEARN_GNN", "drug_opt_endpoint": "physchem_pka_quantitative_ml (currently unavailable); existing ionization_pka_estimated is rule based", "candidate_endpoint": "ionizable-center pKa plus iterative protonation state", "species": "chemical", "assay": "pKa", "unit": "pKa", "output_type": "ionization-center regression", "threshold": "not applicable", "transformation_required": "macro-pKa/endpoint policy validation", "compatibility_status": "ENDPOINT_COMPATIBLE_WITH_NEW_CONTRACT", "notes": "API semantics are technically adequate but training/weight lineage remains incomplete."},
    ]
    dump("stage4e2_endpoint_contract_matrix.json", {"production_contracts_changed": False, "entries": contracts})

    assets = [
        {"candidate_id": "MODEL_CARDIOGENAI_HERG", "acquired": False, "source": "SRC_CARDIOGENAI_E2", "version": "git 2f403b84bb6fb4a44762ef8d118377a2b1529416", "filename": None, "size_bytes": None, "sha256": None, "license": "LICENSE_UNCLEAR", "local_location": None, "download_date": None, "reason": "Early license/checkpoint-manifest gate; no model weights downloaded."},
        {"candidate_id": "MODEL_METABOGNN_CLEARANCE", "acquired": False, "source": "SRC_METABOGNN_E2", "version": "git c442f715356c5554902c1fc1a53521f28fa39172", "filename": None, "size_bytes": None, "sha256": None, "license": "LICENSE_UNCLEAR", "local_location": None, "download_date": None, "reason": "Early license and endpoint-unit gate; no checkpoint archive downloaded."},
        {"candidate_id": "MODEL_PKASOLVER_LITE", "acquired": False, "source": "SRC_PKASOLVER_E2", "version": "git a6ec86e0337474efc9a6797b5a6374e6961a2748", "filename": "best_model_0.pt (source-inspected; not retained)", "size_bytes": 2167218, "sha256": "8190bc02122478812eda0d4fd7a41e39cd9dfce803fbc293e936516845877ab2", "license": "MIT repository distribution", "local_location": None, "download_date": ACCESS_DATE, "reason": "Temporary source clone was removed after identity/hash inspection; no binary is stored in the repository."},
        {"candidate_id": "MODEL_PKALEARN_GNN", "acquired": False, "source": "SRC_PKALEARN_E2", "version": "git ba07b6887a32f027e70510b9565bf61d48cada66", "filename": "train_AAc-1_best.pth (source-inspected; not retained)", "size_bytes": 9864714, "sha256": "8e6fcd437d87e176740104bfc3533106a4e2603af99fb5d3a8bc3ac160145f45", "license": "PASS_WITH_LEGAL_REVIEW", "local_location": None, "download_date": ACCESS_DATE, "reason": "Temporary source clone was removed after identity/hash inspection; no binary is stored in the repository."},
    ]
    dump("stage4e2_asset_manifest.json", {"assets_are_not_committed": True, "assets": assets})

    arm = [
        {"candidate_id": "MODEL_CARDIOGENAI_HERG", "install_status": "NOT_ATTEMPTED_EARLY_GATE", "python": "unverified", "framework": "PyTorch", "compiled_extensions": "unverified", "CPU_load": "NOT_TESTED", "CPU_inference": "NOT_TESTED", "GPU_possible": "documented but unverified", "cold_seconds": None, "warm_seconds": None, "batch_seconds": None, "peak_memory_mb": None, "checkpoint_mb": None, "deterministic": "NOT_TESTED", "errors": ["Exact legal checkpoint and endpoint head were not qualified."], "decision": "LEGAL_REVIEW_REQUIRED"},
        {"candidate_id": "MODEL_METABOGNN_CLEARANCE", "install_status": "NOT_ATTEMPTED_EARLY_GATE", "python": "unverified", "framework": "Graph GNN; exact stack unresolved", "compiled_extensions": "PyG/DGL risk unresolved", "CPU_load": "NOT_TESTED", "CPU_inference": "NOT_TESTED", "GPU_possible": "unverified", "cold_seconds": None, "warm_seconds": None, "batch_seconds": None, "peak_memory_mb": None, "checkpoint_mb": None, "deterministic": "NOT_TESTED", "errors": ["License and endpoint-unit gates failed closed before installation."], "decision": "LEGAL_REVIEW_REQUIRED"},
        {"candidate_id": "MODEL_PKASOLVER_LITE", "install_status": "NOT_INSTALLED_PRODUCTION_SAFE", "python": "production Python has torch 2.8.0+cpu; torch_geometric absent", "framework": "PyTorch + torch_geometric", "compiled_extensions": "torch_geometric installation/build required on ARM64", "CPU_load": "NOT_TESTED", "CPU_inference": "NOT_TESTED", "GPU_possible": "not required by source", "cold_seconds": None, "warm_seconds": None, "batch_seconds": None, "peak_memory_mb": None, "checkpoint_mb": 2.07, "deterministic": "NOT_TESTED", "errors": ["No isolated ARM64 PyG environment was established; production environment deliberately untouched."], "decision": "ARM64_WORKAROUND_REQUIRED"},
        {"candidate_id": "MODEL_PKALEARN_GNN", "install_status": "NOT_INSTALLED_PRODUCTION_SAFE", "python": "production Python has torch 2.8.0+cpu; torch_geometric absent", "framework": "PyTorch + torch_geometric", "compiled_extensions": "torch_geometric installation/build required on ARM64", "CPU_load": "NOT_TESTED", "CPU_inference": "NOT_TESTED", "GPU_possible": "source environment pins pytorch-cuda=12.4 and Windows prefix", "cold_seconds": None, "warm_seconds": None, "batch_seconds": None, "peak_memory_mb": None, "checkpoint_mb": 9.41, "deterministic": "NOT_TESTED", "errors": ["CUDA/Windows-oriented supplied environment is not Xavier evidence; legal/data lineage review incomplete."], "decision": "ARM64_WORKAROUND_REQUIRED"},
    ]
    dump("stage4e2_arm64_runtime_matrix.json", {"host_architecture": "aarch64 / Jetson Xavier target", "production_environment_modified": False, "entries": arm})

    datasets = [
        {"dataset_id": "DATA_BIOGEN_PROSPECTIVE", "source": "SRC_BIOGEN_E2", "raw_n": 3521, "usable_n": 0, "endpoint_counts": {}, "unit_inventory": [], "license": "LICENSE_UNCLEAR", "checksum": None, "overlap_status": "NOT_AUDITABLE_WITHOUT_RAW_STRUCTURES", "processing_script": None, "decision": "NO_GO_DATA_ACCESS", "notes": "Lineage metadata exists but no raw immutable source/terms were available; no data acquired."},
        {"dataset_id": "DATA_EXPANSIONRX", "source": "SRC_EXPANSIONRX_E2", "raw_n": None, "usable_n": 0, "endpoint_counts": {}, "unit_inventory": [], "license": "LICENSE_UNCLEAR", "checksum": None, "overlap_status": "NOT_AUDITABLE_WITHOUT_RAW_STRUCTURES", "processing_script": None, "decision": "LICENSE_REVIEW_REQUIRED", "notes": "Raw structure/value release and usage terms require resolution; protocol-specific Caco-2 A→B must be retained."},
        {"dataset_id": "DATA_LOGD74_1130", "source": "SRC_LOGD74_E2", "raw_n": 1130, "usable_n": 0, "endpoint_counts": {}, "unit_inventory": ["logD at pH 7.4 (reported)"], "license": "LICENSE_UNCLEAR", "checksum": None, "overlap_status": "NOT_AUDITABLE_WITHOUT_RAW_STRUCTURES", "processing_script": None, "decision": "LICENSE_REVIEW_REQUIRED", "notes": "No raw copy acquired because a reusable dataset license was not established."},
    ]
    dump("stage4e2_dataset_manifest.json", {"raw_data_committed": False, "datasets": datasets})
    dump("stage4e2_overlap_audit.json", {"method": "Exact canonical-SMILES overlap is required once legally acquired raw structures exist.", "performed": False, "reason": "No listed dataset passed license/source access gates for raw acquisition.", "datasets": [{"dataset_id": row["dataset_id"], "status": row["overlap_status"], "decision": row["decision"]} for row in datasets]})

    model_decisions = [
        {"candidate_id": "MODEL_CARDIOGENAI_HERG", "decision": "LEGAL_REVIEW_REQUIRED", "failed_gate": "GATE_2_LICENSE", "gate_summary": "Code MIT; checkpoint/data terms and versioned hERG head identity unavailable.", "may_enter_runtime": False},
        {"candidate_id": "MODEL_METABOGNN_CLEARANCE", "decision": "LEGAL_REVIEW_REQUIRED", "failed_gate": "GATE_2_LICENSE", "gate_summary": "Code/weights/data rights unresolved; output unit contract also unresolved.", "may_enter_runtime": False},
        {"candidate_id": "MODEL_PKASOLVER_LITE", "decision": "ARM64_WORKAROUND_REQUIRED", "failed_gate": "GATE_8_ARM64_RUNTIME", "gate_summary": "MIT source-contained lite weights and microstate semantics pass earlier technical gates; PyG ARM64 isolated runtime has not been reproduced.", "may_enter_runtime": False},
        {"candidate_id": "MODEL_PKALEARN_GNN", "decision": "LEGAL_REVIEW_REQUIRED", "failed_gate": "GATE_2_LICENSE", "gate_summary": "MIT code and source-contained weights found, but weight/data lineage requires review; PyG ARM64 runtime is also unproven.", "may_enter_runtime": False},
    ]
    dataset_decisions = [{"dataset_id": row["dataset_id"], "decision": row["decision"], "may_be_used_as_external_benchmark": False} for row in datasets]
    dump("stage4e2_candidate_decisions.json", {"stage": "4E-2", "scientific_superiority_established": False, "production_changed": False, "models": model_decisions, "datasets": dataset_decisions})
    dump("stage4e2_stage4e3_plan.json", {"stage": "4E-3", "ready": False, "pass_model_ids": [], "pass_dataset_ids": [], "blocking_conditions": ["Resolve separate code/weight/data rights for CardioGenAI, MetaboGNN, pKaLearn, and datasets.", "Build and record a pinned isolated ARM64 PyG runtime for pkasolver-lite (and pKaLearn only after lineage review).", "Acquire immutable licensed raw benchmark structures; then perform exact canonical-SMILES overlap exclusions before any metrics."], "conditional_plans": [{"candidate_id": "MODEL_PKASOLVER_LITE", "core_comparator": "ionization_smarts_rules_v1 (rule estimate; not an ML accuracy comparator)", "external_dataset": "licensed site-resolved pKa benchmark to be acquired", "metrics": ["site pKa MAE", "RMSE", "polyprotic completeness", "zwitterion handling"], "overlap_control": "exact canonical-SMILES exclusion against training/source lineage", "split": "scaffold-aware held-out", "bootstrap": "paired 1000-replicate", "ad_analysis": "in-domain/borderline/OOD", "runtime_role_if_successful": "new unregistered quantitative-pKa candidate only"}]})


if __name__ == "__main__":
    main()
