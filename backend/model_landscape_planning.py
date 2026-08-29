"""Stage 4E-1 model-landscape planning; intentionally not a runtime registry."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from .admet_predictor import MODEL_SPECS
from .endpoint_contracts import ENDPOINT_CONTRACTS
from .endpoint_strategy_registry import get_all_strategies


ROOT = Path(__file__).resolve().parents[1]
REVIEW_VERSION = "stage4e1-model-landscape-v1"
ACCESS_DATE = "2026-08-29"


def source(source_id: str, title: str, url: str, notes: str, *, candidate_id: str = "", source_type: str = "OFFICIAL_DOCUMENTATION", year: str = "") -> dict[str, Any]:
    return {"source_id": source_id, "candidate_id": candidate_id, "source_type": source_type, "title": title, "url": url, "authors": "See primary source", "year": year, "accessed_at": ACCESS_DATE, "license_evidence": notes, "notes": notes}


SOURCES = [
    source("SRC_TDC_ADMET", "TDC ADMET Benchmark Group", "https://tdcommons.ai/benchmark/admet_group/overview/", "Official benchmark documentation; scaffold test partition documented.", source_type="OFFICIAL_DOCUMENTATION"),
    source("SRC_TDC_CACO2", "TDC Caco2_Wang leaderboard", "https://tdcommons.ai/benchmark/admet_group/caco2_wang/", "906 compounds, regression MAE, scaffold split; data-license lineage requires review.", source_type="OFFICIAL_DATASET"),
    source("SRC_TDC_HERG", "TDC hERG leaderboard", "https://tdcommons.ai/benchmark/admet_group/20herg/", "648 binary records, scaffold split; current core lineage overlaps Wang hERG compilation.", source_type="OFFICIAL_DATASET"),
    source("SRC_TDC_ADME", "TDC ADME single-instance tasks", "https://tdcommons.ai/single_pred_tasks/adme/", "Official source catalog for P-gp, clearance and related datasets; original-data terms remain separate.", source_type="OFFICIAL_DATASET"),
    source("SRC_TDC_CLEARANCE", "TDC Clearance_Microsome_AZ leaderboard", "https://tdcommons.ai/benchmark/admet_group/18clmicro/", "1,102 records in mL.min-1.g-1; species/scale mismatch must be resolved before use.", source_type="OFFICIAL_DATASET"),
    source("SRC_OPENADMET", "OpenADMET ExpansionRx Challenge Tutorial", "https://github.com/OpenADMET/ExpansionRx-Challenge-Tutorial", "Official challenge repository lists HLM/MLM CLint and Caco-2 Papp A>B/efflux endpoints; dataset terms require review.", source_type="OFFICIAL_REPOSITORY"),
    source("SRC_BAYESHERG", "BayeshERG official repository", "https://github.com/GIST-CSBL/BayeshERG", "MIT code, but stated trained weights/outputs are CC-BY-NC-4.0; non-commercial restriction blocks internal pharmaceutical deployment.", candidate_id="MODEL_BAYESHERG", source_type="OFFICIAL_REPOSITORY", year="2022"),
    source("SRC_CARDIOGENAI", "CardioGenAI official repository", "https://github.com/gregory-kyro/cardiogenai", "MIT repository with model_parameters tree and CPU option; some required files are external Google Drive and weight terms must be separately confirmed.", candidate_id="MODEL_CARDIOGENAI_HERG", source_type="OFFICIAL_REPOSITORY"),
    source("SRC_METABOGNN", "MetaboGNN official repository", "https://github.com/qwon135/MetaboGNN", "Official repository describes HLM/MLM checkpoints via Google Drive and common-test RMSE; license and output-unit compatibility require review.", candidate_id="MODEL_METABOGNN_CLEARANCE", source_type="OFFICIAL_REPOSITORY"),
    source("SRC_PKASOLVER", "pkasolver official repository", "https://github.com/mayrf/pkasolver", "MIT repository; pkasolver-lite model is distributed, while Epik-transfer weights are explicitly not redistributable.", candidate_id="MODEL_PKASOLVER_LITE", source_type="OFFICIAL_REPOSITORY"),
    source("SRC_PKALEARN", "pKaLearn official repository", "https://github.com/MoitessierLab/pKaLearn", "MIT repository describing GNN ionization-center pKa prediction; publication is a preprint and artifact/weight provenance requires verification.", candidate_id="MODEL_PKALEARN_GNN", source_type="OFFICIAL_REPOSITORY", year="2024"),
    source("SRC_LOGD74", "logD7.4 curated dataset", "https://github.com/nanxstats/logd74", "1,130 structure/logD7.4 records with pH definition; repository page does not establish a reusable data license.", candidate_id="DATA_LOGD74_1130", source_type="OFFICIAL_DATASET", year="2015"),
    source("SRC_MMTKPRED", "MMTKPred transporter kinetics repository", "https://github.com/SizheQiu/MMTKPred", "Predicts transporter Vmax/Km, not Drug-OPT inhibitor/substrate classification contracts; endpoint mismatch.", candidate_id="MODEL_MMTKPRED_TRANSPORTER", source_type="OFFICIAL_REPOSITORY"),
    source("SRC_THEMOL", "TheMol official repository", "https://github.com/themolsubmission/TheMol", "Provides foundation/pretraining material but no already-qualified Drug-OPT endpoint head; a foundation encoder alone is not a usable endpoint predictor.", candidate_id="MODEL_THEMOL_FOUNDATION", source_type="OFFICIAL_REPOSITORY"),
]


def candidate(candidate_id: str, endpoint_id: str, model_name: str, family: str, recommendation: str, source_ids: list[str], **extra: Any) -> dict[str, Any]:
    row = {
        "candidate_id": candidate_id, "endpoint_id": endpoint_id, "model_name": model_name,
        "model_family": family, "architecture": "See official source", "authors_organization": "See official source",
        "publication": "Not independently qualification-grade until Stage 4E-2", "publication_year": None,
        "source_repository": "", "checkpoint_source": "", "checkpoint_available": False,
        "code_available": True, "inference_code_available": False, "training_dataset": "Not fully disclosed/verified",
        "training_n": None, "validation_dataset": "Not independently verified", "validation_n": None,
        "external_validation": "Not established for Drug-OPT contract", "endpoint_definition": "Must match Drug-OPT endpoint contract",
        "assay_definition": "Must be verified", "species": "Must be verified", "units": "Must be verified",
        "classification_threshold": None, "reported_metrics": {}, "probability_calibration": "Not established",
        "applicability_domain": "Not established", "input_representation": "See source", "required_dependencies": [],
        "python_requirements": "REQUIRES_BUILD_TEST", "pytorch_requirements": "REQUIRES_BUILD_TEST", "cuda_requirements": "NONE_CONFIRMED",
        "cpu_inference_support": "UNKNOWN", "arm64_feasibility": "REQUIRES_BUILD_TEST", "estimated_ram": "UNKNOWN",
        "estimated_checkpoint_size": "UNKNOWN", "estimated_latency": "UNKNOWN", "license_code": "LICENSE_UNCLEAR",
        "license_checkpoint": "LICENSE_UNCLEAR", "license_dataset": "LICENSE_UNCLEAR", "commercial_use_status": "LICENSE_UNCLEAR",
        "redistribution_status": "LICENSE_UNCLEAR", "training_overlap_risk": "UNKNOWN", "current_CORE_similarity": "UNKNOWN",
        "expected_model_diversity": "MEDIUM", "expected_error_complementarity": "PLAUSIBLE", "integration_complexity": "HIGH",
        "qualification_priority": "TIER_1", "recommended_action": recommendation, "scientific_notes": "Planning only; not installed or registered.",
        "source_ids": source_ids,
    }
    row.update(extra)
    return row


CANDIDATES = [
    candidate("MODEL_CARDIOGENAI_HERG", "safety_herg_blocker_prob", "CardioGenAI hERG discriminator", "SMILES transformer + discriminator", "GO_LICENSE_REVIEW_FIRST", ["SRC_CARDIOGENAI"], architecture="SMILES autoregressive transformer plus cardiac-ion-channel discriminator", source_repository="https://github.com/gregory-kyro/cardiogenai", checkpoint_source="Repository model_parameters plus required external Google Drive files", checkpoint_available="PARTIAL_EXTERNAL", inference_code_available=True, training_dataset="ChEMBL 33, GuacaMol, MOSES, BindingDB for transformer; hERG-head lineage requires verification", training_overlap_risk="UNKNOWN", current_CORE_similarity="MEDIUM", expected_model_diversity="HIGH", expected_error_complementarity="PLAUSIBLE", cpu_inference_support="DOCUMENTED_CPU_OPTION", arm64_feasibility="REQUIRES_BUILD_TEST", license_code="MIT", license_checkpoint="LICENSE_UNCLEAR", license_dataset="LICENSE_UNCLEAR", commercial_use_status="LIKELY_ALLOWED_BUT_REVIEW_REQUIRED", integration_complexity="HIGH", scientific_notes="Potentially diverse hERG secondary, but no promotion claim until exact head/checkpoint, threshold, labels, overlap, and external validation are audited."),
    candidate("MODEL_BAYESHERG", "safety_herg_blocker_prob", "BayeshERG", "Bayesian graph neural network", "NO_GO_LICENSE", ["SRC_BAYESHERG"], architecture="Bayesian GNN", source_repository="https://github.com/GIST-CSBL/BayeshERG", checkpoint_source="Repository trained model", checkpoint_available=True, inference_code_available=True, required_dependencies=["Python 3.6", "DGL", "PyTorch", "RDKit"], cpu_inference_support="LIKELY", arm64_feasibility="REQUIRES_BUILD_TEST", license_code="MIT", license_checkpoint="CC-BY-NC-4.0", license_dataset="LICENSE_UNCLEAR", commercial_use_status="NONCOMMERCIAL_ONLY", redistribution_status="RESTRICTED", training_overlap_risk="UNKNOWN", current_CORE_similarity="MEDIUM", expected_model_diversity="HIGH", expected_error_complementarity="PLAUSIBLE", scientific_notes="Explicit trained-weight non-commercial restriction makes it ineligible for internal pharmaceutical deployment."),
    candidate("MODEL_METABOGNN_CLEARANCE", "hlm_intrinsic_clearance_scaled_log10", "MetaboGNN HLM/MLM", "graph contrastive learning GNN", "GO_LICENSE_REVIEW_FIRST", ["SRC_METABOGNN"], architecture="GraphCL/pretrained GNN with clearance head", source_repository="https://github.com/qwon135/MetaboGNN", checkpoint_source="Official Google Drive checkpoint archive", checkpoint_available=True, inference_code_available=True, training_dataset="Not fully verified from release metadata", validation_dataset="Common test set reported by repository", reported_metrics={"HLM_RMSE": 28.39, "MLM_RMSE": 27.88}, species="Human and mouse only; rat absent", units="UNRESOLVED; cannot assume current scaled log10(mL/min/kg)", cpu_inference_support="UNKNOWN", arm64_feasibility="REQUIRES_BUILD_TEST", license_code="LICENSE_UNCLEAR", license_checkpoint="LICENSE_UNCLEAR", license_dataset="LICENSE_UNCLEAR", commercial_use_status="LICENSE_UNCLEAR", training_overlap_risk="UNKNOWN", current_CORE_similarity="MEDIUM", expected_model_diversity="HIGH", expected_error_complementarity="PLAUSIBLE", integration_complexity="HIGH", scientific_notes="Potential HLM/MLM challenger only after license and unit/endpoint reconciliation; never applies to RLM, dog, or monkey."),
    candidate("MODEL_PKASOLVER_LITE", "physchem_pka_quantitative_ml", "pkasolver-lite", "graph convolutional network", "GO_ARM64_FEASIBILITY_FIRST", ["SRC_PKASOLVER"], architecture="GNN pKa predictor", source_repository="https://github.com/mayrf/pkasolver", checkpoint_source="Model provided with repository/Colab", checkpoint_available=True, inference_code_available=True, training_dataset="Open pkasolver-lite data; transfer-learning Epik model excluded", endpoint_definition="Ionization-center pKa; micro/macro conversion must be qualified", species="Chemical", units="pKa", cpu_inference_support="LIKELY", arm64_feasibility="REQUIRES_BUILD_TEST", license_code="MIT", license_checkpoint="MIT_REPOSITORY_CLAIM", license_dataset="MIT_REPOSITORY_CLAIM", commercial_use_status="CLEAR_FOR_INTERNAL_RESEARCH", training_overlap_risk="MODERATE_OVERLAP_RISK", current_CORE_similarity="LOW", expected_model_diversity="HIGH", expected_error_complementarity="PLAUSIBLE", integration_complexity="MEDIUM", scientific_notes="Quantitative pKa candidate, not a replacement for rule estimate until site, polyprotic, zwitterion, and external validation qualification."),
    candidate("MODEL_PKALEARN_GNN", "physchem_pka_quantitative_ml", "pKaLearn", "ionization-center GNN", "GO_ARM64_FEASIBILITY_FIRST", ["SRC_PKALEARN"], architecture="GNN predicting pKa at ionizable centers", source_repository="https://github.com/MoitessierLab/pKaLearn", checkpoint_source="Repository model directory; exact released artifact checksum required", checkpoint_available="REPOSITORY_CLAIM_REQUIRES_MANIFEST", inference_code_available=True, training_dataset="Repository/paper lineage to be audited", endpoint_definition="Ionization-center pKa", species="Chemical", units="pKa", cpu_inference_support="LIKELY", arm64_feasibility="REQUIRES_BUILD_TEST", license_code="MIT", license_checkpoint="LICENSE_UNCLEAR", license_dataset="LICENSE_UNCLEAR", commercial_use_status="LIKELY_ALLOWED_BUT_REVIEW_REQUIRED", training_overlap_risk="UNKNOWN", current_CORE_similarity="LOW", expected_model_diversity="HIGH", expected_error_complementarity="PLAUSIBLE", integration_complexity="MEDIUM", scientific_notes="Preprint/repository candidate; exact checkpoint, test-set lineage, and micro/macro semantics remain Stage 4E-2 gates."),
    candidate("MODEL_THEMOL_FOUNDATION", "MULTI_ENDPOINT", "TheMol foundation encoder", "3D molecular foundation model", "NO_GO_NO_CHECKPOINT", ["SRC_THEMOL"], source_repository="https://github.com/themolsubmission/TheMol", checkpoint_available="BASE_ONLY", inference_code_available=True, endpoint_definition="No qualified endpoint-specific head", cpu_inference_support="UNKNOWN", arm64_feasibility="GPU_DEPENDENT", license_code="LICENSE_UNCLEAR", commercial_use_status="LICENSE_UNCLEAR", current_CORE_similarity="LOW", expected_model_diversity="HIGH", expected_error_complementarity="NONE", scientific_notes="A foundation encoder without a released, validated endpoint head is not an endpoint predictor and cannot enter qualification."),
    candidate("MODEL_MMTKPRED_TRANSPORTER", "transporter_pgp_substrate_prob", "MMTKPred", "sequence + SMILES XGBoost kinetics", "NO_GO_ENDPOINT_MISMATCH", ["SRC_MMTKPRED"], source_repository="https://github.com/SizheQiu/MMTKPred", checkpoint_available="UNVERIFIED", inference_code_available=True, endpoint_definition="Transporter Vmax/Km, not P-gp substrate or inhibitor classification", units="Vmax/Km", cpu_inference_support="LIKELY", arm64_feasibility="ARM64_PYTHON_ONLY_LIKELY", license_code="LICENSE_UNCLEAR", commercial_use_status="LICENSE_UNCLEAR", expected_error_complementarity="NONE", scientific_notes="Mechanistically different output is not contract-compatible with Drug-OPT transporter classifications."),
]


DATASETS = [
    {"dataset_id":"DATA_BIOGEN_PROSPECTIVE","endpoint_ids":["solubility_aqueous_logs","ppb_human_percent_bound","hlm_intrinsic_clearance_scaled_log10","rlm_intrinsic_clearance_scaled_log10"],"name":"Biogen Public ADME Prospective Benchmark","n":3521,"assay_quality":"Corporate prospective benchmark, already quarantined in Drug-OPT","source_ids":[],"overlap_risk":"LOW_OVERLAP_RISK","license":"LICENSE_REVIEW_REQUIRED","recommended_action":"GO_DATASET_QUALIFICATION_FIRST","use":"Independent validation only; never train a candidate on it."},
    {"dataset_id":"DATA_EXPANSIONRX","endpoint_ids":["permeability_caco2_logpapp","hlm_intrinsic_clearance_scaled_log10","mlm_intrinsic_clearance_scaled_log10"],"name":"ExpansionRx/OpenADMET challenge data","n":None,"assay_quality":"Exact Caco-2 Papp A>B, HLM/MLM CLint endpoints listed; obtain protocol/terms before use","source_ids":["SRC_OPENADMET"],"overlap_risk":"HIGH_OVERLAP_RISK_FOR_CURRENT_OPENADMET_CLEARANCE","license":"LICENSE_UNCLEAR","recommended_action":"GO_LICENSE_REVIEW_FIRST","use":"Potential Caco-2 external dataset only after raw labels, terms, and overlap de-duplication."},
    {"dataset_id":"DATA_TDC_CACO2_WANG","endpoint_ids":["permeability_caco2_logpapp"],"name":"TDC Caco2_Wang","n":906,"assay_quality":"Caco-2 regression, scaffold benchmark; direction/protocol must be reconciled","source_ids":["SRC_TDC_CACO2","SRC_TDC_ADMET"],"overlap_risk":"HIGH_OVERLAP_RISK","license":"LICENSE_REVIEW_REQUIRED","recommended_action":"NO_GO_OVERLAP_RISK","use":"Not independent: current core trained on Wang compilation."},
    {"dataset_id":"DATA_TDC_HERG","endpoint_ids":["safety_herg_blocker_prob"],"name":"TDC hERG","n":648,"assay_quality":"Binary, scaffold benchmark","source_ids":["SRC_TDC_HERG"],"overlap_risk":"HIGH_OVERLAP_RISK","license":"LICENSE_REVIEW_REQUIRED","recommended_action":"NO_GO_OVERLAP_RISK","use":"Not independent from current Wang-derived hERG lineage."},
    {"dataset_id":"DATA_LOGD74_1130","endpoint_ids":["physchem_logd74_quantitative_ml"],"name":"Curated logD7.4","n":1130,"assay_quality":"pH 7.4 distribution coefficient with structures","source_ids":["SRC_LOGD74"],"overlap_risk":"UNKNOWN","license":"LICENSE_UNCLEAR","recommended_action":"GO_LICENSE_REVIEW_FIRST","use":"Potential independent quantitative-logD benchmark after exact license/protocol and overlap audit."},
    {"dataset_id":"DATA_TDC_CLEARANCE_MICROSOME_AZ","endpoint_ids":["microsomal_clearance_generic"],"name":"TDC Clearance_Microsome_AZ","n":1102,"assay_quality":"Regression, mL.min-1.g-1, scaffold split","source_ids":["SRC_TDC_CLEARANCE"],"overlap_risk":"UNKNOWN","license":"LICENSE_REVIEW_REQUIRED","recommended_action":"NO_GO_ENDPOINT_MISMATCH","use":"Not automatically transformable to species-specific scaled Drug-OPT HLM/RLM/MLM contract."},
]


def build_current_baseline() -> dict[str, Any]:
    stage4d7 = json.loads((ROOT / "validation/stage4d7_endpoint_accuracy_matrix.json").read_text())
    decisions = {row["endpoint_name"]: row for row in stage4d7["endpoints"]}
    inventory = json.loads((ROOT / "validation/stage4d0_current_model_inventory.json").read_text())
    active = {row["endpoint"]: row for row in inventory["active_models"]}
    rows=[]
    for name, policy in sorted(get_all_strategies().items()):
        contract=ENDPOINT_CONTRACTS.get(name)
        spec=MODEL_SPECS.get(name, {})
        row={"endpoint_name":name,"endpoint_id":policy.endpoint_id,"endpoint_definition":contract.scientific_definition if contract else "Policy-only endpoint","species":contract.species if contract else "N/A","assay_definition":contract.assay_type if contract else "N/A","output_type":contract.output_type.value if contract else "N/A","unit":contract.canonical_unit if contract else "N/A","current_production_model":policy.primary_model_ids,"model_version":policy.primary_model_versions,"current_shadow_models":policy.shadow_model_ids,"current_strategy":policy.primary_strategy.value,"calibration_status":policy.calibration_status.value,"ad_status":active.get(name,{}).get("applicability_domain_method","Not applicable/unavailable"),"known_validation":spec.get("validation",active.get(name,{}).get("reported_validation",{})),"independent_validation":spec.get("independent_validation",active.get(name,{}).get("independent_validation",{})),"training_n":active.get(name,{}).get("training_n"),"model_family":spec.get("model_family",active.get(name,{}).get("model_family","rule/mechanistic/unavailable")),"stage4d7_decision":decisions[name]["decision_flags"],"scientific_limitation":policy.limitations,"current_model_gap":decisions[name]["evidence_summary"]}
        rows.append(row)
    return {"artifact":"STAGE4E1_CURRENT_MODEL_BASELINE","review_version":REVIEW_VERSION,"production_changed":False,"endpoint_count":len(rows),"reconciliation":{"strategy_policies":len(get_all_strategies()),"runtime_ml_endpoints":len(MODEL_SPECS),"stage4d7_rows":len(decisions),"contradictions":[]},"endpoints":rows}


def build_priority() -> dict[str, Any]:
    priorities=[
        (1,"hERG discrimination + independent dataset","TIER_1",["MODEL_CARDIOGENAI_HERG","MODEL_BAYESHERG"],["DATA_TDC_HERG"],"Current M1 discrimination/calibration limits; current M2 lacks complementarity."),
        (2,"Caco-2 Papp A>B model/data","TIER_1",[],["DATA_EXPANSIONRX","DATA_TDC_CACO2_WANG"],"Current independent cohort is small; Wang dataset is not independent."),
        (3,"Species-specific microsomal clearance","TIER_1",["MODEL_METABOGNN_CLEARANCE"],["DATA_BIOGEN_PROSPECTIVE","DATA_EXPANSIONRX"],"Downstream PK impact and low current prospective generalization; species isolation mandatory."),
        (4,"Quantitative pKa","TIER_3",["MODEL_PKASOLVER_LITE","MODEL_PKALEARN_GNN"],[],"Current output is rule estimate; quantitative candidate requires micro/macro qualification."),
        (5,"Quantitative logD7.4","TIER_3",[],["DATA_LOGD74_1130"],"Current output is derived estimate; pH-specific independent data is more valuable than an unqualified model."),
        (6,"Transporter coverage","TIER_2",["MODEL_MMTKPRED_TRANSPORTER"],[],"Most endpoints unavailable; identified kinetics model fails classification contract."),
    ]
    return {"artifact":"STAGE4E1_MODEL_GAP_PRIORITY","review_version":REVIEW_VERSION,"score_is_planning_only":True,"dimensions":["scientific_impact","current_weakness","coverage_gap","candidate_quality","dataset_quality","complementarity","license","arm64","integration_effort"],"priorities":[{"rank":r,"gap":g,"tier":t,"model_candidates":m,"dataset_candidates":d,"rationale":why} for r,g,t,m,d,why in priorities]}
