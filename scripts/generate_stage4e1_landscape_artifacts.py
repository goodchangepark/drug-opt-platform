"""Generate Stage 4E-1 planning artifacts; no network, model, or DB access."""
from __future__ import annotations
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT)) if str(ROOT) not in sys.path else None
from backend.model_landscape_planning import CANDIDATES, DATASETS, REVIEW_VERSION, SOURCES, build_current_baseline, build_priority

def write(name, value):
    (ROOT/"validation"/name).write_text(json.dumps(value,indent=2,sort_keys=True)+"\n")

def main():
    baseline=build_current_baseline()
    candidates={"artifact":"STAGE4E1_CANDIDATE_MODEL_LANDSCAPE","review_version":REVIEW_VERSION,"installed_or_registered":False,"candidates":CANDIDATES}
    datasets={"artifact":"STAGE4E1_DATASET_LANDSCAPE","review_version":REVIEW_VERSION,"datasets":DATASETS}
    priority=build_priority()
    pilot={"artifact":"STAGE4E1_STAGE4E2_PILOT_PLAN","review_version":REVIEW_VERSION,"no_installation_in_stage4e1":True,"pilot_model_ids":["MODEL_CARDIOGENAI_HERG","MODEL_METABOGNN_CLEARANCE","MODEL_PKASOLVER_LITE","MODEL_PKALEARN_GNN"],"pilot_dataset_ids":["DATA_BIOGEN_PROSPECTIVE","DATA_EXPANSIONRX","DATA_LOGD74_1130"],"universal_gate":["verify_code_checkpoint_dataset_license_separately","download_only_after_license_gate","record_checksum","verify_endpoint_contract_and_units","exact_structure_overlap_exclusion","scaffold_split_and_external_holdout","paired_bootstrap_against_core","AD_and_disagreement_analysis","ARM64_cold_warm_benchmark","stage4d5_promotion_gates"],"plans":[{"candidate_id":c["candidate_id"],"adapter":"New isolated Stage4E2 shadow adapter only","baseline":"Current endpoint CORE","metrics":"Contract-appropriate metrics plus paired bootstrap","overlap_control":"Exact canonical-SMILES exclusion against training and validation sources","xavier_benchmark":"CPU cold/warm latency, RAM, import/build reliability","promotion":"Stage4D5 fail-closed gates; never auto-activate"} for c in CANDIDATES if c["candidate_id"] in {"MODEL_CARDIOGENAI_HERG","MODEL_METABOGNN_CLEARANCE","MODEL_PKASOLVER_LITE","MODEL_PKALEARN_GNN"}]}
    manifest={"artifact":"STAGE4E1_SOURCE_MANIFEST","review_version":REVIEW_VERSION,"sources":SOURCES}
    license_matrix={"artifact":"STAGE4E1_LICENSE_MATRIX","review_version":REVIEW_VERSION,"candidates":[{"candidate_id":c["candidate_id"],"code":c["license_code"],"checkpoint":c["license_checkpoint"],"dataset":c["license_dataset"],"commercial":c["commercial_use_status"]} for c in CANDIDATES]}
    arm={"artifact":"STAGE4E1_ARM64_FEASIBILITY","review_version":REVIEW_VERSION,"candidates":[{"candidate_id":c["candidate_id"],"arm64_feasibility":c["arm64_feasibility"],"cpu":c["cpu_inference_support"],"dependencies":c["required_dependencies"]} for c in CANDIDATES]}
    for name,value in {"stage4e1_current_model_baseline.json":baseline,"stage4e1_candidate_model_landscape.json":candidates,"stage4e1_dataset_landscape.json":datasets,"stage4e1_model_gap_priority.json":priority,"stage4e1_stage4e2_pilot_plan.json":pilot,"stage4e1_source_manifest.json":manifest,"stage4e1_license_matrix.json":license_matrix,"stage4e1_arm64_feasibility.json":arm}.items(): write(name,value)
if __name__=="__main__": main()
