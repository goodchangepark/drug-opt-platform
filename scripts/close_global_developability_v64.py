"""Close Global Developability v1 and freeze pKa/logD data readiness.

No pKa/logD values are synthesized here.  Existing rule/derived outputs are
kept separate from experimental training truth.
"""
from __future__ import annotations

import hashlib
import json
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TODAY = str(date.today())


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def dataset(endpoint, unit, definition, sources, blockers):
    data = {
        "dataset_version": f"{endpoint.lower()}_experimental_v1",
        "created_at": TODAY,
        "endpoint": endpoint,
        "canonical_unit": unit,
        "endpoint_definition": definition,
        "qualified_records": [],
        "qualified_n": 0,
        "acid_n": 0 if endpoint == "PKA" else None,
        "base_n": 0 if endpoint == "PKA" else None,
        "source_families": sources,
        "independence": "NO_QUALIFIED_SOURCE_LINEAGE",
        "rejections": blockers,
        "status": "DATA_LIMITED",
    }
    data["dataset_hash"] = hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()
    return data


def main():
    pka = dataset(
        "PKA", "pKa", "experimental acidic/basic dissociation constant; macro/micro/site semantics retained",
        ["ChEMBL", "PubChem", "pKaLearn", "pKaSolver"],
        [
            {"reason": "MODELLED_OR_RULE_DERIVED", "count": "not accepted", "detail": "Drug-OPT IonizationEngine estimates are not experimental labels."},
            {"reason": "CHECKPOINT_OR_SEMANTICS_UNRESOLVED", "count": "not accepted", "detail": "pKaSolver strict load failed and pKaLearn checkpoint/lineage is unresolved."},
            {"reason": "MULTI_PKA_SEMANTICS_UNRESOLVED", "count": "not accepted", "detail": "No record-level source was available to distinguish macro/basic/acidic/site labels reproducibly."},
        ],
    )
    logd = dataset(
        "LOGD_7_4", "logD", "experimental octanol-water distribution coefficient at explicit pH 7.4",
        ["ChEMBL", "PubChem", "logD74 curated dataset", "OpenADMET-related sources"],
        [
            {"reason": "LICENSE_OR_REUSE_UNRESOLVED", "count": 1130, "detail": "The discovered logD74 candidate dataset has pH-specific records but its reusable data license is unresolved."},
            {"reason": "LOGP_NOT_LOGD74", "count": "all unrelated", "detail": "cLogP/logP is never relabeled as experimental logD7.4."},
            {"reason": "DERIVED_NOT_EXPERIMENTAL", "count": "all Drug-OPT derived outputs", "detail": "Henderson-Hasselbalch outputs remain derived estimates only."},
        ],
    )
    (ROOT / "validation/pka_experimental_dataset_v1.json").write_text(json.dumps(pka, indent=2, sort_keys=True) + "\n")
    (ROOT / "validation/logd74_experimental_dataset_v1.json").write_text(json.dumps(logd, indent=2, sort_keys=True) + "\n")

    matrix = [
        ["SOLUBILITY_GENERIC", "AqSolDB/production route", "PRODUCTION_STABLE", "PRODUCTION_STABLE", "validated core route", "none"],
        ["CACO2_PAPP_AB", "TDC/OpenADMET route", "PRODUCTION_STABLE", "PRODUCTION_STABLE", "validated core route", "none"],
        ["HUMAN_PPB", "existing PPB route", "RESEARCH_CANDIDATE", "INDEPENDENT_VALIDATION_REQUIRED", "DrugBank/reference CV only", "independent human plasma fu"],
        ["HLM_CLINT", "existing HLM route", "PRODUCTION_STABLE", "PRODUCTION_STABLE_WITH_LIMITATIONS", "assay/context limitations", "independent context-rich validation"],
        ["RLM_CLINT", "existing RLM route", "PRODUCTION_STABLE", "PRODUCTION_STABLE_WITH_LIMITATIONS", "species-specific route", "none"],
        ["MLM_CLINT", "existing MLM route", "PRODUCTION_STABLE", "PRODUCTION_STABLE_WITH_LIMITATIONS", "species-specific route", "none"],
        ["PKA", "ionization_smarts_rules_v1", "L1", "DATA_LIMITED", "no qualified experimental dataset", "licensed macro/site pKa data"],
        ["LOGD_7_4", "henderson_hasselbalch_logd_v1", "L1", "DATA_LIMITED", "derived, not experimental ML", "licensed pH 7.4 logD data"],
        ["VDSS", "mechanistic consensus", "L3", "DATA_LIMITED", "strict source-qualified N=50", "independent true Vss data"],
        ["HUMAN_TOTAL_IV_CL", "Random Forest + RDKit descriptors", "RESEARCH_CANDIDATE", "INDEPENDENT_VALIDATION_REQUIRED", "same-source validation only", "independent human IV total CL"],
        ["QUANTITATIVE_CYP_2C19", "none", "MODEL_UNAVAILABLE", "MODEL_UNAVAILABLE", "no validated checkpoint", "validated human quantitative dataset/model"],
        ["QUANTITATIVE_PGP", "none", "MODEL_UNAVAILABLE", "MODEL_UNAVAILABLE", "no validated checkpoint", "validated quantitative transporter data/model"],
        ["QUANTITATIVE_BCRP", "none", "MODEL_UNAVAILABLE", "MODEL_UNAVAILABLE", "no validated checkpoint", "validated quantitative transporter data/model"],
    ]
    status_rows = [
        {"endpoint": a, "production_model": b, "maturity": c, "status": d, "limitation": e, "next_requirement": f,
         "ad": "endpoint-specific AD or explicit fail-closed status"}
        for a, b, c, d, e, f in matrix
    ]
    baseline = {
        "artifact": "GLOBAL_DEVELOPABILITY_V1_BASELINE",
        "created_at": TODAY,
        "reference_library": {"artifact": "REFERENCE_LIBRARY_V1", "n": 1000, "sha256": sha("validation/reference_library_v1_1000.json")},
        "production_engine": {"engine": "drugopt-prediction-engine-v1", "status": "UNCHANGED", "historical_runs_modified": 0},
        "endpoint_registry": {"source": "backend/canonical_endpoints.py", "sha256": sha("backend/canonical_endpoints.py")},
        "evidence_warehouse": {"schema": "canonical evidence/provenance contracts", "source_family_separation": True, "fail_closed": True},
        "model_status_matrix": status_rows,
        "candidate_models": {"PPB_FU": "RESEARCH_CANDIDATE", "HUMAN_TOTAL_IV_CL": "RESEARCH_CANDIDATE", "VDSS": "DATA_LIMITED", "PKA": "DATA_LIMITED", "LOGD_7_4": "DATA_LIMITED"},
        "licensed_data_rebuild_trigger": ["NEW_LICENSED_DATA", "canonical qualification", "dataset version increment", "leakage-safe rebuild", "candidate comparison", "independent validation", "production review"],
        "project_adapter_contract": ["global_prediction", "global_uncertainty", "AD", "model_version", "model_status", "residual-compatible output"],
        "known_limitations": ["PPB/fu independent validation absent", "human total IV CL independent validation absent", "strict Vdss source diversity absent", "pKa experimental data insufficient", "logD7.4 licensed experimental data unavailable"],
        "readiness": "GLOBAL_DEVELOPABILITY_V1_BASELINE_READY_WITH_LIMITATIONS",
    }
    baseline["baseline_hash"] = hashlib.sha256(json.dumps(baseline, sort_keys=True).encode()).hexdigest()
    (ROOT / "validation/global_developability_v1_baseline.json").write_text(json.dumps(baseline, indent=2, sort_keys=True) + "\n")
    result = {
        "artifact": "global_developability_v1_closure_v6_4",
        "created_at": TODAY,
        "pka": {"qualified_experimental_n": 0, "acid_n": 0, "base_n": 0, "best_model": "none", "validation_n": 0, "status": "DATA_LIMITED", "reason": "No qualified experimental label set; rule estimate remains L1."},
        "logd74": {"qualified_experimental_n": 0, "source_candidate_n": 1130, "best_model": "none", "validation_n": 0, "status": "DATA_LIMITED", "reason": "pH-specific candidate reuse/license unresolved; derived logD remains L1."},
        "existing_pk_candidates": {"PPB_FU": "RESEARCH_CANDIDATE", "HUMAN_TOTAL_IV_CL": "RESEARCH_CANDIDATE", "VDSS": "DATA_LIMITED"},
        "global_status": baseline["readiness"],
        "baseline_artifact": "global_developability_v1_baseline.json",
        "engine_version_change": False,
        "production": "UNCHANGED",
    }
    (ROOT / "validation/global_developability_v1_closure_v6_4.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"pka_experimental_n": 0, "logd74_experimental_n": 0, "logd74_candidate_n": 1130, "global_status": result["global_status"]}, indent=2))


if __name__ == "__main__":
    main()
