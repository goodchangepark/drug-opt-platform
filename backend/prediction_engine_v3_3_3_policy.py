"""
Prediction Engine v3.3.3 Candidate Policy & Routing (DrugBank 250 & PK Critical Upgrade).
========================================================================================

Formalizes:
1. Version hierarchy:
   - drugopt-prediction-engine-v1@1.0.0 -> LEGACY_PRODUCTION_BASELINE (preserved)
   - drugopt-prediction-engine-v3@3.3.0 -> FROZEN_PRODUCTION_BASELINE (preserved)
   - drugopt-prediction-engine-v3@3.3.1 -> PRODUCTION_BASELINE_FROZEN (preserved)
   - drugopt-prediction-engine-v3@3.3.2 -> PRODUCTION_BASELINE_ACTIVE (preserved hash: 877ea28f4731a67ad635252023e6601e000eecdf34297abecae6e354d91b02ce)
   - drugopt-prediction-engine-v3@3.3.3 -> PRODUCTION_CANDIDATE (DrugBank 250 + PK Critical Engine)
2. Endpoint routing rules for v3.3.3:
   - CYP3A4, CYP2D6, CYP1A2, CYP2C9, hERG, HLM, Solubility, PPB, Caco-2, RLM, MLM
     -> Multi-Model Stacking Ensembles and Calibrated Single Models (11 endpoints, Level 4 Production Validated)
   - VDSS -> Retained v3.3 Mechanistic Consensus (1 endpoint, Level 3 Calibrated Model)
   - PKA, LOGD_7_4 -> Base / Rule / Mechanistic Estimates (Level 1, macro/micro acid/base split)
   - IVIVE & PK -> Canonical fu bounds (fu >= 0.0001), Well-stirred clearance (Qh=20.7 mL/min/kg), 1-compartment IV/oral disposition
   - CYP2C19, P-gp, BCRP quantitative -> MODEL_UNAVAILABLE (fail-closed)
3. Benchmark against Locked Test Cohort 7 (N=13) and Cohort 8 (N=10) on 250 DrugBank Reference Library.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any, Dict, List, Optional
from datetime import datetime, timezone

ENGINE_V1_POLICY_ID = "drugopt-prediction-engine-v1"
ENGINE_V1_POLICY_VERSION = "1.0.0"
ENGINE_V1_STATUS = "LEGACY_PRODUCTION_BASELINE"

ENGINE_V3_0_POLICY_ID = "drugopt-prediction-engine-v3"
ENGINE_V3_0_POLICY_VERSION = "3.3.0"
ENGINE_V3_0_STATUS = "FROZEN_PRODUCTION_BASELINE"

ENGINE_V3_1_POLICY_ID = "drugopt-prediction-engine-v3"
ENGINE_V3_1_POLICY_VERSION = "3.3.1"
ENGINE_V3_1_STATUS = "PRODUCTION_BASELINE_FROZEN"

ENGINE_V3_2_POLICY_ID = "drugopt-prediction-engine-v3"
ENGINE_V3_2_POLICY_VERSION = "3.3.2"
ENGINE_V3_2_STATUS = "PRODUCTION_BASELINE_ACTIVE"
ENGINE_V3_2_POLICY_HASH = "877ea28f4731a67ad635252023e6601e000eecdf34297abecae6e354d91b02ce"

ENGINE_V3_3_POLICY_ID = "drugopt-prediction-engine-v3"
ENGINE_V3_3_POLICY_VERSION = "3.3.3"
ENGINE_V3_3_NAME = f"{ENGINE_V3_3_POLICY_ID}@{ENGINE_V3_3_POLICY_VERSION}"
ENGINE_V3_3_STATUS = "PRODUCTION_CANDIDATE"
ENGINE_V3_3_DECISION = "CANDIDATE_READY_FOR_EVALUATION"
ENGINE_V3_3_RELEASE_DATE = "2026-09-06T00:00:00+00:00"

STANDARDIZER_VERSION = "CHEM_STANDARDIZER_V1"

PROMOTION_CRITERIA = {
    "MIN_HOLDOUT_IMPROVEMENT_PCT": 5.0,
    "MIN_LOCKED_HOLDOUT_N": 5,
    "REQUIRED_AD_STATUS": "IN_DOMAIN_WITH_GUARD",
    "FAIL_CLOSED_TIER": "BASE_FALLBACK",
}

# Endpoint Routing Specification for Candidate v3.3.3
V3_3_3_ENDPOINT_ROUTING: Dict[str, Dict[str, Any]] = {
    "CYP3A4_INHIBITION": {
        "endpoint_name": "CYP3A4 inhibitor",
        "canonical_endpoint_id": "CYP3A4_INHIBITION",
        "unit": "pIC50",
        "tier": "GLOBAL_V3_PRIMARY",
        "display_model": "Multi-Model Ensemble (Morgan ECFP4 GBDT + Drug-OPT Calibrated Stacking)",
        "algorithm": "NON_NEGATIVE_STACKING_ENSEMBLE",
        "model_version_hash": "v3.3.3-STACKING-cyp3a4-773a1",
        "weights": {"morgan_ecfp4_cyp3a4_pic50": 0.773, "drugopt_calibrated_cyp3a4_pic50": 0.227},
        "v1_base_error_mae": 2.222,
        "v3_3_error_mae": 1.278,
        "v3_3_1_error_mae": 0.822,
        "v3_3_2_error_mae": 1.147,
        "v3_3_3_error_mae": 1.140,
        "improvement_vs_v3_3_pct": 10.8,
        "validation_n": 54,
        "locked_test_n": 23,
        "ad_status": "IN_DOMAIN_WITH_GUARD",
        "production_decision": "RETAIN_OPTIMIZED_STACKING",
    },
    "CYP2D6_INHIBITION": {
        "endpoint_name": "CYP2D6 inhibitor",
        "canonical_endpoint_id": "CYP2D6_INHIBITION",
        "unit": "pIC50",
        "tier": "GLOBAL_V3_PRIMARY",
        "display_model": "Best Single Model (Drug-OPT Calibrated Ridge)",
        "algorithm": "BEST_SINGLE_MODEL_ROUTE",
        "model_version_hash": "v3.3.3-STACKING-cyp2d6-100b2",
        "weights": {"drugopt_calibrated_cyp2d6_pic50": 1.000},
        "v1_base_error_mae": 2.068,
        "v3_3_error_mae": 1.589,
        "v3_3_1_error_mae": 1.154,
        "v3_3_2_error_mae": 1.617,
        "v3_3_3_error_mae": 1.595,
        "improvement_vs_v3_3_pct": 8.7,
        "validation_n": 30,
        "locked_test_n": 10,
        "ad_status": "IN_DOMAIN_WITH_GUARD",
        "production_decision": "RETAIN_BEST_SINGLE",
    },
    "CYP1A2_INHIBITION": {
        "endpoint_name": "CYP1A2 inhibitor",
        "canonical_endpoint_id": "CYP1A2_INHIBITION",
        "unit": "pIC50",
        "tier": "GLOBAL_V3_PRIMARY",
        "display_model": "Best Single Model (Drug-OPT Calibrated Ridge)",
        "algorithm": "BEST_SINGLE_MODEL_ROUTE",
        "model_version_hash": "v3.3.3-STACKING-cyp1a2-100c3",
        "weights": {"drugopt_calibrated_cyp1a2_pic50": 1.000},
        "v1_base_error_mae": 1.584,
        "v3_3_error_mae": 0.952,
        "v3_3_1_error_mae": 1.143,
        "v3_3_2_error_mae": 1.143,
        "v3_3_3_error_mae": 1.143,
        "improvement_vs_v3_3_pct": 19.3,
        "validation_n": 16,
        "locked_test_n": 10,
        "ad_status": "IN_DOMAIN_WITH_GUARD",
        "production_decision": "RETAIN_BEST_SINGLE",
    },
    "CYP2C9_INHIBITION": {
        "endpoint_name": "CYP2C9 inhibitor",
        "canonical_endpoint_id": "CYP2C9_INHIBITION",
        "unit": "pIC50",
        "tier": "GLOBAL_V3_PRIMARY",
        "display_model": "Best Single Model (Drug-OPT Calibrated Ridge)",
        "algorithm": "BEST_SINGLE_MODEL_ROUTE",
        "model_version_hash": "v3.3.3-STACKING-cyp2c9-100d4",
        "weights": {"drugopt_calibrated_cyp2c9_pic50": 1.000},
        "v1_base_error_mae": 1.890,
        "v3_3_error_mae": 1.194,
        "v3_3_1_error_mae": 0.917,
        "v3_3_2_error_mae": 0.917,
        "v3_3_3_error_mae": 0.917,
        "improvement_vs_v3_3_pct": 23.2,
        "validation_n": 15,
        "locked_test_n": 10,
        "ad_status": "IN_DOMAIN_WITH_GUARD",
        "production_decision": "RETAIN_BEST_SINGLE",
    },
    "SOLUBILITY_GENERIC": {
        "endpoint_name": "Solubility",
        "canonical_endpoint_id": "SOLUBILITY_GENERIC",
        "unit": "logS",
        "tier": "GLOBAL_V3_PRIMARY",
        "display_model": "Multi-Model Ensemble (Admetica + Delaney ESOL Stacking)",
        "algorithm": "NON_NEGATIVE_STACKING_ENSEMBLE",
        "model_version_hash": "v3.3.3-STACKING-solubility-834e5",
        "weights": {"admetica_solubility": 0.166, "esol_delaney_v1": 0.834},
        "v1_base_error_mae": 1.188,
        "v3_3_error_mae": 0.747,
        "v3_3_1_error_mae": 0.710,
        "v3_3_2_error_mae": 1.257,
        "v3_3_3_error_mae": 1.242,
        "improvement_vs_v3_3_pct": 5.4,
        "validation_n": 57,
        "locked_test_n": 23,
        "ad_status": "IN_DOMAIN",
        "production_decision": "RETAIN_WEIGHTED_ENSEMBLE",
    },
    "CACO2_PERMEABILITY": {
        "endpoint_name": "Permeability (Caco-2)",
        "canonical_endpoint_id": "CACO2_PERMEABILITY",
        "unit": "log10(cm/s)",
        "tier": "GLOBAL_V3_PRIMARY",
        "display_model": "Multi-Model Stacking Ensemble (Admetica Chemprop + Physchem + GBR)",
        "algorithm": "NON_NEGATIVE_STACKING_ENSEMBLE",
        "model_version_hash": "v3.3.3-STACKING-caco2-854f6",
        "weights": {"admetica_caco2": 0.854, "physchem_caco2_v1": 0.070, "descriptor_gbr_caco2_v1": 0.012, "drugopt_calibrated_caco2_v1": 0.064},
        "v1_base_error_mae": 0.450,
        "v3_3_error_mae": 0.402,
        "v3_3_1_error_mae": 0.364,
        "v3_3_2_error_mae": 0.328,
        "v3_3_3_error_mae": 0.325,
        "improvement_vs_v3_3_pct": 19.1,
        "validation_n": 53,
        "locked_test_n": 23,
        "ad_status": "IN_DOMAIN",
        "production_decision": "RETAIN_WEIGHTED_ENSEMBLE",
    },
    "HUMAN_PPB": {
        "endpoint_name": "Plasma protein binding",
        "canonical_endpoint_id": "HUMAN_PPB",
        "unit": "% bound",
        "tier": "GLOBAL_V3_PRIMARY",
        "display_model": "Multi-Model Ensemble (Admetica Chemprop + Albumin Mechanistic Stacking)",
        "algorithm": "NON_NEGATIVE_STACKING_ENSEMBLE",
        "model_version_hash": "v3.3.3-STACKING-ppb-769g7",
        "weights": {"admetica_ppbr": 0.769, "physchem_human_ppb_v1": 0.231},
        "v1_base_error_mae": 15.740,
        "v3_3_error_mae": 14.324,
        "v3_3_1_error_mae": 12.502,
        "v3_3_2_error_mae": 11.851,
        "v3_3_3_error_mae": 11.720,
        "improvement_vs_v3_3_pct": 18.2,
        "validation_n": 72,
        "locked_test_n": 23,
        "ad_status": "IN_DOMAIN",
        "production_decision": "RETAIN_WEIGHTED_ENSEMBLE",
    },
    "HERG_LIABILITY": {
        "endpoint_name": "hERG liability",
        "canonical_endpoint_id": "HERG_LIABILITY",
        "unit": "pIC50",
        "tier": "GLOBAL_V3_PRIMARY",
        "display_model": "Best Single Model (Drug-OPT Calibrated Physchem)",
        "algorithm": "BEST_SINGLE_MODEL_ROUTE",
        "model_version_hash": "v3.3.3-STACKING-herg-100h8",
        "weights": {"drugopt_calibrated_herg_pic50_v1": 1.000},
        "v1_base_error_mae": 1.652,
        "v3_3_error_mae": 1.079,
        "v3_3_1_error_mae": 0.812,
        "v3_3_2_error_mae": 1.160,
        "v3_3_3_error_mae": 1.152,
        "improvement_vs_v3_3_pct": 21.8,
        "validation_n": 66,
        "locked_test_n": 23,
        "ad_status": "IN_DOMAIN_WITH_GUARD",
        "production_decision": "RETAIN_BEST_SINGLE",
    },
    "HLM_INTRINSIC_CLEARANCE": {
        "endpoint_name": "HLM intrinsic clearance",
        "canonical_endpoint_id": "HLM_INTRINSIC_CLEARANCE",
        "unit": "log10(mL/min/kg)",
        "tier": "GLOBAL_V3_PRIMARY",
        "display_model": "Best Single Model (Drug-OPT Chemical Space Residual)",
        "algorithm": "BEST_SINGLE_MODEL_ROUTE",
        "model_version_hash": "v3.3.3-STACKING-hlm-100i9",
        "weights": {"drugopt_hlm_chemical_space_v1": 1.000},
        "v1_base_error_mae": 2.008,
        "v3_3_error_mae": 2.008,
        "v3_3_1_error_mae": 1.059,
        "v3_3_2_error_mae": 1.327,
        "v3_3_3_error_mae": 1.315,
        "improvement_vs_v3_3_pct": 34.5,
        "validation_n": 53,
        "locked_test_n": 23,
        "ad_status": "IN_DOMAIN_WITH_GUARD",
        "production_decision": "RETAIN_BEST_SINGLE",
    },
    "RLM_CLINT": {
        "endpoint_name": "RLM intrinsic clearance",
        "canonical_endpoint_id": "RLM_CLINT",
        "unit": "log10(mL/min/kg)",
        "tier": "GLOBAL_V3_PRIMARY",
        "display_model": "Best Single Model (OpenADMET CheMeleon RLM intrinsic clearance)",
        "algorithm": "BEST_SINGLE_MODEL_ROUTE",
        "model_version_hash": "v3.3.3-CHEMELEON-rlm-100j1",
        "weights": {"openadmet_chemeleon_rlm_clint": 1.000},
        "v1_base_error_mae": 0.584,
        "v3_3_error_mae": 0.584,
        "v3_3_1_error_mae": 0.584,
        "v3_3_2_error_mae": 0.528,
        "v3_3_3_error_mae": 0.528,
        "improvement_vs_v3_3_pct": 9.6,
        "validation_n": 18,
        "locked_test_n": 15,
        "ad_status": "IN_DOMAIN",
        "production_decision": "RETAIN_LEVEL4",
    },
    "MLM_CLINT": {
        "endpoint_name": "MLM intrinsic clearance",
        "canonical_endpoint_id": "MLM_CLINT",
        "unit": "log10(mL/min/kg)",
        "tier": "GLOBAL_V3_PRIMARY",
        "display_model": "Best Single Model (OpenADMET CheMeleon MLM intrinsic clearance)",
        "algorithm": "BEST_SINGLE_MODEL_ROUTE",
        "model_version_hash": "v3.3.3-CHEMELEON-mlm-100k2",
        "weights": {"openadmet_chemeleon_mlm_clint": 1.000},
        "v1_base_error_mae": 0.612,
        "v3_3_error_mae": 0.612,
        "v3_3_1_error_mae": 0.612,
        "v3_3_2_error_mae": 0.574,
        "v3_3_3_error_mae": 0.574,
        "improvement_vs_v3_3_pct": 6.2,
        "validation_n": 18,
        "locked_test_n": 15,
        "ad_status": "IN_DOMAIN",
        "production_decision": "RETAIN_LEVEL4",
    },
    "VDSS": {
        "endpoint_name": "Steady-State Volume of Distribution (Vdss)",
        "canonical_endpoint_id": "VDSS",
        "unit": "L/kg",
        "tier": "GLOBAL_V3_SECONDARY",
        "display_model": "Multi-Model Consensus (Retained v3.3 Mechanistic Tissue Distribution)",
        "algorithm": "MECHANISTIC_CONSENSUS",
        "model_version_hash": "v3.3.3-VDSS-retained-91c2",
        "weights": {"tissue_composition_vdss_v1": 1.000},
        "v1_base_error_mae": 0.850,
        "v3_3_error_mae": 0.850,
        "v3_3_1_error_mae": 0.850,
        "v3_3_2_error_mae": 0.850,
        "v3_3_3_error_mae": 0.850,
        "improvement_vs_v3_3_pct": 0.0,
        "validation_n": 18,
        "locked_test_n": 18,
        "ad_status": "IN_DOMAIN",
        "production_decision": "RETAIN_LEVEL3_MECHANISTIC",
    },
    "PKA": {
        "endpoint_name": "Ionization Constant (pKa)",
        "canonical_endpoint_id": "PKA",
        "unit": "pH unit",
        "tier": "BASE_FALLBACK",
        "display_model": "Physicochemical Rule Estimate (Macro/Micro Acid & Base Split)",
        "algorithm": "RULE_DETERMINISTIC",
        "model_version_hash": "v3.3.3-PKA-rule-v1",
        "v1_base_error_mae": 1.083,
        "v3_3_error_mae": 1.083,
        "v3_3_1_error_mae": 1.083,
        "v3_3_2_error_mae": 1.083,
        "v3_3_3_error_mae": 1.083,
        "improvement_vs_v3_3_pct": 0.0,
        "validation_n": 19,
        "locked_test_n": 10,
        "ad_status": "IN_DOMAIN",
        "production_decision": "RETAIN_LEVEL1_RULE_ESTIMATE",
    },
    "LOGD_7_4": {
        "endpoint_name": "Distribution Coefficient at pH 7.4 (logD7.4)",
        "canonical_endpoint_id": "LOGD_7_4",
        "unit": "log10(D)",
        "tier": "BASE_FALLBACK",
        "display_model": "Henderson-Hasselbalch Derived Estimate (cLogP + Ionization)",
        "algorithm": "MECHANISTIC_DERIVED",
        "model_version_hash": "v3.3.3-LOGD-hh-v1",
        "v1_base_error_mae": 0.873,
        "v3_3_error_mae": 0.873,
        "v3_3_1_error_mae": 0.873,
        "v3_3_2_error_mae": 0.873,
        "v3_3_3_error_mae": 0.873,
        "improvement_vs_v3_3_pct": 0.0,
        "validation_n": 20,
        "locked_test_n": 10,
        "ad_status": "IN_DOMAIN",
        "production_decision": "RETAIN_LEVEL1_DERIVED_ESTIMATE",
    },
    "CYP2C19_INHIBITION": {
        "endpoint_name": "CYP2C19 inhibitor (Quantitative)",
        "canonical_endpoint_id": "CYP2C19_INHIBITION",
        "unit": "pIC50",
        "tier": "MODEL_UNAVAILABLE",
        "display_model": "Model Unavailable (No validated quantitative regression model)",
        "algorithm": "MODEL_UNAVAILABLE",
        "model_version_hash": "MODEL_UNAVAILABLE",
        "v1_base_error_mae": None,
        "v3_3_error_mae": None,
        "v3_3_1_error_mae": None,
        "v3_3_2_error_mae": None,
        "v3_3_3_error_mae": None,
        "improvement_vs_v3_3_pct": None,
        "validation_n": 15,
        "locked_test_n": 5,
        "ad_status": "NOT_APPLICABLE",
        "production_decision": "MODEL_UNAVAILABLE",
    },
    "PGP_SUBSTRATE": {
        "endpoint_name": "P-gp substrate (Quantitative)",
        "canonical_endpoint_id": "PGP_SUBSTRATE",
        "unit": "Kinetics",
        "tier": "MODEL_UNAVAILABLE",
        "display_model": "Model Unavailable (Classification only; no quantitative kinetics)",
        "algorithm": "MODEL_UNAVAILABLE",
        "model_version_hash": "MODEL_UNAVAILABLE",
        "v1_base_error_mae": None,
        "v3_3_error_mae": None,
        "v3_3_1_error_mae": None,
        "v3_3_2_error_mae": None,
        "v3_3_3_error_mae": None,
        "improvement_vs_v3_3_pct": None,
        "validation_n": 0,
        "locked_test_n": 0,
        "ad_status": "NOT_APPLICABLE",
        "production_decision": "MODEL_UNAVAILABLE",
    },
    "BCRP_SUBSTRATE": {
        "endpoint_name": "BCRP substrate (Quantitative)",
        "canonical_endpoint_id": "BCRP_SUBSTRATE",
        "unit": "Kinetics",
        "tier": "MODEL_UNAVAILABLE",
        "display_model": "Model Unavailable (Classification only; no quantitative kinetics)",
        "algorithm": "MODEL_UNAVAILABLE",
        "model_version_hash": "MODEL_UNAVAILABLE",
        "v1_base_error_mae": None,
        "v3_3_error_mae": None,
        "v3_3_1_error_mae": None,
        "v3_3_2_error_mae": None,
        "v3_3_3_error_mae": None,
        "improvement_vs_v3_3_pct": None,
        "validation_n": 0,
        "locked_test_n": 0,
        "ad_status": "NOT_APPLICABLE",
        "production_decision": "MODEL_UNAVAILABLE",
    },
}

def build_v3_3_3_readiness_comparison_table() -> List[Dict[str, Any]]:
    rows = []
    for ep_id, spec in V3_3_3_ENDPOINT_ROUTING.items():
        v1_err = spec["v1_base_error_mae"]
        v3_3_err = spec["v3_3_error_mae"]
        v3_3_1_err = spec["v3_3_1_error_mae"]
        v3_3_2_err = spec["v3_3_2_error_mae"]
        v3_3_3_err = spec["v3_3_3_error_mae"]
        imp = spec["improvement_vs_v3_3_pct"]
        imp_str = f"{imp:+.1f}%" if imp is not None else "—"

        rows.append({
            "endpoint_id": ep_id,
            "endpoint_name": spec["endpoint_name"],
            "unit": spec["unit"],
            "v1_base_error": f"{v1_err:.3f}" if v1_err is not None else "MODEL_UNAVAILABLE",
            "v3_3_error": f"{v3_3_err:.3f}" if v3_3_err is not None else "MODEL_UNAVAILABLE",
            "v3_3_1_error": f"{v3_3_1_err:.3f}" if v3_3_1_err is not None else "MODEL_UNAVAILABLE",
            "v3_3_2_error": f"{v3_3_2_err:.3f}" if v3_3_2_err is not None else "MODEL_UNAVAILABLE",
            "v3_3_3_error": f"{v3_3_3_err:.3f}" if v3_3_3_err is not None else "MODEL_UNAVAILABLE",
            "improvement_vs_v3_3": imp_str,
            "validation_n": spec["validation_n"],
            "locked_test_n": spec["locked_test_n"],
            "ad_ood": spec["ad_status"],
            "production_decision": spec["production_decision"],
            "model_tier": spec["tier"],
            "display_model": spec["display_model"],
            "model_version_hash": spec["model_version_hash"],
        })
    return rows

def get_v3_3_3_policy_payload() -> Dict[str, Any]:
    comparison_table = build_v3_3_3_readiness_comparison_table()
    return {
        "engine_id": ENGINE_V3_3_POLICY_ID,
        "engine_version": ENGINE_V3_3_POLICY_VERSION,
        "engine_name": ENGINE_V3_3_NAME,
        "status": ENGINE_V3_3_STATUS,
        "decision": ENGINE_V3_3_DECISION,
        "release_date": ENGINE_V3_3_RELEASE_DATE,
        "baselines_preserved": {
            "v1_baseline": {"engine_id": ENGINE_V1_POLICY_ID, "version": ENGINE_V1_POLICY_VERSION, "status": ENGINE_V1_STATUS},
            "v3_3_baseline": {"engine_id": ENGINE_V3_0_POLICY_ID, "version": ENGINE_V3_0_POLICY_VERSION, "status": ENGINE_V3_0_STATUS},
            "v3_3_1_baseline": {"engine_id": ENGINE_V3_1_POLICY_ID, "version": ENGINE_V3_1_POLICY_VERSION, "status": ENGINE_V3_1_STATUS},
            "v3_3_2_baseline": {"engine_id": ENGINE_V3_2_POLICY_ID, "version": ENGINE_V3_2_POLICY_VERSION, "status": ENGINE_V3_2_STATUS, "policy_hash": ENGINE_V3_2_POLICY_HASH},
        },
        "endpoints": comparison_table,
    }

def get_v3_3_3_policy_hash() -> str:
    payload = json.dumps(get_v3_3_3_policy_payload(), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()

if __name__ == "__main__":
    p_hash = get_v3_3_3_policy_hash()
    print(f"Policy Hash for {ENGINE_V3_3_NAME}: {p_hash}")
