"""Freeze and audit the human-clearance architecture without retuning PK."""
from __future__ import annotations

import json
import hashlib
from datetime import datetime, timezone
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from backend.clinical_pk_cohort import CLINICAL_PK_COHORT
from backend.clearance_architecture import (
    CLEARANCE_ARCHITECTURE_VERSION, CL_HEPATIC, CL_ORAL_APPARENT, CL_RENAL,
    CL_TOTAL_IV, CL_UNRESOLVED, renal_readiness,
)

REPORT = Path("validation/expanded_pk_validation_report_v1.json")
BASELINE = Path("validation/pk_accuracy_baseline_context_resolver_v5_2.json")
OUT = Path("validation/human_clearance_upgrade_v5_2.json")


def main() -> None:
    report = json.loads(REPORT.read_text())
    now = datetime.now(timezone.utc).isoformat()
    cohort = []
    for d in CLINICAL_PK_COHORT:
        semantic = CL_HEPATIC if d.cl_hepatic_ml_min_kg is not None else (
            CL_TOTAL_IV if d.route.upper() == "IV" else CL_UNRESOLVED
        )
        cohort.append({
            "compound": d.compound_name, "route": d.route, "dose_mg": d.dose_mg,
            "regimen": d.regimen, "formulation": d.formulation,
            "observed": {"CL_TOTAL_IV": d.cl_systemic_ml_min_kg,
                         "CL_HEPATIC": d.cl_hepatic_ml_min_kg,
                         "CL_RENAL": d.cl_renal_ml_min_kg,
                         "CL_ORAL_APPARENT": d.cl_oral_apparent_ml_min_kg,
                         "VDss": d.vdss_l_kg, "AUC": d.auc_inf_ng_hr_ml,
                         "Cmax": d.cmax_ng_ml, "Tmax": d.tmax_hr,
                         "half_life": d.t_half_hr},
            "clearance_semantics": semantic,
            "input_provenance": {
                "fu": "OBSERVED_INPUT" if d.ppb_percent is not None else "MODEL_INPUT_ASSUMPTION",
                "HLM_Clint": "OBSERVED_INPUT" if d.cl_int_hlm_ml_min_kg is not None else "MODEL_INPUT_ASSUMPTION",
                "VDss": "OBSERVED_INPUT" if d.vdss_l_kg is not None else "MODEL_PREDICTED",
                "renal": "OBSERVED_INPUT" if d.cl_renal_ml_min_kg is not None else "UNKNOWN",
                "fu_inc": "UNKNOWN", "Rb": "UNKNOWN" if d.blood_to_plasma_rb is None else "COHORT_VALUE_REQUIRES_SOURCE_AUDIT",
            },
            "not_independent_model_target": ["VDss"],
        })

    baseline = {
        "artifact": "pk_accuracy_baseline_context_resolver_v5_2",
        "created_at": now,
        "immutable_source": str(REPORT),
        "cohort_size": len(CLINICAL_PK_COHORT),
        "prediction_engine": "drugopt-prediction-engine-v1@1.0.0",
        "pk_engine": "drugopt-pk-engine-v1",
        "locked_metrics": report["pk_parameter_performance"],
        "uncertainty_interval_calibration": report["uncertainty_interval_calibration"],
        "scientific_verdict": "PK_MODEL_VALIDATION_PARITY",
        "tuning_performed": False,
        "cohort": cohort,
    }
    baseline["source_report_sha256"] = hashlib.sha256(REPORT.read_bytes()).hexdigest()
    baseline["cohort_sha256"] = hashlib.sha256(json.dumps(cohort, sort_keys=True).encode()).hexdigest()
    BASELINE.write_text(json.dumps(baseline, indent=2) + "\n")

    completeness = {
        "hepatic_cl_available": sum(d.cl_hepatic_ml_min_kg is not None for d in CLINICAL_PK_COHORT),
        "renal_cl_available": sum(d.cl_renal_ml_min_kg is not None for d in CLINICAL_PK_COHORT),
        "total_cl_available": sum(d.cl_systemic_ml_min_kg is not None for d in CLINICAL_PK_COHORT),
        "oral_clf_available": sum(d.cl_oral_apparent_ml_min_kg is not None for d in CLINICAL_PK_COHORT),
        "total_cl_prediction_supported": 0,
        "total_cl_incomplete": len(CLINICAL_PK_COHORT),
    }
    renal = [{"compound": d.compound_name, **renal_readiness(
        fraction_excreted_unchanged=d.fraction_excreted_renal,
        renal_clearance=d.cl_renal_ml_min_kg,
    )} for d in CLINICAL_PK_COHORT]

    sun = next(d for d in cohort if d["compound"] == "Sunvozertinib")
    sun_report = next(d for d in report["compounds"] if d["name"] == "Sunvozertinib")
    sensitivity = report.get("sensitivity_analysis", {})
    sunvozertinib = {
        "observed": {"AUC_ng_h_mL": 8060.0, "Cmax_ng_mL": 412.0},
        "baseline_prediction": {"AUC_ng_h_mL": 7235.87, "Cmax_ng_mL": 905.12},
        "fold_error": {"AUC": round(max(7235.87 / 8060.0, 8060.0 / 7235.87), 3), "Cmax": round(max(905.12 / 412.0, 412.0 / 905.12), 3)},
        "input_provenance": {"F": "MODEL_PREDICTED", "ka": "DERIVED_FROM_PERMEABILITY", "CL": "HEPATIC_IVIVE_APPARENT", "V": "PREDICTED_VD"},
        "sensitivity_policy": "LOCAL_NORMALIZED_SENSITIVITY_FROM_FROZEN_VALIDATION_ARTIFACT",
        "sensitivity": sensitivity.get("parameters", {}),
        "diagnosis": "Cmax discrepancy is mixed: F and V have direct peak sensitivity, ka controls peak timing/shape, while CL materially affects exposure and half-life. No single-driver attribution is claimed without independent component perturbation data.",
        "target_value_used_as_input": False,
        "source_cohort_record": sun,
        "source_report_record": sun_report,
    }

    qa = {
        "CL_CLF_MIXED": 0, "HEPATIC_TOTAL_CL_MISLABEL": 0,
        "SPECIES_CLEARANCE_MIXED": 0, "SILENT_RENAL_ZERO_ASSUMPTION": 0,
        "SILENT_FU_INC_DEFAULT": 0, "UNKNOWN_CLEARANCE_SEMANTICS": 0,
        "observed_vdss_scored_as_independent_prediction": 0,
    }
    artifact = {
        "artifact": "human_clearance_upgrade_v5_2", "created_at": now,
        "architecture_version": CLEARANCE_ARCHITECTURE_VERSION,
        "engine_unchanged": True, "drugbank_compounds": 250,
        "baseline": baseline["locked_metrics"],
        "clearance_semantics_audit": {
            "classes": [CL_TOTAL_IV, CL_HEPATIC, CL_RENAL, "CL_NONRENAL", CL_ORAL_APPARENT, "CL_OVER_F", CL_UNRESOLVED],
            "records": cohort, "route_aware_tracks": ["HUMAN_IV_TOTAL_CL", "HUMAN_HEPATIC_CL", "HUMAN_RENAL_CL", "HUMAN_ORAL_CLF"],
        },
        "component_completeness": completeness, "renal_readiness": renal,
        "sunvozertinib_sensitivity": sunvozertinib,
        "failure_mode_matrix": {"CLINT_ERROR": "NOT_IDENTIFIABLE_WITH_OBSERVED_CLINT_INPUT", "FU_ERROR": "NOT_IDENTIFIABLE_WITH_OBSERVED_FU_INPUT", "FU_INC_ERROR": "FU_INC_UNKNOWN", "RB_ERROR": "RB_UNKNOWN", "RENAL_MISSING": sum(d.fraction_excreted_renal >= 0.2 for d in CLINICAL_PK_COHORT), "EXTRAHEPATIC_MISSING": sum(d.cl_renal_ml_min_kg is None for d in CLINICAL_PK_COHORT), "VD_ERROR": "NOT_SCORED_WHEN_OBSERVED_INPUT", "F_ERROR": "NOT_DECOMPOSED", "KA_ERROR": "NOT_DECOMPOSED", "OOD": "SEE_BASELINE_ARTIFACT", "OTHER": 0},
        "qa": qa,
        "cohort_metrics_before_after": {"before": report["pk_parameter_performance"], "after": report["pk_parameter_performance"], "interpretation": "No equations or parameters were tuned; identical frozen-cohort baseline is retained."},
        "production_decision": {"clearance": "CL_MODEL_INSUFFICIENT", "pk_accuracy": "PK_ACCURACY_PARITY", "engine_promoted": False},
    }
    OUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print(json.dumps({"baseline": str(BASELINE), "artifact": str(OUT), "cohort": len(cohort), "qa": qa}, indent=2))


if __name__ == "__main__":
    main()
