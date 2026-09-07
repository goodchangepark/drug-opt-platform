"""Leakage-safe predictive hepatic-clearance chain candidate.

This is a candidate PK-layer calibration, not a replacement for the frozen
Prediction Engine.  It accepts model outputs only; observed PPB/HLM values
are intentionally rejected by the interface.
"""
from __future__ import annotations

import math
from typing import Any

PREDICTIVE_HEPATIC_CL_VERSION = "predictive-hepatic-cl-v5.5-candidate"
PPB_DEVELOPMENT_OFFSET_PERCENT = 2.246017699115046
PPB_CALIBRATION_VERSION = "human-ppb-development-offset-v1"
HLM_CALIBRATION_VERSION = "production-hlm-v3.3.2"
HUMAN_QH_ML_MIN_KG = 20.7142857143


def predictive_hepatic_clearance(*, predicted_ppb_percent: float,
                                 predicted_hlm_log10_ml_min_kg: float,
                                 applicability_domain: str,
                                 apply_ppb_development_calibration: bool = True,
                                 qh_ml_min_kg: float = HUMAN_QH_ML_MIN_KG) -> dict[str, Any]:
    """Calculate CLH from predicted PPB and HLM Clint, never observations."""
    for name, value in (("predicted_ppb_percent", predicted_ppb_percent),
                        ("predicted_hlm_log10_ml_min_kg", predicted_hlm_log10_ml_min_kg)):
        if not math.isfinite(float(value)):
            raise ValueError(f"{name} must be finite")
    if applicability_domain not in {"IN_DOMAIN", "BORDERLINE", "OUT_OF_DOMAIN"}:
        raise ValueError("unknown applicability domain")
    ppb = float(predicted_ppb_percent)
    if apply_ppb_development_calibration:
        ppb -= PPB_DEVELOPMENT_OFFSET_PERCENT
    ppb = min(99.0, max(0.0, ppb))
    fu_p = (100.0 - ppb) / 100.0
    # Preserve the production HLM AD guard; this does not use the target CL.
    correction = 0.237 if applicability_domain == "IN_DOMAIN" else (0.5 * 0.237 if applicability_domain == "BORDERLINE" else 0.0)
    hlm_log = float(predicted_hlm_log10_ml_min_kg) - correction
    clint = 10.0 ** hlm_log
    driving = fu_p * clint
    clh = float(qh_ml_min_kg) * driving / (float(qh_ml_min_kg) + driving)
    return {
        "value": clh,
        "unit": "mL/min/kg",
        "prediction_type": "MODEL_PREDICTED_CHAIN",
        "provenance": {
            "chain_version": PREDICTIVE_HEPATIC_CL_VERSION,
            "ppb_source": "MODEL_PREDICTED",
            "ppb_calibration": PPB_CALIBRATION_VERSION if apply_ppb_development_calibration else "NONE",
            "hlm_source": "MODEL_PREDICTED",
            "hlm_calibration": HLM_CALIBRATION_VERSION,
            "ivive": "WELL_STIRRED",
            "physiological_assumptions": {"Qh": float(qh_ml_min_kg), "Qh_unit": "mL/min/kg"},
            "observed_inputs_used": [],
        },
        "predicted_ppb_percent": ppb,
        "predicted_fu_p": fu_p,
        "predicted_hlm_log10_ml_min_kg": hlm_log,
        "predicted_hlm_ml_min_kg": clint,
        "applicability_domain": applicability_domain,
        "confidence": "MEDIUM" if applicability_domain == "IN_DOMAIN" else "LOW",
    }
