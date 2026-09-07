"""Leakage-safe validation-mode contracts for end-to-end PK evaluation."""
from __future__ import annotations

from enum import Enum
from typing import Any


class PKValidationMode(str, Enum):
    OBSERVED_PARAMETER_ASSISTED = "OBSERVED_PARAMETER_ASSISTED"
    HYBRID_PREDICTION = "HYBRID_PREDICTION"
    FULL_END_TO_END_PREDICTION = "FULL_END_TO_END_PREDICTION"


CRITICAL_UPSTREAM = ("pKa", "logD7.4", "fu", "HLM_Clint", "CL", "VDss", "F", "ka")

# These are intentionally validation contracts, not a promise that every
# production route can currently populate them.  Keeping the profiles
# explicit prevents an assisted run from being relabelled as end-to-end.
HYBRID_PROFILES = {
    "HYBRID_CL_TEST": ("fu", "HLM_Clint", "dose", "route"),
    "HYBRID_VDSS_TEST": ("VDss", "dose", "route"),
    "HYBRID_ABSORPTION_TEST": ("F", "ka", "dose", "route"),
    "HYBRID_FULL_PK": CRITICAL_UPSTREAM + ("dose", "route"),
}


def validate_mode_inputs(mode: PKValidationMode | str, inputs: dict[str, Any]) -> dict:
    """Validate provenance gates; missing values become incomplete, never defaults."""
    mode = PKValidationMode(mode)
    provenance = inputs.get("provenance", {}) or {}
    missing = [name for name in CRITICAL_UPSTREAM if inputs.get(name) is None]
    observed = [name for name in CRITICAL_UPSTREAM if provenance.get(name) in {"EXPERIMENTAL", "OBSERVED_INPUT"}]
    forbidden = []
    if mode == PKValidationMode.FULL_END_TO_END_PREDICTION:
        forbidden = observed
    complete = not missing and not forbidden
    return {
        "mode": mode.value,
        "complete": complete,
        "status": "COMPLETE" if complete else "END_TO_END_INCOMPLETE",
        "missing_parameters": missing,
        "observed_inputs": observed,
        "forbidden_observed_inputs": forbidden,
        "target_leakage": any(name in {"CL", "F", "VDss"} for name in forbidden),
    }


def classify_compound_modes(*, assisted_inputs: dict, hybrid_inputs: dict, full_inputs: dict) -> dict:
    return {"assisted": validate_mode_inputs(PKValidationMode.OBSERVED_PARAMETER_ASSISTED, assisted_inputs),
            "hybrid": validate_mode_inputs(PKValidationMode.HYBRID_PREDICTION, hybrid_inputs),
            "full_end_to_end": validate_mode_inputs(PKValidationMode.FULL_END_TO_END_PREDICTION, full_inputs)}


def endpoint_completeness(*, mode: PKValidationMode | str, outputs: dict[str, object],
                          eligible: bool = True) -> dict[str, object]:
    """Return an endpoint-level result without hiding a partial profile.

    A missing Tmax, for example, must not erase an otherwise valid AUC result.
    ``outputs`` is expected to contain only outputs produced by the declared
    mode; observed targets are rejected for hybrid/full validation.
    """
    mode = PKValidationMode(mode)
    result = {}
    for endpoint, value in outputs.items():
        observed = isinstance(value, dict) and value.get("provenance") in {"EXPERIMENTAL", "OBSERVED_INPUT"}
        result[endpoint] = {
            "complete": bool(eligible and value is not None and not (mode != PKValidationMode.OBSERVED_PARAMETER_ASSISTED and observed)),
            "provenance": value.get("provenance") if isinstance(value, dict) else None,
        }
    return {"mode": mode.value, "eligible": eligible, "endpoints": result,
            "complete_n": sum(1 for row in result.values() if row["complete"])}
