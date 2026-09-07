"""Leakage-safe validation-mode contracts for end-to-end PK evaluation."""
from __future__ import annotations

from enum import Enum
from typing import Any


class PKValidationMode(str, Enum):
    OBSERVED_PARAMETER_ASSISTED = "OBSERVED_PARAMETER_ASSISTED"
    HYBRID_PREDICTION = "HYBRID_PREDICTION"
    FULL_END_TO_END_PREDICTION = "FULL_END_TO_END_PREDICTION"


CRITICAL_UPSTREAM = ("pKa", "logD7.4", "fu", "HLM_Clint", "CL", "VDss", "F", "ka")


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
