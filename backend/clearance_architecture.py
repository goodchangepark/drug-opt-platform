"""Explicit clearance semantics and fail-closed systemic clearance assembly.

This is a metadata/assembly layer.  It does not change the frozen prediction
engine or the existing hepatic IVIVE equations.  In particular, a missing
renal or extra-hepatic component is never silently treated as zero.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import re
from typing import Any

CLEARANCE_ARCHITECTURE_VERSION = "drugopt-clearance-architecture-v1"

CL_TOTAL_IV = "CL_TOTAL_IV"
CL_HEPATIC = "CL_HEPATIC"
CL_RENAL = "CL_RENAL"
CL_NONRENAL = "CL_NONRENAL"
CL_ORAL_APPARENT = "CL_ORAL_APPARENT"
CL_OVER_F = "CL_OVER_F"
CL_UNRESOLVED = "CL_UNRESOLVED"

TOTAL_CL_COMPLETE = "TOTAL_CL_COMPLETE"
TOTAL_CL_INCOMPLETE = "TOTAL_CL_INCOMPLETE"


def _text(value: Any) -> str:
    return str(value or "").strip().upper().replace("–", "-")


def classify_clearance_semantics(endpoint: Any = "", route: Any = "", *, context: dict | None = None) -> str:
    """Classify a clearance observation without guessing unresolved semantics."""
    text = _text(endpoint) + " " + _text((context or {}).get("parameter_subtype"))
    route_text = _text(route or (context or {}).get("route"))
    if re.search(r"CL\s*/\s*F|CLF|APPARENT|OVER[_ ]?F", text):
        return CL_ORAL_APPARENT if route_text in {"", "PO", "ORAL", "P.O."} else CL_OVER_F
    if re.search(r"RENAL|URINARY|KIDNEY|EXCRET(?:ED|ION)", text):
        return CL_RENAL
    if re.search(r"NON[- ]?RENAL|EXTRA[- ]?HEPATIC|OTHER", text):
        return CL_NONRENAL
    if re.search(r"HLM|HEPATIC|LIVER|CLINT|IVIVE", text):
        return CL_HEPATIC
    if route_text in {"IV", "INTRAVENOUS", "BOLUS", "INFUSION"}:
        return CL_TOTAL_IV
    return CL_UNRESOLVED


@dataclass(frozen=True)
class ClearanceComponent:
    name: str
    value: float | None
    unit: str
    source: str
    provenance: dict


def assemble_total_clearance(*, hepatic: ClearanceComponent | None = None,
                             renal: ClearanceComponent | None = None,
                             other: ClearanceComponent | None = None) -> dict:
    """Return total CL only when every required component is supported."""
    components = {"CL_H": hepatic, "CL_R": renal, "CL_OTHER": other}
    missing = [name for name, component in components.items()
               if component is None or component.value is None]
    result = {
        "version": CLEARANCE_ARCHITECTURE_VERSION,
        "status": TOTAL_CL_INCOMPLETE if missing else TOTAL_CL_COMPLETE,
        "value": None if missing else sum(float(c.value) for c in components.values()),
        "unit": "mL/min/kg",
        "missing_components": missing,
        "components": {name: (asdict(component) if component else None)
                        for name, component in components.items()},
        "assumptions": [],
    }
    if missing:
        result["reason"] = "Renal and/or other clearance contribution is unresolved; no zero substitution was made."
    return result


def renal_readiness(*, fraction_excreted_unchanged: float | None = None,
                    renal_clearance: float | None = None,
                    urinary_recovery: float | None = None,
                    gfr: float | None = None,
                    secretion: float | None = None,
                    reabsorption: float | None = None) -> dict:
    """Describe renal evidence; never turn absent renal evidence into CL_R=0."""
    fe = fraction_excreted_unchanged
    evidence = any(x is not None for x in (renal_clearance, urinary_recovery, gfr, secretion, reabsorption, fe))
    if not evidence:
        status = "RENAL_UNRESOLVED"
    elif fe is not None and fe >= 0.5:
        status = "RENAL_DOMINANT"
    elif fe is not None and fe < 0.2:
        status = "LOW_RENAL_CONTRIBUTION"
    else:
        status = "MIXED_CLEARANCE"
    return {
        "status": status,
        "evidence_present": evidence,
        "fraction_excreted_unchanged": fe,
        "renal_clearance": renal_clearance,
        "mechanisms": {"filtration": gfr is not None, "active_secretion": secretion is not None,
                        "reabsorption": reabsorption is not None, "unknown": not evidence},
        "quantitative_prediction_supported": renal_clearance is not None,
        "silent_zero_assumption": False,
    }


def clearance_confidence(*, hepatic_available: bool, renal_status: str,
                         fu_inc_known: bool, rb_known: bool, applicability: str = "IN_DOMAIN") -> str:
    if applicability == "OOD":
        return "OOD_CL"
    if not hepatic_available:
        return "MODEL_UNAVAILABLE"
    if renal_status == "RENAL_UNRESOLVED" or not fu_inc_known or not rb_known:
        return "TOTAL_CL_INCOMPLETE"
    return "HIGH_CONFIDENCE_TOTAL_CL"


def conversion_cl_l_h_kg_to_ml_min_kg(value: float) -> float:
    return float(value) * 1000.0 / 60.0


def conversion_cl_ml_min_kg_to_l_h_kg(value: float) -> float:
    return float(value) * 60.0 / 1000.0
