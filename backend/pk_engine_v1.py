"""Fail-closed PK overlay contract.  This is deliberately separate from Engine v1.

It only produces a transparent one-compartment estimate when every required
input is supplied with provenance.  It is not a validated structure model.
"""
from __future__ import annotations

import hashlib
import json
import math
from typing import Any

PK_ENGINE_VERSION = "drugopt-pk-engine-v1"
MECHANISTIC_ESTIMATE = "MECHANISTIC_ESTIMATE"
DERIVED_ESTIMATE = "DERIVED_ESTIMATE"
INSUFFICIENT_INPUT = "INSUFFICIENT_INPUT"
UNAVAILABLE = "UNAVAILABLE"

REQUIRED_ORAL = ("dose_mg_per_kg", "f_fraction", "ka_per_h", "cl_l_per_h_per_kg", "v_l_per_kg")
REQUIRED_IV = ("dose_mg_per_kg", "cl_l_per_h_per_kg", "v_l_per_kg")


def request_fingerprint(context: dict[str, Any]) -> str:
    """Stable fingerprint; a caller may reuse the immutable completed run."""
    return hashlib.sha256(json.dumps(context, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def _complete(inputs: dict[str, Any], required: tuple[str, ...]) -> list[str]:
    missing=[]
    for name in required:
        value=inputs.get(name)
        source=(inputs.get("sources") or {}).get(name)
        if value is None or not source or source in {"DEFAULT", "SILENT_DEFAULT", "UNAVAILABLE"}:
            missing.append(name)
    return missing


def estimate_one_compartment(*, species: str, route: str, inputs: dict[str, Any]) -> dict[str, Any]:
    """Return no quantitative value unless the explicit model contract holds."""
    route=str(route).upper(); required=REQUIRED_IV if route=="IV" else REQUIRED_ORAL
    missing=_complete(inputs, required)
    base={"pk_engine_version":PK_ENGINE_VERSION,"species":species,"route":route,"inputs":inputs,
          "fingerprint":request_fingerprint({"engine":PK_ENGINE_VERSION,"species":species,"route":route,"inputs":inputs})}
    if missing:
        return base|{"status":INSUFFICIENT_INPUT,"prediction_type":INSUFFICIENT_INPUT,"missing_inputs":missing,"outputs":{},"assumptions":[]}
    dose=float(inputs["dose_mg_per_kg"]); cl=float(inputs["cl_l_per_h_per_kg"]); v=float(inputs["v_l_per_kg"])
    if min(dose,cl,v)<=0 or (route!="IV" and (float(inputs["f_fraction"])<=0 or float(inputs["f_fraction"])>1 or float(inputs["ka_per_h"])<=0)):
        return base|{"status":INSUFFICIENT_INPUT,"prediction_type":INSUFFICIENT_INPUT,"missing_inputs":["positive, dimensionally valid PK inputs"],"outputs":{},"assumptions":[]}
    ke=cl/v; outputs={"cl":{"value":cl,"unit":"L/h/kg","prediction_type":MECHANISTIC_ESTIMATE},"v":{"value":v,"unit":"L/kg","prediction_type":MECHANISTIC_ESTIMATE},"t_half":{"value":math.log(2)/ke,"unit":"h","prediction_type":DERIVED_ESTIMATE}}
    if route=="IV":
        outputs["auc_inf"]={"value":dose/cl*1000,"unit":"ng*h/mL","prediction_type":MECHANISTIC_ESTIMATE}
    else:
        f=float(inputs["f_fraction"]); ka=float(inputs["ka_per_h"])
        outputs["auc_inf"]={"value":f*dose/cl*1000,"unit":"ng*h/mL","prediction_type":MECHANISTIC_ESTIMATE}
        if abs(ka-ke)<1e-10:
            tmax=1/ke; cmax=(f*dose/v)*ke*tmax*math.exp(-ke*tmax)*1000
        else:
            tmax=math.log(ka/ke)/(ka-ke); cmax=(f*dose*ka/(v*(ka-ke)))*(math.exp(-ke*tmax)-math.exp(-ka*tmax))*1000
        outputs["tmax"]={"value":tmax,"unit":"h","prediction_type":MECHANISTIC_ESTIMATE}
        outputs["cmax"]={"value":cmax,"unit":"ng/mL","prediction_type":MECHANISTIC_ESTIMATE}
    return base|{"status":"COMPLETE","prediction_type":MECHANISTIC_ESTIMATE,"model_family":"one_compartment_first_order","outputs":outputs,"missing_inputs":[],"assumptions":["Linear PK and one-compartment disposition"]}
