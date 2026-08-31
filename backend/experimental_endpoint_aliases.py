"""Versioned, context-gated public experimental endpoint aliases.

Aliases alone never qualify a record: ``classify_experimental_endpoint``
requires the source assay context needed for each canonical family.
"""
from __future__ import annotations

import re

VERSION = "drugopt-experimental-endpoint-aliases-v1.1"

_RULES = (
    ("SOLUBILITY", ("solubility",), "ADMET", None),
    ("LOGP", ("logp", "partition coefficient"), "ADMET", None),
    ("LOGD", ("logd", "distribution coefficient"), "ADMET", None),
    ("PKA", ("pka", "acid dissociation"), "ADMET", None),
    ("CACO2_PAPP_AB", ("papp", "apparent permeability", "permeability"), "ADMET", "caco2_ab"),
    ("PPB", ("plasma protein binding", "protein binding", "fraction unbound", "fu"), "ADMET", "ppb"),
    ("MICROSOMAL_STABILITY", ("microsomal stability", "microsome", "intrinsic clearance", "clint"), "METABOLISM", "microsome"),
    ("HEPATOCYTE_STABILITY", ("hepatocyte", "hepatocytes"), "METABOLISM", "hepatocyte"),
    ("CYP", ("cyp", "cytochrome p450"), "METABOLISM", "cyp"),
    ("TRANSPORTER", ("p-gp", "pgp", "bcrp", "bsep", "oatp", "oct", "mate"), "ADMET", "transporter"),
    ("HERG", ("herg", "kcn h2", "kcn h2"), "TOXICITY", "herg"),
    ("AMES", ("ames", "salmonella"), "TOXICITY", "ames"),
    ("DILI", ("drug induced liver injury", "dili"), "TOXICITY", "dili"),
    ("PK", ("auc", "cmax", "tmax", "half-life", "half life", "clearance", "volume of distribution", "bioavailability"), "PK", "pk"),
)


def classify_experimental_endpoint(*, label: str, assay_type: str = "", description: str = "", species: str = "", cell_line: str = "", unit: str = "") -> dict:
    """Return a category plus context protection; never infer missing context."""
    text = " ".join(str(v or "") for v in (label, assay_type, description, species, cell_line, unit)).lower()
    normalized = re.sub(r"[^a-z0-9%]+", " ", text)
    for endpoint, aliases, category, context in _RULES:
        if not any(alias in text or alias.replace("-", " ") in normalized for alias in aliases):
            continue
        if context == "caco2_ab":
            if "caco" not in text:
                return {"endpoint": endpoint, "category": category, "qualified": False, "reason": "Caco-2 cell context required"}
            if any(v in text for v in ("b to a", "b→a", "b-a", "efflux")):
                return {"endpoint": "CACO2_PAPP_BA_OR_EFFLUX", "category": category, "qualified": False, "reason": "Not Caco-2 A→B"}
        if context == "microsome" and "microsom" not in text and "clint" not in text:
            return {"endpoint": endpoint, "category": category, "qualified": False, "reason": "Microsomal context required"}
        if context == "ppb" and "plasma" not in text and "protein binding" not in text and "fraction unbound" not in text and " fu" not in text:
            return {"endpoint": endpoint, "category": category, "qualified": False, "reason": "Plasma/protein-binding context required"}
        return {"endpoint": endpoint, "category": category, "qualified": True, "reason": "Context-qualified alias"}
    return {"endpoint": label or "Unclassified", "category": "ACTIVITY", "qualified": False, "reason": "No curated ADMET/PK alias match"}
