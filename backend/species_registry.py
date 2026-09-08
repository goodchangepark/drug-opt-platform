"""Canonical species identities and prediction-support declarations."""
from __future__ import annotations

from typing import Any

SPECIES_REGISTRY = {
    "HUMAN": {
        "display_name": "Human",
        "scientific_name": "Homo sapiens",
        "aliases": ("human", "homo sapiens", "patient", "patients", "healthy volunteer", "healthy volunteers", "clinical", "hlm"),
        "physiology_available": True,
    },
    "RAT": {
        "display_name": "Rat",
        "scientific_name": "Rattus norvegicus",
        "aliases": ("rat", "rats", "rattus norvegicus", "sd rat", "sprague-dawley rat", "sprague dawley rat", "rlm"),
        "physiology_available": True,
    },
    "MOUSE": {
        "display_name": "Mouse",
        "scientific_name": "Mus musculus",
        "aliases": ("mouse", "mice", "mus musculus", "mlm"),
        "physiology_available": True,
    },
    "DOG": {
        "display_name": "Dog",
        "scientific_name": "Canis lupus familiaris",
        "aliases": ("dog", "dogs", "beagle", "canine", "dlm"),
        "physiology_available": True,
    },
    "MONKEY": {
        "display_name": "Monkey",
        "scientific_name": "Non-human primate",
        "aliases": ("monkey", "monkeys", "cynomolgus", "cynomolgus monkey", "nonhuman primate", "non-human primate", "rhesus", "nhp", "cyno"),
        "physiology_available": True,
    },
}


def normalize_species_code(value: Any, context: Any = "", *, allow_other: bool = True) -> str:
    text = f"{value or ''} {context or ''}".lower().replace("–", "-")
    aliases = [
        (alias, code)
        for code, record in SPECIES_REGISTRY.items()
        for alias in record["aliases"]
    ]
    for alias, code in sorted(aliases, key=lambda item: -len(item[0])):
        if alias in text:
            return code
    raw = str(value or "").strip().upper()
    if raw in SPECIES_REGISTRY:
        return raw
    if raw in {"", "UNSPECIFIED", "UNKNOWN", "N/A", "NA", "NONE"}:
        return "UNSPECIFIED"
    if allow_other:
        return "OTHER"
    raise ValueError(f"Unsupported species: {value!r}; choose {', '.join(SPECIES_REGISTRY)}")


def species_registry_payload() -> dict:
    return {
        "registry_version": "DRUGOPT_SPECIES_V1",
        "canonical_species": list(SPECIES_REGISTRY),
        "species": [
            {
                "code": code,
                "display_name": record["display_name"],
                "scientific_name": record["scientific_name"],
                "physiology_available": record["physiology_available"],
            }
            for code, record in SPECIES_REGISTRY.items()
        ],
        "invariants": [
            "Every PK evidence and prediction row declares a species.",
            "Cross-species substitution is never implicit.",
            "Pooled metrics never replace per-species validation metrics.",
        ],
    }
