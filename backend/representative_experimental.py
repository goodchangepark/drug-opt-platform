"""Deterministic experimental representative selection for scientific rows.

Policy Version: drugopt-representative-experimental-v5.2

Context-Specific Representative Hierarchy:
1. Approved/Clinical Human Context (e.g. therapeutic clinical dose/regimen, steady state, human plasma)
2. Validated Human In-Vitro (recombinant kinase / primary cell assay / human hepatocytes/microsomes)
3. Relevant Animal In-Vivo / In-Vitro (Rat, Dog, Mouse, NHP)
4. Context Completeness (species, matrix, route, dose, target context)
5. Origin Priority (INTERNAL_EXPERIMENTAL > EXTERNAL_IMPORTED > AUTO_QUALIFIED_EXTERNAL > EXTERNAL_CANDIDATE > RELATED_EXTERNAL)
6. Semantic Comparability (DIRECTLY_COMPARABLE > COMPARABLE_AFTER_DETERMINISTIC_CONVERSION > RELATED_SAME_SCIENTIFIC_GROUP)
7. Source/Reference Quality
8. Stable Identity Determinism

PREDICTION PROXIMITY IS NEVER CONSULTED.
"""
from __future__ import annotations

REPRESENTATIVE_EXPERIMENTAL_VERSION = "drugopt-representative-experimental-v5.2"

_ORIGIN_PRIORITY = {
    "INTERNAL_EXPERIMENTAL": 0,
    "EXTERNAL_IMPORTED": 1,
    "AUTO_QUALIFIED_EXTERNAL": 2,
    "EXTERNAL_CANDIDATE": 3,
    "RELATED_EXTERNAL": 4,
}

_SEMANTIC_PRIORITY = {
    "DIRECTLY_COMPARABLE": 0,
    "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION": 1,
    "QUALIFIED_DIRECT": 0,
    "QUALIFIED_DETERMINISTIC_CONVERSION": 1,
    "RELATED_SAME_SCIENTIFIC_GROUP": 3,
}


def representative_rank(item: dict) -> tuple:
    """Rank without accepting or consulting any predicted numeric value."""
    context = item.get("context") or {}
    qualification = item.get("qualification_details") or {}
    stages = qualification.get("stages") or {}
    origin = item.get("origin") or item.get("state") or "EXTERNAL_CANDIDATE"
    semantic = item.get("comparability") or item.get("qualification") or ""

    # 1. Clinical / Approved Human Context Hierarchy
    species = str(qualification.get("species") or context.get("species") or item.get("species") or "").upper()
    regimen = str(qualification.get("regimen") or context.get("regimen") or item.get("regimen") or "").upper()
    target_context = str(qualification.get("target_context") or context.get("target") or "").upper()

    # Clinical relevance rank
    # Human steady-state clinical dose: 0
    # Human clinical single-dose: 1
    # Human in-vitro / primary assay: 2
    # Animal in-vivo / in-vitro: 3
    # Unspecified: 5
    if species == "HUMAN":
        if "STEADY" in regimen or "QD" in regimen or "DAILY" in regimen or "CLINICAL" in str(context).upper():
            clinical_rank = 0
        elif "SINGLE" in regimen or "ORAL" in str(context).upper() or "IV" in str(context).upper():
            clinical_rank = 1
        else:
            clinical_rank = 2
    elif species in {"RAT", "DOG", "MOUSE", "MONKEY"}:
        clinical_rank = 3
    else:
        clinical_rank = 5

    # 2. Target Specificity / Pharmacologic Assay Relevance Rank
    # Direct Agonist / Primary Functional Response > Target Specific Mutation > General Agonist/Inhibitor > Allosteric Modulation (PAM/NAM) > Kinome / Binding panel
    ctx_str = f"{str(context)} {str(qualification)} {item.get('raw_endpoint', '')}".upper()
    is_pam_or_nam = any(term in ctx_str for term in ("ALLOSTERIC", "PAM", "NAM", "POTENTIATION"))
    has_primary_functional = any(term in ctx_str for term in ("CAMP ASSAY", "AGONIST ACTIVITY", "KINASE ACTIVITY", "AUTOPHOSPHORYLATION"))
    has_target_mutation = any(m in target_context for m in ("EXON20INS", "CAMP", "MUTANT", "T790M", "L858R"))

    if (has_primary_functional or has_target_mutation) and not is_pam_or_nam:
        target_rank = 0
    elif not is_pam_or_nam and target_context not in {"", "GENERAL", "UNSPECIFIED"}:
        target_rank = 1
    elif is_pam_or_nam:
        target_rank = 2
    else:
        target_rank = 3

    # 3. Context Completeness
    complete_context = bool(stages.get("CONTEXT_QUALIFIED")) or (
        species not in {"", "UNSPECIFIED"} and
        context.get("matrix") not in {None, "", "UNSPECIFIED"}
    )

    # 4. Source Quality
    source_quality = str(context.get("source_quality") or item.get("source_quality") or "D").upper()
    reference = item.get("reference") or {}
    reference_present = bool(reference.get("reference") or reference.get("url") or reference.get("source_record_id"))
    stable_id = str(item.get("display_evidence_group_id") or item.get("independent_experiment_group_id") or item.get("id") or "")

    return (
        clinical_rank,
        target_rank,
        _ORIGIN_PRIORITY.get(origin, 9),
        _SEMANTIC_PRIORITY.get(semantic, 5),
        0 if complete_context else 1,
        source_quality,
        0 if reference_present else 1,
        stable_id,
    )


def select_representative(items: list[dict]) -> tuple[dict | None, str]:
    eligible = [item for item in items if item.get("display", {}).get("value") is not None]
    if not eligible:
        return None, "NO_DISPLAYABLE_NUMERIC_OBSERVATION"
    selected = min(eligible, key=representative_rank)
    return selected, f"{REPRESENTATIVE_EXPERIMENTAL_VERSION}: clinical context, target specificity, origin priority, semantic comparability, completeness, source quality"
