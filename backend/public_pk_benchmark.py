"""Immutable public-PK benchmark curation and baseline helpers (v4.4B).

This module deliberately has no database dependency.  Public benchmark
observations are validation infrastructure, never project evidence.  The
curated v1 set is intentionally limited: it freezes only source-backed,
context-complete regulatory observations and makes no predictive-performance
claim.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from copy import deepcopy
import hashlib
import json
import math
from typing import Any, Iterable

from rdkit import Chem
from rdkit.Chem import inchi

from backend.pk_engine_v1 import PK_ENGINE_VERSION, estimate_one_compartment


BENCHMARK_VERSION = "drugopt-public-pk-benchmark-v1"
BENCHMARK_SCHEMA_VERSION = "drugopt-public-pk-benchmark-schema-v1"
SPLIT_SEED = "drugopt-public-pk-benchmark-v1-4404"
TRACK_A = "PREDICTIVE_VALIDATION"
TRACK_B = "MECHANISTIC_VERIFICATION"
DEVELOPMENT = "DEVELOPMENT"
FINAL_EVALUATION = "FINAL_EVALUATION"


# A source link is provenance for an observation, not another independent
# observation.  Values are transcribed only where the public source gives the
# parameter and context explicitly.  Ranges and context-incomplete claims are
# excluded to the review queue instead of being converted into point values.
SOURCES: dict[str, dict[str, str]] = {
    "FDA_SUNVOZERTINIB_NDA_219839": {
        "source_type": "PRIMARY_REGULATORY",
        "quality_tier": "TIER_1",
        "document": "FDA NDA 219839 multidisciplinary review (ZEGFROVY/sunvozertinib)",
        "url": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2025/219839Orig1s000MultidisciplineR.pdf",
        "section": "Clinical Pharmacology / Drug Exposure at Steady State",
    },
    "DAILYMED_ACETAMINOPHEN_IV": {
        "source_type": "PRIMARY_REGULATORY",
        "quality_tier": "TIER_1",
        "document": "DailyMed acetaminophen injection label",
        "url": "https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=f80b6f72-f280-4b4a-a6a9-e6490b94705b",
        "section": "12.3 Pharmacokinetics, Table 5",
    },
    "DAILYMED_METFORMIN_IR": {
        "source_type": "PRIMARY_REGULATORY",
        "quality_tier": "TIER_1",
        "document": "DailyMed metformin hydrochloride immediate-release label",
        "url": "https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=542cec22-eeae-4704-9bb6-4176640e5ea8",
        "section": "12.3 Pharmacokinetics, Table 3",
    },
    "DAILYMED_OSIMERTINIB": {
        "source_type": "PRIMARY_REGULATORY",
        "quality_tier": "TIER_1",
        "document": "DailyMed TAGRISSO (osimertinib) label",
        "url": "https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=5e81b4a7-b971-45e1-9c31-29cea8c87ce7",
        "section": "12.3 Pharmacokinetics",
    },
    "DAILYMED_MIDAZOLAM": {
        "source_type": "PRIMARY_REGULATORY",
        "quality_tier": "TIER_1",
        "document": "DailyMed midazolam label",
        "url": "https://dailymed.nlm.nih.gov/dailymed/fda/fdaDrugXsl.cfm?setid=90f7b7f6-d0fd-40af-bbfa-70b152e3e27c",
        "section": "Clinical Pharmacology, Absorption",
    },
    "DAILYMED_WARFARIN": {
        "source_type": "PRIMARY_REGULATORY",
        "quality_tier": "TIER_1",
        "document": "DailyMed COUMADIN (warfarin sodium) label",
        "url": "https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=91fa852c-b43d-4a55-983b-74aa6827125d",
        "section": "12.3 Pharmacokinetics, Excretion",
    },
}


COMPOUNDS: dict[str, dict[str, str]] = {
    "acetaminophen": {"canonical_name": "Acetaminophen", "aliases": ["paracetamol"], "smiles": "CC(=O)NC1=CC=C(O)C=C1", "pubchem_cid": "1983", "identity_source": "PubChem CID 1983"},
    "metformin": {"canonical_name": "Metformin", "aliases": ["metformin hydrochloride"], "smiles": "CN(C)C(=N)NC(N)=N", "pubchem_cid": "4091", "identity_source": "PubChem CID 4091"},
    "midazolam": {"canonical_name": "Midazolam", "aliases": [], "smiles": "Cc1ncc2c(n1)N=C(c1ccccc1F)CN2c1ccc(Cl)cc1", "pubchem_cid": "4192", "identity_source": "PubChem CID 4192"},
    "osimertinib": {"canonical_name": "Osimertinib", "aliases": ["TAGRISSO"], "smiles": "C=CC(=O)Nc1cc(Nc2nccc(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C", "pubchem_cid": "71496458", "identity_source": "PubChem CID 71496458"},
    "sunvozertinib": {"canonical_name": "Sunvozertinib", "aliases": ["ZEGFROVY"], "smiles": "COC1=CC(N2CC[C@H](C2)N(C)C)=C(NC(=O)C=C)C=C1NC1=NC=CC(NC2=CC(Cl)=C(F)C=C2C(C)(C)O)=N1", "pubchem_cid": "139377809", "identity_source": "PubChem CID 139377809"},
    "warfarin": {"canonical_name": "Warfarin", "aliases": ["COUMADIN"], "smiles": "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O", "pubchem_cid": "54678486", "identity_source": "PubChem CID 54678486"},
}


def _observation(
    oid: str, compound: str, *, source: str, study: str, parameter: str,
    value: float, unit: str, route: str, dose: float | None, dose_unit: str | None,
    regimen: str, steady_state: bool, population: str, subtype: str = "",
    analyte: str = "PARENT", normalized_value: float | None = None,
    normalized_unit: str | None = None, notes: str = "",
) -> dict[str, Any]:
    source_row = SOURCES[source]
    return {
        "benchmark_observation_id": oid,
        "benchmark_version": BENCHMARK_VERSION,
        "track": TRACK_A,
        "compound_id": compound,
        "canonical_compound": COMPOUNDS[compound]["canonical_name"],
        "species": "HUMAN",
        "strain_or_population": population,
        "route": route,
        "dose": dose,
        "dose_unit": dose_unit,
        "regimen": regimen,
        "single_or_multiple_dose": "MULTIPLE" if steady_state else "SINGLE",
        "steady_state": steady_state,
        "fed_fasted": "UNSPECIFIED",
        "analyte": analyte,
        "parent_or_metabolite": "PARENT" if analyte == "PARENT" else "METABOLITE",
        "matrix": "PLASMA",
        "canonical_parameter": parameter,
        "parameter_subtype": subtype or parameter,
        "raw_value": value,
        "raw_unit": unit,
        "normalized_value": normalized_value if normalized_value is not None else value,
        "normalized_unit": normalized_unit or unit,
        "relation": "=",
        "study_id": study,
        "independence_id": f"{compound}|{study}",
        "source_ids": [source],
        "source_type": source_row["source_type"],
        "quality_tier": source_row["quality_tier"],
        "source_document": source_row["document"],
        "source_url": source_row["url"],
        "section_table": source_row["section"],
        "qualification_status": "BENCHMARK_QUALIFIED",
        "identity_status": "EXACT_STRUCTURE_MATCH",
        "notes": notes,
    }


def curated_observations() -> list[dict[str, Any]]:
    """Return a new copy of the frozen, source-backed v1 observation set."""
    rows = [
        _observation("PKB1-APAP-AUC", "acetaminophen", source="DAILYMED_ACETAMINOPHEN_IV", study="APAP-IV-ADULT-TABLE5", parameter="AUC0_T", value=43, unit="µg*h/mL", normalized_value=43000, normalized_unit="ng*h/mL", route="IV", dose=1000, dose_unit="mg", regimen="single 15-minute infusion", steady_state=False, population="healthy adults", subtype="AUC0-6h"),
        _observation("PKB1-APAP-CMAX", "acetaminophen", source="DAILYMED_ACETAMINOPHEN_IV", study="APAP-IV-ADULT-TABLE5", parameter="CMAX", value=28, unit="µg/mL", normalized_value=28000, normalized_unit="ng/mL", route="IV", dose=1000, dose_unit="mg", regimen="single 15-minute infusion", steady_state=False, population="healthy adults"),
        _observation("PKB1-APAP-THALF", "acetaminophen", source="DAILYMED_ACETAMINOPHEN_IV", study="APAP-IV-ADULT-TABLE5", parameter="T_HALF", value=2.4, unit="h", route="IV", dose=1000, dose_unit="mg", regimen="single 15-minute infusion", steady_state=False, population="healthy adults", subtype="terminal"),
        _observation("PKB1-APAP-CL", "acetaminophen", source="DAILYMED_ACETAMINOPHEN_IV", study="APAP-IV-ADULT-TABLE5", parameter="CL", value=0.27, unit="L/h/kg", normalized_value=4.5, normalized_unit="mL/min/kg", route="IV", dose=1000, dose_unit="mg", regimen="single 15-minute infusion", steady_state=False, population="healthy adults"),
        _observation("PKB1-APAP-VSS", "acetaminophen", source="DAILYMED_ACETAMINOPHEN_IV", study="APAP-IV-ADULT-TABLE5", parameter="VSS", value=0.8, unit="L/kg", route="IV", dose=1000, dose_unit="mg", regimen="single 15-minute infusion", steady_state=False, population="healthy adults"),
        _observation("PKB1-MET-CMAX-500", "metformin", source="DAILYMED_METFORMIN_IR", study="MET-IR-HEALTHY-500-SD", parameter="CMAX", value=1.03, unit="µg/mL", normalized_value=1030, normalized_unit="ng/mL", route="ORAL", dose=500, dose_unit="mg", regimen="single dose", steady_state=False, population="healthy nondiabetic adults"),
        _observation("PKB1-MET-TMAX-500", "metformin", source="DAILYMED_METFORMIN_IR", study="MET-IR-HEALTHY-500-SD", parameter="TMAX", value=2.75, unit="h", route="ORAL", dose=500, dose_unit="mg", regimen="single dose", steady_state=False, population="healthy nondiabetic adults"),
        _observation("PKB1-MET-CMAX-850", "metformin", source="DAILYMED_METFORMIN_IR", study="MET-IR-HEALTHY-850-SD", parameter="CMAX", value=1.60, unit="µg/mL", normalized_value=1600, normalized_unit="ng/mL", route="ORAL", dose=850, dose_unit="mg", regimen="single dose", steady_state=False, population="healthy nondiabetic adults"),
        _observation("PKB1-MET-TMAX-850", "metformin", source="DAILYMED_METFORMIN_IR", study="MET-IR-HEALTHY-850-SD", parameter="TMAX", value=2.64, unit="h", route="ORAL", dose=850, dose_unit="mg", regimen="single dose", steady_state=False, population="healthy nondiabetic adults"),
        _observation("PKB1-MET-VDF", "metformin", source="DAILYMED_METFORMIN_IR", study="MET-IR-VDF-850-SD", parameter="VDF", value=654, unit="L", route="ORAL", dose=850, dose_unit="mg", regimen="single dose", steady_state=False, population="healthy adults", subtype="apparent oral volume"),
        _observation("PKB1-MET-THALF", "metformin", source="DAILYMED_METFORMIN_IR", study="MET-IR-PLASMA-ELIMINATION", parameter="T_HALF", value=6.2, unit="h", route="ORAL", dose=None, dose_unit=None, regimen="clinical oral dosing", steady_state=False, population="adults", subtype="plasma elimination", notes="Dose is not required for half-life eligibility; source retains clinical context."),
        _observation("PKB1-OSIM-TMAX", "osimertinib", source="DAILYMED_OSIMERTINIB", study="OSIM-PO-POP-PK", parameter="TMAX", value=6, unit="h", route="ORAL", dose=80, dose_unit="mg", regimen="once daily", steady_state=True, population="patients", subtype="median"),
        _observation("PKB1-OSIM-VSSF", "osimertinib", source="DAILYMED_OSIMERTINIB", study="OSIM-PO-POP-PK", parameter="VSSF", value=918, unit="L", route="ORAL", dose=80, dose_unit="mg", regimen="once daily", steady_state=True, population="patients", subtype="apparent steady-state volume"),
        _observation("PKB1-OSIM-THALF", "osimertinib", source="DAILYMED_OSIMERTINIB", study="OSIM-PO-POP-PK", parameter="T_HALF", value=48, unit="h", route="ORAL", dose=80, dose_unit="mg", regimen="once daily", steady_state=True, population="patients", subtype="population estimated"),
        _observation("PKB1-OSIM-CLF", "osimertinib", source="DAILYMED_OSIMERTINIB", study="OSIM-PO-POP-PK", parameter="CLF", value=14.3, unit="L/h", route="ORAL", dose=80, dose_unit="mg", regimen="once daily", steady_state=True, population="patients", subtype="apparent oral clearance"),
        _observation("PKB1-MDZ-CMAX", "midazolam", source="DAILYMED_MIDAZOLAM", study="MDZ-HEALTHY-7.5-IM", parameter="CMAX", value=90, unit="ng/mL", route="IM", dose=7.5, dose_unit="mg", regimen="single dose", steady_state=False, population="healthy subjects"),
        _observation("PKB1-MDZ-TMAX", "midazolam", source="DAILYMED_MIDAZOLAM", study="MDZ-HEALTHY-7.5-IM", parameter="TMAX", value=0.5, unit="h", route="IM", dose=7.5, dose_unit="mg", regimen="single dose", steady_state=False, population="healthy subjects"),
        _observation("PKB1-WARF-THALF", "warfarin", source="DAILYMED_WARFARIN", study="WARFARIN-EFFECTIVE-HALFLIFE", parameter="T_HALF", value=40, unit="h", route="ORAL", dose=None, dose_unit=None, regimen="single dose", steady_state=False, population="adults", subtype="effective half-life", notes="Dose is not required for half-life eligibility; terminal and effective half-life remain distinct."),
        _observation("PKB1-SUN-AUC-200", "sunvozertinib", source="FDA_SUNVOZERTINIB_NDA_219839", study="SUN-PO-200-QD-CLINICAL", parameter="AUC_TAU", value=8060, unit="ng*h/mL", route="ORAL", dose=200, dose_unit="mg/day", regimen="once daily", steady_state=True, population="NSCLC patients", subtype="steady-state AUC", notes="FDA review table context retained from NDA219839 MultidisciplineR.pdf:L1623."),
        _observation("PKB1-SUN-CMAX-200", "sunvozertinib", source="FDA_SUNVOZERTINIB_NDA_219839", study="SUN-PO-200-QD-CLINICAL", parameter="CMAX", value=412, unit="ng/mL", route="ORAL", dose=200, dose_unit="mg", regimen="once daily", steady_state=True, population="NSCLC patients", subtype="steady-state Cmax", notes="FDA review table context retained from NDA219839 MultidisciplineR.pdf:L1623."),
        _observation("PKB1-SUN-TMAX-200", "sunvozertinib", source="FDA_SUNVOZERTINIB_NDA_219839", study="SUN-PO-200-QD-CLINICAL", parameter="TMAX", value=7, unit="h", route="ORAL", dose=200, dose_unit="mg", regimen="once daily", steady_state=True, population="NSCLC patients", subtype="median Tmax", notes="FDA review table context retained from NDA219839 MultidisciplineR.pdf:L1623."),
        _observation("PKB1-SUN-AUC", "sunvozertinib", source="FDA_SUNVOZERTINIB_NDA_219839", study="SUN-PO-300-QD-STEADY-STATE", parameter="AUC_TAU", value=12089, unit="ng*h/mL", route="ORAL", dose=300, dose_unit="mg", regimen="once daily", steady_state=True, population="NSCLC patients", subtype="steady-state AUC"),
        _observation("PKB1-SUN-CMAX", "sunvozertinib", source="FDA_SUNVOZERTINIB_NDA_219839", study="SUN-PO-300-QD-STEADY-STATE", parameter="CMAX", value=619, unit="ng/mL", route="ORAL", dose=300, dose_unit="mg", regimen="once daily", steady_state=True, population="NSCLC patients", subtype="steady-state geometric mean"),
        _observation("PKB1-SUN-TMAX", "sunvozertinib", source="FDA_SUNVOZERTINIB_NDA_219839", study="SUN-PO-POP-PK", parameter="TMAX", value=6, unit="h", route="ORAL", dose=300, dose_unit="mg", regimen="once daily", steady_state=True, population="NSCLC patients", subtype="median"),
        _observation("PKB1-SUN-VDF", "sunvozertinib", source="FDA_SUNVOZERTINIB_NDA_219839", study="SUN-PO-POP-PK", parameter="VDF", value=2116, unit="L", route="ORAL", dose=300, dose_unit="mg", regimen="once daily", steady_state=True, population="NSCLC patients", subtype="apparent oral volume"),
        _observation("PKB1-SUN-THALF", "sunvozertinib", source="FDA_SUNVOZERTINIB_NDA_219839", study="SUN-PO-POP-PK", parameter="T_HALF", value=50, unit="h", route="ORAL", dose=300, dose_unit="mg", regimen="once daily", steady_state=True, population="NSCLC patients", subtype="elimination half-life"),
        _observation("PKB1-SUN-CLF", "sunvozertinib", source="FDA_SUNVOZERTINIB_NDA_219839", study="SUN-PO-POP-PK", parameter="CLF", value=29, unit="L/h", route="ORAL", dose=300, dose_unit="mg", regimen="once daily", steady_state=True, population="NSCLC patients", subtype="apparent oral clearance"),
    ]
    return deepcopy(rows)


REVIEW_QUEUE = [
    {"candidate_id": "PKB-RQ-SUN-RAT-F", "compound": "Sunvozertinib", "parameter": "F", "species": "RAT", "reason": "DOSE_MISSING", "raw_context": "Oral bioavailability estimated at 39.6% in rats", "source_id": "FDA_SUNVOZERTINIB_NDA_219839", "action": "Keep out of quantitative benchmark until the source study dose/regimen is independently resolved."},
    {"candidate_id": "PKB-RQ-SUN-DOG-F", "compound": "Sunvozertinib", "parameter": "F", "species": "DOG", "reason": "DOSE_MISSING", "raw_context": "Oral bioavailability estimated at 48.8% in dogs", "source_id": "FDA_SUNVOZERTINIB_NDA_219839", "action": "Keep out of quantitative benchmark until the source study dose/regimen is independently resolved."},
    {"candidate_id": "PKB-RQ-MDZ-VD", "compound": "Midazolam", "parameter": "VD", "species": "HUMAN", "reason": "MEASUREMENT_TYPE_MISSING", "raw_context": "Vd ranged from 1.0 to 3.1 L/kg across six studies", "source_id": "DAILYMED_MIDAZOLAM", "action": "Do not turn a heterogeneous range into a single benchmark target."},
]


def _identity(compound_id: str) -> dict[str, Any]:
    entry = COMPOUNDS[compound_id]
    molecule = Chem.MolFromSmiles(entry["smiles"])
    if molecule is None:
        raise ValueError(f"invalid curated structure for {compound_id}")
    return entry | {"inchikey": inchi.MolToInchiKey(molecule), "canonical_smiles": Chem.MolToSmiles(molecule, canonical=True), "identity_status": "EXACT_STRUCTURE_MATCH"}


def compounds_with_identity() -> list[dict[str, Any]]:
    return [{"compound_id": key, **_identity(key)} for key in sorted(COMPOUNDS)]


def observation_fingerprint(row: dict[str, Any]) -> str:
    fields = ("compound_id", "study_id", "species", "route", "dose", "dose_unit", "regimen", "analyte", "canonical_parameter", "parameter_subtype", "normalized_value", "normalized_unit")
    payload = {field: row.get(field) for field in fields}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def deduplicate_observations(rows: Iterable[dict[str, Any]]) -> tuple[list[dict[str, Any]], int]:
    """Collapse source representations of one scientific observation."""
    output: dict[str, dict[str, Any]] = {}
    duplicates = 0
    for raw in rows:
        row = deepcopy(raw); fingerprint = observation_fingerprint(row)
        row["scientific_fingerprint"] = fingerprint
        previous = output.get(fingerprint)
        if previous is None:
            output[fingerprint] = row
        else:
            previous["source_ids"] = sorted(set(previous["source_ids"]) | set(row["source_ids"]))
            duplicates += 1
    return sorted(output.values(), key=lambda row: row["benchmark_observation_id"]), duplicates


def validate_observation(row: dict[str, Any]) -> list[str]:
    required = ("compound_id", "species", "route", "analyte", "canonical_parameter", "parameter_subtype", "raw_value", "raw_unit", "normalized_value", "normalized_unit", "study_id", "source_ids", "quality_tier", "qualification_status", "independence_id")
    missing = [field for field in required if row.get(field) in (None, "", [])]
    if row.get("canonical_parameter") in {"CMAX", "AUC0_T", "AUC0_INF", "AUC_TAU", "F"} and row.get("dose") is None:
        missing.append("dose")
    if row.get("canonical_parameter") in {"CMAX", "AUC0_T", "AUC0_INF", "AUC_TAU"} and not row.get("regimen"):
        missing.append("regimen")
    if row.get("compound_id") not in COMPOUNDS:
        missing.append("exact_identity")
    return sorted(set(missing))


def dataset_hash(rows: Iterable[dict[str, Any]]) -> str:
    material = [{key: row[key] for key in sorted(row) if key not in {"scientific_fingerprint"}} for row in rows]
    return hashlib.sha256(json.dumps(material, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _split_rank(compound_id: str) -> str:
    return hashlib.sha256(f"{SPLIT_SEED}|{compound_id}".encode()).hexdigest()


def freeze_compound_split(rows: Iterable[dict[str, Any]]) -> dict[str, Any]:
    compounds = sorted({row["compound_id"] for row in rows}, key=_split_rank)
    holdout_count = max(1, round(len(compounds) * 0.25)) if compounds else 0
    holdout = set(compounds[:holdout_count])
    assignment = {compound: ("HOLDOUT" if compound in holdout else DEVELOPMENT) for compound in compounds}
    scaffold_summary: dict[str, list[str]] = {}
    for compound in compounds:
        molecule = Chem.MolFromSmiles(COMPOUNDS[compound]["smiles"])
        # The canonical SMILES is a deterministic and inspectable diversity
        # signature here; no model fitting consumes it.
        scaffold_summary[compound] = [Chem.MolToSmiles(molecule, canonical=True)]
    return {"benchmark_version": BENCHMARK_VERSION, "split_version": "drugopt-public-pk-split-v1", "seed": SPLIT_SEED, "assignment": assignment, "development_compounds": sorted(c for c in compounds if assignment[c] == DEVELOPMENT), "holdout_compounds": sorted(holdout), "compound_overlap": [], "holdout_leakage": False, "scaffold_summary": scaffold_summary, "split_hash": hashlib.sha256(json.dumps(assignment, sort_keys=True).encode()).hexdigest(), "frozen": True}


def records_for_mode(rows: Iterable[dict[str, Any]], split: dict[str, Any], *, mode: str, partition: str) -> list[dict[str, Any]]:
    """Access guard: development callers cannot retrieve held-out targets."""
    if mode == DEVELOPMENT and partition == "HOLDOUT":
        raise PermissionError("Final-holdout experimental targets are locked during development mode")
    if mode not in {DEVELOPMENT, FINAL_EVALUATION}:
        raise ValueError(f"unknown benchmark mode {mode}")
    if partition not in {DEVELOPMENT, "HOLDOUT"}:
        raise ValueError(f"unknown benchmark partition {partition}")
    return [deepcopy(row) for row in rows if split["assignment"][row["compound_id"]] == partition]


def baseline_rows(rows: Iterable[dict[str, Any]], *, split: dict[str, Any], partition: str, mode: str) -> list[dict[str, Any]]:
    """Run the unchanged fail-closed PK overlay without experimental targets."""
    selected = records_for_mode(rows, split, mode=mode, partition=partition)
    output = []
    for observed in selected:
        # Track A contains no endpoint target in the request; current v1 has no
        # structure-only input assembly, so it correctly fails closed.
        result = estimate_one_compartment(species=observed["species"], route=observed["route"], inputs={})
        output.append({"benchmark_observation_id": observed["benchmark_observation_id"], "track": TRACK_A, "partition": partition, "engine": PK_ENGINE_VERSION, "observed_value": observed["normalized_value"], "observed_unit": observed["normalized_unit"], "prediction_available": False, "prediction_type": result["prediction_type"], "input_completeness": "INSUFFICIENT", "status": result["status"], "missing_inputs": result["missing_inputs"], "predicted_value": None, "predicted_unit": None, "error": None, "error_metric_type": None})
    return output


def mechanistic_verification() -> list[dict[str, Any]]:
    """Track B verifies an equation only; it is excluded from predictive metrics."""
    derived = math.log(2) * 0.8 / 0.27
    return [{"track": TRACK_B, "study_id": "APAP-IV-ADULT-TABLE5", "compound_id": "acetaminophen", "species": "HUMAN", "parameter": "T_HALF", "experimental_inputs": {"CL": {"value": 0.27, "unit": "L/h/kg"}, "VSS": {"value": 0.8, "unit": "L/kg"}}, "derived_value": derived, "derived_unit": "h", "observed_value": 2.4, "observed_unit": "h", "difference_h": derived - 2.4, "classification": "MECHANISTIC_VERIFICATION", "counts_as_predictive_validation": False}]


def coverage(rows: Iterable[dict[str, Any]]) -> dict[str, Any]:
    counts: dict[str, dict[str, Any]] = defaultdict(lambda: {"eligible_observations": 0, "independent_compounds": set(), "prediction_attempted": 0, "prediction_available": 0, "insufficient_input": 0, "unavailable": 0})
    for row in rows:
        key = f"{row['species']} {row['canonical_parameter']}"
        item = counts[key]; item["eligible_observations"] += 1; item["independent_compounds"].add(row["compound_id"]); item["prediction_attempted"] += 1; item["insufficient_input"] += 1
    result = {}
    for key, item in sorted(counts.items()):
        compounds = len(item.pop("independent_compounds"))
        item["independent_compounds"] = compounds
        item["coverage_percent"] = 0.0
        item["dataset_status"] = "STRONG_BENCHMARK_CANDIDATE" if compounds >= 30 else ("BENCHMARK_READY_LIMITED" if compounds >= 20 else ("LIMITED_VALIDATION" if compounds >= 10 else "INSUFFICIENT_DATA"))
        result[key] = item
    return result


def benchmark_package() -> dict[str, Any]:
    observations, duplicate_count = deduplicate_observations(curated_observations())
    invalid = {row["benchmark_observation_id"]: validate_observation(row) for row in observations}
    invalid = {key: value for key, value in invalid.items() if value}
    if invalid:
        raise ValueError(f"invalid benchmark observations: {invalid}")
    return {"benchmark_version": BENCHMARK_VERSION, "schema_version": BENCHMARK_SCHEMA_VERSION, "frozen_at": "2026-09-02", "dataset_status": "LIMITED_BENCHMARK", "freeze_status": "SUCCEEDED", "purpose": "Independent public-PK benchmark curation; not project evidence and not a production-validation claim.", "compounds": compounds_with_identity(), "observations": observations, "observation_count": len(observations), "independent_compound_count": len({row['compound_id'] for row in observations}), "independent_study_count": len({row['independence_id'] for row in observations}), "duplicate_source_representations_collapsed": duplicate_count, "dataset_hash": dataset_hash(observations), "tracks": {TRACK_A: "Experimental target is not supplied to the PK engine.", TRACK_B: "Experimental inputs may verify equations but never contribute predictive performance."}, "project_safety": {"stored_in_project_database": False, "can_increase_effective_n": False, "can_increase_maturity": False, "can_train_adapter": False}}
