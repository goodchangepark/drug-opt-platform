"""Record v6.3 independent-evidence qualification decisions.

The acquisition gate is intentionally fail-closed.  A related public dataset
is not an accepted validation record unless its endpoint, species, context,
experimental status, and original-study lineage are reproducible.
"""
from __future__ import annotations

import hashlib
import json
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TODAY = str(date.today())


SOURCES = [
    {
        "source_family": "PKSmart",
        "database": "PKSmart Human_PK_data.csv / External_test_315.csv",
        "status": "SAME_SOURCE_FAMILY_NOT_INDEPENDENT",
        "accepted": 0,
        "reason": "Already supplies the development/consumed stress lineage; mirrored or split files cannot qualify external validation.",
        "lineage": "PKSmart",
    },
    {
        "source_family": "PK-DB",
        "database": "PK-DB REST API",
        "status": "NO_REPRODUCIBLE_TARGET_RECORDS",
        "accepted": 0,
        "reason": "Public API discovery was inspected; no reproducible qualified human IV total-CL or strict Vdss/fu validation records were extracted in this gate.",
        "lineage": "PK-DB",
    },
    {
        "source_family": "TDC",
        "database": "PPBR / Clearance_Hepatocyte benchmark families",
        "status": "RELATED_ENDPOINT_NOT_VALIDATION_TARGET",
        "accepted": 0,
        "reason": "Related PPB or hepatocyte datasets do not establish human plasma fu study lineage or human IV total CL/Vdss semantics required here.",
        "lineage": "TDC-associated source datasets",
    },
    {
        "source_family": "ADMET-MTFT",
        "database": "public ADME benchmark files",
        "status": "FOUND_NOT_LINEAGE_QUALIFIED",
        "accepted": 0,
        "reason": "Human fup is documented, but original-study deduplication and independent validation lineage are not reproducible from the current accessible artifact.",
        "lineage": "ADME benchmark compilation",
    },
    {
        "source_family": "FDA-DailyMed-DrugCentral",
        "database": "regulatory labels and curated reference evidence",
        "status": "EXISTING_REFERENCE_FAMILY",
        "accepted": 0,
        "reason": "Current structured PPB/Vdss evidence is already represented by the DrugBank/reference family; no new non-overlapping study set was established.",
        "lineage": "DrugBank/reference",
    },
    {
        "source_family": "ChEMBL-PubChem",
        "database": "bioassay repositories",
        "status": "SEMANTICS_NOT_TARGET",
        "accepted": 0,
        "reason": "Useful assay evidence is not a substitute for clinical human IV total CL, true IV Vdss, or study-qualified plasma fu validation.",
        "lineage": "mixed/public assay aggregation",
    },
]


def endpoint_artifact(endpoint, definition, desired_n, notes):
    return {
        "dataset_version": f"{endpoint.lower()}_independent_v1",
        "created_at": TODAY,
        "endpoint": endpoint,
        "endpoint_definition": definition,
        "canonical_unit": {"HUMAN_VDSS": "L/kg", "HUMAN_TOTAL_IV_CL": "mL/min/kg", "HUMAN_FU_PLASMA": "fraction unbound"}[endpoint],
        "qualification_policy": [
            "human species",
            "exact endpoint semantics",
            "record-level experimental value and unit",
            "chemical identity and chemical form resolution",
            "original study/publication lineage",
            "source-family independence from PKSmart",
        ],
        "source_families_searched": [x["source_family"] for x in SOURCES],
        "qualified_records": [],
        "qualified_n": 0,
        "unique_compounds": 0,
        "unique_studies": 0,
        "unique_scaffolds": 0,
        "rejections": [],
        "desired_new_n": desired_n,
        "role": "LOCKED_EXTERNAL_VALIDATION",
        "status": "NO_QUALIFIED_INDEPENDENT_RECORDS",
        "notes": notes,
    }


def main():
    prior = json.loads((ROOT / "validation/model_qualification_data_closure_v6_2.json").read_text())
    payloads = {
        "HUMAN_VDSS": endpoint_artifact(
            "HUMAN_VDSS", "human IV-derived true steady-state volume; exclude Vz/Vc/Vd/F", 100,
            "PKSmart values without record-level analyte/formulation/study context remain ambiguous; existing 50 structured records are the only strict set and are not an independent source family.",
        ),
        "HUMAN_TOTAL_IV_CL": endpoint_artifact(
            "HUMAN_TOTAL_IV_CL", "human intravenous total systemic clearance; exclude CL/F and organ-specific clearance", 30,
            "No non-PKSmart source-family cohort with reproducible IV total-CL study lineage was qualified.",
        ),
        "HUMAN_FU_PLASMA": endpoint_artifact(
            "HUMAN_FU_PLASMA", "human plasma fraction unbound or explicitly invertible fraction bound", 50,
            "The current candidate uses DrugBank/reference PPB records; PKSmart fup is not independent PPB validation and related public compilations lack reproducible study lineage in this gate.",
        ),
    }
    for endpoint, artifact in payloads.items():
        artifact["dataset_hash"] = hashlib.sha256(json.dumps(artifact, sort_keys=True).encode()).hexdigest()
        (ROOT / f"validation/{endpoint.lower()}_independent_v1.json").write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

    result = {
        "artifact": "independent_evidence_acquisition_v6_3",
        "created_at": TODAY,
        "prior_qualification_artifact": "model_qualification_data_closure_v6_2.json",
        "reference_library_n": prior["reference_library_n"],
        "source_policy": "source_family and original_study are separate; mirrors are not independent",
        "sources_searched": SOURCES,
        "endpoint_datasets": {k: {x: v for x, v in a.items() if x != "qualified_records"} for k, a in payloads.items()},
        "study_lineage": {
            "newly_accepted_records": 0,
            "independent_studies": 0,
            "duplicate_studies_detected_in_accepted_set": 0,
            "cross_database_duplicates": 0,
            "unknown_lineage_records_accepted": 0,
            "policy": "unresolved lineage is held/rejected, never counted as independent",
        },
        "frozen_candidate_evaluation": {
            "PPB_FU": {"external_n": 0, "status": "RESEARCH_CANDIDATE"},
            "HUMAN_TOTAL_IV_CL": {"external_n": 0, "status": "RESEARCH_CANDIDATE"},
            "HUMAN_VDSS": {"external_n": 0, "status": "DATA_LIMITED"},
        },
        "commercial_vendor_framework": {
            "scoring": {
                "endpoint_match": 20,
                "species_context_match": 15,
                "record_level_metadata": 15,
                "original_study_lineage": 15,
                "independent_source_families": 10,
                "unique_compounds_and_scaffolds": 10,
                "duplicate_rate": 5,
                "license_internal_ml_rights": 10,
            },
            "decision_bands": {"HIGH_VALUE": ">=80", "MODERATE_VALUE": "60-79", "LOW_VALUE": "40-59", "REDUNDANT": "<40 or source-family overlap"},
            "minimums": {
                "HUMAN_VDSS": {"minimum_n": 100, "preferred_n": "150-300", "metadata": ["true Vss/Vdss", "human IV", "L/kg", "analyte", "body weight", "study/publication", "chemical form"], "role": "development plus locked validation"},
                "HUMAN_TOTAL_IV_CL": {"minimum_n": 30, "preferred_n": "50-100", "metadata": ["human IV", "total systemic CL", "mL/min/kg", "analyte", "dose", "formulation", "study/publication", "chemical form"], "role": "locked validation"},
                "HUMAN_FU_PLASMA": {"minimum_n": 50, "preferred_n": "100", "metadata": ["human plasma", "fu or bound fraction", "method", "concentration", "temperature", "analyte", "study/publication"], "role": "locked validation"},
            },
        },
        "global_status": "MORE_FOUNDATION_REQUIRED",
        "global_blocker": "Independent, record-level study-qualified validation data for PPB/fu and total IV CL, plus strict Vdss source diversity, are still absent.",
        "production": "UNCHANGED",
    }
    out = ROOT / "validation/independent_evidence_acquisition_v6_3.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"sources": len(SOURCES), "accepted": {k: 0 for k in payloads}, "global_status": result["global_status"]}, indent=2))


if __name__ == "__main__":
    main()
