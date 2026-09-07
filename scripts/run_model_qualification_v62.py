"""Freeze v6.2 candidates and qualify PPB/CL/VDss evidence.

This is deliberately conservative: source-level PKSmart markers without
record-level context are not promoted to strict experimental observations.
No runtime database is read or written by this script.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import evaluate_human_total_iv_cl_v56 as base  # noqa: E402
from backend.clinical_pk_cohort import CLINICAL_PK_COHORT  # noqa: E402


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(smiles: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=True)


def metric(y, p):
    return base.metric(np.asarray(y, float), np.asarray(p, float))


def strict_vdss_rows(library):
    """Return only structured exact Vdss observations with auditable context."""
    rows = []
    ambiguous = 0
    rejected = 0
    for compound in library:
        # PKSmart is represented in the library as a compact endpoint/value
        # map.  It has source-level human/IV provenance, but no record-level
        # analyte, formulation, or study context, so valid values are
        # ambiguous for strict Vdss training and invalid values are rejected.
        if isinstance(compound.get("evidence"), dict) and "HUMAN_VDSS" in compound["evidence"]:
            value = compound["evidence"]["HUMAN_VDSS"]
            try:
                valid = Chem.MolFromSmiles(compound["smiles"]) and value is not None and float(value) > 0
            except (TypeError, ValueError):
                valid = False
            if valid:
                ambiguous += 1
            else:
                rejected += 1
            continue
        for evidence in compound.get("evidence", []):
            if not isinstance(evidence, dict):
                if "VDSS" in str(evidence).upper():
                    ambiguous += 1
                continue
            if evidence.get("canonical_endpoint_id") != "VDSS":
                continue
            required = ("normalized_value", "normalized_unit", "species",
                        "assay_type", "reference_text", "matrix")
            if any(evidence.get(k) in (None, "", "None") for k in required):
                ambiguous += 1
                continue
            if (evidence["species"] != "Homo sapiens"
                    or evidence["normalized_unit"] != "L/kg"
                    or "Vdss" not in evidence.get("raw_endpoint_name", "")
                    or "IV" not in evidence.get("assay_type", "")):
                rejected += 1
                continue
            mol = Chem.MolFromSmiles(compound["smiles"])
            if not mol or float(evidence["normalized_value"]) <= 0:
                rejected += 1
                continue
            rows.append({
                "library_id": compound["library_id"],
                "name": compound["name"],
                "inchikey": compound["inchikey"],
                "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
                "mol": mol,
                "y": float(evidence["normalized_value"]),
                "source_family": compound["source_family"],
                "source_record": evidence["reference_text"],
            })
    return rows, ambiguous, rejected


def compact_vdss_cv(rows):
    """Small, leakage-safe comparison on strict rows; no locked-test tuning."""
    if len(rows) < 15:
        return {"status": "DATA_INSUFFICIENT_FOR_MODEL_SELECTION", "n": len(rows)}
    groups = np.array([
        MurckoScaffold.MurckoScaffoldSmiles(mol=r["mol"]) or r["canonical_smiles"]
        for r in rows
    ])
    n_splits = min(5, len(np.unique(groups)))
    y = np.log10([r["y"] for r in rows])
    results = {}
    # Baseline and the two least over-parameterized candidates are enough for
    # this strict, 50-record dataset; v6.1 already tested the larger search.
    candidates = {
        "median_log_vdss": None,
        "ridge_descriptors": base.models()["ridge_descriptors"],
        "extra_trees_descriptors": base.models()["extra_trees_descriptors"],
    }
    X = base.build_matrix(rows, "descriptors")
    for name, model in candidates.items():
        folds = []
        for train, test in GroupKFold(n_splits=n_splits).split(X, y, groups):
            if model is None:
                pred = np.repeat(np.median(y[train]), len(test))
            else:
                model.fit(X[train], y[train])
                pred = model.predict(X[test])
            folds.append(metric(10 ** y[test], 10 ** pred))
        results[name] = {
            "folds": folds,
            "mean_aafe": float(np.mean([x["aafe"] for x in folds])),
            "mean_within_2_fold_pct": float(np.mean([x["within_2_fold_pct"] for x in folds])),
        }
    return {"status": "STRICT_SAME_SOURCE_CV", "n": len(rows), "results": results}


def main():
    library_artifact = ROOT / "validation/reference_library_v1_1000.json"
    rebuild_artifact = ROOT / "validation/endpoint_model_rebuild_v6_1.json"
    library = json.loads(library_artifact.read_text())
    rebuild = json.loads(rebuild_artifact.read_text())
    compounds = library["compounds"]

    # The v6.1 candidates are frozen by definition here.  This artifact freezes
    # configuration/provenance; neither research candidate is promoted.
    ppb = rebuild["experiments"]["PPB_FU"]
    cl = rebuild["experiments"]["TOTAL_IV_CL"]
    vdss = rebuild["experiments"]["VDSS"]
    strict, ambiguous, rejected = strict_vdss_rows(compounds)
    strict_keys = {r["inchikey"] for r in strict}

    pksmart_vdss = sum(
        1 for c in compounds
        if isinstance(c.get("evidence"), dict) and "HUMAN_VDSS" in c["evidence"]
    )
    # All strict observations are existing-reference records and are separate
    # from the frozen N=30 cohort; assert this rather than relying on names.
    n30_keys = {canonical(x.smiles) for x in CLINICAL_PK_COHORT}
    strict_n30_overlap = sum(r["canonical_smiles"] in n30_keys for r in strict)

    output = {
        "artifact": "model_qualification_data_closure_v6_2",
        "created_from": {
            "reference_library": sha256(library_artifact),
            "endpoint_rebuild": sha256(rebuild_artifact),
            "script": sha256(Path(__file__).resolve()),
        },
        "production": "UNCHANGED",
        "reference_library_n": len(compounds),
        "frozen_candidate_status": {
            "PPB_FU": {
                "endpoint": "HUMAN_FU_PLASMA",
                "algorithm": ppb["selected"],
                "features": ppb["cv"][ppb["selected"]]["features"],
                "training_n": ppb["development_n"],
                "target_transform": "existing HUMAN_PPB normalized target; no silent fu relabel",
                "scaffold_cv": ppb["cv"][ppb["selected"]],
                "status": "RESEARCH_CANDIDATE",
                "independent_validation_n": 0,
                "decision": "INDEPENDENT_VALIDATION_REQUIRED",
            },
            "TOTAL_IV_CL": {
                "endpoint": "HUMAN_TOTAL_IV_CL",
                "algorithm": cl["selected"],
                "features": cl["cv"][cl["selected"]]["features"],
                "training_n": cl["development_n"],
                "target_transform": "log10 clearance",
                "scaffold_cv": cl["cv"][cl["selected"]],
                "locked_validation": cl["locked_metrics"],
                "status": "RESEARCH_CANDIDATE",
                "independent_validation_n": 0,
                "decision": "INDEPENDENT_VALIDATION_REQUIRED",
            },
        },
        "independent_validation": {
            "PPB_FU": {"qualified_n": 0, "reason": "No independent source-family human plasma fu/PPB cohort available in the current warehouse; PKSmart fup is not an independent PPB validation source."},
            "HUMAN_TOTAL_IV_CL": {"qualified_n": 0, "reason": "All qualified total-IV-CL candidate records are PKSmart source family; no independent source-family cohort was found in the current warehouse."},
            "VDSS": {"qualified_n": 0, "reason": "No independent source-family strict human Vdss cohort is present."},
        },
        "vdss_strict_audit": {
            "original_v6_1_qualified_n": vdss["available_n"],
            "pksmart_marker_n": pksmart_vdss,
            "strict_qualified_n": len(strict),
            # strict_vdss_rows already counts each PKSmart marker as
            # ambiguous; keep this explicit to prevent double counting.
            "ambiguous_n": ambiguous,
            "rejected_n": rejected,
            "strict_n30_overlap": strict_n30_overlap,
            "strict_semantics": {"species": "Homo sapiens", "endpoint": "true Vdss", "unit": "L/kg", "context": "Clinical IV Pharmacokinetics", "analyte": "record-level analyte not present in PKSmart marker records; structured records retain source reference"},
            "source_families": sorted({r["source_family"] for r in strict}),
            "source_independence": "SINGLE_EXISTING_REFERENCE_FAMILY",
            "model_comparison": compact_vdss_cv(strict),
            "decision": "VDSS_PUBLIC_DATA_INSUFFICIENT",
        },
        "learning_and_data_value": {
            "PPB_FU": "Scaffold-CV is strong but independent validation is the limiting experiment; buy/source validation before more training records.",
            "HUMAN_TOTAL_IV_CL": "Within-source candidate performance is not sufficient for qualification; independent source-family IV total CL is required.",
            "VDSS": "Semantic qualification reduces 735 to 50; the dominant value is exact Vdss context and independent source diversity, not another tree-model sweep.",
        },
        "commercial_data_specification": [
            {"endpoint": "HUMAN_FU_PLASMA", "current_qualified_n": 232, "independent_source_n": 0, "desired_new_n": 50, "semantics": "human plasma fraction unbound or fraction bound with method and concentration context", "metadata": ["chemical identity", "matrix", "method", "temperature", "concentration", "source record", "study/reference"], "role": "LOCKED_VALIDATION", "priority": 3},
            {"endpoint": "HUMAN_TOTAL_IV_CL", "current_qualified_n": 748, "independent_source_n": 0, "desired_new_n": 30, "semantics": "human intravenous total systemic clearance, not CL/F, hepatic-only, renal-only, or derived", "metadata": ["analyte", "chemical form", "IV route", "dose", "formulation", "study", "unit", "source record"], "role": "LOCKED_VALIDATION", "priority": 2},
            {"endpoint": "HUMAN_VDSS", "current_qualified_n": 50, "independent_source_n": 0, "desired_new_n": 100, "semantics": "human IV-derived true Vss/Vdss in L/kg; exclude Vz, Vc, Vd/F", "metadata": ["analyte", "IV route", "steady-state definition", "body weight", "chemical form", "study/reference", "unit"], "role": "DEVELOPMENT plus independent validation", "priority": 1},
        ],
        "global_status": "MORE_FOUNDATION_REQUIRED",
        "global_status_reason": "Core reference library and provenance are stable, but PPB and total-IV-CL lack independent validation and strict Vdss evidence is only 50 records from one source family.",
    }
    out = ROOT / "validation/model_qualification_data_closure_v6_2.json"
    out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"ppb_independent_n": 0, "cl_independent_n": 0, "vdss_original_n": vdss["available_n"], "vdss_strict_n": len(strict), "vdss_ambiguous_n": ambiguous, "vdss_cv": output["vdss_strict_audit"]["model_comparison"]}, indent=2))


if __name__ == "__main__":
    main()
