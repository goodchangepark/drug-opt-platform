#!/usr/bin/env python3
"""Leakage-safe v3.3.3 research ladders for hepatic CL, Vdss and total IV CL.

The three contracts are intentionally separate.  Hepatic CL is an assisted
N=30 mechanistic benchmark, Vdss uses only strict true-Vss observations, and
total CL uses the frozen PKSmart within-source cohort.  None can promote a
production route without the independent gates recorded in the artifact.
"""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.clinical_pk_cohort import CLINICAL_PK_COHORT  # noqa: E402
from backend.ionization import IonizationClass, analyze_ionization  # noqa: E402
from backend.pk_parameter_set import calculate_well_stirred_clearance, canonical_fu_from_ppb  # noqa: E402
from scripts import evaluate_human_total_iv_cl_v56 as base  # noqa: E402
from scripts.run_endpoint_rebuild_v61 import clean, rows as library_rows  # noqa: E402
from scripts.run_model_qualification_v62 import strict_vdss_rows  # noqa: E402

OUT = ROOT / "validation/pk_core_hybrid_ladders_v78.json"
LIBRARY_PATH = ROOT / "validation/reference_library_v1_1000.json"
SEED = 78
RDLogger.DisableLog("rdApp.warning")


def metrics(y: np.ndarray, pred: np.ndarray) -> dict:
    y, pred = np.asarray(y, float), np.clip(np.asarray(pred, float), 1e-6, None)
    fold = np.maximum(pred / y, y / pred)
    return {
        "n": int(len(y)), "aafe": float(np.mean(fold)),
        "median_fold_error": float(np.median(fold)),
        "within_2_fold_pct": float(np.mean(fold <= 2) * 100),
        "within_3_fold_pct": float(np.mean(fold <= 3) * 100),
        "mae_log10": float(np.mean(np.abs(np.log10(pred) - np.log10(y)))),
        "bias_log10": float(np.mean(np.log10(pred) - np.log10(y))),
    }


def groups(rows: list[dict]) -> np.ndarray:
    return np.asarray([
        MurckoScaffold.MurckoScaffoldSmiles(mol=row["mol"], includeChirality=False)
        or row.get("canonical_smiles") or row.get("name")
        for row in rows
    ])


def similarity(train: list[dict], test: list[dict]) -> np.ndarray:
    train_fp = [base.FINGERPRINT.GetFingerprint(row["mol"]) for row in train]
    return np.asarray([
        DataStructs.BulkTanimotoSimilarity(base.FINGERPRINT.GetFingerprint(row["mol"]), train_fp)
        for row in test
    ], dtype=float)


def local_average(values: np.ndarray, similarities: np.ndarray, k: int = 8) -> np.ndarray:
    k = min(k, len(values))
    selected = np.argpartition(similarities, -k, axis=1)[:, -k:]
    return np.asarray([
        np.average(values[index], weights=np.maximum(similarities[i, index], 0.05) ** 2)
        for i, index in enumerate(selected)
    ])


def hepatic_ladder() -> dict:
    rows = []
    for drug in CLINICAL_PK_COHORT:
        if not (drug.cl_hepatic_ml_min_kg and drug.ppb_percent is not None and drug.cl_int_hlm_ml_min_kg):
            continue
        mol = Chem.MolFromSmiles(drug.smiles)
        fu = canonical_fu_from_ppb(drug.ppb_percent)
        prior = calculate_well_stirred_clearance(
            drug.cl_int_hlm_ml_min_kg, fu, rb=drug.blood_to_plasma_rb,
        )["cl_h_plasma_ml_min_kg"]
        rows.append({
            "name": drug.compound_name, "mol": mol,
            "y": float(drug.cl_hepatic_ml_min_kg), "fu": float(fu),
            "clint": float(drug.cl_int_hlm_ml_min_kg), "prior": max(float(prior), 1e-6),
            "extraction": drug.extraction_category,
        })
    rows.sort(key=lambda row: row["name"])
    y = np.asarray([row["y"] for row in rows])
    prior = np.asarray([row["prior"] for row in rows])
    desc = np.asarray([base.descriptor_vector(row["mol"]) for row in rows], float)
    dependencies = np.asarray([[np.log10(row["fu"]), np.log10(row["clint"])] for row in rows])
    splits = list(GroupKFold(5).split(desc, groups=groups(rows)))
    names = [
        "well_stirred", "structure_direct", "fu_only", "clint_only", "fu_clint",
        "mechanistic_structure_residual", "mechanistic_fu_residual",
        "mechanistic_clint_residual", "mechanistic_fu_clint_residual",
        "mechanistic_local_residual", "mechanistic_cluster_residual",
    ]
    pred = {name: np.empty(len(rows)) for name in names}
    pred["well_stirred"] = prior.copy()
    max_sim = np.empty(len(rows))
    for fold, (train, test) in enumerate(splits):
        train_rows, test_rows = [rows[i] for i in train], [rows[i] for i in test]
        residual = np.log10(y[train]) - np.log10(prior[train])
        candidates = {
            "structure_direct": (desc, np.log10(y), ExtraTreesRegressor(n_estimators=400, min_samples_leaf=3, max_features=.8, random_state=SEED, n_jobs=2), False),
            "fu_only": (dependencies[:, :1], np.log10(y), make_pipeline(StandardScaler(), Ridge(alpha=3.0)), False),
            "clint_only": (dependencies[:, 1:], np.log10(y), make_pipeline(StandardScaler(), Ridge(alpha=3.0)), False),
            "fu_clint": (dependencies, np.log10(y), ExtraTreesRegressor(n_estimators=300, min_samples_leaf=3, random_state=SEED, n_jobs=2), False),
            "mechanistic_structure_residual": (desc, residual, ExtraTreesRegressor(n_estimators=400, min_samples_leaf=3, max_features=.8, random_state=SEED, n_jobs=2), True),
            "mechanistic_fu_residual": (dependencies[:, :1], residual, make_pipeline(StandardScaler(), Ridge(alpha=3.0)), True),
            "mechanistic_clint_residual": (dependencies[:, 1:], residual, make_pipeline(StandardScaler(), Ridge(alpha=3.0)), True),
            "mechanistic_fu_clint_residual": (dependencies, residual, ExtraTreesRegressor(n_estimators=300, min_samples_leaf=3, random_state=SEED, n_jobs=2), True),
        }
        for name, (x, target, model, residual_mode) in candidates.items():
            model.fit(x[train], target if residual_mode else target[train])
            value = model.predict(x[test])
            pred[name][test] = 10 ** (np.log10(prior[test]) + value) if residual_mode else 10 ** value
        similarities = similarity(train_rows, test_rows)
        max_sim[test] = similarities.max(axis=1)
        correction = local_average(residual, similarities, 8)
        pred["mechanistic_local_residual"][test] = 10 ** (np.log10(prior[test]) + correction)
        scaler = StandardScaler().fit(desc[train])
        train_scaled, test_scaled = scaler.transform(desc[train]), scaler.transform(desc[test])
        cluster = KMeans(n_clusters=min(5, len(train) // 4), n_init=10, random_state=SEED + fold).fit(train_scaled)
        train_labels, test_labels = cluster.labels_, cluster.predict(test_scaled)
        correction_by_cluster = {
            label: float(np.median(residual[train_labels == label])) for label in set(train_labels)
        }
        pred["mechanistic_cluster_residual"][test] = 10 ** (
            np.log10(prior[test]) + np.asarray([correction_by_cluster[label] for label in test_labels])
        )
    results = {}
    for name, values in pred.items():
        results[name] = {
            "pooled_oof": metrics(y, values),
            "by_extraction": {
                label: metrics(y[ix], values[ix])
                for label in ("LOW", "INTERMEDIATE", "HIGH")
                if (ix := np.asarray([i for i, row in enumerate(rows) if row["extraction"] == label])).size
            },
            "by_ad": {
                label: metrics(y[mask], values[mask])
                for label, mask in {
                    "OOD": max_sim < .30, "BORDERLINE": (max_sim >= .30) & (max_sim < .50), "IN_DOMAIN": max_sim >= .50,
                }.items() if mask.any()
            } if name != "well_stirred" else {},
        }
    baseline = results["well_stirred"]["pooled_oof"]
    selected = min(results, key=lambda name: results[name]["pooled_oof"]["aafe"])
    candidate = results[selected]["pooled_oof"]
    return {
        "canonical_contract": {
            "endpoint": "HUMAN_HEPATIC_CL", "species": "HUMAN", "unit": "mL/min/kg",
            "dataset": "curated N=30 clinical hepatic-clearance cohort",
            "dataset_hash": hashlib.sha256(json.dumps([(r["name"], r["y"], r["fu"], r["clint"]) for r in rows]).encode()).hexdigest(),
            "mode": "ASSISTED", "inputs": "observed fu and HLM Clint", "split": "5-fold Murcko scaffold GroupKFold",
        },
        "results": results, "baseline": baseline, "best_candidate": selected,
        "best_metrics": candidate,
        "decision": "RESEARCH_CANDIDATE" if candidate["aafe"] <= baseline["aafe"] * .90 else "REJECT_NO_MATERIAL_GAIN",
        "promotion": False,
        "blockers": ["Assisted rather than full-prediction benchmark", "N=30", "No independent hepatic-clearance cohort"],
    }


def _fu_for_library_compound(compound: dict) -> float | None:
    evidence = compound.get("evidence")
    if isinstance(evidence, dict) and evidence.get("HUMAN_FUP") is not None:
        return float(evidence["HUMAN_FUP"])
    observations = compound.get("observations", evidence if isinstance(evidence, list) else [])
    for record in observations or []:
        if record.get("canonical_endpoint_id") == "HUMAN_PPB" and record.get("normalized_value") is not None:
            bound = float(record["normalized_value"])
            if record.get("normalized_unit") == "% bound" and 0 <= bound < 100:
                return max((100.0 - bound) / 100.0, 1e-4)
    return None


def mechanistic_vdss(row: dict, fu: float) -> tuple[float, list[float]]:
    ion = analyze_ionization(row["canonical_smiles"])
    ion_class = ion.get("ionization_class", IonizationClass.NEUTRAL)
    clogp = float(ion.get("clogp") or 0.0)
    logd = float(ion.get("physiological_state_7_4", {}).get("estimated_logd74", clogp))
    lipo = np.clip(logd, -1.5, 4.5)
    if ion_class == IonizationClass.ACID:
        value = np.clip(.08 + .15 * fu + .05 * fu * 10 ** (.2 * lipo), .05, 1.5)
    elif ion_class == IonizationClass.BASE:
        value = np.clip(.6 + .4 * fu + .30 * fu * 10 ** (.35 * lipo), .2, 30)
    elif ion_class in (IonizationClass.ZWITTERION_POSSIBLE, IonizationClass.AMPHOLYTE):
        value = np.clip(.3 + .2 * fu + .10 * fu * 10 ** (.2 * lipo), .1, 5)
    else:
        value = np.clip(.6 + .4 * fu + .15 * fu * 10 ** (.3 * lipo), .1, 20)
    class_bits = [float(ion_class == label) for label in (IonizationClass.ACID, IonizationClass.BASE, IonizationClass.ZWITTERION_POSSIBLE, IonizationClass.AMPHOLYTE)]
    return float(value), [fu, clogp, logd, *class_bits]


def vdss_ladder(library: dict) -> dict:
    strict, ambiguous, rejected = strict_vdss_rows(library["compounds"])
    by_key = {compound["inchikey"]: compound for compound in library["compounds"]}
    rows = [row for row in strict if _fu_for_library_compound(by_key[row["inchikey"]]) is not None]
    y = np.asarray([row["y"] for row in rows])
    desc = np.asarray([base.descriptor_vector(row["mol"]) for row in rows])
    priors, physiology = zip(*[
        mechanistic_vdss(row, _fu_for_library_compound(by_key[row["inchikey"]])) for row in rows
    ])
    priors, physiology = np.asarray(priors), np.asarray(physiology)
    x_hybrid = np.hstack([desc, physiology])
    splits = list(GroupKFold(min(5, len(set(groups(rows))))).split(desc, groups=groups(rows)))
    names = ["mechanistic_prior", "direct_structure", "direct_physchem_fu", "mechanistic_residual", "local_residual", "cluster_residual"]
    pred = {name: np.empty(len(rows)) for name in names}
    pred["mechanistic_prior"] = priors.copy()
    max_sim = np.empty(len(rows))
    for fold, (train, test) in enumerate(splits):
        residual = np.log10(y[train]) - np.log10(priors[train])
        direct = ExtraTreesRegressor(n_estimators=400, min_samples_leaf=3, max_features=.8, random_state=SEED, n_jobs=2)
        direct.fit(desc[train], np.log10(y[train])); pred["direct_structure"][test] = 10 ** direct.predict(desc[test])
        hybrid = ExtraTreesRegressor(n_estimators=400, min_samples_leaf=3, max_features=.8, random_state=SEED, n_jobs=2)
        hybrid.fit(x_hybrid[train], np.log10(y[train])); pred["direct_physchem_fu"][test] = 10 ** hybrid.predict(x_hybrid[test])
        residual_model = ExtraTreesRegressor(n_estimators=400, min_samples_leaf=3, max_features=.8, random_state=SEED, n_jobs=2)
        residual_model.fit(x_hybrid[train], residual); pred["mechanistic_residual"][test] = 10 ** (np.log10(priors[test]) + residual_model.predict(x_hybrid[test]))
        sims = similarity([rows[i] for i in train], [rows[i] for i in test]); max_sim[test] = sims.max(axis=1)
        pred["local_residual"][test] = 10 ** (np.log10(priors[test]) + local_average(residual, sims, 6))
        scaler = StandardScaler().fit(x_hybrid[train]); scaled_train = scaler.transform(x_hybrid[train])
        cluster = KMeans(n_clusters=min(5, len(train) // 4), n_init=10, random_state=SEED + fold).fit(scaled_train)
        corrections = {label: float(np.median(residual[cluster.labels_ == label])) for label in set(cluster.labels_)}
        pred["cluster_residual"][test] = 10 ** (np.log10(priors[test]) + np.asarray([corrections[label] for label in cluster.predict(scaler.transform(x_hybrid[test]))]))
    results = {name: metrics(y, values) for name, values in pred.items()}
    best = min(results, key=lambda name: results[name]["aafe"])
    uncertainty = np.std(np.column_stack([np.log10(values) for name, values in pred.items() if name != "mechanistic_prior"]), axis=1)
    error = np.abs(np.log10(pred[best]) - np.log10(y))
    return {
        "canonical_contract": {
            "endpoint": "HUMAN_VDSS", "species": "HUMAN", "unit": "L/kg",
            "strict_semantics": "true Vss/Vdss only; Vz and Vd/F excluded", "qualified_n": len(strict),
            "hybrid_complete_case_n": len(rows), "ambiguous_n": ambiguous, "rejected_n": rejected,
            "split": "Murcko scaffold GroupKFold", "mode": "ASSISTED_FU",
        },
        "results": results, "best_candidate": best, "best_metrics": results[best],
        "uncertainty": {
            "ensemble_disagreement_error_spearman": float(spearmanr(uncertainty, error).statistic),
            "chemical_distance_error_spearman": float(spearmanr(1 - max_sim, error).statistic),
            "claim": "RESEARCH_DIAGNOSTIC_NOT_CALIBRATED",
        },
        "decision": "RESEARCH_ONLY_SINGLE_SOURCE", "promotion": False,
        "blockers": ["Only one source family", "No independent strict true-Vss cohort", "fu-assisted benchmark"],
    }


def total_cl_ladder() -> dict:
    rows = clean(library_rows("HUMAN_TOTAL_IV_CL", external_only=True))
    y = np.asarray([row["y"] for row in rows])
    log_y = np.log10(y)
    desc = base.build_matrix(rows, "descriptors")
    splits = list(GroupKFold(5).split(desc, groups=groups(rows)))
    names = ["extra_trees_v61", "random_forest", "hist_gradient", "fingerprint_knn", "local_residual", "cluster_residual"]
    pred = {name: np.empty(len(rows)) for name in names}
    max_sim = np.empty(len(rows))
    for fold, (train, test) in enumerate(splits):
        train_rows, test_rows = [rows[i] for i in train], [rows[i] for i in test]
        extra = ExtraTreesRegressor(n_estimators=400, min_samples_leaf=2, max_features=.8, random_state=SEED, n_jobs=2)
        extra.fit(desc[train], log_y[train]); global_pred = extra.predict(desc[test]); pred["extra_trees_v61"][test] = 10 ** global_pred
        forest = RandomForestRegressor(n_estimators=350, min_samples_leaf=3, max_features=.8, random_state=SEED, n_jobs=2)
        forest.fit(desc[train], log_y[train]); pred["random_forest"][test] = 10 ** forest.predict(desc[test])
        hist = HistGradientBoostingRegressor(loss="absolute_error", max_iter=240, max_leaf_nodes=31, l2_regularization=2, random_state=SEED)
        hist.fit(desc[train], log_y[train]); pred["hist_gradient"][test] = 10 ** hist.predict(desc[test])
        sims = similarity(train_rows, test_rows); max_sim[test] = sims.max(axis=1)
        neighbor = local_average(log_y[train], sims, 10); pred["fingerprint_knn"][test] = 10 ** neighbor
        # Base-model residuals are generated by an inner scaffold split.
        inner_groups = groups(train_rows); inner = GroupKFold(min(4, len(set(inner_groups))))
        inner_oof = np.empty(len(train))
        for inner_train, inner_test in inner.split(desc[train], groups=inner_groups):
            model = ExtraTreesRegressor(n_estimators=220, min_samples_leaf=2, max_features=.8, random_state=SEED, n_jobs=2)
            model.fit(desc[train][inner_train], log_y[train][inner_train]); inner_oof[inner_test] = model.predict(desc[train][inner_test])
        residual = log_y[train] - inner_oof
        pred["local_residual"][test] = 10 ** (global_pred + local_average(residual, sims, 10))
        scaler = StandardScaler().fit(desc[train]); train_scaled = scaler.transform(desc[train])
        cluster = KMeans(n_clusters=16, n_init=10, random_state=SEED + fold).fit(train_scaled)
        corrections = {label: float(np.median(residual[cluster.labels_ == label])) for label in set(cluster.labels_)}
        pred["cluster_residual"][test] = 10 ** (global_pred + np.asarray([corrections[label] for label in cluster.predict(scaler.transform(desc[test]))]))
    results = {name: metrics(y, values) for name, values in pred.items()}
    best = min(results, key=lambda name: results[name]["aafe"])
    error = np.abs(np.log10(pred[best]) - log_y)
    return {
        "canonical_contract": {
            "endpoint": "HUMAN_TOTAL_IV_CL", "species": "HUMAN", "route": "IV", "unit": "mL/min/kg",
            "dataset": "PKSmart compact endpoint markers", "dataset_n": len(rows), "split": "Murcko scaffold GroupKFold",
            "role": "WITHIN_SOURCE_RESEARCH_ONLY", "locked_stress_cohort_used_for_selection": False,
        },
        "results": results, "best_candidate": best, "best_metrics": results[best],
        "ad": {"chemical_distance_error_spearman": float(spearmanr(1 - max_sim, error).statistic), "claim": "RESEARCH_DIAGNOSTIC"},
        "decision": "RESEARCH_ONLY_NO_INDEPENDENT_SOURCE", "promotion": False,
        "blockers": ["Compact source markers lack record-level IV context", "No independent source-family validation"],
    }


def main() -> None:
    library = json.loads(LIBRARY_PATH.read_text(encoding="utf-8"))
    artifact = {
        "artifact": "PK_CORE_HYBRID_LADDERS_V78", "created_at": datetime.now(timezone.utc).isoformat(),
        "engine_release": "drugopt-prediction-engine-v3@3.3.3",
        "hepatic_cl": hepatic_ladder(), "vdss": vdss_ladder(library), "total_iv_cl": total_cl_ladder(),
        "production_decision": "UNCHANGED_ENDPOINT_ROUTING",
    }
    OUT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({name: {"best": result["best_candidate"], "metrics": result["best_metrics"], "decision": result["decision"]} for name, result in artifact.items() if isinstance(result, dict) and "best_candidate" in result}, indent=2))


if __name__ == "__main__":
    main()
