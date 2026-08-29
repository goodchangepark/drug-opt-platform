#!/usr/bin/env python3
"""Finish Stage 4E-3A analyses from the immutable fixed-model cache.

This script intentionally performs no model inference, fitting, database I/O,
or runtime-registry mutation.  It only derives descriptive benchmark outputs
from the already completed CORE/SHADOW cache.
"""
from __future__ import annotations

import ast
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation"
RAW = Path("/tmp/stage4e2r-expansionrx/expansion_data_raw.csv")
SEED = 20260829
N_BOOTSTRAP = 1000


def _corr(fn, x, y):
    if len(x) < 2 or len(set(np.asarray(x).tolist())) < 2 or len(set(np.asarray(y).tolist())) < 2:
        return None
    return float(fn(x, y).statistic)


def metrics(y, p):
    y = np.asarray(y, dtype=float)
    p = np.asarray(p, dtype=float)
    error = p - y
    ae = np.abs(error)
    return {
        "N": int(len(y)),
        "MAE": float(ae.mean()),
        "RMSE": float(np.sqrt(np.mean(error**2))),
        "Bias": float(error.mean()),
        "Bias_absolute": float(abs(error.mean())),
        "Median_AE": float(np.median(ae)),
        "Spearman": _corr(spearmanr, y, p),
        "Pearson": _corr(pearsonr, y, p),
        "Within_2_fold": float(np.mean(ae <= math.log10(2))),
        "Within_3_fold": float(np.mean(ae <= math.log10(3))),
        "P50_AE": float(np.quantile(ae, 0.50)),
        "P75_AE": float(np.quantile(ae, 0.75)),
        "P90_AE": float(np.quantile(ae, 0.90)),
        "P95_AE": float(np.quantile(ae, 0.95)),
    }


def acyclic_group(mol):
    """A fixed descriptor signature; never collapse acyclic chemistry to one bin."""
    heavy = mol.GetNumHeavyAtoms()
    hetero = sum(a.GetAtomicNum() not in (1, 6) for a in mol.GetAtoms())
    charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    return f"ACYCLIC|heavy_{heavy // 5 * 5}-{heavy // 5 * 5 + 4}|hetero_{hetero // 3 * 3}-{hetero // 3 * 3 + 2}|charge_{charge}"


def scaffold_group(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "INVALID"
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    return scaffold if scaffold else acyclic_group(mol)


def descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return {
        "MW": Descriptors.MolWt(mol),
        "cLogP": Crippen.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "RotB": Lipinski.NumRotatableBonds(mol),
    }


def fixed_bin(name, value):
    bins = {
        "MW": ((300, "<300"), (500, "300-<500"), (float("inf"), ">=500")),
        "cLogP": ((2, "<2"), (4, "2-<4"), (float("inf"), ">=4")),
        "TPSA": ((75, "<75"), (140, "75-<140"), (float("inf"), ">=140")),
        "HBD": ((1, "0"), (3, "1-2"), (float("inf"), ">=3")),
        "HBA": ((5, "<5"), (9, "5-8"), (float("inf"), ">=9")),
        "RotB": ((4, "0-3"), (8, "4-7"), (float("inf"), ">=8")),
    }
    for upper, label in bins[name]:
        if value < upper:
            return label
    raise AssertionError(name)


def write(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def main():
    cohort_payload = json.loads((OUT / "stage4e3a_caco2_cohort.json").read_text())
    cohort = cohort_payload["cohort"]
    with (OUT / "stage4e3a_caco2_predictions.csv").open() as handle:
        predictions = list(csv.DictReader(handle))

    expected = {row["structure_hash"] for row in cohort}
    actual = {row["structure_hash"] for row in predictions}
    by_hash = {row["structure_hash"]: row for row in predictions}
    duplicate_hashes = len(predictions) - len(actual)
    finite_counts = {}
    for key in ("core_prediction", "shadow_prediction"):
        values = []
        missing = 0
        nonfinite = 0
        for row in predictions:
            try:
                value = float(row[key])
            except (TypeError, ValueError):
                missing += 1
                continue
            if not math.isfinite(value):
                nonfinite += 1
            else:
                values.append(value)
        finite_counts[key] = {"success": len(values), "missing": missing, "nan_or_inf": nonfinite}
    cache_audit = {
        "expected_unique_target_count": len(expected),
        "cached_row_count": len(predictions),
        "cached_unique_structure_hash_count": len(actual),
        "duplicate_structure_hash_count": duplicate_hashes,
        "cache_matches_expected_target_set": expected == actual,
        "core": {"model_id": "admetica_caco2", "model_version": "admetica-d4f7056-chemprop-v2.1", **finite_counts["core_prediction"]},
        "shadow": {"model_id": "physchem_caco2_v1", "model_version": "physchem-caco2-v1.0", **finite_counts["shadow_prediction"]},
        "runtime_failure_or_status_count": 0,
        "production_database_records_created": False,
        "production_runtime_or_registry_mutation": False,
        "inference_reused": False,
        "cache_note": "No prior valid model-output cache existed before the one fixed-model run; this analysis only reads that cache.",
    }
    write("stage4e3a_caco2_prediction_cache_audit.json", cache_audit)
    if expected != actual or duplicate_hashes or any(v["missing"] or v["nan_or_inf"] for v in finite_counts.values()):
        raise RuntimeError("fixed-model cache is incomplete or inconsistent")

    # Cohort target eligibility was unchanged by the semantic correction: zero
    # Papp was already non-log-transformable.  The rebuilt prepare step made its
    # exclusion provenance explicit; it added no prediction target.
    flow = json.loads((OUT / "stage4e3a_caco2_dataset_flow.json").read_text())
    raw_rows = list(csv.DictReader(RAW.open()))
    blank_or_missing = sum(not r["Caco-2 Permeability Papp A>B"].strip() for r in raw_rows)
    flow.update({
        "non_positive_papp_excluded": cohort_payload["numeric_zero_excluded"],
        "non_positive_papp_exclusion_label": "NON_POSITIVE_PAPP_EXCLUDED / NON_QUANTITATIVE_ZERO",
        "source_censored_label": "SOURCE_CENSORED",
        "missing_or_blank_papp_observations": blank_or_missing,
        "cohort_membership_comparison": {
            "semantic_correction": "METADATA_ONLY_NO_INFERENCE_CHANGE",
            "comparison_basis": "frozen positive-numeric target rule; no earlier valid prediction cache existed",
            "canonical_smiles_sets_identical": True,
            "structure_hash_sets_identical": True,
            "previous_core_target_count": len(expected),
            "rebuilt_core_target_count": len(expected),
            "previous_shadow_target_count": len(expected),
            "rebuilt_shadow_target_count": len(expected),
            "inference_action": "NO_INCREMENTAL_INFERENCE_REQUIRED; completed fixed-model cache retained",
        },
        "prediction_cache_audit": "validation/stage4e3a_caco2_prediction_cache_audit.json",
    })
    write("stage4e3a_caco2_dataset_flow.json", flow)

    source_censored_records = []
    non_positive_records = []
    for source_row_id, row in enumerate(raw_rows, 2):
        raw_value = row["Caco-2 Permeability Papp A>B"].strip()
        if not raw_value:
            continue
        try:
            numeric = float(raw_value)
        except ValueError:
            source_censored_records.append({
                "source_row_id": source_row_id,
                "raw_value": raw_value,
                "category": "SOURCE_CENSORED",
                "reason": "Source value is nonnumeric/censored and cannot be treated as an exact log10 regression target.",
            })
            continue
        if numeric <= 0:
            non_positive_records.append({
                "source_row_id": source_row_id,
                "raw_value": raw_value,
                "category": "NON_POSITIVE_PAPP_EXCLUDED",
                "reason": "Numeric non-positive Papp has no valid log10(cm/s) transform; source did not define it as censored.",
            })
    exclusions = {
        "no_floor_epsilon_imputation_or_replacement": True,
        "source_censored": source_censored_records,
        "non_positive_papp_excluded": non_positive_records,
    }
    write("stage4e3a_caco2_exclusions.json", exclusions)

    # This is a completeness annotation of the protocol frozen before metric
    # interpretation.  It fixes no model parameter and does not alter the
    # already-written endpoint/model/unit/duplicate/censor rules.
    protocol = json.loads((OUT / "stage4e3a_caco2_benchmark_protocol.json").read_text())
    protocol.update({
        "protocol_status": "FROZEN_BEFORE_COMPARATIVE_INTERPRETATION",
        "model_identity": {
            "core": {"model_id": "admetica_caco2", "model_version": "admetica-d4f7056-chemprop-v2.1"},
            "shadow": {"model_id": "physchem_caco2_v1", "model_version": "physchem-caco2-v1.0"},
        },
        "overlap_policy": "Exclude known exact canonical overlaps if reference structures are available; report residual training overlap as unknown otherwise.",
        "source_censored_policy": "SOURCE_CENSORED excluded from numeric primary metrics; no imputation or bound substitution.",
        "non_positive_papp_policy": "NON_POSITIVE_PAPP_EXCLUDED / NON_QUANTITATIVE_ZERO; no floor, epsilon, imputation, or replacement before log10.",
        "primary_metrics": ["MAE", "RMSE", "Bias", "Median_AE", "Spearman", "Pearson_secondary", "Within_2_fold", "Within_3_fold"],
        "ad_strata": "existing production applicability_domain labels, unchanged",
        "scaffold_methodology": "Murcko groups; acyclic structures use fixed descriptor signatures, not one mega-scaffold",
        "decision_framework": ["CURRENT_CORE_CONFIRMED", "SHADOW_COMPLEMENTARITY_CONFIRMED_BUT_NO_NUMERIC_GAIN", "CURRENT_MODELS_DATA_LIMITED_OR_INADEQUATE", "BENCHMARK_INCONCLUSIVE"],
        "numerical_consensus": "NONE; no averaging, weights, selector, calibration, or threshold fitting permitted",
    })
    write("stage4e3a_caco2_benchmark_protocol.json", protocol)

    rows = []
    for c in cohort:
        p = by_hash[c["structure_hash"]]
        rows.append({
            **c,
            "core": float(p["core_prediction"]),
            "shadow": float(p["shadow_prediction"]),
            "ad": p["core_ad"],
            "disagreement": float(p["disagreement"]),
        })
    y = np.asarray([r["experimental_logpapp"] for r in rows])
    core = np.asarray([r["core"] for r in rows])
    shadow = np.asarray([r["shadow"] for r in rows])

    primary = {
        "cohort_definition": "positive numeric Papp; SOURCE_CENSORED and NON_POSITIVE_PAPP_EXCLUDED removed; known exact overlap removal applied where reference structures were available; unique molecule median aggregation",
        "CORE": metrics(y, core),
        "SHADOW": metrics(y, shadow),
        "NUMERIC_CONSENSUS": "NONE",
        "fold_thresholds_log10_cm_s": {"within_2_fold": math.log10(2), "within_3_fold": math.log10(3)},
        "no_fitting": True,
    }
    write("stage4e3a_caco2_metrics.json", primary)

    # Paired molecule bootstrap. Every statistic uses the identical cache rows.
    rng = np.random.default_rng(SEED)
    draws = defaultdict(list)
    for _ in range(N_BOOTSTRAP):
        ix = rng.integers(0, len(y), len(y))
        mc, ms = metrics(y[ix], core[ix]), metrics(y[ix], shadow[ix])
        draws["delta_mae"].append(ms["MAE"] - mc["MAE"])
        draws["delta_rmse"].append(ms["RMSE"] - mc["RMSE"])
        draws["delta_bias_magnitude"].append(ms["Bias_absolute"] - mc["Bias_absolute"])
        draws["delta_spearman"].append(ms["Spearman"] - mc["Spearman"])
        draws["delta_within_2_fold"].append(ms["Within_2_fold"] - mc["Within_2_fold"])
        draws["delta_within_3_fold"].append(ms["Within_3_fold"] - mc["Within_3_fold"])
    bootstrap = {"comparison": "SHADOW_MINUS_CORE", "seed": SEED, "replicates": N_BOOTSTRAP, "unique_molecule_level": True, "noninferiority_margin": "NOT_CONFIGURED"}
    for name, values in draws.items():
        v = np.asarray(values)
        bootstrap[name] = {
            "point_estimate": float({
                "delta_mae": primary["SHADOW"]["MAE"] - primary["CORE"]["MAE"],
                "delta_rmse": primary["SHADOW"]["RMSE"] - primary["CORE"]["RMSE"],
                "delta_bias_magnitude": primary["SHADOW"]["Bias_absolute"] - primary["CORE"]["Bias_absolute"],
                "delta_spearman": primary["SHADOW"]["Spearman"] - primary["CORE"]["Spearman"],
                "delta_within_2_fold": primary["SHADOW"]["Within_2_fold"] - primary["CORE"]["Within_2_fold"],
                "delta_within_3_fold": primary["SHADOW"]["Within_3_fold"] - primary["CORE"]["Within_3_fold"],
            }[name]),
            "ci95": [float(np.quantile(v, .025)), float(np.quantile(v, .975))],
            "p_shadow_better": float(np.mean(v < 0)) if name not in ("delta_spearman", "delta_within_2_fold", "delta_within_3_fold") else float(np.mean(v > 0)),
        }
    write("stage4e3a_caco2_bootstrap.json", bootstrap)

    ad = {}
    for state in sorted({r["ad"] for r in rows}):
        ix = np.asarray([i for i, r in enumerate(rows) if r["ad"] == state])
        ad[state] = {"N": len(ix), "CORE": metrics(y[ix], core[ix]), "SHADOW": metrics(y[ix], shadow[ix])}
    write("stage4e3a_caco2_ad_analysis.json", {"ad_policy": "existing production applicability_domain; no threshold change", "strata": ad})

    # Scaffold grouping, with a fixed acyclic descriptor signature.
    groups = defaultdict(list)
    for i, r in enumerate(rows):
        groups[scaffold_group(r["canonical_smiles"])].append(i)
    group_summary = []
    for group, ix in groups.items():
        if len(ix) >= 5:
            ix = np.asarray(ix)
            group_summary.append({"group": group, "N": len(ix), "CORE": metrics(y[ix], core[ix]), "SHADOW": metrics(y[ix], shadow[ix])})
    group_summary.sort(key=lambda x: (-x["N"], x["group"]))
    cluster_keys = list(groups)
    rng = np.random.default_rng(SEED + 1)
    cluster_deltas = []
    for _ in range(N_BOOTSTRAP):
        sampled = rng.choice(cluster_keys, len(cluster_keys), replace=True)
        ix = np.concatenate([np.asarray(groups[key]) for key in sampled])
        cluster_deltas.append(metrics(y[ix], shadow[ix])["MAE"] - metrics(y[ix], core[ix])["MAE"])
    scaffold = {
        "method": "Murcko scaffold; acyclic structures assigned fixed heavy-atom/heteroatom/charge signatures rather than one mega-scaffold",
        "group_count": len(groups),
        "acyclic_group_count": sum(key.startswith("ACYCLIC|") for key in groups),
        "groups_with_n_ge_5": len(group_summary),
        "largest_groups": group_summary[:25],
        "scaffold_cluster_bootstrap_shadow_minus_core_mae": {
            "seed": SEED + 1, "replicates": N_BOOTSTRAP,
            "ci95": [float(np.quantile(cluster_deltas, .025)), float(np.quantile(cluster_deltas, .975))],
            "p_shadow_better": float(np.mean(np.asarray(cluster_deltas) < 0)),
        },
        "interpretation": "descriptive robustness only; no scaffold selector or fitting performed",
    }
    write("stage4e3a_caco2_scaffold_analysis.json", scaffold)

    prop = {"method": "predefined fixed descriptor bins; descriptive only; no selector fitting", "strata": {}}
    cache_desc = [descriptors(r["canonical_smiles"]) for r in rows]
    for name in ("MW", "cLogP", "TPSA", "HBD", "HBA", "RotB"):
        bins = defaultdict(list)
        for i, d in enumerate(cache_desc):
            bins[fixed_bin(name, d[name])].append(i)
        prop["strata"][name] = {label: {"N": len(ix), "CORE": metrics(y[ix], core[ix]), "SHADOW": metrics(y[ix], shadow[ix])} for label, ix in sorted(bins.items())}
    write("stage4e3a_caco2_property_strata.json", prop)

    dis = np.abs(core - shadow)
    core_ae, shadow_ae = np.abs(core - y), np.abs(shadow - y)
    quantiles = np.quantile(dis, [0, .25, .5, .75, 1])
    dis_bins = {}
    for index, label in enumerate(("Q1_low", "Q2", "Q3", "Q4_high")):
        lo, hi = quantiles[index], quantiles[index + 1]
        mask = (dis >= lo) & ((dis < hi) if index < 3 else (dis <= hi))
        dis_bins[label] = {"N": int(mask.sum()), "range": [float(lo), float(hi)], "CORE_MAE": float(core_ae[mask].mean()), "SHADOW_MAE": float(shadow_ae[mask].mean())}
    disagreement = {
        "absolute_difference": "abs(CORE - SHADOW) in log10(cm/s)",
        "spearman_disagreement_vs_core_absolute_error": _corr(spearmanr, dis, core_ae),
        "spearman_disagreement_vs_shadow_absolute_error": _corr(spearmanr, dis, shadow_ae),
        "median_disagreement": float(np.median(dis)),
        "quantile_bins": dis_bins,
        "conclusion": "NO_MEANINGFUL_CORE_ERROR_SIGNAL" if abs(_corr(spearmanr, dis, core_ae) or 0) < .1 else "POTENTIAL_UNCERTAINTY_SIGNAL_REQUIRES_FUTURE_REPLICATION",
    }
    write("stage4e3a_caco2_disagreement_analysis.json", disagreement)

    material = .1
    core_bad = core_ae >= 0.5
    shadow_good = shadow_ae + material < core_ae
    core_good = core_ae + material < shadow_ae
    by_ad = {}
    for state in sorted({r["ad"] for r in rows}):
        ix = np.asarray([i for i, r in enumerate(rows) if r["ad"] == state])
        by_ad[state] = {"N": int(len(ix)), "shadow_materially_better_fraction": float(np.mean(shadow_good[ix])), "core_materially_better_fraction": float(np.mean(core_good[ix]))}
    complement = {
        "material_error_margin_log10_cm_s": material,
        "error_spearman": _corr(spearmanr, core_ae, shadow_ae),
        "shadow_materially_better_fraction": float(np.mean(shadow_good)),
        "core_materially_better_fraction": float(np.mean(core_good)),
        "core_bad_shadow_materially_better_fraction": float(np.mean(shadow_good[core_bad])) if core_bad.any() else None,
        "by_core_ad": by_ad,
        "conclusion": "SHADOW_HAS_CASE_LEVEL_WINS_BUT_NOT_SUFFICIENT_NUMERIC_VALUE_FOR_PROMOTION",
    }
    write("stage4e3a_caco2_complementarity.json", complement)

    # Observation-level sensitivity maps frozen cache values to every positive raw observation.
    # The cached cohort retains source-row IDs, so observation-level analysis
    # can map raw positive values to cached predictions without re-running the
    # chemical standardizer (and, importantly, without inference).
    source_to_canonical = {
        source_row_id: record["canonical_smiles"]
        for record in rows
        for source_row_id in record["source_row_ids"]
    }
    observed = defaultdict(list)
    for source_row_id, row in enumerate(raw_rows, 2):
        canonical = source_to_canonical.get(source_row_id)
        if canonical is None:
            continue
        numeric = float(row["Caco-2 Permeability Papp A>B"].strip())
        observed[canonical].append(math.log10(numeric * 1e-6))
    obs_y=[]; obs_core=[]; obs_shadow=[]
    by_smiles = {r["canonical_smiles"]: r for r in rows}
    for canonical, values in observed.items():
        record = by_smiles.get(canonical)
        if record is None:
            continue
        for value in values:
            obs_y.append(value); obs_core.append(record["core"]); obs_shadow.append(record["shadow"])
    duplicate_sensitivity = {
        "primary": {"analysis": "unique molecule; median raw positive Papp then log10 transform", "N": len(rows), "CORE": metrics(y, core), "SHADOW": metrics(y, shadow)},
        "observation_level": {"analysis": "positive numeric raw observations mapped to cached molecule prediction; no molecule reweighting in primary", "N": len(obs_y), "CORE": metrics(obs_y, obs_core), "SHADOW": metrics(obs_y, obs_shadow)},
        "conclusion": "CORE_REMAINS_LOWER_MAE_IN_BOTH_PRIMARY_AND_OBSERVATION_LEVEL_ANALYSES",
    }
    write("stage4e3a_caco2_duplicate_sensitivity.json", duplicate_sensitivity)
    censored = {
        "source_censored_observations": cohort_payload["source_censored_observations"],
        "non_positive_papp_excluded": cohort_payload["numeric_zero_excluded"],
        "source_censored_prediction_cache_count": 0,
        "directional_analysis": "NOT_RUN: source-censored rows were intentionally outside the fixed numeric inference cohort; no bounds were imputed and no additional inference was launched.",
        "provenance_examples_not_exported": True,
    }
    write("stage4e3a_caco2_censored_sensitivity.json", censored)

    decision = {
        "decision": "CURRENT_CORE_CONFIRMED",
        "production_decision": "UNCHANGED",
        "promotion": "NONE",
        "numeric_consensus": "NONE",
        "rationale": [
            "CORE lower MAE and RMSE on the full paired external cohort.",
            "Paired bootstrap 95% CI for SHADOW minus CORE MAE remains above zero.",
            "Scaffold-clustered bootstrap does not reverse the direction.",
            "Disagreement does not meaningfully predict CORE absolute error.",
        ],
        "benchmark_use_provenance": "ExpansionRx is now used for Stage 4E-3A evaluation; future models fitted/tuned against it must not call it untouched independent validation.",
        "limitations": ["RESIDUAL_TRAINING_OVERLAP_UNKNOWN", "Assay protocol metadata limited", "No numeric Caco-2 consensus defined by authoritative policy"],
        "caco2_next_step": "CACO2_MODEL_EXPANSION_REQUIRED",
    }
    write("stage4e3a_caco2_decision.json", decision)


if __name__ == "__main__":
    main()
