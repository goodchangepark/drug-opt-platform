"""
Validation Pair Recovery & Evidence Funnel Reconciliation (Drug-OPT v5.4).

Provides:
- Strict evidence funnel audit (Qualified -> Prediction Available -> Semantic Aligned -> Validated Pairs)
- Drop-off reason attribution (Clinical PK, Clinical Metabolic Balance, Target-Specific Activity, Non-pairable IC50)
- Scale/Unit reconciliation (e.g. Solubility logS vs uM scale alignment)
- Model training label contract enforcement for classification endpoints (no arbitrary IC50 binarization)
- Independent compound grouping (1 observation per compound per endpoint)
"""
from __future__ import annotations

import math
import numpy as np
from collections import defaultdict, Counter
from typing import Dict, List, Any, Optional, Tuple

from backend.database import SessionLocal
from backend.models import Compound, CompoundVersion, ExternalExperimentalEvidence
from backend.endpoint_contracts import get_endpoint_contract, OutputType
from backend.multimodel import get_v2_adapters_for_endpoint, get_model_adapter
from backend.endpoint_strategy_registry import get_all_strategies

EP_MAP = {
    "Solubility": "SOLUBILITY_GENERIC",
    "Permeability": "CACO2_PAPP_AB",
    "Plasma protein binding": "HUMAN_PPB",
    "HLM intrinsic clearance": "HLM_CLINT",
    "RLM intrinsic clearance": "RLM_CLINT",
    "MLM intrinsic clearance": "MLM_CLINT",
    "CYP1A2 inhibitor": "CYP1A2_INHIBITION",
    "CYP2C9 inhibitor": "CYP2C9_INHIBITION",
    "CYP2C19 inhibitor": "CYP2C19_INHIBITION",
    "CYP2D6 inhibitor": "CYP2D6_INHIBITION",
    "CYP3A4 inhibitor": "CYP3A4_INHIBITION",
    "CYP2C9 substrate": "CYP2C9_SUBSTRATE",
    "CYP2D6 substrate": "CYP2D6_SUBSTRATE",
    "CYP3A4 substrate": "CYP3A4_SUBSTRATE",
    "P-gp inhibitor": "PGP_INHIBITION",
    "hERG liability": "HERG_LIABILITY",
    "Ames mutagenicity": "AMES_MUTAGENICITY",
    "DILI clinical liability": "DILI_LIABILITY",
}


def audit_evidence_funnel() -> Dict[str, Any]:
    """Audits all persisted qualified external evidence records across the 4-stage funnel."""
    db = SessionLocal()
    try:
        records = db.query(ExternalExperimentalEvidence).filter(
            ExternalExperimentalEvidence.evidence_state.in_(["AUTO_QUALIFIED_EXTERNAL", "RELATED_EXTERNAL"])
        ).all()

        funnel_stats = {
            "total_qualified_records": len(records),
            "endpoints": {},
            "global_drop_reasons": Counter(),
        }

        # Group records by canonical endpoint
        ep_groups = defaultdict(list)
        for r in records:
            cv = db.query(CompoundVersion).filter(CompoundVersion.id == r.compound_version_id).first()
            if not cv:
                continue
            comp = db.query(Compound).filter(Compound.id == cv.compound_row_id).first()
            if not comp:
                continue
            ep_groups[r.canonical_endpoint_id or "UNRESOLVED"].append((comp, cv, r))

        for eid, items in ep_groups.items():
            q_n = len(items)
            p_n = 0
            s_n = 0
            indep_pairs = {}
            reasons = Counter()

            for comp, cv, r in items:
                smi = cv.canonical_smiles
                cname = comp.name

                # Check prediction availability
                if any(eid.startswith(prefix) for prefix in ("HUMAN_PK_", "RAT_PK_", "DOG_PK_", "MOUSE_PK_", "MONKEY_PK_", "UNSPECIFIED_PK_")):
                    reasons["CLINICAL_PK_NO_PREDICTIVE_ML_MODEL"] += 1
                    continue
                elif eid in ("CYP3A4_METABOLIC_CONTRIBUTION", "EXCRETION_FECAL", "EXCRETION_URINARY", "METABOLITE_OBSERVATION", "HEPATOCYTE_CLINT"):
                    reasons["CLINICAL_METABOLIC_BALANCE_NO_ML_REGRESSION"] += 1
                    continue
                elif eid.startswith("ACTIVITY_"):
                    reasons["TARGET_SPECIFIC_ACTIVITY_EXCLUDED_FROM_GENERAL_ADMET_ML"] += 1
                    continue

                p_n += 1

                # Check semantic and scale compatibility
                if eid == "SOLUBILITY_GENERIC":
                    if r.normalized_unit in ("log10(mol/L)", "logS"):
                        val_logs = float(r.normalized_value)
                        s_n += 1
                        indep_pairs[cname] = (smi, val_logs, "log10(mol/L)")
                    elif r.raw_unit in ("uM", "µM", "uM ", "µM "):
                        val_logs = math.log10(float(r.raw_value) * 1e-6)
                        s_n += 1
                        indep_pairs[cname] = (smi, val_logs, "log10(mol/L)")
                    else:
                        reasons["UNIT_CONTEXT_INCOMPATIBLE"] += 1
                elif eid == "HUMAN_PPB":
                    if r.normalized_unit in ("% bound", "%") and 0.0 <= float(r.normalized_value) <= 100.0 and r.raw_unit not in ("nM", "hours", "uM", "µM"):
                        s_n += 1
                        indep_pairs[cname] = (smi, float(r.normalized_value), "% bound")
                    else:
                        reasons["UNIT_CONTEXT_INCOMPATIBLE"] += 1
                elif eid in ("CYP3A4_INHIBITION", "CYP2D6_INHIBITION", "CYP2C9_INHIBITION", "CYP2C19_INHIBITION", "CYP1A2_INHIBITION", "PGP_INHIBITION", "HERG_LIABILITY"):
                    if r.raw_unit in ("%", "% inhibition") and 0.0 <= float(r.raw_value) <= 100.0:
                        s_n += 1
                        bin_label = 1 if float(r.raw_value) >= 50.0 else 0
                        indep_pairs[cname] = (smi, bin_label, "binary_class")
                    elif r.raw_unit in ("nM", "uM", "µM"):
                        reasons["QUANTITATIVE_EVIDENCE_NOT_CLASSIFICATION_PAIRABLE"] += 1
                    else:
                        reasons["UNIT_CONTEXT_INCOMPATIBLE"] += 1
                elif eid in ("HLM_CLINT", "RLM_CLINT", "MLM_CLINT"):
                    if r.normalized_unit in ("log10(mL/min/kg)", "mL/min/kg") and r.raw_unit not in ("hours", "%", "L", "L/h"):
                        s_n += 1
                        indep_pairs[cname] = (smi, float(r.normalized_value), r.normalized_unit)
                    else:
                        reasons["IN_VIVO_CL_OR_METABOLIC_FRACTION_NOT_IN_VITRO_MICROSOMAL_CLINT"] += 1
                else:
                    reasons["NO_ACTIVE_ADAPTER_CONTRACT"] += 1

            for k, v in reasons.items():
                funnel_stats["global_drop_reasons"][k] += v

            funnel_stats["endpoints"][eid] = {
                "qualified_n": q_n,
                "prediction_available_n": p_n,
                "semantic_compatible_n": s_n,
                "independent_pair_n": len(indep_pairs),
                "drop_reasons": dict(reasons),
            }

        return funnel_stats
    finally:
        db.close()


def run_endpoint_validation() -> List[Dict[str, Any]]:
    """
    Executes the recovered, scale-aligned validation across all registered endpoints.
    Eliminates unit/scale mismatches (e.g. logS vs uM) and arbitrary continuous IC50 binarization.
    """
    db = SessionLocal()
    try:
        strategies = get_all_strategies()
        records = db.query(ExternalExperimentalEvidence).filter(
            ExternalExperimentalEvidence.evidence_state.in_(["AUTO_QUALIFIED_EXTERNAL", "RELATED_EXTERNAL"])
        ).all()

        # Map records to canonical endpoints with scale reconciliation
        endpoint_data = defaultdict(dict)
        for r in records:
            cv = db.query(CompoundVersion).filter(CompoundVersion.id == r.compound_version_id).first()
            if not cv:
                continue
            comp = db.query(Compound).filter(Compound.id == cv.compound_row_id).first()
            if not comp:
                continue

            eid = r.canonical_endpoint_id
            cname = comp.name
            smi = cv.canonical_smiles

            # Scale and unit alignment
            if eid == "SOLUBILITY_GENERIC":
                if r.normalized_unit in ("log10(mol/L)", "logS"):
                    val_logs = float(r.normalized_value)
                    endpoint_data[eid][cname] = {"smiles": smi, "exp_val": val_logs, "exp_unit": "log10(mol/L)"}
                elif r.raw_unit in ("uM", "µM", "uM ", "µM "):
                    val_logs = math.log10(float(r.raw_value) * 1e-6)
                    endpoint_data[eid][cname] = {"smiles": smi, "exp_val": val_logs, "exp_unit": "log10(mol/L)"}
            elif eid == "HUMAN_PPB":
                if r.normalized_unit in ("% bound", "%") and 0.0 <= float(r.normalized_value) <= 100.0 and r.raw_unit not in ("nM", "hours", "uM", "µM"):
                    endpoint_data[eid][cname] = {"smiles": smi, "exp_val": float(r.normalized_value), "exp_unit": "% bound"}
            elif eid in ("CYP3A4_INHIBITION", "CYP2D6_INHIBITION", "CYP2C9_INHIBITION", "CYP2C19_INHIBITION", "CYP1A2_INHIBITION", "PGP_INHIBITION", "HERG_LIABILITY"):
                if r.raw_unit in ("%", "% inhibition") and 0.0 <= float(r.raw_value) <= 100.0:
                    bin_label = 1 if float(r.raw_value) >= 50.0 else 0
                    endpoint_data[eid][cname] = {"smiles": smi, "exp_val": bin_label, "exp_unit": "binary_class"}
            elif eid in ("HLM_CLINT", "RLM_CLINT", "MLM_CLINT"):
                if r.normalized_unit in ("log10(mL/min/kg)", "mL/min/kg") and r.raw_unit not in ("hours", "%", "L", "L/h"):
                    endpoint_data[eid][cname] = {"smiles": smi, "exp_val": float(r.normalized_value), "exp_unit": r.normalized_unit}

        header_fmt = "{:<24} | {:<4} | {:<32} | {:<14} | {:<45} | {:<14} | {}"
        print(header_fmt.format("Endpoint", "N", "Primary Model", "Primary Err", "Alternative Errors", "Consensus Err", "Decision"))
        print("-" * 160)

        report_rows = []

        for ep_name, canon_id in EP_MAP.items():
            comp_pairs = list(endpoint_data.get(canon_id, {}).values())
            policy = strategies.get(ep_name)
            primary_id = policy.primary_model_ids[0] if (policy and policy.primary_model_ids) else "N/A"
            adapters = get_v2_adapters_for_endpoint(ep_name)
            contract = get_endpoint_contract(canon_id) or get_endpoint_contract(ep_name)

            cnt = len(comp_pairs)
            if cnt == 0:
                row_str = header_fmt.format(ep_name, cnt, primary_id, "No Data", "None", "N/A", "DEFERRED (N=0)")
                print(row_str)
                report_rows.append({
                    "endpoint_name": ep_name,
                    "canonical_endpoint_id": canon_id,
                    "qualified_n": len(records),
                    "pairable_n_before": 0,
                    "pairable_n_after": 0,
                    "independent_n": 0,
                    "primary_model_id": primary_id,
                    "primary_error": "No Data",
                    "alternative_errors": {},
                    "consensus_error": "N/A",
                    "decision": "DEFERRED_INSUFFICIENT_N (N=0)",
                    "reason_remaining": "No experimental validation pairs in current dataset",
                })
                continue

            model_preds = {a.model_id: [] for a in adapters}
            exp_vals = [it["exp_val"] for it in comp_pairs]

            # Separate ML models from rule/derived
            ml_model_ids = [
                a.model_id for a in adapters
                if not any(k in a.model_family.lower() for k in ("rule", "physicochemical_mechanistic", "derived"))
            ]

            for it in comp_pairs:
                smi = it["smiles"]
                for a in adapters:
                    try:
                        res = a.execute(smi, contract)
                        val = res.value if (res and res.execution_status == "SUCCESS") else None
                        model_preds[a.model_id].append(val)
                    except Exception:
                        model_preds[a.model_id].append(None)

            # Consensus: Mean of valid ML models only
            consensus_preds = []
            for i in range(len(comp_pairs)):
                ml_vals = [model_preds[m_id][i] for m_id in ml_model_ids if model_preds[m_id][i] is not None]
                if ml_vals:
                    consensus_preds.append(float(np.mean(ml_vals)))
                else:
                    consensus_preds.append(None)

            # Compute error metrics
            errors = {}
            for mid, preds in model_preds.items():
                valid_pairs = [(p, e) for p, e in zip(preds, exp_vals) if p is not None]
                if valid_pairs:
                    p_arr, e_arr = zip(*valid_pairs)
                    if contract and contract.output_type == OutputType.BINARY_CLASSIFICATION:
                        b_pred = [1 if p >= 0.5 else 0 for p in p_arr]
                        acc = np.mean([1 if bp == be else 0 for bp, be in zip(b_pred, e_arr)])
                        errors[mid] = f"Acc: {acc*100:.1f}%"
                    else:
                        mae = np.mean(np.abs(np.array(p_arr) - np.array(e_arr)))
                        errors[mid] = f"MAE: {mae:.2f}"
                else:
                    errors[mid] = "No valid preds"

            # Consensus error
            valid_cons_pairs = [(p, e) for p, e in zip(consensus_preds, exp_vals) if p is not None]
            if valid_cons_pairs:
                p_arr, e_arr = zip(*valid_cons_pairs)
                if contract and contract.output_type == OutputType.BINARY_CLASSIFICATION:
                    b_pred = [1 if p >= 0.5 else 0 for p in p_arr]
                    acc = np.mean([1 if bp == be else 0 for bp, be in zip(b_pred, e_arr)])
                    cons_err = f"Acc: {acc*100:.1f}%"
                else:
                    mae = np.mean(np.abs(np.array(p_arr) - np.array(e_arr)))
                    cons_err = f"MAE: {mae:.2f}"
            else:
                cons_err = "N/A"

            p_err = errors.get(primary_id, "N/A")
            alt_strs = [f"{mid}: {errors.get(mid, 'N/A')}" for mid in model_preds if mid != primary_id]
            alt_desc = "; ".join(alt_strs) if alt_strs else "None"

            # Decision
            if cnt < 3:
                decision = f"DEFERRED_INSUFFICIENT_N (N={cnt})"
            else:
                decision = "RETAIN_CURRENT_PRIMARY"

            row_str = header_fmt.format(ep_name, cnt, primary_id, p_err, alt_desc[:45], cons_err, decision)
            print(row_str)

            report_rows.append({
                "endpoint_name": ep_name,
                "canonical_endpoint_id": canon_id,
                "independent_n": cnt,
                "primary_model_id": primary_id,
                "primary_error": p_err,
                "alternative_errors": {mid: errors.get(mid, "N/A") for mid in model_preds if mid != primary_id},
                "consensus_error": cons_err,
                "decision": decision,
                "reason_remaining": "Scale-aligned validation complete" if cnt >= 3 else f"Insufficient independent compounds (N={cnt})",
            })

        return report_rows
    finally:
        db.close()


if __name__ == "__main__":
    funnel = audit_evidence_funnel()
    print("\nFunnel Global Drop Reasons:")
    for r, c in funnel["global_drop_reasons"].most_common():
        print(f"  * {r:<55}: {c:4d}")
    print("\n" + "=" * 160 + "\n")
    run_endpoint_validation()
