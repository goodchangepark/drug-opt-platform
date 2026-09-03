"""
Validation and Multi-Model Comparison across All Endpoints (Drug-OPT v5.3).
"""
import math
import numpy as np
from typing import Dict, List, Any, Optional

from backend.database import SessionLocal
from backend.models import Compound
from backend.endpoint_contracts import get_endpoint_contract, OutputType
from backend.multimodel import (
    get_v2_adapters_for_endpoint,
    get_model_adapter,
)
from backend.endpoint_strategy_registry import get_all_strategies
from backend.endpoint_comparison import build_endpoint_comparison

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


def run_endpoint_validation():
    db = SessionLocal()
    compounds = db.query(Compound).all()
    valid_comps = [(c, c.versions[-1]) for c in compounds if c.versions and c.versions[-1].canonical_smiles]

    exp_by_ep = {}
    for c, v in valid_comps:
        res = build_endpoint_comparison(db, v.id)
        for r in res.get("scientific_rows", []):
            eid = r.get("canonical_endpoint")
            disp = r.get("primary_experimental_display", {})
            val = disp.get("value")
            if val is not None and isinstance(val, (int, float)):
                exp_by_ep.setdefault(eid, []).append({
                    "compound_id": c.id,
                    "name": c.name,
                    "smiles": v.canonical_smiles,
                    "exp_val": float(val),
                    "exp_unit": disp.get("unit"),
                })

    strategies = get_all_strategies()

    header_fmt = "{:<24} | {:<4} | {:<32} | {:<14} | {:<45} | {:<14} | {}"
    print(header_fmt.format("Endpoint", "N", "Primary Model", "Primary Err", "Alternative Errors", "Consensus Err", "Decision"))
    print("-" * 160)

    report_rows = []

    for ep_name, canon_id in EP_MAP.items():
        items = exp_by_ep.get(canon_id, [])
        policy = strategies.get(ep_name)
        primary_id = policy.primary_model_ids[0] if (policy and policy.primary_model_ids) else "N/A"
        adapters = get_v2_adapters_for_endpoint(ep_name)
        contract = get_endpoint_contract(canon_id) or get_endpoint_contract(ep_name)

        cnt = len(items)
        if cnt == 0:
            row_str = header_fmt.format(ep_name, cnt, primary_id, "No Data", "None", "N/A", "DEFERRED (N=0)")
            print(row_str)
            report_rows.append({
                "endpoint_name": ep_name,
                "canonical_endpoint_id": canon_id,
                "independent_n": 0,
                "primary_model_id": primary_id,
                "primary_error": "No Data",
                "alternative_errors": {},
                "consensus_error": "N/A",
                "decision": "DEFERRED_INSUFFICIENT_N (N=0)",
            })
            continue

        model_preds = {a.model_id: [] for a in adapters}
        exp_vals = [it["exp_val"] for it in items]

        # Filter ML models vs Rule/derived for consensus
        ml_model_ids = [
            a.model_id for a in adapters
            if not any(k in a.model_family.lower() for k in ("rule", "physicochemical_mechanistic", "derived"))
        ]

        for it in items:
            smi = it["smiles"]
            for a in adapters:
                try:
                    res = a.execute(smi, contract)
                    val = res.value if (res and res.execution_status == "SUCCESS") else None
                    model_preds[a.model_id].append(val)
                except Exception:
                    model_preds[a.model_id].append(None)

        # Compute Consensus predictions (mean of valid ML models only)
        consensus_preds = []
        for i in range(len(items)):
            ml_vals = [model_preds[m_id][i] for m_id in ml_model_ids if model_preds[m_id][i] is not None]
            if ml_vals:
                consensus_preds.append(float(np.mean(ml_vals)))
            else:
                consensus_preds.append(None)

        # Compute error metrics per model and consensus
        errors = {}
        for mid, preds in model_preds.items():
            valid_pairs = [(p, e) for p, e in zip(preds, exp_vals) if p is not None]
            if valid_pairs:
                p_arr, e_arr = zip(*valid_pairs)
                if contract and contract.output_type == OutputType.BINARY_CLASSIFICATION:
                    b_exp = [1 if e > 50 or e == 1.0 else 0 for e in e_arr]
                    b_pred = [1 if p >= 0.5 else 0 for p in p_arr]
                    acc = np.mean([1 if bp == be else 0 for bp, be in zip(b_pred, b_exp)])
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
                b_exp = [1 if e > 50 or e == 1.0 else 0 for e in e_arr]
                b_pred = [1 if p >= 0.5 else 0 for p in p_arr]
                acc = np.mean([1 if bp == be else 0 for bp, be in zip(b_pred, b_exp)])
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
        })

    db.close()
    return report_rows


if __name__ == "__main__":
    run_endpoint_validation()
