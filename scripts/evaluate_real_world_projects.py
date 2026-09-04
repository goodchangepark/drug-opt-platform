"""
Evaluates Prediction Engine v3.3 vs v3.3.1 on Real-World Internal Drug-OPT Projects:
- Project 1: GLP-1 (small molecule) [4 compounds]
- Project 3: EGFR [7 compounds]
- Project 5: AMYR [4 compounds]
Total = 15 active therapeutic compounds.
Exports: validation/real_world_project_benchmark_v3_3_1.json
"""
from __future__ import annotations

import json
from pathlib import Path
import sqlite3
from typing import Any, Dict, List

from backend.candidate_model_registry import register_candidate_models_to_multimodel
from backend.multimodel import get_model_adapter
from backend.endpoint_contracts import get_endpoint_contract

register_candidate_models_to_multimodel()

con = sqlite3.connect("drug_opt.db")
cur = con.cursor()

PROJECTS = [
    (1, "GLP-1 (small molecule)", 4),
    (3, "EGFR", 7),
    (5, "AMYR (small molecules)", 4),
]

# Stacking weights obtained from DrugBank 150 benchmark:
STACKING_WEIGHTS = {
    "SOLUBILITY": [
        ("admetica_solubility", 0.191),
        ("esol_delaney_v1", 0.723),
        ("rdkit_gbr_solubility_v1", 0.086),
    ],
    "CACO2": [
        ("admetica_caco2", 0.704),
        ("physchem_caco2_v1", 0.296),
    ],
    "PPB": [
        ("admetica_ppbr", 0.728),
        ("physchem_human_ppb_v1", 0.223),
        ("descriptor_gbr_ppb_v1", 0.048),
    ],
    "HLM": [
        ("drugopt_hlm_chemical_space_v1", 1.000),
    ],
    "CYP3A4_PIC50": [
        ("openadmet_chemeleon_cyp3a4_pic50", 0.118),
        ("drugopt_calibrated_cyp3a4_pic50", 0.882),
    ],
    "CYP2D6_PIC50": [
        ("openadmet_chemeleon_cyp2d6_pic50", 0.385),
        ("drugopt_calibrated_cyp2d6_pic50", 0.615),
    ],
    "CYP1A2_PIC50": [
        ("drugopt_calibrated_cyp1a2_pic50", 1.000),
    ],
    "CYP2C9_PIC50": [
        ("drugopt_calibrated_cyp2c9_pic50", 1.000),
    ],
    "HERG_PIC50": [
        ("physchem_gbr_herg_pic50_v1", 1.000),
    ]
}

CONTRACT_MAP = {
    "SOLUBILITY": "Solubility",
    "CACO2": "Permeability",
    "PPB": "Plasma protein binding",
    "HLM": "HLM intrinsic clearance",
    "CYP3A4_PIC50": "CYP3A4 inhibitor",
    "CYP2D6_PIC50": "CYP2D6 inhibitor",
    "CYP1A2_PIC50": "CYP1A2 inhibitor",
    "CYP2C9_PIC50": "CYP2C9 inhibitor",
    "HERG_PIC50": "hERG liability",
}

def predict_v3_3_1_ensemble(smi: str, ep_key: str) -> Dict[str, Any]:
    weights = STACKING_WEIGHTS.get(ep_key, [])
    if not weights:
        return {"value": None, "status": "MODEL_UNAVAILABLE"}
    
    contract_name = CONTRACT_MAP.get(ep_key)
    contract = get_endpoint_contract(contract_name) if contract_name else None
    
    vals = []
    w_sum = 0.0
    for m_id, w in weights:
        ad = get_model_adapter(m_id)
        if ad is None:
            continue
        res = ad.execute(smi, contract)
        if res.execution_status.value == "SUCCESS" and res.value is not None:
            vals.append(float(res.value) * w)
            w_sum += w
    if w_sum > 0:
        return {"value": round(sum(vals) / w_sum, 4), "status": "SUCCESS"}
    return {"value": None, "status": "RUNTIME_ERROR"}

compounds_data = []

for pid, pname, target_n in PROJECTS:
    cur.execute(
        "SELECT c.id, c.name, cv.canonical_smiles FROM compounds c "
        "JOIN compound_versions cv ON c.id = cv.compound_row_id "
        "WHERE c.project_id = ? AND cv.version_number = c.current_version",
        (pid,)
    )
    rows = cur.fetchall()
    print(f"\nEvaluating Project {pid}: {pname} (N={len(rows)})")
    
    for cid, cname, smi in rows:
        comp_res = {
            "compound_id": cid,
            "name": cname,
            "project_id": pid,
            "project_name": pname,
            "smiles": smi,
            "predictions_v3_3_1": {}
        }
        for ep_key in STACKING_WEIGHTS.keys():
            res = predict_v3_3_1_ensemble(smi, ep_key)
            comp_res["predictions_v3_3_1"][ep_key] = res
        compounds_data.append(comp_res)

con.close()

out_file = Path("validation/real_world_project_benchmark_v3_3_1.json")
with open(out_file, "w", encoding="utf-8") as f:
    json.dump({
        "total_compounds": len(compounds_data),
        "compounds": compounds_data
    }, f, indent=2)

print(f"\nSaved Real-World Benchmark for {len(compounds_data)} compounds to {out_file}")
