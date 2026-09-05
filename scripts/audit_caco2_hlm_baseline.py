import json
import sqlite3
import numpy as np
from backend.prediction_maturity import ENDPOINT_MATURITY_REGISTRY

# 1. 50-Endpoint Taxonomy Inventory Audit
total_endpoints = len(ENDPOINT_MATURITY_REGISTRY)
by_type = {
    "DETERMINISTIC": 0,
    "RULE_ESTIMATE": 0,
    "QUANTITATIVE_MODEL": 0,
    "CLASSIFICATION_MODEL": 0,
    "MODEL_UNAVAILABLE": 0
}

for ep_id, ep in ENDPOINT_MATURITY_REGISTRY.items():
    route = ep.get("model_route", "")
    is_mech = ep.get("is_mechanistic", False)
    is_unavail = ep.get("is_unavailable", False)
    lvl = ep.get("maturity_level", 1)
    
    if is_unavail or route == "MODEL_UNAVAILABLE":
        by_type["MODEL_UNAVAILABLE"] += 1
    elif route == "CLASSIFICATION_ONLY" or "CLASSIFICATION" in route:
        by_type["CLASSIFICATION_MODEL"] += 1
    elif is_mech:
        by_type["RULE_ESTIMATE"] += 1
    elif ep_id in ["MW", "CLOGP", "TPSA", "HBD", "HBA", "ROTB", "FSP3", "QED", "FORMAL_CHARGE", "NUM_RINGS", "HEAVY_ATOM_COUNT", "LOGS_ESOL", "PKA_ACID", "PKA_BASE"]:
        by_type["DETERMINISTIC"] += 1
    else:
        by_type["QUANTITATIVE_MODEL"] += 1

print("=== 50-ENDPOINT INVENTORY AUDIT ===")
print(f"Total Canonical Endpoints: {total_endpoints}")
for k, v in by_type.items():
    print(f"  - {k}: {v}")
print(f"Sum of types: {sum(by_type.values())}")
assert sum(by_type.values()) == 50, f"Expected 50 endpoints, got {sum(by_type.values())}"

# 2. Caco-2 Audit
# Load Caco-2 validation cohort from existing benchmark/dataset if available
print("\n=== CACO-2 AUDIT (Log10 cm/s normalized scale) ===")
# Check Caco-2 cohort records
conn = sqlite3.connect("drug_opt.db")
c = conn.cursor()

# Check evidence observations for Caco-2
c.execute("""
    SELECT count(*) FROM external_experimental_evidence 
    WHERE canonical_endpoint_id LIKE '%CACO2%' OR raw_endpoint_name LIKE '%Caco%'
""")
caco2_ev_cnt = c.fetchone()[0]
print(f"Total Caco-2 evidence observations in DB: {caco2_ev_cnt}")

# Check paired predictions on Caco-2
c.execute("""
    SELECT count(*) FROM prediction_endpoint_snapshots
    WHERE endpoint_name LIKE '%CACO2%' OR endpoint_name LIKE '%Caco%'
""")
caco2_snap_cnt = c.fetchone()[0]
print(f"Total Caco-2 snapshots in DB: {caco2_snap_cnt}")

# Check HLM audit
print("\n=== HLM AUDIT (species=human, matrix=liver microsome) ===")
c.execute("""
    SELECT count(*) FROM external_experimental_evidence 
    WHERE canonical_endpoint_id LIKE '%HLM%' OR raw_endpoint_name LIKE '%HLM%' OR raw_endpoint_name LIKE '%microsom%'
""")
hlm_ev_cnt = c.fetchone()[0]
print(f"Total HLM evidence observations in DB: {hlm_ev_cnt}")

conn.close()
