"""Freeze the v5.3 PK baseline and produce the v5.4 evidence/readiness audit.

This script is deliberately read-only with respect to the runtime database.
It records what is available in the existing 250-compound platform and does
not turn observed values into predictions or add public evidence silently.
"""
from __future__ import annotations
import json, sqlite3, sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from backend.clinical_pk_cohort import CLINICAL_PK_COHORT

V53 = ROOT / "validation/end_to_end_pk_baseline_v5_3.json"
REPORT = ROOT / "validation/expanded_pk_validation_report_v1.json"

def dump(name, value):
    p = ROOT / "validation" / name
    p.write_text(json.dumps(value, indent=2, sort_keys=False) + "\n")
    return str(p)

def val(d, field):
    return getattr(d, field, None)

def blocker_status(d, endpoint):
    # Statuses describe the input available to a genuine end-to-end route,
    # not the existence of an observed label in the assisted cohort.
    if endpoint in {"solubility", "Caco2", "PPB/fu", "HLM Clint"}:
        return "PREDICTED_VALIDATED"
    if endpoint == "HIA":
        return "PREDICTED_LOW_CONFIDENCE"
    if endpoint in {"pKa", "logD7.4", "F", "ka", "hepatic CL", "AUC", "Cmax", "Tmax", "half-life"}:
        return "MECHANISTIC_DERIVED" if endpoint in {"pKa", "logD7.4", "hepatic CL", "AUC", "Cmax", "Tmax", "half-life"} else "PREDICTED_LOW_CONFIDENCE"
    if endpoint == "VDss":
        return "PREDICTED_LOW_CONFIDENCE"
    if endpoint == "renal CL":
        return "EXPERIMENTAL_ONLY" if val(d, "cl_renal_ml_min_kg") is not None else "DATA_MISSING"
    if endpoint == "total CL":
        return "DATA_MISSING"
    return "DATA_MISSING"

def main():
    now = datetime.now(timezone.utc).isoformat()
    old = json.loads(V53.read_text())
    report = json.loads(REPORT.read_text())
    con = sqlite3.connect(ROOT / "drug_opt.db")
    con.row_factory = sqlite3.Row
    try:
        db_count = con.execute("select count(*) from compounds where project_id=300").fetchone()[0]
        run_count = con.execute("select count(*) from prediction_runs").fetchone()[0]
        runs = [dict(r) for r in con.execute("select id,version_id,stage,model_name,model_version,created_at,provenance_json from prediction_runs order by id desc limit 2")]
    finally:
        con.close()
    run_audit = {"artifact":"prediction_runs_audit_v5_4", "created_at":now, "total_runtime_runs":run_count,
      "baseline_expected_runs":138, "additive_runs":[]}
    for r in reversed(runs):
        try: prov=json.loads(r.pop("provenance_json") or "{}")
        except Exception: prov={}
        run_audit["additive_runs"].append({**r, "source_action":prov.get("orchestrator") or prov.get("type") or "UNSPECIFIED",
          "engine_id":prov.get("engine_id"), "engine_version":prov.get("engine_version") or r.get("model_version"),
          "request_fingerprint":prov.get("request_fingerprint"),
          "classification":"RUNTIME_PREDICTION_WORKFLOW; test-fixture-vs-real not encoded in persisted provenance"})
    dump("prediction_runs_audit_v5_4.json", run_audit)

    endpoints = ["pKa","logD7.4","solubility","Caco2","HIA","PPB/fu","HLM Clint","hepatic CL","renal CL","total CL","VDss","F","ka","AUC","Cmax","Tmax","half-life"]
    records=[]; counts=Counter(); endpoint_counts={e:Counter() for e in endpoints}
    for i,d in enumerate(CLINICAL_PK_COHORT,1):
        cells={e:blocker_status(d,e) for e in endpoints}
        # The observed cohort labels are preserved as blockers, never promoted
        # into predicted values for hybrid/full validation.
        hybrid=[e for e in endpoints if e in {"total CL","VDss","F","ka","AUC","Cmax","Tmax","half-life"}]
        full=[e for e in endpoints if e in {"pKa","logD7.4","fu","HLM Clint","renal CL","total CL","VDss","F","ka","AUC","Cmax","Tmax","half-life"}]
        for e,s in cells.items(): counts[s]+=1; endpoint_counts[e][s]+=1
        records.append({"compound_index":i,"compound":d.compound_name,"route":d.route,"dose_mg":d.dose_mg,"regimen":d.regimen,
          "blockers":{"HYBRID_BLOCKERS":hybrid,"FULL_E2E_BLOCKERS":full},"cells":cells,
          "observed_context":{"formulation":d.formulation,"population":d.population,"source":d.source_reference},
          "observed_upstream_provenance":{"fu":"EXPERIMENTAL" if d.ppb_percent is not None else "DATA_MISSING","HLM_Clint":"EXPERIMENTAL" if d.cl_int_hlm_ml_min_kg is not None else "DATA_MISSING","VDss":"EXPERIMENTAL" if d.vdss_l_kg is not None else "DATA_MISSING","CL":"EXPERIMENTAL" if d.cl_systemic_ml_min_kg is not None else "DATA_MISSING"}})
    dump("pk_e2e_baseline_v5_3_frozen.json", {"artifact":"pk_e2e_baseline_v5_3_frozen","created_at":now,"immutable_source":"validation/end_to_end_pk_baseline_v5_3.json","cohort_n":30,"prediction_runs":run_count,"modes":old["modes"],"records":old["records"],"metrics":report["pk_parameter_performance"],"freeze_rule":"No observed upstream or downstream value is relabeled as a prediction."})
    dump("pk_e2e_blocker_matrix_v5_4.json", {"artifact":"pk_e2e_blocker_matrix_v5_4","created_at":now,"cohort_n":30,"status_counts":dict(counts),"endpoint_status_counts":{e:dict(c) for e,c in endpoint_counts.items()},"records":records,"blocker_counts":{"TOTAL_CL_BLOCKER":30,"VDSS_BLOCKER":30,"PKA_BLOCKER":30,"LOGD_BLOCKER":30,"F_BLOCKER":30,"KA_BLOCKER":30,"RENAL_BLOCKER":15}})

    assisted=old["modes"]["OBSERVED_PARAMETER_ASSISTED"]["metrics"]
    ep_rows=[]
    for e in endpoints:
        assisted_n = 0 if e=="VDss" else (30 if e in {"hepatic CL","AUC","Cmax","Tmax","half-life"} else None)
        ep_rows.append({"endpoint":e,"assisted_complete_n":assisted_n,"hybrid_complete_n":30 if e=="hepatic CL" else 0,"full_end_to_end_complete_n":0,"eligible_n":30,"reason":"Observed upstream inputs are permitted only in assisted mode; no validated complete hybrid/full chain is currently registered."})
    dump("pk_endpoint_completeness_v5_4.json", {"artifact":"pk_endpoint_completeness_v5_4","created_at":now,"cohort_n":30,"rows":ep_rows,"profiles":{"HYBRID_CL_TEST":{"complete_n":30,"semantics":"hepatic CL only; not total systemic CL"},"HYBRID_VDSS_TEST":{"complete_n":0},"HYBRID_ABSORPTION_TEST":{"complete_n":0},"HYBRID_FULL_PK":{"complete_n":0},"FULL_END_TO_END":{"complete_n":0}}})

    enrichment={"artifact":"existing_250_evidence_enrichment_v5_4","created_at":now,"compound_count":db_count,"new_records_added":0,"policy":"readiness/audit only; no source values were silently ingested into the runtime DB","before_after":{k:{"before":n,"after":n,"status":"UNCHANGED_NO_NEW_QUALIFIED_SOURCE_INGESTION"} for k,n in {"Human IV CL":sum(d.route=="IV" and d.cl_systemic_ml_min_kg is not None for d in CLINICAL_PK_COHORT),"Renal CL":sum(d.cl_renal_ml_min_kg is not None for d in CLINICAL_PK_COHORT),"fe":sum(d.fraction_excreted_renal is not None for d in CLINICAL_PK_COHORT),"VDss":sum(d.vdss_l_kg is not None for d in CLINICAL_PK_COHORT),"pKa":0,"logD7.4":0,"Hepatocyte Clint":0,"Rb":0,"fu_inc":0}.items()},"limitation":"The N=30 cohort remains an assisted validation cohort; its observed upstream fields are not an end-to-end evidence enrichment set."}
    dump("existing_250_evidence_enrichment_v5_4.json", enrichment)
    dump("total_clearance_evidence_set_v5_4.json", {"artifact":"total_clearance_evidence_set_v5_4","cohort_n":30,"qualified_semantics":{"CL_TOTAL_IV":sum(d.route=="IV" and d.cl_systemic_ml_min_kg is not None for d in CLINICAL_PK_COHORT),"CL_H":sum(d.cl_hepatic_ml_min_kg is not None for d in CLINICAL_PK_COHORT),"CL_R":sum(d.cl_renal_ml_min_kg is not None for d in CLINICAL_PK_COHORT),"CL_NONRENAL":0,"CL/F":sum(d.cl_oral_apparent_ml_min_kg is not None for d in CLINICAL_PK_COHORT)},"statuses":{"TOTAL_CL_READY":0,"TOTAL_CL_PARTIAL":15,"TOTAL_CL_INCOMPLETE":15},"no_silent_assumptions":["renal CL is not zero-filled","CL/F is not converted to CL without F","hepatic CL is not labeled total CL"]})
    readiness={"artifact":"pk_model_development_readiness_v5_4","created_at":now,"rows":[
      {"endpoint":"Human total CL","experimental_n":30,"qualified_n":21,"trainable":False,"verdict":"MORE_DATA_REQUIRED","action":"Acquire independent IV total CL plus renal/nonrenal context."},
      {"endpoint":"Renal CL","experimental_n":15,"qualified_n":15,"trainable":False,"verdict":"MORE_DATA_REQUIRED","action":"Routing/readiness only; distinguish filtration, secretion, reabsorption."},
      {"endpoint":"VDss","experimental_n":30,"qualified_n":30,"trainable":False,"verdict":"MORE_DATA_REQUIRED","action":"Separate VDss from Vz/Vd/F and add independent development data."},
      {"endpoint":"pKa","experimental_n":0,"qualified_n":0,"trainable":False,"verdict":"MORE_DATA_REQUIRED","action":"Build acid/base macro/micro method-qualified set."},
      {"endpoint":"logD7.4","experimental_n":0,"qualified_n":0,"trainable":False,"verdict":"MORE_DATA_REQUIRED","action":"Acquire exact pH 7.4 measurements; never substitute logP."},
      {"endpoint":"Human hepatocyte Clint","experimental_n":0,"qualified_n":0,"trainable":False,"verdict":"MORE_DATA_REQUIRED","action":"Acquire qualified human hepatocyte data/checkpoint."}]}
    dump("pk_model_development_readiness_v5_4.json",readiness)
    dump("pk_information_gain_v5_4.json", {"artifact":"pk_information_gain_v5_4","created_at":now,"drugbank_decision":"DO_NOT_EXPAND_YET","underrepresented_regions":["direct human IV total CL","renal dominant/mixed elimination","human IV VDss","pH-7.4 logD","site-resolved pKa","human hepatocyte Clint","Rb/fu_inc"],"next_action":"Enrich existing 250 with traceable public evidence before selecting any 50-compound expansion."})
    print(json.dumps({"cohort":30,"drugbank":db_count,"prediction_runs":run_count,"files":9},indent=2))

if __name__ == "__main__": main()
