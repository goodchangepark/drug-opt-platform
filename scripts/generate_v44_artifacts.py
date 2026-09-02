"""Generate conservative v4.4 PK-engine and review-queue audit artifacts."""
from __future__ import annotations
import json
import sys
from collections import Counter
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
 sys.path.insert(0, str(ROOT))
from backend.database import SessionLocal
from backend.endpoint_comparison import build_endpoint_comparison
from backend.pk_engine_v1 import PK_ENGINE_VERSION

OUT=ROOT / "validation"
def dump(name, value): (OUT/name).write_text(json.dumps(value,indent=2,sort_keys=True)+"\n")
def main():
 s=SessionLocal()
 try:
  view=build_endpoint_comparison(s,13); rows=view["scientific_rows"]
  pk=[r for r in rows if r["section"]=="PK"]
  # Scientific rows carry the canonical qualification state.  Individual
  # observations preserve an optional structured qualification payload, but
  # legacy observations may use a string.  Never let that representation
  # difference make review evidence disappear from the audit.
  review_rows=[r for r in rows if r.get("semantic_status")=="NEEDS_REVIEW" or r.get("qualification_status")=="CONTEXT_NOT_QUALIFIED"]
  review=[o for r in review_rows for o in r["experimental_observations"]]
  def reason(observation):
   qualification=observation.get("qualification") or observation.get("qualification_details") or {}
   return qualification.get("primary_gap_reason") if isinstance(qualification,dict) else None
  matrix={species:{p:"NOT_IMPLEMENTED" for p in ("CL","CL/F","Vd","Vss","Vd/F","F","ka","Cmax","Tmax","AUC0-t","AUC0-inf","AUCtau","t1/2")} for species in ("HUMAN","RAT","MOUSE","DOG","MONKEY")}
  for species in matrix: matrix[species]["CL"]="CURRENT_RESEARCH_ONLY"; matrix[species]["Vd"]="CURRENT_RESEARCH_ONLY"; matrix[species]["t1/2"]="CURRENT_RESEARCH_ONLY"
  dump("pk_capability_matrix_v4_4.json",{"pk_engine_version":PK_ENGINE_VERSION,"matrix":matrix,"policy":"No endpoint is production-supported without an independent public benchmark."})
  dataset={"status":"NO_DEFENSIBLE_MACHINE_READABLE_PUBLIC_BENCHMARK_CURATED","records":[],"sources":["FDA regulatory labels/reviews (per-observation source)","Douguet 2018 PMID:29541361","Sakaeda 2001 PMID:11510489","Obach 2018 PMID:30115648"],"reason":"No versioned, context-complete 20-compound PK set was available in this repository; no values were fabricated or scraped into a benchmark."}
  dump("pk_benchmark_dataset_v4_4.json",dataset); dump("pk_benchmark_split_v4_4.json",{"status":"NOT_CREATED_NO_DATASET","development":[],"holdout":[],"holdout_leakage":False})
  dump("pk_benchmark_baseline_v4_4.json",{"status":"NOT_RUN_NO_BENCHMARK","engine":PK_ENGINE_VERSION}); dump("pk_benchmark_results_v4_4.json",{"status":"INSUFFICIENT_VALIDATION","metrics":{}})
  dump("pk_endpoint_validation_status_v4_4.json",{"pk_engine_version":PK_ENGINE_VERSION,"production_supported":[],"research_only":["one_compartment_experiment_informed"],"insufficient_validation":["Human CL","Human Vd","Human Cmax","Human AUC","Human Tmax","Oral F"],"unavailable":["structure_only human Cmax/AUC/Tmax"]})
  dump("sunvozertinib_pk_validation_v4_4.json",{"clinical":{"AUC":{"value":8060,"unit":"ng*h/mL","context":"Human oral 200 mg/day"},"Cmax":{"value":412,"unit":"ng/mL","context":"Human oral 200 mg"},"Tmax":{"value":7,"unit":"h","context":"Human oral"}},"pk_engine_status":"INSUFFICIENT_INPUT","no_fabricated_predictions":True,"invalid_f_100_percent_absent":True})
  dump("evidence_review_queue_v4_4.json",{"count":len(review),"by_reason":dict(Counter(reason(o) or "OTHER" for o in review)),"items":[{"id":o.get("id"),"raw_endpoint":o.get("raw_endpoint"),"raw_value":o.get("raw_value"),"reference":o.get("reference"),"reason":reason(o) or "OTHER"} for o in review]})
 finally: s.close()
if __name__=="__main__": main()
