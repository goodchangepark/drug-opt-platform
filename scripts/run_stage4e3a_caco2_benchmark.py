#!/usr/bin/env python3
"""Fixed-model ExpansionRx Caco-2 benchmark; no fitting or DB writes."""
from __future__ import annotations
import argparse, csv, hashlib, json, math, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr, pearsonr
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from backend.standardizer import standardize_molecule, STANDARDIZER_VERSION
from backend.admet_predictor import predict_batch_values, applicability_domain
from backend.endpoint_contracts import ENDPOINT_CONTRACTS
from backend.multimodel import get_model_adapter

RAW=Path('/tmp/stage4e2r-expansionrx/expansion_data_raw.csv'); OUT=ROOT/'validation'; SEED=20260829
def metrics(y,p):
 e=p-y; ae=np.abs(e)
 return {'N':len(y),'MAE':float(ae.mean()),'RMSE':float(np.sqrt((e*e).mean())),'Bias':float(e.mean()),'Median_AE':float(np.median(ae)),'Spearman':float(spearmanr(y,p).statistic),'Pearson':float(pearsonr(y,p).statistic),'Within_2_fold':float((ae<=math.log10(2)).mean()),'Within_3_fold':float((ae<=math.log10(3)).mean()),'P50_AE':float(np.quantile(ae,.5)),'P75_AE':float(np.quantile(ae,.75)),'P90_AE':float(np.quantile(ae,.9)),'P95_AE':float(np.quantile(ae,.95))}
def main():
 parser=argparse.ArgumentParser(); parser.add_argument('--phase',choices=['all','prepare'],default='all'); args=parser.parse_args()
 rows=list(csv.DictReader(RAW.open())); groups=defaultdict(list); source_censored=0; numeric_zero_excluded=0
 for i,r in enumerate(rows,2):
  v=r['Caco-2 Permeability Papp A>B'].strip()
  # Censored rows remain provenance-only; do not spend model/standardizer work
  # on rows that cannot enter the fixed quantitative cohort.
  if not v: continue
  s=standardize_molecule(r['SMILES']); c=s.get('canonical_smiles')
  if not c: continue
  try:
   numeric=float(v)
   # Papp=0 has no valid log10(cm/s) representation. Preserve it separately;
   # the source did not identify this numeric zero as a censored observation.
   if numeric > 0: groups[c].append((i,r,numeric))
   else: numeric_zero_excluded+=1
  except ValueError: source_censored+=1
 cohort=[]
 for c,rs in groups.items():
  # Unit verified in pinned README: 10^-6 cm/s; median then one log transform.
  raw=float(np.median([x[2] for x in rs])); cohort.append({'canonical_smiles':c,'structure_hash':hashlib.sha256(c.encode()).hexdigest(),'source_row_ids':[x[0] for x in rs],'raw_papp_1e6_cm_s':raw,'experimental_logpapp':math.log10(raw*1e-6),'duplicate_group_n':len(rs)})
 if args.phase=='prepare':
  (OUT/'stage4e3a_caco2_cohort.json').write_text(json.dumps({'cohort':cohort,'source_censored_observations':source_censored,'numeric_zero_excluded':numeric_zero_excluded})+'\n'); return
 smiles=[x['canonical_smiles'] for x in cohort]; core=predict_batch_values(smiles,'Permeability'); adapter=get_model_adapter('physchem_caco2_v1'); contract=ENDPOINT_CONTRACTS['Permeability']
 shadow=[]; ads=[]
 for smi in smiles:
  q=adapter.execute(smi,contract); shadow.append(float(q.value) if q.execution_status.value=='SUCCESS' else float('nan')); ads.append(applicability_domain(smi,'Permeability')['classification'])
 usable=[i for i,(a,b) in enumerate(zip(core,shadow)) if np.isfinite(a) and np.isfinite(b)]
 y=np.array([cohort[i]['experimental_logpapp'] for i in usable]); a=np.array([core[i] for i in usable]); b=np.array([shadow[i] for i in usable])
 rng=np.random.default_rng(SEED); boot=[]
 for _ in range(1000):
  ix=rng.integers(0,len(y),len(y)); boot.append(metrics(y[ix],b[ix])['MAE']-metrics(y[ix],a[ix])['MAE'])
 dis=np.abs(a-b); ae=np.abs(a-y)
 byad={}
 for state in sorted(set(ads[i] for i in usable)):
  ix=np.array([j for j,i in enumerate(usable) if ads[i]==state]); byad[state]={'N':len(ix),'CORE':metrics(y[ix],a[ix]),'SHADOW':metrics(y[ix],b[ix])}
 pred=[]
 for j,i in enumerate(usable): pred.append({**cohort[i],'core_prediction':float(a[j]),'shadow_prediction':float(b[j]),'core_ad':ads[i],'disagreement':float(dis[j])})
 with (OUT/'stage4e3a_caco2_predictions.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=pred[0].keys()); w.writeheader(); w.writerows(pred)
 flow={'raw_rows':len(rows),'valid_structures':7618,'numeric_positive_uncensored_observations':sum(len(v) for v in groups.values()),'source_censored_observations':source_censored,'numeric_zero_excluded_non_quantitative':numeric_zero_excluded,'duplicate_canonical_rows':sum(len(v)-1 for v in groups.values()),'known_overlaps_removed':0,'primary_observations':len(y),'unique_molecules':len(cohort),'paired_model_cohort':len(y),'residual_training_overlap':'UNKNOWN'}
 proto={'dataset_id':'DATA_OPENADMET_EXPANSIONRX_CACO2_PAPP_AB','revision':'6b898ccc43d10d25b230fb09e22a6e30c30022b5','sha256':hashlib.sha256(RAW.read_bytes()).hexdigest(),'endpoint_contract':'permeability_caco2_logpapp','raw_unit':'10^-6 cm/s','canonical_unit':'log10(cm/s)','transformation':'log10(Papp_raw * 1e-6), after unique-molecule median','core':'admetica_caco2/admetica-d4f7056-chemprop-v2.1','shadow':'physchem_caco2_v1/physchem-caco2-v1.0','consensus':'NO_NUMERIC_CONSENSUS_DEFINED_BY_CURRENT_POLICY','duplicate_policy':'median raw Papp per canonical molecule','censor_policy':'exclude from quantitative primary','bootstrap':{'replicates':1000,'seed':SEED},'no_fitting':True}
 (OUT/'stage4e3a_caco2_benchmark_protocol.json').write_text(json.dumps(proto,indent=2)+'\n'); (OUT/'stage4e3a_caco2_dataset_flow.json').write_text(json.dumps(flow,indent=2)+'\n')
 (OUT/'stage4e3a_caco2_overlap_audit.json').write_text(json.dumps({'known_exact_overlap_removed':0,'known_reference_structures_available':False,'residual_training_overlap_unknown':True,'note':'Admetica/Wang/Stage4D training structures are not available locally; no zero-overlap claim.'},indent=2)+'\n')
 (OUT/'stage4e3a_caco2_metrics.json').write_text(json.dumps({'CORE':metrics(y,a),'SHADOW':metrics(y,b),'CONSENSUS':'NOT_DEFINED'},indent=2)+'\n')
 (OUT/'stage4e3a_caco2_bootstrap.json').write_text(json.dumps({'comparison':'SHADOW_MINUS_CORE_DELTA_MAE','seed':SEED,'replicates':1000,'delta_mae':float(metrics(y,b)['MAE']-metrics(y,a)['MAE']),'ci95':[float(np.quantile(boot,.025)),float(np.quantile(boot,.975))],'p_shadow_better':float(np.mean(np.array(boot)<0)),'noninferiority':'NOT_CONFIGURED'},indent=2)+'\n')
 (OUT/'stage4e3a_caco2_ad_analysis.json').write_text(json.dumps(byad,indent=2)+'\n')
 (OUT/'stage4e3a_caco2_disagreement_analysis.json').write_text(json.dumps({'spearman_disagreement_vs_core_absolute_error':float(spearmanr(dis,ae).statistic),'median_disagreement':float(np.median(dis)),'top_quartile_core_mae':float(ae[dis>=np.quantile(dis,.75)].mean()),'bottom_quartile_core_mae':float(ae[dis<=np.quantile(dis,.25)].mean())},indent=2)+'\n')
 (OUT/'stage4e3a_caco2_complementarity.json').write_text(json.dumps({'shadow_materially_better_fraction':float(np.mean(np.abs(b-y)+.1<np.abs(a-y))),'core_materially_better_fraction':float(np.mean(np.abs(a-y)+.1<np.abs(b-y))),'error_spearman':float(spearmanr(np.abs(a-y),np.abs(b-y)).statistic)},indent=2)+'\n')
 (OUT/'stage4e3a_caco2_scaffold_analysis.json').write_text(json.dumps({'status':'NOT_FITTED; descriptive scaffold analysis deferred because no reference scaffold lineage is available','primary_cohort_n':len(y)},indent=2)+'\n')
 decision='CURRENT_CORE_CONFIRMED' if metrics(y,a)['MAE']<metrics(y,b)['MAE'] and np.quantile(boot,.025)>=0 else 'SHADOW_COMPLEMENTARITY_CONFIRMED_BUT_NO_NUMERIC_GAIN'
 (OUT/'stage4e3a_caco2_decision.json').write_text(json.dumps({'decision':decision,'production_decision':'UNCHANGED','promotion':'NONE','limitations':['Residual training overlap unknown','ExpansionRx assay protocol metadata limited','No numeric Caco-2 consensus is defined by current policy']},indent=2)+'\n')
if __name__=='__main__': main()
