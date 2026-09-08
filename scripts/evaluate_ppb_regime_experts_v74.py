"""Train-fold-only PPB binding-regime experts on the frozen V72 cohort."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import ExtraTreesRegressor, RandomForestClassifier
from sklearn.model_selection import GroupKFold
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from scripts.evaluate_human_total_iv_cl_v56 import descriptor_vector

def metric(y,p):
 f=np.maximum(p/y,y/p); return {"n":int(len(y)),"aafe":float(f.mean()),"within_2_fold_pct":float((f<=2).mean()*100),"within_3_fold_pct":float((f<=3).mean()*100),"mae_log10":float(np.abs(np.log10(p)-np.log10(y)).mean())}
def main():
 lib=json.loads((ROOT/'validation/reference_library_v1_1000.json').read_text())['compounds']; rows=[]
 for x in lib:
  fu=(x.get('evidence') or {}).get('HUMAN_FUP') if isinstance(x.get('evidence'),dict) else None; m=Chem.MolFromSmiles(x.get('smiles',''))
  if m and fu and 0<float(fu)<=1: rows.append((x['inchikey'],m,float(fu)))
 rows.sort(); x=np.array([descriptor_vector(r[1]) for r in rows]); y=np.array([r[2] for r in rows]); ly=np.log10(y)
 groups=np.array([MurckoScaffold.MurckoScaffoldSmiles(mol=r[1]) or r[0] for r in rows]); oof=np.empty(len(y)); base=np.empty(len(y))
 for tr,te in GroupKFold(5).split(x,ly,groups):
  global_=ExtraTreesRegressor(n_estimators=400,min_samples_leaf=2,max_features=.8,random_state=74,n_jobs=-1).fit(x[tr],ly[tr]); base[te]=global_.predict(x[te])
  regime=np.digitize(ly[tr],[-2,-1]) # high binding, middle, freer
  clf=RandomForestClassifier(n_estimators=400,min_samples_leaf=4,class_weight='balanced',random_state=74,n_jobs=-1).fit(x[tr],regime)
  probs=clf.predict_proba(x[te]); pred=np.zeros(len(te))
  for cls in range(3):
   idx=np.where(regime==cls)[0]
   if len(idx)<12: expert=global_
   else: expert=ExtraTreesRegressor(n_estimators=400,min_samples_leaf=2,max_features=.8,random_state=74,n_jobs=-1).fit(x[tr][idx],ly[tr][idx])
   col=list(clf.classes_).index(cls) if cls in clf.classes_ else None
   if col is not None: pred+=probs[:,col]*expert.predict(x[te])
  oof[te]=pred
 pred=np.clip(10**oof,1e-4,1); basep=np.clip(10**base,1e-4,1)
 out={'artifact':'PPB_REGIME_EXPERTS_V74','canonical_benchmark':'PPB_CANONICAL_BENCHMARK_V72','n':len(y),'global_baseline':metric(y,basep),'soft_regime_experts':metric(y,pred),'subgroups':{k:metric(y[s],pred[s]) for k,s in {'fu_lt_0_01':y<.01,'fu_0_01_to_0_1':(y>=.01)&(y<.1),'fu_ge_0_1':y>=.1}.items()},'decision':'RETAIN_RESEARCH_CANDIDATE' if metric(y,pred)['aafe']<4.760840445074518*.95 else 'REJECT_NO_MEANINGFUL_IMPROVEMENT','production_engine':'drugopt-prediction-engine-v3@3.3.2'}
 (ROOT/'validation/ppb_regime_experts_v74.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n'); print(json.dumps(out,indent=2))
if __name__=='__main__': main()
