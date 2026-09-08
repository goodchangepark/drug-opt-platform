"""Evaluate leakage-safe similarity-weighted kNN log-fu on canonical V72."""
import json, sys
from pathlib import Path
import numpy as np
from rdkit import DataStructs, RDLogger
from rdkit.Chem import AllChem
from sklearn.model_selection import GroupKFold
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT)); RDLogger.DisableLog('rdApp.warning')
from scripts.reconcile_ppb_benchmark_v72 import canonical_fu_rows, fold_indices, group_key, library, matrix, metric
def main():
 rows=canonical_fu_rows(library()); y=np.asarray([r['fu'] for r in rows],float); logy=np.log10(y); fps=[AllChem.GetMorganFingerprintAsBitVect(r['mol'],2,nBits=2048) for r in rows]; folds,groups=fold_indices(rows); pred=np.empty(len(rows))
 for tr,te in folds:
  sims=np.asarray([DataStructs.BulkTanimotoSimilarity(fps[i],[fps[j] for j in tr]) for i in te]); k=min(16,len(tr)); top=np.argpartition(sims,-k,axis=1)[:,-k:]
  for a,idx in enumerate(top): pred[te[a]]=np.average(logy[tr][idx],weights=np.maximum(sims[a,idx],.05)**2)
 fu=10**pred; result={'artifact':'PPB_KNN_V80','benchmark':'validation/ppb_canonical_benchmark_v72.json','n':len(y),'metric':metric(y,fu),'subgroups':{k:metric(y[m],fu[m]) for k,m in {'fu_lt_0_001':y<.001,'fu_0_001_to_0_01':(y>=.001)&(y<.01),'fu_0_01_to_0_1':(y>=.01)&(y<.1),'fu_0_1_to_0_5':(y>=.1)&(y<.5),'fu_ge_0_5':y>=.5}.items() if m.any()},'decision':'RETAIN_ONLY_IF_ROBUST_GAIN'}
 Path(ROOT/'validation/ppb_knn_v80.json').write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result['metric'],indent=2))
if __name__=='__main__': main()
