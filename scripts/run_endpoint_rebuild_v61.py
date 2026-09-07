"""Endpoint-specific rebuild experiments over the frozen 1000-compound library."""
from __future__ import annotations
import json, hashlib, sys
from pathlib import Path
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupKFold
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import evaluate_human_total_iv_cl_v56 as base

ROOT=Path(__file__).resolve().parents[1]
LIB=json.loads((ROOT/'validation/reference_library_v1_1000.json').read_text())['compounds']
from backend.clinical_pk_cohort import CLINICAL_PK_COHORT
N30={Chem.MolToSmiles(Chem.MolFromSmiles(x.smiles),canonical=True) for x in CLINICAL_PK_COHORT}
N30M=[Chem.MolFromSmiles(x.smiles) for x in CLINICAL_PK_COHORT]

def rows(target, external_only=False):
    out=[]
    for x in LIB:
        if external_only and x['role'] != 'DIVERSITY_SELECTED_EXTERNAL':
            continue
        e=x['evidence']; y=e.get(target) if isinstance(e,dict) else None
        if y is None and isinstance(e,list):
            for o in e:
                if o.get('canonical_endpoint_id') in {target,'HUMAN_PPB'} and o.get('normalized_value') is not None: y=o['normalized_value']; break
        if y is None: continue
        m=Chem.MolFromSmiles(x['smiles'])
        if m and float(y)>0: out.append({'key':x['inchikey'],'canonical_smiles':Chem.MolToSmiles(m,canonical=True),'mol':m,'y':float(y)})
    return out
def clean(rs):
    nfp=[base.FINGERPRINT.GetFingerprint(m) for m in N30M]
    return [r for r in rs if r['canonical_smiles'] not in N30 and not any(DataStructs.TanimotoSimilarity(base.FINGERPRINT.GetFingerprint(r['mol']),f)>=.999999 for f in nfp)]
def cv(rs):
    X={m:base.build_matrix(rs,m) for m in ('descriptors','fingerprints','combined')}; y=np.log10([r['y'] for r in rs]); g=np.array([MurckoScaffold.MurckoScaffoldSmiles(mol=r['mol']) or r['canonical_smiles'] for r in rs]); ns=min(5,len(np.unique(g))); out={}
    for name,model in base.models().items():
        mode='fingerprints' if 'fingerprints' in name else ('combined' if 'combined' in name else 'descriptors'); fs=[]
        for tr,te in GroupKFold(ns).split(X[mode],y,g):
            model.fit(X[mode][tr],y[tr]); fs.append(base.metric(10**y[te],10**model.predict(X[mode][te])))
        out[name]={'features':mode,'mean_aafe':float(np.mean([x['aafe'] for x in fs])),'mean_within_2':float(np.mean([x['within_2_fold_pct'] for x in fs])),'folds':fs}
    return out
def fit(rs, best, test):
    mode='fingerprints' if 'fingerprints' in best else ('combined' if 'combined' in best else 'descriptors'); m=base.models()[best]; m.fit(base.build_matrix(rs,mode),np.log10([r['y'] for r in rs])); p=10**m.predict(base.build_matrix(test,mode)); return m,mode,p
def learning(rs,best):
    vals={}
    for n in (100,250,500,len(rs)):
        if n>len(rs): continue
        sub=rs[:n]; X=base.build_matrix(sub,'descriptors'); y=np.log10([r['y'] for r in sub]); g=np.array([MurckoScaffold.MurckoScaffoldSmiles(mol=r['mol']) or r['canonical_smiles'] for r in sub]); ns=min(3,len(np.unique(g))); fs=[]
        for tr,te in GroupKFold(ns).split(X,y,g):
            m=base.models()[best]; m.fit(X[tr],y[tr]); fs.append(base.metric(10**y[te],10**m.predict(X[te])))
        vals[str(n)]={'aafe':float(np.mean([x['aafe'] for x in fs])),'within_2_fold_pct':float(np.mean([x['within_2_fold_pct'] for x in fs]))}
    return vals
def main():
    out={'artifact':'endpoint_model_rebuild_v6_1','library':'REFERENCE_LIBRARY_V1','n30_excluded_from_fit':True,'experiments':{},'production':'UNCHANGED'}
    # VDss and total IV CL are sourced from PKSmart evidence in the manifest.
    for target,label in [('HUMAN_VDSS','VDSS'),('HUMAN_TOTAL_IV_CL','TOTAL_IV_CL')]:
        allr=rows(target, external_only=True); dev=clean(allr); cvr=cv(dev); best=min(cvr,key=lambda k:cvr[k]['mean_aafe']);
        locked=[]
        for x in CLINICAL_PK_COHORT:
            if (target=='HUMAN_TOTAL_IV_CL' and x.route=='IV' and x.cl_systemic_ml_min_kg) or (target=='HUMAN_VDSS' and x.vdss_l_kg):
                locked.append({'canonical_smiles':Chem.MolToSmiles(Chem.MolFromSmiles(x.smiles),canonical=True),'mol':Chem.MolFromSmiles(x.smiles),'y':float(x.cl_systemic_ml_min_kg if target=='HUMAN_TOTAL_IV_CL' else x.vdss_l_kg),'name':x.compound_name})
        model,mode,pred=fit(dev,best,locked); lockmetrics=base.metric(np.array([r['y'] for r in locked]),pred)
        out['experiments'][label]={'target':target,'available_n':len(allr),'development_n':len(dev),'model_selection_n':len(dev),'locked_validation_n':len(locked),'cv':cvr,'selected':best,'locked_metrics':lockmetrics,'learning_curve':learning(dev,best),'source_status':'WITHIN_SOURCE_VALIDATION_ONLY','decision':'MODEL_CANDIDATE_ONLY','reason':'No independent source validation; source-qualified endpoint semantics retained separately.'}
    # PPB uses exact canonical HUMAN_PPB records from the existing backbone.
    ppb=rows('HUMAN_PPB'); ppbdev=clean(ppb); ppbcv=cv(ppbdev); ppbbest=min(ppbcv,key=lambda k:ppbcv[k]['mean_aafe'])
    out['experiments']['PPB_FU']={'target':'HUMAN_PPB','available_n':len(ppb),'development_n':len(ppbdev),'model_selection_n':len(ppbdev),'cv':ppbcv,'selected':ppbbest,'learning_curve':learning(ppbdev,ppbbest),'strata':'fu<0.01 / 0.01-0.1 / >0.1 requires source-level fu records; current backbone target is %bound and was not silently converted.','decision':'RESEARCH_CANDIDATE_ONLY','reason':'Cross-source/locked endpoint validation is not available in this cycle; current production PPB route unchanged.'}
    out['model_hashes']={'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (ROOT/'validation/endpoint_model_rebuild_v6_1.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:{'selected':v.get('selected'),'locked_metrics':v.get('locked_metrics'),'development_n':v.get('development_n')} for k,v in out['experiments'].items()},indent=2))
if __name__=='__main__': main()
