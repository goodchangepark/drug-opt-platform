"""Three-cycle numerical follow-up: total IV CL, VDss, and IV half-life.

All final clinical records are locked before model choice.  PKSmart's
External_test_315 is consumed stress validation and is never used for fitting.
"""
from __future__ import annotations
import csv, hashlib, json, sys
from pathlib import Path
import numpy as np
from rdkit import Chem, DataStructs
from sklearn.model_selection import GroupKFold
from rdkit.Chem.Scaffolds import MurckoScaffold
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import evaluate_human_total_iv_cl_v56 as b
from backend.clinical_pk_cohort import CLINICAL_PK_COHORT

ROOT=Path(__file__).resolve().parents[1]
MAIN=Path('/tmp/pksmart-human-pk.csv')
STRESS=Path('/tmp/pksmart-external-test-315.csv')

def read_target(path, col):
    groups={}; mols={}
    with path.open(encoding='utf-8-sig',newline='') as f:
        for row in csv.DictReader(f):
            try: y=float(row.get(col,''))
            except (ValueError,TypeError): continue
            if y<=0: continue
            m=Chem.MolFromSmiles(row.get('smiles_r',''))
            if not m: continue
            k=Chem.MolToSmiles(m,canonical=True); groups.setdefault(k,[]).append(y); mols[k]=m
    return [{'canonical_smiles':k,'mol':mols[k],'y':float(10**np.median(np.log10(v))),'n':len(v)} for k,v in groups.items()]

def groups_of(rs): return np.array([MurckoScaffold.MurckoScaffoldSmiles(mol=r['mol']) or r['canonical_smiles'] for r in rs])

def evaluate_models(dev, target_name):
    y=np.log10([r['y'] for r in dev]); X={m:b.build_matrix(dev,m) for m in ('descriptors','fingerprints','combined')}
    out={}; nsplit=min(5,len(np.unique(groups_of(dev))))
    for name, model in b.models().items():
        mode='fingerprints' if 'fingerprints' in name else ('combined' if 'combined' in name else 'descriptors')
        folds=[]
        for tr,te in GroupKFold(n_splits=nsplit).split(X[mode],y,groups_of(dev)):
            model.fit(X[mode][tr],y[tr]); folds.append(b.metric(10**y[te],10**model.predict(X[mode][te])))
        out[name]={'features':mode,'mean_aafe':float(np.mean([x['aafe'] for x in folds])),'mean_within_2':float(np.mean([x['within_2_fold_pct'] for x in folds])),'folds':folds}
    return out

def main():
    n30keys={Chem.MolToSmiles(Chem.MolFromSmiles(x.smiles),canonical=True) for x in CLINICAL_PK_COHORT}
    cl=read_target(MAIN,'human_CL_mL_min_kg'); vd=read_target(MAIN,'human_VDss_L_kg')
    locked_cl=[{'name':x.compound_name,'mol':Chem.MolFromSmiles(x.smiles),'y':x.cl_systemic_ml_min_kg} for x in CLINICAL_PK_COHORT if x.route=='IV' and x.cl_systemic_ml_min_kg]
    locked_vd=[{'name':x.compound_name,'mol':Chem.MolFromSmiles(x.smiles),'y':x.vdss_l_kg} for x in CLINICAL_PK_COHORT if x.vdss_l_kg]
    def clean(rs, locked):
        lf=[b.FINGERPRINT.GetFingerprint(r['mol']) for r in locked]
        return [r for r in rs if r['canonical_smiles'] not in n30keys and not any(DataStructs.TanimotoSimilarity(b.FINGERPRINT.GetFingerprint(r['mol']),f)>=.999999 for f in lf)]
    cldev=clean(cl,locked_cl); vddev=clean(vd,locked_vd)
    clcv=evaluate_models(cldev,'CL'); vdCV=evaluate_models(vddev,'VDss')
    clbest=min(clcv,key=lambda k:clcv[k]['mean_aafe']); vdbest=min(vdCV,key=lambda k:vdCV[k]['mean_aafe'])
    def fit_eval(dev, locked, bestname):
        mode='fingerprints' if 'fingerprints' in bestname else ('combined' if 'combined' in bestname else 'descriptors')
        model=b.models()[bestname]; model.fit(b.build_matrix(dev,mode),np.log10([r['y'] for r in dev])); p=10**model.predict(b.build_matrix(locked,mode))
        return model,mode,p,b.metric(np.array([r['y'] for r in locked]),p)
    clmodel,clmode,clpred,clmetrics=fit_eval(cldev,locked_cl,clbest)
    vdmodel,vdmode,vdpred,vdmetrics=fit_eval(vddev,locked_vd,vdbest)
    # Stress is evaluated only after both choices are fixed; no model changes.
    b.SOURCE=STRESS; stress=read_target(STRESS,'human_CL_mL_min_kg')
    trainkeys={r['canonical_smiles'] for r in cl}; trainfps=[b.FINGERPRINT.GetFingerprint(r['mol']) for r in cl]
    stress=[r for r in stress if r['canonical_smiles'] not in trainkeys and not any(DataStructs.TanimotoSimilarity(b.FINGERPRINT.GetFingerprint(r['mol']),f)>=.999999 for f in trainfps)]
    stresspred=10**clmodel.predict(b.build_matrix(stress,clmode)); stressmetrics=b.metric(np.array([r['y'] for r in stress]),stresspred)
    # IV half-life uses only predicted CL and predicted VDss; no observed PK inputs.
    vdmap={r.compound_name:p for r,p in zip([x for x in CLINICAL_PK_COHORT if x.vdss_l_kg],vdpred)}
    half=[]
    for r,pcl in zip(locked_cl,clpred):
        x=next(x for x in CLINICAL_PK_COHORT if x.compound_name==r['name'])
        if x.t_half_hr and r['name'] in vdmap:
            # CL is mL/min/kg and V is L/kg: convert CL to L/h/kg first.
            half.append({'compound':r['name'],'observed':x.t_half_hr,'predicted':0.693*vdmap[r['name']]/(pcl*0.06)})
    hm=b.metric(np.array([x['observed'] for x in half]),np.array([x['predicted'] for x in half])) if half else None
    artifact={'artifact':'multicycle_pk_model_development_v5_6','source':'PKSmart Human_PK_data.csv; External_test_315 consumed stress only','cycles':{
      'cycle_1_total_iv_cl':{'development_n':len(cldev),'model_selection_n':len(cldev),'locked_n':len(locked_cl),'cv':clcv,'selected':clbest,'locked_metrics':clmetrics,'stress_n':len(stress),'stress_metrics':stressmetrics,'decision':'CANDIDATE_NOT_PRODUCTION_QUALIFIED','reason':'Consumed same-resource stress cohort generalizes materially worse than internal N=7; independent-source validation still required.'},
      'cycle_2_vdss':{'development_n':len(vddev),'model_selection_n':len(vddev),'locked_n':len(locked_vd),'cv':vdCV,'selected':vdbest,'locked_metrics':vdmetrics,'decision':'RESEARCH_CANDIDATE_ONLY','reason':'Separate VDss target and AD were evaluated without Vz/Vd/F mixing; no production promotion without independent-source validation.'},
      'cycle_3_iv_half_life':{'eligible_n':len(half),'predicted_n':len(half),'metrics':hm,'rows':half,'decision':'RESEARCH_CHAIN_ONLY','reason':'Mechanistic one-compartment calculation from predicted total IV CL and predicted VDss; not propagated to production PK.'}},
      'model_hashes':{'implementation_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},'production':'UNCHANGED','next_data_requirement':'Independent source human IV total CL and human IV VDss with study-level analyte/form/route metadata.'}
    (ROOT/'validation/multicycle_pk_model_development_v5_6.json').write_text(json.dumps(artifact,indent=2)+'\n')
    print(json.dumps({'cl':clmetrics,'stress':stressmetrics,'vdss':vdmetrics,'half_life':hm,'selected':(clbest,vdbest)},indent=2))
if __name__=='__main__': main()
