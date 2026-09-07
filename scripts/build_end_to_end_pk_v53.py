"""Build the three-mode N=30 PK benchmark and evidence readiness map."""
from __future__ import annotations
import json, sqlite3, sys
from datetime import datetime, timezone
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from backend.clinical_pk_cohort import CLINICAL_PK_COHORT

REPORT=ROOT/'validation/expanded_pk_validation_report_v1.json'
OUT=ROOT/'validation/end_to_end_pk_baseline_v5_3.json'
READINESS=ROOT/'validation/pk_model_development_readiness_v5_3.json'
ERROR=ROOT/'validation/end_to_end_pk_error_propagation_v5_3.json'

def main():
    now=datetime.now(timezone.utc).isoformat(); old=json.loads(REPORT.read_text())
    db=sqlite3.connect(ROOT/'drug_opt.db')
    try:
        drugbank=db.execute('select count(*) from compounds where project_id=300').fetchone()[0]
        runs=db.execute('select count(*) from prediction_runs').fetchone()[0]
    finally: db.close()
    modes={
      'OBSERVED_PARAMETER_ASSISTED': {'complete_n':30,'predicted_n':30,'eligible_n':30,'metrics':old['pk_parameter_performance'],'status':'COMPLETE_WITH_OBSERVED_UPSTREAM_INPUTS'},
      'HYBRID_PREDICTION': {'complete_n':0,'predicted_n':0,'eligible_n':30,'metrics':None,'status':'END_TO_END_INCOMPLETE','reason':'No pre-registered hybrid upstream selection was evaluated on this frozen cohort; observed HLM/fu/VDss inputs cannot be relabeled as predicted.'},
      'FULL_END_TO_END_PREDICTION': {'complete_n':0,'predicted_n':0,'eligible_n':30,'metrics':None,'status':'END_TO_END_INCOMPLETE','reason':'No qualified fully predicted pKa, logD, fu, HLM Clint, renal, VDss, F and ka chain was available without observed target inputs.'}
    }
    records=[]; error_rows=[]
    for d in CLINICAL_PK_COHORT:
        records.append({'compound':d.compound_name,'dose_mg':d.dose_mg,'route':d.route,'regimen':d.regimen,'formulation':d.formulation,'population':d.population,'source':d.source_reference,'modes':{'OBSERVED_PARAMETER_ASSISTED':'COMPLETE_WITH_OBSERVED_UPSTREAM_INPUTS','HYBRID_PREDICTION':'END_TO_END_INCOMPLETE','FULL_END_TO_END_PREDICTION':'END_TO_END_INCOMPLETE'},'observed_upstream_parameters':{'fu':d.ppb_percent is not None,'HLM_Clint':d.cl_int_hlm_ml_min_kg is not None,'VDss':d.vdss_l_kg is not None,'CL':d.cl_systemic_ml_min_kg is not None,'renal_CL':d.cl_renal_ml_min_kg is not None,'F':d.f_oral_bioavailability is not None}})
        error_rows.append({'compound':d.compound_name,'stages':{'pKa':'NOT_EVALUATED','logD7.4':'NOT_EVALUATED','solubility':'NOT_EVALUATED','Caco2':'INPUT_TO_ASSISTED_ROUTE','fu':'OBSERVED_INPUT','HLM_Clint':'OBSERVED_INPUT','hepatic_CL':'DERIVED_FROM_OBSERVED_INPUT','renal_CL':'OBSERVED' if d.cl_renal_ml_min_kg is not None else 'MISSING','total_CL':'OBSERVED_TARGET_NOT_PREDICTED','VDss':'OBSERVED_INPUT','F':'OBSERVED' if d.f_oral_bioavailability is not None else 'MODEL/DERIVED','ka':'DERIVED_FROM_PERMEABILITY','AUC':'DOWNSTREAM_OUTPUT','Cmax':'DOWNSTREAM_OUTPUT','Tmax':'DOWNSTREAM_OUTPUT','half_life':'DOWNSTREAM_OUTPUT'},'leakage_safe_for_full_mode':False})
    baseline={'artifact':'end_to_end_pk_baseline_v5_3','created_at':now,'cohort_n':30,'prediction_engine':'drugopt-prediction-engine-v3@3.3.2','pk_engine':'drugopt-pk-engine-v1','frozen_engine_v1':{'id':'drugopt-prediction-engine-v1@1.0.0','hash':'12757ab197b5a70d8ea1754678d9a342ab0b6ea0d82f2896bebb767d686bbdeb'},'drugbank':drugbank,'prediction_runs':runs,'modes':modes,'records':records,'production_changed':False,'target_leakage':False}
    error={'artifact':'end_to_end_pk_error_propagation_v5_3','created_at':now,'cohort_n':30,'stages':['pKa','logD7.4','solubility','Caco2','fu','HLM_Clint','hepatic_CL','renal_CL','total_CL','VDss','F','ka','AUC','Cmax','Tmax','half_life'],'records':error_rows,'statistical_interpretation':'Stage-wise correlations are not identifiable for observed-input stages in this cohort. The artifact preserves that limitation instead of assigning causal error to fu/Clint/VDss. Existing local sensitivity coefficients remain diagnostic only.' ,'top_contributors':[{'stage':'CL','status':'PRIMARY_CURRENT_WEAKNESS','evidence':'hepatic CL AAFE 3.295'},{'stage':'renal_CL','status':'MISSING_OR_INCOMPLETE','count':sum(d.cl_renal_ml_min_kg is None for d in CLINICAL_PK_COHORT)},{'stage':'fu_inc','status':'UNKNOWN_PROVENANCE','count':30},{'stage':'Rb','status':'UNKNOWN_PROVENANCE','count':30},{'stage':'VDss','status':'OBSERVED_INPUT_NOT_SCORED','count':30}]}
    readiness={'artifact':'pk_model_development_readiness_v5_3','created_at':now,'rows':[
      {'endpoint':'Human total CL','experimental_n':30,'qualified_n':30,'context_completeness':'IV/oral semantics heterogeneous; component provenance incomplete','checkpoint_available':False,'trainable':False,'independent_validation_possible':True,'verdict':'MORE_DATA_REQUIRED','action':'Acquire IV CL plus renal/nonrenal component data.'},
      {'endpoint':'Human hepatocyte Clint','experimental_n':0,'qualified_n':0,'context_completeness':'NONE','checkpoint_available':False,'trainable':False,'independent_validation_possible':False,'verdict':'MORE_DATA_REQUIRED','action':'Acquire qualified human hepatocyte data first.'},
      {'endpoint':'Renal CL','experimental_n':15,'qualified_n':15,'context_completeness':'fe/CL coverage partial','checkpoint_available':False,'trainable':False,'independent_validation_possible':True,'verdict':'MORE_DATA_REQUIRED','action':'Routing/readiness only; no universal regression.'},
      {'endpoint':'VDss','experimental_n':30,'qualified_n':30,'context_completeness':'IV-derived vs Vz/Vd/F separation required','checkpoint_available':False,'trainable':False,'independent_validation_possible':True,'verdict':'MORE_DATA_REQUIRED','action':'Collect independent IV VDss and run leakage-safe evaluation.'},
      {'endpoint':'pKa','experimental_n':0,'qualified_n':0,'context_completeness':'site/type/method absent','checkpoint_available':False,'trainable':False,'independent_validation_possible':False,'verdict':'NO_VALIDATED_MODEL','action':'PKA_DATA_INSUFFICIENT.'},
      {'endpoint':'logD7.4','experimental_n':0,'qualified_n':0,'context_completeness':'exact pH-7.4 labels absent','checkpoint_available':False,'trainable':False,'independent_validation_possible':False,'verdict':'NO_VALIDATED_MODEL','action':'LOGD_DATA_INSUFFICIENT.'},
      {'endpoint':'CYP2C19/P-gp/BCRP quantitative','experimental_n':0,'qualified_n':0,'context_completeness':'quantitative assay/checkpoint absent','checkpoint_available':False,'trainable':False,'independent_validation_possible':False,'verdict':'NO_VALIDATED_MODEL','action':'Remain fail-closed.'},
      {'endpoint':'DrugBank 250→300','experimental_n':None,'qualified_n':None,'context_completeness':'targeted gaps not yet populated','checkpoint_available':None,'trainable':None,'independent_validation_possible':None,'verdict':'NOT_PRIORITY','action':'DO_NOT_EXPAND_YET.'}]}
    OUT.write_text(json.dumps(baseline,indent=2)+'\n'); ERROR.write_text(json.dumps(error,indent=2)+'\n'); READINESS.write_text(json.dumps(readiness,indent=2)+'\n')
    print(json.dumps({'outputs':[str(OUT),str(ERROR),str(READINESS)],'cohort_n':30,'drugbank':drugbank,'prediction_runs':runs},indent=2))
if __name__=='__main__': main()
