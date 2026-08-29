import json
from pathlib import Path
v=Path(__file__).resolve().parents[1]/'validation'
def put(n,o):
    (v/n).write_text(json.dumps(o,indent=2)+'\n')
put('stage4e3e_pka_logd_protocol.json',{'artifact':'STAGE4E3E_PKA_LOGD_PROTOCOL','version':'stage4e3e-v1','status':'FROZEN','candidate_limits':{'pka':2,'logd':2},'current':{'pka':'RULE_ESTIMATE','logd':'DERIVED_ESTIMATE'},'rules':['strict checkpoint loading','no retraining or fitting','site/micro/macro semantics preserved','logP is not logD','pH-specific data required','production unchanged']})
put('stage4e3e_current_pka_baseline.json',{'implementation':'IonizationEngine_v1 SMARTS structural rules','evidence_class':'RULE_ESTIMATE','model_id':'ionization_smarts_rules_v1','version':'stage4c4-ionization-v1','semantics':['acidic/basic centers','atom indices','multiple centers','ampholyte/zwitterion classification'],'limitations':['not quantitatively ML validated','typical accuracy ±1–2 pKa units','polyprotic and complex systems high uncertainty']})
put('stage4e3e_current_logd_baseline.json',{'implementation':'RDKit Crippen cLogP + rule/experimental pKa + simplified monoprotic Henderson-Hasselbalch','evidence_class':'DERIVED_ESTIMATE','endpoint':'logD pH 7.4','model_id':'henderson_hasselbalch_logd_v1','version':'stage4c4-ionization-v1','limitations':['not standalone quantitative ML','polyprotic/zwitterion uncertainty','cLogP never relabeled logD']})
put('stage4e3e_pka_candidate_audit.json',{'candidates':[{'id':'MODEL_PKASOLVER_LITE','license':'MIT source/weights','checkpoint':'SOURCE_EMBEDDED_WEIGHTS','strict_load':False,'result':'NO_GO_REPRODUCIBILITY','reason':'legacy PyG state-dict incompatible; partial loading prohibited'},{'id':'MODEL_PKALEARN_GNN','license':'code MIT; weights/data unresolved','checkpoint':'NO_CHECKPOINT_VERIFIED','result':'NO_GO_REPRODUCIBILITY'}]})
put('stage4e3e_pka_checkpoint_manifest.json',{'pkasolver_lite':{'classification':'SOURCE_EMBEDDED_WEIGHTS','strict_load':'FAILED'},'pkaLearn':{'classification':'NO_CHECKPOINT','hash':'NOT_AVAILABLE'}})
put('stage4e3e_pka_arm64.json',{'pkasolver_lite':{'environment':'isolated ARM64 Python3.11 torch 2.8 CPU PyG 2.8.0.post1','build':'successful','load':'failed strict','decision':'NO_GO_REPRODUCIBILITY'},'pkaLearn':{'decision':'NOT_TESTED_UNRESOLVED_CHECKPOINT'}})
put('stage4e3e_pka_external_validation.json',{'status':'NOT_AVAILABLE','reason':'No site-resolved licensed independent cohort.'})
put('stage4e3e_logd_candidate_audit.json',{'candidates':[],'status':'NO_QUALIFIED_QUANTITATIVE_MODEL'})
put('stage4e3e_logd_dataset_audit.json',{'dataset':'DATA_LOGD74_1130','status':'NO_GO_LICENSE_FOR_ENGINE_V1_VALIDATION','reason':'Reuse/license terms unresolved'})
put('stage4e3e_logd_external_validation.json',{'status':'NOT_PERFORMED'})
put('stage4e3e_pka_logd_final_decisions.json',{'engine_v1_status':'CLOSED','pka':{'decision':'PKA_NO_REPRODUCIBLE_QUANTITATIVE_MODEL_RULE_ESTIMATE_FROZEN','reliability':'LOW-MEDIUM'},'logd':{'decision':'LOGD_NO_QUALIFIED_QUANTITATIVE_MODEL_DERIVED_ESTIMATE_FROZEN','reliability':'LOW-MEDIUM'},'production':'UNCHANGED'})
