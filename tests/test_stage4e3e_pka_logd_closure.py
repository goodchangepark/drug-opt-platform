import json
from pathlib import Path
V=Path('validation')
def test_physchem_closure_is_conservative():
    d=json.loads((V/'stage4e3e_pka_logd_final_decisions.json').read_text())
    assert d['engine_v1_status']=='CLOSED'
    assert d['pka']['decision'].startswith('PKA_') and d['logd']['decision'].startswith('LOGD_')
    assert json.loads((V/'stage4e3e_current_pka_baseline.json').read_text())['evidence_class']=='RULE_ESTIMATE'
    assert json.loads((V/'stage4e3e_current_logd_baseline.json').read_text())['evidence_class']=='DERIVED_ESTIMATE'
    assert json.loads((V/'stage4e3e_pka_candidate_audit.json').read_text())['candidates'][0]['strict_load'] is False
