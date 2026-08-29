import json
from pathlib import Path

V=Path('validation')
def test_clearance_closure_artifacts_and_species_isolation():
    d=json.loads((V/'stage4e3d_clearance_final_decisions.json').read_text())
    assert d['engine_v1_status']=='CLOSED'
    assert set(d['species'])=={'hlm_intrinsic_clearance_scaled_log10','rlm_intrinsic_clearance_scaled_log10','mlm_intrinsic_clearance_scaled_log10'}
    assert d['dog']==d['monkey']==d['generic']=='MODEL_UNAVAILABLE'
    p=json.loads((V/'stage4e3d_clearance_protocol.json').read_text())
    assert p['rules']['species_isolation'] and p['rules']['no_species_average']
    assert json.loads((V/'stage4e3d_clearance_dataset_flow.json').read_text())['status']=='NO_DEFENSIBLE_INDEPENDENT_RAW_BENCHMARK'
    assert json.loads((V/'stage4e3d_clearance_predictions.json').read_text())['production_db_untouched']
