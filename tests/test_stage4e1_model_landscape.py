"""Stage 4E-1 artifact consistency; this stage must not alter runtime policy."""
from __future__ import annotations
import json
from pathlib import Path
from backend.endpoint_strategy_registry import StrategyType, get_all_strategies

ROOT=Path(__file__).resolve().parents[1]
def load(name): return json.loads((ROOT/'validation'/name).read_text())

def test_1_all_current_policies_are_represented():
    baseline=load('stage4e1_current_model_baseline.json')
    assert baseline['endpoint_count']==49==len(get_all_strategies())
    assert {row['endpoint_name'] for row in baseline['endpoints']}==set(get_all_strategies())
    assert baseline['reconciliation']['contradictions']==[]

def test_2_landscape_does_not_claim_installation_or_production_change():
    baseline=load('stage4e1_current_model_baseline.json')
    candidates=load('stage4e1_candidate_model_landscape.json')
    assert baseline['production_changed'] is False
    assert candidates['installed_or_registered'] is False

def test_3_candidate_records_have_required_planning_fields():
    required={'candidate_id','endpoint_id','model_name','model_family','checkpoint_available','code_available','inference_code_available','license_code','license_checkpoint','license_dataset','training_overlap_risk','arm64_feasibility','recommended_action','source_ids'}
    for candidate in load('stage4e1_candidate_model_landscape.json')['candidates']:
        assert required <= set(candidate)
        assert candidate['recommended_action']

def test_4_sources_resolve_for_every_candidate_and_dataset():
    source_ids={row['source_id'] for row in load('stage4e1_source_manifest.json')['sources']}
    for row in load('stage4e1_candidate_model_landscape.json')['candidates']+load('stage4e1_dataset_landscape.json')['datasets']:
        assert set(row['source_ids']) <= source_ids

def test_5_go_stage4e2_candidates_would_have_all_required_gates():
    for row in load('stage4e1_candidate_model_landscape.json')['candidates']:
        if row['recommended_action']=='GO_STAGE4E2':
            assert row['source_ids'] and row['checkpoint_available'] is True and row['code_available']
            assert row['license_code'] not in {'LICENSE_UNCLEAR','RESTRICTED'}
            assert row['license_checkpoint'] not in {'LICENSE_UNCLEAR','RESTRICTED'}

def test_6_known_noncommercial_and_endpoint_mismatch_candidates_cannot_be_pilots():
    rows={row['candidate_id']:row for row in load('stage4e1_candidate_model_landscape.json')['candidates']}
    assert rows['MODEL_BAYESHERG']['recommended_action']=='NO_GO_LICENSE'
    assert rows['MODEL_MMTKPRED_TRANSPORTER']['recommended_action']=='NO_GO_ENDPOINT_MISMATCH'

def test_7_model_unavailable_policies_remain_unavailable():
    baseline={row['endpoint_name']:row for row in load('stage4e1_current_model_baseline.json')['endpoints']}
    for name,policy in get_all_strategies().items():
        if policy.primary_strategy==StrategyType.MODEL_UNAVAILABLE:
            assert baseline[name]['current_strategy']=='MODEL_UNAVAILABLE'

def test_8_pilot_plan_references_known_candidates_and_datasets():
    plan=load('stage4e1_stage4e2_pilot_plan.json')
    candidates={row['candidate_id'] for row in load('stage4e1_candidate_model_landscape.json')['candidates']}
    datasets={row['dataset_id'] for row in load('stage4e1_dataset_landscape.json')['datasets']}
    assert set(plan['pilot_model_ids']) <= candidates
    assert set(plan['pilot_dataset_ids']) <= datasets
    assert plan['no_installation_in_stage4e1'] is True

def test_9_high_overlap_tdc_benchmarks_are_not_mislabeled_independent():
    datasets={row['dataset_id']:row for row in load('stage4e1_dataset_landscape.json')['datasets']}
    assert datasets['DATA_TDC_CACO2_WANG']['recommended_action']=='NO_GO_OVERLAP_RISK'
    assert datasets['DATA_TDC_HERG']['recommended_action']=='NO_GO_OVERLAP_RISK'
