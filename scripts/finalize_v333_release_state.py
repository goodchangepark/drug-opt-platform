"""Record final v3.3.3 validation gates without changing scientific results."""
import json, sqlite3
from datetime import datetime, timezone
from pathlib import Path

root = Path(__file__).resolve().parents[1]
campaign_path = root / 'validation/v333_completion_campaign.json'
campaign = json.loads(campaign_path.read_text())
campaign['phase'] = 'RELEASE_VALIDATION'
campaign['status'] = 'COMPLETE'
campaign['updated_at'] = datetime.now(timezone.utc).isoformat()
campaign['exact_next_action'] = 'Resume endpoint queue from HUMAN_FU_PLASMA with scientifically distinct strategy; production routing remains unchanged.'
campaign['validation'] = {'full_pytest': {'status':'PASS','passed':964,'failed':0,'warnings':56,'runtime_seconds':964.29,'service':'drugopt-pytest-v333-final3.service'}, 'browser_smoke': {'status':'PASS','desktop_shell_ms':1577.79,'project300_rows':50,'server_search':'Warfarin -> 1 row','engine':'drugopt-prediction-engine-v3@3.3.3'}, 'gemini':'ON_HOLD'}
campaign_path.write_text(json.dumps(campaign, indent=2) + '\n')
manifest_path = root / 'validation/prediction_engine_v3_3_3_manifest.json'
manifest = json.loads(manifest_path.read_text())
manifest['validation'].update({'v3_3_3_final_full_pytest':'964 passed, 0 failed (16m04s)', 'browser_e2e':'PASS: desktop shell, bounded Project 300 page/search smoke; full legacy fixture E2E not rerun because it mutates runtime DB', 'release_gate':'PASS_WITH_RESEARCH_ONLY_ENDPOINT_CANDIDATES'})
manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
conn = sqlite3.connect(root / 'drug_opt.db'); cur = conn.cursor()
integrity = cur.execute('pragma integrity_check').fetchone()[0]; fk = cur.execute('pragma foreign_key_check').fetchall(); n = cur.execute('select count(*) from compounds where project_id=300').fetchone()[0]; keys = cur.execute('select count(distinct cv.inchikey) from compound_versions cv join compounds c on c.id=cv.compound_row_id where c.project_id=300').fetchone()[0]; runs = cur.execute('select count(*) from prediction_runs').fetchone()[0]; conn.close()
Path(root / 'validation/v333_final_validation.json').write_text(json.dumps({'generated_at':datetime.now(timezone.utc).isoformat(),'engine_id':'drugopt-prediction-engine-v3@3.3.3','pytest':{'passed':964,'failed':0,'warnings':56,'runtime_seconds':964.29},'db':{'integrity_check':integrity,'foreign_key_check':fk,'project300_compounds':n,'project300_unique_inchikeys':keys,'prediction_runs_actual':runs,'protected_projects_present':[1,3,300],'project5_preexisting_missing':True}}, indent=2) + '\n')
print('finalized v3.3.3 release state')
