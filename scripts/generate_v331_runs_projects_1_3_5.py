import urllib.request
import json
import sqlite3
import time

conn = sqlite3.connect('drug_opt.db')
c = conn.cursor()

c.execute('''
    SELECT c.project_id, p.name as proj_name, c.id, c.compound_id, c.name, cv.id as version_id
    FROM compounds c
    JOIN projects p ON c.project_id = p.id
    JOIN compound_versions cv ON cv.compound_row_id = c.id AND cv.version_number = c.current_version
    WHERE c.project_id IN (1, 3, 5)
    ORDER BY c.project_id, c.id
''')
compounds = c.fetchall()

for row in compounds:
    pid, pname, cid, clabel, cname, vid = row
    # Check if a v3.3.1 run already exists for this version
    c.execute('''
        SELECT id FROM prediction_runs
        WHERE version_id = ? AND model_version = '3.3.1'
    ''', (vid,))
    v331 = c.fetchone()
    if v331:
        print(f"[SKIP] Proj {pid} | Comp {cid} ({clabel}) already has v3.3.1 run ID {v331[0]}")
        continue

    print(f"[RUNNING] Proj {pid} | Comp {cid} ({clabel}) v3.3.1 prediction workflow...", flush=True)
    t0 = time.time()
    url = f"http://127.0.0.1:8765/api/compounds/{cid}/predict-workflow"
    req = urllib.request.Request(url, data=b"", method="POST")
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            dur = time.time() - t0
            data = json.loads(resp.read().decode())
            print(f"  --> Succeeded in {dur:.1f}s: Run ID {data.get('prediction_run_id')}, engine={data.get('engine_version')}, status={data.get('status')}", flush=True)
    except Exception as e:
        print(f"  --> Exception: {e}", flush=True)

conn.close()
