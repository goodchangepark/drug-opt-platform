"""Conservatively remove only positively identified runtime test projects."""
import json, sqlite3
from datetime import datetime, timezone
from pathlib import Path
from sqlalchemy import select, func
from backend.database import SessionLocal, engine
from backend.models import Project, Compound, CompoundVersion, ExternalExperimentalEvidence, PredictionRun
from backend.main import _delete_project_tree_rows
from backend.stabilization import classify_project

PROTECTED = {1, 3, 5, 300}
OUT = Path("validation/test_project_cleanup_audit.json")

def counts(s, pid):
    p = s.get(Project, pid)
    if not p: return None
    c = s.scalar(select(func.count(Compound.id)).where(Compound.project_id == pid)) or 0
    cvs = s.scalar(select(func.count(CompoundVersion.id)).join(Compound).where(Compound.project_id == pid)) or 0
    ev = s.scalar(select(func.count(ExternalExperimentalEvidence.id)).join(CompoundVersion).join(Compound).where(Compound.project_id == pid)) or 0
    return {"project_id": pid, "project_name": p.name, "compound_count": c, "version_count": cvs, "evidence_count": ev}

def main():
    s = SessionLocal()
    before_runs = [(r.id, r.model_version, r.stage, r.version_id) for r in s.scalars(select(PredictionRun).order_by(PredictionRun.id))]
    rows=[]; delete_ids=[]
    for p in s.scalars(select(Project).order_by(Project.id)):
        cls, reason = classify_project({"id":p.id,"name":p.name,"target":p.target,"description":p.description})
        row = counts(s,p.id) | {"classification":cls,"classification_evidence":reason,"decision":"PRESERVE"}
        if p.id not in PROTECTED and cls == "CONFIRMED_TEST":
            row["decision"]="DELETE"; delete_ids.append(p.id)
        rows.append(row)
    if PROTECTED - {p.id for p in s.scalars(select(Project))}:
        # Missing protected rows are a pre-existing state; never manufacture data.
        rows.append({"classification":"PROTECTED_MISSING","decision":"BLOCKED","missing":sorted(PROTECTED-{p.id for p in s.scalars(select(Project))})})
    if delete_ids:
        _delete_project_tree_rows(s, sorted(delete_ids)); s.commit()
    after_runs = [(r.id, r.model_version, r.stage, r.version_id) for r in s.scalars(select(PredictionRun).order_by(PredictionRun.id))]
    s.close()
    conn=sqlite3.connect("drug_opt.db"); cur=conn.cursor()
    cur.execute("PRAGMA integrity_check"); integrity=cur.fetchone()[0]
    cur.execute("PRAGMA foreign_key_check"); fk=cur.fetchall()
    cur.execute("SELECT id, count(*) FROM compounds WHERE project_id=300 GROUP BY project_id"); p300=cur.fetchone()
    cur.execute("SELECT count(DISTINCT cv.inchikey) FROM compound_versions cv JOIN compounds c ON c.id=cv.compound_row_id WHERE c.project_id=300")
    unique_keys=cur.fetchone()[0]
    conn.close()
    result={"generated_at":datetime.now(timezone.utc).isoformat(),"protected_project_ids":sorted(PROTECTED),"projects":rows,"deleted_ids":sorted(delete_ids),"prediction_runs_before":len(before_runs),"prediction_runs_after":len(after_runs),"prediction_runs_unchanged":before_runs==after_runs,"integrity_check":integrity,"foreign_key_check":fk,"project300_compounds":p300[1] if p300 else 0,"project300_unique_inchikeys":unique_keys}
    OUT.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps({k:result[k] for k in ("deleted_ids","prediction_runs_before","prediction_runs_after","prediction_runs_unchanged","integrity_check","project300_compounds","project300_unique_inchikeys")},indent=2))
if __name__ == "__main__": main()
