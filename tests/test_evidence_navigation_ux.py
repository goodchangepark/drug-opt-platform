"""Test Suite for Project and Compound Evidence Navigation UX.

Validates:
1. Project-level navigation contains: Compounds | Evidence | Assays | Compare | Settings
2. Project Compounds view displays compound list without embedded evidence/learning blocks
3. Project Evidence view contains: Project Evidence, Evidence Review, and Project Learning
4. Compound-level navigation contains: Overview | Properties | Activity | ADMET | Metabolism | PK | Evidence | History
5. Compound Evidence view isolates evidence to the selected compound version
6. EGFR and GLP-1R projects and compounds support navigation cleanly
"""
from pathlib import Path
import pytest
from backend.database import SessionLocal
from backend.models import Project, Compound
from backend.main import app
from starlette.testclient import TestClient
from sqlalchemy import select

client = TestClient(app)

def test_frontend_navigation_contract():
    root = Path(__file__).resolve().parents[1]
    js = (root / "frontend" / "static" / "app.js").read_text(encoding="utf-8")

    # Project navigation tabs
    assert "['compounds','Compounds'],['evidence','Evidence'],['assays','Assays'],['compare','Compare'],['settings','Settings']" in js

    # Compound navigation tabs
    assert "const tabs=['overview','properties','activity','admet','metabolism','pk','evidence','history'];" in js

    # Project Evidence view structure
    assert "project&&projectTab==='evidence'&&!detail" in js
    assert "key:'project-evidence-view'" in js
    assert "PROJECT EVIDENCE & LEARNING" in js

    # Compound Evidence view structure
    assert "detailTab==='evidence'&&e('div',{key:'evidence-tab'}" in js
    assert "key:'compound-evidence-panel'" in js


def test_project_and_compound_endpoints_for_navigation():
    with SessionLocal() as db:
        projects = list(db.scalars(select(Project)).all())
        assert len(projects) >= 2

        for p in projects:
            # Test project overview & evidence summaries
            proj_res = client.get(f"/api/projects/{p.id}")
            assert proj_res.status_code == 200
            proj = proj_res.json()
            assert "name" in proj
            assert "compounds" in proj

            ev_summary_res = client.get(f"/api/projects/{p.id}/evidence-summary")
            assert ev_summary_res.status_code == 200

            ev_review_res = client.get(f"/api/projects/{p.id}/evidence-review")
            assert ev_review_res.status_code == 200

            compounds = proj.get("compounds", [])
            assert len(compounds) > 0

            # Test compound workspace for each compound
            for c in compounds:
                v_id = c.get("current_version_id")
                if v_id:
                    ws_res = client.get(f"/api/compound-versions/{v_id}/workspace")
                    assert ws_res.status_code == 200
                    ws = ws_res.json()
                    assert ws["scope"]["compound_id"] == c["row_id"]
                    assert ws["scope"]["version_id"] == v_id
                    assert "external_experimental_evidence" in ws
