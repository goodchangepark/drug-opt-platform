"""Tests for Compound-Scoped Evidence, Persistent Prediction & Refined Search v4.8."""
from __future__ import annotations

import json
from pathlib import Path

import pytest
from fastapi.testclient import TestClient
from sqlalchemy import select

from backend.database import SessionLocal
from backend.experimental_refinement import (
    REFINEMENT_POLICY_VERSION,
    refine_scientific_observation,
    reprocess_all_persisted_evidence,
)
from backend.main import app
from backend.models import (
    Compound,
    CompoundVersion,
    ExternalExperimentalEvidence,
    Project,
)
from backend.platform_info import version_history

client = TestClient(app)


def test_v48_version_history():
    history = version_history()
    assert any(h["version"] == "v4.8" for h in history)
    v48 = next(h for h in history if h["version"] == "v4.8")
    assert "Compound-Scoped Evidence" in v48["milestone"]


def test_experimental_refinement_context_hierarchy():
    obs = {
        "raw_endpoint_name": "Activity",
        "raw_value": 0.45,
        "raw_unit": "nM",
        "assay_conditions_json": {
            "table_header": "GLP-1R cAMP EC50 (nM)",
            "footnotes": "Human recombinant GLP-1 receptor expressing CHO cells",
        },
    }
    refined = refine_scientific_observation(obs)
    assert refined["qualification"] == "AUTO_QUALIFIED"
    assert refined["measurement_type"] == "EC50"
    assert refined["species"] == "HUMAN"
    assert refined["canonical_endpoint_id"] == "ACTIVITY_EC50"


def test_experimental_refinement_review_required_reason():
    obs = {
        "raw_endpoint_name": "Clearance",
        "raw_value": 15.2,
        "raw_unit": "",
        "assay_conditions_json": {},
    }
    refined = refine_scientific_observation(obs)
    assert refined["qualification"] == "REVIEW_REQUIRED"
    assert refined["unresolved_reason"] in {"UNIT_MISSING", "SPECIES_MISSING", "ENDPOINT_AMBIGUOUS"}


def test_compound_evidence_isolation():
    with SessionLocal() as db:
        glp1_proj = db.scalar(select(Project).where(Project.name.like("%GLP-1%")))
        if glp1_proj:
            comp_rows = list(db.scalars(
                select(Compound).where(Compound.project_id == glp1_proj.id)
            ).all())
            for c in comp_rows:
                res = client.get(f"/api/compounds/{c.id}")
                if res.status_code == 200 and res.json().get("version"):
                    v_id = res.json()["version"]["id"]
                    ws = client.get(f"/api/compound-versions/{v_id}/workspace")
                    assert ws.status_code == 200
                    ws_data = ws.json()
                    assert ws_data["scope"]["compound_id"] == c.id
                    for ev in ws_data.get("external_experimental_evidence", []):
                        ev_row = db.get(ExternalExperimentalEvidence, ev["id"])
                        ev_ver = db.get(CompoundVersion, ev_row.compound_version_id)
                        assert ev_ver.compound_row_id == c.id


def test_reprocess_all_persisted_evidence():
    with SessionLocal() as db:
        stats = reprocess_all_persisted_evidence(db)
        assert stats["total"] > 0
        assert stats["unusable"] == 0
        assert stats["resolved_count"] > 0
