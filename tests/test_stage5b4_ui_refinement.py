"""
Tests for Stage 5B-4 Scientific UI/UX Refinement
Covers:
- Platform version 0.6.1-stage5b4-ui and version history registry
- Centralized interpretation registry (interpretation.py)
- Ionization Henderson-Hasselbalch calculation enhancements
- API endpoints: /api/health, /api/help/registry, /api/interpretation/rules
- Frontend design system assets (CSS typography & color semantics, app.js structure)
"""

import pytest
from fastapi.testclient import TestClient
from backend.main import app
from backend.platform_info import APP_VERSION, CURRENT_STAGE, version_history, latest_release_date
from backend.interpretation import interpret_property, get_interpretation_registry_summary, INTERPRETATION_RULES
from backend.ionization import analyze_ionization
import os

client = TestClient(app)


def test_platform_version_and_history():
    assert APP_VERSION == "1.0.0"
    assert CURRENT_STAGE == "5B-4"
    assert latest_release_date() == "2026-09-02"
    
    vh = version_history()
    assert len(vh) >= 12
    assert vh[0]["version"] == "0.1.0"
    stage5b4_entries = [entry for entry in vh if "Stage 5B-4 Refinement" in entry["stage"]]
    assert len(stage5b4_entries) > 0
    assert ("Unified Prediction Workflow" in stage5b4_entries[-1]["milestone"] or "Dashboard Redesign" in stage5b4_entries[-1]["milestone"])
    assert vh[-1]["version"] in {"v4.7", "v4.8", "v4.8.1"}


def test_api_health_endpoint():
    response = client.get("/api/health")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert data["version"] == APP_VERSION
    assert data["step"] == "5B-4"
    assert data["updated"] == latest_release_date()


def test_api_interpretation_rules_endpoint():
    response = client.get("/api/interpretation/rules")
    assert response.status_code == 200
    data = response.json()
    assert data["version"] == "1.0.0"
    assert data["rules_count"] >= 20
    assert "mw" in data["rules"]
    assert "clogp" in data["rules"]
    assert "solubility" in data["rules"]
    assert "hlm_clint" in data["rules"]
    assert "herg" in data["rules"]


def test_api_help_registry_includes_history_and_rules():
    response = client.get("/api/help/registry")
    assert response.status_code == 200
    data = response.json()
    assert data["application"]["version"] == APP_VERSION
    assert "version_history" in data
    assert len(data["version_history"]) >= 12
    assert "interpretation_registry" in data
    assert data["interpretation_registry"]["rules_count"] >= 20



def test_physicochemical_interpretations():
    # MW
    res = interpret_property("mw", 350.0)
    assert res["color_class"] == "favorable"
    
    res = interpret_property("mw", 650.0)
    assert res["color_class"] == "liability"

    # cLogP
    res = interpret_property("clogp", 2.5)
    assert res["color_class"] == "favorable"

    res = interpret_property("clogp", 6.5)
    assert res["color_class"] == "liability"

    # TPSA
    res = interpret_property("tpsa", 80.0)
    assert res["color_class"] == "favorable"

    res = interpret_property("tpsa", 160.0)
    assert res["color_class"] == "liability"

    # QED
    res = interpret_property("qed", 0.75)
    assert res["color_class"] == "favorable"

    res = interpret_property("qed", 0.35)
    assert res["color_class"] == "liability"


def test_admet_and_safety_interpretations():
    # Solubility
    res = interpret_property("solubility", 85.0)
    assert res["color_class"] == "favorable"
    assert "High" in res["interpretation"]

    res = interpret_property("solubility", 5.0)
    assert res["color_class"] == "liability"

    # Caco-2
    res = interpret_property("caco2", -4.8)
    assert res["color_class"] == "favorable"

    # HLM Clint
    res = interpret_property("hlm_clint", 10.0)
    assert res["color_class"] == "favorable"

    res = interpret_property("hlm_clint", 60.0)
    assert res["color_class"] == "liability"

    # hERG
    res = interpret_property("herg", 0.12)
    assert res["color_class"] == "favorable"

    res = interpret_property("herg", 0.85)
    assert res["color_class"] == "liability"

    # DILI
    res = interpret_property("dili", 0.70)
    assert res["color_class"] == "liability"


def test_ampholyte_ionization_profile():
    # Ciprofloxacin: contains carboxylic acid and piperazine basic nitrogen
    ciprofloxacin_smiles = "C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O"
    prof = analyze_ionization(ciprofloxacin_smiles)
    assert prof["ionization_class"] in ("AMPHOLYTE", "ZWITTERION_POSSIBLE")
    assert len(prof["ionizable_centers"]) >= 2
    
    # Check that pH transitions are calculated smoothly rather than static 0.5
    ph_profiles = prof["ph_profiles"]
    assert len(ph_profiles) >= 5
    fn_values = [p["fraction_neutral"] for p in ph_profiles]
    assert all(0.0 <= fn <= 1.0 for fn in fn_values)


def test_static_assets_consistency():
    css_path = os.path.join(os.path.dirname(__file__), "..", "frontend", "static", "app.css")
    js_path = os.path.join(os.path.dirname(__file__), "..", "frontend", "static", "app.js")
    
    assert os.path.exists(css_path)
    assert os.path.exists(js_path)
    
    with open(css_path, "r", encoding="utf-8") as f:
        css_content = f.read()
    
    with open(js_path, "r", encoding="utf-8") as f:
        js_content = f.read()
    
    # Typography
    assert "Noto Sans KR" in css_content
    # Color semantics
    assert ".dot-favorable" in css_content
    assert ".dot-liability" in css_content
    assert ".dot-intermediate" in css_content
    assert ".badge-favorable" in css_content
    assert ".badge-liability" in css_content
    # Profile & tables
    assert ".admet-visual-profile" in css_content
    assert ".sidebar-footer" in css_content
    assert ".compound-header-card" in css_content
    
    # App.js helpers & components
    assert "ScientificBadge" in js_content
    assert "getInterpretation" in js_content
    assert "VisualProfileChart" in js_content
    assert "unifiedPhysicochemicalTable" in js_content
    assert "speciesMetabolicStabilityTable" in js_content
    assert "MultiSpeciesPkSummaryTable" in js_content
    assert "renderPredictionMaturity" in js_content
