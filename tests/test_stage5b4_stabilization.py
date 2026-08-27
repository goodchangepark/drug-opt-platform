"""Pre-Stage-5C stabilization contracts."""

import importlib.metadata
import json
import math
from pathlib import Path

from backend.admet_predictor import MODEL_SPECS, predict_endpoint
from backend.database import SessionLocal
from backend.main import app, dashboard_summary, health, help_registry
from backend.stabilization import classify_project


ROOT = Path(__file__).resolve().parents[1]


def test_help_registry_endpoint_uses_runtime_registries_and_versions():
    assert any(getattr(route, "path", None) == "/api/help/registry" for route in app.routes)
    with SessionLocal() as db:
        payload = help_registry(db)
    assert payload["application"]["current_stage"] == "5B-4"
    assert payload["application"]["version"] in ("0.6.0-stage5b4-stable", "0.6.1-stage5b4-ui", "0.6.2-stage5b4-ui", "0.6.3-stage5b4-ui")
    versions = {row["package"]: row["version"] for row in payload["package_inventory"]}
    assert versions["RDKit"] == importlib.metadata.version("rdkit")
    assert versions["Chemprop"] == importlib.metadata.version("chemprop")
    assert versions["PyTorch"] == importlib.metadata.version("torch")
    assert len(payload["models"]) >= len(MODEL_SPECS)
    assert payload["pk_method_registry"]
    assert payload["source"].startswith("RUNTIME_PACKAGE_INVENTORY")


def test_help_and_dashboard_share_capability_source_and_pk_is_ready():
    with SessionLocal() as db:
        help_payload = help_registry(db)
        dashboard = dashboard_summary(db)
    assert help_payload["capability_summary"] == dashboard["capability_summary"]
    groups = {row["title"]: row for row in dashboard["capability_summary"]["groups"]}
    assert groups["PK / DMPK"]["status"] == "READY"
    assert groups["CYP & Transporters"]["status"] == "PARTIAL"
    assert groups["Safety / Toxicology"]["status"] == "PARTIAL"
    assert health()["step"] == "5B-4"


def test_help_frontend_is_registry_driven_and_optimization_navigation_is_guarded():
    source = (ROOT / "frontend/static/app.js").read_text(encoding="utf-8")
    assert "api.get('/help/registry')" in source
    assert "helpRegistry.models" in source
    assert "helpRegistry.package_inventory" in source
    assert "data-testid':'help-registry'" in source
    assert "Current Limitations" in source and "PK / DMPK" in source
    assert "Select a project to begin optimization." in source
    assert "selectedProject&&compounds.length===0" in source
    assert "selectedProject&&compounds.length>0&&!selectedCompound" in source


def test_project_cleanup_classifier_deletes_only_positive_test_markers():
    confirmed = classify_project({"project_name": "Stage 5B-3 Browser Acceptance 20260827-010539",
                                  "target": "Translational PK Acceptance", "compound_count": 1})
    genuine = classify_project({"project_name": "NME", "target": "GLP-1R",
                                "description": "best-in-class", "compound_count": 2})
    ambiguous = classify_project({"project_name": "Discovery", "target": "", "compound_count": 0})
    assert confirmed[0] == "CONFIRMED_TEST"
    assert genuine[0] == "KEEP"
    assert ambiguous[0] == "AMBIGUOUS"


def test_installed_model_runtime_audit_and_live_inference_are_operational():
    audit = json.loads((ROOT / "validation/stage5b4_stabilization_model_audit.json").read_text())
    assert audit["summary"]["status"] == "PASS"
    assert audit["summary"]["passed"] == len(MODEL_SPECS) == 18
    assert all(row["registry_entry"] and row["assets_available"] and row["finite_prediction"]
               and row["endpoint_mapping"] == "PASS" for row in audit["models"])
    live = predict_endpoint("CC(=O)Oc1ccccc1C(=O)O", "Solubility")
    assert live["status"] == "COMPLETE"
    assert live["unit"] == MODEL_SPECS["Solubility"]["unit"]
    assert math.isfinite(live["predicted_value"])


def test_cleanup_artifact_preserves_every_nonconfirmed_project():
    cleanup = json.loads((ROOT / "validation/test_projects_cleanup.json").read_text())
    assert cleanup["deleted_count"] == cleanup["summary"]["CONFIRMED_TEST"]
    assert all(row["deletion_status"] == "DELETED" for row in cleanup["projects"]
               if row["classification"] == "CONFIRMED_TEST")
    assert all(row["deletion_status"] == "PRESERVED" for row in cleanup["projects"]
               if row["classification"] in {"AMBIGUOUS", "KEEP"})
