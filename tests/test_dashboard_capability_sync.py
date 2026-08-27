"""Dashboard capability status synchronization tests."""

from pathlib import Path

from backend.capabilities import build_capability_summary
from backend.database import SessionLocal
from backend.main import dashboard_summary


ROOT = Path(__file__).resolve().parents[1]


def _groups(payload):
    return {group["title"]: group for group in payload["capability_summary"]["groups"]}


def _items(group):
    return {item["label"]: item for item in group["items"]}


def _dashboard_payload():
    db = SessionLocal()
    try:
        return dashboard_summary(db)
    finally:
        db.close()


def test_dashboard_capabilities_follow_live_backend_registries():
    payload = _dashboard_payload()
    groups = _groups(payload)
    registry = {row["endpoint"]: row for row in payload["model_registry"]}

    adme = _items(groups["ADME"])
    for label in ("Solubility", "Caco-2", "PPB", "fu", "HLM", "RLM", "MLM"):
        assert adme[label]["availability"] == "READY"
    assert adme["HLM"]["confidence"] == "LOW"
    assert adme["HLM"]["availability"] != adme["HLM"]["confidence"]

    cyp = _items(groups["CYP & Transporters"])
    assert groups["CYP & Transporters"]["status"] == "PARTIAL"
    for label in (
        "CYP1A2 inhibitor", "CYP2C9 inhibitor", "CYP2C19 inhibitor", "CYP2D6 inhibitor",
        "CYP3A4 inhibitor", "CYP2C9 substrate", "CYP2D6 substrate", "CYP3A4 substrate", "P-gp inhibitor",
    ):
        assert cyp[label]["availability"] == "READY"
        assert cyp[label]["availability"] == registry[cyp[label]["endpoint"]]["availability"]
    for label in ("BCRP inhibitor", "BSEP inhibitor", "OATP1B1 inhibitor", "OCT1 inhibitor", "MATE1 inhibitor"):
        assert cyp[label]["availability"] == "MODEL_UNAVAILABLE"

    safety = _items(groups["Safety / Toxicology"])
    assert groups["Safety / Toxicology"]["status"] == "PARTIAL"
    for label in ("hERG", "Ames", "DILI", "Structural Alerts"):
        assert safety[label]["availability"] == "READY"

    pk = _items(groups["PK / DMPK"])
    assert groups["PK / DMPK"]["status"] == "READY"
    expected_pk = {
        "Experimental PK Data Management", "NCA", "IVIVE / Hepatic Clearance", "Vd / Absorption Foundation",
        "IV Simulation", "PO / SC / IP Simulation", "Cross-Species Scaling", "Human Translational PK",
        "Prospective Prediction Freeze", "Retrospective Validation",
    }
    assert expected_pk == set(pk)
    assert all(item["availability"] == "READY" for item in pk.values())


def test_availability_confidence_and_conformal_status_are_independent():
    payload = _dashboard_payload()
    registry = {row["endpoint"]: row for row in payload["model_registry"]}
    assert registry["HLM intrinsic clearance"]["availability"] == "READY"
    assert registry["HLM intrinsic clearance"]["confidence"] == "LOW"
    assert registry["RLM intrinsic clearance"]["conformal_status"] == "CONFORMAL_UNAVAILABLE"
    assert registry["RLM intrinsic clearance"]["availability"] == "READY"


def test_missing_backend_route_marks_feature_unavailable():
    payload = build_capability_summary([], route_paths=())
    groups = _groups({"capability_summary": payload})
    assert groups["PK / DMPK"]["status"] == "MODEL_UNAVAILABLE"
    assert all(item["availability"] == "MODEL_UNAVAILABLE" for item in groups["PK / DMPK"]["items"])


def test_dashboard_frontend_consumes_summary_and_busts_stale_assets():
    source = (ROOT / "frontend/static/app.js").read_text(encoding="utf-8")
    index = (ROOT / "frontend/static/index.html").read_text(encoding="utf-8")
    dashboard = source[source.index("function MainDashboard()") : source.index("function ProjectWorkspace()")]
    assert "dashboard?.capability_summary?.groups" in dashboard
    assert "const modules=[" not in dashboard
    assert "Confidence: " in dashboard and "Conformal: " in dashboard
    assert "PK / DMPK is planned" not in source
    assert "stage5b4-stable-1" in index
