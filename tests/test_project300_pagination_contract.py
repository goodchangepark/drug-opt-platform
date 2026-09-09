"""Regression coverage for the bounded Global Reference Library list path."""

from pathlib import Path

from backend.database import SessionLocal
from backend.main import get_compound_scientific_summary, get_project


def test_reference_library_api_returns_bounded_lightweight_page():
    with SessionLocal() as db:
        response = get_project(300, db, page=1, page_size=50, search="")

    assert response["list_mode"] == "PAGINATED_SUMMARY"
    assert response["pagination"] == {"page": 1, "page_size": 50, "total": 1000, "total_pages": 20}
    assert len(response["compounds"]) == 50
    assert all(row["version"]["svg"] == "" for row in response["compounds"] if row["version"])
    assert all(row["versions"] == [] for row in response["compounds"])


def test_reference_library_search_is_server_side_and_frontend_uses_bounded_query():
    with SessionLocal() as db:
        response = get_project(300, db, page=1, page_size=50, search="Warfarin")

    assert response["pagination"]["total"] == 1
    assert [row["name"] for row in response["compounds"]] == ["Warfarin"]
    js = (Path(__file__).resolve().parents[1] / "frontend/static/app.js").read_text(encoding="utf-8")
    assert "page_size=50" in js
    assert "project-reference-search" in js
    assert "Server-paginated reference summaries" in js


def test_compound_open_uses_bounded_local_summary_and_lazy_tabs():
    with SessionLocal() as db:
        result = get_project(300, db, page=1, page_size=1, search="Warfarin")
        row_id = result["compounds"][0]["row_id"]
        summary = get_compound_scientific_summary(row_id, db)

    assert summary["payload_scope"] == "COMPOUND_CORE_SUMMARY"
    assert summary["name"] == "Warfarin"
    assert summary["prediction_engine"]["engine_id"] == "drugopt-prediction-engine-v3@3.3.3"
    assert summary["version"]["inchikey"]
    assert "prediction_history" not in summary
    assert "project_learning" not in summary
    assert "external_experimental_evidence" not in summary

    js = (Path(__file__).resolve().parents[1] / "frontend/static/app.js").read_text(encoding="utf-8")
    assert "api.get('/compounds/'+rowId+'/summary')" in js
    assert "['properties','activity','admet','metabolism','pk','evidence','history'].includes(detailTab)" in js
    assert "'/compound-versions/'+versionId+'/scientific-tabs/'+tab" in js
    assert "detailTab==='admet'" in js
    assert "detailTab==='metabolism'" in js
