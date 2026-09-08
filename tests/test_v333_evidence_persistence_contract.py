import pytest
from fastapi import HTTPException

from backend.database import SessionLocal
from backend.main import search_external_experimental_data


def test_legacy_transient_search_is_retired_in_favor_of_persisted_workflow():
    with SessionLocal() as db:
        with pytest.raises(HTTPException) as exc:
            search_external_experimental_data(1, db)
    assert exc.value.status_code == 410
    detail = exc.value.detail
    assert detail["status"] == "TRANSIENT_SEARCH_RETIRED"
    assert detail["replacement"].endswith("/experimental-harvest/preview")
    assert detail["confirmation_required"] is True


def test_normal_ui_uses_only_persisting_search_workflow():
    from pathlib import Path

    source = (Path(__file__).resolve().parents[1] / "frontend/static/app.js").read_text(encoding="utf-8")
    assert "'/compounds/'+detail.row_id+'/experimental-harvest/preview'" in source
    assert "/external-experimental/search" not in source
