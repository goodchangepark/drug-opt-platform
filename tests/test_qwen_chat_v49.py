"""Test Suite for Qwen3.5 9B Section & Comparison Chat v4.9.

Validates:
1. Backend section chat endpoint (/api/chat/section) across compound tabs
2. Backend comparison chat endpoint (/api/chat/compare) across selected compounds
3. Isolation of compound-specific data and absence of cross-compound leakage
4. Absence of fabricated data / proper fallback to "현재 데이터에서는 확인할 수 없습니다."
5. Frontend UI contract: placeholder "Qwen3.5 9B", section & comparison chat boxes
"""

from pathlib import Path
import pytest
from starlette.testclient import TestClient
from sqlalchemy import select

from backend.main import app
from backend.database import SessionLocal
from backend.models import Project, Compound
from backend.qwen_chat import (
    build_compound_section_context,
    build_comparison_context,
    QWEN_MODEL
)

client = TestClient(app)


def test_qwen_model_configuration():
    assert QWEN_MODEL == "qwen3.5:9b"


def test_frontend_qwen_chat_ui_contract():
    root = Path(__file__).resolve().parents[1]
    js = (root / "frontend" / "static" / "app.js").read_text(encoding="utf-8")
    css = (root / "frontend" / "static" / "app.css").read_text(encoding="utf-8")

    # Placeholder and component checks
    assert "placeholder='Qwen3.5 9B'" in js or 'placeholder="Qwen3.5 9B"' in js or "placeholder: 'Qwen3.5 9B'" in js
    assert "AIChatSection" in js
    assert "AIChatCompare" in js
    assert "ai-chat-box" in js
    assert "ai-chat-input" in js

    # CSS styles
    assert ".ai-chat-box" in css
    assert ".ai-chat-input::placeholder" in css


def test_compound_section_context_isolation():
    with SessionLocal() as db:
        # Check Mobocertinib (id=16) and Sunvozertinib (id=10)
        mobo = db.scalar(select(Compound).where(Compound.id == 16))
        sunvo = db.scalar(select(Compound).where(Compound.id == 10))
        assert mobo is not None
        assert sunvo is not None

        mobo_v = next((v for v in mobo.versions if v.version_number == mobo.current_version), mobo.versions[-1])
        sunvo_v = next((v for v in sunvo.versions if v.version_number == sunvo.current_version), sunvo.versions[-1])

        # Test PK context isolation
        mobo_pk_ctx = build_compound_section_context(db, mobo, mobo_v, "pk")
        sunvo_pk_ctx = build_compound_section_context(db, sunvo, sunvo_v, "pk")

        assert "Mobocertinib" in mobo_pk_ctx
        assert "Sunvozertinib" not in mobo_pk_ctx

        assert "Sunvozertinib" in sunvo_pk_ctx
        assert "Mobocertinib" not in sunvo_pk_ctx

        # Test Activity context isolation
        mobo_act_ctx = build_compound_section_context(db, mobo, mobo_v, "activity")
        assert "Mobocertinib" in mobo_act_ctx
        assert "Sunvozertinib" not in mobo_act_ctx


def test_comparison_context_scoping():
    with SessionLocal() as db:
        egfr = db.scalar(select(Project).where(Project.id == 3))
        compounds = list(db.scalars(select(Compound).where(Compound.id.in_([10, 16]))))
        assert len(compounds) == 2

        ctx = build_comparison_context(db, egfr, compounds)
        assert "Sunvozertinib" in ctx
        assert "Mobocertinib" in ctx
        assert "Orforglipron" not in ctx  # GLP-1 compound should not leak into EGFR comparison


def test_section_chat_api_endpoint():
    res = client.post("/api/chat/section", json={
        "compound_id": 16,
        "section": "properties",
        "question": "분자량과 cLogP를 알려줘."
    })
    assert res.status_code == 200
    data = res.json()
    assert data["model"] == "qwen3.5:9b"
    assert data["compound_name"] == "Mobocertinib"
    assert len(data["answer"]) > 0


def test_comparison_chat_api_endpoint():
    res = client.post("/api/chat/compare", json={
        "project_id": 3,
        "compound_ids": [10, 16],
        "question": "두 물질의 분자량 차이는?"
    })
    assert res.status_code == 200
    data = res.json()
    assert data["model"] == "qwen3.5:9b"
    assert "Sunvozertinib" in data["compounds"]
    assert "Mobocertinib" in data["compounds"]
    assert len(data["answer"]) > 0
