from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_shared_unified_endpoint_comparison_has_all_four_states_and_columns():
    js = (ROOT / "frontend/static/app.js").read_text(encoding="utf-8")
    assert "function UnifiedEndpointComparison" in js
    assert "Experimental" in js and "Prediction" in js
    assert "Difference" in js and "Project Learning" in js
    assert "BOTH" in js and "PREDICTION_ONLY" in js and "EXPERIMENTAL_ONLY" in js
    assert "No experimental value yet" in js
    assert "No matching prediction endpoint" in js
    assert "project_adapted_prediction" in js


def test_unified_endpoint_rows_are_deduplicated_and_candidates_remain_importable():
    js = (ROOT / "frontend/static/app.js").read_text(encoding="utf-8")
    assert "display_sources||[row.source||'External']" in js
    assert "Individual observations preserved" in js
    assert "External candidate · not imported" in js
    assert "Raw source values remain preserved" in js
    assert "no aggregate value is used for adaptation" in js


def test_project_learning_panel_exposes_endpoint_specific_learning_and_activation():
    js = (ROOT / "frontend/static/app.js").read_text(encoding="utf-8")
    assert "function projectLearningPanel" in js
    assert "Independent compounds" in js
    assert "Effective N" in js
    assert "Activate" in js
    assert "No qualified project endpoint pairs yet." in js
