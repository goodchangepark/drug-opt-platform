from backend.prediction_maturity import maturity_for_adapter
from pathlib import Path

def get(status, n, decision="ACTIVATED", **kw):
    return maturity_for_adapter(status=status, effective_n=n, activation_decision=decision, **kw)

def test_five_levels_require_validated_endpoint_specific_state():
    assert get("BASE_ONLY", 4).level == 1
    assert get("LIGHT_PROJECT_ADAPTATION", 5).level == 2
    assert get("REGULARIZED_PROJECT_ENSEMBLE", 10).level == 3
    assert get("LOCAL_SERIES_ADAPTATION", 20, representative_series=True).level == 4
    assert get("LOCAL_SERIES_ADAPTATION", 40, representative_series=True, stable_history_count=3).level == 5

def test_n_alone_or_related_evidence_cannot_promote():
    assert get("REGULARIZED_PROJECT_ENSEMBLE", 20, "BASE_RETAINED").level == 1
    assert get("LOCAL_SERIES_ADAPTATION", 50, compatible_evidence_only=False, representative_series=True, stable_history_count=9).level == 1

def test_accessible_stars_are_stable_metadata():
    item = get("REGULARIZED_PROJECT_ENSEMBLE", 10)
    assert item.stars == "★★★☆☆" and "3 of 5" in item.to_dict()["aria_label"]

def test_frontend_uses_accessible_gold_and_muted_star_classes():
    root = Path(__file__).parents[1]
    js, css = (root / "frontend/static/app.js").read_text(), (root / "frontend/static/app.css").read_text()
    assert "maturity-stars" in js and "aria-label" in js
    assert "#F5B700" in css and "maturity-star-empty" in css
