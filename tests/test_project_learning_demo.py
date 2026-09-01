from backend.project_learning_demo import run_synthetic_learning_demo


def test_synthetic_six_compound_learning_progression_is_explicit_and_endpoint_scoped():
    demo = run_synthetic_learning_demo()
    assert len(demo["compounds"]) == 6
    assert [row["status"] for row in demo["tiers"][:4]] == ["BASE_ONLY"] * 4
    assert demo["tiers"][4]["status"] == "LIGHT_PROJECT_ADAPTATION"
    assert demo["candidate"]["requires_explicit_activation"] is True
    assert demo["candidate"]["activation_decision"] == "ACTIVATED"
    assert demo["compound_six_before_experiment"]["experimental"] is None
    assert demo["compound_six_before_experiment"]["maturity"]["level"] == 2
    assert demo["compound_six_before_experiment"]["project_prediction"] != demo["compound_six_before_experiment"]["base_prediction"]
    assert demo["compound_six_after_experiment"]["project_error"] < demo["compound_six_after_experiment"]["base_error"]


def test_synthetic_fixture_does_not_claim_real_project_activation():
    demo = run_synthetic_learning_demo()
    assert demo["candidate"]["requires_explicit_activation"] is True
    assert demo["compound_six_before_experiment"]["prediction_source"] == "Project-adapted Prediction"

