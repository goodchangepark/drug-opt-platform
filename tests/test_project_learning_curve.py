import copy
import json
from pathlib import Path

from backend.project_adaptation_v2 import fit_project_adapter
from backend.project_learning_curve import (
    build_disposable_learning_demo,
    build_learning_curve,
    build_synthetic_validation_dataset,
    maturity_policy,
)


GLOBAL = {"model_a": .60, "model_b": .25, "model_c": .15}


def test_learning_curve_is_endpoint_specific_and_base_until_five():
    events = build_synthetic_validation_dataset("Solubility", 12)
    result = build_learning_curve("Solubility", events, GLOBAL)
    assert result["unique_experiments"] == 12
    points = {point["n"]: point for point in result["aggregate"]}
    assert points[0]["adapted_mae"] is None
    assert points[3]["adapted_mae"] is None
    assert points[5]["adapted_mae"] is not None
    assert points[5]["validation_decisions"] == ["LIGHT_PROJECT_ADAPTATION_CANDIDATE"]
    assert result["engine_hash"] == "12757ab197b5a70d8ea1754678d9a342ab0b6ea0d82f2896bebb767d686bbdeb"


def test_every_holdout_is_excluded_from_training_and_repeated_orderings_exist():
    events = build_synthetic_validation_dataset("PPB", 12)
    result = build_learning_curve("PPB", events, GLOBAL)
    assert len(result["ordering_seeds"]) == 20
    for point in result["primary_ordering"]:
        assert not set(point["training_compounds"]).intersection(point["holdout_compounds"])
        for observation in point["holdout_predictions"]:
            assert observation["compound_version_id"] in point["holdout_compounds"]


def test_structure_and_series_splits_are_reported_without_changing_frozen_predictions():
    events = build_synthetic_validation_dataset("Caco2", 12)
    before = copy.deepcopy([event.frozen_predictions for event in events])
    result = build_learning_curve("Caco2", events, GLOBAL)
    assert {row["mode"] for row in result["split_evaluations"]} == {"random", "scaffold", "series"}
    assert [event.frozen_predictions for event in events] == before
    assert set(result["similarity_analysis"]).issubset({"high", "medium", "low"})


def test_weights_are_a_simplex_and_equal_models_do_not_claim_learning():
    events = build_synthetic_validation_dataset("HLM", 12)
    result = fit_project_adapter("HLM", events[:5], GLOBAL)
    assert min(result.project_weights.values()) >= 0
    assert abs(sum(result.project_weights.values()) - 1) < 1e-12
    equal = [event.__class__(event.evidence_id, event.compound_version_id, event.smiles, "HLM", event.value,
                             {"model_a": event.value, "model_b": event.value, "model_c": event.value}) for event in events[:5]]
    assert fit_project_adapter("HLM", equal, GLOBAL).activation_decision == "BASE_RETAINED"


def test_disposable_new_compound_snapshot_and_reveal_are_separate():
    demo = build_disposable_learning_demo()
    assert demo["experiment_initially_absent"] is True
    assert demo["adapter_activation_decision"] == "ACTIVATED"
    assert demo["project_prediction"] is not None
    assert demo["project_error"] < demo["base_error"]
    assert demo["snapshot_immutable"] is True
    assert demo["prediction_snapshot"]["experiment_used_for_prediction"] is False
    assert demo["prediction_snapshot"]["adapter_version"]
    assert demo["prediction_snapshot"]["training_compounds"] == ["1", "2", "3", "4", "5"]


def test_maturity_policy_requires_validation_and_activation_not_n_alone():
    policy = maturity_policy()
    assert "explicit activation" in policy["rules"]["level_2"]["gate"]
    assert "not N alone" in policy["rules"]["level_3"]["gate"]
    assert "never N alone" in policy["rules"]["level_5"]["gate"]


def test_validation_artifacts_are_versioned_and_include_endpoint_snapshots():
    root = Path(__file__).parents[1]
    curve = json.loads((root / "validation/project_learning_curve_v3_6.json").read_text())
    policy = json.loads((root / "validation/project_adaptation_maturity_policy_v1.json").read_text())
    assert set(curve["endpoints"]) == {"PPB", "Solubility", "Caco2", "HLM"}
    assert all(item["independent_compounds"] == 12 for item in curve["endpoints"].values())
    assert all(len(item["snapshots"]) == 12 for item in curve["endpoints"].values())
    assert policy["calibration_basis"]["validation_artifact"] == "project_learning_curve_v3_6.json"
