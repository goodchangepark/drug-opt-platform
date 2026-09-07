import pytest

from backend.project_learning_contract import (
    adapter_capability, freeze_prospective_prediction, loco_training,
    recommend_prediction,
)


class Event:
    def __init__(self, compound_version_id):
        self.compound_version_id = compound_version_id


def test_loco_removes_every_target_compound_observation():
    events = [Event(1), Event(1), Event(2), Event(3)]
    assert [x.compound_version_id for x in loco_training(events, 1)] == [2, 3]


def test_recommendation_never_uses_unvalidated_project_value():
    result = recommend_prediction(2.0, 1.0, status="PROJECT_RESEARCH_CANDIDATE", global_error=2.0, project_error=1.0)
    assert result["source"] == "GLOBAL" and result["value"] == 2.0
    result = recommend_prediction(2.0, 1.0, status="PROJECT_VALIDATED_CANDIDATE", global_error=2.0, project_error=1.0)
    assert result["source"] == "PROJECT"


def test_prospective_freeze_rejects_known_experiment_and_is_immutable():
    with pytest.raises(ValueError, match="PROSPECTIVE_FREEZE"):
        freeze_prospective_prediction(project_id=3, compound_id=1, endpoint="HLM", global_prediction=1, global_uncertainty=.1, global_ad="IN_DOMAIN", model_version="v1", timestamp="t", experiment_known=True)
    snap = freeze_prospective_prediction(project_id=3, compound_id=1, endpoint="HLM", global_prediction=1, global_uncertainty=.1, global_ad="IN_DOMAIN", model_version="v1", timestamp="t")
    assert snap["immutable"] is True and snap["prospective_status"] == "PROSPECTIVE_PROJECT_PREDICTION"


def test_weak_global_endpoints_cannot_receive_false_precision_adapter():
    assert adapter_capability("VDSS")["eligible"] is False
    assert adapter_capability("Solubility")["eligible"] is True
