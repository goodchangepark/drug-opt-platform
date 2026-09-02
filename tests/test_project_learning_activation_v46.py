"""v4.6 activation/review contract checks."""
import json
from pathlib import Path


def test_v46_activation_artifacts_preserve_explicit_import_and_antileakage():
    root = Path(__file__).resolve().parents[1] / "validation"
    activation = json.loads((root / "egfr_evidence_activation_v4_6.json").read_text())
    flow = json.loads((root / "project_learning_import_flow_v4_6.json").read_text())
    pairing = json.loads((root / "project_learning_pair_generation_v4_6.json").read_text())
    assert activation["automatic_import"] is False
    assert flow["import_batch_persisted"] is True
    assert flow["automatic_import"] is False
    assert pairing["post_experiment_prediction_regenerated"] is False
    assert pairing["historical_prediction_mutated"] is False


def test_v46_candidate_evaluation_requires_threshold_and_explicit_activation():
    root = Path(__file__).resolve().parents[1] / "validation"
    evaluation = json.loads((root / "project_learning_candidate_evaluation_v4_6.json").read_text())
    maturity = json.loads((root / "project_learning_maturity_transition_v4_6.json").read_text())
    assert evaluation["threshold"] == 5
    assert evaluation["activation"] == "explicit user action required"
    assert maturity["searched_candidates_promote_maturity"] is False
    assert maturity["n_only_promotion"] is False
