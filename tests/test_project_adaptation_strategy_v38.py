from backend.project_adaptation_strategy import (
    BASE_ONLY, MULTI_MODEL_ENSEMBLE_ADAPTATION, SINGLE_MODEL_RESIDUAL_CALIBRATION,
    fit_project_adaptation_strategy, resolve_adaptation_strategy,
)
from backend.project_adaptation_v2 import QualifiedEvidencePair


def pairs(values, predictions):
    return [QualifiedEvidencePair(str(i), i, "CCO", "Solubility", value, predictions[i - 1]) for i, value in enumerate(values, 1)]


def test_strategy_resolver_distinguishes_single_and_multi_model_endpoints():
    assert resolve_adaptation_strategy(0).strategy_type == BASE_ONLY
    assert resolve_adaptation_strategy(1, "Solubility").strategy_type == SINGLE_MODEL_RESIDUAL_CALIBRATION
    assert resolve_adaptation_strategy(3, "Solubility").strategy_type == MULTI_MODEL_ENSEMBLE_ADAPTATION
    assert resolve_adaptation_strategy(1, "Ames", "binary_classification").strategy_type != SINGLE_MODEL_RESIDUAL_CALIBRATION


def test_single_model_systematic_bias_can_validate_conservative_calibration():
    evidence = pairs([1.30] * 6, [{"model": 1.0}] * 6)
    result = fit_project_adaptation_strategy("Solubility", evidence, {"model": 1.0})
    assert result.strategy_type == SINGLE_MODEL_RESIDUAL_CALIBRATION
    assert result.activation_decision == "ACTIVATED"
    assert result.shrinkage_factor < 0.55
    assert result.calibration_adjustment > 0
    assert result.adapted_validation_error < result.base_validation_error


def test_single_model_random_signed_residuals_retain_base():
    evidence = pairs([1.2, 0.8, 1.2, 0.8, 1.2, 0.8], [{"model": 1.0}] * 6)
    result = fit_project_adaptation_strategy("Solubility", evidence, {"model": 1.0})
    assert result.strategy_type == SINGLE_MODEL_RESIDUAL_CALIBRATION
    assert result.activation_decision == "BASE_RETAINED"
    assert result.status == "BASE_RETAINED_NO_STABLE_PROJECT_BIAS"
    assert result.calibration_adjustment == 0


def test_multi_model_strategy_keeps_weights_nonnegative_and_normalized():
    evidence = [
        QualifiedEvidencePair(str(i), i, "CCO", "Solubility", 1.0, {"a": 1.0, "b": 1.1, "c": 1.3})
        for i in range(1, 7)
    ]
    result = fit_project_adaptation_strategy("Solubility", evidence, {"a": 1 / 3, "b": 1 / 3, "c": 1 / 3})
    assert result.strategy_type == MULTI_MODEL_ENSEMBLE_ADAPTATION
    assert all(value >= 0 for value in result.project_weights.values())
    assert abs(sum(result.project_weights.values()) - 1) < 1e-9
