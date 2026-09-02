"""Deterministic, conservative interpretation policies for scientific rows."""

SCIENTIFIC_INTERPRETATION_VERSION = "drugopt-scientific-interpretation-policy-v1"
AGREEMENT_POLICY_VERSION = "drugopt-experimental-prediction-agreement-v1"


def interpret_row(*, prediction_available: bool, direct: bool, difference_available: bool) -> dict:
    """Never infer clinical desirability or unvalidated accuracy."""
    return {
        "value_assessment": "CONTEXT_DEPENDENT",
        "value_assessment_color": "neutral",
        "agreement": ("NO_PREDICTION" if not prediction_available else
                       "NOT_NUMERICALLY_COMPARABLE" if not direct else
                       "NOT_CALIBRATED" if not difference_available else "NOT_CALIBRATED"),
        "agreement_color": "neutral",
        "interpretation_policy": SCIENTIFIC_INTERPRETATION_VERSION,
        "agreement_policy": AGREEMENT_POLICY_VERSION,
        "policy_explanation": "No endpoint-specific validated target/agreement threshold is applied.",
    }


def policy_report() -> dict:
    return {
        "version": SCIENTIFIC_INTERPRETATION_VERSION,
        "agreement_version": AGREEMENT_POLICY_VERSION,
        "default_value_assessment": "CONTEXT_DEPENDENT",
        "default_agreement": "NOT_CALIBRATED",
        "colors": {"IN_TARGET": "green", "BORDERLINE": "amber", "OUT_OF_TARGET": "red", "CONTEXT_DEPENDENT": "neutral", "UNKNOWN": "gray", "GOOD_AGREEMENT": "green", "MODERATE_DIFFERENCE": "amber", "LARGE_DIFFERENCE": "red", "NOT_CALIBRATED": "neutral"},
        "universal_thresholds": False,
    }
