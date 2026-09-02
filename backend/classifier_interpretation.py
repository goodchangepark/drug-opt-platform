"""Classifier Interpretation Policy & Contract.

Policy Version: drugopt-classifier-interpretation-v1

Defines researcher-meaningful semantics for binary classification predictions
(CYP3A4 inhibitor, P-gp inhibitor, hERG, Ames, DILI).

Key Rules:
1. Model score is a classifier output (probability / logit proxy), NOT % inhibition or potency.
2. If probability calibration is not established, display "Model score: X.XXX (Calibration: Not established)".
3. Do NOT subtract experimental IC50/Ki from classifier score. Difference is "—".
4. If a deterministic continuous-to-binary mapping is not proven by model training contract,
   report "Related measurement — quantitative agreement not calibrated".
5. Quantitative potency prediction gaps are explicitly noted as QUANTITATIVE_MODEL_GAP.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

CLASSIFIER_INTERPRETATION_POLICY_VERSION = "drugopt-classifier-interpretation-v1"

CLASSIFIER_REGISTRY: dict[str, dict[str, Any]] = {
    "CYP3A4_INHIBITION": {
        "endpoint_name": "CYP3A4 inhibition",
        "model_type": "BINARY_CLASSIFICATION",
        "positive_class": "Inhibitor",
        "negative_class": "Non-inhibitor",
        "decision_threshold": 0.50,
        "calibrated_probability": False,
        "calibration_method": "NONE",
        "training_target_definition": "Binary CYP3A4 inhibition (positive if IC50 <= 10 µM in in-vitro assay)",
        "score_label": "Inhibitor score",
        "unit": "score (0-1)",
        "quantitative_model_available": False,
        "quantitative_gap_code": "QUANTITATIVE_MODEL_GAP",
    },
    "CYP1A2_INHIBITION": {
        "endpoint_name": "CYP1A2 inhibition",
        "model_type": "BINARY_CLASSIFICATION",
        "positive_class": "Inhibitor",
        "negative_class": "Non-inhibitor",
        "decision_threshold": 0.50,
        "calibrated_probability": False,
        "calibration_method": "NONE",
        "training_target_definition": "Binary CYP1A2 inhibition",
        "score_label": "Inhibitor score",
        "unit": "score (0-1)",
        "quantitative_model_available": False,
        "quantitative_gap_code": "QUANTITATIVE_MODEL_GAP",
    },
    "CYP2C19_INHIBITION": {
        "endpoint_name": "CYP2C19 inhibition",
        "model_type": "BINARY_CLASSIFICATION",
        "positive_class": "Inhibitor",
        "negative_class": "Non-inhibitor",
        "decision_threshold": 0.50,
        "calibrated_probability": False,
        "calibration_method": "NONE",
        "training_target_definition": "Binary CYP2C19 inhibition",
        "score_label": "Inhibitor score",
        "unit": "score (0-1)",
        "quantitative_model_available": False,
        "quantitative_gap_code": "QUANTITATIVE_MODEL_GAP",
    },
    "CYP2C9_INHIBITION": {
        "endpoint_name": "CYP2C9 inhibition",
        "model_type": "BINARY_CLASSIFICATION",
        "positive_class": "Inhibitor",
        "negative_class": "Non-inhibitor",
        "decision_threshold": 0.50,
        "calibrated_probability": False,
        "calibration_method": "NONE",
        "training_target_definition": "Binary CYP2C9 inhibition",
        "score_label": "Inhibitor score",
        "unit": "score (0-1)",
        "quantitative_model_available": False,
        "quantitative_gap_code": "QUANTITATIVE_MODEL_GAP",
    },
    "CYP2D6_INHIBITION": {
        "endpoint_name": "CYP2D6 inhibition",
        "model_type": "BINARY_CLASSIFICATION",
        "positive_class": "Inhibitor",
        "negative_class": "Non-inhibitor",
        "decision_threshold": 0.50,
        "calibrated_probability": False,
        "calibration_method": "NONE",
        "training_target_definition": "Binary CYP2D6 inhibition",
        "score_label": "Inhibitor score",
        "unit": "score (0-1)",
        "quantitative_model_available": False,
        "quantitative_gap_code": "QUANTITATIVE_MODEL_GAP",
    },
    "CYP2C9_SUBSTRATE": {
        "endpoint_name": "CYP2C9 substrate",
        "model_type": "BINARY_CLASSIFICATION",
        "positive_class": "Substrate",
        "negative_class": "Non-substrate",
        "decision_threshold": 0.50,
        "calibrated_probability": False,
        "calibration_method": "NONE",
        "training_target_definition": "Binary CYP2C9 substrate",
        "score_label": "Substrate score",
        "unit": "score (0-1)",
        "quantitative_model_available": False,
        "quantitative_gap_code": "NOT_APPLICABLE",
    },
    "CYP2D6_SUBSTRATE": {
        "endpoint_name": "CYP2D6 substrate",
        "model_type": "BINARY_CLASSIFICATION",
        "positive_class": "Substrate",
        "negative_class": "Non-substrate",
        "decision_threshold": 0.50,
        "calibrated_probability": False,
        "calibration_method": "NONE",
        "training_target_definition": "Binary CYP2D6 substrate",
        "score_label": "Substrate score",
        "unit": "score (0-1)",
        "quantitative_model_available": False,
        "quantitative_gap_code": "NOT_APPLICABLE",
    },
    "CYP3A4_SUBSTRATE": {
        "endpoint_name": "CYP3A4 substrate",
        "model_type": "BINARY_CLASSIFICATION",
        "positive_class": "Substrate",
        "negative_class": "Non-substrate",
        "decision_threshold": 0.50,
        "calibrated_probability": False,
        "calibration_method": "NONE",
        "training_target_definition": "Binary CYP3A4 substrate",
        "score_label": "Substrate score",
        "unit": "score (0-1)",
        "quantitative_model_available": False,
        "quantitative_gap_code": "NOT_APPLICABLE",
    },
    "PGP_INHIBITION": {
        "endpoint_name": "P-gp inhibition",
        "model_type": "BINARY_CLASSIFICATION",
        "positive_class": "Inhibitor",
        "negative_class": "Non-inhibitor",
        "decision_threshold": 0.50,
        "calibrated_probability": False,
        "calibration_method": "NONE",
        "training_target_definition": "Binary P-gp / MDR1 inhibition in cell-based transport assays",
        "score_label": "Inhibitor score",
        "unit": "score (0-1)",
        "quantitative_model_available": False,
        "quantitative_gap_code": "QUANTITATIVE_MODEL_GAP",
    },
    "HERG_LIABILITY": {
        "endpoint_name": "hERG liability",
        "model_type": "BINARY_CLASSIFICATION",
        "positive_class": "hERG Blocker",
        "negative_class": "Non-blocker",
        "decision_threshold": 0.50,
        "calibrated_probability": False,
        "calibration_method": "NONE",
        "training_target_definition": "Binary hERG inhibition (IC50 <= 10 µM / 50% block at test conc)",
        "score_label": "Blocker score",
        "unit": "score (0-1)",
        "quantitative_model_available": False,
        "quantitative_gap_code": "QUANTITATIVE_MODEL_GAP",
    },
    "AMES_MUTAGENICITY": {
        "endpoint_name": "Ames mutagenicity",
        "model_type": "BINARY_CLASSIFICATION",
        "positive_class": "Mutagenic",
        "negative_class": "Non-mutagenic",
        "decision_threshold": 0.50,
        "calibrated_probability": False,
        "calibration_method": "NONE",
        "training_target_definition": "Salmonella typhimurium reverse mutation assay (Ames positive/negative)",
        "score_label": "Mutagenicity score",
        "unit": "score (0-1)",
        "quantitative_model_available": False,
        "quantitative_gap_code": "NOT_APPLICABLE",
    },
    "DILI_LIABILITY": {
        "endpoint_name": "DILI liability",
        "model_type": "BINARY_CLASSIFICATION",
        "positive_class": "DILI Positive",
        "negative_class": "DILI Negative",
        "decision_threshold": 0.50,
        "calibrated_probability": False,
        "calibration_method": "NONE",
        "training_target_definition": "Drug-induced liver injury clinical annotation",
        "score_label": "DILI score",
        "unit": "score (0-1)",
        "quantitative_model_available": False,
        "quantitative_gap_code": "NOT_APPLICABLE",
    },
}


def interpret_classifier_prediction(endpoint_id: str, raw_score: float | None) -> dict[str, Any]:
    """Return scientific interpretation for a classifier score."""
    key = str(endpoint_id or "").upper()
    spec = CLASSIFIER_REGISTRY.get(key)
    if not spec:
        # Check aliases
        for k, v in CLASSIFIER_REGISTRY.items():
            if k in key or v["endpoint_name"].upper() in key:
                spec = v
                break

    if not spec:
        return {
            "is_classifier": False,
            "display_prediction": f"{raw_score:.3f}" if raw_score is not None else "—",
            "score": raw_score,
        }

    if raw_score is None:
        return {
            "is_classifier": True,
            "endpoint_name": spec["endpoint_name"],
            "prediction_class": "UNAVAILABLE",
            "display_text": "Unavailable",
            "model_score": None,
            "threshold": spec["decision_threshold"],
            "calibrated": spec["calibrated_probability"],
            "quantitative_gap": spec.get("quantitative_gap_code"),
        }

    is_positive = raw_score >= spec["decision_threshold"]
    predicted_class = spec["positive_class"] if is_positive else spec["negative_class"]

    if spec["calibrated_probability"]:
        display_text = f"{predicted_class} (Probability: {raw_score * 100:.1f}%, Threshold: {spec['decision_threshold']:.2f})"
    else:
        display_text = f"{predicted_class} (Score: {raw_score:.3f}, Threshold: {spec['decision_threshold']:.2f}, Uncalibrated)"

    return {
        "is_classifier": True,
        "policy_version": CLASSIFIER_INTERPRETATION_POLICY_VERSION,
        "endpoint_name": spec["endpoint_name"],
        "predicted_class": predicted_class,
        "display_text": display_text,
        "raw_score": raw_score,
        "decision_threshold": spec["decision_threshold"],
        "calibrated_probability": spec["calibrated_probability"],
        "calibration_status": "CALIBRATED" if spec["calibrated_probability"] else "NOT_ESTABLISHED",
        "training_target_definition": spec["training_target_definition"],
        "score_label": spec["score_label"],
        "quantitative_gap": spec.get("quantitative_gap_code"),
    }


def compare_classifier_with_experiment(endpoint_id: str, raw_score: float | None, exp_val: float | None, exp_unit: str, exp_measurement_type: str = "") -> dict[str, Any]:
    """Provide scientifically honest comparison between classifier prediction and experimental observation."""
    interp = interpret_classifier_prediction(endpoint_id, raw_score)
    if not interp.get("is_classifier"):
        return {"numeric_difference": None, "agreement_status": "NOT_CLASSIFIER", "details": ""}

    if raw_score is None or exp_val is None:
        return {
            "numeric_difference": None,
            "difference_display": "—",
            "agreement_status": "INCOMPLETE_PAIR",
            "details": "Missing prediction or experiment",
        }

    # If experiment is a quantitative potency (IC50, Ki, etc.) and prediction is a classifier probability:
    # Do NOT calculate a numeric difference.
    mtype = str(exp_measurement_type or "").upper()
    u = str(exp_unit or "").lower()
    is_potency = mtype in {"IC50", "KI", "KD"} or any(k in u for k in ["um", "nm", "µm", "mol/l"])

    if is_potency:
        # If experimental IC50 is <= 10 µM, it is qualitatively an in-vitro inhibitor
        # Convert only if consistent with standard definition, but flag agreement clearly
        exp_in_um = exp_val
        if "nm" in u:
            exp_in_um = exp_val / 1000.0
        elif "m" in u and "um" not in u and "nm" not in u and "µm" not in u:
            exp_in_um = exp_val * 1e6

        exp_class = "Inhibitor" if exp_in_um <= 10.0 else "Non-inhibitor"
        pred_class = interp.get("predicted_class", "")

        match = (exp_class == pred_class)
        agreement = "Qualitative Match" if match else "Qualitative Mismatch"

        return {
            "numeric_difference": None,
            "difference_display": "—",
            "agreement_status": agreement,
            "comparison_type": "QUALITATIVE_CLASSIFICATION_COMPARISON",
            "experimental_class": exp_class,
            "predicted_class": pred_class,
            "details": f"Experimental: {mtype or 'Potency'} {exp_val} {exp_unit} ({exp_class}) vs Prediction: {pred_class} (Score: {raw_score:.3f}). Quantitative agreement not calibrated.",
            "quantitative_gap": "QUANTITATIVE_MODEL_GAP",
        }

    return {
        "numeric_difference": None,
        "difference_display": "—",
        "agreement_status": "RELATED_MEASUREMENT",
        "details": "Related measurement — quantitative agreement not calibrated",
        "quantitative_gap": "QUANTITATIVE_MODEL_GAP",
    }
