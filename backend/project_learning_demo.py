"""Disposable, deterministic acceptance fixture for project learning UX.

This module never touches the application database and never changes Engine v1.
It exercises the same project-adaptation fitter used by the production preview
endpoint so the UI contract can be tested without activating a real project.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

from .project_adaptation_v2 import QualifiedEvidencePair, fit_project_adapter


@dataclass(frozen=True)
class LearningDemo:
    compounds: list[dict]
    tiers: list[dict]
    candidate: dict
    compound_six_before_experiment: dict
    compound_six_after_experiment: dict


def _compound(index: int, *, measured: bool = True) -> dict:
    # These are disposable public-like structures used only for fingerprinting.
    smiles = ["CCO", "CCCO", "CCCC", "CCN", "CCCl", "CCBr"][index - 1]
    model_predictions = {"A": 1.8, "B": 1.0, "C": 1.3}
    return {
        "compound_version_id": 9000 + index,
        "smiles": smiles,
        "prediction_created_at": "2026-01-01T00:00:00+00:00",
        "experimental_created_at": "2026-02-01T00:00:00+00:00" if measured else None,
        "model_predictions": model_predictions,
        "experimental_value": 1.0,
    }


def run_synthetic_learning_demo() -> dict:
    """Return a six-compound, in-memory learning progression.

    Five pre-experimental freezes form the candidate dataset.  Activation is
    represented explicitly in the returned fixture; it is never performed by
    this helper or by any real project endpoint.
    """
    compounds = [_compound(index, measured=index <= 5) for index in range(1, 7)]
    evidence = [
        QualifiedEvidencePair(
            evidence_id=f"DEMO-E-{item['compound_version_id']}",
            compound_version_id=item["compound_version_id"],
            smiles=item["smiles"],
            endpoint_id="solubility",
            value=item["experimental_value"],
            frozen_predictions=item["model_predictions"],
            origin="EXPERIMENTAL_INTERNAL",
            source_quality="A",
            comparability_status="DIRECTLY_COMPARABLE",
        )
        for item in compounds[:5]
    ]
    global_weights = {"A": 0.6, "B": 0.2, "C": 0.2}
    tiers = []
    for count in range(1, 6):
        result = fit_project_adapter("solubility", evidence[:count], global_weights)
        tiers.append({
            "n": count,
            "status": result.status if count >= 5 else "BASE_ONLY",
            "activation_decision": result.activation_decision if count >= 5 else "BASE_RETAINED",
            "effective_n": result.effective_n,
            "project_weights": result.project_weights if count >= 5 else result.global_weights,
            "base_validation_error": result.base_validation_error,
            "adapted_validation_error": result.adapted_validation_error,
        })

    candidate_result = fit_project_adapter("solubility", evidence, global_weights)
    candidate = {
        "status": candidate_result.status,
        "activation_decision": candidate_result.activation_decision,
        "requires_explicit_activation": True,
        "base_weights": candidate_result.global_weights,
        "candidate_weights": candidate_result.project_weights,
        "effective_n": candidate_result.effective_n,
        "base_validation_error": candidate_result.base_validation_error,
        "candidate_validation_error": candidate_result.adapted_validation_error,
    }

    sixth = compounds[5]
    base_value = sum(global_weights[key] * sixth["model_predictions"][key] for key in global_weights)
    project_value = sum(candidate_result.project_weights[key] * sixth["model_predictions"][key] for key in global_weights)
    before = {
        "compound_version_id": sixth["compound_version_id"],
        "experimental": None,
        "base_prediction": base_value,
        "project_prediction": project_value,
        "maturity": {"level": 2, "label": "Early Adaptation", "stars": "★★☆☆☆"},
        "prediction_source": "Project-adapted Prediction",
        "prior_qualified_compounds": 5,
    }
    after = before | {
        "experimental": sixth["experimental_value"],
        "base_error": abs(base_value - sixth["experimental_value"]),
        "project_error": abs(project_value - sixth["experimental_value"]),
    }
    return asdict(LearningDemo(compounds, tiers, candidate, before, after))

