"""Canonical, side-effect-free contracts for project learning.

Project learning is an overlay on frozen Global predictions.  These helpers
make the scientific rules explicit without changing the runtime database or
the Global Prediction Engine.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping


PROJECT_LEARNING_VERSION = "drugopt-project-learning-v7.0"
GLOBAL_ENGINE_POLICY = "drugopt-prediction-engine-v1@1.0.0"
GLOBAL_ENGINE_HASH = "12757ab197b5a70d8ea1754678d9a342ab0b6ea0d82f2896bebb767d686bbdeb"


@dataclass(frozen=True)
class ProjectPredictionResidual:
    """Canonical projection of an immutable prediction/experiment pair.

    Residual is defined consistently as experimental minus Global prediction.
    The existing SQL pair row remains the storage record; this is its stable
    scientific contract and does not rewrite historical rows.
    """

    project_id: int
    compound_id: int
    endpoint: str
    global_model_version: str
    global_prediction: float
    global_uncertainty: float | None
    global_ad: str
    experimental_value: float
    residual: float
    absolute_error: float
    fold_error: float | None
    prediction_timestamp: str | None
    experiment_timestamp: str | None
    prospective_status: str
    adapter_version: str
    project_prediction: float | None
    project_residual: float | None

    @classmethod
    def from_pair(cls, row):
        if row.base_prediction is None or row.experimental_value is None:
            raise ValueError("A residual requires numeric Global and experimental values")
        residual = float(row.experimental_value) - float(row.base_prediction)
        base, observed = abs(float(row.base_prediction)), abs(float(row.experimental_value))
        fold = max(base / observed, observed / base) if base > 0 and observed > 0 else None
        snapshot = row.snapshot_json or {}
        return cls(
            project_id=row.project_id,
            compound_id=row.compound_version_id,
            endpoint=row.endpoint_name,
            global_model_version=str(snapshot.get("engine_policy") or GLOBAL_ENGINE_POLICY),
            global_prediction=float(row.base_prediction),
            global_uncertainty=snapshot.get("global_uncertainty"),
            global_ad=str(snapshot.get("global_ad") or snapshot.get("ood_applicability") or "UNKNOWN"),
            experimental_value=float(row.experimental_value),
            residual=residual,
            absolute_error=abs(residual),
            fold_error=fold,
            prediction_timestamp=row.prediction_created_at.isoformat() if row.prediction_created_at else None,
            experiment_timestamp=row.experiment_created_at.isoformat() if row.experiment_created_at else None,
            prospective_status=("PROSPECTIVE_PROJECT_VALIDATION" if row.pair_class == "TRUE_PROSPECTIVE" else "RETROSPECTIVE_OOF"),
            adapter_version=row.adapter_version or "",
            project_prediction=row.project_prediction,
            project_residual=(float(row.experimental_value) - float(row.project_prediction)) if row.project_prediction is not None else None,
        )


def loco_training(events: Iterable, target_compound_id: int | str) -> list:
    """Return adapter events with every target-compound observation removed."""
    target = str(target_compound_id)
    return [event for event in events if str(event.compound_version_id) != target]


def recommend_prediction(global_prediction: float, project_prediction: float | None,
                         *, status: str, global_error: float | None = None,
                         project_error: float | None = None) -> dict:
    """Select Project only after a strict held-out improvement gate."""
    improved = (
        project_prediction is not None
        and status in {"PROJECT_VALIDATED_CANDIDATE", "PROJECT_RECOMMENDED"}
        and global_error is not None and project_error is not None
        and project_error < global_error
    )
    return {
        "value": float(project_prediction if improved else global_prediction),
        "source": "PROJECT" if improved else "GLOBAL",
        "reason": "LOCO_VALIDATED_PROJECT_IMPROVEMENT" if improved else "GLOBAL_RETAINED_UNTIL_LOCO_IMPROVEMENT",
        "project_candidate_present": project_prediction is not None,
    }


def freeze_prospective_prediction(*, project_id: int, compound_id: int,
                                  endpoint: str, global_prediction: float,
                                  global_uncertainty: float | None, global_ad: str,
                                  model_version: str, timestamp: str,
                                  experiment_known: bool = False) -> dict:
    """Create a prediction snapshot; reject a false prospective claim."""
    if experiment_known:
        raise ValueError("PROSPECTIVE_FREEZE_REQUIRES_EXPERIMENT_UNKNOWN")
    return {
        "project_id": project_id, "compound_id": compound_id, "endpoint": endpoint,
        "global_prediction": float(global_prediction), "global_uncertainty": global_uncertainty,
        "global_ad": global_ad, "model_version": model_version,
        "timestamp": timestamp, "prospective_status": "PROSPECTIVE_PROJECT_PREDICTION",
        "immutable": True,
    }


def adapter_capability(endpoint: str, *, global_status: str = "PRODUCTION_STABLE") -> dict:
    """Endpoint gate preventing false precision over weak Global routes."""
    weak = {"VDSS", "pKa", "PKA", "LOGD_7_4", "HUMAN_TOTAL_IV_CL"}
    if endpoint in weak or global_status in {"DATA_LIMITED", "MODEL_UNAVAILABLE"}:
        return {"endpoint": endpoint, "capability": "INSUFFICIENT_GLOBAL_MODEL", "eligible": False}
    return {"endpoint": endpoint, "capability": "ADAPTER_SUPPORTED", "eligible": True}


def project_registry_entry(*, project_id: int, endpoint: str, adapter_version: str,
                           global_model_version: str, dataset_hash: str,
                           experimental_n: int, validation: Mapping,
                           status: str, artifact_hash: str, created_at: str) -> dict:
    return {
        "project_id": project_id, "endpoint": endpoint,
        "adapter_version": adapter_version, "global_model_dependency": global_model_version,
        "training_dataset_hash": dataset_hash, "experimental_N": experimental_n,
        "validation_method": validation.get("method", "NOT_RUN"),
        "validation_metrics": dict(validation), "status": status,
        "artifact_hash": artifact_hash, "created_at": created_at,
    }
