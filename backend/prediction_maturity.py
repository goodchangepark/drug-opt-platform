"""Endpoint-specific maturity metadata for project adaptation, never confidence."""
from __future__ import annotations

from dataclasses import dataclass

LEVEL_LABELS = {1: "Base Prediction", 2: "Early Adaptation", 3: "Project Adapted", 4: "Series Adapted", 5: "Mature Project Prediction"}

@dataclass(frozen=True)
class Maturity:
    level: int
    label: str
    stars: str
    active: bool

    def to_dict(self):
        return {"level": self.level, "label": self.label, "stars": self.stars, "active": self.active,
                "aria_label": f"Prediction maturity {self.level} of 5 — {self.label}",
                "help": "Stars indicate the maturity of project-specific experimental adaptation, not an absolute guarantee of predictive accuracy."}

def maturity_for_adapter(*, status: str, effective_n: float, activation_decision: str, stable_history_count: int = 0, representative_series: bool = False, compatible_evidence_only: bool = True) -> Maturity:
    """N is necessary but never sufficient: a validated active adapter is required."""
    active = activation_decision == "ACTIVATED" and compatible_evidence_only
    level = 1
    if active and effective_n >= 5 and status == "LIGHT_PROJECT_ADAPTATION": level = 2
    if active and effective_n >= 10 and status == "REGULARIZED_PROJECT_ENSEMBLE": level = 3
    if active and effective_n >= 20 and status == "LOCAL_SERIES_ADAPTATION" and representative_series: level = 4
    if level == 4 and effective_n >= 40 and representative_series and stable_history_count >= 3: level = 5
    return Maturity(level, LEVEL_LABELS[level], "★" * level + "☆" * (5-level), level >= 2)
