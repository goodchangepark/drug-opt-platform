"""Generate deterministic Stage 4D-7 review artifacts from authoritative prior evidence."""

from __future__ import annotations

import json
from pathlib import Path
import sys


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from backend.preexperimental_optimization import (
    ROOT,
    STAGE4D7_POLICY_VERSION,
    _artifact,
    build_candidate_results,
    build_endpoint_accuracy_matrix,
    canonical_hash,
)
from backend.endpoint_strategy_registry import get_all_strategies


OUT = ROOT / "validation"


def write(name: str, payload: dict) -> None:
    (OUT / name).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main() -> None:
    matrix = build_endpoint_accuracy_matrix()
    candidates = build_candidate_results()
    source_paths = sorted({
        source["path"] for row in matrix["endpoints"] for source in row["source_artifacts"]
    })
    sources = [_artifact(path) for path in source_paths]
    bootstrap = {
        "artifact": "STAGE4D7_BOOTSTRAP_RESULTS",
        "review_version": STAGE4D7_POLICY_VERSION,
        "new_bootstrap_fitted": False,
        "reason": "Reuses authoritative 1,000-replicate Stage 4D bootstrap analyses; no raw validation data are re-fit in Stage 4D-7.",
        "source_artifacts": [{"path": item["path"], "sha256": item["sha256"]} for item in sources if "bootstrap" in item["path"]],
        "solubility": _artifact("validation/stage4d3a2_m1_bootstrap.json")["data"]["bootstrap_vs_m1"],
        "caco2": _artifact("validation/stage4d2c_bootstrap_comparison.json")["data"]["comparisons"]["Caco-2"]["bootstrap_results"],
    }
    calibration = {
        "artifact": "STAGE4D7_CALIBRATION_RESULTS",
        "review_version": STAGE4D7_POLICY_VERSION,
        "cyp3a4": {
            "source": "validation/stage4d3b1a_calibration.json",
            "ece": _artifact("validation/stage4d3b1a_calibration.json")["data"]["expected_calibration_error"],
            "decision": "RESEARCH_ONLY; Stage 4D-5 endpoint requirements and non-inferiority margin remain unconfigured.",
        },
        "herg": {
            "source": "validation/stage4d3b2a_calibration.json",
            "split": "scaffold_aware: calibration N=546, untouched test N=182",
            "raw": _artifact("validation/stage4d3b2a_calibration.json")["data"]["holdout_calibration_comparison"]["m1_raw"],
            "platt": _artifact("validation/stage4d3b2a_calibration.json")["data"]["holdout_calibration_comparison"]["m1_platt"],
            "decision": "CALIBRATION_RESEARCH_ONLY; raw M1 threshold 0.50 remains production.",
        },
    }
    domain = {
        "artifact": "STAGE4D7_DOMAIN_ANALYSIS",
        "review_version": STAGE4D7_POLICY_VERSION,
        "domain_aware_selector_fitted": False,
        "reason": "No nested, leakage-safe evidence establishes M2 superiority in an AD partition; OOD never triggers an automatic model switch.",
        "solubility": "Stage 4D-3A2 retained M1 globally; no generalizable selector was validated.",
        "caco2": "N≈34 is inadequate for a stable domain selector.",
        "cyp3a4": "Dynamic project/series adaptation previously had NO_ADAPTIVE_VALUE.",
        "herg": {
            "source": "validation/stage4d3b2a_model_metrics.json",
            "out_of_domain": "N=2; insufficient for correction or selector claims.",
            "m2": "Limited discriminatory complementarity; retained calibration-supporting shadow only.",
        },
    }
    decisions = {
        "artifact": "STAGE4D7_PRODUCTION_POLICY_DECISIONS",
        "review_version": STAGE4D7_POLICY_VERSION,
        "activation_authorization_consumed": False,
        "manual_promotion_required": True,
        "reason_no_activation": "Stage 4D-5 endpoint-specific minimum sample and non-inferiority requirements remain unconfigured and fail closed.",
        "previous_policy_version": "stage4d4-endpoint-strategy-v1",
        "new_policy_version": "stage4d4-endpoint-strategy-v1",
        "production_policy_changes": [],
        "rollback": "No activation occurred; existing deterministic Stage 4D-4 rollback metadata remains authoritative.",
        "decisions": [{
            "endpoint_name": row["endpoint_name"],
            "endpoint_id": row["endpoint_id"],
            "decision": row["production_decision"],
            "flags": row["decision_flags"],
        } for row in matrix["endpoints"]],
    }
    baseline = {
        "artifact": "STAGE4D7_PREEXPERIMENTAL_BASELINE",
        "review_version": STAGE4D7_POLICY_VERSION,
        "source_stage4d5_baseline": _artifact("validation/stage4d5_production_baseline.json")["sha256"],
        "production_policy_changes": [],
        "historical_freezes_mutated": False,
        "policies": [{
            "endpoint_name": policy.endpoint_name,
            "endpoint_id": policy.endpoint_id,
            "policy_version": policy.policy_version,
            "primary_strategy": policy.primary_strategy.value,
            "primary_models": policy.primary_model_ids,
            "calibration_status": policy.calibration_status.value,
            "decision_threshold": policy.decision_threshold,
        } for policy in get_all_strategies().values()],
    }
    for name, payload in {
        "stage4d7_endpoint_accuracy_matrix.json": matrix,
        "stage4d7_candidate_strategy_results.json": candidates,
        "stage4d7_bootstrap_results.json": bootstrap,
        "stage4d7_domain_analysis.json": domain,
        "stage4d7_calibration_results.json": calibration,
        "stage4d7_production_policy_decisions.json": decisions,
        "stage4d7_preexperimental_baseline.json": baseline,
    }.items():
        write(name, payload)


if __name__ == "__main__":
    main()
