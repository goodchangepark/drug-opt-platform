"""Write deterministic v4.4B public-PK benchmark artifacts from curated code."""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from backend.public_pk_benchmark import (  # noqa: E402
    DEVELOPMENT, FINAL_EVALUATION, REVIEW_QUEUE, SOURCES, baseline_rows,
    benchmark_package, coverage, freeze_compound_split, mechanistic_verification,
)


def write(name: str, value: object) -> None:
    (ROOT / "validation" / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    package = benchmark_package(); rows = package["observations"]
    split = freeze_compound_split(rows)
    development = baseline_rows(rows, split=split, partition=DEVELOPMENT, mode=DEVELOPMENT)
    # FINAL_EVALUATION is the explicit, audited mode which alone may access
    # final-holdout targets; this baseline is recorded once at dataset freeze.
    holdout = baseline_rows(rows, split=split, partition="HOLDOUT", mode=FINAL_EVALUATION)
    source_counts = Counter(row["source_type"] for row in rows)
    quality_counts = Counter(row["quality_tier"] for row in rows)
    write("public_pk_benchmark_v1.json", package)
    write("public_pk_benchmark_sources_v1.json", {"benchmark_version": package["benchmark_version"], "sources": SOURCES, "observation_source_counts": dict(sorted(source_counts.items())), "pkdb_rest_api_inspected": {"url": "https://pk-db.com/api/v1/swagger/", "status": "INSPECTED_DISCOVERY_SOURCE", "use": "PK-DB documents study/intervention/output context; no PK-DB value is auto-qualified without source-level review."}})
    write("public_pk_benchmark_quality_v1.json", {"benchmark_version": package["benchmark_version"], "quality_tier_counts": dict(sorted(quality_counts.items())), "identity_qualified": len(rows), "context_qualified": len(rows), "benchmark_qualified": len(rows), "exclusion_policy": "Context-incomplete ranges and records with missing dose where dose is endpoint-relevant remain in the review queue."})
    write("pk_benchmark_review_queue_v4_4b.json", {"benchmark_version": package["benchmark_version"], "count": len(REVIEW_QUEUE), "by_reason": dict(sorted(Counter(item["reason"] for item in REVIEW_QUEUE).items())), "items": REVIEW_QUEUE})
    write("pk_benchmark_split_v4_4b.json", split)
    write("pk_benchmark_development_baseline_v4_4b.json", {"benchmark_version": package["benchmark_version"], "engine": "drugopt-pk-engine-v1", "track": "PREDICTIVE_VALIDATION", "partition": DEVELOPMENT, "holdout_targets_accessed": False, "records": development})
    write("pk_benchmark_holdout_baseline_v4_4b.json", {"benchmark_version": package["benchmark_version"], "engine": "drugopt-pk-engine-v1", "track": "PREDICTIVE_VALIDATION", "partition": "HOLDOUT", "mode": FINAL_EVALUATION, "frozen_once": True, "used_for_tuning": False, "records": holdout})
    write("pk_benchmark_coverage_v4_4b.json", {"benchmark_version": package["benchmark_version"], "coverage": coverage(rows), "mechanistic_verification": mechanistic_verification(), "performance_policy": "No predictive accuracy metric is calculated when the unchanged fail-closed engine returns no quantitative Track-A prediction."})


if __name__ == "__main__":
    main()
