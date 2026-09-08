#!/usr/bin/env python3
"""Audit durable accepted/candidate evidence without external network access."""
from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from sqlalchemy import select

from backend.database import SessionLocal
from backend.models import Compound, CompoundVersion, ExternalExperimentalEvidence
from backend.species_registry import normalize_species_code

OUTPUT = ROOT / "validation" / "evidence_enrichment_v333.json"
PROJECT_IDS = (1, 3, 5, 300)
ACCEPTED_STATES = {"EXTERNAL_IMPORTED", "AUTO_QUALIFIED_EXTERNAL"}


def build() -> dict:
    with SessionLocal() as db:
        records = db.execute(
            select(ExternalExperimentalEvidence, CompoundVersion, Compound)
            .join(CompoundVersion, CompoundVersion.id == ExternalExperimentalEvidence.compound_version_id)
            .join(Compound, Compound.id == CompoundVersion.compound_row_id)
            .where(
                Compound.project_id.in_(PROJECT_IDS),
                ExternalExperimentalEvidence.lifecycle_status == "ACTIVE",
            )
        ).all()
    states = Counter()
    endpoints = Counter()
    species = Counter()
    sources = Counter()
    compounds: dict[int, dict] = {}
    coverage = defaultdict(lambda: {"accepted": 0, "candidate": 0, "rejected_or_review": 0})
    for evidence, version, compound in records:
        state = evidence.evidence_state or "UNSPECIFIED"
        states[state] += 1
        endpoint = evidence.canonical_endpoint_id or evidence.raw_endpoint_name
        normalized_species = normalize_species_code(evidence.species)
        endpoints[endpoint] += 1
        species[normalized_species] += 1
        sources[evidence.source_database or "UNSPECIFIED"] += 1
        bucket = "accepted" if state in ACCEPTED_STATES else ("candidate" if state in {"EXTERNAL_CANDIDATE", "RELATED_EXTERNAL"} else "rejected_or_review")
        coverage[(endpoint, normalized_species)][bucket] += 1
        row = compounds.setdefault(compound.id, {
            "compound_id": compound.id,
            "project_id": compound.project_id,
            "name": compound.name,
            "inchikey": version.inchikey,
            "records": 0,
            "accepted": 0,
            "candidate": 0,
        })
        row["records"] += 1
        row["accepted" if state in ACCEPTED_STATES else "candidate"] += 1
    return {
        "artifact": "EVIDENCE_ENRICHMENT_V333",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "scope": {"project_ids": list(PROJECT_IDS), "network_search_performed": False},
        "contract": {
            "search_results_are_persisted": True,
            "accepted_states": sorted(ACCEPTED_STATES),
            "candidate_states": ["EXTERNAL_CANDIDATE", "RELATED_EXTERNAL", "REVIEW_REQUIRED"],
            "displayable_does_not_imply_learning_eligible": True,
        },
        "summary": {
            "compounds_with_evidence": len(compounds),
            "records": len(records),
            "accepted": sum(count for state, count in states.items() if state in ACCEPTED_STATES),
            "candidate_or_review": sum(count for state, count in states.items() if state not in ACCEPTED_STATES),
        },
        "states": dict(sorted(states.items())),
        "sources": dict(sorted(sources.items())),
        "species": dict(sorted(species.items())),
        "coverage": [
            {"endpoint": endpoint, "species": sp, **counts}
            for (endpoint, sp), counts in sorted(coverage.items())
        ],
        "compounds": sorted(compounds.values(), key=lambda row: (row["project_id"], row["compound_id"])),
    }


if __name__ == "__main__":
    payload = build()
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUTPUT.relative_to(ROOT)} ({payload['summary']['records']} records)")
