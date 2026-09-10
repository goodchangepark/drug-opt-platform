"""Compatibility shapes backed exclusively by Stable Core scientific rows."""

from __future__ import annotations

from .scientific_core_service import build_scientific_endpoint_rows


def _legacy_row(row: dict) -> dict:
    experimental = row.get("experimental")
    prediction = row.get("prediction")
    comparison = row.get("comparison") or {}
    observation = None if not experimental else {
        "id": experimental.get("observation_id"),
        "origin": "CANONICAL_EXPERIMENTAL_OBSERVATION",
        "state": experimental.get("status"),
        "raw_value": experimental.get("display_value"),
        "normalized_value": experimental.get("value"),
        "raw_unit": experimental.get("unit"),
        "normalized_unit": experimental.get("unit"),
        "species": row.get("species"),
        "context": row.get("context") or {},
        "reference": {"source": experimental.get("source"), "url": experimental.get("source_url")},
        "canonical_endpoint_id": row.get("canonical_endpoint"),
    }
    legacy_prediction = {
        "available": bool(prediction),
        "display_value": prediction.get("value") if prediction else None,
        "unit": prediction.get("unit") if prediction else "",
        "classification": prediction.get("classification") if prediction else "",
        "model_id": prediction.get("model_id") if prediction else None,
        "model_version": prediction.get("model_version") if prediction else None,
        "model_artifact_hash": prediction.get("model_artifact_hash") if prediction else None,
        "engine_version": prediction.get("engine_version") if prediction else None,
        "prediction_mode": prediction.get("mode") if prediction else None,
        "maturity": row.get("maturity"),
        "availability_status": "CURRENT" if prediction else "MODEL_UNAVAILABLE",
    }
    endpoint_id = row.get("canonical_endpoint")
    return {
        "endpoint_id": endpoint_id,
        "canonical_endpoint": endpoint_id,
        "canonical_endpoint_id": endpoint_id,
        "display_name": row.get("display_name"),
        "section": row.get("category"),
        "group": f"{row.get('species', 'UNSPECIFIED')} PK" if row.get("category") == "PK" else row.get("category"),
        "species": row.get("species"),
        "context": row.get("context") or {},
        "experimental": experimental,
        "experimental_internal": [observation] if observation else [],
        "experimental_external_imported": [],
        "experimental_external_candidates": [],
        "needs_review": [],
        "related_evidence": [],
        "experimental_observations": [observation] if observation else [],
        "experimental_display_value": experimental.get("value") if experimental else None,
        "experimental_display_unit": experimental.get("unit") if experimental else "",
        "representative_observation_id": experimental.get("observation_id") if experimental else None,
        "primary_experimental_display": None if not experimental else {
            "value": experimental.get("value"), "unit": experimental.get("unit"),
            "label": f"{experimental.get('display_value')} {experimental.get('unit', '')}".strip(),
            "observation_count": 1,
        },
        "prediction": legacy_prediction,
        "difference": comparison,
        "comparison": comparison,
        "semantic_status": "DIRECTLY_COMPARABLE" if comparison.get("numeric_pairable") else "NOT_COMPARABLE",
        "quantitative_prediction": None,
        "maturity": row.get("maturity"),
        "model_status": row.get("model_status"),
        "stable_core_row": row,
    }

def build_legacy_comparison_adapter(db, version_id: int) -> dict:
    canonical = build_scientific_endpoint_rows(db, version_id)
    rows = [_legacy_row(row) for row in canonical.get("rows", [])]
    paired = sum(1 for row in rows if row["comparison"].get("numeric_pairable"))
    accepted = sum(1 for row in rows if row.get("experimental"))
    return {
        "contract": "LegacyComparisonAdapter/StableCore-v1.2",
        "authority": canonical.get("contract"),
        "compound_version_id": version_id,
        "current_engine": canonical.get("current_engine"),
        "endpoints": rows,
        "scientific_rows": rows,
        "summary": {
            "comparison_coverage": {
                "accepted_observations": accepted,
                "successfully_paired": paired,
                "direct_pairs": paired,
                "converted_pairs": 0,
                "pairable_evidence_without_pair": 0,
            },
            "qualification": {},
        },
    }
