from backend.evidence_display_dedup import deduplicate_for_display, display_group_key
import json
from pathlib import Path


def row(source, record, value=91.2, *, doi="10.1000/example", concentration="1 uM"):
    return {
        "source": source, "source_record_id": record, "endpoint": "PPB", "canonical_endpoint_id": "ppb_human_percent_bound",
        "measurement_type": "% bound", "species": "Human", "value": value, "unit": "% bound", "normalized_value": value,
        "normalized_unit": "% bound", "raw_relation": "=", "doi": doi, "concentration": concentration,
        "reference": "Primary paper", "routing": {"canonical_endpoint_id": "ppb_human_percent_bound"},
    }


def test_cross_source_same_measurement_collapses_and_aggregates_provenance():
    display, collapsed = deduplicate_for_display([row("ChEMBL", "A"), row("PubChem", "B"), row("Paper", "C")])
    assert len(display) == 1 and collapsed == 2
    assert display[0]["display_source_count"] == 3
    assert set(display[0]["display_sources"]) == {"ChEMBL", "Paper", "PubChem"}


def test_distinct_concentration_or_doi_remains_distinct():
    display, collapsed = deduplicate_for_display([row("Paper", "A"), row("Paper", "B", concentration="10 uM"), row("Paper", "C", doi="10.1000/other")])
    assert len(display) == 3 and collapsed == 0


def test_display_group_is_stable_and_does_not_use_source_record_id():
    assert display_group_key(row("ChEMBL", "A")) == display_group_key(row("PubChem", "B"))


def test_five_drug_display_dedup_artifact_is_present():
    artifact = json.loads((Path(__file__).parents[1] / "validation" / "experimental_evidence_display_dedup_v3_4.json").read_text())
    assert {item["name"] for item in artifact["drugs"]} == {"sunvozertinib", "osimertinib", "midazolam", "warfarin", "metformin"}
    assert all(item["raw_source_records"] >= item["unique_scientific_observations"] for item in artifact["drugs"])


def test_detail_load_uses_selected_compound_for_workspace_and_pair_data():
    js = (Path(__file__).parents[1] / "frontend/static/app.js").read_text()
    assert "loadWorkspace(compound.version.id,compound.row_id)" in js
    assert "prediction-experimental-comparisons" in js
    assert "Activate Project Adapter" in js and "confirm_activation:true" in js
