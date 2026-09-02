"""v4.3 scientific result contract and presentation safeguards."""

from pathlib import Path

from backend.endpoint_comparison import _blank, _scientific_rows
from backend.pk_context import resolve_pk_study_context
from backend.platform_info import APP_VERSION, version_history


def test_pk_table_title_context_is_inherited_with_auditable_sources():
    context = resolve_pk_study_context(
        raw_endpoint="Cmax", raw_value="412", raw_unit="ng/mL", species="",
        context={"table_title": "Pharmacokinetic parameters following a single oral 300 mg dose in healthy volunteers"},
        source_database="FDA / Regulatory", source_record_id="fixture",
    )
    assert context["species"] == "HUMAN"
    assert context["species_source"] == "TABLE_TITLE"
    assert context["route"] == "ORAL"
    assert context["route_source"] == "TABLE_TITLE"
    assert context["dose"] == 300.0
    assert context["dose_source"] == "TABLE_TITLE"
    assert context["regimen"] == "SINGLE_DOSE"


def test_regulatory_clinical_context_can_recover_oral_route_without_claiming_a_prediction():
    context = resolve_pk_study_context(
        raw_endpoint="Cmax", raw_value="412", raw_unit="ng/mL", species="Human",
        context={"conditions": "mean Cmax in patients given 200 mg sunvozertinib"},
        source_database="FDA / Regulatory", source_record_id="NDA219839:review:L1",
    )
    assert context["route"] == "ORAL"
    assert context["route_source"] == "REGULATORY_DOCUMENT_CONTEXT"
    assert context["dose"] == 200.0


def test_pk_percentage_change_is_not_treated_as_an_absolute_parent_pk_measurement():
    context = resolve_pk_study_context(
        raw_endpoint="Cmax", raw_value="36", raw_unit="%", species="Human",
        context={"conditions": "Cmax increased by 36% following coadministration with digoxin"},
    )
    assert context["analyte"] == "OTHER_ANALYTE"
    assert context["measurement_semantics_issue"] == "MEASUREMENT_SEMANTICS_DIFFER"


def test_scientific_rows_convert_model_native_admet_display_units():
    caco = _blank("CACO2_PAPP_AB")
    caco["prediction"] = {"available": True, "display_value": -4.9045887, "unit": "log10(cm/s)", "maturity": {}}
    sol = _blank("SOLUBILITY_GENERIC")
    sol["prediction"] = {"available": True, "display_value": -2.0, "unit": "log10(mol/L)", "maturity": {}}
    rows = _scientific_rows([caco, sol])
    by_id = {row["canonical_endpoint"]: row for row in rows}
    assert by_id["CACO2_PAPP_AB"]["prediction"]["display"]["unit"] == "×10^-6 cm/s"
    assert round(by_id["CACO2_PAPP_AB"]["prediction"]["display"]["value"], 1) == 12.5
    assert by_id["SOLUBILITY_GENERIC"]["prediction"]["display"]["unit"] == "µM"
    assert by_id["SOLUBILITY_GENERIC"]["prediction"]["display"]["value"] == 10000.0


def test_equivalent_prediction_only_dog_vd_rows_have_one_primary_scientific_row():
    rows = []
    for route in ("IV", "IP", "SC"):
        row = _blank(f"DOG_PK_VD_{route}")
        row.update({"section": "PK", "species": "DOG", "route": route, "display_name": "Dog VD"})
        row["prediction"] = {"available": True, "display_value": 1.066, "unit": "L/kg", "maturity": {}}
        rows.append(row)
    scientific = _scientific_rows(rows)
    assert [row["canonical_endpoint"] for row in scientific] == ["DOG_PK_VD_SYSTEMIC"]
    assert scientific[0]["route_contexts"] == ["IV", "IP", "SC"]


def test_pk_labels_are_scientist_facing_and_keep_route_as_context():
    row = _blank("DOG_PK_VDF_ORAL_ORAL")
    row.update({"section": "PK", "species": "DOG", "route": "ORAL", "display_name": "Dog VDF ORAL"})
    row["prediction"] = {"available": True, "display_value": 1.066, "unit": "L/kg", "maturity": {}}
    result = _scientific_rows([row])[0]
    assert result["display_name"] == "Oral Vd/F"
    assert result["route"] == "ORAL"


def test_primary_ppb_display_deduplicates_repeated_public_value_but_keeps_raw_observations():
    row = _blank("HUMAN_PPB")
    row["experimental_external_candidates"] = [
        {"id": index, "raw_endpoint": "Plasma protein binding", "normalized_value": value,
         "normalized_unit": "% bound", "raw_unit": "%", "comparability": "DIRECT",
         "display_evidence_group_id": f"source-{index}"}
        for index, value in enumerate((89.0, 91.46, 91.62, 89.0, 91.46), start=1)
    ]
    result = _scientific_rows([row])[0]
    assert result["primary_experimental_display"]["observation_count"] == 5
    assert result["primary_experimental_display"]["distinct_display_count"] == 3
    assert result["primary_experimental_display"]["value"] == 89.0
    assert result["primary_experimental_display"]["additional_observation_count"] == 4
    assert result["representative_observation_id"] == 1


def test_cyp_heterogeneous_measurements_are_not_aggregated_as_a_numeric_range():
    row = _blank("CYP3A4_INHIBITION")
    row["experimental_external_candidates"] = [
        {"id": 1, "raw_endpoint": "CYP3A4 IC50", "normalized_value": 18.0, "normalized_unit": "µM", "raw_unit": "µM", "comparability": "RELATED"},
        {"id": 2, "raw_endpoint": "CYP3A4 inhibition", "normalized_value": 85.0, "normalized_unit": "%", "raw_unit": "% inhibition", "comparability": "RELATED"},
    ]
    row["prediction"] = {"available": True, "display_value": 0.797, "unit": "probability", "prediction_type": "MODEL", "maturity": {}}
    result = _scientific_rows([row])[0]
    display = result["primary_experimental_display"]
    assert display["heterogeneous"] is True
    assert all("range" not in item for item in display["measurement_types"])


def test_v43_help_and_frontend_use_one_scientific_row_contract():
    root = Path(__file__).resolve().parents[1]
    js = (root / "frontend" / "static" / "app.js").read_text(encoding="utf-8")
    assert APP_VERSION == "1.0.0"
    assert version_history()[-1]["version"] in {"v4.7", "v4.8"}
    assert "scientific_rows" in js
    assert "Other base-prediction endpoints" in js
    assert "Learning curve / leakage-safe validation" not in js


def test_v43_validation_artifact_paths_exist():
    root = Path(__file__).resolve().parents[1]
    for filename in (
        "scientific_results_ux_contract_v4_3.json",
        "sunvozertinib_final_scientific_results_v4_3.json",
        "sunvozertinib_human_clinical_pk_v4_3.json",
        "pk_context_qualification_v4_3.json",
        "qualification_completion_v4_3.json",
    ):
        assert (root / "validation" / filename).is_file()
