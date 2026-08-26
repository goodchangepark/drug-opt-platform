import math
from pathlib import Path

import pytest
from sqlalchemy import create_engine, inspect, select
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

import backend.main as main_module
from backend.admet import (
    ADMETEndpoint, ADMETMeasurement, ADMETModelRegistry, ADMETPrediction, ADMETPredictionRun,
)
from backend.database import Base
from backend.ivive import (
    CANONICAL_CLEARANCE_UNIT, IVIVEInputSet, IVIVEMethodRegistry, IVIVERun, IVIVEUnitError,
    PHYSIOLOGY_VERSION, PhysiologicalParameterOverride, PhysiologicalParameterSet,
    calculate_ivive, calculate_validation_metrics, confidence_ceiling,
    convert_clearance_from_ml_min_kg, convert_clearance_to_ml_min_kg,
    ensure_ivive_schema, extraction_class, fraction_unbound_from_candidate,
    gather_ivive_candidates, normalize_species, observed_iv_clearance,
    resolve_physiology, scale_intrinsic_clearance, well_stirred_clearance,
)
from backend.models import Compound, CompoundVersion, Project
from backend.pk import PKNCAResult, PKStudy


@pytest.fixture()
def db_engine():
    engine = create_engine("sqlite://", connect_args={"check_same_thread": False}, poolclass=StaticPool)
    Base.metadata.create_all(engine)
    ensure_ivive_schema(engine)
    session = sessionmaker(bind=engine, autoflush=False, autocommit=False)()
    try:
        yield session, engine
    finally:
        session.close()


def make_version(db, project_name="IVIVE Project", label="IV-001", version_number=1):
    project = Project(name=project_name, target="IVIVE validation")
    db.add(project); db.flush()
    compound = Compound(project_id=project.id, compound_id=label, name=label, current_version=version_number)
    db.add(compound); db.flush()
    version = CompoundVersion(
        compound_row_id=compound.id, version_number=version_number, original_smiles="CCO",
        canonical_smiles="CCO", isomeric_smiles="CCO", inchikey=f"{label}-{version_number}",
    )
    db.add(version); db.commit()
    return project, compound, version


def physiology(db, project_id, species="Rat"):
    return resolve_physiology(db, project_id, species)


def add_manual(db, project, version, species, endpoint, value, unit, input_type="", source="EXPERIMENTAL", confidence="HIGH"):
    row = IVIVEInputSet(
        project_id=project.id, version_id=version.id, species=species, source_type=source,
        input_type=input_type, input_endpoint=endpoint, input_value=value, unit=unit,
        record_type="Experimental" if source == "EXPERIMENTAL" else "Calculated",
        model_source="test fixture", confidence=confidence, applicability_domain="IN_DOMAIN",
        provenance_json={"fixture": True},
    )
    db.add(row); db.commit(); db.refresh(row)
    return row


def add_admet_measurement(db, project, version, name, species, matrix, value, unit):
    endpoint = db.scalar(select(ADMETEndpoint).where(
        ADMETEndpoint.project_id == project.id, ADMETEndpoint.name == name
    ))
    if not endpoint:
        endpoint = ADMETEndpoint(project_id=project.id, name=name, preferred_unit=unit)
        db.add(endpoint); db.flush()
    row = ADMETMeasurement(
        version_id=version.id, endpoint_id=endpoint.id, species=species, matrix=matrix,
        value=value, unit=unit, method="fixture", source="Experimental fixture",
        provenance_json={"confidence": "HIGH"},
    )
    db.add(row); db.commit(); db.refresh(row)
    return row, endpoint


def add_prediction(db, version, endpoint, species, value, unit, confidence="LOW", output_type="regression"):
    model = ADMETModelRegistry(
        endpoint_name=endpoint.name, model_name=f"{species} quantitative clearance model", model_version="1",
        implementation_status="READY", supported_species=[species], supported_matrix=["liver microsomes"],
        output_unit=unit, source="https://example.test/model", species=species,
        output_type=output_type, is_active=True,
    )
    db.add(model); db.flush()
    run = ADMETPredictionRun(version_id=version.id, inputs_hash=f"run-{version.id}-{model.id}", status="COMPLETE")
    db.add(run); db.flush()
    prediction = ADMETPrediction(
        run_id=run.id, endpoint_id=endpoint.id, version_id=version.id, model_id=model.id,
        predicted_value=value, unit=unit, confidence=confidence, applicability_domain="IN_DOMAIN",
    )
    db.add(prediction); db.commit(); db.refresh(prediction)
    return prediction


def test_schema_architecture_and_versioned_seed(db_engine):
    db, engine = db_engine
    tables = set(inspect(engine).get_table_names())
    assert {"ivive_input_sets", "physiological_parameter_sets", "physiological_parameter_overrides",
            "pk_ivive_method_registry", "ivive_runs"}.issubset(tables)
    assert db.query(PhysiologicalParameterSet).filter_by(version=PHYSIOLOGY_VERSION).count() == 25
    method = db.query(IVIVEMethodRegistry).filter_by(method_key="WELL_STIRRED").one()
    assert "CLh" in method.equation_json["hepatic_clearance"]
    assert "not an ML model" in method.reference_json["model_type"]


@pytest.mark.parametrize("alias,expected", [
    ("mouse", "Mouse"), ("RLM", "Rat"), ("canine", "Dog"), ("cyno", "Monkey"), ("HLM", "Human"),
])
def test_species_normalization(alias, expected):
    assert normalize_species(alias) == expected


def test_raw_microsomal_scaling_exactly_once(db_engine):
    db, _ = db_engine
    project, _, _ = make_version(db)
    result = scale_intrinsic_clearance(10, "µL/min/mg protein", "RAW_MICROSOMAL", physiology(db, project.id))
    assert result["scaled_clint"] == pytest.approx(10 * 47 * 36.6 / 1000)
    assert result["scaling_count"] == 1
    assert "MPPGL" in result["equation"]


def test_raw_microsomal_ml_input_conversion(db_engine):
    db, _ = db_engine
    project, _, _ = make_version(db)
    a = scale_intrinsic_clearance(0.01, "mL/min/mg", "RAW_MICROSOMAL", physiology(db, project.id))
    b = scale_intrinsic_clearance(10, "µL/min/mg", "RAW_MICROSOMAL", physiology(db, project.id))
    assert a["scaled_clint"] == pytest.approx(b["scaled_clint"])


def test_hepatocyte_scaling_uses_hepatocellularity_not_mppgl(db_engine):
    db, _ = db_engine
    project, _, _ = make_version(db)
    result = scale_intrinsic_clearance(10, "µL/min/10^6 cells", "RAW_HEPATOCYTE", physiology(db, project.id))
    assert result["scaled_clint"] == pytest.approx(10 * 128 * 36.6 / 1000)
    assert result["mppgl_used"] is False
    assert "hepatocellularity" in result["equation"]


def test_prescaled_and_log10_clint_never_receive_physiology_scaling(db_engine):
    db, _ = db_engine
    project, _, _ = make_version(db)
    phys = physiology(db, project.id)
    linear = scale_intrinsic_clearance(60, "mL/min/kg", "PRESCALED_CLINT", phys)
    logged = scale_intrinsic_clearance(2, "log10(mL/min/kg)", "PRESCALED_CLINT", phys)
    assert linear["scaled_clint"] == 60 and linear["scaling_count"] == 0
    assert logged["scaled_clint"] == pytest.approx(100) and logged["scaling_count"] == 0
    assert logged["double_scaling_prevented"] is True and logged["mppgl_used"] is False


def test_double_scaling_type_guards_are_hard_failures(db_engine):
    db, _ = db_engine
    project, _, _ = make_version(db)
    phys = physiology(db, project.id)
    with pytest.raises(IVIVEUnitError, match="pre-scaled"):
        scale_intrinsic_clearance(20, "mL/min/kg", "RAW_MICROSOMAL", phys)
    with pytest.raises(IVIVEUnitError, match="raw assay"):
        scale_intrinsic_clearance(20, "µL/min/mg", "PRESCALED_CLINT", phys)


def test_clearance_unit_engine_round_trip():
    assert convert_clearance_to_ml_min_kg(3.6, "L/h/kg") == pytest.approx(60)
    assert convert_clearance_from_ml_min_kg(60, "L/h/kg") == pytest.approx(3.6)
    assert convert_clearance_to_ml_min_kg(60, "mL/min/kg") == 60
    with pytest.raises(IVIVEUnitError):
        convert_clearance_to_ml_min_kg(1, "mg/kg")


def test_well_stirred_synthetic_low_and_high_extraction():
    low = well_stirred_clearance(100, 0.1, 10)
    assert low["clh"] == pytest.approx(100 / 101)
    assert low["extraction_ratio"] == pytest.approx(1 / 101)
    assert low["hepatic_availability"] == pytest.approx(100 / 101)
    assert low["extraction_class"] == "Low"
    high = well_stirred_clearance(100, 1, 1000)
    assert high["clh"] == pytest.approx(1000 / 11)
    assert high["extraction_ratio"] == pytest.approx(10 / 11)
    assert high["hepatic_availability"] == pytest.approx(1 / 11)
    assert high["extraction_class"] == "High"
    assert extraction_class(0.3) == "Intermediate" and extraction_class(0.7) == "Intermediate"


@pytest.mark.parametrize("unit,value,expected", [
    ("% bound", 90, 0.1), ("fraction bound", 0.9, 0.1), ("fu", 0.1, 0.1), ("% unbound", 10, 0.1),
])
def test_fu_p_unit_handling(unit, value, expected):
    candidate = {"endpoint": "PLASMA_PROTEIN_BINDING", "unit": unit, "value": value}
    assert fraction_unbound_from_candidate(candidate) == pytest.approx(expected)


def test_fu_b_with_blood_plasma_ratio_and_plasma_approximation(db_engine):
    db, _ = db_engine
    project, _, version = make_version(db)
    add_manual(db, project, version, "Rat", "CLINT", 100, "mL/min/kg", "PRESCALED_CLINT")
    add_manual(db, project, version, "Rat", "FU_PLASMA", 0.2, "fu")
    add_manual(db, project, version, "Rat", "BLOOD_PLASMA_RATIO", 2, "ratio")
    run = calculate_ivive(db, version, "Rat")
    assert run.outputs_json["fu_p"] == pytest.approx(0.2)
    assert run.outputs_json["fu_b"] == pytest.approx(0.1)
    assert run.outputs_json["binding_basis"] == "BLOOD"
    assert run.confidence == "HIGH"

    project2, _, version2 = make_version(db, "Approximation Project", "IV-002")
    add_manual(db, project2, version2, "Rat", "CLINT", 100, "mL/min/kg", "PRESCALED_CLINT")
    add_manual(db, project2, version2, "Rat", "FU_PLASMA", 0.2, "fu")
    approximate = calculate_ivive(db, version2, "Rat")
    assert approximate.outputs_json["fu_b"] == pytest.approx(0.2)
    assert approximate.outputs_json["binding_basis"] == "PLASMA_APPROXIMATION"
    assert approximate.confidence == "MEDIUM"
    assert any("B/P is unavailable" in text for text in approximate.warnings_json)


def test_species_parameters_are_isolated_and_provenanced(db_engine):
    db, _ = db_engine
    project, _, _ = make_version(db)
    mouse = physiology(db, project.id, "Mouse")
    rat = physiology(db, project.id, "Rat")
    dog = physiology(db, project.id, "Dog")
    monkey = physiology(db, project.id, "Monkey")
    human = physiology(db, project.id, "Human")
    assert [row["HEPATIC_BLOOD_FLOW"]["value"] for row in (mouse, rat, dog, monkey, human)] == pytest.approx([120, 67.6, 30.9, 43.6, 20.7142857143])
    assert rat["MPPGL"]["value"] == 47 and human["MPPGL"]["value"] == 32
    for rows in (mouse, rat, dog, monkey, human):
        assert all(item["version"] == PHYSIOLOGY_VERSION for item in rows.values())
        assert all(item["reference"].get("doi") for item in rows.values())


def test_project_scoped_user_override_does_not_leak(db_engine):
    db, _ = db_engine
    first, _, _ = make_version(db, "Override A", "A")
    second, _, _ = make_version(db, "Override B", "B")
    db.add(PhysiologicalParameterOverride(
        project_id=first.id, species="Rat", parameter="HEPATIC_BLOOD_FLOW", value=75,
        unit="mL/min/kg", source="Study-specific Doppler flow", confidence="HIGH",
        provenance_json={"raw_value": 4.5, "raw_unit": "L/h/kg"},
    )); db.commit()
    first_rows = physiology(db, first.id, "Rat")
    second_rows = physiology(db, second.id, "Rat")
    assert first_rows["HEPATIC_BLOOD_FLOW"]["value"] == 75
    assert first_rows["HEPATIC_BLOOD_FLOW"]["source_label"] == "USER OVERRIDE"
    assert "Study-specific" in first_rows["HEPATIC_BLOOD_FLOW"]["reference"]["source"]
    assert second_rows["HEPATIC_BLOOD_FLOW"]["value"] == 67.6
    assert second_rows["HEPATIC_BLOOD_FLOW"]["source_label"] == "DEFAULT PHYSIOLOGY"


def test_experimental_hepatocyte_then_microsome_precedence_over_prediction(db_engine):
    db, _ = db_engine
    project, _, version = make_version(db)
    _, endpoint = add_admet_measurement(db, project, version, "RLM intrinsic clearance", "Rat", "Rat liver microsomes", 20, "µL/min/mg protein")
    hepatocyte, _ = add_admet_measurement(db, project, version, "Rat hepatocyte intrinsic clearance", "Rat", "Rat hepatocytes", 5, "µL/min/10^6 cells")
    prediction = add_prediction(db, version, endpoint, "Rat", 2, "log10(mL/min/kg)", confidence="LOW")
    candidates = gather_ivive_candidates(db, project.id, version.id, "Rat")["clint"]
    assert candidates[0]["origin"] == "ADMET_MEASUREMENT"
    assert candidates[0]["origin_id"] == hepatocyte.id
    assert candidates[0]["input_type"] == "RAW_HEPATOCYTE" and candidates[0]["selected"] is True
    assert any(row["origin_id"] == prediction.id and row["source_label"] == "PRED" for row in candidates)
    run = calculate_ivive(db, version, "Rat")
    assert run.inputs_snapshot_json["clint"]["origin_id"] == hepatocyte.id


def test_project_calibrated_precedes_external_prediction(db_engine):
    db, _ = db_engine
    project, _, version = make_version(db)
    endpoint = ADMETEndpoint(project_id=project.id, name="RLM intrinsic clearance")
    db.add(endpoint); db.commit()
    add_prediction(db, version, endpoint, "Rat", 2, "log10(mL/min/kg)", confidence="LOW")
    calibrated = add_manual(db, project, version, "Rat", "CLINT", 40, "mL/min/kg", "PRESCALED_CLINT", "PROJECT_CALIBRATED", "MEDIUM")
    selected = gather_ivive_candidates(db, project.id, version.id, "Rat")["clint"][0]
    assert selected["origin_id"] == calibrated.id and selected["source_type"] == "PROJECT_CALIBRATED"


def test_classification_only_result_is_never_quantitative_clint(db_engine):
    db, _ = db_engine
    project, _, version = make_version(db)
    endpoint = ADMETEndpoint(project_id=project.id, name="Microsomal stability class")
    db.add(endpoint); db.flush()
    add_prediction(db, version, endpoint, "Rat", 0.9, "probability", output_type="binary_classification")
    candidates = gather_ivive_candidates(db, project.id, version.id, "Rat")
    assert candidates["clint"] == []


def test_compound_version_and_project_isolation(db_engine):
    db, _ = db_engine
    project = Project(name="Version Isolation"); db.add(project); db.flush()
    compound = Compound(project_id=project.id, compound_id="ISO", current_version=2); db.add(compound); db.flush()
    v1 = CompoundVersion(compound_row_id=compound.id, version_number=1, original_smiles="CC", canonical_smiles="CC", isomeric_smiles="CC", inchikey="ISO-V1")
    v2 = CompoundVersion(compound_row_id=compound.id, version_number=2, original_smiles="CCC", canonical_smiles="CCC", isomeric_smiles="CCC", inchikey="ISO-V2")
    db.add_all([v1, v2]); db.commit()
    add_manual(db, project, v1, "Rat", "CLINT", 25, "mL/min/kg", "PRESCALED_CLINT")
    assert gather_ivive_candidates(db, project.id, v2.id, "Rat")["clint"] == []
    other, _, other_v = make_version(db, "Other Project", "OTHER")
    assert gather_ivive_candidates(db, other.id, other_v.id, "Rat")["clint"] == []


def test_missing_data_policy_stores_unavailable_run_without_fake_values(db_engine):
    db, _ = db_engine
    project, _, version = make_version(db)
    add_manual(db, project, version, "Rat", "CLINT", 10, "µL/min/mg", "RAW_MICROSOMAL")
    run = calculate_ivive(db, version, "Rat")
    assert run.status == "UNAVAILABLE" and run.confidence == "NOT_AVAILABLE"
    assert run.outputs_json["predicted_total_clearance"] is None
    assert run.outputs_json["non_hepatic_clearance"] == "Not modeled"
    assert "clh" not in run.outputs_json
    assert any("PPB/fu,p is unavailable" in warning for warning in run.warnings_json)


def test_low_confidence_prediction_caps_run_confidence(db_engine):
    db, _ = db_engine
    project, _, version = make_version(db)
    endpoint = ADMETEndpoint(project_id=project.id, name="RLM intrinsic clearance")
    db.add(endpoint); db.commit()
    add_prediction(db, version, endpoint, "Rat", math.log10(60), "log10(mL/min/kg)", confidence="LOW")
    add_manual(db, project, version, "Rat", "FU_PLASMA", 0.5, "fu", confidence="HIGH")
    add_manual(db, project, version, "Rat", "BLOOD_PLASMA_RATIO", 1, "ratio", confidence="HIGH")
    run = calculate_ivive(db, version, "Rat")
    assert run.status == "COMPLETE" and run.confidence == "LOW"
    assert confidence_ceiling(["LOW", "HIGH", "HIGH"]) == "LOW"


def test_observed_iv_systemic_clearance_comparison_and_over_100_warning(db_engine):
    db, _ = db_engine
    project, compound, version = make_version(db)
    study = PKStudy(
        project_id=project.id, compound_row_id=compound.id, version_id=version.id,
        study_name="Rat IV", species="Rat", route="IV", dose=1, dose_unit="mg/kg",
    )
    db.add(study); db.flush()
    nca = PKNCAResult(
        pk_study_id=study.id, version_id=version.id, analysis_version=1, is_latest=True,
        cl=10, cl_unit="mL/min/kg", nca_engine="fixture", nca_engine_version="1",
    )
    db.add(nca); db.commit()
    observed = observed_iv_clearance(db, version.id, "Rat")
    assert observed["observed_systemic_cl"] == 10
    assert "total systemic" in observed["limitation"]
    add_manual(db, project, version, "Rat", "CLINT", 1000, "mL/min/kg", "PRESCALED_CLINT")
    add_manual(db, project, version, "Rat", "FU_PLASMA", 1, "fu")
    add_manual(db, project, version, "Rat", "BLOOD_PLASMA_RATIO", 1, "ratio")
    run = calculate_ivive(db, version, "Rat")
    comparison = run.outputs_json["experimental_comparison"]
    assert comparison["estimated_hepatic_contribution"] > 1
    assert comparison["observed_systemic_cl"] == 10
    assert any("exceeds 100%" in warning for warning in run.warnings_json)


def test_literature_hepatocyte_example_matches_published_diclofenac_order(db_engine):
    """Li et al. 2025: 261 µL/min/10^6 cells -> 664.35 mL/min/kg and CLh 2.91.

    Our consensus physiology gives 663.9 and plasma approximation gives 2.87;
    the small difference reflects the paper's fu,inc=0.963 correction, which is
    outside this stage's supported input set and is not silently assumed.
    """
    db, _ = db_engine
    project, _, _ = make_version(db, "Literature", "DICLOFENAC")
    result = scale_intrinsic_clearance(261, "µL/min/10^6 cells", "RAW_HEPATOCYTE", physiology(db, project.id, "Human"))
    assert result["scaled_clint"] == pytest.approx(663.98, rel=2e-3)
    hepatic = well_stirred_clearance(20.7142857143, 0.005, result["scaled_clint"])
    assert hepatic["clh"] == pytest.approx(2.87, rel=0.02)
    assert hepatic["extraction_class"] == "Low"


def test_literature_microsome_examples_cover_low_and_high_extraction(db_engine):
    """Obach's public 29-drug table includes prescaled HLM Clint 2.1 and 189."""
    db, _ = db_engine
    project, _, _ = make_version(db, "Obach", "OBACH")
    phys = physiology(db, project.id, "Human")
    low_clint = scale_intrinsic_clearance(2.1, "mL/min/kg", "PRESCALED_CLINT", phys)
    high_clint = scale_intrinsic_clearance(189, "mL/min/kg", "PRESCALED_CLINT", phys)
    assert well_stirred_clearance(20.7142857143, 1, low_clint["scaled_clint"])["extraction_class"] == "Low"
    assert well_stirred_clearance(20.7142857143, 1, high_clint["scaled_clint"])["extraction_class"] == "High"
    assert low_clint["scaling_count"] == high_clint["scaling_count"] == 0


def test_validation_metrics_fold_error_and_limitations():
    metrics = calculate_validation_metrics([
        {"predicted": 10, "observed": 10}, {"predicted": 10, "observed": 20},
        {"predicted": 30, "observed": 10},
    ])
    assert metrics["fold_errors"] == pytest.approx([1, 2, 3])
    assert metrics["average_absolute_fold_error"] == pytest.approx(6 ** (1 / 3))
    assert metrics["within_2_fold_pct"] == pytest.approx(200 / 3)
    assert metrics["within_3_fold_pct"] == 100
    assert "total systemic" in metrics["limitation"]


def test_no_ivive_run_mutates_existing_prediction(db_engine):
    db, _ = db_engine
    project, _, version = make_version(db)
    endpoint = ADMETEndpoint(project_id=project.id, name="RLM intrinsic clearance")
    db.add(endpoint); db.commit()
    prediction = add_prediction(db, version, endpoint, "Rat", 1.5, "log10(mL/min/kg)", confidence="LOW")
    original = (prediction.predicted_value, prediction.unit, prediction.confidence, prediction.outputs_json)
    add_manual(db, project, version, "Rat", "FU_PLASMA", 0.5, "fu")
    calculate_ivive(db, version, "Rat")
    db.refresh(prediction)
    assert (prediction.predicted_value, prediction.unit, prediction.confidence, prediction.outputs_json) == original
    assert db.query(IVIVERun).count() == 1


def test_project_delete_removes_pk_and_ivive_tree(db_engine):
    db, _ = db_engine
    project, compound, version = make_version(db, "Delete IVIVE", "DELETE")
    add_manual(db, project, version, "Rat", "CLINT", 30, "mL/min/kg", "PRESCALED_CLINT")
    add_manual(db, project, version, "Rat", "FU_PLASMA", 0.5, "fu")
    add_manual(db, project, version, "Rat", "BLOOD_PLASMA_RATIO", 1, "ratio")
    calculate_ivive(db, version, "Rat")
    study = PKStudy(project_id=project.id, compound_row_id=compound.id, version_id=version.id,
                    study_name="Delete IV", species="Rat", route="IV", dose=1, dose_unit="mg/kg")
    db.add(study); db.flush()
    db.add(PKNCAResult(pk_study_id=study.id, version_id=version.id, analysis_version=1,
                       is_latest=True, cl=5, cl_unit="mL/min/kg"))
    db.commit()
    deleted = main_module.delete_project(project.id, {"confirmation_name": project.name}, db)
    assert deleted["deleted_project_ids"] == [project.id]
    assert db.query(IVIVERun).count() == 0 and db.query(IVIVEInputSet).count() == 0
    assert db.query(PKStudy).count() == 0 and db.query(PKNCAResult).count() == 0


def test_ivive_ui_contract_and_no_total_clearance_claim():
    root = Path(__file__).parents[1]
    source = (root / "frontend/static/app.js").read_text()
    css = (root / "frontend/static/app.css").read_text()
    for text in (
        "IVIVE HEPATIC CLEARANCE FOUNDATION", "Inputs · Clint", "Inputs · PPB / fu,p",
        "Species Physiology", "Run IVIVE", "Hepatic Clearance Estimate (CLh)",
        "Predicted Hepatic Availability (Fh)", "Observed IV systemic CL",
        "Estimated hepatic contribution", "Assumptions & Warnings", "IVIVE Provenance & Equations",
        "DEFAULT PHYSIOLOGY", "USER OVERRIDE", "Raw microsomal", "Raw hepatocyte", "Pre-scaled Clint",
    ):
        assert text in source
    assert "Predicted Total Clearance: Not generated" in source
    assert "Renal and other non-hepatic clearance are not modeled" in source
    for selector in (".ivive-section", ".ivive-source-exp", ".ivive-source-pred",
                     ".ivive-source-calc", ".ivive-source-user-override", ".ivive-comparison"):
        assert selector in css


def test_ivive_routes_are_registered_separately_from_admet_registry():
    paths = {route.path for route in main_module.app.routes}
    assert "/api/compound-versions/{version_id}/ivive" in paths
    assert "/api/compound-versions/{version_id}/ivive/run" in paths
    assert "/api/compound-versions/{version_id}/ivive-inputs" in paths
    assert "/api/projects/{project_id}/ivive/physiology-overrides" in paths
    assert "/api/ivive/methods" in paths
