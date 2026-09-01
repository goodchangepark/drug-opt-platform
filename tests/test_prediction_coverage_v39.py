from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session
from sqlalchemy.pool import StaticPool

from backend.admet import ADMETPredictionRun, PredictionEndpointSnapshot, ensure_admet_schema
from backend.canonical_endpoints import (
    PREDICTION_DERIVED, PREDICTION_MECHANISTIC, PREDICTION_RULE,
    canonicalize_prediction_endpoint, normalize_experimental_observation,
    prediction_source_type,
)
from backend.database import Base
from backend.endpoint_comparison import _pk_snapshot_values, persist_pk_prediction_snapshots
from backend.ivive import PKParameterSet, ensure_ivive_schema
from backend.models import Compound, CompoundVersion, Project
from backend.simulation import ensure_simulation_schema


def test_prediction_source_classification_is_explicit():
    assert prediction_source_type(source="Hepatic_IVIVE", endpoint="CL") == PREDICTION_MECHANISTIC
    assert prediction_source_type(source="SyGMa empirical rules", endpoint="soft spots") == PREDICTION_RULE
    assert prediction_source_type(source="Calculated normalized value", endpoint="F", default=PREDICTION_DERIVED) == PREDICTION_DERIVED


def test_pk_aliases_and_units_remain_semantically_canonical():
    assert canonicalize_prediction_endpoint("Apparent oral clearance", species="SD rat", context={"route": "oral"})["canonical_endpoint_id"] == "RAT_PK_CLF_ORAL_ORAL"
    assert normalize_experimental_observation("Cmax", 1, "µg/L", species="Rat", context={"route": "IV"})["normalized_unit"] == "ng/mL"
    assert normalize_experimental_observation("AUC0-inf", 1, "µg*h/L", species="Rat", context={"route": "IV"})["normalized_unit"] == "ng*h/mL"
    assert normalize_experimental_observation("F", 0.4, "fraction", species="Rat", context={"route": "PO"})["normalized_value"] == 40.0


def test_pk_parameter_set_emits_only_real_non_unavailable_outputs():
    pset = PKParameterSet(
        project_id=1, compound_row_id=1, version_id=1, species="Rat", route="IV",
        cl_value=12.0, cl_source_type="HEPATIC_IVIVE", v_value=1.2,
        v_source_type="PREDICTED_VD", v_type="Vss", f_predicted=100.0,
    )
    emitted = list(_pk_snapshot_values(pset))
    # IV is the reference arm for absolute bioavailability, not an F
    # prediction.  It must never be indexed as an oral F endpoint.
    assert {row[0] for row in emitted} == {"CL", "VSS"}
    assert all(row[3] in {PREDICTION_MECHANISTIC, PREDICTION_DERIVED} for row in emitted)
    oral = PKParameterSet(
        project_id=1, compound_row_id=1, version_id=1, species="Rat", route="PO",
        f_predicted=40.0,
        provenance_json={"absorption_info": {"fa_value": .8, "fg_value": .9, "fh_value": .7}},
    )
    assert [row[0] for row in _pk_snapshot_values(oral)] == ["F"]
    unavailable = PKParameterSet(project_id=1, compound_row_id=1, version_id=1, species="Rat", route="PO", cl_source_type="MODEL_UNAVAILABLE")
    assert list(_pk_snapshot_values(unavailable)) == []


def test_persisted_pk_outputs_get_durable_snapshot_identity():
    engine = create_engine("sqlite://", connect_args={"check_same_thread": False}, poolclass=StaticPool)
    Base.metadata.create_all(engine)
    ensure_admet_schema(engine); ensure_ivive_schema(engine); ensure_simulation_schema(engine)
    with Session(engine) as db:
        project = Project(name="coverage fixture")
        db.add(project); db.flush()
        compound = Compound(project_id=project.id, compound_id="C1", name="C1", current_version=1)
        db.add(compound); db.flush()
        version = CompoundVersion(compound_row_id=compound.id, version_number=1, original_smiles="CCO", canonical_smiles="CCO", isomeric_smiles="CCO", inchikey="fixture", svg="", highlighted_svg="")
        db.add(version); db.flush()
        db.add(PKParameterSet(project_id=project.id, compound_row_id=compound.id, version_id=version.id, species="Rat", route="IV", cl_value=12, cl_source_type="HEPATIC_IVIVE"))
        db.flush()
        result = persist_pk_prediction_snapshots(db, version.id)
        db.commit()
        snapshots = db.scalars(select(PredictionEndpointSnapshot).where(PredictionEndpointSnapshot.compound_version_id == version.id)).all()
        assert result["created"] == 1
        assert len(snapshots) == 1
        assert snapshots[0].prediction_type == PREDICTION_MECHANISTIC
        assert snapshots[0].snapshot_json["source_type"] == PREDICTION_MECHANISTIC
        assert db.get(ADMETPredictionRun, snapshots[0].prediction_run_id) is not None


def test_species_and_route_are_not_collapsed_by_prediction_aliases():
    rat = canonicalize_prediction_endpoint("CL", species="Rat", route="IV")
    human = canonicalize_prediction_endpoint("CL", species="Human", route="IV")
    oral = canonicalize_prediction_endpoint("CL/F", species="Rat", route="PO")
    assert rat["canonical_endpoint_id"] == "RAT_PK_CL_IV"
    assert human["canonical_endpoint_id"] == "HUMAN_PK_CL_IV"
    assert oral["canonical_endpoint_id"] == "RAT_PK_CLF_ORAL_ORAL"
    assert len({rat["canonical_endpoint_id"], human["canonical_endpoint_id"], oral["canonical_endpoint_id"]}) == 3


def test_classification_and_mechanistic_outputs_are_not_called_model_values():
    assert prediction_source_type(source="SyGMa", prediction_type="RULE_ESTIMATE") == PREDICTION_RULE
    assert prediction_source_type(source="Stage-5 PK simulation", prediction_type="MECHANISTIC_ESTIMATE") == PREDICTION_MECHANISTIC
