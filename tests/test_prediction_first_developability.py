from __future__ import annotations

from fastapi.testclient import TestClient
from sqlalchemy import create_engine, select
from sqlalchemy.orm import sessionmaker

from backend.current_prediction_publisher import (
    publish_global_current_predictions,
    publish_property_current_predictions,
)
from backend.database import Base
from backend.developability_profile import build_developability_profile
from backend.main import app
from backend.models import Compound, CompoundVersion, Project
from backend.stable_core import CurrentPredictionSnapshot, admit_current_prediction, ensure_stable_core_schema


def _scientific_fixture(tmp_path):
    engine = create_engine(f"sqlite:///{tmp_path / 'developability.sqlite'}")
    Base.metadata.create_all(engine)
    ensure_stable_core_schema(engine)
    db = sessionmaker(bind=engine)()
    project = Project(name="Prediction-first fixture")
    db.add(project); db.flush()
    compound = Compound(project_id=project.id, compound_id="PF-1", name="PF-1", current_version=1)
    db.add(compound); db.flush()
    version = CompoundVersion(
        compound_row_id=compound.id, version_number=1, original_smiles="CCO",
        canonical_smiles="CCO", isomeric_smiles="CCO",
        inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        properties_json={
            "exact_molecular_weight": 46.0419, "clogp": -0.0014, "tpsa": 20.23,
            "hbd": 1, "hba": 1, "rotatable_bonds": 0, "fraction_csp3": 1.0,
            "qed": 0.4068, "formal_charge": 0, "heavy_atom_count": 3,
        },
        calculation_json={"ionization": {"primary_pka": None, "physiological_state_7_4": {"estimated_logd74": -0.0014}}},
    )
    db.add(version); db.commit()
    return db, version


def test_profile_api_is_complete_prediction_first_and_read_only():
    client = TestClient(app)
    before = client.get("/api/compound-versions/11/developability-profile")
    assert before.status_code == 200
    payload = before.json()
    assert payload["authority"] == "ScientificEndpointRow/stable-core-v1"
    assert set(payload["groups"]) == {
        "physchem", "absorption", "distribution", "metabolic_stability",
        "cyp", "transporters", "safety", "metabolism", "pk",
    }
    rows = {row["query_endpoint"]: row for row in payload["availability_catalog"]}
    assert rows["PAMPA_PERMEABILITY"]["status"] == "MODEL_UNAVAILABLE"
    assert rows["PLASMA_STABILITY"]["status"] == "MODEL_UNAVAILABLE"
    assert rows["HUMAN_PK_AUC_ORAL"]["availability"] == "CONTEXT_REQUIRED"
    assert rows["HUMAN_PK_F_ORAL"]["experimental"] is not None
    assert rows["HUMAN_PK_F_ORAL"]["prediction"] is None
    assert rows["HUMAN_PK_F_ORAL"]["status"] == "EXPERIMENTAL_AVAILABLE"
    assert client.get("/api/compound-versions/999999999/developability-profile").status_code == 404


def test_qualified_property_and_global_values_publish_only_through_admission(tmp_path):
    db, version = _scientific_fixture(tmp_path)
    try:
        property_results = publish_property_current_predictions(db, version)
        global_results = publish_global_current_predictions(db, version, {
            "SOLUBILITY_GENERIC": {
                "production_prediction": -2.75, "applicability_domain": "IN_DOMAIN",
                "nearest_neighbor_similarity": 0.72, "prediction_uncertainty": 0.31,
            },
            "CACO2_PAPP_AB": {
                "production_prediction": -5.1, "applicability_domain": "BORDERLINE",
                "nearest_neighbor_similarity": 0.42, "prediction_uncertainty": 0.45,
            },
        })
        db.commit()
        assert any(row["endpoint"] == "MW" and row["status"] == "CALCULATED_AND_PUBLISHED" for row in property_results)
        assert all(row["status"] == "CALCULATED_AND_PUBLISHED" for row in global_results)
        snapshots = list(db.scalars(select(CurrentPredictionSnapshot)))
        assert snapshots and all(admit_current_prediction(db, row).eligible for row in snapshots)
        profile = build_developability_profile(db, version.id)
        rows = {row["query_endpoint"]: row for row in profile["availability_catalog"]}
        assert rows["MW"]["prediction"]["value"] == 46.0419
        assert rows["SOLUBILITY_GENERIC"]["prediction"]["value"] == -2.75
        assert rows["CACO2_PAPP_AB"]["AD"]["classification"] == "BORDERLINE"
        assert rows["PAMPA_PERMEABILITY"]["prediction"] is None
    finally:
        db.close()


def test_frontend_core_tabs_and_mobile_prediction_first_contract():
    source = open("frontend/static/app.js", encoding="utf-8").read()
    css = open("frontend/static/app.css", encoding="utf-8").read()
    assert "const tabs=['overview','properties','admet','metabolism','pk','evidence','history'];" in source
    assert "'Developability Summary'" in source
    assert "['Endpoint','Prediction','Experimental','Difference','Status','AD','Model Maturity']" in source
    assert "profileSection('metabolic_stability','Metabolic Stability (MS)')" in source
    assert "profileSection('transporters','Transporters')" in source
    assert "@media (max-width: 640px)" in css
    assert ".developability-table thead { display:none; }" in css
    workflow = open("backend/main.py", encoding="utf-8").read()
    predict_block = workflow[workflow.index("def run_compound_prediction_workflow"):workflow.index("@app.post(\"/api/compounds/{row_id}/predict-all\"")]
    assert "run_pk_simulation" not in predict_block
    assert "contextual-pk-excluded" in predict_block
