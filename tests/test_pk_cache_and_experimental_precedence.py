"""
Tests for PK Cache Instant Loading, Auto-Recomputation upon Experimental Data Ingestion,
and Default Precedence with Provenance Comments.
"""

import pytest
from fastapi.testclient import TestClient
from backend.main import app, run_compound_prediction_workflow
from backend.database import SessionLocal
from backend.models import Compound, CompoundVersion, Project
from backend.pk import PKStudy, PKObservation, PKNCAResult
from backend.ivive import PKParameterSet, get_multi_species_pk_profile, get_pk_foundation_profile, refresh_pk_and_ivive_for_version

client = TestClient(app)


def test_pk_cache_instant_loading_and_multi_species():
    # 1. Create a project and compound
    p_res = client.post("/api/projects", json={"name": "PK Cache Test", "target_name": "Kinase"})
    assert p_res.status_code == 201
    project_id = p_res.json()["id"]

    try:
        c_res = client.post(f"/api/projects/{project_id}/compounds", json={
            "compound_id": "CMPD-CACHE-01",
            "smiles": "CC(=O)Nc1ccc(O)cc1"  # Acetaminophen
        })
        assert c_res.status_code == 201
        c_data = c_res.json()
        compound_row_id = c_data["row_id"]
        version_id = c_data["version"]["id"]

        # 2. Run prediction workflow (which pre-calculates and caches PK parameter sets)
        with SessionLocal() as db:
            run_compound_prediction_workflow(compound_row_id, db)

        # 3. Verify cached retrieval for pk-multi-species is immediate and non-empty
        pk_multi_res = client.get(f"/api/compound-versions/{version_id}/pk-multi-species")
        assert pk_multi_res.status_code == 200
        multi_data = pk_multi_res.json()
        assert "species_profiles" in multi_data
        assert "Rat" in multi_data["species_profiles"]
        assert "Human" in multi_data["species_profiles"]
        assert "Dog" in multi_data["species_profiles"]
        assert "Monkey" in multi_data["species_profiles"]
        assert "Mouse" in multi_data["species_profiles"]

        # Rat CL and V should be computed from IVIVE
        rat_prof = multi_data["species_profiles"]["Rat"]
        assert rat_prof["cl"]["value"] is not None
        assert rat_prof["v"]["value"] is not None
    finally:
        client.delete(f"/api/projects/{project_id}")


def test_experimental_data_auto_recalculation_and_default_precedence():
    # 1. Create project & compound
    p_res = client.post("/api/projects", json={"name": "Exp Precedence Test", "target_name": "Protease"})
    assert p_res.status_code == 201
    project_id = p_res.json()["id"]

    try:
        c_res = client.post(f"/api/projects/{project_id}/compounds", json={
            "compound_id": "CMPD-EXP-01",
            "smiles": "c1ccccc1NC(=O)C"
        })
        assert c_res.status_code == 201
        c_data = c_res.json()
        compound_row_id = c_data["row_id"]
        version_id = c_data["version"]["id"]

        # 2. Run prediction
        with SessionLocal() as db:
            run_compound_prediction_workflow(compound_row_id, db)

        # Initial baseline multi-species
        init_multi = client.get(f"/api/compound-versions/{version_id}/pk-multi-species").json()
        assert init_multi["species_profiles"]["Rat"]["is_experimental"] is False

        # 3. Add an In Vivo IV PK study with NCA results (Measured CL = 15.4 mL/min/kg, Vz = 2.1 L/kg)
        study_res = client.post(f"/api/compounds/{compound_row_id}/pk-studies", json={
            "study_name": "Rat IV 5mpk Study",
            "species": "Rat",
            "route": "IV",
            "dose": 5.0,
            "dose_unit": "mg/kg",
            "matrix": "Plasma"
        })
        assert study_res.status_code == 201
        study_id = study_res.json()["id"]

        # Add observations
        obs_res = client.post(f"/api/pk-studies/{study_id}/observations", json=[
            {"time_raw": 0.08, "time_unit": "h", "concentration_raw": 5000.0, "concentration_unit": "ng/mL"},
            {"time_raw": 0.5, "time_unit": "h", "concentration_raw": 3200.0, "concentration_unit": "ng/mL"},
            {"time_raw": 1.0, "time_unit": "h", "concentration_raw": 1800.0, "concentration_unit": "ng/mL"},
            {"time_raw": 2.0, "time_unit": "h", "concentration_raw": 800.0, "concentration_unit": "ng/mL"},
            {"time_raw": 4.0, "time_unit": "h", "concentration_raw": 150.0, "concentration_unit": "ng/mL"},
        ])
        assert obs_res.status_code == 200

        # Run NCA
        nca_res = client.post(f"/api/pk-studies/{study_id}/run-nca", json={})
        assert nca_res.status_code == 200
        nca_data = nca_res.json()
        assert nca_data["cl"] is not None
        measured_cl = nca_data["cl"]

        # 4. Check that multi-species PK and Foundation profile automatically updated with experimental values by default!
        updated_multi = client.get(f"/api/compound-versions/{version_id}/pk-multi-species").json()
        rat_prof = updated_multi["species_profiles"]["Rat"]
        assert rat_prof["is_experimental"] is True
        assert rat_prof["experimental_notes"] is not None
        assert "실험값 반영" in rat_prof["experimental_notes"]

        # Check PK parameter foundation
        found_res = client.get(f"/api/compound-versions/{version_id}/pk-foundation?species=Rat").json()
        iv_pset = found_res["route_parameter_sets"]["IV"]
        assert iv_pset["cl_source_type"] == "EXPERIMENTAL_NCA"
        assert abs(iv_pset["cl_value"] - measured_cl) < 1e-2
        assert found_res["distribution"]["v_source_type"] in ("EXPERIMENTAL_VZ", "EXPERIMENTAL_VSS")
    finally:
        client.delete(f"/api/projects/{project_id}")
