import pytest
from fastapi.testclient import TestClient
from backend.main import app
from backend.database import SessionLocal
from backend.models import Project, Compound, CompoundVersion
from backend.ivive import get_multi_species_pk_profile
from backend.simulation import run_pk_simulation, PKSimulationRequest

client = TestClient(app)


def test_dog_and_monkey_metabolism_unavailable_status():
    """Verify Dog and Monkey microsomal models are explicitly MODEL_UNAVAILABLE (Case C)."""
    with SessionLocal() as db:
        v = db.query(CompoundVersion).first()
        assert v is not None

        # Test multi-species PK profile Dog and Monkey CLint status
        res = get_multi_species_pk_profile(db, v.id)
        assert "Dog" in res["species_profiles"]
        assert "Monkey" in res["species_profiles"]

        dog = res["species_profiles"]["Dog"]
        monkey = res["species_profiles"]["Monkey"]

        # CL must be None and source MODEL_UNAVAILABLE without fabrication
        assert dog["cl"]["value"] is None
        assert dog["cl"]["source"] == "MODEL_UNAVAILABLE"
        assert monkey["cl"]["value"] is None
        assert monkey["cl"]["source"] == "MODEL_UNAVAILABLE"


def test_multi_species_pk_assembly_all_five_species():
    """Verify Mouse, Rat, Dog, Monkey, and Human are assembled with independent provenance."""
    with SessionLocal() as db:
        v = db.query(CompoundVersion).first()
        res = get_multi_species_pk_profile(db, v.id)

        sp_map = res["species_profiles"]
        assert set(sp_map.keys()) == {"Mouse", "Rat", "Dog", "Monkey", "Human"}

        # Human, Rat, Mouse must have valid IVIVE clearance
        assert sp_map["Human"]["readiness"] == "READY"
        assert sp_map["Human"]["cl"]["value"] is not None
        assert sp_map["Human"]["cl"]["source"] == "HEPATIC_IVIVE"
        assert sp_map["Human"]["v"]["value"] is not None
        assert sp_map["Human"]["t_half_hours"] is not None

        assert sp_map["Rat"]["readiness"] == "READY"
        assert sp_map["Rat"]["cl"]["value"] is not None
        assert sp_map["Rat"]["cl"]["source"] == "HEPATIC_IVIVE"

        assert sp_map["Mouse"]["readiness"] == "READY"
        assert sp_map["Mouse"]["cl"]["value"] is not None
        assert sp_map["Mouse"]["cl"]["source"] == "HEPATIC_IVIVE"


def test_pk_simulation_normalized_dose_and_routes():
    """Verify PK simulation runs for IV and PO with normalized 1.0 mg/kg dose."""
    with SessionLocal() as db:
        v = db.query(CompoundVersion).first()

        # IV simulation at 1.0 mg/kg
        req_iv = PKSimulationRequest(
            species="Human",
            route="IV",
            administration_type="IV_BOLUS",
            dose=1.0,
            dose_unit="mg/kg",
            infusion_duration_hours=0.0,
            dosing_frequency="Single Dose",
            dose_interval_hours=24.0,
            num_doses=1,
            model_type="ONE_COMPARTMENT"
        )
        sim_iv = run_pk_simulation(db, v.id, req_iv)
        assert sim_iv.output_metrics is not None
        assert sim_iv.output_metrics["cmax_ng_ml"] > 0
        assert sim_iv.output_metrics["auc_inf_analytical_ng_h_ml"] > 0
        assert len(sim_iv.time_series) > 50

        # PO simulation with F and ka override
        req_po = PKSimulationRequest(
            species="Rat",
            route="PO",
            administration_type="EXTRAVASCULAR_1COMP",
            dose=1.0,
            dose_unit="mg/kg",
            infusion_duration_hours=0.0,
            dosing_frequency="Single Dose",
            dose_interval_hours=24.0,
            num_doses=1,
            model_type="ONE_COMPARTMENT",
            user_f_override=75.0,
            user_ka_override=1.2
        )
        sim_po = run_pk_simulation(db, v.id, req_po)
        assert sim_po.output_metrics is not None
        assert sim_po.output_metrics["cmax_ng_ml"] > 0
        assert sim_po.output_metrics["tmax_hours"] > 0
        assert sim_po.output_metrics["auc_inf_analytical_ng_h_ml"] > 0
        assert len(sim_po.time_series) > 50


def test_two_compound_comparison_api_expanded_fields():
    """Verify two-compound comparison API includes ADME, Metabolism, and PK metrics."""
    with SessionLocal() as db:
        proj = Project(name="Comparison Test Project Dedicated", target="EGFR", molecule_type="Small Molecule")
        db.add(proj)
        db.commit()
        db.refresh(proj)
        
        c1 = Compound(project_id=proj.id, compound_id="CMPD_COMP_1", name="Gefitinib", status="STRUCTURE_READY")
        c2 = Compound(project_id=proj.id, compound_id="CMPD_COMP_2", name="Erlotinib", status="STRUCTURE_READY")
        db.add_all([c1, c2])
        db.commit()
        
        v1 = CompoundVersion(
            compound_row_id=c1.id,
            version_number=1,
            original_smiles="COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
            canonical_smiles="COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
            isomeric_smiles="COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
            inchi="InChI=1S/C22H24ClFN4O3",
            inchikey="XGALLCVXEZPNRV-UHFFFAOYSA-N"
        )
        v2 = CompoundVersion(
            compound_row_id=c2.id,
            version_number=1,
            original_smiles="COCCOc1cc2c(cc1OCCOC)ncnc2Nc1cccc(c1)C#C",
            canonical_smiles="COCCOc1cc2c(cc1OCCOC)ncnc2Nc1cccc(c1)C#C",
            isomeric_smiles="COCCOc1cc2c(cc1OCCOC)ncnc2Nc1cccc(c1)C#C",
            inchi="InChI=1S/C22H23N3O4",
            inchikey="AAHAUSGZVGZGBF-UHFFFAOYSA-N"
        )
        db.add_all([v1, v2])
        db.commit()
        
        from backend.main import run_compound_prediction_workflow
        run_compound_prediction_workflow(c1.id, db)
        run_compound_prediction_workflow(c2.id, db)
        
        c1_id, c2_id, proj_id = c1.id, c2.id, proj.id

    try:
        resp = client.get(f"/api/projects/{proj_id}/compare?ids={c1_id},{c2_id}")
        assert resp.status_code == 200
        data = resp.json()

        # Check metrics
        metrics = data["metrics"]
        # Properties
        assert "MW" in metrics and "cLogP" in metrics and "TPSA" in metrics and "QED" in metrics
        # ADME
        assert "Solubility" in metrics and "Caco-2" in metrics and "PPB" in metrics and "fu" in metrics
        # Metabolism
        assert "HLM" in metrics and "RLM" in metrics and "MLM" in metrics and "DLM" in metrics and "CyLM" in metrics
        assert "CYP3A4 Inh" in metrics and "Soft Spots" in metrics
        # PK
        assert "Mouse CL (IV)" in metrics and "Rat CL (IV)" in metrics and "Human CL (IVIVE)" in metrics
        assert "Human Vd (pred)" in metrics and "Human t1/2 (pred)" in metrics and "Human AUC (1mg/kg IV)" in metrics

        # Check compound entries
        for cmpd in data["compounds"]:
            assert cmpd["MW"] is not None
            assert cmpd["cLogP"] is not None
            assert cmpd["sources"]["DLM"] == "MODEL_UNAVAILABLE"
            assert cmpd["sources"]["CyLM"] == "MODEL_UNAVAILABLE"
            assert cmpd["Rat CL (IV)"] is not None
            assert cmpd["Human CL (IVIVE)"] is not None
    finally:
        with SessionLocal() as db:
            p_del = db.get(Project, proj_id)
            if p_del:
                db.delete(p_del)
                db.commit()
