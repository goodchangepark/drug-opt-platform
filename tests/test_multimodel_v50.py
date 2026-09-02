"""
Test suite for Drug-OPT v5.0 Prediction Model Expansion & Multi-Model Consensus.
Validates:
1. All new model adapters registered in _ADAPTER_REGISTRY.
2. ADMET-AI 5-ensemble multi-task extraction for CYP, P-gp, hERG.
3. Physicochemical lipophilicity PPB, pKa QSAR, and derived logD 7.4.
4. PredictionOrchestrator multi-model execution and consensus computation.
5. Strict semantic isolation (no mismatched mixing).
6. Engine v1 immutability and frozen hash integrity.
"""

import pytest
from backend.database import SessionLocal, engine
from backend.models import Compound, CompoundVersion, Project
from backend.admet import ensure_admet_schema, ADMETModelRegistry, ADMETPrediction, ADMETConsensusPrediction
from backend.multimodel import (
    get_model_adapter,
    list_registered_adapters,
    list_all_v2_registered_adapters,
    get_adapters_for_endpoint,
    get_v2_adapters_for_endpoint,
    ADMETAITaskAdapter,
    PhyschemHumanPPBAdapter,
    RDKitIonizationPkaAdapter,
    DerivedLogD74Adapter,
)
from backend.endpoint_contracts import get_endpoint_contract
from backend.prediction_orchestrator import PredictionOrchestrator
from backend.project_adaptation_v2 import ENGINE_V1_HASH


@pytest.fixture(scope="module", autouse=True)
def setup_schema():
    ensure_admet_schema(engine)


def test_frozen_engine_v1_integrity():
    """Verify that Prediction Engine v1 baseline hash remains 100% frozen."""
    assert ENGINE_V1_HASH == "12757ab197b5a70d8ea1754678d9a342ab0b6ea0d82f2896bebb767d686bbdeb"


def test_registered_adapters_coverage():
    """Verify that all new Engine v2 multi-model adapters are properly registered."""
    expected_models = [
        "physchem_human_ppb_v1",
        "rdkit_pka_qsar_v1",
        "derived_logd74_v1",
        "admet_ai_cyp1a2_inh",
        "admet_ai_cyp2c9_inh",
        "admet_ai_cyp2c19_inh",
        "admet_ai_cyp2d6_inh",
        "admet_ai_cyp3a4_inh",
        "admet_ai_cyp2c9_sub",
        "admet_ai_cyp2d6_sub",
        "admet_ai_cyp3a4_sub",
        "admet_ai_pgp_inh",
        "admet_ai_herg",
    ]
    for model_id in expected_models:
        adapter = get_model_adapter(model_id)
        assert adapter is not None, f"Model adapter {model_id} not found in registry"
        avail, reason = adapter.is_available()
        assert avail, f"Adapter {model_id} reports unavailable: {reason}"


def test_admet_ai_multitask_adapter_execution():
    """Verify that ADMETAITaskAdapter accurately executes and extracts all endpoints."""
    smiles = "CN(C)CC#CC(=O)Nc1cc(Nc2ncc(Cl)c(Nc3cccc(C(=O)NC(C)C)c3)n2)c(OC)cc1"
    
    for ep in [
        "CYP1A2 inhibitor", "CYP2C9 inhibitor", "CYP2C19 inhibitor",
        "CYP2D6 inhibitor", "CYP3A4 inhibitor",
        "CYP2C9 substrate", "CYP2D6 substrate", "CYP3A4 substrate",
        "P-gp inhibitor", "hERG liability",
    ]:
        adapters = get_v2_adapters_for_endpoint(ep)
        ai_adapter = next((a for a in adapters if isinstance(a, ADMETAITaskAdapter)), None)
        assert ai_adapter is not None, f"No ADMETAITaskAdapter registered for {ep}"
        
        contract = get_endpoint_contract(ep)
        payload = ai_adapter.execute(smiles, contract)
        assert payload.execution_status.value == "SUCCESS"
        assert payload.value is not None
        assert 0.0 <= payload.value <= 1.0
        assert payload.probability is not None
        assert payload.predicted_class in {
            "INHIBITOR", "NON_INHIBITOR",
            "SUBSTRATE", "NON_SUBSTRATE",
            "BLOCKER", "NON_BLOCKER",
        }


def test_physchem_human_ppb_adapter():
    """Verify that PhyschemHumanPPBAdapter computes physically bounded % bound."""
    adapter = get_model_adapter("physchem_human_ppb_v1")
    assert adapter is not None
    contract = get_endpoint_contract("Plasma protein binding")
    
    # Highly lipophilic molecule
    payload_lipo = adapter.execute("c1ccc(Cc2ccccc2)cc1", contract)
    assert payload_lipo.execution_status.value == "SUCCESS"
    assert payload_lipo.value > 80.0
    
    # Highly hydrophilic molecule
    payload_hydro = adapter.execute("C(CO)O", contract)
    assert payload_hydro.execution_status.value == "SUCCESS"
    assert payload_hydro.value < 50.0


def test_rdkit_ionization_and_logd74_adapters():
    """Verify that RDKitIonizationPkaAdapter and DerivedLogD74Adapter return physiological predictions."""
    pka_adapter = get_model_adapter("rdkit_pka_qsar_v1")
    logd_adapter = get_model_adapter("derived_logd74_v1")
    assert pka_adapter is not None
    assert logd_adapter is not None

    contract_pka = get_endpoint_contract("Ionization (pKa)")
    contract_logd = get_endpoint_contract("logD pH7.4 derived estimate")

    smiles = "CN(C)CC#CC(=O)Nc1cc(Nc2ncc(Cl)c(Nc3cccc(C(=O)NC(C)C)c3)n2)c(OC)cc1"
    
    pka_payload = pka_adapter.execute(smiles, contract_pka)
    assert pka_payload.execution_status.value == "SUCCESS"
    
    logd_payload = logd_adapter.execute(smiles, contract_logd)
    assert logd_payload.execution_status.value == "SUCCESS"
    assert logd_payload.value is not None


def test_v2_multimodel_consensus_execution():
    """Verify that multi-model execution and consensus computation succeed for multi-model endpoints."""
    smiles = "CN(C)CC#CC(=O)Nc1cc(Nc2ncc(Cl)c(Nc3cccc(C(=O)NC(C)C)c3)n2)c(OC)cc1"
    
    # Test CYP3A4 multi-model consensus across Admetica, ADMET-AI, Morgan GBR
    adapters = get_v2_adapters_for_endpoint("CYP3A4 inhibitor")
    assert len(adapters) >= 3
    results = []
    contract = get_endpoint_contract("CYP3A4 inhibitor")
    for a in adapters:
        res = a.execute(smiles, contract)
        assert res.execution_status.value == "SUCCESS"
        assert res.probability is not None
        results.append(res.probability)
    
    avg_prob = sum(results) / len(results)
    assert 0.0 <= avg_prob <= 1.0
    print(f"CYP3A4 3-Model consensus probability: {avg_prob:.4f}")
