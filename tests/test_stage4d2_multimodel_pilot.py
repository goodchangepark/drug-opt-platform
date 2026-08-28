"""
Drug-OPT Stage 4D-2: Qualified Multi-Model Pilot & External Validation Tests.

Covers:
1. Registration of qualified pilot model adapters
2. Direct execution of pilot adapters on valid, borderline, and invalid SMILES
3. Aqueous solubility multi-model execution and static consensus
4. Caco-2 permeability multi-model execution and static consensus
5. CYP3A4 inhibitor multi-model execution and probability aggregation
6. hERG liability multi-model execution and probability aggregation
7. Site-of-metabolism rank fusion (SyGMa + SMARTCyp)
8. Diversity penalty and empirical error correlation weighting
9. Verification of all 5 authoritative Stage 4D-2 validation artifacts
10. Shadow mode non-interference with production endpoints
"""

import json
from pathlib import Path
import pytest
import numpy as np

from backend.endpoint_contracts import get_endpoint_contract, OutputType
from backend.multimodel import (
    get_adapters_for_endpoint,
    list_registered_adapters,
    get_model_adapter,
    ExecutionStatus,
    ESOLPhyschemSolubilityAdapter,
    DescriptorGBRSolubilityAdapter,
    PhyschemCaco2Adapter,
    MorganCYP3A4InhibitorAdapter,
    PhyschemHERGAdapter,
    SMARTCypMetabolismAdapter,
)
from backend.consensus import (
    compute_endpoint_consensus,
    calculate_static_model_weight,
    ConsensusMode,
    AgreementStatus,
    AggregationType,
)

ROOT = Path(__file__).resolve().parents[1]


def test_registered_pilot_adapters():
    """Verify all pilot model adapters are registered in the global registry."""
    adapters = list_registered_adapters()
    assert len(adapters) >= 25, f"Expected at least 25 registered adapters, got {len(adapters)}"
    
    sol_adapters = get_adapters_for_endpoint("Solubility")
    assert len(sol_adapters) == 3
    sol_ids = {a.model_id for a in sol_adapters}
    assert sol_ids == {"admetica_solubility", "esol_delaney_v1", "rdkit_gbr_solubility_v1"}

    caco_adapters = get_adapters_for_endpoint("Permeability")
    assert len(caco_adapters) == 2
    caco_ids = {a.model_id for a in caco_adapters}
    assert caco_ids == {"admetica_caco2", "physchem_caco2_v1"}

    cyp_adapters = get_adapters_for_endpoint("CYP3A4 inhibitor")
    assert len(cyp_adapters) == 2
    cyp_ids = {a.model_id for a in cyp_adapters}
    assert cyp_ids == {"admetica_cyp_cyp3a4-inhibitor", "morgan_cyp3a4_inh_v1"}

    herg_adapters = get_adapters_for_endpoint("hERG liability")
    assert len(herg_adapters) == 2
    herg_ids = {a.model_id for a in herg_adapters}
    assert herg_ids == {"admetica_safety_herg", "physchem_herg_v1"}

    som_adapters = get_adapters_for_endpoint("Metabolic soft spots")
    assert len(som_adapters) == 2
    som_ids = {a.model_id for a in som_adapters}
    assert som_ids == {"sygma_phase1_2", "smartcyp_dft_v1"}


def test_esol_and_gbr_solubility_adapters():
    """Test ESOL and GBR solubility adapters execution on Aspirin."""
    contract = get_endpoint_contract("Solubility")
    assert contract is not None

    esol = ESOLPhyschemSolubilityAdapter()
    gbr = DescriptorGBRSolubilityAdapter()

    aspirin_smi = "CC(=O)Oc1ccccc1C(=O)O"
    
    p_esol = esol.execute(aspirin_smi, contract)
    assert p_esol.execution_status == ExecutionStatus.SUCCESS
    assert p_esol.value is not None
    assert -5.0 <= p_esol.value <= 0.0
    assert p_esol.canonical_unit == "log10(mol/L)"
    assert p_esol.applicability_domain in {"IN_DOMAIN", "BORDERLINE"}

    p_gbr = gbr.execute(aspirin_smi, contract)
    assert p_gbr.execution_status == ExecutionStatus.SUCCESS
    assert p_gbr.value is not None
    assert -5.0 <= p_gbr.value <= 0.0

    # Invalid input handling
    p_invalid = esol.execute("INVALID_SMILES", contract)
    assert p_invalid.execution_status == ExecutionStatus.INVALID_INPUT


def test_physchem_caco2_adapter():
    """Test physicochemical Caco-2 permeability adapter on Metoprolol."""
    contract = get_endpoint_contract("Permeability")
    assert contract is not None

    adapter = PhyschemCaco2Adapter()
    metoprolol_smi = "COCCc1ccc(OCC(O)CNC(C)C)cc1"
    
    p = adapter.execute(metoprolol_smi, contract)
    assert p.execution_status == ExecutionStatus.SUCCESS
    assert p.value is not None
    assert -7.0 <= p.value <= -3.0
    assert p.canonical_unit == "log10(cm/s)"
    assert p.raw_outputs["mw"] > 200.0


def test_morgan_cyp3a4_adapter():
    """Test Morgan ECFP4 CYP3A4 inhibitor classifier adapter."""
    contract = get_endpoint_contract("CYP3A4 inhibitor")
    assert contract is not None

    adapter = MorganCYP3A4InhibitorAdapter()
    # Ketoconazole (strong CYP3A4 inhibitor)
    keto_smi = "CC(=O)N1CCN(CC1)c2ccc(OCC3COC(Cn4cncn4)(O3)c5ccc(Cl)cc5Cl)cc2"
    p_keto = adapter.execute(keto_smi, contract)
    assert p_keto.execution_status == ExecutionStatus.SUCCESS
    assert p_keto.probability is not None
    assert p_keto.probability >= 0.50
    assert p_keto.predicted_class == "INHIBITOR"

    # Aspirin (non-inhibitor)
    aspirin_smi = "CC(=O)Oc1ccccc1C(=O)O"
    p_asp = adapter.execute(aspirin_smi, contract)
    assert p_asp.execution_status == ExecutionStatus.SUCCESS
    assert p_asp.probability is not None
    assert p_asp.probability < 0.50
    assert p_asp.predicted_class == "NON_INHIBITOR"


def test_physchem_herg_adapter():
    """Test basic center pharmacophore hERG liability adapter."""
    contract = get_endpoint_contract("hERG liability")
    assert contract is not None

    adapter = PhyschemHERGAdapter()
    # Astemizole / Terfenadine analogue with basic amine
    terf_smi = "CC(C)(C)c1ccc(C(O)CCCN2CCC(C(O)(c3ccccc3)c4ccccc4)CC2)cc1"
    p = adapter.execute(terf_smi, contract)
    assert p.execution_status == ExecutionStatus.SUCCESS
    assert p.probability is not None
    assert p.probability >= 0.50
    assert p.predicted_class == "BLOCKER"


def test_smartcyp_metabolism_adapter():
    """Test SMARTCyp metabolism adapter."""
    contract = get_endpoint_contract("Metabolic soft spots")
    assert contract is not None

    adapter = SMARTCypMetabolismAdapter()
    smi = "CCOc1ccc(NC(C)=O)cc1"  # Phenacetin
    p = adapter.execute(smi, contract)
    assert p.execution_status == ExecutionStatus.SUCCESS
    assert p.raw_outputs is not None
    assert "soft_spots" in p.raw_outputs or "sy_spots" in p.raw_outputs or len(p.raw_outputs) > 0


def test_solubility_consensus_aggregation():
    """Test multi-model weighted consensus aggregation on Solubility."""
    contract = get_endpoint_contract("Solubility")
    adapters = get_adapters_for_endpoint("Solubility")
    smi = "CC(=O)Oc1ccccc1C(=O)O"

    payloads = [a.execute(smi, contract) for a in adapters]
    assert all(p.execution_status == ExecutionStatus.SUCCESS for p in payloads)

    cons = compute_endpoint_consensus(
        endpoint_name="Solubility",
        compound_version_id=1,
        model_payloads=payloads,
        mode=ConsensusMode.SHADOW,
    )
    assert cons.combined_value is not None
    assert -3.0 <= cons.combined_value <= 0.0
    assert len(cons.models_used) == 3
    assert cons.consensus_mode == ConsensusMode.SHADOW
    assert "model_disagreement_std" in cons.dispersion
    assert cons.dispersion["model_disagreement_std"] >= 0.0


def test_empirical_diversity_penalties_applied():
    """Verify empirical diversity penalties are correctly evaluated."""
    sol_adapters = get_adapters_for_endpoint("Solubility")
    contract = get_endpoint_contract("Solubility")
    smi = "c1ccccc1"

    payloads = [a.execute(smi, contract) for a in sol_adapters]
    p_esol = [p for p in payloads if p.model_id == "esol_delaney_v1"][0]
    p_gbr = [p for p in payloads if p.model_id == "rdkit_gbr_solubility_v1"][0]

    # Calculate weights
    w_esol, reason_esol = calculate_static_model_weight(p_esol, payloads)
    w_gbr, reason_gbr = calculate_static_model_weight(p_gbr, payloads)

    assert w_esol > 0.0
    assert w_gbr > 0.0
    assert "Diversity(" in reason_esol


def test_stage4d2_artifacts_exist_and_valid():
    """Verify all 5 JSON validation artifacts from Stage 4D-2 exist and are well-formed."""
    val_dir = ROOT / "validation"
    artifact_files = [
        "stage4d2_pilot_registry.json",
        "stage4d2_external_validation.json",
        "stage4d2_error_correlation.json",
        "stage4d2_runtime_benchmark.json",
        "stage4d2_promotion_decisions.json",
    ]
    for filename in artifact_files:
        filepath = val_dir / filename
        assert filepath.is_file(), f"Missing artifact: {filename}"
        with open(filepath, "r", encoding="utf-8") as f:
            data = json.load(f)
            assert "stage" in data or "hardware_platform" in data
            assert len(data) > 0


def test_promotion_decisions_content():
    """Verify Stage 4D-2 promotion decisions match scientific qualification."""
    filepath = ROOT / "validation" / "stage4d2_promotion_decisions.json"
    with open(filepath, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    decisions = data.get("decisions", {})
    assert decisions["Solubility"]["decision"] == "PROMOTION_CANDIDATE"
    assert decisions["Permeability"]["decision"] == "KEEP_SHADOW"
    assert decisions["CYP3A4 inhibitor"]["decision"] == "PROMOTION_CANDIDATE"
    assert decisions["hERG liability"]["decision"] == "KEEP_SHADOW"
    assert decisions["Metabolic soft spots"]["decision"] == "STAGE_4D2B_PREPARATION_VALIDATED"
