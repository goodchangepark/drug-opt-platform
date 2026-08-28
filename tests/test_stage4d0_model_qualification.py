"""
Unit tests for Stage 4D-0: Model Qualification & Ensemble Foundation.
"""

import json
from pathlib import Path
import pytest

from backend.endpoint_contracts import (
    ENDPOINT_CONTRACTS,
    EndpointContract,
    EndpointCategory,
    OutputType,
    Directionality,
    ExecutionTier,
    ARM64Status,
    AdapterPredictionResult,
    BaseModelAdapter,
    check_ensemble_compatibility,
    execute_fault_tolerant_ensemble,
    get_endpoint_contract,
)


def test_endpoint_contracts_registry_completeness():
    """Verify all required primary and translational endpoints exist in registry."""
    required_endpoints = [
        "Solubility",
        "Permeability",
        "Plasma protein binding",
        "HLM intrinsic clearance",
        "RLM intrinsic clearance",
        "MLM intrinsic clearance",
        "CYP1A2 inhibitor",
        "CYP2C9 inhibitor",
        "CYP2C19 inhibitor",
        "CYP2D6 inhibitor",
        "CYP3A4 inhibitor",
        "CYP2C9 substrate",
        "CYP2D6 substrate",
        "CYP3A4 substrate",
        "P-gp inhibitor",
        "hERG liability",
        "Ames mutagenicity",
        "DILI clinical liability",
        "Metabolic soft spots",
        "PK Systemic Clearance",
        "PK Volume of Distribution",
        "PK Bioavailability",
    ]
    for name in required_endpoints:
        contract = get_endpoint_contract(name)
        assert contract is not None, f"Missing authoritative contract for '{name}'"
        assert len(contract.endpoint_id) > 0
        assert len(contract.canonical_unit) > 0
        assert len(contract.scientific_definition) > 0


def test_ensemble_compatibility_success():
    """Verify identical contracts pass compatibility check."""
    contract = get_endpoint_contract("Solubility")
    is_compat, reason = check_ensemble_compatibility(contract, contract)
    assert is_compat is True
    assert "Ensemble compatible" in reason


def test_ensemble_compatibility_species_mismatch():
    """Verify cross-species mixing (HLM vs RLM) is rejected."""
    hlm = get_endpoint_contract("HLM intrinsic clearance")
    rlm = get_endpoint_contract("RLM intrinsic clearance")
    is_compat, reason = check_ensemble_compatibility(hlm, rlm)
    assert is_compat is False
    assert "Incompatible" in reason


def test_ensemble_compatibility_role_mismatch():
    """Verify inhibitor vs substrate mixing is rejected."""
    cyp_inh = get_endpoint_contract("CYP3A4 inhibitor")
    cyp_sub = get_endpoint_contract("CYP3A4 substrate")
    is_compat, reason = check_ensemble_compatibility(cyp_inh, cyp_sub)
    assert is_compat is False
    assert "Incompatible" in reason


def test_ensemble_compatibility_unit_mismatch():
    """Verify mismatched canonical units fail."""
    sol = get_endpoint_contract("Solubility")
    ppb = get_endpoint_contract("Plasma protein binding")
    is_compat, reason = check_ensemble_compatibility(sol, ppb)
    assert is_compat is False
    assert "Incompatible" in reason


class MockSuccessAdapter(BaseModelAdapter):
    def __init__(self, model_id: str, value: float):
        self.model_id = model_id
        self.model_version = "v1.0"
        self.model_family = "mock_family"
        self.supported_endpoints = {"solubility_aqueous_logs"}
        self.execution_tier = ExecutionTier.TIER_1_LOCAL_FAST
        self.arm64_status = ARM64Status.RUNS_LOCAL_ARM64
        self.value = value

    def is_available(self):
        return True, ""

    def predict(self, canonical_smiles, endpoint_contract, compound_metadata=None):
        return AdapterPredictionResult(
            model_id=self.model_id,
            model_version=self.model_version,
            model_family=self.model_family,
            endpoint_id=endpoint_contract.endpoint_id,
            canonical_unit=endpoint_contract.canonical_unit,
            value=self.value,
            applicability_domain="IN_DOMAIN",
            confidence="HIGH",
        )


class MockFailingAdapter(BaseModelAdapter):
    def __init__(self, model_id: str):
        self.model_id = model_id
        self.model_version = "v1.0"
        self.model_family = "failing_family"
        self.supported_endpoints = {"solubility_aqueous_logs"}
        self.execution_tier = ExecutionTier.TIER_1_LOCAL_FAST
        self.arm64_status = ARM64Status.RUNS_LOCAL_ARM64

    def is_available(self):
        return True, ""

    def predict(self, canonical_smiles, endpoint_contract, compound_metadata=None):
        raise RuntimeError("Synthetic simulated hardware numerical exception")


def test_fault_tolerant_failure_isolation():
    """Verify that 1 failing model does not crash the ensemble."""
    adapter_1 = MockSuccessAdapter("model_1", -2.5)
    adapter_2 = MockSuccessAdapter("model_2", -2.3)
    adapter_failing = MockFailingAdapter("model_3")

    contract = get_endpoint_contract("Solubility")
    result = execute_fault_tolerant_ensemble(
        adapters=[adapter_1, adapter_failing, adapter_2],
        canonical_smiles="CC(=O)Oc1ccccc1C(=O)O",
        contract=contract,
    )

    assert result.is_valid is True
    assert result.member_count == 2
    assert len(result.successful_predictions) == 2
    assert len(result.failed_models) == 1
    assert result.failed_models[0]["model_id"] == "model_3"
    assert "Synthetic simulated" in result.failed_models[0]["reason"]
    assert "SUCCESS (2/3 models passed)" in result.status_summary


def test_stage4d0_validation_artifacts_exist_and_valid():
    """Verify all 6 machine-readable JSON validation files exist and are valid JSON."""
    validation_dir = Path("validation")
    required_files = [
        "stage4d0_current_model_inventory.json",
        "stage4d0_candidate_model_qualification.json",
        "stage4d0_dataset_lineage.json",
        "stage4d0_model_diversity.json",
        "stage4d0_arm64_compatibility.json",
        "stage4d0_endpoint_contracts.json",
    ]

    for fname in required_files:
        fpath = validation_dir / fname
        assert fpath.is_file(), f"Missing required validation artifact: {fname}"
        with open(fpath, encoding="utf-8") as f:
            data = json.load(f)
            assert isinstance(data, (dict, list)), f"Invalid JSON structure in {fname}"


def test_stage4d0_documentation_files_exist():
    """Verify all 5 required markdown documentation files exist in docs/."""
    docs_dir = Path("docs")
    required_docs = [
        "stage4d0-model-qualification.md",
        "stage4d0-endpoint-contracts.md",
        "stage4d0-ensemble-architecture.md",
        "stage4d0-license-audit.md",
        "stage4d0-candidate-models.md",
    ]

    for fname in required_docs:
        fpath = docs_dir / fname
        assert fpath.is_file(), f"Missing required doc file: {fname}"
        content = fpath.read_text(encoding="utf-8")
        assert len(content) > 500, f"Doc file {fname} is unexpectedly short"
