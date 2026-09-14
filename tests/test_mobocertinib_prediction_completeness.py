from __future__ import annotations

from backend.current_production_executor import (
    EXECUTION_CONTRACT,
    execute_current_production_quantitative,
)
from backend.model_artifact_authority import (
    artifact_bundle_sha256,
    model_artifact_registration,
    validation_registration_error,
)
from backend.predict_all_core_contract import (
    CONTEXT_REQUIRED_ENDPOINTS,
    CURRENT_DATA_CEILING_ENDPOINTS,
    MECHANISTIC_EXECUTION_ENDPOINTS,
    MODEL_NOT_REGISTERED_ENDPOINTS,
    MODEL_UNAVAILABLE_ENDPOINTS,
    PREDICT_ALL_EXECUTION_ENDPOINTS,
    PUBLISHABLE_CURRENT_ENDPOINTS,
)
from backend.prediction_orchestrator import PREDICT_ALL_MODEL_UNAVAILABLE


MOBOCERTINIB = "C=CC(=O)Nc1cc(Nc2ncc(C(=O)OC(C)C)c(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C"


def test_59_endpoint_contract_is_complete_and_disjoint():
    groups = (
        PUBLISHABLE_CURRENT_ENDPOINTS,
        MECHANISTIC_EXECUTION_ENDPOINTS,
        MODEL_NOT_REGISTERED_ENDPOINTS,
        MODEL_UNAVAILABLE_ENDPOINTS,
        CURRENT_DATA_CEILING_ENDPOINTS,
        CONTEXT_REQUIRED_ENDPOINTS,
    )
    assert sum(len(group) for group in groups) == 59
    assert len(set().union(*groups)) == 59
    assert len(PUBLISHABLE_CURRENT_ENDPOINTS) == 32
    assert len(PREDICT_ALL_EXECUTION_ENDPOINTS) == 36
    assert "Ames mutagenicity" in PREDICT_ALL_MODEL_UNAVAILABLE


def test_mobocertinib_executes_every_frozen_quantitative_route():
    outputs = execute_current_production_quantitative(MOBOCERTINIB)
    assert len(outputs) == 9
    for endpoint, output in outputs.items():
        assert output["execution_contract"] == EXECUTION_CONTRACT
        assert output["execution_status"] == "SUCCESS", endpoint
        assert output["production_prediction"] is not None
        assert output["components"]
        assert all(row["status"] == "SUCCESS" for row in output["components"])


def test_publishable_routes_have_truthful_installed_artifact_registration():
    # The other publishable endpoints are deterministic properties or existing
    # orchestrator classifiers; every one must still resolve an immutable file
    # bundle through the same fail-closed authority.
    for endpoint in sorted(PUBLISHABLE_CURRENT_ENDPOINTS):
        registration = model_artifact_registration(endpoint)
        assert registration is not None, endpoint
        assert validation_registration_error(registration) is None, endpoint
        assert artifact_bundle_sha256(registration), endpoint
        assert all(path.is_file() for path in registration.resolved_paths()), endpoint


def test_hia_and_bbb_are_not_promoted_from_route_labels_alone():
    assert model_artifact_registration("HIA") is None
    assert model_artifact_registration("BBB_PENETRATION") is None
