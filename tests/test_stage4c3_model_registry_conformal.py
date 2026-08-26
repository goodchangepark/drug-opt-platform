"""Targeted tests for Stage 4C-3: Model Registry Recovery, Optimization UI Fix, & Conformal Uncertainty."""

import pytest
from fastapi.testclient import TestClient

from backend.admet_predictor import MODEL_SPECS, model_files_available, predict_endpoint
from backend.conformal import (CONFORMAL_CALIBRATION_REGISTRY, compute_calibrated_uncertainty,
                               evaluate_conformal_calibration_coverage)
from backend.main import app


@pytest.fixture
def client():
    return TestClient(app)


def test_model_files_availability():
    """Verify all 18 installed endpoints return available: True."""
    installed_endpoints = [
        "Solubility", "Permeability", "Plasma protein binding",
        "HLM intrinsic clearance", "RLM intrinsic clearance", "MLM intrinsic clearance",
        "CYP1A2 inhibitor", "CYP2C9 inhibitor", "CYP2C19 inhibitor", "CYP2D6 inhibitor", "CYP3A4 inhibitor",
        "CYP2C9 substrate", "CYP2D6 substrate", "CYP3A4 substrate",
        "P-gp inhibitor", "hERG liability", "Ames mutagenicity", "DILI clinical liability",
    ]
    for endpoint in installed_endpoints:
        available, reason = model_files_available(endpoint)
        assert available is True, f"Endpoint {endpoint} failed availability check: {reason}"


def test_real_inference_execution():
    """Run real inference through production code path for key installed models."""
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    for endpoint in ["Solubility", "Plasma protein binding", "hERG liability", "CYP3A4 inhibitor"]:
        res = predict_endpoint(aspirin, endpoint)
        assert res["status"] == "COMPLETE"
        assert "predicted_value" in res
        assert "calibrated_uncertainty" in res
        assert "data_provenance" in res["calibrated_uncertainty"]
        assert "calibration_quality" in res["calibrated_uncertainty"]


def test_conformal_uncertainty_regression():
    """Verify calibrated 90% conformal prediction interval for quantitative endpoints."""
    domain_in = {"classification": "IN_DOMAIN", "similarity": 0.85}
    res = compute_calibrated_uncertainty("HLM intrinsic clearance", 1.50, domain_in, nominal_level="0.90")
    assert res["data_provenance"] == "EXTERNAL"
    assert res["calibration_quality"] == "UNDERCOVERED"
    assert res["display_label"] == "90% Conformal Prediction Interval"
    assert res["lower_bound"] == 0.452
    assert res["upper_bound"] == 2.548
    assert res["interval_width"] == 2.096
    assert res["empirical_coverage"] == 0.796


def test_conformal_uncertainty_classification():
    """Verify conformal prediction set for binary classification endpoints."""
    domain_in = {"classification": "IN_DOMAIN", "similarity": 0.80}
    res = compute_calibrated_uncertainty("hERG liability", 0.9999, domain_in)
    assert res["data_provenance"] == "EXTERNAL"
    assert res["calibration_quality"] == "UNDERCOVERED"
    assert res["prediction_set"] == ["POSITIVE"]

    res_uncertain = compute_calibrated_uncertainty("hERG liability", 0.50, domain_in)
    assert res_uncertain["data_provenance"] == "EXTERNAL"
    assert res_uncertain["prediction_set"] == ["POSITIVE", "NEGATIVE"]
    assert res_uncertain["is_uncertain_set"] is True


def test_ood_conformal_safety_rule():
    """Verify OOD warning is attached to conformal output when compound is OUT_OF_DOMAIN."""
    domain_ood = {"classification": "OUT_OF_DOMAIN", "similarity": 0.15}
    res = compute_calibrated_uncertainty("HLM intrinsic clearance", 1.50, domain_ood)
    assert any("OUT OF DOMAIN — CONFORMAL COVERAGE MAY NOT GENERALIZE" in w for w in res["warnings"])


def test_conformal_unavailable_fallback():
    """Verify uncalibrated endpoints return quality UNAVAILABLE with reason."""
    res = compute_calibrated_uncertainty("Unknown Endpoint", 1.0)
    assert res["calibration_quality"] == "UNAVAILABLE"
    assert "Independent external calibration set not provided" in res["reason"]


def test_dashboard_model_registry_api(client):
    """Verify GET /api/dashboard returns model_registry with correct READY status."""
    response = client.get("/api/dashboard")
    assert response.status_code == 200
    data = response.json()
    assert "model_registry" in data
    registry = {row["endpoint"]: row for row in data["model_registry"]}
    assert registry["Solubility"]["status"] == "READY"
    assert registry["Plasma protein binding"]["status"] == "READY"
    assert registry["hERG liability"]["status"] == "READY"
