import pytest

from backend.predictive_hepatic_clearance import predictive_hepatic_clearance


def test_predictive_chain_uses_only_model_outputs_and_records_provenance():
    result = predictive_hepatic_clearance(
        predicted_ppb_percent=90.0,
        predicted_hlm_log10_ml_min_kg=1.0,
        applicability_domain="IN_DOMAIN",
    )
    assert result["prediction_type"] == "MODEL_PREDICTED_CHAIN"
    assert result["provenance"]["observed_inputs_used"] == []
    assert result["provenance"]["ppb_source"] == "MODEL_PREDICTED"
    assert result["value"] > 0


def test_predictive_chain_rejects_unknown_domain():
    with pytest.raises(ValueError, match="unknown applicability domain"):
        predictive_hepatic_clearance(predicted_ppb_percent=90, predicted_hlm_log10_ml_min_kg=1, applicability_domain="UNKNOWN")
