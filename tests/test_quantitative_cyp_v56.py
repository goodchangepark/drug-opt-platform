"""Unit and integration tests for Quantitative CYP Model v5.6."""
import pytest
from backend.openadmet_cyp import (
    predict_chemeleon_cyp_pic50,
    pic50_to_ic50_um,
    pic50_to_ic50_nm,
    ic50_nm_to_pic50,
    ic50_um_to_pic50,
    OPENADMET_CYP_BENCHMARKS,
    SUPPORTED_CYP_ISOFORMS,
)
from backend.multimodel import get_v2_adapters_for_endpoint
from backend.endpoint_model_validation import build_dmpk_quantitative_expansion_report
from backend.database import SessionLocal
from backend.models import Compound
from backend.endpoint_comparison import build_endpoint_comparison


def test_pic50_unit_conversions():
    """Verify exact logarithmic conversions between pIC50 and molar/micromolar/nanomolar."""
    # 6.0 pIC50 = 1.0 uM = 1000 nM
    assert pic50_to_ic50_um(6.0) == pytest.approx(1.0)
    assert pic50_to_ic50_nm(6.0) == pytest.approx(1000.0)
    assert ic50_nm_to_pic50(1000.0) == pytest.approx(6.0)
    assert ic50_um_to_pic50(1.0) == pytest.approx(6.0)

    # 9.0 pIC50 = 1 nM = 0.001 uM
    assert pic50_to_ic50_nm(9.0) == pytest.approx(1.0)
    assert pic50_to_ic50_um(9.0) == pytest.approx(0.001)


def test_chemeleon_cyp_predictions_and_isoforms():
    """Verify CheMeleon predicts CYP1A2, CYP2C9, CYP2D6, CYP3A4 and excludes CYP2C19."""
    smi = "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C"  # Imatinib
    for iso in ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]:
        pred = predict_chemeleon_cyp_pic50(smi, iso)
        assert 3.0 <= pred.pic50 <= 10.0
        assert pred.ic50_um > 0
        assert pred.ic50_nm > 0
        assert pred.applicability_domain in ("IN_DOMAIN", "BORDERLINE", "OUT_OF_DOMAIN")
        assert pred.provenance["status"] == "CANDIDATE_EXTERNAL_MODEL"

    with pytest.raises(ValueError, match="not supported"):
        predict_chemeleon_cyp_pic50(smi, "CYP2C19")


def test_chemeleon_adapters_registered_as_candidate_v2():
    """Verify OpenADMET CheMeleon adapters are registered in v2 expansion registry."""
    adapters_3a4 = get_v2_adapters_for_endpoint("CYP3A4 quantitative inhibition")
    mids = [a.model_id for a in adapters_3a4]
    assert "openadmet_chemeleon_cyp3a4_pic50" in mids

    adapters_1a2 = get_v2_adapters_for_endpoint("CYP1A2 quantitative inhibition")
    mids_1a2 = [a.model_id for a in adapters_1a2]
    assert "openadmet_chemeleon_cyp1a2_pic50" in mids_1a2


def test_scientific_rows_separate_quantitative_and_classifier():
    """Verify endpoint comparison displays separate quantitative pIC50 and classifier without mixing."""
    db = SessionLocal()
    try:
        comp = db.query(Compound).filter(Compound.name.ilike("%Orforglipron%")).first()
        assert comp is not None
        res = build_endpoint_comparison(db, comp.versions[-1].id)

        cyp_rows = [r for r in res.get("scientific_rows", []) if "CYP" in r["canonical_endpoint"]]
        assert len(cyp_rows) > 0

        for r in cyp_rows:
            if r["canonical_endpoint"] in ("CYP1A2_INHIBITION", "CYP2C9_INHIBITION", "CYP2D6_INHIBITION", "CYP3A4_INHIBITION"):
                qp = r.get("quantitative_prediction")
                assert qp is not None
                assert "pic50" in qp
                assert "ic50_um" in qp
                assert qp["status"] == "CANDIDATE_EXTERNAL_MODEL"

            if r["canonical_endpoint"] == "CYP2C19_INHIBITION":
                assert r.get("quantitative_prediction") is None
    finally:
        db.close()


def test_dmpk_expansion_report_v56():
    """Verify quantitative expansion table contains CheMeleon candidate models and scaffold metrics."""
    report = build_dmpk_quantitative_expansion_report()
    cyp3a4 = next((r for r in report if r["endpoint"] == "CYP3A4 quantitative inhibition"), None)
    assert cyp3a4 is not None
    assert cyp3a4["n"] >= 1
    assert cyp3a4["quantitative_model"] == "OpenADMET CheMeleon CYP3A4 pIC50"
    assert cyp3a4["status"] == "CANDIDATE_EXTERNAL_MODEL_EVALUATED"

    cyp2c19 = next((r for r in report if r["endpoint"] == "CYP2C19 quantitative inhibition"), None)
    assert cyp2c19 is not None
    assert cyp2c19["quantitative_model"] == "MODEL_UNAVAILABLE"
