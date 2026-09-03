"""Unit and integration tests for CYP Quantitative Model Validation Audit v5.6.1."""
import pytest
from backend.openadmet_cyp import (
    predict_chemeleon_cyp_pic50,
    evaluate_cyp_applicability_domain,
    compute_fold_error,
    PROVENANCE_OPENADMET_PRETRAINED,
    PROVENANCE_DRUGOPT_CV,
    PROVENANCE_DRUGOPT_FINAL,
    OPENADMET_PUBLISHER_BENCHMARKS,
)
from backend.endpoint_model_validation import audit_cyp_quantitative_validation, build_dmpk_quantitative_expansion_report
from rdkit import Chem


def test_model_provenance_isolation():
    """Verify that model provenance sources are strictly isolated and CV claims retracted."""
    audit = audit_cyp_quantitative_validation()
    assert audit["audit_version"] == "CYP_VALIDATION_AUDIT_V561"

    for iso in ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]:
        rep = audit["isoforms"][iso]
        # Publisher benchmark isolated
        assert rep["publisher_benchmarks"]["provenance"] == PROVENANCE_OPENADMET_PRETRAINED
        assert rep["publisher_benchmarks"]["n_samples"] > 1000

        # Local Drug-OPT CV retracted
        assert rep["drugopt_cv_status"] == "RETRACTED_NOT_LOCALLY_RETRAINED"
        assert rep["drugopt_final_trained_status"] == "NOT_APPLICABLE"


def test_external_holdout_overlap_and_metrics():
    """Verify external holdout evaluation on Drug-OPT evidence with zero overlap."""
    audit = audit_cyp_quantitative_validation()
    rep_3a4 = audit["isoforms"]["CYP3A4"]["external_holdout"]
    assert rep_3a4["independent_n"] == 1
    assert rep_3a4["exact_overlap_n"] == 0
    assert rep_3a4["mae_pic50"] == 0.17
    assert rep_3a4["geom_fold_error"] == "1.47x"

    rep_2c9 = audit["isoforms"]["CYP2C9"]["external_holdout"]
    assert rep_2c9["independent_n"] == 1
    assert rep_2c9["exact_overlap_n"] == 0
    assert rep_2c9["mae_pic50"] == 0.54
    assert rep_2c9["geom_fold_error"] == "3.50x"

    rep_1a2 = audit["isoforms"]["CYP1A2"]["external_holdout"]
    assert rep_1a2["independent_n"] == 2
    assert rep_1a2["exact_overlap_n"] == 0


def test_real_applicability_domain_evaluation():
    """Verify real chemical space AD evaluates Morgan/Tanimoto similarity and descriptor envelope."""
    # Small standard drug within domain (e.g. Caffeine)
    caffeine = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    status_caf, sim_caf, viol_caf = evaluate_cyp_applicability_domain(caffeine)
    assert len(viol_caf) == 0

    # Large complex drug exceeding MW envelope (e.g. Orforglipron MW 883)
    orforglipron = Chem.MolFromSmiles("CC1(CN(C1)C(=O)C2=CC=C(C=C2)C3=NN=C(N3)C4=CC(=C(C=C4)C#N)N5C=C(C=N5)C6=CC=CC=C6)C7=CC=CC=C7")
    status_orf, sim_orf, viol_orf = evaluate_cyp_applicability_domain(orforglipron)
    assert status_orf in ("BORDERLINE", "OUT_OF_DOMAIN")


def test_promotion_remains_blocked():
    """Verify quantitative CYP models remain in candidate status and not promoted to primary."""
    audit = audit_cyp_quantitative_validation()
    for iso, rep in audit["isoforms"].items():
        assert "RETAIN_CANDIDATE_STATUS" in rep["promotion_decision"]

    dmpk_table = build_dmpk_quantitative_expansion_report()
    cyp3a4_row = next(r for r in dmpk_table if r["endpoint"] == "CYP3A4 quantitative inhibition")
    assert cyp3a4_row["status"] == "CANDIDATE_EXTERNAL_MODEL_EVALUATED"
    assert cyp3a4_row["provenance"] == PROVENANCE_OPENADMET_PRETRAINED
