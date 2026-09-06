"""
Comprehensive Regression & Scientific Integrity Test Suite:
PK Critical Parameter Upgrade & Candidate Governance.
===========================================================

Verifies all 14 mandatory Section 16 directives:
1. acid/base pKa separation
2. macro/micro pKa semantics
3. logP != logD
4. pH context
5. VDss != Vz != Vd/F
6. PPB / fu conversion (all 4 modes)
7. Human / Rat / Mouse clearance separation
8. Classifier cannot become quantitative
9. IVIVE assumption provenance (MODEL_INPUT_ASSUMPTION)
10. PK simulation requires required inputs (fails closed)
11. Uncertainty propagation for extreme binding
12. OOD downstream confidence downgrade
13. Representative evidence ranking regression (IC50 > % > categorical)
14. Mobocertinib P-gp IC50 36.1 µM regression protection
"""
import math
import pytest
from typing import Dict, Any

from backend.ionization import analyze_ionization, IonizationClass
from backend.canonical_endpoints import (
    REGISTRY,
    normalize_experimental_observation,
)
from backend.unit_normalization import (
    convert_ppb,
    convert_clearance,
    convert_volume,
)
from backend.pk_parameter_set import (
    SPECIES_PHYSIOLOGY,
    ASSUMPTION_LABEL,
    SOURCE_TYPE_EXPERIMENTAL,
    SOURCE_TYPE_VALIDATED_MODEL,
    SOURCE_TYPE_MECHANISTIC_DERIVED,
    SOURCE_TYPE_RULE_ESTIMATE,
    SOURCE_TYPE_USER_INPUT,
    SOURCE_TYPE_MODEL_UNAVAILABLE,
    canonical_fu_from_ppb,
    calculate_well_stirred_clearance,
    simulate_one_compartment_disposition,
    propagate_pk_uncertainty,
    build_pk_parameter_set,
)
from backend.representative_experimental import (
    representative_rank,
    select_representative,
)


def test_acid_base_pka_separation():
    """Directive 1: Acidic and basic pKa must be explicitly separated, not collapsed."""
    # Ciprofloxacin: Ampholyte with acidic carboxyl (~6.09) and basic piperazine (~8.62)
    cipro_smiles = "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O"
    res = analyze_ionization(cipro_smiles)

    assert res["ionization_class"] in {IonizationClass.ZWITTERION_POSSIBLE, IonizationClass.AMPHOLYTE}
    assert res["pka_acid"] is not None
    assert res["pka_base"] is not None
    assert res["macro_pka"]["pka_acidic"] is not None
    assert res["macro_pka"]["pka_basic"] is not None
    assert res["macro_pka"]["separation_status"] == "EXPLICIT_ACID_BASE_SEPARATED"

    # Verify canonical endpoint normalization routes separately
    norm_acid = normalize_experimental_observation("acidic pKa", 4.2)
    assert norm_acid["canonical_endpoint_id"] == "PKA_ACID"
    assert norm_acid["display_name"] == "Acidic pKa"

    norm_base = normalize_experimental_observation("basic pKa", 9.4)
    assert norm_base["canonical_endpoint_id"] == "PKA_BASE"
    assert norm_base["display_name"] == "Basic pKa"


def test_macro_micro_pka_semantics():
    """Directive 2: Macro pKa and micro pKa must be distinguished with atom-level centers."""
    # Propranolol: secondary amine base
    prop_smiles = "CC(C)NCC(O)COc1cccc2ccccc12"
    res = analyze_ionization(prop_smiles)

    assert "macro_pka" in res
    assert "micro_pka" in res
    assert res["micro_pka"]["macro_micro_distinct"] is True
    assert len(res["micro_pka"]["basic_centers"]) >= 1

    center = res["micro_pka"]["basic_centers"][0]
    assert "atom_index" in center
    assert "motif_name" in center
    assert "estimated_micro_pka" in center
    assert center["estimated_micro_pka"] >= 9.0


def test_logp_not_equal_logd():
    """Directive 3: logP (neutral partition) must not be substituted for logD (pH-dependent)."""
    # Ibuprofen: carboxylic acid (cLogP ~3.07, logD7.4 ~0.45 due to extensive ionization)
    ibu_smiles = "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"
    res = analyze_ionization(ibu_smiles)

    clogp = res["clogp"]
    logd74 = res["physiological_state_7_4"]["estimated_logd74"]

    assert clogp != logd74
    # For an acid, logD7.4 must be substantially lower than cLogP (ionized carboxylate partitions less into octanol)
    assert logd74 < clogp
    assert "clogp_vs_logd_distinction" in res["physiological_state_7_4"]


def test_ph_context():
    """Directive 4: pH context must define physiological gradients with defined tolerance."""
    res = analyze_ionization("CC(=O)NC1=CC=C(O)C=C1")  # Acetaminophen
    profiles = res["ph_profiles"]
    ph_values = [p["ph"] for p in profiles]

    assert 1.2 in ph_values  # Fasted stomach
    assert 4.5 in ph_values  # Duodenum
    assert 6.5 in ph_values  # Jejunum
    assert 7.4 in ph_values  # Blood / Plasma

    assert res["physiological_state_7_4"]["target_ph"] == 7.4
    assert res["physiological_state_7_4"]["ph_tolerance"] == 0.2


def test_vdss_not_vz_not_vdf():
    """Directive 5: VDss must not be mixed with Vz or oral apparent volume Vd/F."""
    vdss_def = REGISTRY.get("VDSS")
    assert vdss_def is not None
    assert vdss_def.canonical_unit == "L/kg"
    assert vdss_def.domain == "distribution"
    assert "steady-state" in vdss_def.scientific_definition.lower()

    # Vd/F is distinct oral apparent volume
    vdf_def = REGISTRY.get("HUMAN_PK_VDF_ORAL")
    assert vdf_def is not None
    assert vdf_def.canonical_endpoint_id != "VDSS"
    assert "apparent" in vdf_def.scientific_definition.lower() or "vdf" in vdf_def.canonical_endpoint_id.lower()


def test_ppb_fu_conversion_all_modes():
    """Directive 6: Deterministic conversion across % bound, fraction bound, % unbound, fu."""
    # 1. % bound -> fraction bound
    r1 = convert_ppb(98.0, "% bound", "fraction bound")
    assert math.isclose(r1.normalized_value, 0.98, rel_tol=1e-4)

    # 2. fraction bound -> % unbound
    r2 = convert_ppb(0.98, "fraction bound", "% unbound")
    assert math.isclose(r2.normalized_value, 2.0, rel_tol=1e-4)

    # 3. % unbound -> fraction unbound
    r3 = convert_ppb(2.0, "% unbound", "fu")
    assert math.isclose(r3.normalized_value, 0.02, rel_tol=1e-4)

    # 4. fraction unbound -> % bound
    r4 = convert_ppb(0.02, "fu", "% bound")
    assert math.isclose(r4.normalized_value, 98.0, rel_tol=1e-4)

    # 5. Extreme binding non-silent clipping check
    r5 = convert_ppb(99.9999, "% bound", "fu")
    assert r5.normalized_value == 0.0001
    assert "bounded at physical lower limit" in r5.conversion_formula


def test_human_rat_mouse_clearance_separation():
    """Directive 7: Multi-species clearance must preserve species separation."""
    cl_int = 30.0  # mL/min/kg
    fu = 0.10

    h_res = calculate_well_stirred_clearance(cl_int, fu, species="HUMAN")
    r_res = calculate_well_stirred_clearance(cl_int, fu, species="RAT")
    m_res = calculate_well_stirred_clearance(cl_int, fu, species="MOUSE")

    # Species physiological parameters must differ
    assert h_res["species"] == "HUMAN"
    assert r_res["species"] == "RAT"
    assert m_res["species"] == "MOUSE"

    assert h_res["qh_ml_min_kg"] == SPECIES_PHYSIOLOGY["HUMAN"]["qh_ml_min_kg"]
    assert r_res["qh_ml_min_kg"] == SPECIES_PHYSIOLOGY["RAT"]["qh_ml_min_kg"]
    assert m_res["qh_ml_min_kg"] == SPECIES_PHYSIOLOGY["MOUSE"]["qh_ml_min_kg"]

    # Scaling factors must differ
    assert h_res["scaling_factor"] != r_res["scaling_factor"]
    assert r_res["scaling_factor"] != m_res["scaling_factor"]


def test_classifier_cannot_become_quantitative():
    """Directive 8: Transporter and CYP classifiers cannot be converted to quantitative IC50."""
    cyp2c19 = REGISTRY.get("CYP2C19_INHIBITION")
    assert cyp2c19 is not None
    # Must fail closed if quantitative regression is not installed
    from backend.prediction_engine_v3_3_3_policy import V3_3_3_ENDPOINT_ROUTING
    assert V3_3_3_ENDPOINT_ROUTING["CYP2C19_INHIBITION"]["tier"] == "MODEL_UNAVAILABLE"


def test_ivive_assumption_provenance():
    """Directive 9: Assumed physiological values must be labeled MODEL_INPUT_ASSUMPTION."""
    res = calculate_well_stirred_clearance(cl_int_input=20.0, fu=0.05, species="HUMAN")
    assumptions = res["assumptions"]
    assert len(assumptions) >= 4
    for a in assumptions:
        assert ASSUMPTION_LABEL in a, f"Assumption '{a}' must contain '{ASSUMPTION_LABEL}'"


def test_pk_simulation_requires_required_inputs():
    """Directive 10: PK simulation must fail closed if explicit required inputs are missing."""
    # Missing dose
    with pytest.raises(ValueError, match="dose_mg"):
        simulate_one_compartment_disposition(dose_mg=0.0, route="IV", cl_plasma_ml_min_kg=5.0, vdss_l_kg=1.0)

    # Missing clearance
    with pytest.raises(ValueError, match="clearance"):
        simulate_one_compartment_disposition(dose_mg=100.0, route="IV", cl_plasma_ml_min_kg=0.0, vdss_l_kg=1.0)

    # Missing volume
    with pytest.raises(ValueError, match="volume"):
        simulate_one_compartment_disposition(dose_mg=100.0, route="IV", cl_plasma_ml_min_kg=5.0, vdss_l_kg=0.0)

    # Oral without bioavailability
    with pytest.raises(ValueError, match="bioavailability"):
        simulate_one_compartment_disposition(dose_mg=100.0, route="ORAL", cl_plasma_ml_min_kg=5.0, vdss_l_kg=1.0, f_oral=0.0)


def test_uncertainty_propagation_and_extreme_binding():
    """Directive 11: Low fu values (< 0.02) must propagate expanded clearance uncertainty."""
    # Normal binding (fu = 0.10)
    norm_ci = propagate_pk_uncertainty(cl_mean=5.0, vd_mean=1.0, fu_val=0.10, n_draws=2000, seed=42)
    # Extreme binding (fu = 0.005)
    ext_ci = propagate_pk_uncertainty(cl_mean=5.0, vd_mean=1.0, fu_val=0.005, n_draws=2000, seed=42)

    # 90% confidence interval width for clearance must be wider under extreme binding
    norm_width = norm_ci["cl_plasma_ml_min_kg"]["p95"] - norm_ci["cl_plasma_ml_min_kg"]["p05"]
    ext_width = ext_ci["cl_plasma_ml_min_kg"]["p95"] - ext_ci["cl_plasma_ml_min_kg"]["p05"]

    assert ext_width > norm_width, f"Expected {ext_width} > {norm_width} for extreme binding uncertainty"


def test_ood_downstream_confidence_downgrade():
    """Directive 12: OOD upstream inputs must downgrade downstream PK confidence."""
    # In-domain upstream
    pset_in = build_pk_parameter_set(
        smiles="CC(=O)Nc1ccc(O)cc1",
        dose_mg=500.0,
        route="ORAL",
        upstream_ad={"solubility": "IN_DOMAIN", "hlm": "IN_DOMAIN"}
    )
    # Out-of-domain upstream
    pset_ood = build_pk_parameter_set(
        smiles="CC(=O)Nc1ccc(O)cc1",
        dose_mg=500.0,
        route="ORAL",
        upstream_ad={"solubility": "OOD", "hlm": "IN_DOMAIN"}
    )

    assert pset_ood.confidence == "LOW_DOWNGRADED_OOD"
    assert any("out-of-domain" in a.lower() for a in pset_ood.assumptions)


def test_representative_evidence_ranking_regression():
    """Directive 13: representative_rank() must prioritize IC50/Ki/EC50 > % > categorical."""
    cand_ic50 = {
        "id": 1,
        "endpoint": "PGP_INHIBITION",
        "value": 36.1,
        "unit": "uM",
        "measurement_type": "IC50",
        "species": "Human",
        "context": "Caco-2 transport assay",
        "reference": "FDA NDA 219839",
    }
    cand_pct = {
        "id": 2,
        "endpoint": "PGP_INHIBITION",
        "value": 16.0,
        "unit": "%",
        "measurement_type": "% inhibition",
        "species": "Human",
        "context": "FDA clinical drug interaction narrative mention",
        "reference": "FDA Drug Label",
    }

    rank_ic50 = representative_rank(cand_ic50)
    rank_pct = representative_rank(cand_pct)

    # Lower tuple is better; IC50 must strictly outrank % inhibition
    assert rank_ic50 < rank_pct
    best, _ = select_representative([cand_pct, cand_ic50])
    assert best["id"] == 1
    assert best["value"] == 36.1


def test_mobocertinib_pgp_ic50_regression_protection():
    """Directive 14: Mobocertinib P-gp 36.1 uM IC50 must be selected over 16% narrative."""
    records = [
        {
            "id": 1076,
            "raw_endpoint_name": "P-gp inhibition",
            "raw_value": 16.0,
            "raw_unit": "%",
            "context": "Clinical PK drug interaction label text",
            "reference": "FDA NDA Label 2021",
        },
        {
            "id": 1366,
            "raw_endpoint_name": "P-gp IC50",
            "raw_value": 36.1,
            "raw_unit": "uM",
            "context": "In vitro Caco-2 bidirectional transport study",
            "reference": "FDA Multidiscipline Review 2021",
        }
    ]
    selected, _ = select_representative(records)
    assert selected is not None
    assert selected["id"] == 1366
    assert selected["raw_value"] == 36.1
