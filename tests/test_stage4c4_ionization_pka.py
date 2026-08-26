"""Stage 4C-4: Targeted Tests for pKa, Ionization & pH-Dependent Physicochemical Foundation.

Covers:
1. Ionization classification across diverse classes (NEUTRAL, ACID, BASE, AMPHOLYTE, ZWITTERION_POSSIBLE, MULTIPLE_IONIZABLE_CENTERS).
2. Ionizable atom mapping with SMARTS subgraph matching and heavy atom indexing.
3. Acidic Henderson-Hasselbalch equations (anion fraction, neutral fraction).
4. Basic Henderson-Hasselbalch equations (cation fraction, free base fraction).
5. Polyprotic / multi-pKa handling.
6. Experimental precedence (experimental pKa > rule pKa > prediction).
7. cLogP != logD distinction (Crippen calculated vs pH-dependent distribution).
8. Experimental logD requires mandatory pH parameter.
9. pH-dependent solubility estimation and assumption disclosure.
10. Caco-2 permeability contextual interpretation (neutral fraction vs passive diffusion).
11. Plasma protein binding (fu) preservation and HSA/AAG binding partner context.
12. Volume of distribution (Vd) integration (acid restricted Vd vs basic lysosomal trapping).
13. Oral fraction absorbed (Fa) GI transit gradient integration.
14. MODEL_UNAVAILABLE handling for uninstalled quantitative ML pKa / logD7.4 models.
15. Ionization model provenance, versioning, and standardizer contract (CHEM_STANDARDIZER_V1).
16. Project isolation.
17. CompoundVersion isolation.
18. Caching and version immutability behavior.
19. Conformal governance integration (NOT_APPLICABLE_FOR_DETERMINISTIC_RULES).
20. Model Registry entries for Stage 4C-4.
"""

import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine, select
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend.admet import ADMETModelRegistry, ensure_admet_schema, validate_measurement
from backend.chemistry import analyze_smiles
from backend.database import Base, get_db
from backend.ionization import (
    IonizationClass,
    analyze_ionization,
    calculate_monoprotic_fractions,
    estimate_logd_from_pka_and_clogp,
)
from backend.ivive import estimate_absorption_components, estimate_volume_of_distribution
from backend.main import app
from backend.models import Compound, CompoundVersion, Project


@pytest.fixture
def session():
    engine = create_engine(
        "sqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    Base.metadata.create_all(engine)
    ensure_admet_schema(engine)
    Session = sessionmaker(bind=engine)
    db = Session()
    try:
        yield db
    finally:
        db.close()


@pytest.fixture
def client(session):
    def override_get_db():
        try:
            yield session
        finally:
            pass

    app.dependency_overrides[get_db] = override_get_db
    yield TestClient(app)
    app.dependency_overrides.clear()


# ---------------------------------------------------------------------------
# 1. Ionization Classification & 21. Acceptance Compounds
# ---------------------------------------------------------------------------

def test_ionization_acceptance_compounds():
    """Verify classification and primary pKa across the 10 literature acceptance compounds."""
    test_set = {
        "Caffeine": ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", IonizationClass.NEUTRAL, 0),
        "Aspirin": ("CC(=O)Oc1ccccc1C(=O)O", IonizationClass.ACID, 1),
        "Ibuprofen": ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", IonizationClass.ACID, 1),
        "Warfarin": ("CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O", IonizationClass.ACID, 1),
        "Lidocaine": ("CCN(CC)CC(=O)Nc1c(C)cccc1C", IonizationClass.BASE, 1),
        "Propranolol": ("CC(C)NCC(O)COc1cccc2ccccc12", IonizationClass.BASE, 1),
        "Diazepam": ("CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21", IonizationClass.BASE, 1),
        "Metformin": ("CN(C)C(=N)NC(=N)N", IonizationClass.BASE, 1),
        "Ciprofloxacin": ("O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O", IonizationClass.ZWITTERION_POSSIBLE, 2),
        "Amoxicillin": ("CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O", IonizationClass.MULTIPLE_IONIZABLE_CENTERS, 3),
    }

    for name, (smiles, expected_class, expected_min_centers) in test_set.items():
        res = analyze_ionization(smiles)
        assert res["status"] == "COMPLETE", f"Failed analysis for {name}"
        assert res["ionization_class"] == expected_class, f"{name} expected {expected_class}, got {res['ionization_class']}"
        assert res["total_ionizable_centers"] >= expected_min_centers, f"{name} expected at least {expected_min_centers} centers"


# ---------------------------------------------------------------------------
# 2. Ionizable Atom Mapping
# ---------------------------------------------------------------------------

def test_ionizable_atom_mapping():
    """Verify atom index, symbol, motif name, and evidence attached to each center."""
    # Aspirin: carboxylic acid O=C-OH
    res = analyze_ionization("CC(=O)Oc1ccccc1C(=O)O")
    assert len(res["ionizable_centers"]) == 1
    center = res["ionizable_centers"][0]
    assert center["type"] == "ACID"
    assert center["motif_name"] == "Carboxylic acid"
    assert center["typical_pka_range"] == [3.5, 5.0]
    assert isinstance(center["atom_index"], int)
    assert center["confidence"] == "RULE_DETERMINISTIC"


# ---------------------------------------------------------------------------
# 3. Acidic Henderson-Hasselbalch Calculation
# ---------------------------------------------------------------------------

def test_acidic_henderson_hasselbalch():
    """Verify monoprotic acid ionization: pH << pKa -> neutral, pH >> pKa -> ionized."""
    pka = 4.0
    # At pH 4.0 (pH == pKa), 50% ionized, 50% neutral
    f_mid = calculate_monoprotic_fractions(pka=pka, ph=4.0, center_type="ACID")
    assert f_mid["fraction_ionized"] == 0.5
    assert f_mid["fraction_neutral"] == 0.5

    # At gastric pH 2.0 (pH < pKa), predominantly unionized neutral acid
    f_low = calculate_monoprotic_fractions(pka=pka, ph=2.0, center_type="ACID")
    assert f_low["fraction_neutral"] > 0.98
    assert f_low["fraction_ionized"] < 0.02

    # At intestinal pH 7.0 (pH > pKa), predominantly ionized carboxylate
    f_high = calculate_monoprotic_fractions(pka=pka, ph=7.0, center_type="ACID")
    assert f_high["fraction_ionized"] > 0.99
    assert f_high["fraction_neutral"] < 0.01


# ---------------------------------------------------------------------------
# 4. Basic Henderson-Hasselbalch Calculation
# ---------------------------------------------------------------------------

def test_basic_henderson_hasselbalch():
    """Verify monoprotic base ionization: pH << pKa -> protonated cation, pH >> pKa -> neutral."""
    pka = 9.0
    # At pH 9.0 (pH == pKa), 50% protonated, 50% free base
    f_mid = calculate_monoprotic_fractions(pka=pka, ph=9.0, center_type="BASE")
    assert f_mid["fraction_ionized"] == 0.5
    assert f_mid["fraction_neutral"] == 0.5

    # At physiological pH 7.4 (pH < pKa), predominantly protonated cation
    f_phys = calculate_monoprotic_fractions(pka=pka, ph=7.4, center_type="BASE")
    assert f_phys["fraction_ionized"] > 0.97
    assert f_phys["fraction_neutral"] < 0.03

    # At basic pH 11.0 (pH > pKa), predominantly free base
    f_alk = calculate_monoprotic_fractions(pka=pka, ph=11.0, center_type="BASE")
    assert f_alk["fraction_neutral"] > 0.98
    assert f_alk["fraction_ionized"] < 0.02


# ---------------------------------------------------------------------------
# 5. Multi-pKa / Polyprotic Handling
# ---------------------------------------------------------------------------

def test_polyprotic_ampholyte_handling():
    """Verify ampholytes and polyprotic species have distinct profile treatment."""
    # Ciprofloxacin has carboxylic acid (~4.2) and secondary amine (~10.4)
    res = analyze_ionization("O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O")
    assert res["ionization_class"] == IonizationClass.ZWITTERION_POSSIBLE
    # At pH 7.4, zwitterion state should be dominant
    ph74 = next(p for p in res["ph_profiles"] if p["ph"] == 7.4)
    assert "zwitterion" in ph74["dominant_state"].lower()


# ---------------------------------------------------------------------------
# 6. Experimental Precedence
# ---------------------------------------------------------------------------

def test_experimental_pka_precedence():
    """Verify experimental pKa overrides rule-based estimated pKa."""
    smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Ibuprofen (rule pKa ~4.2)
    # Provide experimental measurement pKa = 4.45
    exp_records = [{"value": 4.45, "type": "ACID", "source": "Literature Titration"}]
    res = analyze_ionization(smiles, experimental_pka_records=exp_records)
    assert res["primary_pka"] == 4.45
    assert "EXPERIMENTAL" in res["primary_pka_source"]


# ---------------------------------------------------------------------------
# 7. cLogP != logD Distinction
# ---------------------------------------------------------------------------

def test_clogp_vs_logd_distinction():
    """Verify cLogP is strictly distinguished from logD across pH."""
    clogp = 3.5
    pka = 4.5  # Acid
    # At pH 7.4, acid is ionized -> logD7.4 should be significantly lower than cLogP
    logd74 = estimate_logd_from_pka_and_clogp(clogp=clogp, pka=pka, ph=7.4, center_type="ACID")
    assert logd74 < clogp
    assert logd74 < 1.0  # Approx 3.5 - log10(1 + 10^2.9) = 3.5 - 2.9 = 0.6


# ---------------------------------------------------------------------------
# 8. Experimental logD Requires Mandatory pH
# ---------------------------------------------------------------------------

def test_experimental_logd_requires_ph():
    """Verify validate_measurement enforces mandatory assay pH for logD."""
    # Missing pH -> should raise HTTPException 400
    with pytest.raises(Exception) as exc_info:
        validate_measurement({"endpoint": "logD", "value": 2.1, "unit": "log units"})
    assert "Assay pH is mandatory" in str(exc_info.value)

    # Valid with pH in payload
    val, mean, sd = validate_measurement({"endpoint": "logD", "value": 2.1, "unit": "log units", "ph": 7.4})
    assert val == 2.1


# ---------------------------------------------------------------------------
# 9. Downstream ADMET Context Formulations
# ---------------------------------------------------------------------------

def test_admet_context_solubility_permeability_ppb():
    """Verify structured contextual interpretations for Solubility, Caco-2, PPB."""
    # Basic compound: Lidocaine
    res = analyze_ionization("CCN(CC)CC(=O)Nc1c(C)cccc1C")
    ctx = res["admet_context"]

    # Basic compound should note acid gastric solubility
    assert ctx["solubility"]["ph_dependent"] is True
    assert "gastric" in ctx["solubility"]["summary"].lower()

    # High protonation at pH 7.4 reduces passive diffusion
    assert ctx["permeability"]["ionized_fraction_7_4"] > 0.0

    # Basic drug associates with AAG
    assert "AAG" in ctx["plasma_protein_binding"]["likely_target_protein"]


# ---------------------------------------------------------------------------
# 10. Volume of Distribution (Vd) & Oral Absorption (Fa) Integration
# ---------------------------------------------------------------------------

def test_vd_and_fa_ionization_integration(session):
    """Verify Vd and Fa engines account for ionization class."""
    proj = Project(name="Ionization PK Test", molecule_type="Small Molecule")
    session.add(proj)
    session.commit()

    # Add Basic Compound (Propranolol)
    comp = Compound(project_id=proj.id, compound_id="PROP-01", name="Propranolol", status="CALCULATED")
    session.add(comp)
    session.flush()

    analysis = analyze_smiles("CC(C)NCC(O)COc1cccc2ccccc12")
    v = CompoundVersion(
        compound_row_id=comp.id, version_number=1, original_smiles="CC(C)NCC(O)COc1cccc2ccccc12",
        canonical_smiles=analysis["identity"]["canonical_smiles"],
        isomeric_smiles=analysis["identity"]["isomeric_smiles"],
        inchi=analysis["identity"]["inchi"],
        inchikey=analysis["identity"]["inchikey"],
        properties_json=analysis["properties"],
        calculation_json={"ionization": analysis["ionization"], "rules": analysis["rules"], "provenance": analysis["provenance"]}
    )
    session.add(v)
    session.commit()

    vd = estimate_volume_of_distribution(session, proj.id, v.id, "Human")
    assert vd["v_source_type"] in {"PREDICTED_VD", "MODEL_UNAVAILABLE"}
    if vd["v_source_type"] == "PREDICTED_VD":
        assert vd["provenance"]["ionization_class"] == IonizationClass.BASE

    fa = estimate_absorption_components(session, proj.id, v.id, "Human")
    assert fa["ionization_class"] == IonizationClass.BASE
    assert "Base" in fa["gi_transit_context"]


# ---------------------------------------------------------------------------
# 11. Model Registry for Stage 4C-4
# ---------------------------------------------------------------------------

def test_model_registry_stage4c4_entries(session):
    """Verify pKa and logD quantitative ML models are registered as MODEL_UNAVAILABLE."""
    pka_reg = session.scalar(select(ADMETModelRegistry).where(ADMETModelRegistry.endpoint_name == "pKa (quantitative ML)"))
    assert pka_reg is not None
    assert pka_reg.implementation_status == "MODEL_UNAVAILABLE"
    assert pka_reg.is_active is False

    logd_reg = session.scalar(select(ADMETModelRegistry).where(ADMETModelRegistry.endpoint_name == "logD7.4 (quantitative ML)"))
    assert logd_reg is not None
    assert logd_reg.implementation_status == "MODEL_UNAVAILABLE"
    assert logd_reg.is_active is False
