from pathlib import Path

from backend.experimental_endpoint_aliases import classify_experimental_endpoint
from backend.schemas import CompoundCreate
from backend.database import Base
from backend.models import Compound, Project
from backend.main import create_compound, compound_out
from sqlalchemy import create_engine
from sqlalchemy.orm import Session


def test_empty_cas_schema_normalizes_to_optional_value():
    assert CompoundCreate(compound_id="X", smiles="CCO").cas_number is None
    assert CompoundCreate(compound_id="X", smiles="CCO", cas_number="").cas_number == ""


def test_casless_compound_creation_and_reload_are_supported():
    engine = create_engine("sqlite:///:memory:")
    Base.metadata.create_all(engine)
    with Session(engine) as db:
        project = Project(name="CAS optional", molecule_type="Small Molecule")
        db.add(project); db.commit(); db.refresh(project)
        saved = create_compound(project.id, CompoundCreate(compound_id="NO-CAS", smiles="CCO", cas_number="   "), db)
        assert saved["cas_number"] is None
        assert compound_out(db.get(Compound, saved["row_id"]))["cas_number"] is None


def test_endpoint_aliases_require_assay_context():
    accepted = classify_experimental_endpoint(label="Apparent permeability", description="Caco-2 Papp A to B")
    rejected = classify_experimental_endpoint(label="Apparent permeability", description="MDCK efflux B to A")
    assert accepted["endpoint"] == "CACO2_PAPP_AB" and accepted["qualified"]
    assert not rejected["qualified"]


def test_shared_maturity_component_uses_five_inline_svg_stars():
    js = (Path(__file__).parents[1] / "frontend/static/app.js").read_text()
    css = (Path(__file__).parents[1] / "frontend/static/app.css").read_text()
    assert "renderPredictionMaturity" in js
    assert "[...Array(5)]" in js and "e('svg'" in js
    assert "normalized=Math.max(1,Math.min(5,Number(level)||1))" in js
    assert "#F5B700" in css and "min-width: 88px" in css


def test_external_observations_route_to_endpoint_sections_not_compound_information():
    js = (Path(__file__).parents[1] / "frontend/static/app.js").read_text()
    assert "routedEvidenceSection('ACTIVITY'" in js
    assert "routedEvidenceSection('ADMET'" in js
    assert "routedEvidenceSection('METABOLISM'" in js
    assert "routedEvidenceSection('PK'" in js
    assert "Imported external observations are displayed in their canonical" in js
    assert "Imported external observations (" not in js


def test_harvester_categories_and_regulatory_adapter_are_exposed(monkeypatch):
    from backend import experimental_harvester as h
    monkeypatch.setattr(h, "configured_adapters", lambda: [h.RegulatoryAdapter()])
    monkeypatch.setattr(h.RegulatoryAdapter, "harvest", lambda _self, _identity: [
        h._record("FDA / Regulatory", "SPL:1", "Cmax", "12 ng/mL", evidence_category="PK", reference_status="REFERENCE_RESOLVED_REGULATORY")
    ])
    result = h.harvest_public_evidence(h.PublicIdentity(name="Reference drug"))
    assert result["summary"]["categories"]["PK"] == 1
    assert result["source_counts"]["FDA / Regulatory"] == {"found": 1, "qualified": 1}
