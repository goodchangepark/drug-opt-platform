"""v4.5 project evidence persistence, representative, and display contracts."""
from pathlib import Path

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from backend.database import Base
from backend.endpoint_comparison import _blank, _scientific_rows, build_endpoint_comparison
from backend.models import Compound, CompoundVersion, ExternalExperimentalEvidence, Project, ensure_ui_schema
from backend.representative_experimental import REPRESENTATIVE_EXPERIMENTAL_VERSION, select_representative


def _db():
    engine = create_engine("sqlite://", connect_args={"check_same_thread": False}, poolclass=StaticPool)
    Base.metadata.create_all(engine); ensure_ui_schema(engine)
    return sessionmaker(bind=engine, autoflush=False, autocommit=False)()


def test_structure_preserving_version_restores_persisted_search_evidence():
    db = _db(); project = Project(name="EGFR scope"); db.add(project); db.flush()
    compound = Compound(project_id=project.id, compound_id="E1", name="Alias revision", current_version=2)
    db.add(compound); db.flush()
    old = CompoundVersion(compound_row_id=compound.id, version_number=1, original_smiles="CCO", canonical_smiles="CCO", isomeric_smiles="CCO", inchikey="SAME")
    current = CompoundVersion(compound_row_id=compound.id, version_number=2, original_smiles="CCO", canonical_smiles="CCO", isomeric_smiles="CCO", inchikey="SAME")
    db.add_all([old, current]); db.flush()
    db.add(ExternalExperimentalEvidence(compound_version_id=old.id, provenance_key="p1", raw_endpoint_name="Plasma protein binding", raw_value="90", raw_unit="% bound", species="Human", source_database="FDA / Regulatory", source_record_id="doc-1", reference_text="Label table", identity_match_status="EXACT_STRUCTURE_MATCH", endpoint_match_status="CANONICAL", canonical_endpoint_id="HUMAN_PPB", normalized_value="90", normalized_unit="% bound", comparability_status="DIRECTLY_COMPARABLE", evidence_state="EXTERNAL_CANDIDATE", lifecycle_status="ACTIVE"))
    db.commit()
    view = build_endpoint_comparison(db, current.id)
    ppb = next(row for row in view["scientific_rows"] if row["canonical_endpoint"] == "HUMAN_PPB")
    assert ppb["representative_observation_id"] is not None
    assert ppb["experimental_display_value"] == 90.0


def test_representative_origin_priority_and_no_prediction_proximity():
    items = [
        {"id": 1, "origin": "EXTERNAL_CANDIDATE", "display": {"value": 99.9, "unit": "% bound"}, "comparability": "DIRECTLY_COMPARABLE"},
        {"id": 2, "origin": "INTERNAL_EXPERIMENTAL", "display": {"value": 70.0, "unit": "% bound"}, "comparability": "DIRECTLY_COMPARABLE"},
    ]
    chosen, reason = select_representative(items)
    assert chosen["id"] == 2
    assert REPRESENTATIVE_EXPERIMENTAL_VERSION in reason


def test_direct_caco_display_uses_same_scientist_facing_unit():
    row = _blank("CACO2_PAPP_AB")
    row["experimental_external_candidates"] = [{"id": 1, "origin": "EXTERNAL_CANDIDATE", "raw_endpoint": "Caco-2 Papp A→B", "normalized_value": -4.89279, "normalized_unit": "log10(cm/s)", "comparability": "DIRECTLY_COMPARABLE"}]
    row["prediction"] = {"available": True, "display_value": -4.9045887, "unit": "log10(cm/s)", "maturity": {"level": 1}}
    row["comparison"] = {"status": "DIRECTLY_COMPARABLE"}
    result = _scientific_rows([row])[0]
    assert result["semantic_status"] == "DIRECTLY_COMPARABLE"
    assert result["experimental_display_unit"] == "×10^-6 cm/s"
    assert result["prediction_display_unit"] == "×10^-6 cm/s"
    assert result["difference_display_value"] is not None


def test_shared_maturity_star_color_contract_and_selected_import_ui():
    root = Path(__file__).resolve().parents[1]
    css = (root / "frontend/static/app.css").read_text()
    js = (root / "frontend/static/app.js").read_text()
    assert ".maturity-star-filled { color: #F5B700; }" in css
    assert "renderPredictionMaturity(prediction.maturity.level" in js
    assert "Import Selected Qualified Evidence" in js
    assert "Review Qualified Evidence" in js
    assert "Qualified Evidence by Compound" in js
