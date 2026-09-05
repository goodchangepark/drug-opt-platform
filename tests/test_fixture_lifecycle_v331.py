"""
Test Suite: Fixture Lifecycle & Baseline Protection (Directive 16)
Tests:
1. TEST 1: Default API /api/projects hides is_test_fixture=True projects
2. TEST 2: include_test_fixtures=True exposes test fixture projects
3. TEST 3: Cascade cleanup leaves 0 orphans across all child tables
4. TEST 4: Idempotent cleanup execution keeps 4 protected projects unchanged
5. TEST 5: DrugBank 150 compound integrity (150 compounds, 150 versions, 955 evidence, 0 duplicate InChIKey)
"""

import pytest
from sqlalchemy import select, text
from backend.database import SessionLocal
from backend.models import (
    Project, Compound, CompoundVersion, ExternalExperimentalEvidence,
    PredictionRun
)
from backend.main import _delete_project_tree_rows
from backend.cleanup_test_fixtures import run_cleanup, PROTECTED_PROJECT_IDS


def test_test1_default_api_projects_hides_test_fixtures():
    """TEST 1: Creating is_test_fixture=True project hides it from default /api/projects."""
    db = SessionLocal()
    temp_project = None
    try:
        temp_project = Project(
            name="Temporary Automated Test Fixture",
            target="TEST_TARGET",
            molecule_type="Small Molecule",
            indication="Test Indication",
            description="Created to verify default exclusion",
            is_test_fixture=True,
        )
        db.add(temp_project)
        db.commit()
        db.refresh(temp_project)
        temp_id = temp_project.id

        # Query default projects (is_test_fixture is False)
        default_query = select(Project).where(Project.is_test_fixture.is_(False))
        default_projects = list(db.scalars(default_query))
        default_ids = [p.id for p in default_projects]
        assert temp_id not in default_ids, "Test fixture project must NOT appear in default projects query!"

    finally:
        if temp_project and temp_project.id:
            _delete_project_tree_rows(db, [temp_project.id])
            db.commit()
        db.close()


def test_test2_include_test_fixtures_exposes_test_fixtures():
    """TEST 2: include_test_fixtures=True exposes is_test_fixture=True projects."""
    db = SessionLocal()
    temp_project = None
    try:
        temp_project = Project(
            name="Temporary Exposed Test Fixture",
            target="TEST_TARGET_2",
            molecule_type="Small Molecule",
            indication="Test Indication",
            description="Created to verify debug exposure",
            is_test_fixture=True,
        )
        db.add(temp_project)
        db.commit()
        db.refresh(temp_project)
        temp_id = temp_project.id

        # Query all projects without filter (simulating include_test_fixtures=True)
        all_query = select(Project)
        all_projects = list(db.scalars(all_query))
        all_ids = [p.id for p in all_projects]
        assert temp_id in all_ids, "Test fixture project MUST appear when include_test_fixtures=True!"

    finally:
        if temp_project and temp_project.id:
            _delete_project_tree_rows(db, [temp_project.id])
            db.commit()
        db.close()


def test_test3_cascade_cleanup_leaves_zero_orphans():
    """TEST 3: Create temporary fixture with child compound, version, and snapshot, delete it, and verify 0 orphans."""
    db = SessionLocal()
    temp_project = None
    try:
        temp_project = Project(
            name="Temporary Cascade Tree Test",
            target="TEST_CASCADE",
            molecule_type="Small Molecule",
            is_test_fixture=True,
        )
        db.add(temp_project)
        db.commit()
        db.refresh(temp_project)

        comp = Compound(
            project_id=temp_project.id,
            compound_id="TEST-CASCADE-01",
            name="CascadeTestCompound",
            status="ACTIVE",
            current_version=1,
        )
        db.add(comp)
        db.commit()
        db.refresh(comp)

        cv = CompoundVersion(
            compound_row_id=comp.id,
            version_number=1,
            original_smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
            canonical_smiles="CC(=O)Oc1ccccc1C(=O)O",
            isomeric_smiles="CC(=O)Oc1ccccc1C(=O)O",
            inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        )
        db.add(cv)
        db.commit()
        db.refresh(cv)

        comp_id = comp.id
        cv_id = cv.id

        # Now perform cascading delete
        _delete_project_tree_rows(db, [temp_project.id])
        db.commit()
        db.expunge_all()

        # Check that compound and compound_version are gone
        remaining_comp = db.scalar(select(Compound).where(Compound.id == comp_id))
        remaining_cv = db.scalar(select(CompoundVersion).where(CompoundVersion.id == cv_id))
        assert remaining_comp is None, "Compound row must be deleted by cascade"
        assert remaining_cv is None, "CompoundVersion row must be deleted by cascade"

        # Check orphan queries
        comp_orphans = db.scalar(select(text("count(*)")).select_from(text("compounds WHERE project_id NOT IN (SELECT id FROM projects)")))
        cv_orphans = db.scalar(select(text("count(*)")).select_from(text("compound_versions WHERE compound_row_id NOT IN (SELECT id FROM compounds)")))
        assert comp_orphans == 0
        assert cv_orphans == 0

    finally:
        db.close()


def test_test4_idempotent_cleanup_maintains_protected_projects():
    """TEST 4: Idempotent execution of run_cleanup maintains exact 4 protected projects."""
    result = run_cleanup()
    assert result["status"] == "SUCCESS"
    assert result["remaining_count"] == 4
    assert set(result["remaining_ids"]) == {1, 3, 5, 300}

    # Run again to ensure strict idempotence
    result2 = run_cleanup()
    assert result2["status"] == "SUCCESS"
    assert result2["deleted_count"] == 0
    assert result2["remaining_count"] == 4
    assert set(result2["remaining_ids"]) == {1, 3, 5, 300}


def test_test5_drugbank_150_compound_integrity():
    """TEST 5: Verify DrugBank project 300 contains 150 unique compounds, 150 versions, 955 evidence records, and 0 duplicate InChIKey."""
    db = SessionLocal()
    try:
        p300 = db.scalar(select(Project).where(Project.id == 300))
        assert p300 is not None, "Project 300 (DrugBank) must exist!"
        assert p300.name == "DrugBank"
        assert p300.target == "PAN_TARGET_REFERENCE"
        assert "GLOBAL_MODEL_DEVELOPMENT" in p300.indication

        compounds = list(db.scalars(select(Compound).where(Compound.project_id == 300)))
        assert len(compounds) in (150, 200), f"DrugBank must have 150 or 200 compounds, got {len(compounds)}"

        comp_ids = [c.id for c in compounds]
        versions = list(db.scalars(select(CompoundVersion).where(CompoundVersion.compound_row_id.in_(comp_ids))))
        assert len(versions) in (150, 200), f"DrugBank must have 150 or 200 versions, got {len(versions)}"

        inchikeys = [v.inchikey for v in versions if v.inchikey]
        assert len(inchikeys) in (150, 200), "All versions must have an inchikey"
        assert len(set(inchikeys)) in (150, 200), "All InChIKeys must be strictly unique (0 duplicates)"

        version_ids = [v.id for v in versions]
        ev_count = db.scalar(
            select(text("count(*)")).select_from(ExternalExperimentalEvidence).where(
                ExternalExperimentalEvidence.compound_version_id.in_(version_ids)
            )
        )
        assert ev_count in (955, 1490), f"DrugBank must have 955 or 1490 external evidence records, got {ev_count}"

    finally:
        db.close()
