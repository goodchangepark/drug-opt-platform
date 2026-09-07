import json
import sqlite3
from pathlib import Path


def test_completion_board_has_exact_production_version_and_terminal_policy():
    data = json.loads(Path("validation/endpoint_completion_board.json").read_text())
    assert data["production_engine_version"] == "drugopt-prediction-engine-v3@3.3.2"
    assert data["reference_library_n"] == 1000
    assert all(row["status"] in data["status_policy"] for row in data["rows"])


def test_runtime_project300_has_identity_safe_1000_compounds():
    con = sqlite3.connect("drug_opt.db")
    n = con.execute("select count(*) from compounds where project_id=300").fetchone()[0]
    keys = [x[0] for x in con.execute("select cv.inchikey from compound_versions cv join compounds c on c.id=cv.compound_row_id where c.project_id=300 and cv.inchikey is not null")]
    assert n == 1000
    assert len(keys) == len(set(keys)) == 1000


def test_reference_library_identity_manifest_is_preserved():
    data = json.loads(Path("validation/project300_identity_manifest_before_v70.json").read_text())
    assert data["project_id"] == 300
    assert data["compound_n"] == 250
