#!/usr/bin/env python3
"""Run public-example directional sanity checks for Stage 4A.

The fixtures are deliberately not presented as independent model validation.
They verify evidence plumbing and transformation direction only; no reaction is run.
"""

import json
import sys
from pathlib import Path

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT))

from backend.admet import ensure_admet_schema
from backend.database import Base
from backend.main import create_admet_measurement, create_compound, create_optimization_run, create_project
from backend.metabolism import ensure_metabolism_schema
from backend.optimization import ensure_optimization_schema
from backend.schemas import CompoundCreate, ProjectCreate


def validate():
    dataset = json.loads((ROOT / "validation/stage4a_acceptance_examples.json").read_text())
    engine = create_engine("sqlite://", connect_args={"check_same_thread": False}, poolclass=StaticPool)
    Base.metadata.create_all(engine)
    ensure_admet_schema(engine); ensure_metabolism_schema(engine); ensure_optimization_schema(engine)
    db = sessionmaker(bind=engine, autoflush=False, autocommit=False)()
    results = []
    try:
        for index, example in enumerate(dataset["examples"], 1):
            project = create_project(ProjectCreate(name=f"Acceptance {index}"), db)
            parent = create_compound(project.id, CompoundCreate(compound_id=example["id"], smiles=example["smiles"]), db)
            version_id = parent["version"]["id"]
            if "HLM unstable" in example["fixture_evidence"]:
                create_admet_measurement(project.id, {
                    "version_id": version_id, "endpoint": "HLM intrinsic clearance", "species": "Human",
                    "matrix": "HLM", "value": 2.2, "unit": "log10(mL/min/kg)",
                    "method": "directional acceptance fixture", "source": example["reference"],
                }, db)
            run = create_optimization_run(project.id, {
                "parent_version_id": version_id, "objectives": example["objectives"],
                "constraints": {"clogp_max": 4.0},
            }, db)
            ranked_ids = [row["id"] for row in run["recommended_transformations"]]
            expected = example["expected_transformation_id"]
            results.append({
                "id": example["id"], "expected": expected, "observed_rank": ranked_ids.index(expected) + 1 if expected in ranked_ids else None,
                "pass": expected in ranked_ids, "analog_generation": run["analog_generation"],
            })
    finally:
        db.close()
    return {"scope": dataset["scope"], "passed": sum(row["pass"] for row in results), "total": len(results), "results": results}


if __name__ == "__main__":
    result = validate()
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["passed"] == result["total"] else 1)
