#!/usr/bin/env python3
"""Run Stage 4B public-direction acceptance with installed deterministic models."""

import json
import sys
from pathlib import Path

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT))

from backend.activity_models import MatchedMolecularPair
from backend.admet import ensure_admet_schema
from backend.database import Base
from backend.main import add_measurement, create_admet_measurement, create_assay, create_compound, create_optimization_run, create_project
from backend.metabolism import ensure_metabolism_schema
from backend.optimization import ensure_optimization_schema
from backend.proposal import OptimizationProposalRun, ensure_proposal_schema
from backend.proposal_engine import execute_proposal_run
from backend.schemas import CompoundCreate, ProjectCreate


def validate():
    dataset = json.loads((ROOT / "validation/stage4b_acceptance_examples.json").read_text())
    engine = create_engine("sqlite://", connect_args={"check_same_thread": False}, poolclass=StaticPool)
    Base.metadata.create_all(engine)
    ensure_admet_schema(engine); ensure_metabolism_schema(engine); ensure_optimization_schema(engine); ensure_proposal_schema(engine)
    db = sessionmaker(bind=engine, autoflush=False, autocommit=False)()
    results = []
    try:
        for example in dataset["examples"]:
            project = create_project(ProjectCreate(name=example["id"]), db)
            parent = create_compound(project.id, CompoundCreate(compound_id="PARENT", name=example["parent"], smiles=example["smiles"]), db)
            assay_id = None
            if example["id"] == "PHENYL_PYRIDYL_PROJECT_MMP":
                analog = create_compound(project.id, CompoundCreate(compound_id="KNOWN-ANALOG", smiles=example["known_analog_smiles"]), db)
                assay = create_assay(project.id, {"name": "Public-direction MMP IC50", "measurement_type": "IC50", "unit": "nM"}, db)
                assay_id = assay["id"]
                for version_id in (parent["version"]["id"], analog["version"]["id"]):
                    add_measurement(assay_id, {"version_id": version_id, "value": 20, "unit": "nM", "source": "Acceptance fixture"}, db)
                pair = MatchedMolecularPair(
                    assay_id=assay_id, version_a_id=parent["version"]["id"], version_b_id=analog["version"]["id"],
                    similarity=0.65, delta_pactivity=0.0, transformation_smiles=example["smiles"] + ">>" + example["known_analog_smiles"],
                    is_cliff=False, provenance_json={"source": example["reference"], "fixture_activity": True},
                )
                db.add(pair); db.commit()
            elif example.get("experimental_fixture"):
                create_admet_measurement(project.id, {
                    "version_id": parent["version"]["id"], "endpoint": "HLM intrinsic clearance", "species": "Human",
                    "matrix": "HLM", "value": 2.2, "unit": "log10(mL/min/kg)",
                    "method": "public-direction acceptance", "source": example["reference"],
                }, db)
            optimization = create_optimization_run(project.id, {
                "parent_version_id": parent["version"]["id"], "assay_id": assay_id,
                "objectives": example["objectives"], "constraints": {"similarity_min": 0.45},
            }, db)
            proposal = OptimizationProposalRun(
                project_id=project.id, optimization_run_id=optimization["id"], parent_version_id=parent["version"]["id"],
                hard_constraints_json={}, settings_json={"max_raw_candidates": 12, "allow_double_transforms": False},
            )
            db.add(proposal); db.commit(); execute_proposal_run(proposal.id, session=db); db.refresh(proposal)
            db.expire_all(); proposal = db.get(OptimizationProposalRun, proposal.id)
            expected = example["expected_transformation"]
            matching = [candidate for candidate in proposal.candidates if any((row.transformation_id.startswith(expected) if expected == "MMP_PROJECT_OBSERVED" else row.transformation_id == expected) for row in candidate.transformations)]
            accepted = [candidate for candidate in matching if candidate.status in {"ACCEPTED", "TOP_10"}]
            results.append({
                "id": example["id"], "raw_count": proposal.raw_candidate_count,
                "accepted_count": proposal.accepted_count, "rejected_count": proposal.rejected_count,
                "expected_transformation": expected, "in_pool": bool(matching),
                "accepted": bool(accepted), "selected_top10": any(candidate.selected_top10 for candidate in accepted),
                "best_rank": min([candidate.rankings[-1].rank for candidate in accepted if candidate.rankings] or [None]),
                "matching_rejections": sorted({
                    reason.code for candidate in matching for reason in candidate.rejection_reasons
                }),
                "accepted_transformations": sorted({
                    row.transformation_id for candidate in proposal.candidates
                    if candidate.status in {"ACCEPTED", "TOP_10"}
                    for row in candidate.transformations
                }),
                "claim_limit": example["claim_limit"],
            })
    finally:
        db.close()
    passed = sum(row["in_pool"] and row["accepted"] and row["selected_top10"] for row in results)
    return {"scope": dataset["scope"], "passed": passed, "total": len(results), "results": results}


if __name__ == "__main__":
    output = validate()
    print(json.dumps(output, indent=2))
    raise SystemExit(0 if output["passed"] == output["total"] else 1)
