#!/usr/bin/env python3
"""Reproduce the Stage 3D known-drug directionality sanity check."""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.metabolic_soft_spot import (ENGINE_NAME, ENGINE_VERSION,
                                         PUBLISHER_VALIDATION,
                                         predict_soft_spots)


REFERENCES = [
    {
        "compound": "Phenacetin", "smiles": "CCOc1ccc(NC(C)=O)cc1",
        "transformation": "O-dealkylation", "atom_index": 2,
        "known_metabolite": "CC(=O)Nc1ccc(O)cc1", "known_metabolite_name": "Paracetamol",
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6225321/",
    },
    {
        "compound": "Diazepam", "smiles": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
        "transformation": "N-dealkylation", "atom_index": 1,
        "known_metabolite": "O=C1CN=C(c2ccccc2)c2cc(Cl)ccc2N1", "known_metabolite_name": "Nordazepam",
        "source": "https://pubmed.ncbi.nlm.nih.gov/2903030/",
    },
    {
        "compound": "Acetanilide", "smiles": "CC(=O)Nc1ccccc1",
        "transformation": "Aromatic hydroxylation", "atom_index": 7,
        "known_metabolite": "CC(=O)Nc1ccc(O)cc1", "known_metabolite_name": "Paracetamol",
        "source": "https://pubmed.ncbi.nlm.nih.gov/1226/",
    },
    {
        "compound": "Celecoxib", "smiles": "NS(=O)(=O)c1ccc(-n2nc(C(F)(F)F)cc2-c2ccc(C)cc2)cc1",
        "transformation": "Benzylic oxidation", "atom_index": 21,
        "known_metabolite": "NS(=O)(=O)c1ccc(-n2nc(C(F)(F)F)cc2-c2ccc(CO)cc2)cc1",
        "known_metabolite_name": "Hydroxycelecoxib",
        "source": "https://pubmed.ncbi.nlm.nih.gov/10681375/",
    },
    {
        "compound": "Procaine", "smiles": "CCN(CC)CCOC(=O)c1ccc(N)cc1",
        "transformation": "Ester hydrolysis", "atom_index": 8,
        "known_metabolite": "Nc1ccc(C(=O)O)cc1", "known_metabolite_name": "p-Aminobenzoic acid",
        "source": "https://doi.org/10.1016/0009-8981(63)90199-2",
    },
]


def validate():
    rows = []
    for reference in REFERENCES:
        prediction = predict_soft_spots(reference["smiles"])
        matching_spots = [
            spot for spot in prediction["spots"]
            if spot["transformation"] == reference["transformation"]
            and spot["atom_index"] == reference["atom_index"]
        ]
        rank = min((spot["rank"] for spot in matching_spots), default=None)
        metabolite_match = any(
            item["canonical_smiles"] == reference["known_metabolite"]
            and item["transformation"] == reference["transformation"]
            for item in prediction["metabolites"]
        )
        rows.append({
            **reference, "predicted_rank": rank,
            "atom_site_recalled": rank is not None,
            "known_metabolite_generated": metabolite_match,
            "top_3_directionally_consistent": rank is not None and rank <= 3 and metabolite_match,
        })
    count = len(rows)
    metrics = {
        "n": count,
        "top_1_accuracy": sum(row["predicted_rank"] == 1 and row["known_metabolite_generated"] for row in rows) / count,
        "top_2_accuracy": sum(row["predicted_rank"] is not None and row["predicted_rank"] <= 2 and row["known_metabolite_generated"] for row in rows) / count,
        "top_3_accuracy": sum(row["predicted_rank"] is not None and row["predicted_rank"] <= 3 and row["known_metabolite_generated"] for row in rows) / count,
        "atom_level_recall": sum(row["atom_site_recalled"] for row in rows) / count,
    }
    return {
        "validation_type": "Known-drug directional sanity check; not an independent benchmark",
        "engine": ENGINE_NAME, "engine_version": ENGINE_VERSION,
        "publisher_reported_validation": PUBLISHER_VALIDATION,
        "training_overlap_audit": {
            "status": "NOT_ASSESSABLE",
            "reason": "SyGMa's historical MDL Metabolite source database is discontinued and its training compound list is not distributed.",
        },
        "metric_definition": "A hit requires the documented transformation at the curated zero-based atom and generation of the documented canonical metabolite.",
        "metrics": metrics, "references": rows,
    }


if __name__ == "__main__":
    result = validate()
    destination = ROOT / "models/sygma/validation/known_drug_sanity.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
