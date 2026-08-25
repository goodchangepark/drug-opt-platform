"""Reproduce the Stage 3E P-gp inhibitor directionality sanity check.

These well-known inhibitors are not an independent validation set: two are exact
training structures and tariquidar has a very close training analogue.  The
script reports that overlap explicitly and does not calculate validation metrics.
"""

from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.admet_predictor import predict_endpoint  # noqa: E402


REFERENCES = {
    "Verapamil": "CC(C)C(CCCN(C)CCC1=CC(=C(C=C1)OC)OC)(C#N)C2=CC(=C(C=C2)OC)OC",
    "Tariquidar": "COC1=C(C=C2CN(CCC2=C1)CCC3=CC=C(C=C3)NC(=O)C4=CC(=C(C=C4NC(=O)C5=CC6=CC=CC=C6N=C5)OC)OC)OC",
    "Zosuquidar": "C1CN(CCN1CC(COC2=CC=CC3=C2C=CC=N3)O)C4C5=CC=CC=C5C6C(C6(F)F)C7=CC=CC=C47",
}


def validate() -> dict:
    rows = []
    for name, smiles in REFERENCES.items():
        result = predict_endpoint(smiles, "P-gp inhibitor")
        rows.append({
            "compound": name,
            "known_direction": "P-gp inhibitor",
            "probability": round(result["probability"], 6),
            "prediction": result["classification"],
            "domain": result["applicability_domain"]["classification"],
            "nearest_training_similarity": result["applicability_domain"]["nearest_training_similarity"],
            "direction_correct": result["classification"] == "INHIBITOR",
        })
    return {
        "endpoint": "Human P-gp/ABCB1 inhibitor classification",
        "independent_validation": "NOT_AVAILABLE",
        "metrics": None,
        "reason": "Reference compounds overlap or closely match the Broccatelli training chemistry; this is a directionality sanity check only.",
        "reference_structure_source": "PubChem PUG REST, retrieved 2026-08-25",
        "results": rows,
    }


if __name__ == "__main__":
    print(json.dumps(validate(), indent=2))
