#!/usr/bin/env python3
"""Build the explicit per-species PK support and validation matrix."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation/v333_species_pk_matrix.json"


def cell(status: str, model: str, validation_n: int = 0, performance: dict | None = None, modes: list[str] | None = None) -> dict:
    return {
        "status": status,
        "model": model,
        "validation_n": validation_n,
        "performance": performance,
        "prediction_modes": modes or [],
    }


def main() -> None:
    ppb = json.loads((ROOT / "validation/ppb_canonical_benchmark_v72.json").read_text())["canonical_contract"]
    micro = json.loads((ROOT / "validation/microsomal_family_audit_v77.json").read_text())
    pk = json.loads((ROOT / "validation/pk_core_hybrid_ladders_v78.json").read_text())
    matrix = {
        "HUMAN": {
            "fu_ppb": cell("RESEARCH_CANDIDATE", "canonical Extra Trees log(fu)", ppb["dataset_n"], json.loads((ROOT / "validation/ppb_canonical_benchmark_v72.json").read_text())["candidates"]["extra_trees_logfu"]["pooled_oof"], ["FULL_PREDICTION"]),
            "microsomal_clint": cell("PRODUCTION_STABLE", "OpenADMET CheMeleon HLM", micro["results"]["HLM"]["independent_validation_n"], micro["results"]["HLM"]["models"]["production_chemeleon"], ["FULL_PREDICTION"]),
            "hepatic_cl": cell("RESEARCH_CANDIDATE", "well-stirred + fu/Clint residual", 30, pk["hepatic_cl"]["best_metrics"], ["ASSISTED"]),
            "total_iv_cl": cell("CURRENT_DATA_CEILING", "within-source Random Forest research", pk["total_iv_cl"]["canonical_contract"]["dataset_n"], pk["total_iv_cl"]["best_metrics"], ["FULL_PREDICTION"]),
            "vdss": cell("CURRENT_DATA_CEILING", "strict direct-structure research", pk["vdss"]["canonical_contract"]["hybrid_complete_case_n"], pk["vdss"]["best_metrics"], ["ASSISTED", "FULL_PREDICTION"]),
        },
        "RAT": {
            "fu_ppb": cell("MODEL_UNAVAILABLE", "no qualified rat fu route"),
            "microsomal_clint": cell("PRODUCTION_STABLE", "OpenADMET CheMeleon RLM", micro["results"]["RLM"]["independent_validation_n"], micro["results"]["RLM"]["models"]["production_chemeleon"], ["FULL_PREDICTION"]),
            "hepatic_cl": cell("MECHANISTIC_ONLY", "rat physiology well-stirred IVIVE", 0, modes=["ASSISTED"]),
            "total_iv_cl": cell("MODEL_UNAVAILABLE", "no qualified rat systemic-CL model"),
            "vdss": cell("MECHANISTIC_ONLY", "ionization/fu distribution estimate", 0, modes=["ASSISTED"]),
        },
        "MOUSE": {
            "fu_ppb": cell("MODEL_UNAVAILABLE", "no qualified mouse fu route"),
            "microsomal_clint": cell("PRODUCTION_STABLE_NO_INDEPENDENT_VALIDATION", "OpenADMET CheMeleon MLM", micro["MLM"]["independent_validation_n"], modes=["FULL_PREDICTION"]),
            "hepatic_cl": cell("MECHANISTIC_ONLY", "mouse physiology well-stirred IVIVE", 0, modes=["ASSISTED"]),
            "total_iv_cl": cell("MODEL_UNAVAILABLE", "no qualified mouse systemic-CL model"),
            "vdss": cell("MECHANISTIC_ONLY", "ionization/fu distribution estimate", 0, modes=["ASSISTED"]),
        },
        "DOG": {
            "fu_ppb": cell("MODEL_UNAVAILABLE", "no qualified dog fu route"),
            "microsomal_clint": cell("MODEL_UNAVAILABLE", "DLM checkpoint unavailable"),
            "hepatic_cl": cell("CONTEXT_REQUIRED", "dog physiology well-stirred IVIVE", 0, modes=["ASSISTED"]),
            "total_iv_cl": cell("MODEL_UNAVAILABLE", "no qualified dog systemic-CL model"),
            "vdss": cell("CONTEXT_REQUIRED", "requires species-matched fu", 0, modes=["ASSISTED"]),
        },
        "MONKEY": {
            "fu_ppb": cell("MODEL_UNAVAILABLE", "no qualified monkey fu route"),
            "microsomal_clint": cell("MODEL_UNAVAILABLE", "CyLM checkpoint unavailable"),
            "hepatic_cl": cell("CONTEXT_REQUIRED", "monkey physiology well-stirred IVIVE", 0, modes=["ASSISTED"]),
            "total_iv_cl": cell("MODEL_UNAVAILABLE", "no qualified monkey systemic-CL model"),
            "vdss": cell("CONTEXT_REQUIRED", "requires species-matched fu", 0, modes=["ASSISTED"]),
        },
    }
    for species in matrix.values():
        species.update({
            "half_life": cell("CONTEXT_REQUIRED", "requires species-matched predicted CL and Vdss"),
            "F": cell("CONTEXT_REQUIRED", "requires route/formulation/regimen and absorption/first-pass inputs"),
            "ka": cell("CONTEXT_REQUIRED", "requires route/formulation or qualified Tmax"),
            "AUC": cell("CONTEXT_REQUIRED", "requires dose/route/F/CL and regimen"),
            "Cmax": cell("CONTEXT_REQUIRED", "requires dose/route/formulation/F/CL/V/ka and regimen"),
        })
    output = {
        "artifact": "V333_SPECIES_PK_VALIDATION_MATRIX",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "engine_id": "drugopt-prediction-engine-v3@3.3.3",
        "species": matrix,
        "mode_contract": {
            "ASSISTED": "one or more same-compound experimental upstream inputs; never called prediction-only validation",
            "HYBRID": "mix of frozen predicted and explicitly supplied upstream inputs",
            "FULL_PREDICTION": "all upstream values are model predictions or deterministic structure-derived inputs",
        },
        "invariants": [
            "Metrics are reported per species, never pooled to hide an unsupported species.",
            "Human fu is never substituted for another species.",
            "Dog and monkey physiology is available, but physiology alone is not an endpoint model.",
            "Context-dependent exposure endpoints remain unavailable without dose/route/formulation/regimen.",
        ],
    }
    OUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
