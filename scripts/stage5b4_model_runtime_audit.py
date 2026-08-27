#!/usr/bin/env python3
"""Run every installed scientific prediction checkpoint on ARM64 CPU."""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

from sqlalchemy import select

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.admet import ADMETModelRegistry
from backend.admet_predictor import (ADMET_AI_ROOT, MODEL_ROOT, MODEL_SPECS,
                                     OPENADMET_ROOT, model_files_available,
                                     predict_endpoint)
from backend.chemistry import analyze_smiles
from backend.database import SessionLocal
from backend.metabolic_soft_spot import predict_soft_spots
from backend.standardizer import standardize_molecule


OUTPUT = ROOT / "validation" / "stage5b4_stabilization_model_audit.json"
SMILES = "CC(=O)Oc1ccccc1C(=O)O"


def asset_paths(endpoint: str) -> list[Path]:
    spec = MODEL_SPECS[endpoint]
    if spec["model_family"] == "openadmet_clearance":
        return [OPENADMET_ROOT / name for name in ("model.pth", "X_train.csv", "y_train.csv", f"ad_index_{spec['index_key']}.npz")]
    if spec["model_family"] == "admet_ai_ensemble":
        training = ADMET_AI_ROOT / "training" / spec["index_key"]
        return [*(ADMET_AI_ROOT / "classification" / f"model_{index}.pt" for index in range(5)),
                training / "training.csv", training / "ad_index.npz"]
    root = MODEL_ROOT / spec["model_key"]
    return [root / "model_v2_1.pt", root / "training.csv", root / "ad_index.npz"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    with SessionLocal() as db:
        registry = {row.endpoint_name: row for row in db.scalars(select(ADMETModelRegistry))}
    rows = []
    failures = []
    for endpoint, spec in MODEL_SPECS.items():
        entry = registry.get(endpoint)
        available, reason = model_files_available(endpoint)
        assets = asset_paths(endpoint)
        try:
            result = predict_endpoint(SMILES, endpoint)
            value = result.get("predicted_value")
            finite = isinstance(value, (int, float)) and math.isfinite(value)
            passed = bool(entry and available and result.get("status") == "COMPLETE" and finite
                          and result.get("unit") == spec["unit"]
                          and entry.model_version == spec["model_version"])
            error = "" if passed else "registry, asset, inference, unit or version mismatch"
        except Exception as exc:  # audit must retain every endpoint result
            result, finite, passed, error = {}, False, False, f"{type(exc).__name__}: {exc}"
        checkpoint_assets = [path for path in assets if path.suffix in {".pt", ".pth"}]
        row = {
            "endpoint": endpoint, "model": spec["display_name"], "version": spec["model_version"],
            "registry_entry": bool(entry), "registry_status": entry.implementation_status if entry else "MISSING",
            "assets_available": available, "asset_error": reason,
            "checkpoint_assets": [{"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
                                  for path in checkpoint_assets if path.is_file()],
            "manifest_status": "AVAILABLE" if any("SHA256SUMS" in str(path) for path in assets) else "CHECKSUM_RECORDED_BY_AUDIT",
            "loader": "PASS" if result else "FAIL", "cpu_inference": "PASS" if passed else "FAIL",
            "finite_prediction": finite, "predicted_value": result.get("predicted_value"),
            "unit": result.get("unit", spec["unit"]), "confidence": result.get("confidence"),
            "applicability_domain": (result.get("applicability_domain") or {}).get("classification"),
            "endpoint_mapping": "PASS" if result.get("unit") == spec["unit"] else "FAIL",
            "error": error,
        }
        rows.append(row)
        if not passed:
            failures.append(endpoint)

    standard = standardize_molecule(SMILES)
    chemistry = analyze_smiles(SMILES)
    metabolism = predict_soft_spots(SMILES, max_spots=8)
    sygma_pass = bool(metabolism.get("engine") and metabolism.get("spots") and metabolism.get("metabolites"))
    supporting = {
        "RDKit": {"status": "PASS", "canonical_smiles": standard["canonical_smiles"],
                  "finite_molecular_weight": math.isfinite(float(chemistry["properties"]["molecular_weight"]))},
        "SyGMa": {"status": "PASS" if sygma_pass else "FAIL",
                  "availability": "LIMITED", "method": "RULE_BASED",
                  "soft_spot_count": len(metabolism.get("spots", [])),
                  "metabolite_count": len(metabolism.get("metabolites", []))},
    }
    supporting_failures = [name for name, row in supporting.items() if row["status"] != "PASS"]
    payload = {
        "audit": "Stage 5B-4 stabilization model runtime audit", "timestamp": datetime.now(timezone.utc).isoformat(),
        "architecture": platform.machine(), "processor": "CPU", "smiles": SMILES,
        "summary": {"total_models": len(rows), "passed": len(rows) - len(failures),
                    "failed": len(failures), "supporting_engines_passed": len(supporting) - len(supporting_failures),
                    "status": "PASS" if not failures and not supporting_failures else "FAIL"},
        "models": rows, "supporting_engines": supporting,
        "note": "This is an execution sanity check, not a scientific performance validation.",
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload["summary"], sort_keys=True))
    return 1 if failures or supporting_failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
