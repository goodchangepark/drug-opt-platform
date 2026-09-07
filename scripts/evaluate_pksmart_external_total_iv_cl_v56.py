"""Locked evaluation on PKSmart's separately published External_test_315 file.

The model is trained only on the already-frozen PKSmart Human_PK development
set. Exact Morgan duplicates are excluded from the external holdout. The
result is a holdout check, not a claim of independent-source validation.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from rdkit import DataStructs

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import scripts.evaluate_human_total_iv_cl_v56 as base


def main() -> None:
    external = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/pksmart-external-test-315.csv")
    base.SOURCE = Path("/tmp/pksmart-human-pk.csv")
    training, source_qc = base.rows_from_source()
    from backend.clinical_pk_cohort import CLINICAL_PK_COHORT
    n30 = {base.canonical(x.smiles) for x in CLINICAL_PK_COHORT}
    n30_mols = [base.Chem.MolFromSmiles(x.smiles) for x in CLINICAL_PK_COHORT if x.route == "IV" and x.cl_systemic_ml_min_kg is not None]
    n30_fps = [base.FINGERPRINT.GetFingerprint(m) for m in n30_mols if m]
    train_fps = [base.FINGERPRINT.GetFingerprint(r["mol"]) for r in training]
    development = [r for r in training if r["canonical_smiles"] not in n30]
    base.SOURCE = external
    ext, external_qc = base.rows_from_source()
    # External_test_315 is a separate file, but remains a PKSmart-derived
    # holdout. Remove any exact structure overlap with training or N30.
    training_keys = {r["canonical_smiles"] for r in training}
    overlap_training = []
    evaluation = []
    for r in ext:
        fp = base.FINGERPRINT.GetFingerprint(r["mol"])
        same_train = r["canonical_smiles"] in training_keys or any(DataStructs.TanimotoSimilarity(fp, x) >= 0.999999 for x in train_fps)
        same_n30 = r["canonical_smiles"] in n30 or any(DataStructs.TanimotoSimilarity(fp, x) >= 0.999999 for x in n30_fps)
        if same_train or same_n30:
            overlap_training.append(r["canonical_smiles"])
        else:
            evaluation.append(r)
    artifact = json.loads((ROOT / "validation/human_total_iv_cl_v5_6.json").read_text())
    selected = artifact["cv"]["selected"]
    mode = artifact["cv"]["results"][selected]["mode"]
    model = base.models()[selected]
    x_train = base.build_matrix(development, mode)
    y_train = np.log10(np.asarray([r["y"] for r in development]))
    model.fit(x_train, y_train)
    pred = 10 ** model.predict(base.build_matrix(evaluation, mode))
    observed = np.asarray([r["y"] for r in evaluation])
    metrics = base.metric(observed, pred)
    train_fps = [base.FINGERPRINT.GetFingerprint(r["mol"]) for r in development]
    ad_rows = []
    for r, p in zip(evaluation, pred):
        sim = max(DataStructs.BulkTanimotoSimilarity(base.FINGERPRINT.GetFingerprint(r["mol"]), train_fps), default=0.0)
        ad = "IN_DOMAIN" if sim >= 0.45 else ("BORDERLINE" if sim >= 0.30 else "OOD")
        ad_rows.append({"ad": ad, "max_train_tanimoto": float(sim), "observed": r["y"], "predicted": float(p)})
    by_ad = {}
    for label in ("IN_DOMAIN", "BORDERLINE", "OOD"):
        rows = [r for r in ad_rows if r["ad"] == label]
        if rows:
            by_ad[label] = base.metric(np.asarray([r["observed"] for r in rows]), np.asarray([r["predicted"] for r in rows]))
    out = {
        "artifact": "pksmart_external_total_iv_cl_v5_6",
        "source": {"name": "PKSmart External_test_315.csv", "path": str(external), "sha256": hashlib.sha256(external.read_bytes()).hexdigest(), "qc": external_qc, "qualification": "Human IV total CL fields; separate PKSmart holdout file, not an independent-source cohort."},
        "training": {"source": "PKSmart Human_PK_data.csv", "development_n": len(development), "model_selection_n": len(development), "model": selected, "features": mode},
        "overlap_guard": {"exact_or_morgan_overlap_excluded": len(overlap_training), "n30_overlap_excluded": sum(1 for r in ext if r["canonical_smiles"] in n30)},
        "evaluation_n": len(evaluation),
        "metrics": metrics,
        "ad": by_ad,
        "production_qualification": "NOT_SUFFICIENT_FOR_INDEPENDENT_EXTERNAL_VALIDATION",
        "reason": "The holdout is published by the same PKSmart resource family; obtain a genuinely independent source before production qualification.",
    }
    (ROOT / "validation/pksmart_external_total_iv_cl_v5_6.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"evaluation_n": len(evaluation), "excluded": len(overlap_training), "metrics": metrics, "ad": by_ad}, indent=2))


if __name__ == "__main__":
    main()
