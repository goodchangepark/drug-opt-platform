#!/usr/bin/env python3
"""Non-destructive Stage 4E-2R ExpansionRx intake summary.

The raw CSV stays outside Git.  This script never changes it and deliberately
keeps censored Papp observations separate from numeric values.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from backend.standardizer import STANDARDIZER_VERSION, standardize_molecule


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("raw_csv", type=Path)
    parser.add_argument("summary_json", type=Path)
    args = parser.parse_args()
    raw_bytes = args.raw_csv.read_bytes()
    rows = list(csv.DictReader(raw_bytes.decode().splitlines()))
    if "SMILES" not in rows[0] or "Caco-2 Permeability Papp A>B" not in rows[0]:
        raise ValueError("ExpansionRx CSV does not have the expected Caco-2 A→B columns")
    valid = numeric = censored = invalid = duplicate_rows = 0
    canonical_with_value: set[str] = set()
    seen: set[str] = set()
    for row in rows:
        result = standardize_molecule(row["SMILES"])
        canonical = result.get("canonical_smiles")
        if not canonical:
            invalid += 1
            continue
        valid += 1
        value = row["Caco-2 Permeability Papp A>B"].strip()
        if not value:
            continue
        if canonical in seen:
            duplicate_rows += 1
        seen.add(canonical)
        canonical_with_value.add(canonical)
        try:
            float(value)
            numeric += 1
        except ValueError:
            censored += 1
    summary = {
        "raw_n": len(rows),
        "valid_structure_n": valid,
        "invalid_structure_n": invalid,
        "caco2_papp_ab_reported_n": numeric + censored,
        "caco2_papp_ab_numeric_n": numeric,
        "caco2_papp_ab_censored_n": censored,
        "unique_canonical_with_caco2": len(canonical_with_value),
        "within_dataset_canonical_duplicate_rows": duplicate_rows,
        "standardizer_version": STANDARDIZER_VERSION,
        "raw_sha256": hashlib.sha256(raw_bytes).hexdigest(),
    }
    args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
