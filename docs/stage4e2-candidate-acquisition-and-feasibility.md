# Stage 4E-2 — Candidate Acquisition & Runtime Feasibility Qualification

Stage 4E-2 is a fail-closed technical and legal gate, not a benchmark or a
model promotion. It leaves the Stage 4D-7 production registry, endpoint
strategies, thresholds, calibration, runtime orchestration, and visible
predictions unchanged.

## Result

None of the Stage 4E-1 queue is eligible for Stage 4E-3 yet. No candidate is
registered in Save & Predict, including as a new shadow model. The resulting
recommendation is `MODEL_ACQUISITION_BLOCKED` until the listed gates are
resolved.

| Candidate | Result | Gate |
| --- | --- | --- |
| CardioGenAI hERG | `LEGAL_REVIEW_REQUIRED` | Separate checkpoint/data terms and exact versioned hERG-head identity are absent. |
| MetaboGNN HLM/MLM | `LEGAL_REVIEW_REQUIRED` | Code, checkpoint, and data rights are unclear; output units cannot yet be mapped to the species-specific contracts. |
| pkasolver-lite | `ARM64_WORKAROUND_REQUIRED` | MIT source-contained lite weights and microstate semantics passed source review, but an isolated ARM64 PyG runtime remains to be reproduced. |
| pKaLearn | `LEGAL_REVIEW_REQUIRED` | MIT code is clear, but weight/data lineage requires review and its supplied CUDA/Windows-oriented environment is not Xavier evidence. |

The temporary source snapshots used for the pKa source/weight identity audit
were removed from `/tmp` after recording repository commits and selected
weight SHA256 values. No checkpoint, raw dataset, isolated environment, or
candidate binary is retained in this repository.

## Dataset gate

Biogen prospective metadata reports N=3,521 but does not itself grant raw
data access or reuse rights. ExpansionRx has endpoint listings but no resolved
raw-release terms. The logD7.4 repository does not establish a reusable data
license. Accordingly no raw rows were acquired or standardized, and no exact
canonical-SMILES overlap analysis was claimed.

The Stage 4E-3 plan is deliberately blocked: after legal source access is
resolved, raw data must be immutable, structures standardized by the current
standardizer, and exact source-lineage overlap exclusions performed before
any scientific metric, bootstrap, calibration, complementarity, or
applicability-domain analysis.

## Technical safety

The production virtual environment was inspected only: it has CPU PyTorch and
RDKit but no `torch_geometric` or DGL. No package was installed or upgraded.
Thus an ARM64 result of `NOT_TESTED` is evidence of an intentionally preserved
production boundary, not a claim that a candidate works on Xavier. Any later
PyG test must use a pinned, isolated environment and must record CPU loading,
single/batch latency, memory, determinism, and build provenance.

All machine-readable decisions, source evidence, source/asset identities,
endpoint mappings, license classifications, ARM64 status, dataset status, and
the conditional Stage 4E-3 plan are in `validation/stage4e2_*.json`.
