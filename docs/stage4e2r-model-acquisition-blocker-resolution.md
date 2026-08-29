# Stage 4E-2R — Model Acquisition Blocker Resolution

Stage 4E-2R resolves only acquisition, legal, endpoint, and technical gates.
It does not benchmark models, fit weights, register a shadow, or change any
production prediction policy.

## Resolved asset

`DATA_OPENADMET_EXPANSIONRX_CACO2_PAPP_AB` is the sole Stage 4E-3-eligible
asset. The official OpenADMET post-challenge release is CC-BY-4.0, has DOI
`10.57967/hf/9687`, and was acquired at immutable revision
`6b898ccc43d10d25b230fb09e22a6e30c30022b5`. The raw CSV is outside Git with
SHA256 `f674ec74cca1146bc386f832a32d4b8d921d3c312f92cb436cc005901c724a3c`.

The intake contains 7,618 valid standardized structures; 3,773 numeric and
33 censored Caco-2 Papp A→B observations. Censoring is preserved, not
converted to point estimates. It remains an external benchmark only: exact
canonical-SMILES exclusions, units/log-transform reconciliation, duplicate
handling, and assay-context review are mandatory before Stage 4E-3 metrics.

## Unresolved models

CardioGenAI and MetaboGNN remain blocked by missing separate asset rights and
endpoint/checkpoint provenance. pKaLearn remains blocked by data/weight lineage
and unproven ARM64 execution. pkasolver-lite was tested more deeply: an
isolated aarch64 environment successfully installed modern PyG, but the
source-contained checkpoint fails strict loading because its legacy GIN
state-dict layout does not match the current released PyG architecture. This
is a reproducibility failure, not a reason to relax checkpoint loading.

ADMET-AI v2 was identified as a narrow hERG replacement with MIT package and
bundled weights, but its pinned dependency resolution selects a CUDA-heavy
Torch stack on the target; no CPU runtime smoke was completed. It is not a
Stage 4E-3 candidate.

## Entry status

The result is `PARTIAL_READY_FOR_STAGE_4E3`: one qualified external Caco-2
dataset is available, but no eligible new model is available. Stage 4E-3 must
not start until a model and its independent evaluation path meet every gate.
Production remains unchanged.
