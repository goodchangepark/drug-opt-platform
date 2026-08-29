# Stage 4E-1 — Model Landscape & Qualification Planning

## Scope

This is a planning-only stage. No package, checkpoint, dataset, model
registry entry, prediction policy, calibration, threshold, shadow execution,
or production value changed. The Stage 4D-7 runtime remains authoritative.

## Current production and gaps

The complete reconciled 49-policy baseline is
`validation/stage4e1_current_model_baseline.json`. The highest-value gaps are
hERG discrimination, independent Caco-2 Papp evidence, species-specific
microsomal clearance validation, quantitative pKa/logD coverage, and missing
transporter endpoints.

Current production models are not roadmap candidates. Current shadows remain
research evidence only: ESOL for solubility, physchem Caco-2, Morgan CYP3A4,
and physchem hERG. None is activated here.

## Narrow qualification roadmap

The candidate landscape distinguishes code, checkpoint, and training-data
licenses. Missing license is not permission. It also separates endpoint
compatibility from architecture novelty and marks likely training overlap.

| Candidate | Potential value | Current action |
|---|---|---|
| CardioGenAI hERG discriminator | Diverse hERG secondary hypothesis | License/checkpoint/label lineage review first |
| MetaboGNN HLM/MLM | Different GNN clearance hypothesis | License and unit-contract review first |
| pkasolver-lite | Quantitative pKa coverage | ARM64 feasibility first |
| pKaLearn | Independent ionization-center GNN hypothesis | ARM64/checkpoint provenance first |
| BayeshERG | hERG BNN | No-go: released weights are non-commercial |
| MMTKPred | Transporter kinetics | No-go: Vmax/Km is not transporter classification |

The exact Stage 4E-2 pilot list is in
`validation/stage4e1_stage4e2_pilot_plan.json`. It is a qualification queue,
not an installation list. Every pilot must first clear license, checksum,
endpoint contract, training-overlap, scaffold holdout, paired bootstrap,
applicability-domain, and Xavier CPU benchmark gates.

## Dataset-first work

Datasets can be more valuable than a new model. TDC Caco2_Wang and TDC hERG
are useful public benchmarks but are not independent from Drug-OPT's
Wang-derived training lineages, so neither is approved for promotion evidence.
The existing quarantined Biogen prospective benchmark is the highest-value
independent validation resource. ExpansionRx has promising exact Caco-2/HLM/
MLM endpoint coverage but requires terms, raw-label, and overlap review.
The 1,130-compound logD7.4 dataset requires license review before acquisition.

## Foundation and web-only tools

A molecular foundation encoder without a released, validated endpoint head is
not a Drug-OPT predictor. It remains a watchlist item. Web-only services are
not candidates for local production dependency.

## Promotion boundary

Stage 4E-2 must use the Stage 4D-5 qualification gates. A candidate becomes a
shadow only after reproducible acquisition and qualification; no landscape
entry authorizes activation. Historical freezes and current production remain
unchanged.
