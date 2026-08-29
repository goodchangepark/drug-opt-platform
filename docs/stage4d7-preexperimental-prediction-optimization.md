# Stage 4D-7 — Pre-Experimental Prediction Optimization

## Scope and decision

Stage 4D-7 audits the global prediction made before an endpoint experiment is
available. It uses only the existing, versioned Stage 4D validation artifacts.
It does not use project observations, chemical-series adaptation, a compound's
later measurement, ALENIGLIPRON, retraining, or new downloaded models.

The review found no candidate eligible for activation. The Stage 4D-5 policy
requires endpoint-specific minimum sample/class balance requirements and a
versioned non-inferiority margin. Both remain intentionally unconfigured;
therefore they fail closed. Production policy stays
`stage4d4-endpoint-strategy-v1`, and all existing rollback metadata remains
valid.

## Data separation and anti-leakage controls

- Existing Stage 4D-3A2 solubility bootstrap, Stage 4D-3B1A CYP3A4
  attribution/calibration, Stage 4D-3B2A hERG scaffold holdout, and Stage
  4D-2C Caco-2 bootstrap are cited by content hash.
- No candidate was fit on, then reported against, the same observations in
  this stage. No new calibration, threshold, selector, or blend weight was
  fitted.
- hERG's historical calibration split remains 546 calibration / 182 untouched
  scaffold-aware holdout. Calibration probability quality is reported
  separately from discrimination.
- Applicability domain does not trigger model switching. In particular, hERG
  OOD N=2 cannot justify an OOD correction.
- Qualification freezes are append-only and were not mutated.

## Endpoint decisions

The machine-readable master matrix is
`validation/stage4d7_endpoint_accuracy_matrix.json`.

| Endpoint | Evidence reviewed | Decision |
|---|---|---|
| Solubility | N=250, 1,000 paired bootstrap | M1 retained; ESOL shadow retained; no domain selector |
| Caco-2 | N≈34, bootstrap CI crosses zero | M1 retained; data-limited accuracy ceiling |
| CYP3A4 inhibitor | M1/M2 fixed blend attribution | M1 retained; blend research only; dynamic adaptation no-go |
| hERG | N=728, scaffold calibration holdout | raw M1/0.50 retained; Platt research only; M2 shadow only |
| PPB / HLM / RLM / MLM | no promotion-grade independent correction evidence | single core retained, species isolated |
| Other CYP / P-gp / Ames / DILI | no endpoint-specific promotion-grade calibration evidence | raw single core retained |
| SoM / metabolites | rank fusion / rules | retained; not an ML-consensus promotion target |
| PK / NCA / IVIVE / allometry | mechanistic methods | mechanistic qualification only, no ML consensus |
| unavailable transporters, dog/monkey clearance, quantitative pKa/logD | no qualified model | MODEL_UNAVAILABLE retained |

## Known accuracy ceilings

- hERG: prior/label shift, assay heterogeneity, limited base discrimination,
  and an inadequate secondary model.
- Caco-2: very small evidence cohort does not support a robust selector.
- PPB and microsomal clearance: no independent evidence adequate for a
  bias-correction claim; HLM/RLM/MLM remain strictly separate.
- Transporters, dog/monkey clearance, and quantitative pKa/logD have model
  coverage gaps, not values that may be fabricated.

## Runtime and history

The Stage 4D-6 canonical orchestrator already executes the exact current
policy: CORE provides the visible production result and authorized shadows are
stored separately in immutable freeze provenance. Stage 4D-7 makes no runtime
policy change. Historical Stage 4D-6 freezes, including ALENIGLIPRON's OOD
record, remain reproducible under their original policy and model identities.

Future work is controlled internal compound validation to configure
endpoint-specific qualification thresholds, or a separately governed model
expansion/qualification stage for documented model gaps. Project/series
adaptation is out of scope.
