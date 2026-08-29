# Stage 4D-5 Drift Monitoring

Drift monitoring is deterministic and advisory. A warning produces `REVIEW_REQUIRED`; its automatic action is always `NONE`.

The warning taxonomy is:

- `PERFORMANCE_DRIFT`
- `CALIBRATION_DRIFT`
- `DOMAIN_SHIFT`
- `PRIOR_SHIFT`
- `INSUFFICIENT_DATA`
- `MODEL_VERSION_CHANGED`
- `ENDPOINT_MISMATCH`

Versioned monitoring thresholds are engineering review triggers, not scientific acceptance or regulatory limits. They compare current and frozen reference metrics, calibration error, out-of-domain rate, chemical-distance summary, label prevalence, model identity, and endpoint contract identity.

Warnings can never retrain a model, change a decision threshold, change ensemble weights, activate calibration, or switch the production strategy. Any follow-up qualification must create reproducible evidence and pass the full promotion gate.

Metrics are maintained separately from model confidence, applicability domain, model disagreement, calibration quality, and validation status. A disagreement signal is not a confidence interval.

The authoritative policy and explicit forbidden actions are in `validation/stage4d5_drift_policy.json`.

