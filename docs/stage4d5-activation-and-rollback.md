# Stage 4D-5 Activation and Rollback

Stage 4D-5 provides an internal service-level capability for controlled future activation. It does not expose an unauthenticated public mutation endpoint and does not activate any current shadow strategy.

## Candidate freeze

At `PRODUCTION_CANDIDATE`, model IDs, versions, checkpoint hashes, weights, threshold, calibration, policy, standardizer, and endpoint contract are identified by one canonical specification hash. Any change requires a new candidate ID and candidate version. Reusing an ID with a changed specification is rejected.

## Manual activation

Activation is allowed only from `PRODUCTION_CANDIDATE` and requires a structured manual authorization containing authorization ID, reviewer identity, timestamp, qualification record hash, and reason. The service captures the exact previous active candidate ID. The active production output is not switched by drift detection or by passing a gate alone.

## Deterministic rollback

Rollback moves the active candidate to `ROLLED_BACK` and restores the recorded previous active specification. It performs no recomputation and never searches for the “latest good” model.

Rollback reasons are `PERFORMANCE_REGRESSION`, `CALIBRATION_FAILURE`, `RUNTIME_FAILURE`, `DATA_QUALITY_ISSUE`, `POLICY_ERROR`, or `MANUAL_ROLLBACK`. Every lifecycle event is append-only and includes the state-machine version.

The frozen Stage 4D-4 rollback reference is `validation/stage4d5_production_baseline.json`. It records all 27 active policies, exact models/versions, thresholds, calibration states, policy versions, and existing rollback provenance. `production_behavior_changed` remains false.

## Synthetic dry run

`validation/stage4d5_dry_run.json` exercises the entire lifecycle using only labeled TEST data: active and shadow freezes, later results, compatibility, reproducible metrics, gate, candidate transition, simulated manual activation, simulated drift, and deterministic rollback. It does not read or mutate research data or production state.
