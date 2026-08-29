# Stage 4D-4 Promotion and Rollback Policy

Stage 4D-4 defines governance only. It does not promote Solubility adaptive weighting, Caco-2 consensus, the CYP3A4 fixed blend, hERG calibration, or hERG M2. Any such change belongs to a separate Stage 4D-5 validation and change-control cycle.

## Lifecycle state machine

The permitted lifecycle is:

```text
SHADOW → VALIDATED → PRODUCTION_CANDIDATE → ACTIVE
```

`DEFERRED` identifies unavailable strategies that cannot enter the lifecycle until a qualified endpoint-specific implementation exists. `FROZEN` can hold a policy when evidence or operational controls are incomplete. State transitions do not occur automatically and are never driven by one metric.

## Mandatory promotion gates

Every promotion requires all of the following:

1. Endpoint contract compatibility, including assay, species, role, unit, transformation, and output semantics.
2. Held-out validation with calibration/tuning data isolated from the final test set.
3. Documented absence of structure, scaffold, source, assay, and temporal leakage.
4. Meaningful improvement across more than one decision metric; a calibration-only gain cannot be represented as improved discrimination.
5. Stable calibration on an untouched holdout when probabilities or thresholds are involved.
6. Robustness across scaffold series and scientifically relevant subgroups.
7. Immutable model IDs, versions, checkpoints, preprocessing, and endpoint-contract identity.
8. Complete validation artifact and deterministic rollback metadata.

Endpoint-specific requirements in the master matrix are additive to these gates.

## Deterministic rollback record

Every `ACTIVE` primary policy contains:

- current `policy_version`;
- `previous_policy_version`;
- exact rollback target;
- rollback strategy;
- rollback model IDs and versions;
- promotion reason;
- validation artifact path;
- model/version provenance mapping.

Rollback selects the recorded policy and immutable model/version identities. It must not choose “latest,” discover a model dynamically, recompute weights from current feedback, or substitute another species/endpoint.

## Endpoint-specific rollback safety

- Solubility rolls back to the existing Admetica M1 core, never to an adaptive shadow.
- Caco-2 rolls back to the existing Admetica M1 core.
- CYP endpoints roll back by isoform and inhibitor/substrate role; CYP3A4 does not provide a fallback for another CYP endpoint.
- hERG rolls back to raw Admetica M1 at the unchanged production threshold. Platt calibration and M2 remain outside the rollback target until independently promoted.
- HLM, RLM, and MLM rollback targets remain species-specific.
- PK rollback preserves the mechanistic evidence hierarchy rather than averaging methods.
- `MODEL_UNAVAILABLE` endpoints have no active strategy to roll back and continue returning explicit unavailability.

## Rollback triggers

A promoted strategy must be rolled back when contract incompatibility, provenance drift, data leakage, calibration instability, subgroup regression, model/version mismatch, or an operational integrity failure is confirmed. The rollback event must preserve the failing policy version and evidence for audit; it must not delete research or project data.
