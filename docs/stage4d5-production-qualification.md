# Stage 4D-5 Production Qualification

Stage 4D-5 adds deterministic qualification governance without changing any visible prediction. The Stage 4D-4 endpoint strategy registry remains the production source of truth. Qualification code records evidence and evaluates whether a frozen candidate may be proposed for later manual activation; it does not select or execute production strategies.

The authoritative lifecycle is versioned as `stage4d5-lifecycle-v1`:

`RESEARCH_ONLY → SHADOW → VALIDATED → PRODUCTION_CANDIDATE → ACTIVE`

`ACTIVE → RETIRED` and `ACTIVE → ROLLED_BACK` are also supported. `SHADOW → ACTIVE` is illegal. Activation requires a `ManualPromotionAuthorization`; no public mutation endpoint exists.

## Evidence record

The append-only qualification store separates four records:

- frozen prediction, including compound version, endpoint contract, exact strategy, model/checkpoint identities, standardizer, policy, value/probability, and freeze time;
- linked experiment, including assay, species, unit, date, availability time, source, quality, and protocol metadata;
- structured qualification decision, including dataset snapshot hash, validation type/timing, metrics, limitations, decision, rollback target, and provenance;
- lifecycle event, including versioned transition, reason, manual authorization ID, and rollback reason.

SQLAlchemy rejects updates and deletes for these evidence rows. Derived monitoring statistics can therefore be rebuilt from the same source hashes.

## Promotion gate

Promotion is conjunctive. Contract compatibility, frozen identity, absence of leakage, independent evidence, endpoint-specific sample requirements, non-inferiority, documented benefit, calibration, subgroup robustness, safety trade-offs, rollback, reproducibility, and manual activation must all pass. One metric cannot promote a strategy.

Minimum sample requirements and practical-equivalence margins are endpoint- and policy-version specific. Stage 4D-5 does not invent values. Unapproved values remain unconfigured and fail closed as `INSUFFICIENT_EVIDENCE`.

## Candidate tracks

- Solubility: M1 remains active; adaptive M1/M2 remains shadow. M3 remains excluded.
- Caco-2: M1 remains active; alternative consensus remains shadow with insufficient evidence.
- CYP3A4 inhibition: M1 remains active; the fixed 95.78/4.22 blend remains shadow. Dynamic adaptation stays closed as `NO_ADAPTIVE_VALUE`.
- hERG: raw M1 remains active; Platt-calibrated M1 is a shadow qualification candidate. The secondary model is supporting-only and is not a discriminative blend member.

All current details are generated in `validation/stage4d5_qualification_policy.json`. The runtime registry validates that every Stage 4D-4 endpoint has exactly one qualification classification.

## Model-unavailable and mechanistic endpoints

`MODEL_UNAVAILABLE` endpoints are `EXCLUDED` and cannot create candidates. PK, NCA, IVIVE, allometry, simulation, deterministic derivations, and rule methods use `METHOD_QUALIFICATION`, not blind ML ensemble promotion semantics. Mechanisms are evaluated independently and are never raw-score averaged.

## Read-only API

- `GET /api/qualification/strategies`
- `GET /api/qualification/endpoint/{endpoint_id}`
- `GET /api/qualification/candidates`
- `GET /api/qualification/drift`

There is intentionally no public POST/PATCH/DELETE activation route.

