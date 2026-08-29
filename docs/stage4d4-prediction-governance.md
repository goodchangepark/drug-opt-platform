# Stage 4D-4 Prediction Governance

Stage 4D-4 is a policy finalization release. It adds a read-only scientific strategy registry and does not alter the prediction values returned by any existing production path. The authoritative implementation is `backend/endpoint_strategy_registry.py`; the machine-readable projection is `validation/stage4d4_endpoint_strategy_matrix.json`.

## Policy boundaries

Each endpoint owns one explicit production strategy. Evidence from one endpoint is never generalized to another endpoint, CYP isoform, substrate/inhibitor role, safety target, assay, or species. In particular, CYP3A4 results do not govern other CYP endpoints, and hERG evidence does not govern Ames or DILI.

The registry keeps five scientific concepts separate:

1. **Model confidence** describes what the model itself reports and is not automatically a calibrated probability.
2. **Applicability domain** describes chemical-domain support independently of confidence.
3. **Model disagreement** is a `MODEL_DISAGREEMENT_SIGNAL`; it is not a confidence interval.
4. **Calibration quality** follows the lifecycle `RAW` → `CALIBRATION_RESEARCH` → `CALIBRATION_VALIDATED` → `CALIBRATION_PRODUCTION`.
5. **Validation status** states the available external evidence and its limitations.

No aggregate “confidence score” may collapse these concepts.

## Production and shadow isolation

Primary policies in `ACTIVE` state codify the behavior already present before Stage 4D-4. A shadow strategy has lifecycle state `SHADOW`, cannot overwrite primary values, and cannot become active through the registry API.

- Solubility adaptive weighting remains `ADAPTIVE_RESEARCH_SHADOW`; M3 remains excluded from adaptive execution.
- Caco-2 static consensus remains shadow because promotion evidence is insufficient.
- CYP3A4 fixed global blending remains research shadow. Dynamic adaptive weighting has no demonstrated value beyond the fixed global prior.
- hERG Platt-calibrated M1 remains calibration research. The calibration-selected blend was 100/0, so M2 is calibration-supporting/shadow-only and is not a discriminative production blend.
- Ames and DILI remain independent single-core policies without extrapolation from hERG.

## Non-ML and unavailable endpoints

Metabolic soft spots use reciprocal rank-fusion semantics for compatible SyGMa and SMARTCyp ranks. Raw model scores are never averaged.

PK NCA, IVIVE, allometry, and simulation retain an evidence hierarchy because their assumptions are mechanistically distinct. They never enter an ML consensus.

The ionization pKa path is `RULE_ESTIMATE`. Conditional logD7.4 from cLogP, pKa, and Henderson-Hasselbalch assumptions is `DERIVED_ESTIMATE`. Neither is a validated quantitative ML prediction. The separate quantitative-ML registry slots remain `MODEL_UNAVAILABLE`.

Every `MODEL_UNAVAILABLE` endpoint has no primary model identity, cannot execute a production prediction, and must return `MODEL_UNAVAILABLE` rather than reuse another endpoint or synthesize a value.

## API contract

`GET /api/model-strategy-registry` is read-only and returns registry version, policy date, consistency violations, and endpoint records containing primary/shadow models, calibration state, adaptive and consensus states, validation status, limitations, lifecycle state, and rollback metadata.

The route neither calls prediction code nor modifies model-registry, database, or project state. Existing APIs and response schemas are unchanged.

## Change control

Registry or endpoint-contract changes must regenerate the matrix with:

```bash
.venv/bin/python scripts/generate_stage4d4_strategy_artifacts.py
```

Tests compare the generated artifact to the runtime registry so hand-edited policy drift is rejected.
