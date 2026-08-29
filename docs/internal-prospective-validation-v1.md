# Engine v1 Internal Prospective Validation

## Policy

**Policy ID:** `drugopt-prediction-engine-v1`  
**Policy Version:** `drugopt-prediction-engine-v1@1.0.0`  
**Policy Hash:** `12757ab197b5a70d8ea1754678d9a342ab0b6ea0d82f2896bebb767d686bbdeb`  
**Standardizer:** `CHEM_STANDARDIZER_V1`  
**Stage:** 6 / 6 Pre-AI Engine Completion Tasks

---

## Scientific Principle

```
STRUCTURE
  → ENGINE V1 PREDICTION
  → IMMUTABLE FREEZE
  → EXPERIMENT
  → EXPERIMENTAL RESULT
  → COMPARISON
```

**Never:**

```
STRUCTURE
  → EXPERIMENTAL RESULT
  → MODEL TUNING
  → "PROSPECTIVE" PREDICTION
```

This distinction is mandatory and enforced at the code level.

---

## Absolute Scientific Rules

The following are **forbidden** during this validation stage:

- Retrain, fine-tune, or recalibrate any model
- Adjust thresholds, fit bias correction, fit ensemble weights
- Modify AD thresholds or confidence rules
- Change endpoint strategy or replace models
- Activate shadow models in production
- Use internal validation data to improve current predictions
- Implement project, series, or local adaptation
- Start AI Scientist

Internal data are **validation data** in this stage — not training data.

---

## Campaign

**Campaign ID:** `IVC-engine-v1-2026-08-29`  
**Protocol ID:** `internal-validation-v1-protocol-2026-08-29`  
**Framework Status:** `READY`  
**Scientific Status:** `COLLECTING`

---

## Cohort

| Compound | Series | Project | InChIKey |
|----------|--------|---------|----------|
| ORFORGLIPRON | GLP1-SM-PYRIDINONE | GLP1-SM | USUWIEBBBWHKNI-KHIFEHGGSA-N |
| ALENIGLIPRON | GLP1-SM-PYRIDINONE | GLP1-SM | CPOJUYUGONJVPZ-WIXASUBBSA-N |
| ELECOGLIPRON | GLP1-SM-PYRIDINONE | GLP1-SM | JMKBTILBGROESC-WNJJXGMVSA-N |

**Total:** 3 compounds, 1 chemical series, 1 project

> **Note:** 3 compounds is below the preferred 20–50 compound target. Data collection continues. The validation framework is fully operational and will produce analysis results as additional compounds and experimental data are added.

---

## True Prospective vs Blinded Retrospective Distinction

Every paired observation is classified into one of three evidence classes:

### TRUE_PROSPECTIVE
- **Definition:** `freeze_timestamp < result_available_at`
- **Use in primary metrics:** Yes
- **Meaning:** Prediction was frozen before the experimental result was available to the team.

### BLINDED_RETROSPECTIVE
- **Definition:** Experimental result existed historically but was hidden during prediction freeze; blinding is documented.
- **Use in primary metrics:** Yes (secondary supporting evidence)
- **Requirement:** Blinding must be explicitly documented.

### HISTORICAL_VISIBLE
- **Definition:** Experimental result was available and potentially visible during Engine v1 model development.
- **Use in primary metrics:** No
- **Handling:** Kept separate as retrospective context only. Must **NOT** be described as prospective validation.

> **Critical:** These three classes must **never** be merged. Only TRUE_PROSPECTIVE and documented BLINDED_RETROSPECTIVE constitute valid prospective evidence.

---

## Endpoint Contracts

Engine v1 covers 49 endpoints across the following strategies:

| Strategy | Count |
|---------|-------|
| SINGLE_CORE_MODEL | 18 |
| MODEL_UNAVAILABLE | 22 |
| RULE_ESTIMATE | 1 |
| DERIVED_ESTIMATE | 2 |
| RULE_BASED | 2 |
| MECHANISTIC_NO_CONSENSUS | 4 |

### Key Assay Compatibility Rules

**Solubility (`solubility_aqueous_logs`)**
- Unit: `log10(mol/L)` or convertible from `mol/L`
- Kinetic vs thermodynamic distinction recorded
- pH 7.4 preferred

**Caco-2 (`permeability_caco2_logpapp`)**
- Assay direction: **A→B required**
- Unit: `log10(cm/s)` or `cm/s` (convertible)
- Only Papp A→B compared

**PPB (`ppb_human_percent_bound`)**
- Species: **human only**
- Unit: `% bound` or `fu` (conversion: `% bound = (1-fu)*100`)
- Strict species isolation — rat PPB not compared

**HLM (`hlm_intrinsic_clearance_scaled_log10`)**
- Species: **human only**
- Assay type: HLM microsomal
- Hepatocyte data NOT comparable

**hERG (`safety_herg_blocker_prob`)**
- IC50 threshold must be specified (10 µM standard)
- Different thresholds → ASSAY_CONTEXT_LIMITED

---

## Prediction Freeze

- **Total freezes:** 18 (ALENIGLIPRON, 18 production endpoints)
- **Policy hash verified:** Yes
- **Experimental data hidden before prediction:** Yes (enforced by code)
- **Freeze immutability:** Append-only (never overwritten)
- **Upstream freeze table:** `qualification_prediction_freezes`
- **Validation freeze table:** `internal_validation_prediction_freezes`

### Historical Freeze Verification (ALENIGLIPRON)

These values are immutable and verified in tests:

| Endpoint | Frozen Value |
|----------|-------------|
| `solubility_aqueous_logs` | −4.287727355957031 |
| `permeability_caco2_logpapp` | −5.135347366333008 |
| `cyp3a4_inhibitor_prob` | 0.9331268668174744 |
| `safety_herg_blocker_prob` | 0.9903563857078552 |

---

## Experimental Import

### Import Template

`validation/internal_validation_v1_experiment_import_template.csv`

### Import Script

```bash
source .venv/bin/activate

# Validate CSV (dry run)
python3 scripts/import_validation_experiments.py \
    --csv your_experimental_data.csv --dry-run

# Import
python3 scripts/import_validation_experiments.py \
    --csv your_experimental_data.csv
```

### Mandatory CSV Fields

| Field | Description |
|-------|-------------|
| `compound_id` | ORFORGLIPRON / ALENIGLIPRON / ELECOGLIPRON |
| `compound_version_id` | 1 / 2 / 3 |
| `endpoint_id` | Engine v1 endpoint ID |
| `raw_value` | Numeric measured value (never pre-log-transformed) |
| `raw_unit` | Exact unit string |
| `qualifier` | `=`, `<`, `>`, `BLQ`, `ULOQ` |
| `species` | human / rat / mouse |
| `assay_type` | Assay description |
| `assay_direction` | A→B (for Caco-2) |
| `result_available_at` | ISO 8601 — critical for TRUE_PROSPECTIVE |

---

## Replicates

- All replicates retained in `internal_validation_experimental_records`
- Aggregation (median / majority vote) applied at analysis time
- Raw replicates never cherry-picked

---

## Censored Values

- `qualifier: "<"`, `">"`, `"BLQ"`, `"ULOQ"` → `censor_flag=True`
- Censored observations excluded from primary quantitative metrics
- Used in directional/censored sensitivity analysis only
- **Never replaced with threshold or epsilon**

---

## Non-Positive Log Handling

- Experimental raw values ≤ 0 cannot be log-transformed
- Status: `NON_POSITIVE_EXCLUDED_FROM_LOG_METRIC`
- These observations are excluded from log-scale primary metrics
- **Never silently floored**

---

## Metrics

### Primary Regression Metrics
`N | MAE | RMSE | Bias | Median_AE | Spearman | P75_AE | P90_AE | P95_AE`

### Fold Accuracy (log10 endpoints)
`Within_2fold (log10(2) = 0.301) | Within_3fold (log10(3) = 0.477)`

### Primary Classification Metrics
`MCC | Balanced_Accuracy | AUROC | AUPRC | Sensitivity | Specificity | Brier | LogLoss | ECE`

### Bootstrap
- N replicates: 1,000
- Seed: 42 (fixed for reproducibility)
- Level: compound
- Minimum N: 10 (smaller → `data_insufficiency_flag=True`)

### Minimum N Rules
- Regression metrics: N ≥ 5
- Classification: N ≥ 5 per class for AUROC, N ≥ 10 for MCC
- Bootstrap: N ≥ 10
- Series Spearman: N ≥ 3 compounds per series

---

## Applicability Domain Validation

For each endpoint, observations stratified by AD class:
- `IN_DOMAIN`
- `BORDERLINE`
- `OUT_OF_DOMAIN`

**Key question:** Does error actually worsen as AD becomes weaker?

GLP-1 compounds are predominantly `OUT_OF_DOMAIN` or `BORDERLINE` across all endpoints — this is itself a key finding confirming the engine's conservative AD assessment for large, complex GLP-1 receptor modulators.

---

## Reliability Validation

Group predictions by Engine v1 reliability class and compare observed accuracy.

**Key question:** Does current reliability classification correlate with observed accuracy?

---

## Shadow Disagreement

Engine v1 GLP-1 predictions use `SINGLE_CORE_MODEL` strategy. No authorized shadow models were active for these compounds. Shadow disagreement analysis will apply when shadow-enabled endpoints receive internal data.

**Key question (for future endpoints):** Does `|CORE - SHADOW|` predict CORE error?

---

## Scaffold / Series Analysis

### Overall Compound-Level Results
Reported when N ≥ 1 endpoint pairs available.

### Within-Series Spearman
Reported when ≥ 3 compounds in the same series share an endpoint. Minimum N for within-series Spearman analysis.

### Series Summary: GLP1-SM-PYRIDINONE
All 3 enrolled compounds belong to the same Murcko scaffold family — limiting diversity but enabling within-series ranking analysis once experimental data arrives.

---

## Current Results (Data Collection Phase)

> **No experimental data has been imported yet.**  
> The validation framework is fully operational. Results will be populated as assays complete.

### Evidence Classification
| Type | N |
|------|---|
| TRUE_PROSPECTIVE | 0 |
| BLINDED_RETROSPECTIVE | 0 |
| HISTORICAL_VISIBLE | 0 |

### Endpoint Coverage

All 49 Engine v1 endpoints are documented. 22 are `MODEL_UNAVAILABLE` (coverage gaps recorded for AI Scientist). 18 production endpoints have frozen predictions for ALENIGLIPRON.

---

## Known Engine v1 Limitations

| Endpoint | Known Issue |
|----------|-------------|
| Caco-2 | External MAE 0.5695 log10; Spearman 0.041. Limited ranking ability. |
| hERG | AUROC 0.6669; Specificity 0.113 at threshold 0.50. Raw probability miscalibrated. |
| HLM/RLM | LOW-MEDIUM evidence. No independent benchmark. |
| MLM | INSUFFICIENT_EVIDENCE. No independent benchmark. |
| pKa | RULE_ESTIMATE. ±1–2 pKa unit uncertainty. |
| logD | DERIVED_ESTIMATE from cLogP + pKa. Not quantitative ML. |

---

## Coverage Gaps (MODEL_UNAVAILABLE)

22 endpoints have no production model in Engine v1. Internal experimental data for these endpoints:
- Recorded as `COVERAGE_GAP` evidence
- Not compared quantitatively  
- Preserved for future Engine v2 AI Scientist

Key gaps: BBB penetration, CYP substrate panels (1A2, 2C19), transporter panel (BCRP, BSEP, MATE1/2K, OATP1B1/B3, OCT1/2, P-gp substrate), quantitative pKa ML, quantitative logD ML.

---

## Final Scientific Decision

**Current decision:** `INTERNAL_VALIDATION_NOT_STARTED_AWAITING_EXPERIMENTAL_DATA`

This will be updated to one of the following as data arrives:

| Code | Decision |
|------|----------|
| A | `INTERNAL_VALIDATION_SUPPORTS_ENGINE_V1_BASELINE` |
| B | `INTERNAL_VALIDATION_MIXED_ENGINE_V1_RETAINED_WITH_LIMITATIONS` |
| C | `INTERNAL_VALIDATION_IDENTIFIES_ENDPOINT_REVIEW_REQUIRED` |
| D | `INTERNAL_VALIDATION_INSUFFICIENT_DATA_CONTINUE_COLLECTION` |
| E | `INTERNAL_VALIDATION_NOT_STARTED_AWAITING_EXPERIMENTAL_DATA` ← current |

> C or D does **NOT** automatically modify Engine v1.

---

## Engine v1 Production Decision

**UNCHANGED**

Engine v1 remains frozen per policy `drugopt-prediction-engine-v1@1.0.0`.  
Any scientifically necessary modification would require a future Engine v1.1 / Engine v2 process.

---

## Future Project / Series Adaptation Readiness

**Status:** `DATA_LIMITED`

3 compounds from 1 series and 1 project is insufficient for meaningful project- or series-level adaptation. 20–50 compounds across multiple series would be required.

**No adaptation implemented in this stage.**

---

## AI Scientist Readiness

**Status:** `PARTIAL`

AI Scientist can receive:
- ✅ Prediction value
- ✅ Prediction timestamp (freeze)
- ✅ Evidence class (MODEL_PREDICTION, RULE_ESTIMATE, etc.)
- ✅ Applicability domain
- ✅ Reliability
- ✅ Limitations
- ✅ MODEL_UNAVAILABLE documented
- ⏳ Experimental result (pending data)
- ⏳ Paired validation performance (pending data)
- ❌ Disagreement (no shadow models for these compounds)

> AI Scientist must **not** invent missing numeric predictions.

---

## Artifacts

All artifacts in `validation/internal_validation_v1_*.json`:

| Artifact | Status |
|----------|--------|
| `protocol.json` | ✅ Frozen |
| `campaign.json` | ✅ Active |
| `dataset_flow.json` | ✅ Current |
| `endpoint_contracts.json` | ✅ All 49 endpoints |
| `prediction_freezes.json` | ✅ 18 freezes |
| `experimental_manifest.json` | ✅ Awaiting data |
| `pairing_audit.json` | ✅ Framework ready |
| `metrics.json` | ⏳ Awaiting data |
| `bootstrap.json` | ⏳ Awaiting data |
| `ad_analysis.json` | ⏳ Awaiting data |
| `reliability_analysis.json` | ⏳ Awaiting data |
| `shadow_disagreement.json` | ✅ N/A for GLP-1 |
| `scaffold_series_analysis.json` | ⏳ Awaiting data |
| `final_decision.json` | ✅ Decision E |

---

## Progress

```
5 / 6 COMPLETE — TASK 6 VALIDATION CAMPAIGN ACTIVE
```

- ✅ Validation protocol frozen
- ✅ Campaign created
- ✅ Cohort enrolled (3 GLP-1 compounds)
- ✅ Prediction freezes registered (18 freezes)
- ✅ Blinding enforcement active
- ✅ All artifacts generated
- ✅ All validation tests passing
- ⏳ Experimental data collection ongoing
- ⏳ Paired analysis pending

**Next step:** `CONTINUE_INTERNAL_EXPERIMENTAL_DATA_COLLECTION`

To mark 6/6 scientifically complete: actual paired prediction/experimental observations  
from multiple compounds, multiple endpoints, and multiple chemical series are required.
