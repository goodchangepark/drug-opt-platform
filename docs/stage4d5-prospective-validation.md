# Stage 4D-5 Prospective Validation

A result is prospective only when the prediction was frozen before the experimental result became available. Re-running a model after result availability is deterministically classified `POST_HOC_PREDICTION` and cannot count toward prospective qualification.

## Validation taxonomy

The machine-readable taxonomy distinguishes `TRAINING_INTERNAL`, `CROSS_VALIDATION`, `EXTERNAL_RETROSPECTIVE`, `PSEUDO_PROSPECTIVE`, `PROSPECTIVE_INTERNAL`, `PROSPECTIVE_EXTERNAL`, `CLINICAL_RETROSPECTIVE`, and `CLINICAL_PROSPECTIVE`. They carry explicit evidence ranks and are never treated as interchangeable. Clinical labels describe evidence context only; Stage 4D-5 makes no clinical or regulatory validation claim.

## Freeze before result

Each immutable freeze contains:

- compound version and optional project/chemical-series scope;
- endpoint ID and endpoint contract version;
- candidate specification hash, strategy, model versions, and checkpoint hashes;
- value and/or probability, unit, applicability-domain state, standardizer, policy, provenance, and UTC timestamp.

The linked experiment includes endpoint contract, value, unit, assay type, species, experiment date, result-availability time, source, quality, and available protocol metadata.

Compatibility is checked before evaluation. Outcomes are `QUALIFICATION_ELIGIBLE`, `LIMITED`, `INCOMPATIBLE`, `MISSING_METADATA`, `NO_FROZEN_PREDICTION`, or `POST_HOC_PREDICTION`. Only `QUALIFICATION_ELIGIBLE` observations count.

## Performance rebuild

Eligible records are isolated by endpoint, model version, and policy version. Monitoring supports `GLOBAL`, `PROJECT`, `CHEMICAL_SERIES`, `MODEL_VERSION`, and `POLICY_VERSION` scopes.

Regression endpoints compute N, MAE, RMSE, bias, Spearman correlation, and within-two-/three-fold rates. Binary classification endpoints compute N, MCC, balanced accuracy, AUROC, AUPRC, Brier score, log loss, sensitivity, specificity, and ECE. Ranking and mechanistic endpoints do not receive inappropriate scalar ML metrics.

The source-record hash list is sorted before hashing. Clearing derived TEST statistics and rebuilding from the same immutable observations produces the same metrics and source snapshot hash.

## Priority prospective plans

hERG requires frozen raw-M1 and calibrated-M1 probabilities linked to compatible experimental hERG labels, class balance, sensitivity, specificity, MCC, balanced accuracy, Brier, log loss, ECE, applicability, and subgroup analysis. No results currently exist and none are fabricated.

Solubility requires future compatible project measurements comparing the adaptive shadow with active M1 on the same prospective compounds. Caco-2 similarly accumulates compatible future A→B Caco-2 measurements. CYP3A4 collects evidence for the fixed blend only; dynamic adaptation stays closed until a new explicit scientific decision reopens it.

