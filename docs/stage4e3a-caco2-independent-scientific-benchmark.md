# Stage 4E-3A — Caco-2 Independent Scientific Benchmark

Stage 4E-3A evaluated the already deployed Caco-2 CORE and the already
authorized SHADOW model against the pinned ExpansionRx Caco-2 Papp A→B
release. It fitted no weights, calibration, selector, threshold, or numerical
consensus, and did not change production policy or runtime execution.

The source is the OpenADMET ExpansionRx release (revision
`6b898ccc43d10d25b230fb09e22a6e30c30022b5`, SHA256
`f674ec74cca1146bc386f832a32d4b8d921d3c312f92cb436cc005901c724a3c`),
licensed CC-BY-4.0. The endpoint is experimental Caco-2 Papp A→B. Raw values
in `10^-6 cm/s` were harmonized as `log10(Papp × 1e-6)` in `log10(cm/s)`.

The fixed models were `admetica_caco2`
(`admetica-d4f7056-chemprop-v2.1`) as CORE and `physchem_caco2_v1`
(`physchem-caco2-v1.0`) as SHADOW. The authoritative policy has no numeric
Caco-2 consensus; none was created.

## Cohort and overlap controls

The primary cohort is unique-molecule level: positive numeric Papp values are
aggregated by a median on the raw scale before the single log transform. It
contains 3,498 unique molecules / paired model results from 3,771 positive
numeric observations. The intake separately preserves 33 `SOURCE_CENSORED`
observations and excludes two numeric-zero records as
`NON_POSITIVE_PAPP_EXCLUDED / NON_QUANTITATIVE_ZERO`; neither class was
imputed or used as an exact primary regression target.
Their source-row IDs and distinct reasons are retained in
`validation/stage4e3a_caco2_exclusions.json`.

No canonical Caco-2 training/reference structure collection was locally
available for Wang/TDC/Admetica/previous validation lineage, so known exact
overlap removal was zero and `RESIDUAL_TRAINING_OVERLAP_UNKNOWN` remains an
explicit limitation. ExpansionRx is now benchmark-use data; future models
tuned against it may not describe it as untouched external validation.

## External results

On the paired primary cohort, CORE had MAE 0.5695 and RMSE 0.7457; SHADOW had
MAE 0.7047 and RMSE 0.9202. CORE was within two-fold for 34.0% and within
three-fold for 52.8%, compared with 31.0% and 45.5% for SHADOW. The paired
1,000-replicate bootstrap for SHADOW minus CORE MAE was +0.1352 (95% CI
+0.1183 to +0.1513; P(SHADOW better) = 0). Scaffold-clustered resampling
also retained the direction of worse SHADOW MAE.

CORE performance remained lower-error in the observation-level duplicate
sensitivity analysis. The current AD labels were stratified unchanged. Model
disagreement had Spearman 0.022 with CORE absolute error, so this external
cohort does not establish it as a meaningful CORE-failure warning signal.
SHADOW had case-level wins, but not enough aggregate or robust complementary
numeric value to support a different production number.

## Decision

`CURRENT_CORE_CONFIRMED`. Production remains `SINGLE_CORE_MODEL` with
`admetica_caco2`. No shadow activation, promotion, calibration, consensus, or
threshold change occurred. The observed result is
`OBSERVED_EXPANSIONRX_EXTERNAL_PERFORMANCE`, not universal model accuracy.
The next Caco-2 scientific need remains a genuinely qualified model expansion
with independent validation controls.

Machine-readable provenance and all detailed analyses are in the
`validation/stage4e3a_caco2_*` artifacts.
