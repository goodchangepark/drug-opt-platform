# Stage 3 Step 3C model selection

Research frozen on 2026-08-25. Inhibitor and substrate status were treated as separate
endpoint/assay contracts. Scientific provenance and external behavior were weighted above the
mere presence of a downloadable checkpoint.

## Candidate review

| Candidate | CYP coverage | Validation/provenance | License and ARM64 | Decision |
| --- | --- | --- | --- | --- |
| Admetica commit `d4f70569901c189f39fa37871e2aeabeef3adc83` | Individual Chemprop checkpoints for five inhibitors and five substrates | Exact released training CSVs and publisher metrics are available. Five inhibitors and CYP2C9/2D6/3A4 substrates have reported metrics; CYP1A2/2C19 substrate checkpoints are omitted from the publisher model table | MIT repository; official Chemprop conversion and CPU inference verified on aarch64 | Selected for five inhibitors and three documented substrates. CYP1A2/2C19 substrate disabled. |
| ADMET-AI v2 | TDC-derived five inhibitor and three substrate tasks | Re-trained v2 differs from the published/v1 web model; useful package, but version-to-publication alignment is less direct for this installation | MIT; CPU compatible | Not selected; it largely covers the same TDC tasks without improving endpoint provenance here. |
| TDC | Five Veith inhibitor datasets; CarbonMangels CYP2C9/2D6/3A4 substrate tasks | Public data/benchmark interface, not one canonical pretrained deployment checkpoint | Apache-2.0 framework; dataset terms remain source-specific; ARM64 feasible | Used to corroborate endpoint availability and definitions, not as a deployed model. |
| DeepChem / generic Chemprop examples | Frameworks capable of training CYP classifiers | No single public, endpoint-versioned checkpoint selected for all requested contracts | MIT; ARM64 feasible | Framework only; not selected as model evidence. |

Sources: [Admetica](https://github.com/datagrok-ai/admetica),
[ADMET-AI](https://github.com/swansonk14/admet_ai),
[TDC ADME tasks](https://tdcommons.ai/single_pred_tasks/adme/),
[PubChem AID 1851](https://pubchem.ncbi.nlm.nih.gov/bioassay/1851), and
[Carbon-Mangels et al.](https://doi.org/10.1002/minf.201100069).

The packaged model version is `admetica-d4f7056-cyp-chemprop-v2.1`. Original Chemprop 2.0
checkpoints were converted, not retrained, with Chemprop 2.2.4's supported v2.0-to-v2.1 converter.
The aarch64 host ran every checkpoint with PyTorch CPU inference. The Admetica repository and
checkpoints are distributed under MIT, which permits commercial/internal R&D use under its terms;
the PubChem/literature training sources retain their own terms, which should be reviewed separately
before model redistribution or product commercialization.

## Implemented endpoints and publisher validation

All inhibitor models are binary classifiers of the human recombinant CYP functional response in
PubChem AID 1851. The assay uses a CYP-specific pro-luciferin substrate at Km, NADPH, a 60-minute
enzyme incubation, and luminescence readout. AID 1851 defines active at `AC50 <= 10 µM` from a
40 µM to 0.24 nM test range. Model probability is the positive-class probability. It is not IC50,
AC50, Km, or a turnover rate and is never converted to one.

| Endpoint | Reported training N | Specificity | Sensitivity | Accuracy | Balanced accuracy |
| --- | ---: | ---: | ---: | ---: | ---: |
| CYP1A2 inhibitor | 13,239 | 0.873 | 0.866 | 0.870 | 0.869 |
| CYP2C9 inhibitor | 12,881 | 0.830 | 0.819 | 0.826 | 0.824 |
| CYP2C19 inhibitor | 13,427 | 0.819 | 0.830 | 0.824 | 0.825 |
| CYP2D6 inhibitor | 13,898 dataset table/released data; publisher model-table size cell 11,127 | 0.866 | 0.751 | 0.843 | 0.808 |
| CYP3A4 inhibitor | 12,997 | 0.815 | 0.842 | 0.826 | 0.829 |

The root overview's conservative CYP2C9 balanced accuracy of 0.824 is used; the repository's
metabolism detail page gives a conflicting 0.890. Both that conflict and the CYP2D6 size-cell
conflict are exposed in model-registry provenance rather than silently resolved.

Substrate models are binary molecular substrate/non-substrate classifiers based on heterogeneous
literature compilations rather than one harmonized kinetic assay. They do not predict turnover,
Km, Vmax, or clearance.

| Endpoint | Training source / reported N | Specificity | Sensitivity | Accuracy | Balanced accuracy |
| --- | --- | ---: | ---: | ---: | ---: |
| CYP2C9 substrate | Carbon-Mangels compilation / 899 | 0.728 | 0.757 | 0.738 | 0.742 |
| CYP2D6 substrate | Carbon-Mangels and Zaretzki compilations / 941 | 0.749 | 0.769 | 0.753 | 0.759 |
| CYP3A4 substrate | Carbon-Mangels and Zaretzki compilations / 1,149 | 0.569 | 0.779 | 0.718 | 0.674 |
| CYP1A2 substrate | `MODEL_UNAVAILABLE` | — | — | — | — |
| CYP2C19 substrate | `MODEL_UNAVAILABLE` | — | — | — | — |

Although the upstream repository contains CYP1A2/2C19 substrate checkpoint files, it provides no
publisher validation row and insufficiently clear dataset/assay provenance for those endpoints.
They are registered as unavailable instead of producing unsupported results.

## AD, confidence, experimental evidence, and liability

Each installed endpoint has its own training index. AD combines nearest-training Morgan radius-2,
2048-bit Tanimoto similarity, `1 - similarity` chemical-space distance, and the MW, cLogP, TPSA,
HBD, HBA, and rotatable-bond training envelope. `IN_DOMAIN` requires similarity at least 0.40 and
all descriptors in range; `BORDERLINE` requires similarity at least 0.25 and no more than one
descriptor outside; otherwise the result is `OUT_OF_DOMAIN`. This remains a transparent heuristic,
not calibrated uncertainty.

Confidence never uses the predicted probability. `MEDIUM` would require in-domain chemistry,
publisher balanced accuracy at least 0.80, and a two-class independent set of at least 30 compounds
with balanced accuracy at least 0.70. Otherwise it is `LOW`; unavailable endpoints are effectively
`UNKNOWN` because no prediction exists. The current external inhibitor results do not reach that
gate, so all installed CYP outputs are conservatively `LOW`. No `HIGH` value is generated.

Experimental binary classification (`class`/`0/1`) is compared only with its matching isoform and
role. Quantitative inhibition evidence such as IC50 is displayed with the classification prediction
as `NOT_NUMERICALLY_COMPARABLE`; absolute/relative errors remain null. The fixed decision threshold
is 0.5. A positive inhibitor classification emits `Potential CYP… inhibition concern`, with the
threshold and screening-only basis stored in provenance. No LLM and no quantitative potency claim
is involved.

## Independent validation and limitations

`scripts/validate_stage3c_models.py` uses probabilities at the fixed 0.5 threshold and removes exact
canonical-isomeric-SMILES overlap with each packaged training set:

| Endpoint | Independent set / n | Removed overlap | AUROC | AUPRC | Balanced accuracy | Sensitivity | Specificity | MCC |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| CYP2C9 inhibitor | ChEMBL 30 / 464 | 0 | 0.5851 | 0.3813 | 0.5717 | 0.5111 | 0.6322 | 0.1324 |
| CYP2D6 inhibitor | ChEMBL 30 / 639 | 0 | 0.5999 | 0.4164 | 0.5657 | 0.5436 | 0.5878 | 0.1216 |
| CYP3A4 inhibitor | ChEMBL 30 / 788 | 0 | 0.6533 | 0.4471 | 0.6096 | 0.6527 | 0.5665 | 0.2015 |

External inhibitor generalization is materially weaker than the publisher validation and is
therefore exposed in model details and confidence policy. CYP3A4 substrate was checked against 24
FDA-approved tyrosine kinase inhibitors labeled as substrates: two exact training overlaps were
removed and 21/22 remaining positives were classified positive (sensitivity 0.9545). This is only a
directionality sanity check; with no negative class, AUROC, AUPRC, specificity, balanced accuracy,
and MCC are not identifiable and are not reported.

No qualified independent sets were secured for CYP1A2/CYP2C19 inhibitors or CYP2C9/CYP2D6
substrates in this step. Exact-overlap exclusion does not eliminate close analogue, series, assay,
or source overlap. Checkpoints are deterministic and publisher probability calibration is unknown.
The PubChem functional signal can be affected by substrate turnover or assay interference, while
literature substrate labels are protocol-heterogeneous. These models support screening hypotheses,
not clinical DDI decisions.
