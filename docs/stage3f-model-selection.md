# Stage 3 Step 3F safety model selection and integration audit

Research was frozen on 2026-08-26. Safety outputs are endpoint-specific binary screening
classifications. No probability is converted to IC50, clinical causality, QT prolongation, or a
single ADMET score. CPU loading and inference were reproduced on this aarch64 host with PyTorch
2.8.0 and Chemprop 2.2.4.

## Model decisions

| Endpoint | Selected model | Definition and training data | Released validation | License decision |
| --- | --- | --- | --- | --- |
| hERG | Admetica commit `d4f70569901c189f39fa37871e2aeabeef3adc83`, converted Chemprop v2.1 checkpoint | Human hERG/KCNH2 blocker-liability classification; Wang et al. heterogeneous aggregation, 22,249 reported / 22,248 valid structures | specificity 0.811, sensitivity 0.897, accuracy 0.885, balanced accuracy 0.854 | MIT code/checkpoint permits commercial/internal R&D; upstream data terms are source-specific |
| Ames | ADMET-AI v2.0.1 commit `c65bf0418e19c65d7228f9e40da5d0152aade756`, five Chemprop-v2 models | Salmonella bacterial reverse-mutation aggregate from four studies, TDC/Xu et al., N=7,255 | ADMET-AI-v2 five-fold held-out AUROC 0.8816, AUPRC 0.8958 | MIT code/checkpoints; TDC page has ambiguous upstream dataset notice (`Not Specified` plus CC BY link), so dataset redistribution/commercial terms need separate review |
| DILI | Same ADMET-AI v2 ensemble, task 13 | FDA/NCTR human clinical drug-level DILI association, N=475; not an in-vitro hepatotoxicity endpoint | ADMET-AI-v2 five-fold held-out AUROC 0.8815, AUPRC 0.8777 | Same license limitation as Ames |

Sources: [Admetica](https://github.com/datagrok-ai/admetica),
[ADMET-AI v2](https://github.com/swansonk14/admet_ai),
[TDC toxicity endpoint definitions](https://tdcommons.ai/single_pred_tasks/tox/),
[Wang et al. hERG data paper](https://doi.org/10.1021/acs.jcim.5b00695),
[Xu et al. Ames paper](https://doi.org/10.1021/ci300400a), and
[Xu et al. DILI paper](https://doi.org/10.1021/acs.jcim.5b00238).

The hERG source does not preserve one assay-mode field per row. The deployed endpoint therefore
does not claim to be purely binding or purely functional patch-clamp inhibition. Ames likewise
does not harmonize strain, metabolic activation, dose, and protocol. DILI is specifically a
clinical association label and is never merged with cell toxicity or mechanistic liver assays.

## Independent evidence and confidence

The OpenADMET ChEMBL-37 human hERG IC50 aggregate contained 7,977 structures. After removing
7,249 exact canonical-SMILES overlaps with the packaged training set, 728 remained (489 positive
at median IC50 <=10 µM, 239 negative). At threshold 0.5 the selected checkpoint produced AUROC
0.6669, AUPRC 0.7854, balanced accuracy 0.5442, sensitivity 0.9755, specificity 0.1130, and MCC
0.1844. Exact removal cannot exclude source, series, or assay-lineage overlap, and low specificity
is scientifically material. hERG confidence is therefore always `LOW`.

No source-lineage-independent, reusable public structure/label set qualified for Ames or DILI.
Their known-compound results are sanity checks only, not independent metrics. Ensemble disagreement
is retained but not treated as calibrated uncertainty. All three endpoints are consequently capped
at `LOW` confidence. AD is calculated independently from each endpoint's transparent structure
index using nearest Morgan/Tanimoto similarity, `1-similarity` chemical-space distance, and the
MW/cLogP/TPSA/HBD/HBA/rotatable-bond descriptor envelope. This is heuristic AD, not a calibrated
uncertainty model.

The public raw Ames table used for the transparent similarity index contains 7,278 rows, whereas
ADMET-AI v2 reports 7,255 valid training labels. The checkpoint does not expose exact fold
membership, so this index is an approximate source-space reference rather than proof that a query
is inside the exact training fold. DILI has 475 source rows matching the release count.

## Sanity and integration audit

The public-reference acceptance dataset checks hERG blocker dofetilide and lower-liability atenolol,
Ames-positive nitrofurantoin and negative caffeine, and DILI-concern troglitazone and lower-risk
buspirone. All six directions match at threshold 0.5. Exact/possible training overlap is recorded
per row, so these results establish pipeline behavior only. The acceptance CSV contains public
compounds and no proprietary structures.

Every serialized prediction exposes record type, model/name/version/source, endpoint/unit/species,
training dataset, license, validation, AD, confidence, timestamp, and CompoundVersion. Compatible
experimental classification has display precedence while predictions remain stored. Quantitative
hERG IC50 and classification probability stay numerically incomparable. Deterministic Strengths,
Concerns, and Unknown lists preserve each evidence source and confidence; no overall score exists.

Mitochondrial toxicity, generic cytotoxicity, skin sensitization, BBB penetration, and composite
CNS liability remain explicit `MODEL_UNAVAILABLE` records. BBB is not relabeled as toxicity.
Existing unavailable transporter and CYP endpoints remain visible with reasons.
