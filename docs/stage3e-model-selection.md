# Stage 3 Step 3E transporter model selection

Research frozen on 2026-08-25. Substrate and inhibitor labels were treated as different
endpoint contracts, and all enabled predictions are human-specific. A dataset or web predictor
was not treated as a redistributable pretrained checkpoint.

## Candidate review and decision

| Candidate | Transporter coverage | Public checkpoint / validation | License / ARM64 | Decision |
| --- | --- | --- | --- | --- |
| Admetica commit `d4f70569901c189f39fa37871e2aeabeef3adc83` | Separate P-gp inhibitor and substrate Chemprop files | Inhibitor has documented training counts and metrics. Substrate has a checkpoint/data but no reported metrics for that exact Chemprop model | MIT repository/checkpoints; Chemprop 2.2.4 CPU inference verified on aarch64 | Human P-gp inhibitor selected; substrate disabled |
| ADMET-AI v2 | TDC Pgp_Broccatelli, which is an inhibitor task | Bundled weights, but no requested transporter expansion beyond P-gp inhibitor | MIT; ARM64-capable | Not selected because it does not improve endpoint coverage/provenance |
| TDC | Pgp_Broccatelli inhibitor benchmark/data | Public benchmark/data rather than one endpoint-versioned pretrained checkpoint | Apache-2.0 framework; source data terms remain separate | Provenance corroboration only |
| vNN-ADMET | Distinct P-gp substrate model (822 compounds; 422 substrate/400 non-substrate) | Reported restricted-AD accuracy/specificity/sensitivity 0.80/0.80/0.80, but the download requires an account and no reusable local checkpoint/license was qualified | Web deployment; local ARM64 reproducibility not established | `MODEL_UNAVAILABLE` |
| ADMETlab 3.0 | BCRP inhibitor N=2,799; BSEP inhibitor N=763; OATP1B1 inhibitor N=2,372; OATP1B3 inhibitor N=2,228 | Web/API models; no public local checkpoint with clear redistribution license was identified | Local ARM64 reproducibility not established | `MODEL_UNAVAILABLE` |
| admetSAR 3.0 | Public endpoint collections include P-gp, BCRP, BSEP, OCT1/2, and MATE1 data | Web predictor/data resources; no endpoint-versioned local checkpoints and redistribution terms qualified | Local ARM64 reproducibility not established | `MODEL_UNAVAILABLE` |
| DeepChem / Chemprop / TDC training recipes | Frameworks can train classifiers | A framework is not a pretrained endpoint model, and retraining was outside this pretrained-model selection | Framework licenses are permissive; ARM64 feasible | Not used as scientific model evidence |

Sources: [Admetica](https://github.com/datagrok-ai/admetica),
[Admetica absorption model card](https://github.com/datagrok-ai/admetica/blob/main/ADMET/absorption/absorption.md),
[Broccatelli et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC3069647/),
[TDC Pgp_Broccatelli benchmark](https://tdcommons.ai/benchmark/admet_group/03pgp/),
[vNN-ADMET available models](https://vnnadmet.bhsai.org/vnnadmet/availablemodels.xhtml),
[ADMETlab 3.0 endpoint inventory](https://admetlab3.scbdd.com/diversity/), and
[admetSAR 3.0 resources](https://lmmd.ecust.edu.cn/admetsar3/resource/ADME.php).

## Active endpoint

`P-gp inhibitor` uses `admetica-d4f7056-pgp-inhibitor-chemprop-v2.1`. It is a binary
human P-glycoprotein/ABCB1 functional inhibitor classifier, never a substrate model. The
training compilation aggregates more than 60 literature sources with heterogeneous cell systems,
probe substrates, concentrations, and assay conditions. Broccatelli labels used IC50 <= 15 µM or
more than 25-30% inhibition for inhibitors and IC50 >= 100 µM or less than 10-12% inhibition for
non-inhibitors when those measurements were available. The model therefore does not represent one
standardized regulatory assay.

The reported source set contains 1,275 compounds (666 inhibitor, 609 non-inhibitor); the packaged
curated file contains 1,227 valid structures. Admetica reports specificity 0.916, sensitivity
0.863, accuracy 0.888, and balanced accuracy 0.889. Output is the checkpoint's positive-class
probability at a fixed 0.5 classification threshold. Calibration was not published, so it is not
described as calibrated and is never converted to Ki, IC50, efflux ratio, or any quantitative unit.

Original Chemprop 2.0 weights were converted, not retrained, with Chemprop 2.2.4's supported
v2.0-to-v2.1 converter. PyTorch CPU loading/inference was verified on the aarch64 host. Admetica is
MIT licensed, which permits internal and commercial R&D under the license; upstream literature/data
terms still require separate review for redistribution or commercial product use.

## Explicit unavailable endpoints

The registry independently records `MODEL_UNAVAILABLE` for P-gp substrate, BCRP substrate, BCRP
inhibitor, BSEP inhibitor, OATP1B1 inhibitor, OATP1B3 inhibitor, OCT1 inhibitor, OCT2 inhibitor,
MATE1 inhibitor, and MATE2-K inhibitor. The active human P-gp inhibitor result is never reused for
a substrate, another transporter, or another species. OCT/MATE role-unspecified predictions are not
fabricated; the registered roles correspond to the investigated public inhibitor resources.

## Applicability domain and confidence

The endpoint has its own 1,227-structure training index. AD combines nearest-training Morgan
radius-2, 2048-bit Tanimoto similarity, `1 - similarity` chemical-space distance, and a descriptor
envelope for MW, cLogP, TPSA, HBD, HBA, and rotatable bonds. `IN_DOMAIN` requires similarity >= 0.40
and every descriptor in range; `BORDERLINE` requires similarity >= 0.25 and at most one descriptor
outside; all other compounds are `OUT_OF_DOMAIN`. This is a transparent heuristic, not calibrated
model uncertainty.

Probability and confidence are separate. `MEDIUM` would require in-domain chemistry, publisher
balanced accuracy >= 0.80, and a two-class independent set of at least 30 compounds with balanced
accuracy >= 0.70. No rigorously independent public set qualified: accessible alternatives share
Broccatelli/Chen source lineage or lack reusable structures. Confidence is consequently capped at
`LOW`, regardless of probability. Unavailable endpoints have no prediction and effectively
`UNKNOWN` confidence.

## Experimental evidence, flags, and validation

Experimental records retain transporter, substrate/inhibitor role, assay, value/class, unit,
species, source, and notes through the existing ADMET measurement contract. Only matching human
P-gp inhibitor binary class evidence is compared to the predicted class. Quantitative IC50/Ki or
efflux evidence remains visible as `NOT_NUMERICALLY_COMPARABLE`; no absolute or relative numerical
error is invented. Non-human and substrate evidence is excluded. A positive result emits the
deterministic `Potential P-gp / ABCB1 inhibition concern` flag with its 0.5 rule and LOW confidence;
no LLM and no overall candidate ranking is used.

`scripts/validate_stage3e_model.py` reproduces a directionality sanity check for verapamil,
tariquidar, and zosuquidar. All three were classified inhibitor (probabilities 0.999378, 0.999809,
and 0.999983), but verapamil and zosuquidar are exact training structures and tariquidar has nearest
training similarity 0.9014. These are therefore training-lineage sanity results, not independent
validation; AUROC, AUPRC, balanced accuracy, sensitivity, specificity, and MCC are deliberately not
reported. This overlap and the absence of an independent set are persisted in the validation
artifact and model registry rather than hidden.
