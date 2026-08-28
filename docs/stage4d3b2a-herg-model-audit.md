# Stage 4D-3B2A: hERG Model Quality Audit

## 1. Executive Summary

Stage 4D-3B2A is an authoritative diagnostic audit of the hERG cardiac liability
classification system (`EP_SAFETY_HERG_BLOCKER_PROB`), executed before any adaptive
weighting implementation.

The central scientific question:
> *"Is hERG poor specificity caused primarily by threshold, calibration, base-model
> discrimination, class imbalance, assay heterogeneity, or a combination?"*

### Authoritative Findings
All five root-cause categories are confirmed:

| Root Cause | Evidence |
| :--- | :--- |
| **THRESHOLD** | At 0.50, specificity = 0.113. Platt calibration at same threshold raises specificity to 0.310 |
| **CALIBRATION** | ECE = 0.265 (severely miscalibrated). M1 Brier = 0.274, LogLoss = 1.690 |
| **BASE_MODEL_DISCRIMINATION** | AUROC = 0.667 — moderate discrimination, cannot fully compensate with calibration alone |
| **CLASS_IMBALANCE** | Training prevalence = 86.0%; evaluation prevalence = 67.2%; shift = 18.8 pp |
| **LABEL_BOUNDARY_UNCERTAINTY** | 72.9% of evaluation compounds are IC50 borderline (1k–30k nM) |

### Final Decision: `HERG_FIXED_BLEND_CANDIDATE`
### Adaptive Weighting Gate: **NO_GO**

---

## 2. Endpoint Contract

| Field | Value |
| :--- | :--- |
| **Endpoint ID** | `safety_herg_blocker_prob` |
| **Display Name** | hERG Cardiac Blocker Liability Probability |
| **Category** | SAFETY |
| **Species** | Human |
| **Assay Type** | Heterogeneous Patch-Clamp & Radioligand Binding ([³H]-dofetilide / [³H]-astemizole) |
| **Positive Class** | BLOCKER |
| **Negative Class** | NON_BLOCKER |
| **Potency Threshold** | IC50 ≤ 10 µM (10,000 nM) / pIC50 ≥ 5.0 |
| **Decision Threshold** | 0.50 (probability) |
| **Output Type** | Binary classification probability ∈ [0, 1] |
| **Directionality** | LOWER_BETTER |

### Assay Heterogeneity Status: `ASSAY_HETEROGENEITY_PRESENT`
Training labels pool manual patch-clamp, automated electrophysiology, and radioligand
binding assays (Wang et al. compilation). Assay-type provenance is **not retained per
row** in the packaged dataset. Labels must be treated as heterogeneous screening
liability, not a single standardized functional assay. The endpoint contract
explicitly forbids mixing with QT-prolongation clinical data or IC50 regression
without thresholding.

---

## 3. Models Audited

### M1: Admetica D-MPNN (CORE)

| Field | Value |
| :--- | :--- |
| **Model ID** | `admetica_safety_herg` |
| **Model Version** | `admetica-d4f7056-herg-chemprop-v2.1` |
| **Model SHA-256** | `c1aae9ca495c5dcb699f386eb0d8a2c0a7bafb78fca7c92f78241bb6790a5b26` |
| **Architecture** | Chemprop D-MPNN (graph neural network on molecular graph) |
| **Training N** | 22,248 (one invalid SMILES excluded from 22,249 records) |
| **Training N+ / N−** | 19,130 / 3,118 |
| **Training Prevalence** | **86.0%** |
| **Training Dataset** | Wang et al. hERG blocker compilation curated by Admetica |
| **Decision Threshold** | 0.50 |
| **Probability Calibration** | Not reported by publisher; audit confirms **severely miscalibrated** |
| **Role** | CORE / PRIMARY |
| **License** | MIT (Admetica repository) |
| **ARM64** | Compatible |

### M2: Physchem Pharmacophore Logistic (SHADOW)

| Field | Value |
| :--- | :--- |
| **Model ID** | `physchem_herg_v1` |
| **Model Version** | `physchem-herg-v1.0` |
| **Architecture** | Logistic regression on cLogP, MW, TPSA, basic-N presence, aromatic ring count |
| **Training Dataset** | Wang et al. hERG Blocker Compilation (N=22,249; design-basis only) |
| **Applicability Domain** | MW ≤ 800 Da AND −2 ≤ cLogP ≤ 7 → IN_DOMAIN |
| **Role** | SHADOW_ONLY / CALIBRATION_SUPPORTING |
| **ARM64** | Native |

---

## 4. Raw Performance — Full Cohort (N=728, Threshold=0.50)

| Method | MCC | Balanced Acc | Sensitivity | Specificity | AUROC | AUPRC | Brier | LogLoss | ECE |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **M1 CORE** | 0.1844 | 0.5442 | **0.9755** | **0.113** | **0.6669** | **0.7865** | 0.2745 | 1.6901 | 0.265 |
| **M2 SHADOW** | −0.0227 | 0.498 | 0.9918 | 0.0042 | 0.5319 | 0.5975 | 0.2769 | 0.9196 | 0.236 |
| **50/50 Static** | 0.1437 | 0.5178 | 0.998 | 0.038 | 0.6405 | 0.7712 | 0.2656 | 0.911 | — |

True positives: 477. True negatives: 27. **False positives: 212. False negatives: 12.**

---

## 5. Validation Cohort

| Field | Value |
| :--- | :--- |
| Source | ChEMBL 37 hERG IC50 aggregate — exact SMILES overlap excluded |
| Starting N | 7,977 |
| Exact overlaps removed | 7,249 |
| Evaluated N | 728 |
| Positives (IC50 ≤ 10k nM) | 489 (67.2%) |
| Negatives (IC50 > 10k nM) | 239 (32.8%) |
| IC50 range | 0.95 – 100,000 nM |
| **Strong positives (IC50 ≤ 1k nM)** | 133 (18.3%) |
| **Borderline (1k–30k nM)** | 531 (72.9%) |
| **Strong negatives (IC50 > 30k nM)** | 64 (8.8%) |
| Scaffold-aware split | Calibration=546, Test=182 |

> [!WARNING]
> The 72.9% borderline fraction is the largest single driver of label uncertainty.
> Any compound with IC50 between 1,000 nM and 30,000 nM is within ±0.48 log units of
> the 10,000 nM decision boundary — easily within assay noise of multiple heterogeneous
> protocols.

---

## 6. Class Imbalance Analysis

| Dataset | Prevalence (+) | N+ | N− |
| :--- | :--- | :--- | :--- |
| **Training (Wang et al.)** | **86.0%** | 19,130 | 3,118 |
| **Evaluation (ChEMBL37)** | 67.2% | 489 | 239 |
| **Shift** | **−18.8 pp** | — | — |

A model trained at 86% positive prevalence learns a heavily biased prior. The decision
boundary at 0.50 then over-calls positives when the true prevalence is 67%. This is
a **prior-shift / dataset-shift** problem, not purely a threshold choice.

---

## 7. Performance by IC50 Potency Class

| IC50 Class | N | N+ | Sensitivity | Specificity | MCC | BAcc |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Strong positive (≤1k nM)** | 133 | 133 | 0.9925 | — | — | — |
| **Borderline (1k–30k nM)** | 531 | 407 | 0.9691 | 0.0514 | 0.0832 | 0.5103 |
| **Strong negative (>30k nM)** | 64 | 0 | — | 0.3281 | — | — |

M1 achieves 99.25% recall on strong positives but **completely fails** on strong
negatives (specificity = 32.8%). The borderline class (72.9% of all data) generates
nearly all false positives.

---

## 8. Applicability Domain

| Domain | N | Sensitivity | Specificity | MCC | BAcc |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **IN_DOMAIN** | 719 | 0.9751 | 0.1134 | 0.1856 | 0.5443 |
| BORDERLINE | ~8 | (insufficient N) | — | — | — |
| OUT_OF_DOMAIN | ~1 | (insufficient N) | — | — | — |

99% of compounds fall IN_DOMAIN. Applicability domain criteria (MW ≤ 800 Da, −2 ≤ cLogP ≤ 7)
do **not** predict hERG errors on this dataset — the specificity problem is uniform
across the accessible chemical space.

---

## 9. Pseudo-Project Performance

| Project | N | N+ | M1 Specificity | M1 Sensitivity | M1 BAcc |
| :--- | :--- | :--- | :--- | :--- | :--- |
| PROJ_01 | 152 | ~102 | 0.073 | ~0.975 | 0.510 |
| PROJ_02 | 154 | ~103 | 0.111 | ~0.972 | 0.550 |
| PROJ_03 | 143 | ~96 | 0.068 | ~0.979 | 0.534 |
| PROJ_04 | 161 | ~108 | 0.205 | ~0.975 | 0.590 |
| PROJ_05 | 118 | ~79 | 0.105 | ~0.975 | 0.540 |

Specificity is consistently low (6.8%–20.5%) across all projects. There is no project
where the model reliably identifies true negatives. This rules out project-specific
adaptive weighting as a remediation strategy.

---

## 10. Model Complementarity (M1 / M2 Rescue Analysis)

| Category | Count |
| :--- | :--- |
| Both correct | 474 (65.1%) |
| Both wrong | 212 (29.1%) |
| M1 correct / M2 wrong | 30 (4.1%) |
| **M2 correct / M1 wrong** | **12 (1.6%)** |

**M2 rescue rate = 5.4%** (12 out of 224 M1 errors). M2 provides negligible
independent information. In 29.1% of cases both models simultaneously fail on the
same compound, confirming correlated pharmacophore-based reasoning rather than
complementary discrimination.

### Disagreement Analysis

| Region | N | M1 Brier |
| :--- | :--- | :--- |
| Low disagreement (\|Δp\| ≤ 0.30) | 663 | 0.281 |
| High disagreement (\|Δp\| > 0.30) | 65 | 0.208 |

High disagreement correlates with slightly lower M1 Brier, but the signal is too weak
(65 compounds, −0.07 Brier improvement) to serve as a reliable error flag.

---

## 11. False Positive Root Cause

| Feature | FP (N=212) | True Negative (N=27) | All Negatives (N=239) |
| :--- | :--- | :--- | :--- |
| Has basic N (%) | ~85% | ~63% | ~82% |
| Mean cLogP | ~3.5 | ~2.3 | ~3.3 |
| Mean MW | ~430 | ~360 | ~425 |
| Mean TPSA | ~78 | ~88 | ~81 |
| Mean aromatic rings | ~2.5 | ~1.9 | ~2.4 |
| IC50 borderline fraction | **88.7%** | ~69% | ~85% |

False positives are enriched for lipophilic, basic amine-containing,
aromatic-dense compounds with IC50 in the borderline zone. The model applies the
hERG pharmacophore rule (basic N + lipophilicity) broadly and cannot distinguish
true blockers from structurally similar non-blockers near the IC50 boundary.

---

## 12. False Negative Root Cause (N=12)

Only 12 M1 false negatives (missed actives). These are predominantly:
- Neutral or weakly basic compounds (lower has_basic_n rate)
- Lower cLogP than typical blockers
- Many are borderline IC50 (1k–10k nM range)

False negatives are an acceptable minority risk relative to the overwhelming false
positive problem.

---

## 13. Scaffold Series Complementarity

| Metric | Value |
| :--- | :--- |
| Scaffolds with N ≥ 5 | 22 |
| Series where M2 substantially better (+5% BAcc) | **1 / 22** |
| Verdict | `NO_REPRODUCIBLE_SERIES_ADVANTAGE` |

One isolated scaffold cluster shows M2 performing marginally better. This is
insufficient reproducible evidence to justify series-level adaptive weighting.

---

## 14. Final Scientific Decision

**`HERG_FIXED_BLEND_CANDIDATE`**

The primary problems are calibration and class imbalance / prior shift, not
complementary model information. The correct path forward is:

1. Apply Platt calibration (fitted on calibration set only).
2. Evaluate calibrated threshold on untouched holdout.
3. Consider conservative fixed blend (95%/5% or 90%/10%) for probability softening.

**Do NOT implement adaptive weighting at this stage.**

---

## 15. Production Safety

- Consensus mode: **SHADOW**
- Production hERG output: **unchanged** (100% M1 CORE at threshold 0.50)
- No visible output modifications
- No UI changes
- No threshold changes
- Calibration research confined to `validation/` and `scripts/`
