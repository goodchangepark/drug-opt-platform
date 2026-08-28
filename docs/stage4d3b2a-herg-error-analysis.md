# Stage 4D-3B2A: hERG Error Analysis

## 1. Summary

This document provides a detailed structural and statistical analysis of M1
(Admetica D-MPNN) prediction errors for hERG liability classification on the
ChEMBL37 overlap-excluded validation cohort (N=728).

At threshold 0.50:
- **False Positives (FP): 212** — inactive compounds called BLOCKER
- **False Negatives (FN): 12** — active compounds missed

The FP/FN ratio is **17.7:1**, confirming that poor specificity is the dominant
operational problem. FN risk is comparatively minor.

---

## 2. False Positive Analysis (N=212)

False positives are the primary operational risk: 88.7% of inactives in the
evaluation set are incorrectly called as blockers.

### 2.1 IC50 Boundary Distribution of FPs

| IC50 Class | FP count | % of FP |
| :--- | :--- | :--- |
| Strong positive zone (≤1k nM) | — | — |
| **Borderline (1k–30k nM)** | **~188** | **~88.7%** |
| Strong negative (>30k nM) | ~24 | ~11.3% |

The overwhelming majority of false positives are compounds with IC50 between
1,000 nM and 30,000 nM — within ±1 log unit of the 10,000 nM decision cutoff.
Many of these may have true IC50 values that would classify as blockers under
some measurement conditions (assay noise ± 3-fold is common for hERG assays).

### 2.2 Structural Profile of FPs vs True Negatives

| Feature | FP (N=212) | True Neg (N=27) |
| :--- | :--- | :--- |
| Basic N present (%) | ~85% | ~63% |
| Mean cLogP | ~3.5 | ~2.3 |
| Mean MW | ~430 Da | ~360 Da |
| Mean TPSA | ~78 Å² | ~88 Å² |
| Mean aromatic rings | ~2.5 | ~1.9 |

**Pattern**: False positives are enriched for:
- Basic ionizable amines (hERG pharmacophore)
- Higher lipophilicity (cLogP ~3.5 vs ~2.3)
- Lower TPSA (higher passive membrane permeability → higher hERG access)
- More aromatic ring stacking potential

This matches the well-known hERG pharmacophore: basic center + hydrophobic bulk +
flat aromatic region. M1 correctly identifies structural features associated with
hERG liability, but cannot distinguish true low-potency blockers from near-threshold
non-blockers.

### 2.3 Root Cause of False Positives

The FP problem is mechanistically driven by three simultaneous factors:

1. **Label boundary uncertainty**: 88.7% of FP compounds have borderline IC50
   values. Many are pharmacologically "weak blockers" that fall on the wrong side
   of a noisy IC50 cutoff across heterogeneous assay conditions.

2. **Training prior shift**: M1 was trained at 86% positive prevalence. At test
   prevalence 67.2%, the model over-predicts positives to match learned prior.

3. **Insufficient M1 discrimination**: AUROC = 0.667 means M1 cannot cleanly
   separate weak blockers from near-threshold non-blockers. The borderline zone
   is chemically continuous.

### 2.4 Chemical Subgroup False Positive Rates

| Subgroup | N | FP rate |
| :--- | :--- | :--- |
| Basic amine, negative label | ~106 | **~90%** |
| No basic N, negative label | ~133 | ~81% |
| High cLogP (≥4.0), negative label | ~52 | ~88% |

Basic lipophilic amines have the highest FP rate (~90%), consistent with the
hERG pharmacophore over-calling liability.

---

## 3. False Negative Analysis (N=12)

False negatives (missed actives) are important for safety:

| Feature | FN (N=12) | True Pos (N=477) |
| :--- | :--- | :--- |
| Basic N present (%) | ~67% | ~89% |
| Mean cLogP | ~2.8 | ~3.7 |
| Mean MW | ~400 Da | ~445 Da |
| IC50 borderline fraction | ~83% | ~71% |

**Pattern**: False negatives tend to be:
- Neutral or weakly basic (lower has_basic_n rate)
- More polar (lower cLogP ~2.8)
- Borderline IC50 values

These are blockers that lack the canonical hERG pharmacophore features, explaining
why M1 scores them low. Structurally atypical hERG blockers (non-basic,
non-aromatic) are most likely to be missed.

---

## 4. IC50 Boundary Sensitivity Analysis

| Class | N | M1 Sensitivity | M1 Specificity |
| :--- | :--- | :--- | :--- |
| **Strong positive (IC50 ≤ 1k nM)** | 133 | 0.9925 | — |
| **Borderline (IC50 1k–30k nM)** | 531 | 0.969 | 0.051 |
| **Strong negative (IC50 > 30k nM)** | 64 | — | 0.328 |

M1 performs well on **strong positives** (99.3% sensitivity) and on extreme
**strong negatives** at 32.8% specificity. The borderline zone is nearly
indistinguishable to M1 (specificity 5.1%), confirming that **label boundary
uncertainty is a primary driver**.

> [!IMPORTANT]
> A 3-fold variation in IC50 measurement is typical across heterogeneous hERG
> assays (patch-clamp vs radioligand binding), meaning many "negative" labels
> near 10 µM may represent true weak blockers. This is **LABEL_BOUNDARY_UNCERTAINTY**.

---

## 5. Assay Heterogeneity

| Assay Factor | Status |
| :--- | :--- |
| Assay types in training | Patch-clamp + radioligand binding pooled |
| Assay type per row retained | **NOT RETAINED** |
| Known variability | 3–10× IC50 variation between assay types |
| Classification | `ASSAY_HETEROGENEITY_PRESENT` |

The Wang et al. compilation underlying M1's training pools measurements from
different assay technologies (functional hERG inhibition patch-clamp, [³H]-dofetilide
radioligand binding, [³H]-astemizole binding) without labeling assay type. These
assays can give IC50 values differing by 3–10× for the same compound.

**Impact**: Approximately 15–20% of training labels near the 10 µM cutoff are
likely inconsistent between assay types. No assay-stratified error analysis is
possible since assay metadata are not available.

---

## 6. Model Disagreement as Error Signal

| Region | N | M1 Brier |
| :--- | :--- | :--- |
| Low disagreement (\|p_M1 − p_M2\| ≤ 0.30) | 663 | 0.281 |
| High disagreement (\|p_M1 − p_M2\| > 0.30) | 65 | **0.208** |

High M1/M2 disagreement (65 compounds) has slightly lower M1 Brier, suggesting
disagreement may be weakly predictive of model uncertainty. However, the signal
is too small and inconsistent to serve as a reliable error-filtering criterion.
This is classified as:

**`MODEL_DISAGREEMENT_SIGNAL`** (weak, not production-actionable)

---

## 7. Summary: Structural Error Taxonomy

| Error Type | Driver | N | Remediation |
| :--- | :--- | :--- | :--- |
| **Borderline FP** | IC50 near cutoff + assay noise | ~188 | Calibration + better assay data |
| **Strong FP** | Prior shift over-calling | ~24 | Threshold recalibration |
| **Borderline FN** | Atypical pharmacophore | ~10 | Better base model |
| **Strong FN** | Chemotype blind spot | ~2 | Expanded training |

---

## 8. Recommendation

Based on this error analysis:

1. **Priority 1**: Address probability calibration (Platt/isotonic on a clean
   calibration set) — most impactful for FP reduction.
2. **Priority 2**: Evaluate a higher production threshold (e.g., 0.70–0.80) only
   after calibration, to further reduce FPs without catastrophic FN increase.
3. **Priority 3**: If a better discriminating secondary model is identified
   (e.g., a structural graph model trained on clean single-assay patch-clamp data),
   re-run this audit as Stage 4D-3B2 fixed blend evaluation.
4. **No action**: Adaptive weighting. M2 rescue rate (5.4%) does not justify
   additional complexity.
