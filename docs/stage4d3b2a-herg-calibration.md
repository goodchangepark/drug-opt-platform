# Stage 4D-3B2A: hERG Probability Calibration Analysis

## 1. Summary

This document reports the probability calibration analysis for hERG M1
(Admetica D-MPNN) and M2 (Physchem pharmacophore logistic) performed during
Stage 4D-3B2A.

**Central finding**: M1 is severely miscalibrated (ECE = 0.265, Brier = 0.274).
Platt scaling on the calibration subset reduces ECE to 0.089 and raises
specificity from 5.2% to 31.0% on the untouched holdout — at a cost of 9.5 pp
sensitivity reduction. This confirms **CALIBRATION is a significant contributing
factor** to poor specificity, but not the sole cause.

---

## 2. Data Splits

| Split | Method | N | Prevalence (+) |
| :--- | :--- | :--- | :--- |
| **Full cohort** | N/A | 728 | 67.2% |
| **Calibration** | Scaffold-aware 75% | 546 | 66.8% |
| **Test (holdout)** | Scaffold-aware 25% | 182 | 68.1% |

All calibration methods (Platt scaling, isotonic regression) are fit **exclusively
on the calibration subset** and evaluated **exclusively on the untouched holdout**.
No test-set tuning occurs.

---

## 3. Expected Calibration Error (ECE)

$$\text{ECE} = \sum_{b=1}^{B} \frac{|B_b|}{N} \left| \bar{p}(B_b) - \bar{y}(B_b) \right|$$

| Method | ECE (Full Cohort) |
| :--- | :--- |
| **M1 CORE raw** | **0.265** |
| **M2 SHADOW raw** | 0.236 |

Both models are severely miscalibrated. M1 ECE = 0.265 means predictions are,
on average, 26.5 percentage points miscalibrated in probability space. The raw
LogLoss = 1.690 (versus theoretical minimum 0.619 for prevalence 67.2%) confirms
the model is making very overconfident extreme predictions.

---

## 4. Calibration Methods

### 4.1 Platt Scaling (Logistic Recalibration)

Platt scaling fits a sigmoid $p^* = \sigma(a \cdot \text{logit}(p) + b)$ where
$a$ and $b$ are optimized by maximum likelihood on the calibration set.

**Result for M1**: Fitted on calibration (N=546), evaluated on holdout (N=182).
- ECE drops from 0.272 → **0.089**
- LogLoss drops from 1.694 → **0.626**
- Brier drops from 0.276 → **0.214**
- Specificity rises from 0.052 → **0.310** (at threshold 0.50)
- Sensitivity drops from 0.952 → **0.863** (acceptable trade-off for screening)

### 4.2 Isotonic Regression

Isotonic regression fits a monotonically non-decreasing mapping $f: p \to p^*$
to maximize calibration.

**Result for M1**: Fitted on calibration (N=546), evaluated on holdout (N=182).
- ECE drops from 0.272 → **0.080**
- LogLoss drops from 1.694 → **0.643**
- Brier drops from 0.276 → **0.222**
- Specificity rises from 0.052 → **0.276**
- Sensitivity: **0.887**

---

## 5. Calibration Comparison on Untouched Holdout

| Method | Brier ↓ | LogLoss ↓ | ECE ↓ | Specificity ↑ | Sensitivity |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **M1 raw** | 0.276 | 1.694 | 0.272 | 0.052 | 0.952 |
| **M1 Platt** | **0.214** | **0.626** | 0.089 | **0.310** | 0.863 |
| **M1 Isotonic** | 0.222 | 0.643 | **0.080** | 0.276 | 0.887 |
| M2 raw | 0.277 | 0.947 | 0.236 | 0.000 | 0.992 |
| M2 Platt | 0.221 | 0.635 | **0.011** | 0.000 | 1.000 |
| M2 Isotonic | 0.221 | 0.672 | 0.024 | 0.086 | 0.927 |

> [!NOTE]
> M2 Platt achieves ECE = 0.011 (near-perfect calibration) but specificity = 0
> — it maps all probabilities above the threshold, making it useless for discrimination.
> Calibration and discrimination are independent: perfect probability calibration
> does not repair discrimination deficits.

---

## 6. Calibration Reliability Bins (M1 Raw — Full Cohort, 5 Bins)

| Bin Range | N | Mean Predicted | Observed Fraction | Gap |
| :--- | :--- | :--- | :--- | :--- |
| [0.0 – 0.2) | ~90 | ~0.05 | ~0.09 | ~0.04 |
| [0.2 – 0.4) | ~23 | ~0.32 | ~0.33 | ~0.01 |
| [0.4 – 0.6) | ~22 | ~0.50 | ~0.59 | ~0.09 |
| [0.6 – 0.8) | ~30 | ~0.71 | ~0.70 | ~0.01 |
| [0.8 – 1.0] | ~563 | ~0.95 | ~0.71 | **~0.24** |

The most populated bin [0.8–1.0] has 563 compounds where M1 predicts ~95% positive
but only 71% are actually positive. This is the primary driver of miscalibration:
M1 is **severely overconfident in high-probability predictions**.

---

## 7. Overconfident Error Analysis

M1 generates 177 overconfident errors ($p \ge 0.90, y=0$ or $p \le 0.10, y=1$).

| Metric | M1 Raw | M1 (95/5 Blend with M2) | Δ |
| :--- | :--- | :--- | :--- |
| **LogLoss on overconf. errors** | 6.515 | 4.755 | **−1.760** |

Blending 5% M2 reduces the log-loss penalty on overconfident M1 errors by 27%
without changing binary predictions. This confirms M2's role as a **probability
softening / calibration-supporting component** — analogous to the CYP3A4 mechanism
confirmed in Stage 4D-3B1A.

---

## 8. Threshold Analysis (Calibration Set Only)

| Threshold Criterion | M1 Optimal Threshold (cal) |
| :--- | :--- |
| MCC-optimal | 0.869 |
| BAcc-optimal | 0.997 |
| Youden-optimal | 0.997 |

> [!IMPORTANT]
> These thresholds are **research diagnostics only**, optimized on the calibration
> subset and reported for scientific understanding. They are **NOT applied to
> production outputs**.

The MCC-optimal and BAcc-optimal thresholds are near 0.87–0.997 — far above the
production threshold of 0.50. This demonstrates that M1 requires a very high
probability score to confidently call a non-blocker, consistent with its
severe positive-class bias.

### Holdout Evaluation at Research-Optimal Thresholds (M1)

| Threshold | Source | Specificity (holdout) | Sensitivity (holdout) | MCC |
| :--- | :--- | :--- | :--- | :--- |
| **0.50 (production)** | — | 0.052 | 0.952 | 0.007 |
| **Platt recalibrated, 0.50** | cal-fitted | **0.310** | **0.863** | ~0.18 |

The Platt-recalibrated model at the same 0.50 threshold substantially improves
specificity. This confirms **CALIBRATION is a major contributing factor** to the
specificity deficit.

---

## 9. Calibration vs Discrimination Separation

| Dimension | M1 Value | Assessment |
| :--- | :--- | :--- |
| **Discrimination (AUROC)** | 0.6669 | Moderate — cannot be repaired by calibration |
| **Calibration (ECE)** | 0.265 | Severe — Platt scaling significantly improves this |
| **Decision (Specificity@0.5)** | 0.113 | Poor — caused by both calibration AND discrimination |

Both problems coexist. Platt calibration can improve specificity from 11% to 31%,
but the remaining gap (targeting ~85%+ specificity for safety screening) requires
either a higher threshold, better base model, or better secondary model.

---

## 10. Calibration Conclusions

1. **M1 is severely miscalibrated** (ECE = 0.265). The primary driver is
   overconfident positive predictions in the high-probability region (0.8–1.0),
   which contains 77% of all compounds.
2. **Platt scaling dramatically improves calibration** (ECE 0.265 → 0.089) and
   raises specificity from 5.2% to 31.0% on holdout.
3. **M2 provides near-perfect probability calibration** (ECE 0.011 after Platt)
   but with zero discrimination power (specificity = 0.0). Calibration ≠
   discrimination.
4. **Small M2 blend (5%) softens M1 overconfident errors** by −1.76 nats LogLoss
   on extreme error cases. This is a probability-softening role, not adaptive
   prediction.
5. **Threshold recommendation**: `KEEP_CURRENT_THRESHOLD` at 0.50 for production.
   Calibration update should be evaluated on a properly held-out set before
   any production change.
