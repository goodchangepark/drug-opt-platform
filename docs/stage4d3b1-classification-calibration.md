# Stage 4D-3B1: Probability Calibration & Classification Performance Audit

## 1. Overview
Evaluating binary classification models for ADMET ensembling requires separating **probability calibration** (Brier score, Bounded Log Loss, reliability) from **decision thresholding** (MCC, Balanced Accuracy, Sensitivity, Specificity).

This document presents the calibration analysis of qualified CYP3A4 models on the authoritative $N=250$ cohort.

---

## 2. Performance Comparison Table ($N=250$)

| Method | Role / Status | MCC | Balanced Accuracy | Brier Score | Bounded Log Loss | AUROC | AUPRC | Sensitivity | Specificity |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **$M_1$ (Admetica D-MPNN)** | **CORE** | **0.8334** | **0.9166** | **0.0726** | 0.2646 | **0.9572** | **0.9217** | **0.9426** | **0.8906** |
| **$M_2$ (Morgan GBDT)** | SHADOW_ONLY | 0.4365 | 0.7043 | 0.2056 | 0.6076 | 0.7627 | 0.6983 | 0.8852 | 0.5234 |
| **Static Consensus** | SHADOW (Stage 4D-1) | 0.7974 | 0.8973 | 0.1046 | 0.3583 | 0.9458 | 0.9012 | 0.9508 | 0.8438 |
| **Adaptive Full Ensembling** | **SHADOW (Stage 4D-3B1)** | **0.8334** | **0.9166** | **0.0728** | **0.2465** | **0.9572** | **0.9217** | **0.9426** | **0.8906** |
| **Base Rate (Dataset Prevalence)** | Trivial Control | 0.0000 | 0.5000 | 0.2499 | 0.6929 | 0.5000 | 0.4880 | 0.0000 | 1.0000 |

---

## 3. Probability Calibration Findings

### 3.1. Brier Score & Bounded Log Loss
- **$M_1$ vs $M_2$ Quality Gap**:
  - $M_1$ exhibits high calibration fidelity ($\text{Brier} = 0.0726$).
  - $M_2$ suffers from poor probability separation ($\text{Brier} = 0.2056, \text{LogLoss} = 0.6076$), producing overly confident false positives in inactive compounds.
- **Adaptive Recovery**:
  - Unweighted static consensus inflates Brier loss to $0.1046$ (+44% error increase over $M_1$).
  - Adaptive ensembling preserves $M_1$'s sharp calibration ($\text{Brier} = 0.0728$) while slightly reducing Log Loss ($\text{LogLoss} = 0.2465$ vs $0.2646$ for $M_1$) by smoothing boundary probabilities.

### 3.2. Class Balance Safeguards
When a medicinal chemistry project or series contains only positive or only negative experimental records ($N \ge 5, \text{unique}(y) = 1$):
- Threshold-dependent metrics (such as MCC or Balanced Accuracy) become mathematically undefined or ill-conditioned.
- The adaptive engine detects this condition, flags `CLASS_BALANCE_LIMITED`, and restricts learning exclusively to well-bounded Brier probability loss with conservative Bayesian shrinkage toward the global prior.

---

## 4. Model Disagreement Signal Analysis ($|p_{M1} - p_{M2}|$)

To determine whether the inter-model disagreement magnitude $\Delta = |p_{M1} - p_{M2}|$ serves as a reliable surrogate for prediction uncertainty:

| Disagreement Bin | Range | Sample Count ($N$) | $M_1$ Brier Error Rate | Adaptive Brier Error Rate |
| :--- | :--- | :--- | :--- | :--- |
| **Low Disagreement** | $\Delta < 0.20$ | 165 | **0.0271** | **0.0272** |
| **Moderate Disagreement** | $0.20 \le \Delta < 0.40$ | 44 | **0.0734** | **0.0740** |
| **High Disagreement** | $\Delta \ge 0.40$ | 41 | **0.2547** | **0.2553** |

### Key Takeaway:
Model disagreement $\Delta$ exhibits strong positive monotonicity with true observation error rate ($0.0271 \to 0.0734 \to 0.2547$). High disagreement ($\Delta \ge 0.40$) flags compounds with ~10x higher error probability, confirming its utility as an objective `MODEL_DISAGREEMENT_SIGNAL` metadata flag for drug discovery teams.
