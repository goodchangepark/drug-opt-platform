# Stage 4D-3A: Hierarchical Bayesian Shrinkage Formulation

## 1. Theoretical Motivation
In small-molecule drug discovery, experimental laboratory data accumulates incrementally within medicinal chemistry series. 

Standard statistical learning systems face two failure modes:
1. **Under-adaptation**: Treating all series identically according to global training priors, ignoring genuine local chemotype-specific SAR trends.
2. **Over-adaptation (Overfitting / Noise Capture)**: Dramatically over-reacting to 1 or 2 noisy or outlying experimental observations, destroying predictive generalization on subsequent analogs.

The **Hierarchical Bayesian Shrinkage Engine** in Stage 4D-3A solves this tradeoff through tiered empirical Bayes shrinkage:
$$\lambda = \frac{N_{\text{effective}}}{N_{\text{effective}} + N_{\text{prior}}}$$

When $N_{\text{effective}} \ll N_{\text{prior}}$, $\lambda \to 0$ and the posterior weight collapses smoothly onto the parent prior. When $N_{\text{effective}} \gg N_{\text{prior}}$, $\lambda \to 1$ and the posterior weight reflects empirical series/local performance.

---

## 2. Mathematical Definition of Shrinkage Parameters

| Hierarchical Level | Evidence Granularity | Prior Sample Size ($N_{\text{prior}}$) | Fallback Behavior When Sparse |
| :--- | :--- | :--- | :--- |
| **Level 1: Global** | Cross-series benchmark ($N=250$) | $\infty$ (Fixed Baseline) | Stage 4D-2C Qualification Benchmark |
| **Level 2: Project** | All project experiments | $N_{\text{prior, project}} = 10.0$ | Smoothly retreats to Global Prior |
| **Level 3: Series** | Bemis-Murcko core scaffold | $N_{\text{prior, series}} = 5.0$ | Smoothly retreats to Project Posterior |
| **Level 4: Local** | Morgan $R=2, 2048\text{b}, T \ge 0.40$ | $N_{\text{prior, local}} = 3.0$ | Smoothly retreats to Series Posterior |

---

## 3. Distance-Weighted Local Neighborhood
Unlike strict categorical groupings, local chemical space is continuous. The effective local sample size $N_{\text{eff, local}}$ and local MAE are computed using continuous similarity weights:
$$T_j = \text{Tanimoto}(F_{\text{query}}, F_j), \quad \forall j \text{ s.t. } T_j \ge 0.40$$
$$N_{\text{eff, local}} = \sum_{j \in \text{neighbors}} T_j^2$$
$$\text{MAE}_{i, \text{local}} = \frac{\sum_{j} T_j \cdot |\hat{y}_{ij} - y_j|}{\sum_j T_j}$$

This ensures that a single high-similarity analog ($T=0.90$) exerts substantially greater weight than several distant analogs ($T=0.42$), while remaining mathematically bounded.

---

## 4. Minimum Floor & Applicability Domain (AD) Safeguards
Even when empirical evidence heavily favors one model in a specific chemotype, models are prevented from complete weight collapse by the minimum floor parameter:
$$\epsilon = 0.02$$

Furthermore, if a model's feature vector or molecular descriptors violate the model's Applicability Domain, its weight is penalized:
- $\gamma_{\text{AD}} = 1.0$ (`IN_DOMAIN`)
- $\gamma_{\text{AD}} = 0.5$ (`BORDERLINE`)
- $\gamma_{\text{AD}} = 0.1$ (`OUT_OF_DOMAIN`)
