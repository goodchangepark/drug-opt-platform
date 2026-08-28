# Stage 4D-2: Model Error Correlation & Diversity Analysis

## 1. Mathematical Rationale for Error Diversity

The foundational requirement of any statistical ensemble is **error independence**:
$$\text{Var}(e_{\text{ensemble}}) = \sum_i w_i^2 \text{Var}(e_i) + 2 \sum_{i < j} w_i w_j \text{Cov}(e_i, e_j)$$

If two models make identical errors ($\text{Cov}(e_i, e_j) \approx \text{Var}(e)$), combining them provides zero variance reduction and merely doubles computational cost. Conversely, when errors are uncorrelated ($\text{Cov}(e_i, e_j) \approx 0$), ensemble variance drops by up to $1/K$.

In Stage 4D-1, diversity penalties used a placeholder assumption of 0.55. In Stage 4D-2, we replace this assumption with **empirical pairwise error correlation matrices** derived from held-out validation cohorts.

---

## 2. Empirical Error Correlation Matrices

### 2.1. Aqueous Solubility Residuals ($e_i = y_i - \hat{y}_i$)
Evaluated on $N=250$ held-out Delaney evaluation compounds:

| Model Pair | Architecture 1 vs Architecture 2 | Pearson $r(e_1, e_2)$ | Diversity Factor ($\gamma = 1 - r^2$) | Interpretation |
|---|---|---|---|---|
| **$M_1$ vs $M_2$** | D-MPNN vs Delaney Physical Model | **0.3861** | **0.85** | **High Independence**: Graph neural net and linear physical descriptors exhibit low residual collinearity. |
| **$M_1$ vs $M_3$** | D-MPNN vs RDKit Topological GBR | **0.4498** | **0.80** | **Moderate Independence**: Topological surface descriptors complement message-passing latent features. |
| **$M_2$ vs $M_3$** | Delaney Physical vs RDKit GBR | **0.8726** | **0.25** | **High Collinearity**: Both models rely heavily on 2D atom counts and cLogP; heavy diversity penalty applied. |

### 2.2. Caco-2 Permeability Residuals
Evaluated on $N=34$ Pham-The external validation set:

| Model Pair | Architecture 1 vs Architecture 2 | Pearson $r(e_1, e_2)$ | Diversity Factor | Interpretation |
|---|---|---|---|---|
| **$M_1$ vs $M_2$** | D-MPNN vs Polar Surface Model | **0.5179** | **0.73** | **Moderate Independence**: Mechanistic TPSA and charge constraints provide an independent bounding floor. |

### 2.3. CYP3A4 Inhibitor Prediction Errors ($|y_i - p_i|$)
Evaluated on $N=788$ ChEMBL 30 validation compounds:

| Model Pair | Architecture 1 vs Architecture 2 | Pearson $r(e_1, e_2)$ | Diversity Factor | Interpretation |
|---|---|---|---|---|
| **$M_1$ vs $M_2$** | D-MPNN vs Morgan Pharmacophore GBR | **0.2072** | **0.95** | **Very High Independence**: ECFP4 substructure rules and deep message passing capture complementary structural alerts. |

### 2.4. hERG Liability Prediction Errors
Evaluated on $N=728$ ChEMBL 37 validation compounds:

| Model Pair | Architecture 1 vs Architecture 2 | Pearson $r(e_1, e_2)$ | Diversity Factor | Interpretation |
|---|---|---|---|---|
| **$M_1$ vs $M_2$** | D-MPNN vs Basic Center Logistic | **0.8915** | **0.20** | **High Collinearity**: Both models strongly predict hERG liability primarily from basic nitrogen presence and lipophilicity. |

---

## 3. Empirical Model Disagreement vs Prediction Error

A critical question investigated in Stage 4D-2 is:
> *"Does multi-model disagreement reliably predict when the consensus prediction is wrong?"*

### Results:
- **Aqueous Solubility**: $\text{Spearman } \rho(\sigma_w, |y - \hat{y}|) = \mathbf{+0.4699}$ ($p < 10^{-14}$).
- **Caco-2 Permeability**: $\text{Spearman } \rho(\sigma_w, |y - \hat{y}|) = \mathbf{+0.3552}$ ($p < 0.04$).

### Conclusion:
Model disagreement standard deviation ($\sigma_w$) is a statistically verified, quantitative proxy for model uncertainty. When $\sigma_w \le 0.30$ (High Agreement), consensus predictions are highly reliable. When $\sigma_w > 0.60$ (Low Agreement), absolute prediction error is on average 2.1x higher, providing medicinal chemists with transparent warning flags.
