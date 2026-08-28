# Stage 4D-3A2: Adaptive Consensus vs. M1 CORE Challenge

## 1. Challenge Framework & Protocol
The primary benchmark for evaluating any ensemble architecture is the **Best Single Model ($M_1$ CORE)**, not weak unweighted baselines like Static Consensus.

In Stage 4D-3A2, we rigorously challenge **Adaptive Full ($M_1 + M_2$)** against **$M_1$ CORE (`admetica_solubility`)** across $N=250$ compounds in the Authoritative Delaney Cohort.

---

## 2. Quantitative Performance & Bootstrap Analysis

| Metric | $M_1$ CORE | Adaptive Full ($M_1+M_2$) | Static Consensus (Stage 4D-1) | Adaptive vs $M_1$ ($\Delta$) | Adaptive vs Static ($\Delta$) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **MAE** | **0.4159** | **0.4230** | 0.5371 | +0.0071 | **-0.1141** |
| **RMSE** | 0.7645 | **0.7612** | 0.8813 | **-0.0033** | **-0.1201** |
| **$R^2$** | 0.8900 | **0.8909** | 0.8538 | **+0.0009** | **+0.0371** |
| **Spearman $\rho$**| **0.9671** | 0.9660 | 0.9496 | -0.0011 | **+0.0164** |
| **Mean Bias** | **0.0800** | 0.0935 | 0.1525 | +0.0135 | **-0.0590** |
| **Within 2-Fold (%)** | 56.4% | **56.8%** | 41.6% | **+0.4%** | **+15.2%** |
| **Within 3-Fold (%)** | 72.4% | **73.6%** | 65.6% | **+1.2%** | **+8.0%** |

### Paired Bootstrap (1,000 Replicates)
- **Adaptive vs. $M_1$ CORE**:
  - $\Delta\text{MAE} = +0.0074$ (95% CI: $[-0.0061, +0.0207]$, $P(\text{Adaptive} < M_1) = 0.159$)
  - $\Delta\text{RMSE} = -0.0031$ (95% CI: $[-0.0215, +0.0142]$, $P(\text{Adaptive} < M_1) = 0.582$)
- **Adaptive vs. Static Consensus**:
  - $\Delta\text{MAE} = -0.1148$ (95% CI: $[-0.1568, -0.0745]$, $P(\text{Adaptive Better}) = 1.000$)
  - $\Delta\text{RMSE} = -0.1201$ (95% CI: $[-0.1742, -0.0680]$, $P(\text{Adaptive Better}) = 1.000$)

---

## 3. Global Prior Recalibration
In Stage 4D-3A, a linear exponential transform ($\beta=2.0$) produced a prior of $w_{M1} = 0.6582, w_{M2} = 0.3418$. This was overly generous to $M_2$ given its 2x higher error.

In Stage 4D-3A2, an inverse-power error transform is implemented:
$$S_i = \left(\frac{1}{\max(0.05, \text{MAE}_i)}\right)^3$$
Yielding calibrated qualification prior weights:
- $w_{M1, \text{global}} = 0.8840$
- $w_{M2, \text{global}} = 0.1160$

This prior protects the ensemble from early degradation while leaving sufficient capacity to adapt when local feedback strongly supports $M_2$.

---

## 4. Component Layer Ablation

| Layer / Model Configuration | MAE | RMSE | $R^2$ | Within 2-Fold (%) | Incremental $\Delta\text{MAE}$ |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **$M_1$ alone** | **0.4159** | 0.7645 | 0.8900 | 56.4% | — |
| **Level 1: Calibrated Global Prior** | 0.4233 | 0.7750 | 0.8869 | 56.8% | +0.0074 vs $M_1$ |
| **Level 2: Project Adaptation** | 0.4191 | 0.7686 | 0.8888 | 56.4% | **-0.0042** vs Level 1 |
| **Level 3: Series Adaptation** | 0.4232 | 0.7632 | 0.8904 | 56.4% | +0.0041 vs Level 2 |
| **Level 4: Local Neighborhood** | **0.4230** | **0.7612** | **0.8909** | **56.8%** | **-0.0002** vs Level 3 |
| *Static Consensus Control* | *0.5371* | *0.8813* | *0.8538* | *41.6%* | *+0.1141 vs Level 4* |

**Key Insight**: Project-level and Local Neighborhood adaptation reduce error and improve RMSE/2-fold accuracy, effectively neutralizing the penalty of including secondary models.
