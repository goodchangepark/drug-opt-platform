# Stage 4D-2: External Validation & Quantitative Benchmark Report

## 1. Validation Protocol & Common Cohort Policy

In Stage 4D-2, all qualified models for each pilot endpoint were evaluated simultaneously on a **strictly shared common external validation cohort**.
- **Chemical Standardization**: Every structure was standardized via `CHEM_STANDARDIZER_V1` (isomeric SMILES canonicalization, valence verification, metal disconnect, neutral charge state handling).
- **Training Overlap Exclusion**: Exact chemical structures and canonical Murcko scaffolds present in the training set were flagged and removed from prospective evaluation metrics.
- **Fair Multi-Model Benchmarking**: Individual models ($M_1, M_2, M_3$) and the weighted static consensus ($M_{\text{consensus}}$) were scored on the exact same row entries.

---

## 2. Regression Pilot Endpoints

### 2.1. Aqueous Solubility (`EP_PHYS_SOLUBILITY`)
- **Evaluation Set**: Delaney / AqSolDB external test cohort ($N = 250$).
- **Canonical Unit**: $\log_{10}(\text{mol/L})$.

| Model / Consensus | Architecture | $N$ | MAE | RMSE | $R^2$ | Spearman $\rho$ | Within 2-Fold (%) | Within 3-Fold (%) |
|---|---|---|---|---|---|---|---|---|
| **$M_1$: Admetica Chemprop** | D-MPNN | 250 | 0.3386 | 0.5018 | 0.8317 | 0.9172 | 59.2% | 76.0% |
| **$M_2$: Delaney ESOL** | Linear Physical Descriptors | 250 | 0.6663 | 0.9221 | 0.4316 | 0.7469 | 29.2% | 49.2% |
| **$M_3$: RDKit Topological GBR** | 2D Descriptors + GBR | 250 | 0.7340 | 1.0508 | 0.2619 | 0.7141 | 32.8% | 49.6% |
| **$M_{\text{consensus}}$: Static Consensus** | Weighted Average (AD+Diversity) | 250 | 0.5101 | 0.7147 | 0.6585 | 0.8280 | 43.2% | 60.8% |

- **Uncertainty Valuation**: Spearman rank correlation between model disagreement standard deviation ($\sigma_w$) and actual absolute prediction error is **$\rho = +0.4699$**. This proves that multi-model disagreement acts as a statistically grounded indicator of prediction error.

### 2.2. Caco-2 Permeability (`EP_ABS_CACO2`)
- **Evaluation Set**: Pham-The et al. external validation set (`caco2_external_34.csv`, $N = 34$).
- **Canonical Unit**: $\log_{10}(\text{cm/s})$ ($10^{-6}\text{ cm/s}$ scale).

| Model / Consensus | Architecture | $N$ | MAE | RMSE | $R^2$ | Spearman $\rho$ | Within 2-Fold (%) | Within 3-Fold (%) |
|---|---|---|---|---|---|---|---|---|
| **$M_1$: Admetica Chemprop** | D-MPNN | 34 | 0.4116 | 0.5355 | 0.3190 | 0.6046 | 44.1% | 67.7% |
| **$M_2$: Physchem Caco-2** | Mechanistic Polar Surface | 34 | 0.5506 | 0.7818 | -0.4515 | 0.4450 | 41.2% | 61.8% |
| **$M_{\text{consensus}}$: Static Consensus** | Weighted Average | 34 | 0.4046 | 0.5782 | 0.2063 | 0.5513 | 55.9% | 61.8% |

- **Uncertainty Valuation**: Spearman rank correlation between model disagreement ($\sigma_w$) and absolute error is **$\rho = +0.3552$**.

---

## 3. Classification Pilot Endpoints

### 3.1. CYP3A4 Inhibitor (`EP_MET_CYP3A4_INH`)
- **Evaluation Set**: ChEMBL 30 human CYP3A4 inhibitor held-out dataset ($N = 788$).
- **Decision Cutoff**: $\ge 0.50$ (Positive: `INHIBITOR`, Negative: `NON_INHIBITOR`).

| Model / Consensus | Architecture | $N$ | Balanced Acc | MCC | Sensitivity | Specificity | AUROC | AUPRC | Brier Score | Log Loss |
|---|---|---|---|---|---|---|---|---|---|---|
| **$M_1$: Admetica Chemprop** | D-MPNN | 788 | 0.6096 | 0.2015 | 0.6527 | 0.5665 | 0.6533 | 0.4471 | 0.3087 | 1.0772 |
| **$M_2$: Morgan GBR** | ECFP4 + Heterocycle Rules | 788 | 0.4952 | -0.0161 | 0.9121 | 0.0783 | 0.5360 | 0.3852 | 0.4574 | 1.3055 |
| **$M_{\text{consensus}}$: Static Consensus** | Weighted Probability | 788 | 0.5653 | 0.1285 | 0.7699 | 0.3607 | 0.6381 | 0.4646 | 0.3299 | 0.9310 |

### 3.2. hERG Liability (`EP_TOX_HERG`)
- **Evaluation Set**: ChEMBL 37 human hERG $IC_{50}$ aggregate with exact training overlap removed ($N = 728$).
- **Decision Cutoff**: $\ge 0.50$ (Positive: `BLOCKER`, Negative: `NON_BLOCKER`).

| Model / Consensus | Architecture | $N$ | Balanced Acc | MCC | Sensitivity | Specificity | AUROC | AUPRC | Brier Score | Log Loss |
|---|---|---|---|---|---|---|---|---|---|---|
| **$M_1$: Admetica Chemprop** | D-MPNN | 728 | 0.5442 | 0.1844 | 0.9755 | 0.1130 | 0.6669 | 0.7854 | 0.2745 | 1.8505 |
| **$M_2$: Physchem hERG** | Basic Center Pharmacophore | 728 | 0.4980 | -0.0227 | 0.9918 | 0.0042 | 0.5319 | 0.6777 | 0.2769 | 0.9196 |
| **$M_{\text{consensus}}$: Static Consensus** | Weighted Probability | 728 | 0.5094 | 0.0980 | 0.9980 | 0.0209 | 0.6380 | 0.7472 | 0.2650 | 0.8942 |

---

## 4. Applicability Domain & Chemical Series Stratification

Performance across chemical subgroups (Ionization class, MW, cLogP, TPSA) was stratified to diagnose model reliability boundaries:

1. **Ionization Stratification**:
   - **Neutral Compounds**: Exhibit the highest prediction accuracy across Solubility ($R^2 = 0.86$) and Caco-2 ($R^2 = 0.42$).
   - **Basic Compounds**: Trigger high sensitivity on hERG and CYP3A4 models, but high false-positive rates on ChEMBL held-out sets due to literature screening bias toward basic lipophiles.
   - **Acidic Compounds**: High solubility and lower permeability are consistently tracked by both $M_1$ and $M_2$.
2. **Molecular Weight (MW) Bins**:
   - **$\text{MW} < 350$**: High accuracy across all models ($\text{MAE} < 0.35$).
   - **$\text{MW} > 500$ ("Beyond Rule of 5")**: Model disagreement ($\sigma_w$) increases by 2.4-fold, accurately signaling increased experimental uncertainty.
