# Stage 4D-2C: Leave-One-Model-Out Contribution & Ensemble Governance

## 1. Leave-One-Model-Out (LOO) Contribution Analysis

To determine the exact mathematical impact of each individual model within an ensemble, systematic leave-one-out and sub-combination evaluations were conducted.

### 1.1. Aqueous Solubility Combinations ($N=250$)

| Combination | Models Included | Aggregation Mode | MAE | RMSE | $R^2$ | Spearman $\rho$ | Contribution Classification |
|---|---|---|---|---|---|---|---|
| **$M_1$ alone** | `admetica_solubility` | Single Model | **0.3386** | **0.5018** | **0.8317** | **0.9172** | **`CORE`** (Best Global Performance) |
| **$M_2$ alone** | `esol_delaney_v1` | Single Model | 0.6663 | 0.9221 | 0.4316 | 0.7469 | `SHADOW_ONLY` (Physical Baseline) |
| **$M_3$ alone** | `rdkit_gbr_solubility_v1` | Single Model | 0.7340 | 1.0508 | 0.2619 | 0.7141 | `EXCLUDED_FROM_CONSENSUS` |
| **$M_1 + M_2$** | Admetica + ESOL | Weighted Mean | 0.4153 | 0.5855 | 0.7708 | 0.8792 | Degrades $M_1$ by +0.077 MAE |
| **$M_1 + M_3$** | Admetica + GBR | Weighted Mean | 0.4479 | 0.6402 | 0.7260 | 0.8654 | Degrades $M_1$ by +0.109 MAE |
| **$M_2 + M_3$** | ESOL + GBR | Weighted Mean | 0.6850 | 0.9547 | 0.3907 | 0.7329 | Poor baseline combination |
| **$M_1 + M_2 + M_3$** | All 3 Models | Weighted Mean | 0.3931 | 0.5513 | 0.7968 | 0.8900 | Degrades $M_1$ by +0.055 MAE |
| **$M_1 + M_2 + M_3$** | All 3 Models | Median | 0.5707 | 0.8151 | 0.5558 | 0.7850 | Severe degradation (+0.232 MAE) |

**Mathematical Insight on Median Aggregation**:
When an ensemble contains one high-accuracy neural model ($M_1$) and two less accurate descriptor models ($M_2, M_3$) that both skew in the same direction, **median aggregation picks the central value between the two weak models**, disastrously worsening performance ($\text{MAE} = 0.5707$). Median aggregation is therefore prohibited for asymmetric accuracy ensembles.

---

## 2. Model Ensemble Contribution Taxonomy

Drug-OPT establishes a clear distinction between **Model Scientific Qualification** and **Ensemble Contribution Status**:

| Model Ensemble Status | Definition | Platform Policy |
|---|---|---|
| **`CORE`** | Primary high-performing model that sets the gold standard for the endpoint. | Emits default visible predictions in shadow mode; primary anchor in any weighted combination. |
| **`SUPPORTING`** | High-performing alternative model that adds demonstrated value across multiple subgroups. | Active participant in weighted combinations when promoted. |
| **`SHADOW_ONLY`** | Qualified model providing independent physical/mechanistic insight or error signals, but not beating $M_1$ globally. | Executed and recorded in shadow provenance; evaluated in Stage 4D-3 adaptive weighting research. |
| **`EXCLUDED_FROM_CONSENSUS`** | Model exhibits high collinearity with another member or degrades ensemble accuracy without orthogonal signal. | Retained in audit registry but excluded from active consensus calculation. |

### Summary of Model Role Assignments

| Endpoint | Model ID | Lineage / Family | Qualification Status | Ensemble Contribution Status |
|---|---|---|---|---|
| **Solubility** | `admetica_solubility`<br>`esol_delaney_v1`<br>`rdkit_gbr_solubility_v1` | D-MPNN<br>Delaney Descriptor<br>2D Topological GBR | QUALIFIED<br>QUALIFIED<br>QUALIFIED | **`CORE`**<br>**`SHADOW_ONLY`**<br>**`EXCLUDED_FROM_CONSENSUS`** |
| **Caco-2** | `admetica_caco2`<br>`physchem_caco2_v1` | D-MPNN<br>Polar Surface Mechanistic | QUALIFIED<br>QUALIFIED | **`CORE`**<br>**`SHADOW_ONLY`** |
| **CYP3A4** | `admetica_cyp_cyp3a4-inhibitor`<br>`morgan_cyp3a4_inh_v1` | D-MPNN<br>Morgan ECFP4 Classifier | QUALIFIED<br>QUALIFIED | **`CORE`**<br>**`SHADOW_ONLY`** |
| **hERG** | `admetica_safety_herg`<br>`physchem_herg_v1` | D-MPNN<br>Basic Amine Logistic | QUALIFIED<br>QUALIFIED | **`CORE`**<br>**`SHADOW_ONLY`** |
| **SoM** | `sygma_phase1_2`<br>`smartcyp_dft_v1` | Rule Engine<br>DFT Lookup | QUALIFIED<br>QUALIFIED | **`CORE`**<br>**`SUPPORTING`** ($RRF$) |

---

## 3. Stage 4D-3 Readiness Evaluation

Stage 4D-3 (**Adaptive Experimental Weighting**) is designed to answer:
> *"Can project-specific and series-specific laboratory feedback learn when to dynamically shift weights from $M_1$ to $M_2$?"*

### Readiness Criteria Gate:
1. **Multiple Qualified Models Available**: **PASS** (25 adapters registered across 18 endpoints).
2. **Demonstrated Model Performance Heterogeneity**: **PASS** (Delaney ESOL beats D-MPNN on acyclic/aliphatic scaffolds; Morgan captures distinct substructure alerts).
3. **Project Series Data Structures Available**: **PASS** (Bemis-Murcko series grouping validated).
4. **Prediction Freeze & Shadow Mode Stability**: **PASS** (Zero visible UI regressions; complete shadow provenance).

**Conclusion**: Drug-OPT is scientifically and architecturally **READY FOR STAGE 4D-3 RESEARCH** when authorized by the user.
