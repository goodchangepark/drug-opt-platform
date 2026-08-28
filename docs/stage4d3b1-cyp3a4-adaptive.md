# Stage 4D-3B1: Hierarchical Experimental Adaptive Weighting for CYP3A4 Inhibition

## 1. Executive Summary & Mission
Stage 4D-3B1 extends Drug-OPT's 4-level Bayesian adaptive weighting architecture from continuous regression (Aqueous Solubility) to **binary classification endpoints**, piloting on **CYP3A4 Inhibition** (`EP_MET_CYP3A4_INH`, canonical unit: `probability`).

The central research question:
> *"Can project/series/local experimental CYP3A4 inhibition feedback improve future classification probability and decision performance over the current CORE model without data leakage?"*

### Primary Findings:
1. **Adaptive Consensus Decisively Prevents Static Consensus Degradation**:
   - Unweighted static consensus incurs severe penalty ($\text{Brier} = 0.1046, \text{LogLoss} = 0.3583, \text{MCC} = 0.7974$) compared to $M_1$ CORE ($\text{Brier} = 0.0726, \text{LogLoss} = 0.2646, \text{MCC} = 0.8334$).
   - Adaptive hierarchical weighting restores performance ($\text{Brier} = 0.0728, \text{LogLoss} = 0.2465, \text{MCC} = 0.8334$), achieving lower bounded Log Loss than $M_1$ alone.
2. **Comparison with $M_1$ CORE**:
   - $M_1$ (Admetica D-MPNN) is a strong foundation classifier ($\text{Balanced Accuracy} = 91.66\%, \text{Sensitivity} = 94.26\%, \text{Specificity} = 89.06\%$).
   - $M_2$ (Morgan GBDT) suffers from low specificity ($52.34\%$) and high Brier loss ($0.2056$).
   - The adaptive consensus engine safely constrains $M_2$'s global prior weight to $w_{M2} \approx 0.0422$, preventing corruption of $M_1$'s predictions while permitting local series-specific adaptations.
3. **Scientific Decision**:
   - **`ARCHITECTURE_VALID_BUT_NO_ACCURACY_GAIN`** globally on benchmark evaluation; **`CONDITIONAL_ADAPTIVE_VALUE`** on specific subseries (e.g. Basic Amines, Heteroaromatics).
   - Consensus mode remains strictly **`SHADOW`**; visible production predictions remain 100% untouched $M_1$ CORE.

---

## 2. Mathematical Formulation for Classification Adaptation

### 2.1. Probability Bounding & Loss Metrics
For compound $j$ and model $i$, predictions output estimated posterior probabilities $p_{ij} \in [0, 1]$. To prevent undefined logarithmic penalties $\log(0)$ and extreme weight divergence, probabilities are bounded:
$$\tilde{p}_{ij} = \text{clip}(p_{ij}, \epsilon, 1 - \epsilon), \quad \epsilon = 10^{-4}$$

Given revealed experimental binary truth $y_j \in \{0, 1\}$:
- **Observation Brier Error**:
  $$e_{ij} = (\tilde{p}_{ij} - y_j)^2 \in [0, 1]$$
- **Bounded Logarithmic Loss**:
  $$\mathcal{L}_{ij} = - \left[ y_j \ln(\tilde{p}_{ij}) + (1 - y_j) \ln(1 - \tilde{p}_{ij}) \right]$$

### 2.2. Error Score Transformation
To transform mean empirical Brier loss $\overline{e}_i$ into an unnormalized positive performance score $S_i$:
$$S_i = \left( \frac{1}{\max(0.01, \overline{e}_i)} \right)^\beta, \quad \beta = 3.0$$

### 2.3. 4-Level Evidence Hierarchy & Empirical Bayes Shrinkage
The evidence hierarchy shrinks sparse local observations toward broader, well-calibrated priors:

1. **Level 1: GLOBAL Prior**:
   - Evaluated on external qualification benchmark ($N=250$):
     $$\overline{e}_{M1, \text{glob}} = 0.0726 \implies S_{M1} = 2,613.8$$
     $$\overline{e}_{M2, \text{glob}} = 0.2056 \implies S_{M2} = 115.1$$
     $$w_{M1, \text{glob}} = \frac{2613.8}{2613.8 + 115.1} = 0.9578, \quad w_{M2, \text{glob}} = 0.0422$$

2. **Level 2: PROJECT Level**:
   - Empirical mean Brier loss across project events $\overline{e}_{i, \text{proj}}$.
   - Shrinkage factor:
     $$\lambda_{\text{proj}} = \frac{N_{\text{proj}}}{N_{\text{proj}} + N_{\text{prior, proj}}}, \quad N_{\text{prior, proj}} = 10.0$$
     $$w_{i, \text{proj\_post}} = (1 - \lambda_{\text{proj}}) w_{i, \text{glob}} + \lambda_{\text{proj}} w_{i, \text{proj}}$$

3. **Level 3: SERIES Level (Bemis-Murcko / Functional Cluster)**:
   - Scaffolds identified via canonical Bemis-Murcko rings; acyclic molecules partitioned into functional clusters (`[acyclic_Alcohol]`, `[acyclic_Amine]`, etc.).
   - Shrinkage factor:
     $$\lambda_{\text{ser}} = \frac{N_{\text{ser}}}{N_{\text{ser}} + N_{\text{prior, ser}}}, \quad N_{\text{prior, ser}} = 5.0$$
     $$w_{i, \text{ser\_post}} = (1 - \lambda_{\text{ser}}) w_{i, \text{proj\_post}} + \lambda_{\text{ser}} w_{i, \text{ser}}$$

4. **Level 4: LOCAL Neighborhood (Morgan Fingerprint Similarity)**:
   - Tanimoto distance $T_j \ge 0.40$ with radius=2, 2048 bits.
   - Effective sample size $N_{\text{local, eff}} = \sum_j T_j^2$.
   - Shrinkage factor:
     $$\lambda_{\text{loc}} = \frac{N_{\text{local, eff}}}{N_{\text{local, eff}} + 3.0}$$
     $$w_{i, \text{loc\_post}} = (1 - \lambda_{\text{loc}}) w_{i, \text{ser\_post}} + \lambda_{\text{loc}} w_{i, \text{loc}}$$

5. **Level 5: Applicability Domain & Weight Floor**:
   - Applicability factor $\gamma_{\text{AD}} \in \{1.0 \text{ (In Domain)}, 0.5 \text{ (Borderline)}, 0.1 \text{ (OOD)}\}$.
   - Minimum weight floor $\epsilon_{\text{floor}} = 0.02$.
   - Final normalized effective weights:
     $$w_{i, \text{eff}} = \frac{\max(\epsilon_{\text{floor}}, w_{i, \text{loc\_post}} \cdot \gamma_{\text{AD}, i})}{\sum_k \max(\epsilon_{\text{floor}}, w_{k, \text{loc\_post}} \cdot \gamma_{\text{AD}, k})}$$

---

## 3. Production Safety & Governance
- **Shadow Mode Only**: Consensus predictions are written to `ADAPTIVE_SHADOW_CONSENSUS` and never displayed as default production predictions in UI.
- **Strict Endpoint Contract**: Incompatible experimental records (e.g. CYP3A4 substrate turnover, $K_i$, time-dependent inhibition, mechanism-based inactivation, induction) are rejected with `EXPERIMENT_NOT_ADAPTATION_COMPATIBLE`.
- **Zero Retrospective Leakage**: Future feedback events after prediction timestamp $t$ are strictly excluded from history during sequential forward replay.
