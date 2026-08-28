# Stage 4D-1: Static Consensus & Aggregation Policy

## 1. Scope & Foundation Purpose

Stage 4D-1 implements an initial **static consensus and aggregation framework** informed strictly by model qualification, applicability domain, and representation diversity.

This policy does **NOT** implement adaptive experimental weighting (which is deferred to Stage 4D-2). Instead, it establishes the deterministic mathematical rules for combining multiple model predictions while preventing artificial consensus distortion.

---

## 2. Endpoint-Specific Aggregation Types

A single generic averaging function is scientifically invalid across diverse endpoints. Drug-OPT enforces four distinct aggregation paradigms:

```
Endpoint Type               Aggregation Strategy          Dispersion Metric
-----------------------------------------------------------------------------------------
Continuous Regression       REGRESSION_WEIGHTED           Weighted Std Dev (Model Disagreement)
Binary Classification       CLASSIFICATION_WEIGHTED       Vote Pattern & Majority Consensus
Metabolic Soft Spots (SoM)  RANK_FUSION                   Reciprocal Rank Fusion (RRF)
Pharmacokinetics (NCA/ODE)  NO_CONSENSUS                  Excluded from ADMET ensemble
```

### 2.1 Continuous Regression Aggregation
For continuous endpoints (Solubility, Caco-2, PPB, HLM, RLM, MLM):
* **Weighted Consensus Mean**:
  $$\bar{x}_w = \sum_{i=1}^M w_i x_i \quad \text{where } \sum_{i=1}^M w_i = 1$$
* **Model Disagreement (Weighted Standard Deviation)**:
  $$\sigma_w = \sqrt{\sum_{i=1}^M w_i (x_i - \bar{x}_w)^2}$$
  *(Note: $\sigma_w$ represents Model Disagreement across different neural network architectures; it is NOT a statistical confidence interval).*
* **Spread / Range**: $[\min(x_i), \max(x_i)]$.

### 2.2 Binary Classification Aggregation
For probability endpoints (CYP inhibitors, CYP substrates, P-gp, hERG, Ames, DILI):
* **Weighted Consensus Probability**:
  $$\bar{p}_w = \sum_{i=1}^M w_i p_i$$
* **Consensus Classification**:
  $$\text{Class} = \begin{cases} \text{POSITIVE}, & \text{if } \bar{p}_w \ge \theta_{\text{threshold}} \\ \text{NEGATIVE}, & \text{otherwise} \end{cases}$$
* **Model Vote Pattern**: Recorded as an explicit string (e.g. `M1:POSITIVE, M2:POSITIVE, M3:NEGATIVE`).

### 2.3 Site-of-Metabolism (SoM) Rank Fusion
For atom-level metabolic soft spots, numeric value averaging is strictly prohibited. Atom positions are aggregated using **Reciprocal Rank Fusion (RRF)**:
$$\text{RRF\_Score}(a) = \sum_{m \in \text{models}} \frac{w_m}{60 + \text{rank}_m(a)}$$

---

## 3. Weight Formulation & Diversity Penalty

The static weight $w_i$ of each model is computed transparently:

$$w_i = \text{BaseQuality}_i \times \text{ApplicabilityDomain}_i \times \text{Confidence}_i \times \text{DiversityPenalty}_i$$

1. **Base Quality Multiplier**:
   * `QUALIFIED`: $1.00$
   * `QUALIFIED_WITH_LIMITATIONS`: $0.85$
   * `RESEARCH_ONLY` / `REJECTED` / `UNAVAILABLE`: $0.00$ (excluded)
2. **Applicability Domain Multiplier**:
   * `IN_DOMAIN`: $1.00$
   * `BORDERLINE`: $0.70$
   * `OUT_OF_DOMAIN`: $0.10$ (prevents OOD predictions from dominating consensus)
3. **Confidence Multiplier**:
   * `HIGH`: $1.00$
   * `MEDIUM`: $0.85$
   * `LOW`: $0.65$
4. **Diversity Penalty Factor**:
   * If two models share the exact same training dataset (e.g. AqSolDB, AID 1851) and inductive bias (Chemprop D-MPNN), they receive a **$0.55\times$ diversity factor** so their combined vote does not double-count identical evidence.

---

## 4. Failure Renormalization

If a model $M_k$ produces a `RUNTIME_ERROR` or `MODEL_UNAVAILABLE`, the remaining valid models are renormalized:

$$w_{i,\text{effective}} = \frac{w_i}{\sum_{j \in \text{valid}} w_j}$$

Both `original_weights` and `effective_weights` are stored immutably in the `admet_consensus_predictions` table.

---

## 5. Model Agreement Classification

Model agreement is classified into deterministic categories:

* **`HIGH_AGREEMENT`**: Regression $\sigma_w \le 0.30$ (within ~2-fold error window) OR 100% unanimous classification vote.
* **`MODERATE_AGREEMENT`**: $0.30 < \sigma_w \le 0.60$ OR majority classification vote ($\ge 70\%$).
* **`LOW_AGREEMENT`**: $\sigma_w > 0.60$ OR split classification vote ($< 70\%$).
* **`SINGLE_MODEL`**: Only 1 qualified model active.
