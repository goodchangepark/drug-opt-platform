# Stage 4D-0: Future Multi-Model Ensemble Architecture & Governance

**Platform Version:** `0.6.3-stage5b4-ui`  
**Stage:** `4D-0 Foundation Design`  
**Status:** Pre-Implementation Specification & Engineering Contracts

---

## 1. Architecture Overview & Execution Tiers

The Drug-OPT multi-model ensemble system is designed for deterministic, fault-tolerant execution on resource-constrained ARM64 hardware (NVIDIA Jetson AGX Xavier). Models are segregated into three distinct execution tiers based on computational complexity and latency:

```
+-----------------------------------------------------------------------------------+
|                           Drug-OPT Ensemble Orchestrator                          |
+-----------------------------------------------------------------------------------+
                                        |
      +---------------------------------+---------------------------------+
      |                                 |                                 |
      v                                 v                                 v
[ Tier 1: Local Fast ]        [ Tier 2: Local Heavy ]        [ Tier 3: Optional Ext ]
• Latency: < 5 ms             • Latency: 10 - 50 ms          • Latency: > 200 ms
• CPU Footprint: < 100 MB     • CPU Footprint: 200 - 500 MB  • Distributed Workers
• Tree Models (CatBoost/GBDT) • Deep MPNNs (Chemprop v2)     • Async Cloud Service
• RDKit 2D / Morgan FP        • 5-Fold Ensemble Inference    • Sandbox Tools
• Deterministic Physics       • Local ARM64 Torch Runtime    • Optional / Fallback
```

---

## 2. Multi-Level Weighting Scope Design

Future ensemble weights will NOT be static global constants. The weighting architecture establishes four hierarchical scopes with Bayesian shrinkage:

```
[Level 1: Global Prior]
  └── Default benchmark weights based on public held-out cross-validation.
        │
        ▼
[Level 2: Project Scope]
  └── Weights adjusted by project-specific experimental assay comparisons.
        │
        ▼
[Level 3: Chemical Series / Scaffold Scope]
  └── Bemis-Murcko framework and Maximum Common Substructure (MCS) clustering.
      Weights favor models with proven accuracy on the target core scaffold.
        │
        ▼
[Level 4: Local Chemical Neighborhood Scope]
  └── Nearest neighbor Tanimoto similarity ($k$-NN in Morgan fingerprint space).
      Weights dynamically upweight models whose training set covers the immediate local analog space.
```

### Effective Sample Size & Bayesian Shrinkage Policy
To prevent overfitting on small project datasets ($N < 10$), project-level and scaffold-level weights must shrink towards the global prior via a regularized shrinkage factor $\alpha$:

$$\alpha = \frac{N_{\text{eff}}}{N_{\text{eff}} + N_{\text{prior}}}$$

$$w_{\text{final}} = \alpha \cdot w_{\text{local}} + (1 - \alpha) \cdot w_{\text{global}}$$

Where $N_{\text{prior}}$ defaults to 15 compounds, ensuring smooth, conservative transitions as project experimental data accumulates.

---

## 3. Error Metrics & Loss Formulations

Ensemble weighting and model comparisons must adhere to the following mathematically rigorous metric contracts:

### 3.1 Continuous Regression Endpoints (Log-Scale)
* **Log-Space Mean Absolute Error (MAE):**
  $$\text{MAE}_{\log} = \frac{1}{N} \sum_{i=1}^N |\log_{10}(y_i) - \log_{10}(\hat{y}_i)|$$
* **Root Mean Squared Error (RMSE):**
  $$\text{RMSE}_{\log} = \sqrt{\frac{1}{N} \sum_{i=1}^N (\log_{10}(y_i) - \log_{10}(\hat{y}_i))^2}$$
* **Fold-Error Acceptance Bounds:**
  $$\text{Fold Error} = 10^{|\log_{10}(y_i) - \log_{10}(\hat{y}_i)|}$$
  * $2\times$ Fold Error Window: $|\log_{10}(y) - \log_{10}(\hat{y})| \le \log_{10}(2) \approx 0.301$
  * $3\times$ Fold Error Window: $|\log_{10}(y) - \log_{10}(\hat{y})| \le \log_{10}(3) \approx 0.477$

### 3.2 Binary Classification Endpoints
* **Balanced Accuracy:**
  $$\text{BA} = \frac{1}{2} \left( \frac{\text{TP}}{\text{TP} + \text{FN}} + \frac{\text{TN}}{\text{TN} + \text{FP}} \right)$$
* **Matthews Correlation Coefficient (MCC):**
  $$\text{MCC} = \frac{\text{TP} \cdot \text{TN} - \text{FP} \cdot \text{FN}}{\sqrt{(\text{TP}+\text{FP})(\text{TP}+\text{FN})(\text{TN}+\text{FP})(\text{TN}+\text{FN})}}$$
* **Brier Calibration Score:**
  $$\text{Brier} = \frac{1}{N} \sum_{i=1}^N (p_i - y_i)^2$$

---

## 4. Prospective Prediction Freeze Contract

To eliminate hindsight bias and ensure scientific auditability:
1. **Freeze Trigger:** When a compound version is saved or an optimization run generated, all model predictions, individual member outputs, applicability domains, and conformal intervals are saved in an immutable database snapshot.
2. **Immutability:** Once frozen, a prediction record is NEVER overwritten when new experimental data is uploaded.
3. **Comparative Provenance:** When experimental assay data arrives, the platform logs the exact prediction error against the frozen prospective baseline.

---

## 5. Fault-Tolerant Failure Isolation Contract

The ensemble orchestrator enforces strict fault isolation. If one model fails during inference (due to numerical instability, missing descriptor generation, or checkpoint read timeout), the remaining models continue execution cleanly:

```python
from backend.endpoint_contracts import execute_fault_tolerant_ensemble, get_endpoint_contract

# If Model 1 and Model 2 succeed, but Model 3 raises an unhandled exception:
result = execute_fault_tolerant_ensemble(
    adapters=[adapter_1, adapter_2, adapter_3],
    canonical_smiles="CC(=O)Oc1ccccc1C(=O)O",
    contract=get_endpoint_contract("Solubility"),
)

# Execution completes successfully:
assert result.is_valid is True
assert result.member_count == 2
assert len(result.failed_models) == 1
```

This guarantees 100% platform uptime even under intermittent single-model failures.
