# Stage 4D-3A: Hierarchical Experimental Adaptive Weighting Engine

## 1. Overview & Scientific Purpose
Stage 4D-3A introduces Drug-OPT's first hierarchical adaptive experimental weighting system for **Aqueous Solubility** (`EP_PHYS_SOLUBILITY`).

In Stage 4D-2C, static consensus was audited and demonstrated to suffer from negative contribution when combining a superior primary model ($M_1$: `admetica_solubility`, $\text{MAE}=0.3386$) with weaker static models ($M_2$: `esol_delaney_v1`, $\text{MAE}=0.6663$). Static consensus unconditionally degraded accuracy from $\text{MAE}=0.3386$ to $\text{MAE}=0.3931$ across the full Delaney cohort.

Stage 4D-3A replaces static weighting with a **4-level hierarchical Bayesian shrinkage engine** that:
1. Respects global validation priors when project-specific data is sparse.
2. Adaptively learns empirical performance at the **Project**, **Chemical Series (Bemis-Murcko scaffold)**, and **Local Chemical Neighborhood (Morgan Tanimoto)** levels.
3. Prevents retrospective leakage by strictly processing historical experiments recorded prior to prediction timestamps.
4. Operates strictly in **SHADOW mode** (`consensus_mode = "SHADOW"`), ensuring zero alterations to visible primary production predictions and preserving the 100% UI design freeze.

---

## 2. 4-Level Hierarchical Evidence Architecture

```
GLOBAL PRIOR (Stage 4D-2C Benchmark Validation)
    ↓
PROJECT EVIDENCE (Project-specific experimental observations, N_prior = 10)
    ↓
SERIES EVIDENCE (Bemis-Murcko Scaffold Series, N_prior = 5)
    ↓
LOCAL NEIGHBORHOOD EVIDENCE (Morgan radius=2 2048-bit T >= 0.40, N_prior = 3)
    ↓
APPLICABILITY DOMAIN & WEIGHT FLOOR (gamma_AD in {1.0, 0.5, 0.1}, eps = 0.02)
    ↓
ADAPTIVE SHADOW CONSENSUS PREDICTION
```

### Level 1: Global Validation Prior
Global error metrics are derived from frozen benchmark qualification:
- $M_1$ (`admetica_solubility`): $\text{MAE} = 0.3386$
- $M_2$ (`esol_delaney_v1`): $\text{MAE} = 0.6663$
- $M_3$ (`rdkit_gbr_solubility_v1`): $\text{MAE} = 0.7340$

Performance scores are computed via exponential error transformation:
$$S_{i, \text{global}} = \exp(-\beta \cdot \text{MAE}_{i, \text{global}})$$
with $\beta = 2.0$, yielding normalized global prior weights:
- $w_{M1, \text{global}} = 0.6582$
- $w_{M2, \text{global}} = 0.3418$

### Level 2: Project-Level Shrinkage
For a project with $N_{\text{project}}$ valid prior experimental measurements:
$$\lambda_{\text{project}} = \frac{N_{\text{project}}}{N_{\text{project}} + N_{\text{prior, project}}}, \quad N_{\text{prior, project}} = 10.0$$
$$\mathbf{w}_{\text{proj\_post}} = (1 - \lambda_{\text{project}}) \mathbf{w}_{\text{global}} + \lambda_{\text{project}} \mathbf{w}_{\text{proj\_emp}}$$

### Level 3: Series-Level Shrinkage
Molecules are partitioned into chemical series using canonical Bemis-Murcko scaffolds (with acyclic structures labeled `[acyclic]`):
$$\lambda_{\text{series}} = \frac{N_{\text{series}}}{N_{\text{series}} + N_{\text{prior, series}}}, \quad N_{\text{prior, series}} = 5.0$$
$$\mathbf{w}_{\text{ser\_post}} = (1 - \lambda_{\text{series}}) \mathbf{w}_{\text{proj\_post}} + \lambda_{\text{series}} \mathbf{w}_{\text{ser\_emp}}$$

### Level 4: Local Neighborhood Shrinkage
Local structural analogs are retrieved using Morgan fingerprints (radius=2, 2048-bit) with Tanimoto similarity $T_j \ge 0.40$:
$$N_{\text{eff, local}} = \sum_{j \in \text{neighbors}} T_j^2$$
$$\text{MAE}_{i, \text{local}} = \frac{\sum_j T_j \cdot |\hat{y}_{ij} - y_j|}{\sum_j T_j}$$
$$\lambda_{\text{local}} = \frac{N_{\text{eff, local}}}{N_{\text{eff, local}} + N_{\text{prior, local}}}, \quad N_{\text{prior, local}} = 3.0$$
$$\mathbf{w}_{\text{loc\_post}} = (1 - \lambda_{\text{local}}) \mathbf{w}_{\text{ser\_post}} + \lambda_{\text{local}} \mathbf{w}_{\text{loc\_emp}}$$

### Level 5: Applicability Domain Scaling & Weight Floor
- Applicability domain multiplier $\gamma_{\text{AD}} \in \{1.0 \text{ (IN\_DOMAIN)}, 0.5 \text{ (BORDERLINE)}, 0.1 \text{ (OUT\_OF\_DOMAIN)}\}$
- Minimum weight floor $\epsilon = 0.02$
$$w_{i, \text{final}} = \text{Normalize}(\max(\epsilon, w_{i, \text{loc\_post}} \cdot \gamma_{\text{AD}, i}))$$

---

## 3. Provenance & Reason Codes
Every adaptive prediction emits deterministic machine-readable reason codes:
- `GLOBAL_PRIOR_DOMINANT`: $N_{\text{project}} < 5$; prediction relies primarily on qualified global benchmark weights.
- `PROJECT_EVIDENCE_ACTIVE`: $N_{\text{project}} \ge 5$; project-level experimental feedback is actively influencing weights.
- `SERIES_M2_OUTPERFORMS_M1`: Series experimental observations demonstrate lower MAE for $M_2$ over $M_1$.
- `SERIES_M1_OUTPERFORMS_M2`: Series experimental observations demonstrate $M_1$ superiority.
- `LOCAL_NEIGHBORHOOD_ACTIVE`: Similar chemical analogs ($T \ge 0.40, N_{\text{eff}} \ge 1.0$) are modulating local weights.
- `INSUFFICIENT_LOCAL_DATA`: No chemical analogs meet the $T \ge 0.40$ threshold; local level gracefully falls back to series/project posteriors.
- `M3_ADAPTIVE_EXCLUDED`: $M_3$ (RDKit GBR) is excluded from production adaptation due to verified negative contribution.
- `UNSTABLE_ADAPTIVE_WEIGHTS`: Flagged if step-to-step weight jump exceeds $0.35$.
- `NO_FROZEN_PREDICTION`: Recorded if experimental measurement has no prior frozen model predictions.

---

## 4. API Endpoints
- `GET /api/compound-versions/{version_id}/adaptive-provenance`:
  Returns complete hierarchical weight breakdown across Global, Project, Series, and Local levels, sample counts, scaffold series ID, effective weights, and reason codes.
