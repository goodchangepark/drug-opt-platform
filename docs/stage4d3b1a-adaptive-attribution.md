# Stage 4D-3B1A: CYP3A4 Adaptive Attribution Audit

## 1. Executive Summary & Scientific Mission

Stage 4D-3B1A conducts an authoritative attribution audit of the hierarchical adaptive weighting architecture on the **CYP3A4 Inhibition** classification endpoint (`EP_MET_CYP3A4_INH`).

The central scientific question:
> *"Does the observed CYP3A4 performance improvement in Stage 4D-3B1 arise from genuine dynamic experimental adaptation (across project, series, and local chemical neighborhoods), or is it merely an artifact of applying a conservative fixed global blend ($w_{M1} \approx 0.9578, w_{M2} \approx 0.0422$)?"*

### Authoritative Findings:
1. **Dynamic Experimental Feedback Adds Zero Predictive Value Beyond Fixed Global Blend**:
   - **Fixed Global Prior (No Adaptation)**: $\text{Brier} = 0.0726$, $\text{LogLoss} = 0.2460$, $\text{MCC} = 0.8334$, $\text{Balanced Acc} = 91.66\%$.
   - **Full Adaptive (Global + Project + Series + Local)**: $\text{Brier} = 0.0735$, $\text{LogLoss} = 0.2535$, $\text{MCC} = 0.8259$, $\text{Balanced Acc} = 91.27\%$.
   - Paired Bootstrap (1,000 resamples) demonstrates that Full Adaptive is worse than Fixed Global Prior with $P(\text{Adaptive better}) = 0.006$ ($0.6\%$) on Brier score and $P(\text{Adaptive better}) = 0.004$ ($0.4\%$) on Log Loss.
2. **Log Loss Improvement Explained by Probability Softening, Not Dynamic Learning**:
   - $M_1$ (Admetica D-MPNN) raw Log Loss is $0.2646$ due to large $-\ln(\epsilon)$ penalties on a tiny set of 9 overconfident mispredictions ($p \ge 0.80$ with $y=0$ or $p \le 0.20$ with $y=1$).
   - Blending with a fixed $4.2\%$ fraction of $M_2$ (Morgan GBDT) softens boundary probabilities (e.g. $0.0000 \to 0.0003$), lowering Log Loss to $0.2460$ without changing binary classification thresholds.
   - Dynamic updating actually introduces small stochastic weight fluctuations that slightly inflate Log Loss to $0.2535$.
3. **Re-evaluation and Correction of Stage 4D-3B1 Subgroup Claims**:
   - The Stage 4D-3B1 hypothesis of conditional dynamic value for *Basic Amines* and *Neutral Heteroaromatics* is **not supported by empirical evidence**.
   - Subgroup comparisons show Fixed Global Prior matches or outperforms Full Adaptive across every chemical subtype.
4. **Final Scientific Verdict**:
   - Scientific Decision: **`FIXED_GLOBAL_BLEND_SUFFICIENT`**.
   - $M_2$ Role: **`CALIBRATION_SUPPORTING`** (and `SHADOW_ONLY`).
   - hERG Gate: **`GO_HERG_CALIBRATION_AUDIT_FIRST`**.

---

## 2. 7-Strategy Component Ablation ($N=250$)

All 7 strategies were evaluated on the identical frozen authoritative evaluation cohort ($N=250$, 122 positives / 128 negatives):

| Strategy | Role / Scope | MCC | Balanced Accuracy | Brier Score | Bounded Log Loss | AUROC | AUPRC | Sensitivity | Specificity |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **1. $M_1$ CORE (Admetica D-MPNN)** | Primary Production Model | **0.8334** | **0.9166** | **0.0726** | 0.2646 | 0.9572 | 0.9217 | **0.9426** | **0.8906** |
| **2. $M_2$ SHADOW (Morgan GBDT)** | Secondary Shadow Model | 0.4365 | 0.7043 | 0.2056 | 0.6076 | 0.7627 | 0.6983 | 0.8852 | 0.5234 |
| **3. 50/50 Static Consensus** | Unweighted Stage 4D-1 Blend | 0.7974 | 0.8973 | 0.1046 | 0.3583 | 0.9483 | 0.9312 | 0.9508 | 0.8438 |
| **4. Fixed Global Prior** | $w_{M1}=0.9578, w_{M2}=0.0422$ (No Adaptation) | **0.8334** | **0.9166** | **0.0726** | **0.2460** | **0.9594** | **0.9483** | **0.9426** | **0.8906** |
| **5. Global + Project** | Project evidence ($\lambda_{\text{proj}} > 0$) | 0.8259 | 0.9127 | 0.0735 | 0.2532 | 0.9576 | 0.9416 | **0.9426** | 0.8828 |
| **6. Global + Project + Series** | Project + Series evidence ($\lambda_{\text{ser}} > 0$) | 0.8259 | 0.9127 | 0.0735 | 0.2533 | 0.9576 | 0.9416 | **0.9426** | 0.8828 |
| **7. Full Adaptive** | All 4 Levels (Global+Project+Series+Local) | 0.8259 | 0.9127 | 0.0735 | 0.2535 | 0.9575 | 0.9411 | **0.9426** | 0.8828 |

---

## 3. Paired Bootstrap Analysis (1,000 Resamples)

### 3.1. Primary Test: Full Adaptive vs Fixed Global Prior

$$\Delta = \text{Metric}_{\text{Full Adaptive}} - \text{Metric}_{\text{Fixed Global Prior}}$$

| Metric | Mean $\Delta$ | Median $\Delta$ | 95% Confidence Interval | $P(\text{Adaptive Better})$ | Verdict |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Brier Score** | $+0.0009$ | $+0.0009$ | $[+0.0002, +0.0018]$ | $0.0060$ ($0.6\%$) | **Fixed Global Superior** |
| **Bounded Log Loss** | $+0.0075$ | $+0.0071$ | $[+0.0022, +0.0151]$ | $0.0040$ ($0.4\%$) | **Fixed Global Superior** |
| **MCC** | $-0.0074$ | $-0.0075$ | $[-0.0231, 0.0000]$ | $0.0000$ ($0.0\%$) | **Fixed Global Superior** |
| **Balanced Accuracy** | $-0.0038$ | $-0.0038$ | $[-0.0123, 0.0000]$ | $0.0000$ ($0.0\%$) | **Fixed Global Superior** |

### 3.2. Secondary Comparators
- **Fixed Global Prior vs $M_1$ CORE**:
  - $\Delta \text{Brier} = 0.0000$ ($[-0.0005, +0.0005]$).
  - $\Delta \text{LogLoss} = -0.0186$ ($[-0.0673, +0.0028]$), $P(\text{Fixed Better}) = 0.963$.
  - $\Delta \text{MCC} = 0.0000$, $\Delta \text{Balanced Acc} = 0.0000$.

---

## 4. Weight Movement Attribution & Trajectory Analysis

Tracking effective weight movements across all 250 prospective sequential decisions:

- **Absolute Shift from Global Prior** $|\Delta w| = |w_{M1, \text{eff}} - w_{M1, \text{glob}}|$:
  - **Median**: $0.0118$
  - **75th Percentile (P75)**: $0.0270$
  - **90th Percentile (P90)**: $0.0816$
  - **Maximum**: $0.1843$
  - **Mean**: $0.0262$
  - Compounds within $\pm 0.01$ of Global Prior: $44.8\%$
  - Compounds within $\pm 0.05$ of Global Prior: $81.6\%$
- **Effective Model Weights**:
  - $w_{M1, \text{eff}}$: Mean $= 0.9376$, Range $= [0.7735, 0.9800]$
  - $w_{M2, \text{eff}}$: Mean $= 0.0624$, Range $= [0.0200, 0.2265]$

Because the weight shift is constrained by the shrinkage prior ($N_{\text{prior, proj}}=10.0, N_{\text{prior, ser}}=5.0$) and the initial global prior is heavily skewed ($0.9578 / 0.0422$), the algorithm operates effectively as a static blend with minor noisy perturbations.

---

## 5. Dynamic Value by Project & Scaffold Series

### 5.1. Pseudo-Project Attribution
| Project | $N$ | Positive Frac. | $M_1$ Brier | Fixed Global Brier | Adaptive Brier | Fixed Global LogLoss | Adaptive LogLoss | Fixed MCC | Adaptive MCC | Classification |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **PROJ_01** | 50 | 0.460 | 0.0715 | 0.0716 | 0.0726 | 0.2398 | 0.2444 | 0.8803 | 0.8803 | `EQUIVALENT` |
| **PROJ_02** | 50 | 0.400 | 0.0659 | 0.0666 | 0.0669 | 0.2229 | 0.2268 | 0.8498 | 0.8498 | `EQUIVALENT` |
| **PROJ_03** | 50 | 0.460 | 0.0732 | 0.0725 | 0.0726 | 0.2351 | 0.2392 | 0.7987 | 0.7987 | `EQUIVALENT` |
| **PROJ_04** | 50 | 0.580 | 0.0955 | 0.0956 | 0.0983 | 0.3248 | 0.3368 | 0.7537 | 0.7114 | `GLOBAL_PRIOR_BETTER` |
| **PROJ_05** | 50 | 0.540 | 0.0569 | 0.0566 | 0.0573 | 0.2073 | 0.2203 | 0.8797 | 0.8797 | `GLOBAL_PRIOR_BETTER` |

Zero projects exhibited `ADAPTIVE_BETTER`. In 2 out of 5 projects, Fixed Global Prior was clearly better due to avoiding negative transfer from small-sample feedback.

### 5.2. Chemical Scaffold Series Attribution
Across all populated Bemis-Murcko scaffold clusters ($N \ge 4$), zero series demonstrated statistically significant or reproducible dynamic adaptive gains over the Fixed Global Prior.

---

## 6. Re-Audit of Stage 4D-3B1 Subgroup Claims

In the Stage 4D-3B1 report, it was hypothesized that dynamic adaptation provided conditional value for basic amines and heteroaromatic ring systems. The matched re-audit reveals:

| Subgroup | $N$ | Pos. Frac. | $M_1$ Brier | Fixed Global Brier | Adaptive Brier | $M_1$ LogLoss | Fixed Global LogLoss | Adaptive LogLoss | Claim Audit Verdict |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Basic Amine (+)** | 106 | 0.340 | 0.0672 | **0.0667** | 0.0682 | 0.2899 | **0.2356** | 0.2484 | `NOT SUPPORTED` (Fixed Global Best) |
| **Basic Amine (-)** | 144 | 0.597 | **0.0765** | 0.0770 | 0.0775 | **0.2459** | 0.2536 | 0.2572 | `NOT SUPPORTED` (Fixed Global $\approx$ M1) |
| **Heteroaromatic (+)** | 158 | 0.468 | **0.0755** | 0.0756 | 0.0762 | **0.2496** | 0.2521 | 0.2567 | `NOT SUPPORTED` (Fixed Global $\approx$ M1) |
| **Heteroaromatic (-)** | 92 | 0.522 | 0.0674 | **0.0675** | 0.0689 | 0.2903 | **0.2355** | 0.2481 | `NOT SUPPORTED` (Fixed Global Best) |
| **Neutral Heteroaromatics**| 148 | 0.486 | **0.0675** | 0.0679 | 0.0688 | **0.2271** | 0.2318 | 0.2378 | `NOT SUPPORTED` (Fixed Global $\approx$ M1) |

**Conclusion**: The earlier claim of conditional adaptive value was an artifact of comparing the dynamic model against unweighted static consensus or raw $M_1$ rather than against the matched Fixed Global Prior.

---

## 7. Prospective Learning Curve (Feedback Volume Sensitivity)

Evaluating future observations after accumulating $K$ prior experimental labels:

| Cumulative Prior Labels ($K$) | Remaining Evaluation Compounds | Fixed Global Brier | Adaptive Brier | Fixed Global LogLoss | Adaptive LogLoss | Fixed MCC | Adaptive MCC |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **0** | 250 | **0.0726** | 0.0735 | **0.2460** | 0.2535 | **0.8334** | 0.8259 |
| **5** | 245 | **0.0740** | 0.0749 | **0.2498** | 0.2575 | **0.8300** | 0.8224 |
| **10** | 240 | **0.0728** | 0.0738 | **0.2471** | 0.2550 | **0.8330** | 0.8254 |
| **20** | 230 | **0.0714** | 0.0724 | **0.2442** | 0.2525 | **0.8378** | 0.8302 |
| **30** | 220 | **0.0726** | 0.0736 | **0.2475** | 0.2556 | **0.8336** | 0.8259 |
| **50** | 200 | **0.0728** | 0.0738 | **0.2475** | 0.2558 | **0.8334** | 0.8256 |

**Finding**: Performance does not improve with increasing prior feedback. The gap between Fixed Global and Adaptive remains constant ($\Delta \text{Brier} \approx +0.0010, \Delta \text{LogLoss} \approx +0.0080$).

---

## 8. Negative Control Audit & Interpretation

- **Real Adaptive Brier**: $0.0735$
- **Shuffled Feedback Brier**: $0.1014$
- **Fixed Global Prior Brier**: $0.0726$

### Methodological Interpretation:
- Real Adaptive significantly outperforms Shuffled Feedback ($0.0735$ vs $0.1014$), proving that the algorithm is sensitive to and correctly updates from sequential feedback.
- However, Real Adaptive fails to outperform Fixed Global Prior ($0.0735$ vs $0.0726$).
- Therefore, the negative control validates code execution and signal pipeline integrity, but scientific attribution shows that the dynamic feedback does not add value beyond the conservative fixed global blend.

---

## 9. Final Decision & Governance

1. **Scientific Decision**: **`FIXED_GLOBAL_BLEND_SUFFICIENT`**.
2. **Model Status**:
   - $M_1$ (Admetica D-MPNN): **`CORE`** (Primary predictive engine).
   - $M_2$ (Morgan GBDT): **`CALIBRATION_SUPPORTING`** (Shadow only; serves only as a probability regularizer).
3. **Consensus Mode**: Strictly **`SHADOW`**. Production predictions remain 100% $M_1$ CORE.
4. **hERG Liability Gate**: **`GO_HERG_CALIBRATION_AUDIT_FIRST`**.
   - Prior validation of hERG models showed marked specificity deficits in Morgan GBDT models.
   - Extending adaptive weighting directly to hERG without first assessing fixed global calibration would obscure underlying model deficiencies.
   - We recommend executing a calibration audit of single models and fixed conservative mixtures before considering any adaptive infrastructure on hERG.
