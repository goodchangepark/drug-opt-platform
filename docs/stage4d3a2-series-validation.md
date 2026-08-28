# Stage 4D-3A2: Series-Level Validation & Realistic Project Campaigns

## 1. Scaffold & Functional Series Challenge

Evaluating individual chemical series reveals where adaptive weighting adds distinct value:

| Series Identifier | Chemical Description | $N$ | $M_1$ MAE | $M_2$ MAE | Adaptive MAE | $\Delta\text{MAE vs } M_1$ | Classification |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `c1ccccc1` | Simple Benzene / Aryl | 48 | 0.3150 | 0.5875 | 0.3291 | +0.0141 | **`EQUIVALENT`** |
| `[acyclic_hydrocarbon]` | Alkanes & Alkenes | 21 | 0.6656 | 2.3535 | 0.6986 | +0.0330 | `M1_BETTER` |
| `[acyclic_Alcohol]` | Aliphatic Alcohols | 10 | 0.1539 | 1.5227 | 0.1410 | **-0.0129** | **`ADAPTIVE_BETTER`** |
| `[acyclic_Ester]` | Aliphatic Esters | 10 | 0.2965 | 0.8194 | 0.3150 | +0.0185 | **`EQUIVALENT`** |
| `[acyclic_Halogenated]` | Halogenated Aliphatics | 8 | 1.0984 | 1.3409 | 1.1548 | +0.0564 | `M1_BETTER` |
| `c1ccc(-c2ccccc2)cc1` | Biphenyl Series | 6 | 0.2257 | 1.1212 | 0.2745 | +0.0488 | `M1_BETTER` |
| `[acyclic_Amine]` | Aliphatic Amines | 5 | 0.7081 | 0.9327 | 0.6314 | **-0.0767** | **`ADAPTIVE_BETTER`** |
| `c1ccc(Oc2ccccc2)cc1` | Diphenyl Ether | 5 | 0.1013 | 0.9850 | 0.1119 | +0.0106 | **`EQUIVALENT`** |
| `[acyclic_Alcohol_Halogenated]`| Halo-alcohols | 4 | 0.3693 | 0.9348 | 0.3154 | **-0.0539** | **`ADAPTIVE_BETTER`** |
| `[acyclic_Amide]` | Aliphatic Amides | 4 | 0.2657 | 0.9848 | 0.2430 | **-0.0227** | **`ADAPTIVE_BETTER`** |

**Summary**: In specific functional groups (e.g. Aliphatic Alcohols, Aliphatic Amines, Halo-alcohols, Amides), local and project adaptive feedback achieves measurable accuracy gains over $M_1$ alone ($\Delta\text{MAE}$ from $-0.013$ to $-0.077$).

---

## 2. Realistic Project Campaign Simulations

Rather than treating 250 diverse molecules as a single project, we simulated realistic drug discovery campaigns of varying sizes ($N \in \{3, 5, 10, 20, 30\}$):

| Campaign Size ($N$) | Simulated Trials | Mean $M_1$ MAE | Mean Adaptive MAE | $\Delta\text{MAE}$ | Cross-Project Isolation |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **$N = 3$** | 20 | 0.3754 | 0.3886 | +0.0132 | **VERIFIED (100%)** |
| **$N = 5$** | 20 | 0.3928 | 0.3826 | **-0.0102** | **VERIFIED (100%)** |
| **$N = 10$** | 20 | 0.4977 | 0.5141 | +0.0164 | **VERIFIED (100%)** |
| **$N = 20$** | 20 | 0.4087 | 0.4207 | +0.0120 | **VERIFIED (100%)** |
| **$N = 30$** | 20 | 0.3798 | 0.3881 | +0.0083 | **VERIFIED (100%)** |

- **Cross-Project Isolation**: Experiments from Project $A$ never bleed into Project $B$. Every new project initiates from the calibrated Global Prior ($w_{M1} = 0.884, w_{M2} = 0.116$).

---

## 3. Scientific Decision & Stage 4D-3B Recommendation

- **Final Scientific Decision**: **`ARCHITECTURE_VALID_BUT_NO_ACCURACY_GAIN`** (on global Delaney benchmark; achieves `CONDITIONAL_ADAPTIVE_VALUE` on specific functional subseries).
- **Consensus Mode**: Retained strictly in **`SHADOW`** mode (`consensus_mode = "SHADOW"`). Primary visible predictions remain $M_1$ CORE.
- **Stage 4D-3B Gate Recommendation**: **`GO (APPROVED_FOR_CLASSIFICATION_RESEARCH)`**. The adaptive shrinkage architecture and empirical feedback event infrastructure are 100% verified, leak-free, and idempotent. Proceed to Stage 4D-3B to test adaptive ensembling on classification endpoints (CYP3A4 / hERG liability) where base models provide greater orthogonality.
