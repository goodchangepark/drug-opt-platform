# Stage 4D-2: Consensus Promotion Criteria & Final Status

## 1. Objective Promotion Gate

To promote an endpoint from **Shadow Consensus** (`consensus_mode = "SHADOW"`) to **Visible Production Consensus** (`consensus_mode = "VISIBLE"` or `"CONSENSUS_ACTIVE"`), the candidate ensemble must satisfy **ALL** of the following quantitative conditions on a qualified external benchmark:

1. **Non-Degradation Criterion**: The consensus ensemble MAE/RMSE (for regression) or Balanced Accuracy/AUROC (for classification) must be equal to or better than the production model ($M_1$).
2. **Statistically Significant Uncertainty Correlation**: Model disagreement standard deviation ($\sigma_w$) must correlate positively ($\rho > +0.25$) with actual prediction error.
3. **Error Diversity Verification**: Pairwise residual correlation between models must satisfy $r < 0.70$ (ensuring genuine error independence).
4. **ARM64 Real-Time Latency**: Total multi-model execution latency per compound must remain $< 1.0\text{ s}$ on Xavier ARM64 CPU.
5. **Memory & Stability Envelope**: Peak memory consumption during batch inference must remain $< 2.0\text{ GB}$ RAM.

---

## 2. Endpoint Promotion Decisions

| Endpoint | Active Qualified Models | Decision | Scientific Justification |
|---|---|---|---|
| **Aqueous Solubility** | `admetica_solubility`<br>`esol_delaney_v1`<br>`rdkit_gbr_solubility_v1` | **`PROMOTION_CANDIDATE`** | Strong error diversity ($r = 0.386$ between $M_1$ and $M_2$), statistically significant uncertainty correlation ($\rho = +0.4699$), and low inference overhead ($< 15\text{ ms}$). Ready for Stage 4D-3 adaptive weighting. |
| **Caco-2 Permeability** | `admetica_caco2`<br>`physchem_caco2_v1` | **`KEEP_SHADOW`** | Model 1 performs adequately ($\text{MAE} = 0.41$). Consensus stabilizes predictions ($\text{MAE} = 0.40$), but small external sample size ($N=34$) warrants retaining shadow mode until larger prospective data arrives. |
| **CYP3A4 Inhibitor** | `admetica_cyp_cyp3a4-inhibitor`<br>`morgan_cyp3a4_inh_v1` | **`PROMOTION_CANDIDATE`** | Complementary sensitivity ($0.65$ vs $0.91$) reduces false-negative risk. Error correlation is very low ($r = 0.207$). Ready for Stage 4D-3. |
| **hERG Liability** | `admetica_safety_herg`<br>`physchem_herg_v1` | **`KEEP_SHADOW`** | Severe class imbalance and low specificity ($0.11$ and $0.004$) on held-out ChEMBL set requires retaining shadow mode until experimental patch-clamp calibration is integrated. |
| **Metabolic Soft Spots** | `sygma_phase1_2`<br>`smartcyp_dft_v1` | **`STAGE_4D2B_PREPARATION_VALIDATED`** | Reciprocal Rank Fusion ($RRF$) architecture validated. Proceed to Stage 4D-2B prospective testing. |

---

## 3. Xavier ARM64 Hardware Benchmark Summary

Performance measured natively on the NVIDIA Jetson Xavier ARM64 8-core Carmel CPU:

- **Cold Model Initialization**: $1,408.28\text{ ms}$
- **Warm Single Compound Multi-Model Execution**: $484.43\text{ ms}$
- **Batch 10 Compounds Multi-Model Execution**: $3,277.50\text{ ms}$ ($\sim 327\text{ ms/compound}$)
- **Deterministic Cache Hit Lookup**: $0.0906\text{ ms}$
- **Peak Process RAM**: $1,427.44\text{ MB}$ ($< 1.5\text{ GB}$)
- **Thermal & Concurrency Policy**: Sequential PyTorch inference with controlled CPU thread pools (`torch.set_num_threads(4)`) guarantees zero thermal throttling and stable execution on Xavier ARM64.
