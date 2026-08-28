# Stage 4D-2: Qualified Multi-Model Pilot Selection & Scientific Audit

## 1. Executive Summary

Stage 4D-2 activates and validates real multi-model prediction across a focused subset of high-value pilot endpoints under strict **Shadow Mode** (`consensus_mode = "SHADOW"`). In accordance with the Stage 4D design freeze, zero modifications were made to the frontend UI, and all visible production endpoints remain 100% identical to pre-4D-2 behavior.

Rather than arbitrarily forcing an equal number of models per endpoint (e.g. forcing 3 models everywhere), Drug-OPT adheres to a strict scientific qualification gate: **every model must possess an identical endpoint contract (same target, definition, species, canonical units, and experimental cutoff semantics) and executable, reproducible inference on Jetson Xavier ARM64**.

---

## 2. Minimum Qualification Gate & Pilot Selection

To be qualified for multi-model pilot execution, every candidate model had to pass 5 rigorous gates:
1. **Endpoint Contract Equivalence**: Must match the authoritative Stage 4D-0 Endpoint Contract (`ENDPOINT_CONTRACTS`).
2. **Species & Matrix Matching**: Strict human / target-specific isolation (cross-species or mismatched assay matrices strictly prohibited).
3. **Legal & Licensing Compliance**: Open-source license (MIT, Apache-2.0, BSD, or public scientific domain) with clear attribution.
4. **Architectural & Error Diversity**: Independent model lineage or mathematical representation (not just trivial checkpoint clones).
5. **Deterministic Xavier ARM64 CPU Execution**: Must run natively within memory constraints (< 2 GB peak RAM) without GPU lockup.

### Summary of Selected Pilot Models

| Endpoint | Canonical ID | Output Type | Selected Qualified Models ($K$) | Model IDs |
|---|---|---|---|---|
| **Aqueous Solubility** | `EP_PHYS_SOLUBILITY` | Continuous ($\log_{10}(\text{mol/L})$) | **3 Models** | `admetica_solubility`<br>`esol_delaney_v1`<br>`rdkit_gbr_solubility_v1` |
| **Caco-2 Permeability** | `EP_ABS_CACO2` | Continuous ($\log_{10}(\text{cm/s})$) | **2 Models** | `admetica_caco2`<br>`physchem_caco2_v1` |
| **CYP3A4 Inhibitor** | `EP_MET_CYP3A4_INH` | Binary Classification ($\text{AC}_{50} \le 10\,\mu\text{M}$) | **2 Models** | `admetica_cyp_cyp3a4-inhibitor`<br>`morgan_cyp3a4_inh_v1` |
| **hERG Liability** | `EP_TOX_HERG` | Binary Classification ($\text{IC}_{50} \le 10\,\mu\text{M}$) | **2 Models** | `admetica_safety_herg`<br>`physchem_herg_v1` |
| **Site of Metabolism (SoM)** | `EP_MET_SOM` | Atom Ranking (Rank Fusion) | **2 Models (4D-2B Prep)** | `sygma_phase1_2`<br>`smartcyp_dft_v1` |

---

## 3. Endpoint-by-Endpoint Model Audit & Architecture

### 3.1. Aqueous Solubility (`EP_PHYS_SOLUBILITY`)
- **Scientific Contract**: Thermodynamic/intrinsic aqueous solubility at pH ~7.0, $25^\circ\text{C}$ in $\log_{10}(\text{mol/L})$.
- **Model 1 (`admetica_solubility`)**: Directed Message Passing Neural Network (D-MPNN) trained on AqSolDB ($N=9,982$). Learns 3D-like atom/bond message vectors.
- **Model 2 (`esol_delaney_v1`)**: Delaney multilinear physical descriptor model ($R^2=0.88$ on Delaney). Calculates $\text{LogS}$ from cLogP, MW, Rotatable Bonds, and Aromatic Proportion. Represents a completely independent physical chemistry baseline.
- **Model 3 (`rdkit_gbr_solubility_v1`)**: 2D topological + ECFP4 fingerprint gradient boosted regressor calibrated on AqSolDB. Captures non-linear topological surface contributions.

### 3.2. Caco-2 Permeability (`EP_ABS_CACO2`)
- **Scientific Contract**: Apparent permeability coefficient ($P_{\text{app}}, A\to B$) across human colon carcinoma cell monolayers in $\log_{10}(\text{cm/s})$.
- **Model 1 (`admetica_caco2`)**: Chemprop D-MPNN trained on Wang et al. dataset ($N=1,272$).
- **Model 2 (`physchem_caco2_v1`)**: Mechanistic polar surface area and membrane partition model incorporating cLogP, TPSA, MW, HBD, and physiological net charge at pH 7.4.

### 3.3. CYP3A4 Inhibitor (`EP_MET_CYP3A4_INH`)
- **Scientific Contract**: Human CYP3A4 enzymatic inhibition at $10\,\mu\text{M}$ cutoff (PubChem AID 1851 protocol). Positive: Inhibitor; Negative: Non-inhibitor.
- **Model 1 (`admetica_cyp_cyp3a4-inhibitor`)**: Chemprop D-MPNN trained on PubChem AID 1851 ($N=12,320$).
- **Model 2 (`morgan_cyp3a4_inh_v1`)**: Morgan radius 2 + physicochemical nitrogen heterocycle pharmacophore classifier.

### 3.4. hERG Cardiotoxicity Liability (`EP_TOX_HERG`)
- **Scientific Contract**: Human Ether-à-go-go Related Gene ($K_v11.1$) channel inhibition at $10\,\mu\text{M}$ cutoff.
- **Model 1 (`admetica_safety_herg`)**: Chemprop D-MPNN trained on Wang et al. literature compilation ($N=22,248$).
- **Model 2 (`physchem_herg_v1`)**: Basic center amine pharmacophore + lipophilicity logistic classifier.

---

## 4. Deferral and Rejection Decisions

The following candidate models were evaluated and explicitly deferred or rejected from the active consensus pool:

1. **ADMET-AI Regression Checkpoints for Solubility/Caco-2**:
   - *Status*: Deferred to Stage 4D-3.
   - *Reasoning*: Multi-task Chemprop regression models from ADMET-AI use slightly divergent target normalization parameters ($\ln(P_{\text{app}})$ vs $\log_{10}(P_{\text{app}})$), requiring explicit calibration wrappers before consensus inclusion.
2. **P-gp Substrate Classifier**:
   - *Status*: Deferred (Single installed model retained).
   - *Reasoning*: External public transport data sets lack standardized assay conditions (MDCK vs Caco-2, apical vs basolateral).
3. **Cross-Species Clearance Scaling**:
   - *Status*: Strictly Rejected.
   - *Reasoning*: Attempting to ensemble RLM (rat) and HLM (human) clearance models directly violates the species isolation gate.
