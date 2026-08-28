# Stage 4D-0: Authoritative Scientific Endpoint Contracts & Semantic Taxonomy

**Platform Version:** `0.6.3-stage5b4-ui`  
**Stage:** `4D-0 (Model Qualification & Ensemble Foundation)`  
**Scope:** Strict physical/biological definitions, units, transformations, species isolation, experimental comparison rules, and ensemble compatibility gates.

---

## 1. Executive Summary & Foundational Principles

In high-throughput ADMET modeling and translational pharmacokinetics, model ensembles frequently fail not from algorithmic defects, but from **semantic conflation**—the erroneous aggregation of models trained on fundamentally different physical endpoints under a common colloquial name.

The Drug-OPT platform establishes a zero-tolerance policy against semantic conflation based on four immutable principles:
1. **Scientific Semantics First:** Two models are NEVER ensemble-compatible merely because they share a string label (e.g., "solubility", "clearance", "hERG").
2. **Species & Matrix Isolation:** Interspecies parameter substitution (e.g., using mouse MLM to impute human HLM) and biological matrix substitution (e.g., isolated albumin for whole plasma) are strictly prohibited in prediction ensembling.
3. **Unit & Scale Invariance:** All models participating in an ensemble must emit values in, or be strictly convertible to, the endpoint's authoritative **Canonical Unit** and log-transformation scale.
4. **Directionality & Role Separation:** Inhibitor and Substrate classifications are physiologically distinct and must never be merged or averaged.

---

## 2. Authoritative Endpoint Taxonomy

The table below summarizes the authoritative contracts governing all primary endpoints in Drug-OPT:

| Endpoint | Category | Species | Canonical Unit | Transformation | Directionality | Output Type |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Solubility** | Physicochemical | In Vitro | $\log_{10}(\text{mol/L})$ | Identity | Higher is better | Continuous Regression |
| **Permeability** | Absorption | Human Caco-2 | $\log_{10}(\text{cm/s})$ | Identity | Higher is better | Continuous Regression |
| **Plasma Protein Binding** | Distribution | Human Plasma | $\%$ bound | Identity ($f_u = \frac{100 - \text{bound}}{100}$) | Neutral | Continuous Regression |
| **HLM Intrinsic Clearance** | Metabolism | Human | $\log_{10}(\text{mL/min/kg})$ | Identity (Scaled In Vivo) | Lower is better | Continuous Regression |
| **RLM Intrinsic Clearance** | Metabolism | Rat | $\log_{10}(\text{mL/min/kg})$ | Identity (Scaled In Vivo) | Lower is better | Continuous Regression |
| **MLM Intrinsic Clearance** | Metabolism | Mouse | $\log_{10}(\text{mL/min/kg})$ | Identity (Scaled In Vivo) | Lower is better | Continuous Regression |
| **CYP1A2 Inhibitor** | CYP Panel | Human Recombinant | Probability $[0, 1]$ | $P(\text{AC}_{50} \le 10\,\mu\text{M})$ | Lower is better | Binary Classification |
| **CYP2C9 Inhibitor** | CYP Panel | Human Recombinant | Probability $[0, 1]$ | $P(\text{AC}_{50} \le 10\,\mu\text{M})$ | Lower is better | Binary Classification |
| **CYP2C19 Inhibitor** | CYP Panel | Human Recombinant | Probability $[0, 1]$ | $P(\text{AC}_{50} \le 10\,\mu\text{M})$ | Lower is better | Binary Classification |
| **CYP2D6 Inhibitor** | CYP Panel | Human Recombinant | Probability $[0, 1]$ | $P(\text{AC}_{50} \le 10\,\mu\text{M})$ | Lower is better | Binary Classification |
| **CYP3A4 Inhibitor** | CYP Panel | Human Recombinant | Probability $[0, 1]$ | $P(\text{AC}_{50} \le 10\,\mu\text{M})$ | Lower is better | Binary Classification |
| **CYP2C9 Substrate** | CYP Panel | Human Recombinant | Probability $[0, 1]$ | $P(\text{Turnover Substrate})$ | Neutral | Binary Classification |
| **CYP2D6 Substrate** | CYP Panel | Human Recombinant | Probability $[0, 1]$ | $P(\text{Turnover Substrate})$ | Neutral | Binary Classification |
| **CYP3A4 Substrate** | CYP Panel | Human Recombinant | Probability $[0, 1]$ | $P(\text{Turnover Substrate})$ | Neutral | Binary Classification |
| **P-gp Inhibitor** | Transporter | Human P-gp/ABCB1 | Probability $[0, 1]$ | $P(\text{Inhibitor at } \le 15\,\mu\text{M})$ | Lower is better | Binary Classification |
| **hERG Liability** | Safety | Human $K_v11.1$ | Probability $[0, 1]$ | $P(\text{Blocker at } \le 10\,\mu\text{M})$ | Lower is better | Binary Classification |
| **Ames Mutagenicity** | Safety | *S. typhimurium* | Probability $[0, 1]$ | $P(\text{Bacterial Mutation})$ | Lower is better | Binary Classification |
| **DILI Clinical Concern** | Safety | Human Clinical | Probability $[0, 1]$ | $P(\text{FDA-LTKB Concern})$ | Lower is better | Binary Classification |
| **Metabolic Soft Spots** | Metabolism | Human | Atom Ranking | Rank-ordered atom indices | Neutral | Atom-Level Ranking |
| **PK Systemic Clearance** | Pharmacokinetics | In Vivo Species | $\text{mL/min/kg}$ | IV NCA $CL$ | Lower is better | Mechanistic / NCA |
| **PK Volume ($V_{ss}$)** | Pharmacokinetics | In Vivo Species | $\text{L/kg}$ | IV NCA $V_{ss}$ | Neutral | Mechanistic / NCA |
| **PK Bioavailability ($F$)**| Pharmacokinetics | In Vivo Species | $\%$ | Matched $\frac{\text{AUC}_{po} \cdot \text{Dose}_{iv}}{\text{AUC}_{iv} \cdot \text{Dose}_{po}} \times 100$ | Higher is better | In Vivo Experimental |

---

## 3. Detailed Endpoint Definitions & Boundary Contracts

### 3.1 Aqueous Solubility (LogS)
* **Authoritative Definition:** Logarithmic molar aqueous solubility, $\text{LogS} = \log_{10}(S\,[\text{mol/L}])$, across general aqueous media at $25^\circ\text{C}$.
* **Boundary Exclusions:**
  * Must NOT be equated to intrinsic thermodynamic solubility ($S_0$, uncharged species) without explicit Henderson-Hasselbalch ionization correction.
  * Kinetic solubility in high DMSO (>2%) cannot be ensembled directly with thermodynamic shake-flask solubility without variance penalties.
* **Conversion from Mass Solubility:**
  $$\text{LogS} = \log_{10}\left(\frac{S_{\text{mg/mL}} \times 1000}{\text{MW}}\right)$$

### 3.2 Caco-2 Permeability (LogPapp)
* **Authoritative Definition:** Apparent apical-to-basolateral ($A \rightarrow B$) membrane permeability in human colon adenocarcinoma (Caco-2) cell monolayers, $\text{LogPapp} = \log_{10}(P_{\text{app}}\,[\text{cm/s}])$.
* **Boundary Exclusions:**
  * PAMPA (Parallel Artificial Membrane Permeability Assay), MDCK, and BBB PAMPA are distinct biological matrices and MUST NOT be ensembled as Caco-2.
  * Efflux Ratio ($\text{ER} = P_{\text{app}, B\rightarrow A} / P_{\text{app}, A\rightarrow B}$) is a transporter metric, never an apparent permeability rate.

### 3.3 Plasma Protein Binding (PPB) & Fraction Unbound ($f_u$)
* **Authoritative Definition:** Percentage of compound bound to plasma proteins in pooled human plasma at equilibrium dialysis, bounded in $[0.0, 100.0]\%$.
* **Derived Unbound Fraction ($f_{u,p}$):**
  $$f_{u,p} = \max\left(0.0001, \frac{100.0 - \text{Percent Bound}}{100.0}\right)$$
* **Boundary Exclusions:**
  * Animal plasma (rat, mouse, dog) has distinct albumin/AAG binding affinities; cross-species PPB ensembling is prohibited.
  * Isolated Bovine Serum Albumin (BSA) or Alpha-1-Acid Glycoprotein (AAG) binding must not substitute for whole plasma PPB.

### 3.4 Hepatic Microsomal Clearance (HLM, RLM, MLM)
* **Authoritative Definition:** Species-specific hepatic microsomal intrinsic clearance scaled to in vivo body weight clearance and expressed in $\log_{10}(\text{mL/min/kg})$.
* **Scaling Equations from In Vitro Substrate Depletion ($Cl_{\text{int,micr}}\ [\mu\text{L/min/mg}]$):**
  $$\text{Scaled } Cl_{\text{int}}\ [\text{mL/min/kg}] = \frac{Cl_{\text{int,micr}}\ [\mu\text{L/min/mg}] \times \text{MPPGL} \times \text{Liver Weight}}{1000}$$
  * Human: $\text{MPPGL} = 45.0\,\text{mg/g}$, $\text{Liver Weight} = 25.7\,\text{g/kg}$
  * Rat: $\text{MPPGL} = 45.0\,\text{mg/g}$, $\text{Liver Weight} = 40.0\,\text{g/kg}$
  * Mouse: $\text{MPPGL} = 45.0\,\text{mg/g}$, $\text{Liver Weight} = 87.5\,\text{g/kg}$
* **Boundary Exclusions:**
  * In vitro half-life ($t_{1/2}\ [\text{min}]$) and raw $\mu\text{L/min/mg}$ MUST be scaled before ensembling.
  * Cross-species clearance substitution (e.g., using rat RLM for human HLM) is strictly rejected.

### 3.5 CYP450 Panel: Inhibitors vs. Substrates
* **Authoritative Definition:** Binary classification class probabilities ($P \in [0.0, 1.0]$).
  * **Inhibitors:** Functional inhibition in recombinant human CYP dealkylation assays with active threshold $\text{AC}_{50} \le 10\,\mu\text{M}$ (PubChem AID 1851 protocol).
  * **Substrates:** Qualitative enzymatic metabolic turnover substrate.
* **Boundary Exclusions:**
  * Inhibitor and Substrate models are physiologically orthogonal and MUST NEVER be combined.
  * Probability scores must not be converted to quantitative $\text{IC}_{50}$ or $K_i$ values without explicit dose-response regression models.

### 3.6 Safety Panel: hERG, Ames, DILI
* **hERG Cardiac Liability:** Blocker liability probability ($P(\text{blocker})$) thresholded at $\text{IC}_{50} \le 10\,\mu\text{M}$ ($p\text{IC}_{50} \ge 5.0$). Must NOT be interpreted as clinical QT prolongation or patch-clamp current kinetics.
* **Ames Mutagenicity:** Bacterial reverse mutation probability in *S. typhimurium* strains $\pm\text{S9}$ activation (OECD 471). Must NOT be conflated with mammalian micronucleus tests.
* **DILI Clinical Liability:** Clinical drug-induced liver injury association from the FDA-NCTR LTKB database. Must NOT be ensembled with in vitro hepatocyte cytotoxicity or isolated mitochondrial toxicity.

---

## 4. Deterministic Compatibility Verification Gate

The `backend.endpoint_contracts.check_ensemble_compatibility` function enforces an automated gate before any candidate model may join an ensemble:

```python
from backend.endpoint_contracts import check_ensemble_compatibility, ENDPOINT_CONTRACTS

contract_hlm = ENDPOINT_CONTRACTS["HLM intrinsic clearance"]
contract_caco2 = ENDPOINT_CONTRACTS["Permeability"]

is_valid, reason = check_ensemble_compatibility(contract_hlm, contract_caco2)
# Output: (False, "Incompatible endpoint IDs: 'hlm_intrinsic_clearance_scaled_log10' vs 'permeability_caco2_logpapp'")
```

All future multi-model ensemble aggregations must strictly pass this compatibility gate.
