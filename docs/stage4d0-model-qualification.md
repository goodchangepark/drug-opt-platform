# Stage 4D-0: Comprehensive Scientific Model Qualification Audit

**Platform Version:** `0.6.3-stage5b4-ui`  
**Execution Environment:** NVIDIA Jetson AGX Xavier (32GB ARM64, Ubuntu 22.04 LTS, Python 3.11.15)  
**Milestone Tag:** `stage4d0-model-qualification-complete`

---

## 1. Executive Summary & Qualification Framework

Stage 4D-0 establishes the scientific foundation required for future multi-model adaptive ensemble prediction in Drug-OPT. This stage executes a strict qualification audit across all currently installed models and open-source candidate model families.

### Qualification Assessment Framework (Criteria A–H):
1. **[Criterion A] Scientific Compatibility:** Exact alignment with the Drug-OPT authoritative endpoint contract (endpoint semantics, species, matrix, and canonical unit).
2. **[Criterion B] Reproducibility:** Availability of deterministic open-source model code, fixed random seeds, and reproducible training/inference pipelines.
3. **[Criterion C] Data Provenance:** Fully documented training datasets with public origin and explicit chemical structure curation.
4. **[Criterion D] Validation Quality:** Rigorous held-out cross-validation and, where accessible, prospective non-overlapping benchmark validation.
5. **[Criterion E] Leakage Risk:** Explicit tracking of training-set overlap against internal and prospective evaluation benchmarks.
6. **[Criterion F] Model Diversity:** Architectural and representational independence (e.g. Graph Neural Networks vs Gradient Boosted Decision Trees vs Mechanistic Physics).
7. **[Criterion G] Legal & Licensing:** Unambiguous open-source licenses permitting redistribution and local deployment without commercial or academic restrictions.
8. **[Criterion H] Technical Compatibility:** Verified, error-free CPU execution on NVIDIA Jetson AGX Xavier ARM64 under Python 3.11 with low latency (<50ms) and reasonable memory (<500MB).

---

## 2. Qualification Status Definitions

* **`QUALIFIED`**: Meets all scientific, licensing, validation, and ARM64 requirements. Approved as a primary ensemble candidate.
* **`QUALIFIED_WITH_LIMITATIONS`**: Scientifically compatible and technically functional on ARM64, but has known constraints (e.g., shared training data lineage or moderate class imbalance) requiring variance penalties or rank fusion.
* **`RESEARCH_ONLY`**: Suitable for sandbox R&D exploration but rejected for company production due to licensing restrictions (e.g. academic-only terms) or unstable dependencies.
* **`REJECTED`**: Scientifically incompatible (e.g., binary TDC HLM for continuous clearance), proprietary, or reliant on unsupported native libraries on Xavier ARM64.
* **`UNAVAILABLE`**: Explicit placeholder in the Drug-OPT Model Registry to prevent fabricated or out-of-domain predictions.

---

## 3. Current Installed Model Audit (18 Active Models)

The Drug-OPT platform currently operates 18 active machine-learning checkpoints and 2 deterministic rule/physics engines:

```
[Installed Model Inventory Summary]
├── Physicochemical:
│   ├── Solubility: Admetica Chemprop Solubility (QUALIFIED, MIT)
│   └── Ionization Engine: Henderson-Hasselbalch + Micro-state pKa (QUALIFIED, Internal)
├── Absorption & Distribution:
│   ├── Caco-2 Permeability: Admetica Chemprop Caco-2 (QUALIFIED, MIT)
│   └── Human PPB: Admetica Chemprop Human PPB (QUALIFIED, MIT)
├── Hepatic Microsomal Clearance (Multi-Species):
│   ├── HLM (Human): OpenADMET CheMeleon MPNN (QUALIFIED, Apache-2.0)
│   ├── RLM (Rat): OpenADMET CheMeleon MPNN (QUALIFIED, Apache-2.0)
│   └── MLM (Mouse): OpenADMET CheMeleon MPNN (QUALIFIED, Apache-2.0)
├── CYP450 Panel (Inhibitors & Substrates):
│   ├── CYP1A2 Inhibitor: Admetica Chemprop (QUALIFIED, MIT)
│   ├── CYP2C9 Inhibitor: Admetica Chemprop (QUALIFIED, MIT)
│   ├── CYP2C19 Inhibitor: Admetica Chemprop (QUALIFIED, MIT)
│   ├── CYP2D6 Inhibitor: Admetica Chemprop (QUALIFIED, MIT)
│   ├── CYP3A4 Inhibitor: Admetica Chemprop (QUALIFIED, MIT)
│   ├── CYP2C9 Substrate: Admetica Chemprop (QUALIFIED, MIT)
│   ├── CYP2D6 Substrate: Admetica Chemprop (QUALIFIED, MIT)
│   └── CYP3A4 Substrate: Admetica Chemprop (QUALIFIED, MIT)
├── Transporters:
│   └── P-gp Inhibitor: Admetica Chemprop (QUALIFIED, MIT)
├── Safety & Toxicology:
│   ├── hERG Liability: Admetica Chemprop (QUALIFIED_WITH_LIMITATIONS, MIT)
│   ├── Ames Mutagenicity: ADMET-AI Chemprop 5-Ensemble (QUALIFIED, MIT)
│   └── DILI Clinical Liability: ADMET-AI Chemprop 5-Ensemble (QUALIFIED, MIT)
└── Metabolism Rule Engine:
    └── SyGMa: Phase I / Phase II SMARTS Metabolite Generator (QUALIFIED, GPL-3.0)
```

Every installed checkpoint has been verified on Xavier ARM64 with SHA-256 integrity checks in `validation/stage4d0_current_model_inventory.json`.

---

## 4. Pilot Endpoint Recommendations for Future Ensembles

To ensure maximum scientific impact without introducing spurious complexity, Stage 4D-0 qualifies four distinct pilot endpoints for initial multi-model ensemble integration:

### 1. Aqueous Solubility ($\text{LogS}$)
* **Candidate 1:** `Admetica Chemprop Solubility` (D-MPNN Molecular Graph) — *Installed Primary*
* **Candidate 2:** `ADMET-AI Chemprop Ensemble` (5-Fold D-MPNN + 200 RDKit Descriptors) — *Qualified (Dataset Overlap)*
* **Candidate 3:** `AqSolPred GBDT / CatBoost` (2D Descriptors + Morgan FP) — *Qualified (High Architectural Diversity)*
* **Scientific Justification:** High public data availability (AqSolDB), low measurement noise compared to biological assays, and genuine architectural diversity between Graph Neural Networks and Tree-based models.

### 2. Caco-2 Membrane Permeability ($\text{LogPapp}$)
* **Candidate 1:** `Admetica Chemprop Caco-2` (D-MPNN) — *Installed Primary*
* **Candidate 2:** `ADMET-AI Chemprop Ensemble Caco-2` (5-Fold D-MPNN) — *Qualified (Dataset Overlap)*
* **Candidate 3:** `TDC Benchmark CatBoost Caco-2` (Tabular Tree) — *Qualified (High Diversity)*
* **Scientific Justification:** Critical driver for IVIVE absorption ($F_a$) and gastrointestinal transit modeling.

### 3. CYP3A4 Inhibition Probability
* **Candidate 1:** `Admetica Chemprop CYP3A4 inhibitor` (D-MPNN) — *Installed Primary*
* **Candidate 2:** `ADMET-AI Chemprop Ensemble CYP3A4` (5-Fold D-MPNN) — *Qualified (Dataset Overlap)*
* **Candidate 3:** `CatBoost GBDT CYP3A4 Classifier` (ECFP4 + 2D Descriptors) — *Qualified (High Diversity)*
* **Scientific Justification:** Major drug-drug interaction (DDI) clearance bottleneck accounting for >50% of marketed drug oxidative metabolism.

### 4. Metabolic Site-of-Metabolism (SoM)
* **Candidate 1:** `SyGMa SMARTS Reaction Engine` (Rule-Based Substructure Matching) — *Installed Primary*
* **Candidate 2:** `SMARTCyp Quantum-Chemical Model` (DFT Activation Energies $E_{act}$) — *Qualified (Mechanistic)*
* **Scientific Justification:** Combining rule-based reaction libraries with physical quantum mechanical transition-state energetics via Reciprocal Rank Fusion (RRF) produces superior soft-spot identification.

---

## 5. Special Scientific Handling Rules

### 5.1 Single Qualified Model for HLM Clearance
* **Policy:** For HLM, RLM, and MLM scaled microsomal clearance ($\log_{10}(\text{mL/min/kg})$), only **OpenADMET CheMeleon** currently qualifies.
* **No Artificial Model Inflation:** We strictly refuse to force arbitrary model counts by forcing incompatible models (such as TDC binary classification models) into continuous clearance ensembles. Operating with a single well-qualified model is scientifically sound; forcing pseudo-diversity is prohibited.

### 5.2 Metabolism Site-of-Metabolism Fusion
* **Policy:** SoM outputs are rank-ordered atom indices, not scalar quantities.
* **Aggregation Method:** Multi-model SoM integration must use **Reciprocal Rank Fusion (RRF)** or **Borda Count**, NEVER arithmetic averaging.

### 5.3 Pharmacokinetics & IVIVE Hierarchy
* **Policy:** In PK parameter determination, mechanistic and experimental evidence take absolute precedence:
  $$\text{Matched In Vivo NCA} \succ \text{Physiological IVIVE (Well-Stirred)} \succ \text{Allometric Interspecies Scaling} \succ \text{Direct QSAR ML}$$
