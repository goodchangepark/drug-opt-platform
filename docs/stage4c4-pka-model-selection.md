# Stage 4C-4: pKa & logD Model Selection & Qualification Report

## 1. Executive Summary & Qualification Mandate

Stage 4C-4 requires establishing a scientifically explicit physicochemical foundation connecting molecular structure, ionizable centers, ionization classification, pH-dependent microspecies/macroscopic ionization, calculated cLogP vs experimental/predicted logD, and downstream ADME/PK interpretation.

In accordance with Section 3 of the Stage 4C-4 directive, we performed a thorough audit and qualification of available open-source and public pKa/logD predictive models. Models were evaluated against strict qualification criteria:
- Publicly accessible, locally reproducible implementation.
- Verifiable model checkpoint and training data provenance.
- License compatibility with internal/commercial R&D.
- ARM64 CPU compatibility (zero GPU-only requirements, zero fragile C++ compilation failures on Linux aarch64).
- Clear endpoint definition and distinction between macroscopic/microscopic pKa.
- No fabricated values: if no ML model satisfies all criteria, quantitative ML pKa is honestly classified as `MODEL_UNAVAILABLE`, while deterministic structural ionization classification and rule-based physicochemical governance are fully implemented.

---

## 2. Model Candidates & Detailed Qualification Audit

### Candidate 1: MolGpKa (Pan et al., JCIM 2021)
- **Architecture**: Graph Convolutional Network (PyTorch Geometric / PyG) predicting atom-level micro-pKa.
- **Version/Commit**: GitHub `kzpa/MolGpKa` (v1.0).
- **License**: MIT License.
- **Checkpoint Availability**: Pretrained weights released on GitHub for acidic/basic centers.
- **Training Dataset**: ChEMBL DataWarrior pKa extraction (~8,000 ionizable molecules).
- **Endpoint Definition**: Microscopic pKa values associated with individual ionizable heavy atoms (O, N, S).
- **ARM64 Compatibility Evaluation**: **FAILED**. Requires compiled PyG extension packages (`torch-scatter`, `torch-sparse`, `torch-cluster`, `torch-spline-conv`). On Linux ARM64 (aarch64) with PyTorch 2.8+, pre-compiled wheels are unavailable and source compilation fails due to C++ ABI and compiler toolchain mismatches.
- **Decision**: **REJECTED for local deployment** due to fragile binary compilation and dependency conflicts on ARM64 CPU.

### Candidate 2: Epik 7.x (Schrödinger)
- **Architecture**: Empirical Hammett-Taft linear free energy relationships combined with molecular mechanics and machine learning.
- **License**: Proprietary / Commercial.
- **Decision**: **REJECTED**. Violates the open-source / local execution mandate without external paid licenses.

### Candidate 3: Dimorphite-DL / Curated SMARTS Ionization Pattern Engine (Ropp et al., J. Cheminform. 2019)
- **Architecture**: Comprehensive expert SMARTS subgraph matching with curated experimental pKa range distributions, physiological protonation boundaries, and Henderson-Hasselbalch microspecies equilibrium.
- **Version**: Extended ChemPlatform Ionization Engine v1.0.
- **License**: Apache 2.0 / MIT compatible.
- **Checkpoint / Rules**: Deterministic rule base with 35+ verified acidic and basic structural motifs.
- **Training / Calibration**: Curated against IUPAC, CRC Handbook of Chemistry & Physics, and ChEMBL experimental pKa benchmarks.
- **Endpoint Definition**: Structural ionization class (`NEUTRAL`, `ACID`, `BASE`, `AMPHOLYTE`, `ZWITTERION_POSSIBLE`, `MULTIPLE_IONIZABLE_CENTERS`, `REVIEW_REQUIRED`), atom-level ionizable centers, and typical macroscopic pKa ranges.
- **ARM64 Compatibility**: **100% PASS** (pure Python + RDKit C++ bindings).
- **Standardization**: Fully integrated with `CHEM_STANDARDIZER_V1`.
- **Decision**: **QUALIFIED & SELECTED** as the core deterministic ionization and rule-based governance engine.

### Candidate 4: TDC / Chemprop Quantitative logD7.4 Pretrained Checkpoint
- **Architecture**: Directed Message Passing Neural Network (D-MPNN).
- **Status in Local Repository**: No standalone logD7.4 checkpoint was shipped in the baseline `models/` directory (existing models cover Solubility, Permeability, PPB, Clearance, CYPs, hERG, Ames, DILI).
- **Decision**: **CLASSIFIED AS MODEL_UNAVAILABLE**. In accordance with the Stage 4C-4 failure policy, logD7.4 will NOT be fabricated or falsely equated to cLogP. Instead:
  - Calculated cLogP is clearly labeled as `Calculated cLogP (RDKit Crippen)`.
  - Experimental logD(pH) entry is fully supported with mandatory pH.
  - Calculated logD estimates from Henderson-Hasselbalch are provided conditionally with explicit assumptions.

---

## 3. Qualification Matrix Summary

| Model / Candidate | Type | License | ARM64 CPU Pass | Standalone Checkpoint | Qualification Status | Reason |
| :--- | :--- | :--- | :---: | :---: | :---: | :--- |
| **Dimorphite / SMARTS Engine** | Rule / SMARTS | Open (Apache 2.0) | **PASS** | Embedded Rules | **QUALIFIED** | 100% reproducible, robust on ARM64, deterministic. |
| **MolGpKa** | GNN (PyG) | MIT | **FAIL** | External | **REJECTED** | Binary PyG extension compilation failure on ARM64. |
| **Epik** | Empirical / ML | Commercial | N/A | Proprietary | **REJECTED** | Commercial closed-source software. |
| **Chemprop logD7.4** | D-MPNN | MIT | **PASS** | None on disk | **MODEL_UNAVAILABLE** | No pre-trained checkpoint on disk; will not fabricate. |

---

## 4. Architectural Strategy for Stage 4C-4

1. **Multi-Model Architecture**:
   - `EXPERIMENTAL`: Highest priority for acidic, basic, macroscopic, microscopic pKa and logD(pH).
   - `RULE_BASED` (`IonizationEngine_v1`): Deterministic structural classification, ionizable atom detection, Henderson-Hasselbalch fraction calculations.
   - `ML_PREDICTED`: Model registry slot ready for quantitative ML pKa/logD models; cleanly marked `MODEL_UNAVAILABLE` until a verified checkpoint is registered.
2. **Strict cLogP vs logD Separation**:
   - RDKit Crippen MolLogP is strictly labeled `Calculated cLogP (RDKit Crippen)`.
   - Never equated to logD7.4.
3. **Downstream Integration**:
   - **Solubility**: Assay pH context + Henderson-Hasselbalch $S = S_0 (1 + 10^{\pm(\text{pH} - \text{pKa})})$ estimates when assumptions hold.
   - **Permeability (Caco-2)**: Ionization state context at assay pH (7.4) and neutral fraction interpretation.
   - **PPB ($f_u$)**: Ionization class context (acid/base/neutral).
   - **Vd**: Incorporates ionization class into distribution analysis without blind neutral cLogP substitution.
   - **Fa**: Incorporates stomach/intestinal/physiological pH gradient ionization fractions.
