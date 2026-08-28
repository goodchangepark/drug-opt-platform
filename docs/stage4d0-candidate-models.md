# Stage 4D-0: Detailed Candidate Model Catalog & Technical Analysis

**Platform Version:** `0.6.3-stage5b4-ui`  
**Hardware Target:** NVIDIA Jetson AGX Xavier (32GB ARM64)

---

## 1. Candidate Catalog by Endpoint

### 1. Aqueous Solubility ($\text{LogS}$)
* **Candidate 1: Admetica Chemprop Solubility**
  * *Architecture:* Chemprop Directed Message Passing Neural Network (D-MPNN, $d_h=2048$, depth 3).
  * *Features:* 2D Molecular Graph + 6 RDKit molecular descriptors (MW, cLogP, TPSA, HBD, HBA, RotB).
  * *Training Set:* AqSolDB curated aqueous solubility (9,982 compounds).
  * *Performance:* $\text{MAE} = 0.714$, $\text{RMSE} = 1.089$, $R^2 = 0.788$, Spearman $\rho = 0.897$.
  * *Verdict:* **`QUALIFIED`** (Installed Primary Model).
* **Candidate 2: ADMET-AI Chemprop 5-Fold Ensemble Solubility**
  * *Architecture:* 5-model D-MPNN ensemble with 200 RDKit 2D features.
  * *Training Set:* TDC AqSolDB (9,982 compounds).
  * *Performance:* $\text{MAE} = 0.767$, $\text{RMSE} = 1.052$.
  * *Verdict:* **`QUALIFIED_WITH_LIMITATIONS`** (High dataset overlap with Admetica; downweighted correlation in ensemble).
* **Candidate 3: AqSolPred / CatBoost Tree Regressor**
  * *Architecture:* Gradient Boosted Decision Trees (CatBoost / LightGBM).
  * *Features:* 200 RDKit 2D Descriptors + 1024-bit Morgan Fingerprints (ECFP4).
  * *Training Set:* AqSolDB (9,982 compounds).
  * *Performance:* $\text{MAE} = 0.82$, $R^2 = 0.75$.
  * *Verdict:* **`QUALIFIED`** (High architectural diversity; orthogonal error profile).

---

### 2. Caco-2 Apparent Permeability ($\text{LogPapp}$)
* **Candidate 1: Admetica Chemprop Caco-2**
  * *Architecture:* D-MPNN (2048-dim message passing).
  * *Training Set:* Wang et al. Caco-2 compiled set (910 compounds).
  * *Performance:* $\text{MAE} = 0.317$, $R^2 = 0.701$, External 34-compound benchmark $\text{MAE} = 0.412$.
  * *Verdict:* **`QUALIFIED`** (Installed Primary Model).
* **Candidate 2: ADMET-AI Chemprop 5-Fold Ensemble Caco-2**
  * *Architecture:* 5-model D-MPNN + 200 RDKit features.
  * *Training Set:* TDC Wang Caco-2 (910 compounds).
  * *Performance:* $\text{MAE} = 0.341$, $R^2 = 0.680$.
  * *Verdict:* **`QUALIFIED_WITH_LIMITATIONS`** (Shared training set).
* **Candidate 3: TDC Benchmark CatBoost Caco-2**
  * *Architecture:* CatBoost Regressor on ECFP4 + Descriptors.
  * *Performance:* $\text{MAE} = 0.380$, $R^2 = 0.640$.
  * *Verdict:* **`QUALIFIED`** (High architectural diversity).

---

### 3. Human Plasma Protein Binding (PPB)
* **Candidate 1: Admetica Chemprop Human PPB**
  * *Architecture:* D-MPNN (2048-dim).
  * *Training Set:* AstraZeneca / ChEMBL 3301361 (2,790 compounds).
  * *Performance:* $\text{MAE} = 6.92\%$ bound, $R^2 = 0.609$, Biogen held-out benchmark $\text{MAE} = 14.6\%$.
  * *Verdict:* **`QUALIFIED`** (Installed Primary Model).
* **Candidate 2: ADMET-AI Chemprop PPBR**
  * *Architecture:* 5-model D-MPNN.
  * *Training Set:* TDC PPBR (1,614 compounds subset).
  * *Performance:* $\text{MAE} = 8.50\%$ bound.
  * *Verdict:* **`QUALIFIED_WITH_LIMITATIONS`** (Partial dataset overlap).
* **Candidate 3: Lombardo Mechanistic Physicochemical Model**
  * *Architecture:* Deterministic ionization and lipophilicity equations for albumin/AAG binding.
  * *Performance:* Average fold error $1.8\times$ across 1,058 compounds.
  * *Verdict:* **`QUALIFIED_WITH_LIMITATIONS`** (Valuable physical guardrail for lipophilic bases/acids).

---

### 4. Hepatic Microsomal Clearance (HLM, RLM, MLM)
* **Candidate 1: OpenADMET CheMeleon Microsomal Clearance**
  * *Architecture:* CheMeleon MPNN Multi-Task Regressor ($d_h=2048$, depth 3).
  * *Training Set:* ChEMBL 35, ASAP-Polaris, ExpansionRx (5,086 HLM, 670 RLM, 5,086 MLM).
  * *Output:* Species-specific scaled in vivo clearance in $\log_{10}(\text{mL/min/kg})$.
  * *Verdict:* **`QUALIFIED`** (Sole qualified continuous scaled clearance model).
* **Candidate 2: TDC / ADMET-AI HLM Classification**
  * *Architecture:* D-MPNN binary classifier (threshold at $T_{1/2} = 30\,\text{min}$).
  * *Verdict:* **`REJECTED`** (Scientifically incompatible: binary class output cannot ensemble with continuous $\log_{10}(\text{mL/min/kg})$).

---

### 5. CYP3A4 Functional Inhibition
* **Candidate 1: Admetica Chemprop CYP3A4 Inhibitor**
  * *Architecture:* D-MPNN binary classifier.
  * *Training Set:* PubChem AID 1851 (12,997 compounds).
  * *Performance:* Balanced Accuracy $= 0.829$, Sensitivity $= 0.842$, Specificity $= 0.815$.
  * *Verdict:* **`QUALIFIED`** (Installed Primary Model).
* **Candidate 2: ADMET-AI Chemprop 5-Fold CYP3A4**
  * *Architecture:* 5-model D-MPNN ensemble.
  * *Training Set:* TDC Veith / AID 1851 (12,328 compounds).
  * *Performance:* $\text{AUROC} = 0.902$, $\text{AUPRC} = 0.895$.
  * *Verdict:* **`QUALIFIED_WITH_LIMITATIONS`** (Shared AID 1851 lineage).
* **Candidate 3: CatBoost GBDT CYP3A4 Classifier**
  * *Architecture:* CatBoost on ECFP4 + 2D Descriptors.
  * *Performance:* $\text{AUROC} = 0.880$, Balanced Accuracy $= 0.810$.
  * *Verdict:* **`QUALIFIED`** (High architectural diversity).

---

### 6. Metabolic Site-of-Metabolism (SoM)
* **Candidate 1: SyGMa Phase I & Phase II Rule Engine**
  * *Architecture:* Expert SMARTS reaction matching + atom vulnerability scoring.
  * *Performance:* Top-3 soft-spot recovery $82\%$ on known drug set.
  * *Verdict:* **`QUALIFIED`** (Installed Primary Engine).
* **Candidate 2: SMARTCyp Quantum Mechanical DFT Model**
  * *Architecture:* Density functional theory (B3LYP/6-31G*) activation barrier energies ($E_{\text{act}}$) + topological steric factors.
  * *Performance:* Top-2 accuracy $76-86\%$ across CYP 3A4, 2D6, 2C9.
  * *Verdict:* **`QUALIFIED_WITH_LIMITATIONS`** (Exceptional physical diversity; must be aggregated via Reciprocal Rank Fusion).
* **Candidate 3: DeepMetab Graph Neural Network**
  * *Architecture:* Attentive Graph Neural Network on BioPrint/DrugBank.
  * *Verdict:* **`RESEARCH_ONLY`** (**REJECTED for Production** due to non-commercial academic license).

---

### 7. Quantitative pKa & LogD Models (Secondary Evaluation)
* **IonizationEngine_v1:** Deterministic Henderson-Hasselbalch + RDKit SMARTS micro-states (**`QUALIFIED`**).
* **MolGpKa:** Graph convolutional network for atomic pKa (**`REJECTED`** due to DGL/PyG ARM64 compilation failures).
* **Schrodinger Epik / ACD/Percepta:** Commercial pKa/logD engines (**`REJECTED`** due to paid proprietary licensing).
