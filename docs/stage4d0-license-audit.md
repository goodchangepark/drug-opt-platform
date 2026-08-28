# Stage 4D-0: Comprehensive Intellectual Property & License Audit

**Platform Version:** `0.6.3-stage5b4-ui`  
**Audit Scope:** Source Code Licenses, Model Checkpoint/Weight Terms, and Training Dataset Lineage Terms.

---

## 1. Tripartite License Governance

Legal qualification of scientific AI models requires auditing three distinct legal layers:
1. **Layer 1: Software Code License:** The legal license governing the model architecture, feature extraction code, and runtime inference scripts.
2. **Layer 2: Pretrained Weights / Checkpoint License:** The specific terms attached to the trained neural network or tree weights file.
3. **Layer 3: Upstream Training Dataset Terms:** The data usage rights and attribution obligations attached to the underlying assay database.

A model is approved for Drug-OPT production ONLY if all three layers permit local redistribution, company internal use, and reproducible offline execution.

---

## 2. Detailed License Audit Matrix

| Model Family | Code License | Checkpoint License | Training Data Terms | Qualification Status | Legal Verdict |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Admetica (Datagrok-AI)** | MIT | MIT | Source-specific public | **QUALIFIED** | Fully permissive; approved for company production. |
| **OpenADMET (CheMeleon)** | Apache-2.0 | Apache-2.0 | ChEMBL (CC BY-SA 3.0), ASAP | **QUALIFIED** | Permissive open-source; approved. |
| **ADMET-AI (Swanson et al.)**| MIT | MIT | TDC (CC BY 4.0 / Public) | **QUALIFIED** | Permissive code & weights; dataset attribution required. |
| **SyGMa (Ridder et al.)** | GPL-3.0 | GPL-3.0 | Expert Literature Rules | **QUALIFIED** | Copyleft; code runs as an isolated subprocess/module. |
| **SMARTCyp (Rodet et al.)** | LGPL-3.0 | LGPL-3.0 | Quantum DFT Calculations | **QUALIFIED** | Permissive dynamically linked module; approved. |
| **AqSolPred / GBDT Models** | MIT / BSD | MIT / BSD | AqSolDB (CC BY 4.0) | **QUALIFIED** | Permissive; approved. |
| **DeepMetab (Wu et al.)** | Academic-Only | Non-Commercial | Proprietary BioPrint + DrugBank| **RESEARCH_ONLY** | **REJECTED for Production**; strictly non-commercial license. |
| **SwissADME / pkCSM** | Proprietary Web | Closed Source | Literature | **REJECTED** | No local offline weights; Terms of Service prohibit API scraping. |
| **Schrodinger / Epik / QikProp**| Proprietary | Proprietary | Proprietary Pharma | **REJECTED** | Commercial paid license required; cannot distribute. |

---

## 3. Copyleft & Commercial Isolation Policies

### 3.1 SyGMa & SMARTCyp Isolation
* **SyGMa** is distributed under the GNU General Public License v3.0 (GPL-3.0).
* **SMARTCyp** is distributed under the GNU Lesser General Public License v3.0 (LGPL-3.0).
* **Compliance Architecture:** In Drug-OPT, SyGMa is invoked strictly as an independent service via well-defined API boundaries, preserving full license compliance without contaminating the proprietary drug discovery application codebase.

### 3.2 Dataset Attribution Obligations
* **AqSolDB & TDC:** Licensed under Creative Commons Attribution 4.0 International (CC BY 4.0).
* **Compliance:** The platform maintains explicit provenance metadata and references in `/api/help/registry` and all exported validation manifests.
