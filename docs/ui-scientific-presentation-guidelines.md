# Drug-OPT UI/UX Scientific Presentation Guidelines

**App Version:** `0.6.1-stage5b4-ui`  
**Date:** 2026-08-27  
**Scope:** Medicinal Chemistry, ADMET, Metabolism, and Translational PK Visualization

---

## 1. Core Principles

1. **Less Clutter, More Interpretation:**
   - Avoid dumping unformatted numeric columns without scientific meaning.
   - Every property and predicted endpoint MUST be paired with a reference range, clinical threshold, or standard interpretation.

2. **Deterministic & Rule-Based Grounding:**
   - Scientific evaluations are computed deterministically via the Centralized Interpretation Registry (`backend/interpretation.py`).
   - Interpretations reference established literature (Lipinski 2001, Veber 2002, Lovering 2009, Bickerton 2012, FDA BCS, FDA DDI Guidelines).

3. **Strict Color Semantics:**
   - **Blue (`#1a56db` / `.badge-favorable`, `.dot-favorable`):** Favorable profile, compliant with drug-likeness, low liability, or safe zone.
   - **Gray (`#6b7280` / `.badge-intermediate`, `.dot-intermediate`):** Intermediate, moderate, uncertain, or borderline values.
   - **Red (`#e02424` / `.badge-liability`, `.dot-liability`):** Unfavorable profile, drug-likeness rule violation, high cardiac/liver toxicity liability, or high clearance.
   - **Accessibility Rule:** Color must NEVER be used alone; always accompany status dots and badges with explicit descriptive text.

4. **Typography & Layout Hierarchy:**
   - **Font Stack:** `"Noto Sans KR", "Noto Sans CJK KR", -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif`.
   - Clear contrast hierarchy: H1 (22px bold), H2 (17px semibold), H3 (14.5px semibold), Table headers (11px uppercase tracking), Monospace for chemical formulas and numerical values (`ui-monospace`, `SFMono-Regular`, `Consolas`).

---

## 2. Tabular Layout Standards

### A. Physicochemical Table (Properties Tab)
Single 4-column aligned table:
1. `Property`: Plain language property name with mathematical abbreviation.
2. `Calculated Value`: Formatted numeric value with explicit unit.
3. `Drug-Likeness / Reference Range`: Standard medicinal chemistry threshold (e.g. `≤ 500 g/mol`, `≤ 140 Å²`, `≥ 0.42`).
4. `Assessment`: Semantic status badge (`FAVORABLE`, `INTERMEDIATE`, `UNFAVORABLE`) with rule description.

### B. ADMET Prediction Table (ADMET Tab)
Standard 4-column matrix:
1. `Evaluation / Endpoint`: Target developability property.
2. `Prediction (Value & Unit)`: Consensus prediction and individual model outputs (M1, M2 chips).
3. `Scientific Interpretation`: Deterministic textual summary and reference threshold.
4. `Model System`: Model family, algorithm version, and training provenance.

### C. Cross-Species Metabolic Stability Table (Metabolism Tab)
Cross-species microsomal stability comparison across **Human**, **Rat**, and **Mouse**:
- Reports Clint in `µL/min/mg protein` and physiological `mL/min/kg`.
- Clear hepatic extraction risk category (`Low Extraction`, `Moderate Extraction`, `High Extraction`).

### D. Multi-Species PK Summary Table (PK Tab)
Comparative pharmacokinetic matrix across **Mouse**, **Rat**, **Dog**, **Monkey**, and **Human**:
- Displays in vivo NCA parameters (Clearance, Vss, Half-life, AUC, Bioavailability) side-by-side.
- Reports human translational clearance hierarchy: `Experimental > Allometry > IVIVE`.

---

## 3. Visual Profile Chart (`VisualProfileChart`)

The ADMET visual profile chart displays 8 key developability dimensions normalized on a 0–100 qualitative scale:
- `Solubility (Aqueous)`
- `Permeability (Caco-2)`
- `Plasma Free Fraction (fu)`
- `Microsomal Stability (HLM)`
- `Cardiac Safety (hERG Neg)`
- `Liver Safety (DILI Neg)`
- `Mutagenic Safety (Ames Neg)`
- `Drug-Likeness (QED)`

Bars are colored according to favorable (blue), intermediate (gray), or liability (red) status to allow medicinal chemists to rapidly identify developability bottlenecks.
