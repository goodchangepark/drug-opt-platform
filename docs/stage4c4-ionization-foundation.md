# Stage 4C-4: Ionization, pKa & pH-Dependent Physicochemical Foundation

## 1. Overview & Theoretical Framework

Ionization is a primary determinant of drug-like behavior, governing:
- **Aqueous solubility** ($S_{\text{pH}} \gg S_0$ when ionized).
- **Passive membrane permeability** (predominantly uncharged species cross lipid bilayers by passive diffusion according to the pH-partition hypothesis).
- **Plasma protein binding** (acidic drugs bind primarily to Human Serum Albumin; basic drugs bind strongly to $\alpha_1$-acid glycoprotein).
- **Volume of distribution** ($V_d$) (basic drugs undergo lysosomal trapping in acidic intracellular organelles, expanding $V_d$; acidic drugs are restricted to plasma and extracellular water, yielding lower $V_d$).
- **Oral fraction absorbed** ($F_a$) (differential ionization across the stomach pH 1.2–2.0, proximal intestine pH 4.5–6.5, and distal intestine/blood pH 7.4).

Stage 4C-4 establishes a deterministic, scientifically rigorous ionization layer that connects molecular structure directly to downstream ADME/PK interpretation without fabricating unavailable quantitative values.

---

## 2. Ionization Classification & Structural Motif Engine

Molecules are deterministically classified into one of seven structural ionization classes based on rigorous SMARTS pattern matching:

1. `NEUTRAL`: No ionizable acidic or basic centers with physiological relevance ($1.0 \le \text{pH} \le 14.0$).
2. `ACID`: Contains one or more acidic centers (e.g. carboxylic acid, sulfonic acid, phosphonic acid, tetrazole, acidic sulfonamide) with zero basic centers.
3. `BASE`: Contains one or more basic centers (e.g. aliphatic amine, amidine, guanidine, basic pyridine, imidazole) with zero acidic centers.
4. `AMPHOLYTE`: Contains at least one acidic center and at least one basic center where macroscopic or microscopic charges are separated.
5. `ZWITTERION_POSSIBLE`: A subset of ampholytes where simultaneous ionization to a net-neutral zwitterion ($+1/-1$) occurs at physiological pH (e.g. alpha-amino acids like amoxicillin, fluoroquinolones like ciprofloxacin).
6. `MULTIPLE_IONIZABLE_CENTERS`: Contains 3 or more ionizable centers (complex polyprotic species).
7. `REVIEW_REQUIRED`: Unresolved charge patterns, unusual quaternary centers, or structural anomalies requiring manual chemist inspection.

### SMARTS Motif Rules & Ionizable Centers
- **Acidic Motifs**:
  - Carboxylic Acid: `[CX3](=O)[OX2H1]` ($\text{pKa} \approx 3.5 - 5.0$)
  - Sulfonic Acid: `[SX4](=O)(=O)[OX2H1]` ($\text{pKa} \approx 0.5 - 2.0$)
  - Phosphonic / Phosphoric Acid: `[PX4](=O)([OX2H1])[OX2H1]` ($\text{pKa}_1 \approx 1.5 - 2.5, \text{pKa}_2 \approx 6.5 - 7.5$)
  - 1H-Tetrazole: `c1nnn[nH]1` ($\text{pKa} \approx 4.5 - 5.5$)
  - Primary / Secondary Sulfonamide (N-H): `[SX4](=O)(=O)[NX3;H1,H2][#6]` ($\text{pKa} \approx 5.5 - 8.5$)
  - Acidic Phenol: `[OX2H1]c1ccccc1` ($\text{pKa} \approx 8.5 - 10.5$)
  - Imide / Cyclic Barbiturate: `[CX3](=O)[NX3H1][CX3](=O)` ($\text{pKa} \approx 7.0 - 9.0$)
- **Basic Motifs**:
  - Primary Aliphatic Amine: `[NX3;H2;!$(NC=O);!$(NS=O);!$(N=O);!$(nc)]` ($\text{pKa} \approx 9.5 - 10.5$)
  - Secondary Aliphatic Amine: `[NX3;H1;!$(NC=O);!$(NS=O);!$(nc)]([#6])[#6]` ($\text{pKa} \approx 9.5 - 10.8$)
  - Tertiary Aliphatic Amine: `[NX3;H0;!$(NC=O);!$(NS=O);!$(nc)]([#6])([#6])[#6]` ($\text{pKa} \approx 8.5 - 10.2$)
  - Amidine: `[NX3;H1,H2][CX3]=[NX2]` ($\text{pKa} \approx 11.0 - 12.5$)
  - Guanidine: `[NX3;H1,H2][CX3](=[NX2])[NX3;H1,H2]` ($\text{pKa} \approx 12.0 - 13.5$)
  - Pyridine / Diazine: `[nX2;r6]` ($\text{pKa} \approx 4.0 - 6.0$)
  - Imidazole: `[nX2;r5]1cc[nH]1` ($\text{pKa} \approx 6.0 - 7.5$)
  - Basic Aniline: `[NX3;H1,H2]c1ccccc1` ($\text{pKa} \approx 3.5 - 5.0$)
- **Explicit Non-Basic Exclusions**:
  - Amides (`NC=O`), Carbamates (`OC(=O)N`), Ureas (`NC(=O)N`), Pyrrole-like aromatic nitrogens (`[nH]`), Nitro groups (`N(=O)=O`), Cyano nitrogens (`C#N`).

---

## 3. Mathematical Equations for pH-Dependent Ionization

### 1. Henderson-Hasselbalch Monoprotic Ionization
For a monoprotic acid ($HA \rightleftharpoons H^+ + A^-$):
$$\text{Fraction Ionized } (f_{\text{ionized}}) = \frac{1}{1 + 10^{\text{pKa} - \text{pH}}}$$
$$\text{Fraction Neutral } (f_{\text{neutral}}) = 1 - f_{\text{ionized}} = \frac{1}{1 + 10^{\text{pH} - \text{pKa}}}$$

For a monoprotic base ($BH^+ \rightleftharpoons H^+ + B$):
$$\text{Fraction Protonated / Ionized } (f_{\text{ionized}}) = \frac{1}{1 + 10^{\text{pH} - \text{pKa}}}$$
$$\text{Fraction Neutral / Free Base } (f_{\text{neutral}}) = 1 - f_{\text{ionized}} = \frac{1}{1 + 10^{\text{pKa} - \text{pH}}}$$

### 2. cLogP vs logD
Calculated $\text{cLogP}$ (via RDKit Crippen) measures the intrinsic partition coefficient of the neutral uncharged molecule between 1-octanol and water:
$$P = \frac{[HA]_{\text{octanol}}}{[HA]_{\text{water}}}$$

Distribution coefficient $\log D_{\text{pH}}$ accounts for the partition of both neutral and ionized species at a given pH:
$$D_{\text{pH}} = \frac{[HA]_{\text{octanol}}}{[HA]_{\text{water}} + [A^-]_{\text{water}}}$$

For a monoprotic acid:
$$\log D_{\text{pH}} \approx \text{cLogP} - \log_{10}\left(1 + 10^{\text{pH} - \text{pKa}}\right)$$

For a monoprotic base:
$$\log D_{\text{pH}} \approx \text{cLogP} - \log_{10}\left(1 + 10^{\text{pKa} - \text{pH}}\right)$$

### 3. pH-Dependent Solubility Estimation
When intrinsic solubility $S_0$ is known:
- For acids: $S_{\text{pH}} = S_0 \times \left(1 + 10^{\text{pH} - \text{pKa}}\right)$
- For bases: $S_{\text{pH}} = S_0 \times \left(1 + 10^{\text{pKa} - \text{pH}}\right)$

---

## 4. Downstream ADME & PK Integration

### 1. Aqueous Solubility
- Generic ML solubility is reported without asserting an unverified assay pH.
- If assay pH is provided in experimental data, the ionization state is reported in context.
- Estimated pH-solubility profiles are labeled as `CALCULATED pH-DEPENDENT SOLUBILITY ESTIMATE` with all underlying assumptions disclosed.

### 2. Membrane Permeability (Caco-2)
- Caco-2 predictions are contextualized with the neutral fraction $f_{\text{neutral}}$ at assay pH (pH 7.4).
- High neutral fraction $\implies$ passive transcellular diffusion favored.
- High ionized fraction $\implies$ passive permeability may be limited; paracellular or transporter-mediated transport may dominate.

### 3. Plasma Protein Binding ($f_u$)
- The fundamental relationship $f_u = 1 - \text{fraction bound}$ is preserved.
- Ionization class (acid vs base vs neutral) is reported to assist in identifying likely binding partners (Albumin vs AAG).

### 4. Volume of Distribution ($V_d$)
- Empirical $V_d$ estimations incorporate the ionization class:
  - Bases: Higher tissue partition and lysosomal trapping potential.
  - Acids: Higher plasma albumin retention and lower tissue penetration.
  - cLogP is never silently substituted for $\log D_{7.4}$ without noting the neutral-state assumption.

### 5. Oral Absorption ($F_a$)
- Contextualizes fraction absorbed with ionization behavior across the gastrointestinal pH gradient:
  - Stomach (pH 1.2–2.0)
  - Duodenum / Proximal Jejunum (pH 4.5–6.5)
  - Distal Ileum / Blood (pH 7.4)

---

## 5. Model Registry & Conformal Governance
- All ionization rules and pKa calculation endpoints adhere to the Stage 4C-3B decoupled governance schema:
  - `Data Provenance`: `EXTERNAL`, `INTERNAL`, `TRAINING_OVERLAP_UNKNOWN`, `UNAVAILABLE`.
  - `Calibration Quality`: `VALIDATED`, `UNDERCOVERED`, `INSUFFICIENT_N`, `UNAVAILABLE`.
- Rule-based classifications are explicitly distinguished from ML predictions.
