# Stage 3D model and tool selection

## Decision

Stage 3D uses the public SyGMa 1.1.0 empirical reaction rules as complementary chemistry evidence and RDKit 2025.03.1 for atom mapping, structure generation, sanitization, canonicalization, duplicate removal, and depiction. It does **not** emit an atom-level ML probability. All generated structures are labeled `PREDICTED METABOLITE HYPOTHESIS` and all rule-only soft spots have `LOW` confidence.

No evaluated atom-level pretrained model qualified as the primary installed model:

| Candidate | Phase / CYP scope | Training / reported validation | License | ARM64, dependencies, reproducibility and deployment result |
|---|---|---|---|---|
| FAME 3 / FAME3R | Phase I and II, atom level; FAME 3 is not a per-isoform CYP claim | MetaQSAR: 1,733 train and 434 test molecules | FAME3R source is MIT, but required FAME 3 model data are restricted to non-commercial research | Source can run on CPU/RDKit-class environments, but the weight license blocks the intended internal/commercial deployment regardless of architecture. Not installed. |
| XenoSite | Phase I CYP/HLM atom level; published CYP settings include nine isoforms | More than 680 substrates; publisher reports about 83–89% top-2 across CYP settings | Web service/paper available; no deployable checkpoint license identified | No maintained reproducible local checkpoint; the historical stack used MOPAC plus OpenBabel/OpenEye-era tooling and was not reproducibly qualified on ARM64. Not installed. |
| ATTNSOM | Phase I SOM, atom level; CYP-focused attention model | 679 Zaretzki substrates plus 120 public AstraZeneca exact-SOM compounds | MIT source and data | PyTorch/PyTorch-Geometric training code is portable in principle, but GPU-oriented retraining would be required and no published pretrained checkpoint exists. Not installed. |
| SMARTCyp 2.4.2 | Phase I CYP2C9/2D6/3A4 atom ranking | Energy/accessibility model; published external CYP2C9 top-1/2/3 roughly 42/58/67% | LGPL-3.0 CDK port | Java/CDK can be architecture-neutral with an ARM JVM, but the public port identifies itself as not tested and this environment has no validated Java/Maven runtime for it. Not installed. |
| GLORYx | Phase I/II metabolite generation; atom ranking inherited from FAME 3 | Publisher reports 77% recall and AUC 0.79 | GPL-3.0 source; FAME 3 model data retain separate restrictions | Python/RDKit plus separately obtained FAME 3 files; reproducibility and deployability are blocked by the missing restricted model data. Not installed. |
| BioTransformer | Phase I/II reaction prediction, including enzyme-associated transformations; not primarily an atom-level SOM model | Knowledge/rule-based multi-step metabolite prediction | LGPL-3.0/GPL components as distributed by the project | JVM/CDK tooling is generally ARM64-capable with a compatible Java runtime, but it does not satisfy the requested primary atom-level pretrained model role. Not selected. |
| SyGMa 1.1.0 | 148 Phase I and 27 Phase II empirical reaction rules; not CYP-isoform-specific | Publisher: 68% metabolite recall on 175 parents; 45% within top 10. Single-step CYP: 84% reaction recall; 66% within top 3. | GPL-3.0 | Python/RDKit CPU path runs reproducibly on the current ARM64 environment. Rule files are packaged locally. Selected as rule evidence only. |

The packaged SyGMa rule tables are copied unmodified from release 1.1.0 and retain GPL-3.0. Only the supported Stage 3D transformation families are executed.

## Scientific interpretation

The stored score is the maximum SyGMa empirical reaction-rule occurrence prior among matching rules. It is neither a calibrated probability nor a model confidence. Site confidence remains `LOW` because no qualified atom-level model checkpoint or rigorous applicability-domain calculation is available.

CYP substrate predictions from Stage 3C and experimental/predicted HLM/RLM evidence from Stage 3B are retained as compound-level supporting evidence. They never assign a CYP isoform to a soft-spot atom. Raw microsomal experimental results are preferred for display even when scaling data are absent; in that case no prediction-unit conversion or stability threshold is applied.

Phase II glucuronidation and sulfation are rule hypotheses only. They do not model enzyme isoform, tissue exposure, competing pathways, cofactor availability, or sequential kinetics.

## Validation boundary

The reproducible artifact in `models/sygma/validation/known_drug_sanity.json` is a five-drug directional sanity check, not an independent benchmark. SyGMa's historical MDL Metabolite source database is discontinued and the training compound list is not distributed, so overlap cannot be audited. Publisher metrics and local sanity metrics are therefore reported separately.
