# Stage 4B Analog Proposal Engine

Stage 4B turns executable Stage 4A strategies into explicit, scoped molecular
edits. It does not use random generation, an LLM, or PK calculations. Every
candidate records its parent `CompoundVersion`, immutable `OptimizationRun`,
transformation provenance, generated structure, and prediction snapshots.

## Job and data architecture

`OptimizationProposalRun` owns `OptimizationCandidate` rows. Each candidate has
one or two `CandidateTransformation` records, endpoint-specific
`CandidatePredictionSnapshot` records, versioned `CandidateRanking` snapshots,
and zero or more `CandidateRejectionReason` records. Job states are `PENDING`,
`GENERATING`, `FILTERING`, `PREDICTING`, `RANKING`, `COMPLETED`, and `FAILED`.
Candidate-local failures are isolated.

Runs snapshot the parent and strategy run IDs, transformation-library version,
model versions, endpoint weights, hard constraints, settings, seed, and
timestamps. The seed is recorded although the current generator is fully
deterministic.

## Generation and scope

The implementation executes 13 local graph edits: targeted fluorination,
steric shielding, methyl deletion, benzylic-linker contraction, aromatic
carbon-to-nitrogen replacement, N/O-alpha fluorination, phenyl-to-heteroaryl,
aliphatic carbon-to-oxygen, alkyl trimming, polar substituent introduction,
ring heteroatom replacement, and linker oxygen replacement. Four Stage 4A
entries remain `STRATEGY_ONLY`: context-dependent basic-center attenuation,
whole-ring saturation, amide/sulfonamide bioisosterism, and generic alert
removal.

Only stored modifiable/soft-spot atoms and the exact source atoms of ranked
Stage 4A strategies are eligible. Protected or `HIGH-RISK TO MODIFY` parent
atoms are hard exclusions unless Stage 4A manual override removed protection.
Single changes are generated first. Two edits are allowed only for distinct
liabilities, after parent-level protected-region, stereochemistry, and
similarity validation. Three changes are never generated.

Project-observed MMP structures rank first and preserve linked experimental
activity/ADMET evidence. Other priorities are project SAR/soft spots, curated
rules, then generic low-confidence strategies. The engine permits fewer than
50 raw candidates when scoped chemistry does not justify filler structures.

## Staged filtering and rescoring

1. RDKit sanitization, valence/aromaticity, fragmentation, charge, canonical
   duplicate, stereochemistry, protected-region, Morgan similarity, and MCS
   coverage checks.
2. Stage 1 properties/alerts and hard MW, cLogP, TPSA, similarity, and new-alert
   gates.
3. Selected project assay QSAR or nearest-neighbour prediction, with project
   experimental evidence taking precedence for an existing MMP analog.
4. All installed Stage 3 models: solubility, Caco-2, PPB/fu, HLM/RLM/MLM, CYP,
   P-gp inhibitor, hERG, Ames, and DILI. `MODEL_UNAVAILABLE` endpoints are
   recorded but omitted from scoring.
5. SyGMa/RDKit soft-spot reanalysis and parent/candidate Top-3 liability change.

The SA score from Ertl and Schuffenhauer plus ring size, stereocenter change,
and Bertz-complexity delta forms a `LOW`, `MODERATE`, or `HIGH SYNTHETIC
COMPLEXITY` surrogate. It is explicitly not a synthesis-success probability or
retrosynthetic analysis.

## Hard gates and uncertainty

Configured MW, cLogP, TPSA, similarity, no-new-alert, absolute potency,
potency-fold, LogS, and Caco-2 limits become candidate gates. OUT-OF-DOMAIN
activity is strongly penalized instead of trusted for automatic potency
rejection. A LOW-confidence hERG/CYP/P-gp classification can add a ranking risk
but cannot by itself reject a candidate. Experimental evidence and compatible
project MMP/SAR stay above external prediction.

## Ranking formula

Available objectives are normalized to `[0, 1]`. Unavailable endpoints do not
enter either numerator or denominator. The disclosed score is:

`100 × max(0, Σ(wᵢ × qᵢ × cᵢ × dᵢ) / Σ(usable wᵢ) − Palerts − Psynthetic − POOD)`

where `q` is objective quality, `c` is evidence-confidence factor (experimental
or HIGH 1.0, MEDIUM 0.75, LOW 0.5, UNKNOWN 0.35), and `d` is applicability
factor (IN 1.0, BORDERLINE 0.7, OUT 0.3, UNKNOWN 0.5). Stage 4A objectives boost
their endpoint weights; explicit endpoint weights override defaults.

The score does not replace non-dominated sorting. Pareto fronts are calculated
over every mutually available objective. Final selection greedily preserves
distinct transformation hypotheses, sites, and Morgan-fingerprint diversity.
Information Value is HIGH for a distinct hypothesis on Pareto front 1 or a
chemically diverse discriminator, MEDIUM for partially novel information, and
LOW for redundant nearby chemistry.

## Limitations

- Local graph edits are hypotheses, not reaction conditions or guaranteed
  isolable compounds.
- SA score and complexity are coarse synthesis surrogates.
- Prediction uncertainty, heterogeneous training assays, and applicability
  domains remain endpoint-specific; rescoring does not make them experimental.
- MCS preservation is deterministic but is not a pharmacophore model.
- Soft-spot ranks are empirical rule priors, not calibrated atom probabilities.
- Diversity and Information Value are deterministic experiment-design
  heuristics, not measures of biological success.
- The public acceptance cases check expected direction and pipeline behavior;
  they are not prospective validation.
