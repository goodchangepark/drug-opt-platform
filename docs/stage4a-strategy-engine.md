# Stage 4A Hit Optimization Strategy Engine

Stage 4A assembles existing activity, calculated properties, ADMET, metabolism,
SAR, MMP, and activity-cliff evidence into a deterministic strategy. It ranks
medicinal-chemistry transformations only. It does not execute reaction SMARTS,
generate analogs, calculate an overall candidate score, run PK, or call an LLM.

## Persistent architecture

`OptimizationRun` belongs to one project and one parent `CompoundVersion`, with
an optional active assay. It stores multiple objectives, constraints, endpoint
weights, manual overrides, immutable result snapshots, engine version, status,
and timestamps. Project and CompoundVersion ownership are checked when a run is
created. Runs are not reused across versions.

The explicit evidence hierarchy is:

1. Experimental
2. Project-specific validated model/SAR
3. External validated quantitative model
4. External classification model
5. Rule-based hypothesis

Calculated RDKit properties are recorded as calculated evidence. Experimental
ADMET is preferred whenever endpoint, role, species, and unit are compatible;
the prediction remains preserved in the evidence snapshot. A LOW-confidence
external classification is supporting-only unless corroborated. Probability is
never treated as quantitative IC50 or clearance.

## Regions and transformations

Activity-cliff regions where the parent is experimentally more potent are
marked `HIGH-RISK TO MODIFY`. Tolerated project MMP regions, metabolic soft
spots, structural alerts, high-lipophilicity aromatic regions, and terminal
substituent heuristics can be modifiable evidence. If protection data is
insufficient, the engine records `UNKNOWN` instead of asserting a pharmacophore.

The versioned transformation library covers fluorination, steric shielding,
methyl or benzylic-CH removal, heteroatom replacement, N/O-dealkylation
blocking, phenyl-to-heteroaryl and ring/linker bioisosteres, alkyl reduction,
polar-group introduction, aromaticity/Fsp3 modulation, basic-center attenuation,
and structural-alert removal. Each entry records reaction SMARTS, motif,
expected direction, possible risk, reference, and version. Project-observed MMP
evidence ranks above generic rules. For activity-tolerated project MMPs, compatible
paired experimental/predicted solubility, Caco-2, HLM/RLM, hERG, and P-gp evidence
is retained separately; an observed favorable direction raises the MMP rank without
claiming causality. Protected overlap raises potency risk.

Manual `protect_atoms`, `allow_atoms`, `exclude_transformations`, and
`prioritize_transformations` overrides are stored on the run and cause a
deterministic rerank. They do not edit the parent structure.

## Sources and limitations

- MMP algorithm: Hussain & Rea, DOI `10.1021/ci900450m`.
- Bioisostere observations: SwissBioisostere, DOI `10.1093/nar/gks1059`.
- Fluorination/metabolic strategy context: Johnson et al., DOI
  `10.1021/acs.jmedchem.9b01877`.
- hERG/lipophilicity direction: Waring & Johnstone, DOI
  `10.1016/j.bmcl.2006.12.061`.

Reaction SMARTS describe a strategic edit class, not a validated product.
Direction and potency transfer remain context dependent. MCS-derived changed
atoms are a project-SAR heuristic rather than pharmacophore inference. A stored
similarity threshold and other future-candidate constraints are not evaluated
until a proposal stage exists. The acceptance examples are pipeline sanity
checks and are not independent accuracy validation.
