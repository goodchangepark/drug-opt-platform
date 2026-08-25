# Packaged SyGMa reaction rules

`phase1.txt` and `phase2.txt` are the unmodified rule files distributed with
SyGMa 1.1.0 from <https://github.com/3D-e-Chem/sygma>.

They are used directly by the Stage 3D RDKit transformation engine so a fresh
installation does not depend on SyGMa's legacy source-distribution build step.
The files contain 148 Phase I and 27 Phase II SMIRKS rules with empirical
occurrence priors derived for the published SyGMa method. Those priors are not
treated as calibrated atom probabilities.

SyGMa and these rule files are licensed under GPL-3.0. See
`LICENSE-GPL-3.0`. The platform's RDKit execution code records the source,
version, license, matched rule, original SMIRKS, and occurrence-count comment
in each prediction.
