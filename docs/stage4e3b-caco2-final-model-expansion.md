# Stage 4E-3B — Caco-2 Final Model Expansion & Closure

This was the final dedicated Caco-2 expansion attempt before Engine-v1
freeze. Search criteria were frozen before candidate screening and limited to
three serious primary-source candidates. No candidate was installed,
registered, benchmarked, or promoted.

The Zenodo GNN-MTL model publishes a small Chemprop checkpoint and a Caco-2
Papp output, but the record provides copyright rather than a clear internal
research weight license, does not publish a SHA256, and does not unambiguously
define the A→B assay direction. MolGrad/Ersilia exposes passive permeability,
but its model is AGPL-3.0-only and the deployment is AMD64. PharmPapp provides
peptide-oriented KNIME archives without a direct small-molecule Papp A→B
contract or clear release licensing. These failures are recorded fail-closed
in the machine-readable artifacts.

Because no candidate passed legal, checkpoint, and direct endpoint gates,
ExpansionRx was not reused for candidate selection and no candidate inference
or bootstrap was run. Stage 4E-3A remains the authoritative fixed CORE/SHADOW
external evidence: CORE is retained with MAE 0.5695, RMSE 0.7457, and Spearman
0.0410; the SHADOW was worse by a paired MAE delta of +0.1352.

Final closure: `CACO2_NO_QUALIFIED_REPLACEMENT_FOUND_CORE_FROZEN`.
Production remains `SINGLE_CORE_MODEL` using
`admetica_caco2` version `admetica-d4f7056-chemprop-v2.1`. No further
dedicated Caco-2 optimization is planned before Engine-v1 freeze. The
reliability state remains limited / low-medium, with residual training overlap
unknown and assay heterogeneity as explicit limitations.
