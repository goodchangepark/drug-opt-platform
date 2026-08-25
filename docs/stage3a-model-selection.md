# Stage 3 Step 3A model selection

Research frozen on 2026-08-25. Scientific endpoint fit and auditability were weighted above
installation convenience.

## Candidate review

| Candidate | Endpoint/data/output | Validation evidence | License/use | ARM64 and reproducibility | Decision |
| --- | --- | --- | --- | --- | --- |
| ADMET-AI v2 | Chemprop v2 models trained on TDC endpoints, including AqSolDB and Caco2_Wang | TDC Caco2 leaderboard lists Chemprop/RDKit MAE 0.330 and Chemprop MAE 0.344, but the listed `admet_ai_v2` submission is 6.328, suggesting an unresolved transform/submission mismatch | MIT; commercial/internal use allowed | CPU supported; current v2 differs from the paper/web-server v1 | Not selected: endpoint transform/version provenance is presently too easy to misapply |
| TDC | Benchmark/data framework: Solubility_AqSolDB (LogS) and Caco2_Wang (LogPapp) | Scaffold-split MAE benchmark; it does not itself provide one canonical deployable checkpoint | Apache-2.0 code; individual dataset provenance still applies | ARM64-feasible Python framework; training choices must be reproduced locally | Data/benchmark reference, not selected as a pretrained model |
| DeepChem | Delaney aqueous-solubility training example | Example trains a new model; no canonical production checkpoint/version for this platform | MIT | ARM64 possible, but checkpoint and split must be created | Not selected: framework/example rather than endpoint checkpoint |
| Chemprop pretrained/examples | General D-MPNN framework; solubility quickstart data is an example | Requires local training/split choice; no canonical Caco-2 + solubility pair | MIT | ARM64 CPU works | Runtime selected, but not a model source |
| Admetica | Chemprop regression checkpoints; AqSolDB LogS and Wang-derived Caco-2 LogPapp | Publisher: Solubility MAE 0.714/RMSE 1.089/R² 0.788/Spearman 0.897; Caco-2 MAE 0.317/RMSE 0.415/R² 0.701/Spearman 0.832. External 34-compound Caco-2 check is reproducible | MIT repository/checkpoints; commercial/internal R&D allowed under MIT terms | Verified on aarch64 with PyTorch 2.8 CPU and Chemprop 2.2.4; source commit and hashes packaged | Selected |

Sources: [Admetica](https://github.com/datagrok-ai/admetica),
[ADMET-AI](https://github.com/swansonk14/admet_ai),
[TDC Caco2 benchmark](https://tdcommons.ai/benchmark/admet_group/01caco2/),
[AqSolDB paper](https://doi.org/10.1038/s41597-019-0151-1),
[Wang et al.](https://doi.org/10.1021/acs.jcim.5b00642), and
[Chemprop documentation](https://chemprop.readthedocs.io/).

## Implemented endpoint definitions

- Solubility: `LogS = log10(S [mol/L])` from the heterogeneous AqSolDB aggregate. It is not
  asserted to be intrinsic solubility or solubility at a requested pH. No pKa/pH-dependent value
  is generated.
- Permeability: `LogPapp = log10(Papp [cm/s])` for Caco-2 cell monolayers. The distributed
  aggregate retains neither A→B/B→A direction nor detailed assay conditions, so no direction is
  claimed. PAMPA, MDCK, efflux ratio, and intrinsic passive permeability are explicitly excluded.

The checkpoint version is `admetica-d4f7056-chemprop-v2.1`: Admetica commit `d4f7056`, with the
original Chemprop 2.0 checkpoint converted using Chemprop 2.2.4's official v2.0→v2.1 converter.

## Applicability domain and confidence

For each endpoint, the system compares a Morgan radius-2/2048-bit fingerprint with every usable
training structure and reports maximum Tanimoto similarity and `1 - similarity` chemical-space
distance. MW, cLogP, TPSA, HBD, HBA, and rotatable-bond values are also checked against the full
training range.

- `IN_DOMAIN`: similarity ≥ 0.40 and every descriptor within range.
- `BORDERLINE`: similarity ≥ 0.25 and at most one descriptor outside range.
- `OUT_OF_DOMAIN`: otherwise.

These thresholds are transparent heuristics, not a native/calibrated model AD. A single checkpoint
does not support ensemble disagreement, so `HIGH` confidence is never emitted. `MEDIUM` is limited
to `IN_DOMAIN`; `BORDERLINE` and `OUT_OF_DOMAIN` are `LOW`. Uncertainty remains null with an explicit
reason.

## Scientific validation

`scripts/validate_stage3a_models.py` reproduces:

- Caco-2: 34 public Pham-The et al. structures filtered by the upstream workflow against Admetica
  training structures. Reproduced MAE 0.411552, RMSE 0.535530, R² 0.319010 in log10(cm/s), with the
  high-vs-low reference direction correct.
- Solubility: ethanol 0.9922 vs 4,4'-dichlorobiphenyl -6.9126 LogS, the expected direction. Both
  structures occur in AqSolDB, so this is only an execution/directional sanity check, not an
  independent validation. The broad public training aggregate prevented claiming independence.

Known limitations include heterogeneous source assays, incomplete Caco-2 protocol/direction
metadata, unknown publisher validation split details, no calibrated predictive uncertainty, and a
heuristic rather than model-native AD.
