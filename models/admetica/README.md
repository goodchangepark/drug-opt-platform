# Stage 3A model assets

Only the aqueous-solubility and Caco-2 checkpoints are integrated. The assets came from the
[Datagrok Admetica repository](https://github.com/datagrok-ai/admetica) at commit
`d4f70569901c189f39fa37871e2aeabeef3adc83` under its MIT license. Original Chemprop 2.0
checkpoints were deterministically converted with Chemprop 2.2.4's
`v2_0_to_v2_1` converter so they load on Python 3.11 / ARM64.

| Endpoint | Model SHA-256 | Training CSV SHA-256 |
| --- | --- | --- |
| Caco-2 | `87bedd7ea3b314557f803b6ec3f3e7726dba0cdcd19128fd232b2d10511bc8d1` | `685364e2d7607caa88b37d447e4910ec17aa33964f9947bb12881241fecc335f` |
| Solubility | `1e41f9e30687954657d3ac19fd579febc1d66452bfa718c9951b87da115dd6f7` | `beffa08a63196e97e9e5d25a99300121a6314f71a23fc79716d6028f4f89ff94` |

AD-index SHA-256: Caco-2 `96e1001c24bd43c3c47802c2f2f87cbfb225f9feff71f5a0926a8574f4d27478`;
Solubility `7c62d63d0e345b119f1ea14200feb6d14b35a3f07da58ec58ca10594135a3056`.

The external Caco-2 reference file has SHA-256
`2ac729c783a8b9995d99e2794f36bc3bf6f9d693faa3e1194ceb12b842899889`. It contains 34
public Pham-The et al. structures that the Admetica comparison workflow filtered against the
Admetica training structures. Run `python scripts/validate_stage3a_models.py` to reproduce the
scientific checks.

The training CSVs are packaged only to calculate a transparent, heuristic applicability domain:
nearest Morgan/Tanimoto similarity plus a molecular-descriptor range envelope. It is not a
model-calibrated uncertainty estimate.
