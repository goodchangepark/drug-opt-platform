# OpenADMET microsomal-clearance model asset

- Model: `openadmet/microsomal-clearance-chemeleon-v1`
- Source revision: `e13549384674e70a536097b7175932e36d5ff271`
- Artifact SHA-256: `dd143760b58af67a9c698759830c62b7daabaedbd17b96556b1c508d927a9ffc`
- Model license: Apache-2.0 (as declared by the Hugging Face model card)
- Framework source: OpenADMET Models, MIT license
- Outputs: `LOG_CLint_HLM`, `LOG_CLint_RLM`, `LOG_CLint_MLM`, each scaled
  `log10(mL/min/kg)` and species-isolated.

`X_train.csv` and `y_train.csv` are retained to make the applicability-domain
indexes reproducible. The released checkpoint was trained on all ExpansionRx
train and test data; benchmark plots in its card were made with an analogous
checkpoint that excluded the ExpansionRx test set.

The Biogen prospective public validation CSV is redistributed from
`https://github.com/biogen/ADME` at commit
`b00df003de117ce9e5b381afd886095c5f2af2d5` under that repository's MIT
license. Its SHA-256 is
`2cfabc2667740c224487876c33b23124159ef43294e0f9e4d926cb6276c95a3b`.
