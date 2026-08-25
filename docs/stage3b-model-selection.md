# Stage 3 Step 3B model selection

Research frozen on 2026-08-25. Endpoint/species fidelity, data provenance, validation evidence,
and honest uncertainty were weighted above installation convenience.

## Candidate review

| Candidate | Endpoint and training data | Output / validation | License / ARM64 / decision |
| --- | --- | --- | --- |
| Admetica PPBR | Human plasma protein binding rate from AstraZeneca ChEMBL `CHEMBL3301361`, 2,790 rows | Percent bound regression. Publisher MAE 6.919 percentage points, RMSE 11.294, R² 0.609, Spearman 0.762 | MIT; commercial/internal R&D permitted by the license. Chemprop 2.2.4 CPU inference verified on aarch64. **Selected for human PPB.** |
| OpenADMET CheMeleon microsomal-clearance v1 | Curated ChEMBL 35, ASAP-Polaris and ExpansionRx; 5,086 HLM, 670 RLM and 5,086 MLM non-missing labels in the packaged training table | Three species-specific regression tasks: `LOG_CLint_HLM/RLM/MLM`, each scaled `log10(mL/min/kg)`. The exact released all-data checkpoint has plots but no numeric held-out metrics in its card | Model Apache-2.0; OpenADMET framework MIT. Strict checkpoint load and CPU inference verified on aarch64. **Selected for HLM/RLM/MLM.** |
| Admetica CL-Micro / TDC Clearance_Microsome_AZ | Human microsomal clearance, AstraZeneca/Di et al., ChEMBL `CHEMBL3301370`, 1,102 rows; HLM at 37 °C, range under 3 to over 150 µL/min/mg | Quantitative `mL/min/g` (numerically equivalent to µL/min/mg); publisher MAE 26.715, RMSE 39.201, R² 0.216, Spearman 0.576 | MIT/Apache-2.0 framework context; ARM64 feasible. Not selected because HLM-only and weaker validation than desired; retained as endpoint/assay reference. |
| ADME@NCATS RLM | Rat liver microsome public subset; stable if t½ >30 min and unstable at ≤30 min | Classification; public web model/data, not a qualified quantitative downloadable checkpoint for this installation | Not selected because checkpoint/license provenance was not sufficient and quantitative species-specific OpenADMET was available. |
| Biogen Computational-ADME models | Proprietary discovery training split with public prospective test structures and HLM/RLM/hPPB endpoints | Repository includes training code/data but no single versioned pretrained checkpoint chosen here | MIT. Used only as an independent public validation set, never as the deployed model. |

Sources: [Admetica](https://github.com/datagrok-ai/admetica),
[OpenADMET model card](https://huggingface.co/openadmet/microsomal-clearance-chemeleon-v1),
[OpenADMET models and datasets](https://openadmet.org/),
[TDC microsomal clearance](https://tdcommons.ai/single_pred_tasks/adme/#microsomal-clearance),
[ChEMBL CHEMBL3301370](https://www.ebi.ac.uk/chembl/explore/assay/CHEMBL3301370),
[ADME@NCATS](https://opendata.ncats.nih.gov/adme/), and
[Biogen Computational-ADME](https://github.com/biogen/ADME).

## Implemented definitions and transformations

- Human PPB is the model's original `% bound`. `fraction_bound = %bound / 100`,
  `fu_fraction = 1 - fraction_bound`, and `fu_percent = 100 - %bound`. The original prediction and
  calculated values are all retained. Non-human PPB is unavailable; human output is never reused.
- HLM, RLM and MLM are independent output heads. Their output is physiologically scaled intrinsic
  clearance `log10(mL/min/kg)`, not raw microsomal `µL/min/mg protein`, hepatocyte
  `µL/min/10^6 cells`, hepatic clearance, or t½.
- Raw experimental microsomal clearance is converted only when provenance supplies both MPPGL and
  liver weight. The recorded equation is
  `scaled mL/min/kg = raw µL/min/mg × MPPGL mg/g × liver weight g/kg / 1000`.
  No t½ is calculated without the incubation volume/protein parameters required by its assay.
- Dog and monkey clearance remain `MODEL_UNAVAILABLE`; no species substitution is performed.

## Applicability domain, confidence, and stability summary

Each endpoint uses only structures with a label for that endpoint. AD combines nearest-training
Morgan radius-2/2048-bit Tanimoto similarity, `1 - similarity` chemical-space distance, and the MW,
cLogP, TPSA, HBD, HBA, and rotatable-bond training envelope. Thresholds remain the transparent
Stage 3A heuristic: IN_DOMAIN at similarity ≥0.40 with all descriptors in range; BORDERLINE at
similarity ≥0.25 with no more than one descriptor outside; otherwise OUT_OF_DOMAIN.

PPB may be MEDIUM only in-domain. The independent clearance results were weak, so HLM/RLM/MLM are
always LOW even if in-domain. No model emits HIGH and no ensemble uncertainty is fabricated.

Stable/moderate/unstable is an operational research summary based on the species-specific 25th and
75th percentiles of the released model's training labels. The exact thresholds and this basis are
returned in every result. `METABOLIC STABILITY CONCERN` is stored for the unstable category. These
data-derived thresholds are not represented as universal assay standards and do not rank compounds.

## Independent and scientific validation

`scripts/validate_stage3b_models.py` excludes exact canonical-SMILES overlap and evaluates the
Biogen public prospective set:

| Endpoint | n | MAE | RMSE | R² | Spearman |
| --- | ---: | ---: | ---: | ---: | ---: |
| Human PPB, percentage points bound | 185 | 14.6194 | 21.7950 | 0.4389 | 0.6105 |
| HLM, log10(mL/min/kg) | 3,078 | 0.6259 | 0.7616 | -0.4911 | 0.3700 |
| RLM, log10(mL/min/kg) | 3,045 | 0.6263 | 0.7716 | -0.0577 | 0.4248 |

The negative/weak clearance R² values are retained rather than hidden and directly constrain
confidence. The named HLM directionality check uses published rifampicin (2.84), isoniazid (13.9),
and ethionamide (77.1 µL/min/mg) values. After explicit human scaling, predictions were 1.9224,
0.4359, and 2.2676 log10(mL/min/kg), respectively, versus experimental 0.5215, 1.2112, and 1.9552;
three-point Spearman was 0.5. Ethionamide remained high, but rifampicin was seriously overpredicted.
This is a sanity check, not an accuracy estimate, and reinforces the LOW clearance confidence.

The Biogen canonical-overlap exclusion cannot rule out close analogue/series overlap. MLM has no
compatible independent endpoint in that set. The OpenADMET released checkpoint includes all
ExpansionRx train and test records; model-card benchmark plots concern an analogous excluded-test
checkpoint, so those plots are not claimed as validation of the packaged all-data checkpoint.
