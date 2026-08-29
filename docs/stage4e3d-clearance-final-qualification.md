# Stage 4E-3D — Clearance Final Qualification

HLM, RLM and MLM remain species-isolated `SINGLE_CORE_MODEL` predictions using the OpenADMET CheMeleon checkpoint `openadmet-microsomal-clearance-chemeleon-v1-e135493`. Contracts are scaled `log10(mL/min/kg)` and require explicit microsomal assay context.

The pinned ExpansionRx file is structurally available (7,618 rows; HLM/RLM/MLM fields), but the current checkpoint training lineage explicitly includes ExpansionRx. It is therefore not an independent benchmark. The previously referenced Biogen prospective cohort has no reproducible raw licensed structures and values in the available source chain. No new inference, bias correction, calibration, species averaging, or production mutation was performed.

Historical reported metrics are retained as non-reproducible metadata: HLM N=3078, MAE 0.6259, RMSE 0.7616, Spearman 0.3700; RLM N=3045, MAE 0.6263, RMSE 0.7716, Spearman 0.4248. MLM has no compatible independent endpoint. Bias, AD and scaffold metrics cannot be freshly assessed without an independent raw cohort. Dog, monkey and generic microsomal clearance remain unavailable.

Clearance dedicated optimization is CLOSED for Engine v1. Production policy is unchanged; limitations and residual overlap uncertainty are preserved for future review.
