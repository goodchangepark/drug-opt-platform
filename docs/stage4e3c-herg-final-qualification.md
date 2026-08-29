# Stage 4E-3C — hERG Final Qualification & Closure

This was the final dedicated hERG optimization attempt before Engine-v1
freeze. The exact production policy remains raw Admetica M1,
`admetica_safety_herg` version `admetica-d4f7056-herg-chemprop-v2.1`, at
threshold 0.50. The existing PhysChem model remains calibration-supporting
shadow-only. No threshold, ensemble, calibration, or runtime policy changed.

The authoritative ChEMBL37 overlap-excluded cohort has N=728 (489 positive,
239 negative; prevalence 0.6717). Raw M1 achieved MCC 0.1844, balanced
accuracy 0.5442, AUROC 0.6669, sensitivity 0.9755, specificity 0.113,
Brier 0.2745, LogLoss 1.6901, and ECE 0.2651. This is usable but low
discrimination with severe class/prior and assay-heterogeneity limitations.

The historical Platt audit used a scaffold-aware 546/182 calibration/holdout
split. On the untouched holdout, Platt improved Brier 0.2763→0.2144,
LogLoss 1.6942→0.6259, and ECE 0.272→0.0887. AUROC remained essentially
unchanged (0.6001→0.599), so this is probability-calibration evidence, not
new discrimination. Because it is one historical holdout with no new
prospective calibration cohort and threshold semantics were not independently
qualified, the final state is `CALIBRATION_RESEARCH_ONLY`.

Three serious primary-source secondary candidates were screened: AttenhERG,
NCATS herg-ml, and ADMET-AI v2. AttenhERG exposes repository model files but
does not separately establish weight/data rights or a versioned checksum.
NCATS requires legacy Git-LFS/conda assets and explicitly warns that the
archived repository should not be used in production. ADMET-AI v2 has likely
TDC lineage overlap and unresolved independent checkpoint/ARM64 gates. No
candidate was downloaded, benchmarked, or registered.

Final closure is `HERG_NO_QUALIFIED_REPLACEMENT_RAW_M1_FROZEN`. hERG
optimization is CLOSED after Stage 4E-3C for Engine-v1. No further dedicated
hERG optimization is planned before Engine-v1 freeze. Future work proceeds
only through the governed next endpoint stage and later final policy review.
