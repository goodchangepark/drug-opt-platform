# Stage 4E-3E — pKa / logD Final Qualification & Closure

The production pKa implementation remains `RULE_ESTIMATE` (`ionization_smarts_rules_v1`, `stage4c4-ionization-v1`). It identifies acidic/basic SMARTS centers, atom indices, multiple-center and ampholyte/zwitterion classes, and representative literature estimates; it is not quantitatively ML validated.

The production logD pH 7.4 implementation remains `DERIVED_ESTIMATE` (`henderson_hasselbalch_logd_v1`), combining RDKit Crippen cLogP with rule/experimental pKa under simplified monoprotic Henderson–Hasselbalch assumptions. cLogP is never relabeled logD.

pkasolver-lite was not accepted: isolated ARM64 dependencies built, but strict checkpoint loading failed because the legacy PyG state-dict layout is incompatible. Partial loading, ignored keys, and random initialization are prohibited. pKaLearn has unresolved checkpoint and weight/data lineage. The pH 7.4 logD dataset lead remains license-unresolved and no quantitative logD checkpoint was qualified.

Final decisions: pKa `PKA_NO_REPRODUCIBLE_QUANTITATIVE_MODEL_RULE_ESTIMATE_FROZEN`; logD `LOGD_NO_QUALIFIED_QUANTITATIVE_MODEL_DERIVED_ESTIMATE_FROZEN`. No production code, policy, prediction, or database records changed. pKa and logD dedicated optimization is CLOSED for Engine v1.
