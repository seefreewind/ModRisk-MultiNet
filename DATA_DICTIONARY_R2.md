# R2 derived data dictionary

The current pair-level dataset is `results/phase3_classification_repair_R2/PAIR_CLASSIFICATION_REBUILT.tsv` (2,485 unordered disease pairs). `CHANGE_BY_RESIDUAL_CROSSTAB.tsv` is its two-dimensional summary. These files are derived from frozen estimates; the R2 repair did not rerun genetic models. Restricted third-party GWAS summary statistics are not redistributed.

| Field(s) | Meaning |
|---|---|
| `pair_id`, `trait_1`, `trait_2` | Unique pair identifier and the two disease IDs. Disease labels are mapped in `config/diseases_gemini.tsv` and `config/FIGURE_DISEASE_LABEL_MAP.tsv`. |
| `baseline_rg`, `baseline_se`, `baseline_p`, `baseline_fdr` | Baseline signed genetic correlation, jackknife standard error, nominal P value and FDR. |
| `joint_rg`, `joint_se`, `joint_p`, `residual_fdr` | Six-trait jointly conditioned residual genetic correlation and its uncertainty/significance. |
| `signed_difference`, `difference_se`, `difference_p`, `difference_fdr` | Signed `baseline_rg − joint_rg` difference, uncertainty and inferential test; FDR across 2,485 pair tests. |
| `baseline_abs_rg`, `joint_abs_rg`, `absolute_attenuation` | Absolute correlation magnitudes and their signed difference `|baseline_rg| − |joint_rg|`; positive is lower magnitude, negative is higher magnitude. |
| `proportional_attenuation` | Descriptive ratio where the pre-specified denominator is valid; missing values are not zero. |
| `residual_ci_lower`, `residual_ci_upper` | 95% interval for residual `rg`. |
| `change_class` | `ATTENUATED` (900), `ENHANCED` (119), or stored code `NO_SIGNIFICANT_CHANGE` (1,466). The last code means **no FDR-significant signed difference**, not equality. |
| `residual_state`, `residual_significant`, `near_null_classified` | Independent residual classification: 1,100 significant, 605 descriptive near-null, 780 residual non-significant other. |
| `sign_reversal_flag`, `significant_reversal_flag` | Separate reversal indicators; the stricter flag additionally requires significant signed difference and residual `|Z|≥2`. |
| `equivalent_005`, `equivalent_010` | Entire residual 95% CI inside ±0.05 or ±0.10. Near-null classification does not imply equivalence. |
| `original_frozen_class`, `original_near_null_flag` | **Deprecated legacy classification**, retained only to audit the prior derived labels. Do not use in R2 figures or results. |

Difference FDR and residual FDR are distinct Benjamini–Hochberg families. See `config/CLASSIFICATION_RULES_R2.md` for the exact decision rules.

Figure mapping: Figure 2 uses `results/phase3_joint/FIGURE2_HEATMAP_DATA.tsv`, `DISEASE_LEVEL_JOINT_SUMMARY.tsv`, R2 change counts and `results/phase2_single_factor/NETWORK_ATTENUATION_SUMMARY.tsv`; Figure 3 uses the rebuilt pair table; Figure 4 uses that table and `figures/bmc_genomics_r2/FIGURE4_R2_DISEASE_SUMMARY.tsv`. The current figure scripts are `scripts/manuscript/render_bmc_genomics_submission_figure2.py`, `render_bmc_genomics_submission_figures34.py`, and `render_bmc_genomics_submission_figure4_compact.py`.

No repository-wide reuse licence has been assigned in this release. Source-provider terms continue to apply to third-party data.
