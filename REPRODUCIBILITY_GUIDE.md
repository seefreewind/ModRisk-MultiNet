# Reproducibility guide

1. Obtain the disease and exposure GWAS inputs listed in `SOURCE_DATA_MANIFEST.tsv`; third-party data are not bundled.
2. Recreate a compatible Python/R environment using `SOFTWARE_ENVIRONMENT.tsv`.
3. Set `MODRISK_ROOT` to the project directory containing the source-data layout referenced by the scripts. If R packages are installed in a project-specific library, set `MODRISK_R_LIB` to that directory; otherwise the normal R library paths are used.
4. Run the audit and analysis scripts in `code/` in their documented order: network audit, pair-classification repair, allele harmonisation, independent UVMR, robust MVMR, table generation and figure rendering.
5. Compare generated outputs with `derived_results/` and verify `SHA256SUMS.txt`.
6. Do not interpret files labelled `INVALID_ARCHIVED_WITHDRAWN` as current evidence.

The final analysis uses the allele-harmonised two-exposure BMI-smoking MVMR only. Exact historical composite-score cut-offs for exploratory pair selection and the historical MR-AHC bootstrap seed/count were not recoverable and are disclosed as unavailable.
