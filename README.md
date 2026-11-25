# ModRisk-MultiNet reproducibility archive v1.0.3

This public release contains project-authored code and derived, non-restricted outputs for the BMC Genomics major revision of *Genome-wide mapping of shared and residual genetic architecture of multimorbidity after conditioning on six modifiable risk-factor traits*.

The current package adds a prespecified frequency audit of all 54 palindromic variants. Thirty-seven variants were resolved using effect-allele frequency and 17 remained frequency-ambiguous. Those 17 variants were excluded before the strict-primary MVMR analysis, leaving 406 variants (405 for four outcomes in which rs705217 was unavailable). The archive also includes complete instrument-set A–D comparisons, strict-primary IVW and robust estimators, covariance-aware heterogeneity statistics, 6,086 leave-one-SNP-out fits, method-specific confidence intervals, the final 2,485-pair classification file and the complete historical 54-variant MR-AHC provenance. MR-AHC remains withdrawn and is included only as provenance; no clustering result is treated as current evidence.

The archive does not redistribute third-party GWAS summary-statistic files. Source data must be obtained from the resources in `SOURCE_DATA_MANIFEST.tsv` and placed in local paths before code is rerun. Source-file paths and checksums in the audit tables identify the inputs used; they do not grant redistribution rights. The 17 frequency-ambiguous palindromic variants are retained only in the disclosed historical 423-variant sensitivity model and are excluded from the strict-primary model.

Manuscripts, response letters and submission documents are excluded. Version v1.0.2 is archived at DOI [10.5281/zenodo.24528817](https://doi.org/10.5281/zenodo.24528817). Version v1.0.3 changes only the project licence to CC BY 4.0 so that the GitHub release and Zenodo rights metadata are consistent; scientific code, data and results are unchanged.

## Licence

Project-authored code, documentation and derived outputs are available under the [Creative Commons Attribution 4.0 International License](LICENSE). Third-party datasets and software remain subject to their original terms and are not relicensed here.

## Key entry points

- `derived_results/PAIR_LEVEL_CLASSIFICATION.csv`: final pair-level data.
- `derived_results/UVMR_FULL_RESULTS.tsv`: independent outcome-specific UVMR estimates.
- `derived_results/MVMR_STRICT_FULL_RESULTS.tsv`: strict-primary IVW, weighted-LAD median, Q-minimisation and directional MVMR-Egger results.
- `derived_results/MVMR_INSTRUMENT_SET_COMPARISON.tsv`: predefined instrument-set A–D comparison.
- `derived_results/MVMR_STRICT_HETEROGENEITY.tsv` and `MVMR_STRICT_LEAVE_ONE_OUT_*.tsv`: covariance-aware Q and influence diagnostics.
- `derived_results/MVMR_STRICT_INSTRUMENT_MANIFEST.tsv`: strict-primary set membership and per-outcome availability.
- `derived_results/MVMR_423_INSTRUMENT_PROVENANCE.tsv`: complete 423-variant source chain.
- `derived_results/PALINDROMIC_SNP_EAF_AUDIT.tsv`: frequency and orientation audit for all 54 palindromic variants.
- `derived_results/MVMR_RESTRICTED_412_*.tsv`: historical implementation of the BMI-or-smoking genome-wide-significant sensitivity set; the current A–D comparison files are authoritative.
- `code/`: analysis, covariance, conditional-F, CI validation and figure scripts.
- `audits/`: scientific qualification and withdrawal records.
- `references/`: reference metadata, claim-support and first-citation-order audits.
