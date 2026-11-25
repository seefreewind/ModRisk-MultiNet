# Changelog

## v1.0.2 — released (2026-11-25)

- Corrected method-specific confidence intervals without changing beta, standard error, P value or FDR q value.
- Added complete provenance for the 423-variant six-exposure-union MVMR matrix and documented the 11 variants not genome-wide significant for BMI or smoking initiation.
- Added the restricted 412-variant MVMR-IVW sensitivity analysis.
- Added an independent frequency audit for all 54 palindromic variants: 37 frequency-supported and 17 frequency-ambiguous.
- Excluded the 17 frequency-ambiguous palindromic variants from the strict-primary MVMR set before outcome modelling, leaving 406 variants (405 for four outcomes).
- Added complete results for strict-primary IVW, weighted-LAD median, Q-minimisation and two directional MVMR-Egger codings.
- Added instrument-set A–D comparisons, covariance-aware Q statistics and 6,086 successful leave-one-SNP-out fits.
- Reproduced the smoking-initiation effective sample size as 619,093 from the frozen LDSC input and corrected the unreproducible historical value of 620,841.
- Added an explicit baseline-by-residual significance cross-tabulation (1,079 both significant, 456 baseline-only, 21 residual-only and 929 neither).
- Replaced the unsupported claim that all palindromic variants were resolved and removed the unsupported maximum-EAF-difference statement.
- Added the final 2,485-pair classification file, disease-code mapping and validation report.
- Expanded the historical withdrawn MR-AHC provenance to the complete 54-variant set without reinstating MR-AHC results.
- Preserved the frozen primary network estimates while updating the exploratory MR layer to the frequency-audited strict-primary analysis.

## v1.0.1 — published historical archive

Published at DOI 10.5281/zenodo.22927913. Superseded by v1.0.2 only after the new version is publicly released.
