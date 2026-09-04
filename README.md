# ModRisk-MultiNet public code and derived-data repository

This repository supports the Genome Medicine manuscript **Modifiable risk factors differentially attenuate the shared genetic architecture of multimorbidity**. It contains the frozen configurations, reproducibility scripts and derived machine-readable tables needed to reproduce the headline numerical checks and inspect the displayed results.

## Scientific overview

The project maps baseline genetic correlations across 71 diseases (2,485 disease pairs), evaluates six modifiable risk axes with single-factor and joint partial genetic correlation, and performs pre-registered BMI + smoking multivariable MR on selected disease pairs. MR-AHC is used to test whether smoking-associated variant heterogeneity forms stable clusters across four outcomes. Phase-4C corrected results are authoritative; the archived negative orientation is excluded.

## Data sources

- GEMINI v1.1 disease GWAS summary statistics: Zenodo record 19890575, subject to source licence.
- Exposure sources: Yengo 2018 BMI; GSCAN 2019 smoking initiation; Neale UK Biobank 4080 SBP; GCST90013994 ApoB; MAGIC Manning 2012 fasting glucose; Dashti 2019 sleep duration.
- LDSC LD reference: 1000 Genomes EUR `eur_w_ld_chr`.
- MR LD reference: MRC-IEU 1000 Genomes phase 3 EUR panel.

Restricted third-party raw GWAS files are intentionally not included and must be obtained from their original repositories under their access conditions.

## Software

The frozen analysis used R 4.4.3, partialLDSC 0.2.0, MVMR 0.4.8, TwoSampleMR 0.7.7 and plink2 2.0.0-a.7.3. Figure rendering uses Python with matplotlib and pandas.

## Reproducing manuscript numbers

From the repository root, run:

```bash
python3 scripts/manuscript/validate_manuscript_numbers.py
```

The script checks the frozen pair counts, FDR counts, primary MVMR family, pair-level classifications and MR-AHC status. It should report 22 assertions and `ALL CHECKS PASSED`.

To render the six draft figures after placing the package in the project environment, run:

```bash
python3 scripts/manuscript/render_figures_gm.py
```

## Directory structure

- `config/` — frozen analysis and terminology configuration.
- `results/` — derived, machine-readable tables used by the manuscript and figures.
- `scripts/manuscript/` — numerical validation and Python figure-rendering scripts.
- `software/` — partialLDSC implementation notes.
Manuscript documents are not included in this repository.

## Expected outputs

The validation script reproduces the headline denominators and counts: 71 diseases, 2,485 pairs, 1,535 baseline FDR-significant pairs, 1,019 significant joint attenuations, 1,100 residual-significant pairs, 30 primary MR effects, BMI 0/15, SMK 15/15, 12/12 concordant SMK pairs and no robust MR-AHC cluster solution.

## Release status and citation

This repository should be cited together with the associated manuscript and, after archival, the versioned Zenodo DOI for the tagged release.
