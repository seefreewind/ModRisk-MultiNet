# ModRisk-MultiNet: public code and derived data

This repository supports the BMC Genomics research article **Genome-wide mapping of shared and residual genetic architecture of multimorbidity after conditioning on six modifiable risk-factor traits**. Manuscript documents and restricted third-party GWAS files are not redistributed here.

The fixed network contains 71 diseases and 2,485 disease pairs. Baseline genetic correlations were FDR-significant for 1,535 pairs. Joint conditioning produced FDR-significant signed baseline-to-joint differences for **1,019** pairs: **900** showed decreased and **119** showed increased absolute genetic-correlation magnitude. Another 1,466 pairs had no FDR-significant signed difference. Across the network, **1,100** pairs retained FDR-significant residual genetic correlations. Change class and residual significance are independent dimensions; 495 attenuated, 47 enhanced, and 558 no-significant-signed-change pairs retained a significant residual correlation.

The R2 repair changed **derived pair-level labels**, not GWAS, LDSC, partialLDSC, joint-conditioning, MVMR, or MR-AHC estimates. See [the correction note](R2_CLASSIFICATION_CORRECTION_NOTE.md), [classification rules](config/CLASSIFICATION_RULES_R2.md), and [data dictionary](DATA_DICTIONARY_R2.md). The machine-readable R2 pair table and cross-tab are in `results/phase3_classification_repair_R2/`.

## Reproducibility

Run `python3 scripts/manuscript/validate_manuscript_numbers.py` from the repository root to verify the headline R2 counts and arithmetic. `scripts/manuscript/render_bmc_genomics_r2_figures.py` is the frozen Figure 3/4 R2 renderer; `scripts/manuscript/render_bmc_genomics_submission_figures34.py` changes only the Figure 3 display label to “No FDR-significant signed change” for submission. `scripts/manuscript/render_bmc_genomics_submission_figure4_compact.py` stacks the unchanged R2 Figure 4 panels at journal page width. `scripts/manuscript/render_bmc_genomics_submission_figure2.py` corrects Figure 2's wording without changing its values or disease order. `scripts/manuscript/render_figures_gm.py` is retained only as a **deprecated legacy classification** renderer and must not be used for current manuscript figures.

The disease GWAS data are from GEMINI v1.1 ([source archive](https://doi.org/10.5281/zenodo.19890575)); exposure sources and LD resources are specified in the frozen configurations. Source licences govern access to third-party GWAS files.

## Release status

The corrected R2 code and derived data are archived in [Zenodo version `v1.0.1`, DOI 10.5281/zenodo.22927913](https://doi.org/10.5281/zenodo.22927913), corresponding to the [GitHub `v1.0.1` release](https://github.com/seefreewind/ModRisk-MultiNet/releases/tag/v1.0.1). The previously published [Zenodo record 10.5281/zenodo.22297319](https://doi.org/10.5281/zenodo.22297319) corresponds to the historical GitHub `1` tag and **does not contain this R2 correction**. The GitHub `v1.0.0` release is also historical. Cite the version-specific R2 DOI for the corrected classification.
