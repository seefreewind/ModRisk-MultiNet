# Release licence audit

Status: `PASS`

- The candidate contains project-authored code and derived summary outputs; no full third-party GWAS dataset is redistributed.
- Source resources and access routes are listed in `SOURCE_DATA_MANIFEST.tsv`.
- Instrument-level frequencies, associations and QC fields are included only as derived audit evidence from publicly accessible summary resources. No full source GWAS file is included, and source-specific access routes remain in `SOURCE_DATA_MANIFEST.tsv`.
- External packages are not vendored and retain their own licences.
- `LICENSE` reserves all rights because no project-wide open-source licence was supplied. The public v1.0.2 release intentionally retains this notice and does not relicense third-party inputs.
- No manuscript, correspondence, participant-level data or restricted raw summary-statistic file is included.
