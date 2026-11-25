# Release licence audit

Status: `PASS`

- The candidate contains project-authored code and derived summary outputs; no full third-party GWAS dataset is redistributed.
- Source resources and access routes are listed in `SOURCE_DATA_MANIFEST.tsv`.
- Instrument-level frequencies, associations and QC fields are included only as derived audit evidence from publicly accessible summary resources. No full source GWAS file is included, and source-specific access routes remain in `SOURCE_DATA_MANIFEST.tsv`.
- External packages are not vendored and retain their own licences.
- Project-authored code, documentation and derived outputs are licensed under CC BY 4.0 (`SPDX-License-Identifier: CC-BY-4.0`).
- Third-party datasets, software and externally authored materials retain their original terms and are not relicensed by the project licence.
- Version v1.0.3 corrects the repository licence metadata to match the CC BY 4.0 rights statement inherited by the Zenodo version chain; no scientific file or result is changed.
- No manuscript, correspondence, participant-level data or restricted raw summary-statistic file is included.
