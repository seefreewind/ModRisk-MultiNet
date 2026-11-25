# Smoking effective-sample-size reproducibility audit

## Frozen source

- Source: `ModRisk-MultiNet-data/raw/exposures/Smoking/SmokingInitiation.txt.gz`
- SHA-256: `cab4a185e3ab63989e9f6b3a3cc64c4af920c5f1e2e880907802e9f8a817b3b9`
- Munged LDSC file: `ModRisk-MultiNet-data/munged/SMK_GSCAN2019.sumstats.gz`
- Munged-file SHA-256: `0453313b70806f4e58ca2c31acfb3a85929a2b29dc027e4d6ac83ac2cca3f4ea`

## Recalculation

The frozen munging pipeline removed strand-ambiguous variants, deduplicated rsIDs, and projected the source to HapMap3. The resulting LDSC input contains 1,185,282 variants. Recalculating the arithmetic mean of its per-variant `N` field gives 619,093.3457, reported as 619,093 in the manuscript.

For comparison, the unfiltered public source contains 11,733,344 rows and has a mean `EFFECTIVE_N` of 555,845.4618. The previously reported value 620,841 is not reproducible from either the complete downloaded file or the frozen LDSC input and has therefore been corrected. The maximum cohort-sum `N` in the public source remains 632,802.

## Decision

Use **mean per-SNP effective n = 619,093 across the 1,185,282 variants retained for LDSC**. Retain **maximum cohort-sum n = 632,802** as a distinct source-file descriptor.
