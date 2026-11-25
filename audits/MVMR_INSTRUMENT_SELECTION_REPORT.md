# MVMR instrument selection report

## Reconstructed selection chain

The final 423-variant matrix was reconstructed exactly from the historical six-exposure genome-wide-significant union. The union was jointly LD-clumped, then restricted to variants with complete BMI and smoking-initiation associations and a per-variant F statistic above 10 for at least one of the two exposures. Complete BMI and smoking data alone yielded 481 variants; the F-statistic gate yielded the final 423.

The previous shorthand description as a BMI-smoking union was incorrect. The correct description is a six-exposure clumped union followed by BMI-smoking completeness and the recorded per-variant strength gate.

## Counts

| metric                 |   observed |   expected | status   |
|:-----------------------|-----------:|-----------:|:---------|
| BMI_GWS                |        404 |        404 | PASS     |
| SMK_GWS                |         16 |         16 | PASS     |
| both_GWS               |          8 |          8 | PASS     |
| either_GWS             |        412 |        412 | PASS     |
| neither_GWS            |         11 |         11 | PASS     |
| complete_before_F_gate |        481 |        481 | PASS     |
| final                  |        423 |        423 | PASS     |

Among the final 423 variants, 404 were genome-wide significant for BMI, 16 for smoking initiation, and 8 for both. Their union contained 412 variants. Eleven variants were not genome-wide significant for either BMI or smoking initiation and entered through another exposure axis before satisfying the BMI-or-smoking F-statistic gate. These 11 variants are listed separately and are not silently labelled BMI-smoking genome-wide-significant instruments.

## Source checksums

- `results/phase4/MR_INSTRUMENT_UNION_CLUMPED.tsv`: `5c9fb3d201f4751c68d4235802438af3f3559b9638e67f63e39d1d6360ab5708`
- `results/phase4/MVMR_EXPOSURE_MATRIX.tsv`: `ac1a8c485f029687b806eec86923b9323af3150ad30f9b9ab9c2c5b5ce58b06b`
- `results/phase4_repair/PRIMARY_MVMR_CORRECTED_INPUT.tsv.gz`: `1ded94f2734a060cb4676f5653b821b5cc9f69b75056b227b5744c576436e03b`

## Methods wording

Candidate instruments were formed as the union of genome-wide-significant variants across the six exposure GWASs and jointly LD-clumped at r² < 0.001 within 10 Mb. The BMI-smoking MVMR matrix then retained variants with complete associations for both exposures and a per-variant F statistic above 10 for at least one exposure, yielding 423 variants. Of these, 412 were genome-wide significant for BMI or smoking initiation and 11 entered through another exposure axis.
