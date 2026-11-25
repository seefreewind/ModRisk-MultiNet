# Final pair CSV validation

Status: **PASS**

The final Additional file 2 contains 2,485 unique unordered disease pairs. All locked counts, derived-field formulas, medians, denominator rules, 3×3 change-by-residual cells, reversal definitions, equivalence flags and disease-code mappings were independently recalculated from the exported CSV.

## Validation checks

| metric                                                    |       observed |   expected |   tolerance | status   |
|:----------------------------------------------------------|---------------:|-----------:|------------:|:---------|
| rows                                                      | 2485           |  2485      |       0     | PASS     |
| unique_pair_id                                            | 2485           |  2485      |       0     | PASS     |
| baseline_FDR_significant                                  | 1535           |  1535      |       0     | PASS     |
| signed_difference_FDR_significant                         | 1019           |  1019      |       0     | PASS     |
| attenuated                                                |  900           |   900      |       0     | PASS     |
| enhanced                                                  |  119           |   119      |       0     | PASS     |
| no_significant_change                                     | 1466           |  1466      |       0     | PASS     |
| residual_FDR_significant                                  | 1100           |  1100      |       0     | PASS     |
| near_null                                                 |  605           |   605      |       0     | PASS     |
| other_residual_non_significant                            |  780           |   780      |       0     | PASS     |
| sign_reversals                                            |  135           |   135      |       0     | PASS     |
| significant_reversals                                     |   66           |    66      |       0     | PASS     |
| legacy_residual_Z_gate                                    |    1           |     1      |       0     | PASS     |
| CI_within_0.05                                            |    2           |     2      |       0     | PASS     |
| CI_within_0.10                                            |  162           |   162      |       0     | PASS     |
| median_absolute_attenuation                               |    0.0205814   |     0.0206 |       5e-05 | PASS     |
| proportional_attenuation_denominator                      | 1317           |  1317      |       0     | PASS     |
| median_proportional_attenuation                           |    0.249158    |     0.2492 |       5e-05 | PASS     |
| signed_difference_formula_max_error                       |    1.80411e-16 |     0      |       1e-12 | PASS     |
| absolute_attenuation_formula_max_error                    |    1.80411e-16 |     0      |       1e-12 | PASS     |
| residual_CI_lower_formula_max_error                       |    0           |     0      |       1e-12 | PASS     |
| residual_CI_upper_formula_max_error                       |    0           |     0      |       1e-12 | PASS     |
| significant_reversal_formula_matches                      | 2485           |  2485      |       0     | PASS     |
| 3x3_ATTENUATED__RESIDUAL_SIGNIFICANT                      |  495           |   495      |       0     | PASS     |
| 3x3_ATTENUATED__NEAR_NULL_CLASSIFIED                      |  155           |   155      |       0     | PASS     |
| 3x3_ATTENUATED__RESIDUAL_NON_SIGNIFICANT_OTHER            |  250           |   250      |       0     | PASS     |
| 3x3_ENHANCED__RESIDUAL_SIGNIFICANT                        |   47           |    47      |       0     | PASS     |
| 3x3_ENHANCED__NEAR_NULL_CLASSIFIED                        |   26           |    26      |       0     | PASS     |
| 3x3_ENHANCED__RESIDUAL_NON_SIGNIFICANT_OTHER              |   46           |    46      |       0     | PASS     |
| 3x3_NO_SIGNIFICANT_CHANGE__RESIDUAL_SIGNIFICANT           |  558           |   558      |       0     | PASS     |
| 3x3_NO_SIGNIFICANT_CHANGE__NEAR_NULL_CLASSIFIED           |  424           |   424      |       0     | PASS     |
| 3x3_NO_SIGNIFICANT_CHANGE__RESIDUAL_NON_SIGNIFICANT_OTHER |  484           |   484      |       0     | PASS     |

## Disease-code mapping

The mapping table contains 71 analysis codes. Three frozen aliases map to distinct GEMINI resource codes: alc_problems→alcohol_problems, deaf→hearing_loss and schizo→SSD. All other analysis codes retain the same resource code.

Validated file SHA256: `08dd31840323479b53dfaa76ca37217ffc549022328548f97b8039e41c4a2798`.
