# R2 pair-classification rules — CLASSIFICATION_RULES_FROZEN_R2

The frozen machine-readable change classes are `ATTENUATED` (900), `ENHANCED` (119), and `NO_SIGNIFICANT_CHANGE` (1,466). The last code is displayed in prose as **no FDR-significant signed change** (`NO_FDR_SIGNIFICANT_SIGNED_CHANGE`); it does not assert a zero effect. Difference FDR is evaluated on signed `baseline rg − joint rg` in a 2,485-test family. Significant pairs with decreased `|rg|` are attenuated; those with increased `|rg|` are enhanced.

Residual states are independently `RESIDUAL_SIGNIFICANT` (1,100, joint residual FDR<0.05), `NEAR_NULL_CLASSIFIED` (605, residual FDR≥0.05, `|rg|<0.05`, `|Z|<2`), and `RESIDUAL_NON_SIGNIFICANT_OTHER` (780). Residual FDR uses a separate 2,485-test family. The cross-tab of residual-significant pairs is 495 attenuated, 47 enhanced, and 558 without FDR-significant signed change.

Sign reversal is a separate flag (135 pairs); 66 had significant signed differences. The original protocol's stricter significant-reversal flag additionally requires residual `|Z|≥2` (1 pair). CI-based residual equivalence requires the whole 95% CI within ±0.05 (2 pairs); ±0.10 is a sensitivity margin (162 pairs). Near-null is descriptive, not equivalence.

The historical `machine_class` in frozen estimate tables is a **deprecated legacy classification**, not the R2 scientific result.
