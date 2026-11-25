# MR confidence-interval correction report

## Scope and result

All 240 UVMR effect rows, all 150 MVMR effect rows, and 30 unique directional MVMR-Egger intercepts were audited against the probability distribution used to calculate each reported P value. Beta estimates, standard errors, P values and FDR values were unchanged.

UVMR IVW and weighted-median P values use the standard normal distribution. UVMR MR-Egger uses a Student t distribution with n−2 degrees of freedom. Weighted mode uses a Student t distribution with n−1 degrees of freedom. The historical export applied 1.96 to every method, causing four threshold inconsistencies. Corrected method-specific 95% confidence intervals remove all four.

MVMR-IVW uses Student t with n−2 degrees of freedom, and directional MVMR-Egger effects and intercepts use Student t with n−3 degrees of freedom. Multivariable weighted LAD median and Q-minimisation use normal-reference intervals. No MVMR row changed nominal significance classification after the CI correction.

## Four corrected UVMR threshold inconsistencies

| instrument_set            | method        | exposure           | outcome              |   n_instruments |   reported_P |   reported_CI_lower |   reported_CI_upper |   corrected_CI_lower |   corrected_CI_upper |
|:--------------------------|:--------------|:-------------------|:---------------------|----------------:|-------------:|--------------------:|--------------------:|---------------------:|---------------------:|
| UNION_SUBSET_UVMR         | MR-Egger      | Smoking initiation | coronary_heart       |              16 |    0.0538867 |          0.086002   |             2.41975 |         -0.0240095   |              2.52976 |
| UNION_SUBSET_UVMR         | Weighted mode | Smoking initiation | gastritis_duodenitis |              16 |    0.0671944 |          0.00656701 |             1.96077 |         -0.0789041   |              2.04624 |
| UNION_SUBSET_UVMR         | MR-Egger      | Smoking initiation | gout                 |              16 |    0.0639267 |          0.0412377  |             3.17716 |         -0.106588    |              3.32499 |
| INDEPENDENT_EXPOSURE_UVMR | MR-Egger      | Smoking initiation | coronary_heart       |              84 |    0.0502414 |          0.00996171 |             1.44683 |         -0.000784972 |              1.45758 |

Rows with corrected numeric CI endpoints: 390. Rows with a P/CI significance mismatch after correction: 0.
