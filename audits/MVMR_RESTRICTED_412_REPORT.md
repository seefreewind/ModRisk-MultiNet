# Restricted 412-variant MVMR sensitivity

The final 423-variant matrix contains 11 variants that were not genome-wide significant for BMI or smoking initiation. A bounded sensitivity analysis excluded those 11 variants and repeated only conditional F statistics and covariance-aware MVMR-IVW with the same allele coding, SNP-specific covariance matrices, outcome files, multiplicity family and estimator.

The restricted analysis used 412 variants for 11 outcomes and 411 variants for four outcomes with the existing rs705217 missingness pattern. Conditional F statistics were 25.2470 for BMI and 3.8261 for smoking initiation.

Direction was concordant with the full model for 30/30 exposure-outcome estimates. Nominal significance status was concordant for 29/30 and FDR significance status for 30/30.

## Priority outcomes

Key: <outcome, exposure>
            outcome           exposure full_n  full_beta    full_SE
             <char>             <char>  <int>      <num>      <num>
 1:            COPD                BMI    423  0.2694841 0.03412228
 2:            COPD Smoking initiation    423  1.5972754 0.13555687
 3:              HF                BMI    422  0.5371892 0.02879071
 4:              HF Smoking initiation    422  0.1087939 0.11203245
 5:              OA                BMI    422  0.4133918 0.02263202
 6:              OA Smoking initiation    422  0.1593800 0.08609680
 7:  coronary_heart                BMI    422  0.3182609 0.02953682
 8:  coronary_heart Smoking initiation    422  0.2953722 0.11918685
 9: diabetes_type_2                BMI    422  1.0188487 0.05020025
10: diabetes_type_2 Smoking initiation    422 -0.2132456 0.20045654
          full_P   full_FDR_q restricted_n restricted_beta restricted_SE
           <num>        <num>        <int>           <num>         <num>
 1: 2.498229e-14 6.813353e-14          412       0.2723569    0.03434318
 2: 6.808933e-28 2.553350e-27          412       1.5872839    0.13749752
 3: 5.022793e-57 3.767094e-56          411       0.5376394    0.02841769
 4: 3.320610e-01 3.831473e-01          411       0.0776846    0.11123282
 5: 2.777977e-55 1.388989e-54          411       0.4125643    0.02238278
 6: 6.484658e-02 8.842716e-02          411       0.1359931    0.08596981
 7: 4.560666e-24 1.520222e-23          411       0.3199573    0.02869551
 8: 1.359568e-02 2.039353e-02          411       0.2688976    0.11720428
 9: 2.555115e-64 3.832673e-63          411       1.0070386    0.04736413
10: 2.880308e-01 3.456369e-01          411      -0.1765182    0.19184729
    restricted_P restricted_FDR_q beta_difference absolute_beta_difference
           <num>            <num>           <num>                    <num>
 1: 2.095212e-14     5.714213e-14    0.0028728283             0.0028728283
 2: 6.924068e-27     2.596526e-26   -0.0099915469             0.0099915469
 3: 8.401456e-58     5.040874e-57    0.0004501642             0.0004501642
 4: 4.853262e-01     5.392514e-01   -0.0311093032             0.0311093032
 5: 1.153633e-55     5.768164e-55   -0.0008274892             0.0008274892
 6: 1.144522e-01     1.560712e-01   -0.0233869530             0.0233869530
 7: 2.169163e-25     7.230544e-25    0.0016964268             0.0016964268
 8: 2.228131e-02     3.342196e-02   -0.0264746219             0.0264746219
 9: 4.151319e-68     6.226979e-67   -0.0118100741             0.0118100741
10: 3.580642e-01     4.155225e-01    0.0367273627             0.0367273627
    direction_concordant nominal_status_concordant FDR_status_concordant
                  <lgcl>                    <lgcl>                <lgcl>
 1:                 TRUE                      TRUE                  TRUE
 2:                 TRUE                      TRUE                  TRUE
 3:                 TRUE                      TRUE                  TRUE
 4:                 TRUE                      TRUE                  TRUE
 5:                 TRUE                      TRUE                  TRUE
 6:                 TRUE                      TRUE                  TRUE
 7:                 TRUE                      TRUE                  TRUE
 8:                 TRUE                      TRUE                  TRUE
 9:                 TRUE                      TRUE                  TRUE
10:                 TRUE                      TRUE                  TRUE

This sensitivity analysis does not replace the pre-specified 423-variant model. It isolates the influence of the 11 variants that entered through another exposure axis and is reported as a provenance-focused robustness check.
