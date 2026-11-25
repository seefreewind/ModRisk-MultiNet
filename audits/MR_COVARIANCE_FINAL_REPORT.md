# MR covariance final report

## Result

The final 423 SNP-specific covariance matrices are symmetric and positive semidefinite; the smallest eigenvalue is 2.45641e-06.
The final lambda is 0.0915189134758081 from 159637 null SNPs.

## Allele-sign audit

The final analysis is a joint common-allele recoding. Relative to common BMI coding, both exposure associations and the outcome association are multiplied by -1 for swapped rows when represented on the common smoking allele. This is S=-I for the two exposure effects, so S Sigma S equals Sigma. A single-exposure-only flip would instead require the covariance off-diagonal to change sign; that counterexample is not the final analysis.

## Conditional F reconciliation

Using identical final data and covariance, MVMR 0.4.8 gives 25.41967/3.843062 and MVMR 0.4.9 gives 25.47990/3.852169 for BMI/smoking initiation.
The small change is therefore reproduced as a package-version effect, with no change in strength classification. Four outcomes use 422 SNPs for effect estimation because rs705217 lacks an outcome association; the reported conditional F values describe the complete 423-SNP exposure model and are not outcome-specific causal estimates.

## Boundary

Exposure-exposure sampling covariance does not resolve unknown exposure-outcome sample overlap. Smoking initiation remains conditionally weak.
