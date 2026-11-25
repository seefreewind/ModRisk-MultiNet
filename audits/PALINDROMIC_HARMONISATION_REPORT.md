# Palindromic SNP harmonisation report

## Audit result

All 54 palindromic variants were traced to GRCh37 coordinates, raw BMI and smoking alleles and frequencies, harmonised exposure alleles, and outcome-alignment evidence. Under a conservative ambiguity rule that treats an aligned MAF of 0.42 or greater in either exposure as intermediate-frequency, 37 variants were resolved by EAF and 17 remained frequency-ambiguous.

For 2 of the 54 variants, the opposite strand orientation produced a smaller EAF discrepancy than the retained orientation. Both occur among the frequency-ambiguous variants. The previously reported maximum aligned EAF difference of 0.0383 is a descriptive maximum only; it does not by itself prove strand resolution and must not be used as a blanket harmonisation claim.

Outcome associations are aligned to the smoking effect allele in the frozen corrected matrix. This verifies internal outcome coding but cannot independently resolve a palindromic exposure strand when the exposure EAFs are near 0.5.

## Decision rule

A palindrome was labelled RESOLVED_BY_EAF only when both aligned exposure frequencies were informative (MAF < 0.42), the retained orientation had a smaller EAF discrepancy than the opposite orientation, and the retained-orientation difference was no greater than 0.10. All other palindromes were labelled UNRESOLVED_FREQUENCY_AMBIGUOUS.

## Manuscript consequence

The sentence claiming that all 54 palindromic variants were resolved by aligned allele frequencies should be replaced. The revised text should report that 54 palindromic variants were audited, 37 had frequency support for the retained orientation, and 17 remained frequency-ambiguous. The unsupported 0.0383 resolution claim should be removed. The ambiguity must be disclosed as a limitation and the final manifest must carry the per-SNP decisions.
