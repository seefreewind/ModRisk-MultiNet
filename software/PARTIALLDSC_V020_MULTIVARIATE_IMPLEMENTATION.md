# partialLDSC v0.2.0 — MULTIVARIATE IMPLEMENTATION AUDIT

**Date:** 2026-08-28
**Source:** `.work/partialLDSC_v020/R/partial_ldsc.R` (804 lines)
**Version:** 0.2.0 (DESCRIPTION). git metadata damaged (`.git` empty dir from storage relocation);
source tree contents preserved byte-identically (relocation integrity 87/87 PASS).

## 1. Multivariate function

`partial_ldsc(conditions, confounders, condition.names, confounder.names, ld, n.blocks, log.name)`

- **confounders**: character vector of munged GWAS paths — **supports ≥1 confounders natively**.
  (lines 66-70: `confounders = normalizePath(confounders)`; `n.confounders = length(confounders)`)
- All traits (conditions + confounders) enter a single joint LDSC covariance estimate `S` (dimension 72 for 71 conditions + 1 confounder; dimension 77 for 71 conditions + 6 confounders).
- Multi-confounder capability confirmed: partial block formulas index `(n.conditions+1):n.traits` as the confounder block.

## 2. Schur-complement implementation (point estimate)

Lines 583-592:

```r
partial.cov[j,k] <- cov[j,k] - t(cov[j, confounder_cols]) %*%
                    solve(cov[confounder_cols, confounder_cols]) %*%
                    t(t(cov[k, confounder_cols]))
```

This is exactly `Sigma_D|E = Sigma_DD − Sigma_DE solve(Sigma_EE) Sigma_ED` — the Schur complement.

## 3. Matrix inversion

- Plain **`solve()`** (LAPACK) for both the LDSC regression (line 323, 335, 461) and the confounder-block inversion (line 588, 629).
- **NO** `nearPD`, **NO** ridge, **NO** pseudoinverse/`ginv`/`pinv`, **NO** matrix smoothing, **NO** eigen-repair. Confirmed by grep over full source.
- Consequence: if `Sigma_EE` is rank-deficient the package will error or produce Inf/NaN — it does NOT silently repair. This is the intended fail-open behavior for validation.

## 4. Block-jackknife uncertainty propagation

Two-layer propagation, both via leave-one-out pseudo-values:

- **Unadjusted S/V**: per block, `solve(xtx.delete, xty.delete)` → `delete.values` → pseudo-values `V.hold = n.blocks*reg − (n.blocks−1)*delete.values` (lines 335-343); `jackknife.cov = cov(pseudo)/n.blocks` (line 343).
- **Partial cov V**: the SAME Schur formula is applied per block to the delete-values block (`partial.V.delete`, lines 626-630), then pseudo-values `partial.V.hold` (line 632), then `partial.V = cov(partial.V.hold)/n.blocks` (line 641).
- **Partial rg V**: per-block partial rg `partial.V_Stand.delete = partial.V.delete[,vjk]/sqrt(partial.V.delete[,vjj]*partial.V.delete[,vkk])` (line 671), pseudo-values line 672, `partial.V_Stand = cov(...)/n.blocks` (line 678). SE = `sqrt(diag(partial.V_Stand))` (line 684).
- **rg–partial_rg covariance** (needed for diff test): `cov(V_Stand.hold[,vraw], partial.V_Stand.hold[,vpartial])/n.blocks` (line 721).

**No naive delta approximation:** SEs come from jackknife pseudo-value replicates, not from a delta-method linearization of the Schur complement.

## 5. Standardisation

- Unadjusted: `S_Stand = S * tcrossprod(1/sqrt(diag(S)))` (standard cov2cor form).
- Partial: `partial.S_Stand = partial.cov * tcrossprod(1/sqrt(diag(partial.cov)))` (line 647-648), guarded by `all(diag(partial.cov) > 0)` (line 644) — if any partial variance ≤ 0, NO standardized results are produced (fail-open, matches Phase-3A TASK 8 requirement to flag invalid conditional variances).

## 6. Difference test

Lines 771-773:

```r
diff.T = (rg − partial_rg) / sqrt(rg.SE² + partial_rg.SE² − 2·rg_cov)
diff.P = 2·pnorm(−|diff.T|)
```

`rg_cov` is the jackknife covariance between raw and partial rg (line 721).

## 7. Returned objects

`list(res_diff, S, V, S_Stand, V_Stand, partial.S, partial.V, partial.S_Stand, partial.V_Stand, I)`

## 8. Validation-relevant properties

| Property | Finding |
|---|---|
| Multi-confounder input | YES (`confounders` vector) |
| Schur complement | YES (line 586) |
| Inversion method | `solve()` only |
| nearPD / ridge / pseudo-inverse / smoothing | NONE |
| Jackknife through partial step | YES (pseudo-values) |
| Jackknife through rg standardisation | YES |
| Fail on non-PD confounder block | implicit (solve error) — no silent repair |
| Fail on partial variance ≤ 0 | yes (`all(diag(partial.cov)>0)` guard) |
| Delta-method SEs | NOT used |

## 9. Conclusion

The package implements the exact Schur complement with full block-jackknife propagation,
no silent matrix repair. Independent reference implementation (TASK 3) can therefore be
compared against `partial.S`, `partial.V`, `partial.S_Stand`, `partial.V_Stand` directly.