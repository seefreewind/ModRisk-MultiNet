#!/usr/bin/env Rscript
r_lib <- Sys.getenv("MODRISK_R_LIB", unset = "")
if (nzchar(r_lib)) .libPaths(c(r_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(MVMR)
})

root <- normalizePath(Sys.getenv("MODRISK_ROOT", unset = "."), mustWork = TRUE)
out <- file.path(root, "pre_submission_technical_closeout", "01_instrument_provenance")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

dat <- fread(file.path(root, "results/phase4_repair/PRIMARY_MVMR_CORRECTED_INPUT.tsv.gz"))
pvals <- fread(file.path(root, "results/phase4/MVMR_EXPOSURE_MATRIX.tsv"), select = c("SNP", "BMI_P", "SMK_P"))
dat <- merge(dat, pvals, by = "SNP", all.x = TRUE, sort = FALSE)
dat <- dat[match(fread(file.path(root, "results/phase4_repair/PRIMARY_MVMR_CORRECTED_INPUT.tsv.gz"), select = "SNP")$SNP, SNP)]
gencov <- readRDS(file.path(root, "results/phase4_repair2/PRIMARY_BMI_SMK_GENCOV_LAMBDA.rds"))
gencov <- gencov[match(dat$SNP, names(gencov))]

keep <- (is.finite(dat$BMI_P) & dat$BMI_P < 5e-8) | (is.finite(dat$SMK_P) & dat$SMK_P < 5e-8)
stopifnot(sum(keep) == 412L)
dat <- dat[keep]
gencov <- gencov[keep]

same <- dat$BMI_EA == dat$SMK_EA & dat$BMI_OA == dat$SMK_OA
swap <- dat$BMI_EA == dat$SMK_OA & dat$BMI_OA == dat$SMK_EA
stopifnot(all(same | swap))
dat[swap, BMI_beta := -BMI_beta]

outcomes <- c(
  "COPD", "HF", "OA", "alcohol_problems", "carpal_tunnel", "coronary_heart",
  "diabetes_type_2", "fibromyalgia", "gastritis_duodenitis", "gout",
  "hypertension", "insomnia", "periph_neuro", "sleep_apnoea", "valve_nonrheumatic"
)

make_input <- function(o) {
  by <- dat[[paste0(o, "_beta")]]
  sy <- dat[[paste0(o, "_SE")]]
  ok <- is.finite(dat$BMI_beta) & is.finite(dat$BMI_SE) & dat$BMI_SE > 0 &
    is.finite(dat$SMK_beta) & is.finite(dat$SMK_SE) & dat$SMK_SE > 0 &
    is.finite(by) & is.finite(sy) & sy > 0
  list(
    input = format_mvmr(
      cbind(dat$BMI_beta[ok], dat$SMK_beta[ok]), by[ok],
      cbind(dat$BMI_SE[ok], dat$SMK_SE[ok]), sy[ok], dat$SNP[ok]
    ),
    gencov = gencov[ok],
    ok = ok
  )
}

run_one <- function(o) {
  z <- make_input(o)
  invisible(capture.output(fit <- ivw_mvmr(z$input, gencov = z$gencov)))
  n <- nrow(z$input)
  crit <- qt(0.975, df = n - 2L)
  data.table(
    outcome = o, n_instruments = n,
    BMI_beta = fit["exposure1", "Estimate"],
    BMI_SE = fit["exposure1", "Std. Error"],
    BMI_P = fit["exposure1", "Pr(>|t|)"],
    BMI_CI_lower = fit["exposure1", "Estimate"] - crit * fit["exposure1", "Std. Error"],
    BMI_CI_upper = fit["exposure1", "Estimate"] + crit * fit["exposure1", "Std. Error"],
    SMK_beta = fit["exposure2", "Estimate"],
    SMK_SE = fit["exposure2", "Std. Error"],
    SMK_P = fit["exposure2", "Pr(>|t|)"],
    SMK_CI_lower = fit["exposure2", "Estimate"] - crit * fit["exposure2", "Std. Error"],
    SMK_CI_upper = fit["exposure2", "Estimate"] + crit * fit["exposure2", "Std. Error"],
    distribution = "Student_t", df = n - 2L
  )
}

res <- rbindlist(lapply(outcomes, run_one))
q <- p.adjust(c(res$BMI_P, res$SMK_P), method = "BH")
res[, `:=`(BMI_FDR_q = q[seq_len(.N)], SMK_FDR_q = q[.N + seq_len(.N)])]
fwrite(res, file.path(out, "MVMR_RESTRICTED_412_IVW.tsv"), sep = "\t")

z <- make_input("COPD")
invisible(capture.output(strength <- strength_mvmr(z$input, gencov = z$gencov)))
strength_out <- data.table(
  exposure = c("BMI", "Smoking initiation"),
  conditional_F = as.numeric(strength[1, c("exposure1", "exposure2")]),
  threshold_10 = ifelse(as.numeric(strength[1, c("exposure1", "exposure2")]) >= 10, "PASS", "FAIL"),
  n_instruments = nrow(z$input),
  selection = "423-set restricted to BMI or smoking-initiation P<5e-8"
)
fwrite(strength_out, file.path(out, "MVMR_RESTRICTED_412_CONDITIONAL_F.tsv"), sep = "\t")

full <- fread(file.path(root, "final_revision_closeout/10_submission_package/07_Additional_file_5_MVMR_Full_Results.tsv"))
full <- full[method == "MVMR-IVW", .(outcome, exposure, full_n = n_instruments, full_beta = beta, full_SE = SE, full_P = P, full_FDR_q = FDR_q)]
restricted_long <- rbindlist(list(
  res[, .(outcome, exposure = "BMI", restricted_n = n_instruments, restricted_beta = BMI_beta, restricted_SE = BMI_SE, restricted_P = BMI_P, restricted_FDR_q = BMI_FDR_q)],
  res[, .(outcome, exposure = "Smoking initiation", restricted_n = n_instruments, restricted_beta = SMK_beta, restricted_SE = SMK_SE, restricted_P = SMK_P, restricted_FDR_q = SMK_FDR_q)]
))
comparison <- merge(full, restricted_long, by = c("outcome", "exposure"), all = TRUE)
comparison[, `:=`(
  beta_difference = restricted_beta - full_beta,
  absolute_beta_difference = abs(restricted_beta - full_beta),
  direction_concordant = sign(restricted_beta) == sign(full_beta),
  nominal_status_concordant = (restricted_P < .05) == (full_P < .05),
  FDR_status_concordant = (restricted_FDR_q < .05) == (full_FDR_q < .05)
)]
fwrite(comparison, file.path(out, "MVMR_RESTRICTED_412_COMPARISON.tsv"), sep = "\t")

key_outcomes <- c("HF", "coronary_heart", "diabetes_type_2", "COPD", "OA")
key <- comparison[outcome %chin% key_outcomes]
report <- c(
  "# Restricted 412-variant MVMR sensitivity", "",
  "The final 423-variant matrix contains 11 variants that were not genome-wide significant for BMI or smoking initiation. A bounded sensitivity analysis excluded those 11 variants and repeated only conditional F statistics and covariance-aware MVMR-IVW with the same allele coding, SNP-specific covariance matrices, outcome files, multiplicity family and estimator.", "",
  sprintf("The restricted analysis used 412 variants for 11 outcomes and 411 variants for four outcomes with the existing rs705217 missingness pattern. Conditional F statistics were %.4f for BMI and %.4f for smoking initiation.", strength_out$conditional_F[1], strength_out$conditional_F[2]), "",
  sprintf("Direction was concordant with the full model for %d/30 exposure-outcome estimates. Nominal significance status was concordant for %d/30 and FDR significance status for %d/30.", sum(comparison$direction_concordant), sum(comparison$nominal_status_concordant), sum(comparison$FDR_status_concordant)), "",
  "## Priority outcomes", "",
  paste(capture.output(print(key)), collapse = "\n"), "",
  "This sensitivity analysis does not replace the pre-specified 423-variant model. It isolates the influence of the 11 variants that entered through another exposure axis and is reported as a provenance-focused robustness check."
)
writeLines(report, file.path(out, "MVMR_RESTRICTED_412_REPORT.md"))

cat("MVMR_RESTRICTED_412_COMPLETE\n")
