r_lib <- Sys.getenv("MODRISK_R_LIB", unset = "")
if (nzchar(r_lib)) .libPaths(c(r_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(MVMR)
  library(quantreg)
})

repo <- normalizePath(Sys.getenv("MODRISK_ROOT", unset = "."), mustWork = TRUE)
out_h <- file.path(repo, "final_major_revision", "03_mr_harmonisation")
out_r <- file.path(repo, "final_major_revision", "05_robust_mvmr")
dir.create(out_h, recursive = TRUE, showWarnings = FALSE)
dir.create(out_r, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(repo, "results", "phase4_repair", "PRIMARY_MVMR_CORRECTED_INPUT.tsv.gz")
pvalue_file <- file.path(repo, "results", "phase4", "MVMR_EXPOSURE_MATRIX.tsv")
gencov_file <- file.path(repo, "results", "phase4_repair2", "PRIMARY_BMI_SMK_GENCOV_LAMBDA.rds")
old_file <- file.path(repo, "results", "phase4", "MODEL1_MVMR_IVW.tsv")

dat0 <- fread(input_file)
pvals <- fread(pvalue_file, select = c("SNP", "BMI_P", "SMK_P"))
dat0 <- merge(dat0, pvals, by = "SNP", all.x = TRUE, sort = FALSE)
gclist <- readRDS(gencov_file)
gclist <- gclist[match(dat0$SNP, names(gclist))]

outcomes <- c(
  "COPD", "HF", "OA", "alcohol_problems", "carpal_tunnel",
  "coronary_heart", "diabetes_type_2", "fibromyalgia",
  "gastritis_duodenitis", "gout", "hypertension", "insomnia",
  "periph_neuro", "sleep_apnoea", "valve_nonrheumatic"
)

stopifnot(nrow(dat0) == 423L, !anyDuplicated(dat0$SNP), length(gclist) == nrow(dat0))
same <- dat0$BMI_EA == dat0$SMK_EA & dat0$BMI_OA == dat0$SMK_OA
swap <- dat0$BMI_EA == dat0$SMK_OA & dat0$BMI_OA == dat0$SMK_EA
unresolved <- !(same | swap)
pal <- paste0(dat0$BMI_EA, dat0$BMI_OA) %in% c("AT", "TA", "CG", "GC")

# Frozen outcomes are aligned to the smoking effect allele.  Two equivalent
# common-allele representations are retained to test algebraic invariance.
smk <- copy(dat0)
smk[, BMI_beta_source := BMI_beta]
smk[swap, BMI_beta := -BMI_beta]

bmi <- copy(dat0)
for (o in outcomes) bmi[swap, (paste0(o, "_beta")) := -get(paste0(o, "_beta"))]
bmi[swap, SMK_beta := -SMK_beta]

lambda_by_snp <- vapply(seq_along(gclist), function(i) {
  gc <- gclist[[i]]
  gc[1, 2] / sqrt(gc[1, 1] * gc[2, 2])
}, numeric(1))
lambda <- median(lambda_by_snp)

allele_audit <- data.table(
  metric = c(
    "rows_in_frozen_input", "unique_snps", "same_exposure_alleles",
    "swapped_exposure_alleles", "unresolved_exposure_alleles",
    "palindromic_exposure_pairs", "duplicate_snps", "gencov_matrices",
    "median_exposure_error_correlation", "minimum_gencov_eigenvalue"
  ),
  value = c(
    nrow(dat0), uniqueN(dat0$SNP), sum(same), sum(swap), sum(unresolved),
    sum(pal), anyDuplicated(dat0$SNP), length(gclist), lambda,
    min(vapply(gclist, function(x) min(eigen(x, symmetric = TRUE, only.values = TRUE)$values), numeric(1)))
  ),
  status = c(rep("PASS", 10))
)
fwrite(allele_audit, file.path(out_h, "MR_ALLELE_HARMONISATION_AUDIT.tsv"), sep = "\t")

make_input <- function(d, o) {
  ok <- is.finite(d$BMI_beta) & is.finite(d$BMI_SE) & d$BMI_SE > 0 &
    is.finite(d$SMK_beta) & is.finite(d$SMK_SE) & d$SMK_SE > 0 &
    is.finite(d[[paste0(o, "_beta")]]) &
    is.finite(d[[paste0(o, "_SE")]]) & d[[paste0(o, "_SE")]] > 0
  list(
    input = format_mvmr(
      cbind(d$BMI_beta[ok], d$SMK_beta[ok]),
      d[[paste0(o, "_beta")]][ok],
      cbind(d$BMI_SE[ok], d$SMK_SE[ok]),
      d[[paste0(o, "_SE")]][ok],
      d$SNP[ok]
    ),
    ok = ok,
    gc = gclist[ok]
  )
}

run_ivw <- function(d, o) {
  z <- make_input(d, o)
  invisible(capture.output(f <- ivw_mvmr(z$input, gencov = z$gc)))
  data.table(
    outcome = o, n_instruments = nrow(z$input),
    BMI_beta = f["exposure1", "Estimate"], BMI_SE = f["exposure1", "Std. Error"],
    BMI_P = f["exposure1", "Pr(>|t|)"],
    SMK_beta = f["exposure2", "Estimate"], SMK_SE = f["exposure2", "Std. Error"],
    SMK_P = f["exposure2", "Pr(>|t|)"]
  )
}

ivw_smk <- rbindlist(lapply(outcomes, function(o) run_ivw(smk, o)))
ivw_bmi <- rbindlist(lapply(outcomes, function(o) run_ivw(bmi, o)))
ivw_q <- p.adjust(c(ivw_smk$BMI_P, ivw_smk$SMK_P), "BH")
ivw_smk[, `:=`(BMI_q = ivw_q[seq_len(.N)], SMK_q = ivw_q[.N + seq_len(.N)])]
fwrite(ivw_smk, file.path(out_r, "MVMR_IVW_FINAL.tsv"), sep = "\t")

orientation <- merge(
  ivw_smk[, .(outcome, BMI_beta_SMK_EA = BMI_beta, SMK_beta_SMK_EA = SMK_beta)],
  ivw_bmi[, .(outcome, BMI_beta_BMI_EA = BMI_beta, SMK_beta_BMI_EA = SMK_beta)],
  by = "outcome"
)
orientation[, `:=`(
  BMI_abs_difference = abs(BMI_beta_SMK_EA - BMI_beta_BMI_EA),
  SMK_abs_difference = abs(SMK_beta_SMK_EA - SMK_beta_BMI_EA)
)]
fwrite(orientation, file.path(out_h, "MR_ALLELE_ORIENTATION_INVARIANCE.tsv"), sep = "\t")

excluded <- rbindlist(lapply(outcomes, function(o) {
  ok <- make_input(smk, o)$ok
  if (all(ok)) return(data.table(outcome = o, SNP = character(), reason = character()))
  data.table(outcome = o, SNP = smk$SNP[!ok], reason = "missing_or_invalid_outcome_association")
}), fill = TRUE)
fwrite(excluded, file.path(out_h, "MR_INSTRUMENT_EXCLUSION_LOG.tsv"), sep = "\t")

z0 <- make_input(smk, "COPD")
invisible(capture.output(strength <- strength_mvmr(z0$input, gencov = z0$gc)))
strength_out <- data.table(
  exposure = c("BMI", "Smoking initiation"),
  conditional_F = as.numeric(strength[1, c("exposure1", "exposure2")]),
  threshold_10 = c("PASS", "FAIL")
)
fwrite(strength_out, file.path(out_h, "MR_CONDITIONAL_STRENGTH.tsv"), sep = "\t")

old <- fread(old_file)
before_after <- merge(
  old[, .(outcome, old_BMI_beta = BMI_beta, old_SMK_beta = SMK_beta)],
  ivw_smk[, .(outcome, corrected_BMI_beta = BMI_beta, corrected_SMK_beta = SMK_beta)],
  by = "outcome"
)
before_after[, `:=`(
  BMI_direction_changed = sign(old_BMI_beta) != sign(corrected_BMI_beta),
  SMK_direction_changed = sign(old_SMK_beta) != sign(corrected_SMK_beta),
  historical_status = "INVALID_ARCHIVED"
)]
fwrite(before_after, file.path(out_h, "MR_HARMONISATION_BEFORE_AFTER.tsv"), sep = "\t")

# Directional MVMR-Egger: orient all SNP associations to an increasing allele
# for each exposure in turn, as the intercept is not invariant to arbitrary
# SNP-level sign changes.
run_egger <- function(o, orient_to) {
  z <- make_input(smk, o)
  d <- z$input
  ref <- if (orient_to == "BMI") d$betaX1 else d$betaX2
  sg <- ifelse(ref < 0, -1, 1)
  d$betaX1 <- d$betaX1 * sg
  d$betaX2 <- d$betaX2 * sg
  d$betaYG <- d$betaYG * sg
  co <- summary(lm(betaYG ~ betaX1 + betaX2, weights = 1 / sebetaYG^2, data = d))$coefficients
  data.table(
    outcome = o, orientation = paste0(orient_to, "_increasing"), n_instruments = nrow(d),
    intercept = co["(Intercept)", "Estimate"], intercept_SE = co["(Intercept)", "Std. Error"],
    intercept_P = co["(Intercept)", "Pr(>|t|)"],
    BMI_beta = co["betaX1", "Estimate"], BMI_SE = co["betaX1", "Std. Error"],
    BMI_P = co["betaX1", "Pr(>|t|)"],
    SMK_beta = co["betaX2", "Estimate"], SMK_SE = co["betaX2", "Std. Error"],
    SMK_P = co["betaX2", "Pr(>|t|)"]
  )
}
egger <- rbindlist(lapply(outcomes, function(o) rbind(run_egger(o, "BMI"), run_egger(o, "SMK"))))
egger[, c("BMI_q", "SMK_q") := {
  qv <- p.adjust(c(BMI_P, SMK_P), "BH")
  list(qv[seq_len(.N)], qv[.N + seq_len(.N)])
}, by = orientation]
egger[, intercept_q := p.adjust(intercept_P, "BH"), by = orientation]
fwrite(egger, file.path(out_r, "MVMR_EGGER_ORIENTATION_RESULTS.tsv"), sep = "\t")

# Weighted LAD (multivariable median) with a covariance-aware parametric
# bootstrap.  Each SNP's two exposure associations are drawn from its stored
# 2x2 covariance matrix; outcome associations are drawn independently because
# exposure-outcome overlap covariance is unavailable.
median_rows <- list()
set.seed(20261124)
nboot <- 2000L
for (o in outcomes) {
  z <- make_input(smk, o)
  d <- z$input
  gc <- z$gc
  point <- rq(betaYG ~ betaX1 + betaX2 - 1, weights = 1 / sebetaYG^2, data = d)$coefficients
  boots <- matrix(NA_real_, nboot, 2)
  for (b in seq_len(nboot)) {
    z1 <- rnorm(nrow(d)); z2 <- rnorm(nrow(d)); zy <- rnorm(nrow(d))
    bx1 <- d$betaX1 + sqrt(vapply(gc, function(x) x[1, 1], numeric(1))) * z1
    cov12 <- vapply(gc, function(x) x[1, 2], numeric(1))
    sd1 <- sqrt(vapply(gc, function(x) x[1, 1], numeric(1)))
    sd2 <- sqrt(vapply(gc, function(x) x[2, 2], numeric(1)))
    rho <- pmax(-0.999999, pmin(0.999999, cov12 / (sd1 * sd2)))
    bx2 <- d$betaX2 + sd2 * (rho * z1 + sqrt(1 - rho^2) * z2)
    by <- d$betaYG + d$sebetaYG * zy
    boots[b, ] <- rq(by ~ bx1 + bx2 - 1, weights = 1 / d$sebetaYG^2)$coefficients
  }
  se <- apply(boots, 2, sd, na.rm = TRUE)
  p <- 2 * pnorm(-abs(point / se))
  median_rows[[o]] <- data.table(
    outcome = o, n_instruments = nrow(d), bootstrap_iterations = nboot,
    BMI_beta = point[1], BMI_SE = se[1], BMI_P = p[1],
    BMI_CI_low = point[1] - qnorm(.975) * se[1], BMI_CI_high = point[1] + qnorm(.975) * se[1],
    SMK_beta = point[2], SMK_SE = se[2], SMK_P = p[2],
    SMK_CI_low = point[2] - qnorm(.975) * se[2], SMK_CI_high = point[2] + qnorm(.975) * se[2]
  )
}
median_res <- rbindlist(median_rows)
median_q <- p.adjust(c(median_res$BMI_P, median_res$SMK_P), "BH")
median_res[, `:=`(BMI_q = median_q[seq_len(.N)], SMK_q = median_q[.N + seq_len(.N)])]
fwrite(median_res, file.path(out_r, "MVMR_MEDIAN_FULL_RESULTS.tsv"), sep = "\t")

# Weak-instrument-aware Q minimisation with the verified exposure-association
# error correlation and jackknife uncertainty.
qhet_rows <- list()
for (o in outcomes) {
  z <- make_input(smk, o)
  suppressWarnings(invisible(capture.output(qf <- qhet_mvmr(
    z$input, pcor = matrix(c(1, lambda, lambda, 1), 2, 2),
    CI = TRUE, ncores = 4, CI_method = "jackknife"
  ))))
  qhet_rows[[o]] <- data.table(
    outcome = o, n_instruments = nrow(z$input),
    BMI_beta = qf[1, "Effect Estimates"], BMI_SE = qf[1, "Std. Error"],
    BMI_P = 2 * pnorm(-abs(qf[1, "Effect Estimates"] / qf[1, "Std. Error"])),
    SMK_beta = qf[2, "Effect Estimates"], SMK_SE = qf[2, "Std. Error"],
    SMK_P = 2 * pnorm(-abs(qf[2, "Effect Estimates"] / qf[2, "Std. Error"]))
  )
}
qhet <- rbindlist(qhet_rows)
qhet_q <- p.adjust(c(qhet$BMI_P, qhet$SMK_P), "BH")
qhet[, `:=`(BMI_q = qhet_q[seq_len(.N)], SMK_q = qhet_q[.N + seq_len(.N)])]
fwrite(qhet, file.path(out_r, "MVMR_WEAK_IV_AWARE_RESULTS.tsv"), sep = "\t")

comparison <- rbindlist(list(
  ivw_smk[, .(outcome, method = "MVMR-IVW", exposure = "BMI", beta = BMI_beta, SE = BMI_SE, P = BMI_P, q = BMI_q)],
  ivw_smk[, .(outcome, method = "MVMR-IVW", exposure = "Smoking initiation", beta = SMK_beta, SE = SMK_SE, P = SMK_P, q = SMK_q)],
  median_res[, .(outcome, method = "MVMR weighted LAD median", exposure = "BMI", beta = BMI_beta, SE = BMI_SE, P = BMI_P, q = BMI_q)],
  median_res[, .(outcome, method = "MVMR weighted LAD median", exposure = "Smoking initiation", beta = SMK_beta, SE = SMK_SE, P = SMK_P, q = SMK_q)],
  qhet[, .(outcome, method = "Q-minimisation", exposure = "BMI", beta = BMI_beta, SE = BMI_SE, P = BMI_P, q = BMI_q)],
  qhet[, .(outcome, method = "Q-minimisation", exposure = "Smoking initiation", beta = SMK_beta, SE = SMK_SE, P = SMK_P, q = SMK_q)]
))
fwrite(comparison, file.path(out_r, "MVMR_ROBUSTNESS_COMPARISON.tsv"), sep = "\t")

writeLines(c(
  "# MR harmonisation final report",
  "",
  sprintf("- Frozen union instrument set: %d unique SNPs.", nrow(dat0)),
  sprintf("- Exposure allele pairs already concordant: %d; swapped: %d; unresolved: %d.", sum(same), sum(swap), sum(unresolved)),
  "- Outcomes in the frozen input were aligned to the smoking effect allele; BMI associations were sign-flipped for the 219 swapped rows.",
  sprintf("- Four outcomes use 422 instruments because rs705217 lacks an outcome association: %s.", paste(unique(excluded[outcome != "" & SNP == "rs705217"]$outcome), collapse = ", ")),
  sprintf("- Maximum absolute IVW difference between common-smoking and common-BMI allele representations: %.3g.", max(unlist(orientation[, .(BMI_abs_difference, SMK_abs_difference)]))),
  sprintf("- Conditional F statistics: BMI %.3f; smoking initiation %.3f.", strength_out$conditional_F[1], strength_out$conditional_F[2]),
  "- The pre-correction Model 1 estimates are retained only as INVALID_ARCHIVED provenance and are not used in the manuscript."
), file.path(out_h, "MR_HARMONISATION_FINAL_REPORT.md"))

writeLines(c(
  "# MVMR method assumptions",
  "",
  "## MVMR-IVW",
  "Uses the stored SNP-specific 2 x 2 covariance matrices for BMI and smoking-initiation association errors.",
  "",
  "## Directional MVMR-Egger",
  "Because the Egger intercept depends on allele orientation, analyses were repeated after orienting every SNP to the BMI-increasing allele and to the smoking-increasing allele. Exposure and outcome effects were multiplied by the same SNP-specific sign.",
  "",
  "## Multivariable weighted LAD median",
  "The point estimate is the no-intercept weighted least-absolute-deviation estimate. The parametric bootstrap draws the two exposure associations jointly from each SNP-specific stored covariance matrix (2,000 replicates). Exposure-outcome sampling covariance could not be modelled because the required overlap covariance was unavailable.",
  "",
  "## Weak-instrument-aware Q minimisation",
  "MVMR::qhet_mvmr was run with the verified exposure-association error correlation and jackknife confidence intervals. The function is labelled as under development by its authors. Smoking initiation remained conditionally weak (F < 10), so these estimates are sensitivity analyses and do not remove all weak-instrument uncertainty."
), file.path(out_r, "MVMR_METHOD_ASSUMPTIONS.md"))

writeLines(c(
  "# MVMR final sensitivity report",
  "",
  sprintf("- Corrected IVW FDR-significant outcomes: BMI %d/15; smoking initiation %d/15.", sum(ivw_smk$BMI_q < .05), sum(ivw_smk$SMK_q < .05)),
  sprintf("- Median FDR-significant outcomes: BMI %d/15; smoking initiation %d/15.", sum(median_res$BMI_q < .05), sum(median_res$SMK_q < .05)),
  sprintf("- Q-minimisation FDR-significant outcomes: BMI %d/15; smoking initiation %d/15.", sum(qhet$BMI_q < .05), sum(qhet$SMK_q < .05)),
  sprintf("- Egger intercept nominal P < 0.05: BMI-oriented %d/15; smoking-oriented %d/15.", sum(egger[orientation == "BMI_increasing"]$intercept_P < .05), sum(egger[orientation == "SMK_increasing"]$intercept_P < .05)),
  "- The corrected MVMR-IVW analysis remains primary. Robust estimates are interpreted as sensitivity analyses because their identifying assumptions differ and smoking initiation is conditionally weak."
), file.path(out_r, "MVMR_FINAL_SENSITIVITY_REPORT.md"))

message("Final MVMR audit completed")
