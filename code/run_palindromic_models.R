#!/usr/bin/env Rscript
r_lib <- Sys.getenv("MODRISK_R_LIB", unset = "")
if (nzchar(r_lib)) .libPaths(c(r_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(MVMR)
})

root <- normalizePath(Sys.getenv("MODRISK_ROOT", unset = "."), mustWork = TRUE)
work <- file.path(root, "final_manuscript_error_correction")
out_pal <- file.path(work, "01_palindromic_validation")
out_mr <- file.path(work, "02_mr_sensitivity")
dir.create(out_pal, recursive = TRUE, showWarnings = FALSE)
dir.create(out_mr, recursive = TRUE, showWarnings = FALSE)

dat0 <- fread(file.path(root, "results/phase4_repair/PRIMARY_MVMR_CORRECTED_INPUT.tsv.gz"))
pvals <- fread(file.path(root, "results/phase4/MVMR_EXPOSURE_MATRIX.tsv"),
               select = c("SNP", "BMI_P", "SMK_P"))
ord <- dat0$SNP
dat0 <- merge(dat0, pvals, by = "SNP", all.x = TRUE, sort = FALSE)
dat0 <- dat0[match(ord, SNP)]
gencov0 <- readRDS(file.path(root, "results/phase4_repair2/PRIMARY_BMI_SMK_GENCOV_LAMBDA.rds"))
gencov0 <- gencov0[match(dat0$SNP, names(gencov0))]

audit <- fread(file.path(root, "pre_submission_technical_closeout/03_palindromic_audit/PALINDROMIC_SNP_EAF_AUDIT.tsv"))
unresolved <- audit[decision == "UNRESOLVED_FREQUENCY_AMBIGUOUS", SNP]
resolved <- audit[decision == "RESOLVED_BY_EAF", SNP]
stopifnot(nrow(dat0) == 423L, length(unresolved) == 17L,
          length(resolved) == 37L, !anyDuplicated(dat0$SNP))

same <- dat0$BMI_EA == dat0$SMK_EA & dat0$BMI_OA == dat0$SMK_OA
swap <- dat0$BMI_EA == dat0$SMK_OA & dat0$BMI_OA == dat0$SMK_EA
stopifnot(all(same | swap))
dat0[swap, BMI_beta := -BMI_beta]

model_b <- (is.finite(dat0$BMI_P) & dat0$BMI_P < 5e-8) |
  (is.finite(dat0$SMK_P) & dat0$SMK_P < 5e-8)
model_c <- !(dat0$SNP %chin% unresolved)
sets <- list(
  A_FULL_423 = rep(TRUE, nrow(dat0)),
  B_EXPOSURE_GWS_412 = model_b,
  C_STRICT_406 = model_c,
  D_EXPOSURE_GWS_STRICT_395 = model_b & model_c
)
expected <- c(423L, 412L, 406L, 395L)
stopifnot(identical(as.integer(vapply(sets, sum, integer(1))), expected))

outcomes <- c(
  "COPD", "HF", "OA", "alcohol_problems", "carpal_tunnel", "coronary_heart",
  "diabetes_type_2", "fibromyalgia", "gastritis_duodenitis", "gout",
  "hypertension", "insomnia", "periph_neuro", "sleep_apnoea", "valve_nonrheumatic"
)

make_input <- function(mask, outcome) {
  d <- dat0[mask]
  gc <- gencov0[mask]
  by <- d[[paste0(outcome, "_beta")]]
  sy <- d[[paste0(outcome, "_SE")]]
  ok <- is.finite(d$BMI_beta) & is.finite(d$BMI_SE) & d$BMI_SE > 0 &
    is.finite(d$SMK_beta) & is.finite(d$SMK_SE) & d$SMK_SE > 0 &
    is.finite(by) & is.finite(sy) & sy > 0
  list(
    input = format_mvmr(
      cbind(d$BMI_beta[ok], d$SMK_beta[ok]), by[ok],
      cbind(d$BMI_SE[ok], d$SMK_SE[ok]), sy[ok], d$SNP[ok]
    ),
    gc = gc[ok],
    included = d$SNP[ok]
  )
}

run_one <- function(model, mask, outcome) {
  z <- make_input(mask, outcome)
  invisible(capture.output(fit <- ivw_mvmr(z$input, gencov = z$gc)))
  invisible(capture.output(qfit <- pleiotropy_mvmr(z$input, gencov = z$gc, estimator = "ivw")))
  invisible(capture.output(strength <- strength_mvmr(z$input, gencov = z$gc)))
  n <- nrow(z$input)
  crit <- qt(0.975, df = n - 2L)
  rbind(
    data.table(model = model, outcome = outcome, exposure = "BMI",
               n_instruments = n, beta = fit["exposure1", "Estimate"],
               SE = fit["exposure1", "Std. Error"], P = fit["exposure1", "Pr(>|t|)"],
               CI_lower = fit["exposure1", "Estimate"] - crit * fit["exposure1", "Std. Error"],
               CI_upper = fit["exposure1", "Estimate"] + crit * fit["exposure1", "Std. Error"],
               conditional_F = as.numeric(strength[1, "exposure1"]),
               Q = qfit$Qstat, Q_df = n - 2L, Q_P = qfit$Qpval),
    data.table(model = model, outcome = outcome, exposure = "Smoking initiation",
               n_instruments = n, beta = fit["exposure2", "Estimate"],
               SE = fit["exposure2", "Std. Error"], P = fit["exposure2", "Pr(>|t|)"],
               CI_lower = fit["exposure2", "Estimate"] - crit * fit["exposure2", "Std. Error"],
               CI_upper = fit["exposure2", "Estimate"] + crit * fit["exposure2", "Std. Error"],
               conditional_F = as.numeric(strength[1, "exposure2"]),
               Q = qfit$Qstat, Q_df = n - 2L, Q_P = qfit$Qpval)
  )
}

rows <- list()
for (model in names(sets)) {
  for (o in outcomes) rows[[paste(model, o)]] <- run_one(model, sets[[model]], o)
}
res <- rbindlist(rows)
res[, FDR_q := p.adjust(P, method = "BH"), by = model]
qtab <- unique(res[, .(model, outcome, Q, Q_df, Q_P)])
qtab[, Q_FDR := p.adjust(Q_P, method = "BH"), by = model]
res <- merge(res, qtab[, .(model, outcome, Q_FDR)], by = c("model", "outcome"), all.x = TRUE)

base <- res[model == "A_FULL_423",
            .(outcome, exposure, A_beta = beta, A_P = P, A_FDR_q = FDR_q)]
res <- merge(res, base, by = c("outcome", "exposure"), all.x = TRUE)
res[, `:=`(
  beta_difference_from_A = beta - A_beta,
  absolute_beta_difference_from_A = abs(beta - A_beta),
  direction_concordant_with_A = sign(beta) == sign(A_beta),
  nominal_status_concordant_with_A = (P < 0.05) == (A_P < 0.05),
  FDR_status_concordant_with_A = (FDR_q < 0.05) == (A_FDR_q < 0.05)
)]
setorder(res, model, outcome, exposure)
fwrite(res, file.path(out_mr, "MVMR_PALINDROME_SENSITIVITY_RESULTS.tsv"), sep = "\t")

set_summary <- rbindlist(lapply(names(sets), function(nm) {
  mask <- sets[[nm]]
  data.table(
    model = nm,
    selected_instruments = sum(mask),
    unresolved_palindromic_retained = sum(dat0$SNP[mask] %chin% unresolved),
    exposure_GWS_retained = sum(model_b & mask),
    min_outcome_instruments = min(res[model == nm, n_instruments]),
    max_outcome_instruments = max(res[model == nm, n_instruments]),
    BMI_min_conditional_F = min(res[model == nm & exposure == "BMI", conditional_F]),
    Smoking_min_conditional_F = min(res[model == nm & exposure == "Smoking initiation", conditional_F]),
    directions_vs_A = sprintf("%d/30", sum(res[model == nm, direction_concordant_with_A])),
    nominal_status_vs_A = sprintf("%d/30", sum(res[model == nm, nominal_status_concordant_with_A])),
    FDR_status_vs_A = sprintf("%d/30", sum(res[model == nm, FDR_status_concordant_with_A]))
  )
}))
fwrite(set_summary, file.path(out_mr, "MVMR_INSTRUMENT_SET_COMPARISON.tsv"), sep = "\t")

manifest <- rbindlist(lapply(names(sets), function(nm) {
  data.table(
    SNP = dat0$SNP,
    model = nm,
    included = sets[[nm]],
    BMI_or_SMK_GWS = model_b,
    palindromic_audit_status = fifelse(dat0$SNP %chin% unresolved,
                                      "UNRESOLVED_FREQUENCY_AMBIGUOUS",
                                      fifelse(dat0$SNP %chin% resolved, "RESOLVED_BY_EAF", "NON_PALINDROMIC")),
    exclusion_reason = fifelse(!sets[[nm]] & dat0$SNP %chin% unresolved,
                              "frequency-ambiguous palindromic SNP excluded",
                              fifelse(!sets[[nm]] & !model_b,
                                      "not genome-wide significant for BMI or smoking initiation", ""))
  )
}))
fwrite(manifest, file.path(out_pal, "MVMR_STRICT_INSTRUMENT_MANIFEST.tsv"), sep = "\t")

audit_out <- copy(audit)
audit_out[, final_category := fifelse(
  decision == "RESOLVED_BY_EAF", "resolved by effect-allele frequency",
  fifelse(decision == "UNRESOLVED_FREQUENCY_AMBIGUOUS", "excluded as frequency-ambiguous", "requires review")
)]
audit_out[, strict_primary_included := !(SNP %chin% unresolved)]
fwrite(audit_out, file.path(out_pal, "PALINDROMIC_54_VARIANT_AUDIT.tsv"), sep = "\t")

disp <- audit_out[decision == "UNRESOLVED_FREQUENCY_AMBIGUOUS",
                  .(SNP, BMI_aligned_EAF, SMK_aligned_EAF, BMI_MAF, SMK_MAF,
                    current_orientation_EAF_difference, opposite_orientation_EAF_difference,
                    disposition = "EXCLUDED_FROM_STRICT_PRIMARY_MVMR",
                    reason = decision_reason)]
disp[, BMI_or_SMK_GWS := model_b[match(SNP, dat0$SNP)]]
fwrite(disp, file.path(out_pal, "UNRESOLVED_VARIANT_DISPOSITION.tsv"), sep = "\t")

key <- c("HF", "coronary_heart", "diabetes_type_2", "hypertension", "sleep_apnoea")
keytab <- res[outcome %chin% key,
              .(model, outcome, exposure, n_instruments, beta, SE, P, FDR_q,
                conditional_F, beta_difference_from_A, direction_concordant_with_A,
                FDR_status_concordant_with_A)]
fwrite(keytab, file.path(out_mr, "MVMR_KEY_OUTCOME_COMPARISON.tsv"), sep = "\t")

strict <- set_summary[model == "C_STRICT_406"]
report <- c(
  "# Final palindromic-variant decision",
  "",
  "The decision rule was frozen before model inspection. Thirty-seven palindromic SNPs were resolved by aligned effect-allele frequencies. Seventeen SNPs had intermediate allele frequencies that did not support reliable strand orientation and were excluded without imputation.",
  "",
  sprintf("Model A retained 423 variants; Model B retained 412; Model C retained 406; Model D retained 395. All 17 unresolved variants were members of Model B."),
  "",
  sprintf("For strict Model C, direction agreed with Model A for %s estimates, nominal-significance status for %s, and FDR status for %s.", strict$directions_vs_A, strict$nominal_status_vs_A, strict$FDR_status_vs_A),
  "",
  sprintf("The minimum outcome-specific conditional F statistics in Model C were %.3f for BMI and %.3f for smoking initiation.", strict$BMI_min_conditional_F, strict$Smoking_min_conditional_F),
  "",
  "Because unresolved palindromic variants cannot be validated from the available frequency evidence, Model C is designated the strict primary instrument set. Model A is retained only as a prespecified sensitivity analysis."
)
writeLines(report, file.path(out_pal, "PALINDROMIC_FINAL_DECISION.md"))

cat("PALINDROMIC_MODELS_COMPLETE\n")
