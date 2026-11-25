r_lib <- Sys.getenv("MODRISK_R_LIB", unset = "")
if (nzchar(r_lib)) .libPaths(c(r_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table); library(MVMR)})

root <- normalizePath(Sys.getenv("MODRISK_ROOT", unset = "."), mustWork = TRUE)
out <- file.path(root, "final_revision_closeout", "02_mr_diagnostics")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

d0 <- fread(file.path(root, "results/phase4_repair/PRIMARY_MVMR_CORRECTED_INPUT.tsv.gz"))
g <- readRDS(file.path(root, "results/phase4_repair2/PRIMARY_BMI_SMK_GENCOV_LAMBDA.rds"))
g <- g[match(d0$SNP, names(g))]
outcomes <- c("COPD", "HF", "OA", "alcohol_problems", "carpal_tunnel",
              "coronary_heart", "diabetes_type_2", "fibromyalgia",
              "gastritis_duodenitis", "gout", "hypertension", "insomnia",
              "periph_neuro", "sleep_apnoea", "valve_nonrheumatic")
swap <- d0$BMI_EA == d0$SMK_OA & d0$BMI_OA == d0$SMK_EA
same <- d0$BMI_EA == d0$SMK_EA & d0$BMI_OA == d0$SMK_OA
stopifnot(sum(swap) == 219L, sum(same) == 204L, all(swap | same))

# Common smoking-effect-allele representation. For swapped rows BMI, smoking
# and outcome effects are all on the same smoking allele; relative to common BMI
# coding this is a joint -1 recoding, so the exposure covariance matrix is unchanged.
d <- copy(d0)
d[swap, BMI_beta := -BMI_beta]

make_input <- function(o) {
  ok <- is.finite(d$BMI_beta) & is.finite(d$BMI_SE) & d$BMI_SE > 0 &
        is.finite(d$SMK_beta) & is.finite(d$SMK_SE) & d$SMK_SE > 0 &
        is.finite(d[[paste0(o, "_beta")]]) &
        is.finite(d[[paste0(o, "_SE")]]) & d[[paste0(o, "_SE")]] > 0
  list(
    input = format_mvmr(cbind(d$BMI_beta[ok], d$SMK_beta[ok]),
                        d[[paste0(o, "_beta")]][ok],
                        cbind(d$BMI_SE[ok], d$SMK_SE[ok]),
                        d[[paste0(o, "_SE")]][ok], d$SNP[ok]),
    gc = g[ok]
  )
}

het_rows <- list(); loo_rows <- list(); loo_full <- list(); conv <- list()
for (o in outcomes) {
  z <- make_input(o); x <- z$input
  fit_ok <- TRUE; fit_msg <- ""
  fit <- tryCatch({
    invisible(capture.output(v <- ivw_mvmr(x, gencov = z$gc))); v
  }, error = function(e) { fit_ok <<- FALSE; fit_msg <<- conditionMessage(e); NULL })
  if (!fit_ok) stop(sprintf("IVW failed for %s: %s", o, fit_msg))
  q_ok <- TRUE; q_msg <- ""
  q <- tryCatch({
    invisible(capture.output(v <- pleiotropy_mvmr(x, gencov = z$gc, estimator = "ivw"))); v
  }, error = function(e) { q_ok <<- FALSE; q_msg <<- conditionMessage(e); NULL })
  if (!q_ok) stop(sprintf("Q failed for %s: %s", o, q_msg))

  full_bmi <- as.numeric(fit["exposure1", "Estimate"])
  full_smk <- as.numeric(fit["exposure2", "Estimate"])
  loo <- rbindlist(lapply(seq_len(nrow(x)), function(j) {
    ok <- TRUE; msg <- ""
    lf <- tryCatch({
      invisible(capture.output(v <- ivw_mvmr(x[-j, ], gencov = z$gc[-j]))); v
    }, error = function(e) { ok <<- FALSE; msg <<- conditionMessage(e); NULL })
    data.table(
      outcome = o, omitted_SNP = x$SNP[j], success = ok,
      BMI_beta = if (ok) as.numeric(lf["exposure1", "Estimate"]) else NA_real_,
      SMK_beta = if (ok) as.numeric(lf["exposure2", "Estimate"]) else NA_real_,
      message = msg
    )
  }))
  loo_full[[o]] <- loo
  good <- loo[success == TRUE]
  bmi_rev <- good[sign(BMI_beta) != sign(full_bmi)]
  smk_rev <- good[sign(SMK_beta) != sign(full_smk)]
  dbmi <- abs(good$BMI_beta - full_bmi); dsmk <- abs(good$SMK_beta - full_smk)
  infl <- unique(c(good$omitted_SNP[which(dbmi == max(dbmi))],
                   good$omitted_SNP[which(dsmk == max(dsmk))],
                   bmi_rev$omitted_SNP, smk_rev$omitted_SNP))
  loo_rows[[o]] <- data.table(
    outcome = o, n_instruments = nrow(x),
    BMI_full_beta = full_bmi, BMI_LOO_min_beta = min(good$BMI_beta),
    BMI_LOO_max_beta = max(good$BMI_beta),
    BMI_LOO_sign_reversal_count = nrow(bmi_rev),
    Smoking_full_beta = full_smk, Smoking_LOO_min_beta = min(good$SMK_beta),
    Smoking_LOO_max_beta = max(good$SMK_beta),
    Smoking_LOO_sign_reversal_count = nrow(smk_rev),
    influential_variant_ids = paste(infl, collapse = ";"),
    convergence_status = if (nrow(good) == nrow(x)) "SUCCESS" else "PARTIAL"
  )
  het_rows[[o]] <- data.table(
    outcome = o, n_instruments = nrow(x),
    heterogeneity_method = "MVMR::pleiotropy_mvmr IVW Q with SNP-specific exposure-association covariance",
    Q = q$Qstat, Q_df = nrow(x) - 2L, Q_P = q$Qpval
  )
  conv[[o]] <- data.table(
    estimator = "MVMR-IVW/Q/LOO", outcome = o, planned_units = nrow(x) + 2L,
    successful_units = nrow(good) + as.integer(fit_ok) + as.integer(q_ok),
    failed_units = nrow(x) - nrow(good), convergence_status = if (nrow(good) == nrow(x)) "SUCCESS" else "PARTIAL",
    message = paste(fit_msg, q_msg)
  )
}

het <- rbindlist(het_rows)
het[, Q_FDR := p.adjust(Q_P, method = "BH")]
het[, `:=`(nominal_significant = Q_P < 0.05, FDR_significant = Q_FDR < 0.05)]
loo <- rbindlist(loo_rows)
fwrite(het, file.path(out, "MVMR_HETEROGENEITY_FINAL.tsv"), sep = "\t")
fwrite(loo, file.path(out, "MVMR_LEAVE_ONE_OUT_FINAL.tsv"), sep = "\t")
fwrite(rbindlist(loo_full), file.path(out, "MVMR_LEAVE_ONE_OUT_PER_SNP.tsv"), sep = "\t")
fwrite(rbindlist(conv), file.path(out, "MR_CONVERGENCE_DIAGNOSTICS_IVW_Q_LOO.tsv"), sep = "\t")

summary <- c(
  "# Final MVMR heterogeneity and leave-one-SNP-out diagnostics",
  "",
  sprintf("- Q heterogeneity was nominally significant for %d/15 outcomes and BH-FDR significant for %d/15 outcomes.", sum(het$nominal_significant), sum(het$FDR_significant)),
  sprintf("- BMI leave-one-SNP-out sign reversals: %d omitted-SNP events across %d outcomes.", sum(loo$BMI_LOO_sign_reversal_count), sum(loo$BMI_LOO_sign_reversal_count > 0)),
  sprintf("- Smoking-initiation leave-one-SNP-out sign reversals: %d omitted-SNP events across %d outcomes.", sum(loo$Smoking_LOO_sign_reversal_count), sum(loo$Smoking_LOO_sign_reversal_count > 0)),
  "- Q was calculated with MVMR 0.4.9 pleiotropy_mvmr using the final 423/422-SNP harmonised inputs and SNP-specific exposure-association covariance matrices. Degrees of freedom are n instruments minus two exposure coefficients.",
  "- Q heterogeneity is an indirect diagnostic of model misfit or variant-level heterogeneity and does not by itself prove horizontal pleiotropy.",
  "- Leave-one-out sign changes are influence diagnostics, not reverse-direction or bidirectional MR tests. No reverse-direction or bidirectional MR was performed."
)
writeLines(summary, file.path(out, "MVMR_DIAGNOSTICS_SUMMARY.md"))
cat("MR diagnostics complete\n")
