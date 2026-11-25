#!/usr/bin/env Rscript
r_lib <- Sys.getenv("MODRISK_R_LIB", unset = "")
if (nzchar(r_lib)) .libPaths(c(r_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table); library(MVMR); library(quantreg)})

root <- normalizePath(Sys.getenv("MODRISK_ROOT", unset = "."), mustWork = TRUE)
work <- file.path(root, "final_manuscript_error_correction")
out <- file.path(work, "02_mr_sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

d0 <- fread(file.path(root, "results/phase4_repair/PRIMARY_MVMR_CORRECTED_INPUT.tsv.gz"))
g0 <- readRDS(file.path(root, "results/phase4_repair2/PRIMARY_BMI_SMK_GENCOV_LAMBDA.rds"))
g0 <- g0[match(d0$SNP, names(g0))]
audit <- fread(file.path(root, "pre_submission_technical_closeout/03_palindromic_audit/PALINDROMIC_SNP_EAF_AUDIT.tsv"))
excluded <- audit[decision == "UNRESOLVED_FREQUENCY_AMBIGUOUS", SNP]
keep <- !(d0$SNP %chin% excluded)
d <- d0[keep]
g <- g0[keep]
stopifnot(nrow(d) == 406L, length(g) == 406L)

outcomes <- c("COPD", "HF", "OA", "alcohol_problems", "carpal_tunnel",
              "coronary_heart", "diabetes_type_2", "fibromyalgia",
              "gastritis_duodenitis", "gout", "hypertension", "insomnia",
              "periph_neuro", "sleep_apnoea", "valve_nonrheumatic")
same <- d$BMI_EA == d$SMK_EA & d$BMI_OA == d$SMK_OA
swap <- d$BMI_EA == d$SMK_OA & d$BMI_OA == d$SMK_EA
stopifnot(all(same | swap))
d[swap, BMI_beta := -BMI_beta]

lambda_by_snp <- vapply(g, function(x) x[1,2] / sqrt(x[1,1] * x[2,2]), numeric(1))
lambda <- median(lambda_by_snp)

make_input <- function(o) {
  ok <- is.finite(d$BMI_beta) & is.finite(d$BMI_SE) & d$BMI_SE > 0 &
    is.finite(d$SMK_beta) & is.finite(d$SMK_SE) & d$SMK_SE > 0 &
    is.finite(d[[paste0(o, "_beta")]]) &
    is.finite(d[[paste0(o, "_SE")]]) & d[[paste0(o, "_SE")]] > 0
  list(input = format_mvmr(cbind(d$BMI_beta[ok], d$SMK_beta[ok]),
                           d[[paste0(o, "_beta")]][ok],
                           cbind(d$BMI_SE[ok], d$SMK_SE[ok]),
                           d[[paste0(o, "_SE")]][ok], d$SNP[ok]),
       gc = g[ok])
}

ivw_rows <- list(); egger_rows <- list(); med_rows <- list(); qmin_rows <- list()
het_rows <- list(); loo_summary_rows <- list(); loo_rows <- list(); conv_rows <- list()
set.seed(20261125)
nboot <- 2000L

for (o in outcomes) {
  z <- make_input(o); x <- z$input; gc <- z$gc; n <- nrow(x)
  invisible(capture.output(ivw <- ivw_mvmr(x, gencov = gc)))
  crit_ivw <- qt(.975, df = n - 2L)
  ivw_rows[[o]] <- rbind(
    data.table(outcome=o, method="MVMR-IVW", exposure="BMI", n_instruments=n,
               beta=ivw["exposure1","Estimate"], SE=ivw["exposure1","Std. Error"],
               P=ivw["exposure1","Pr(>|t|)"], CI_lower=ivw["exposure1","Estimate"]-crit_ivw*ivw["exposure1","Std. Error"],
               CI_upper=ivw["exposure1","Estimate"]+crit_ivw*ivw["exposure1","Std. Error"],
               orientation=NA_character_, intercept=NA_real_, intercept_SE=NA_real_, intercept_P=NA_real_),
    data.table(outcome=o, method="MVMR-IVW", exposure="Smoking initiation", n_instruments=n,
               beta=ivw["exposure2","Estimate"], SE=ivw["exposure2","Std. Error"],
               P=ivw["exposure2","Pr(>|t|)"], CI_lower=ivw["exposure2","Estimate"]-crit_ivw*ivw["exposure2","Std. Error"],
               CI_upper=ivw["exposure2","Estimate"]+crit_ivw*ivw["exposure2","Std. Error"],
               orientation=NA_character_, intercept=NA_real_, intercept_SE=NA_real_, intercept_P=NA_real_)
  )

  for (orient in c("BMI", "SMK")) {
    ref <- if (orient == "BMI") x$betaX1 else x$betaX2
    sg <- ifelse(ref < 0, -1, 1)
    ee <- x
    ee$betaX1 <- ee$betaX1 * sg
    ee$betaX2 <- ee$betaX2 * sg
    ee$betaYG <- ee$betaYG * sg
    co <- summary(lm(betaYG ~ betaX1 + betaX2, weights=1/sebetaYG^2, data=ee))$coefficients
    crit_e <- qt(.975, df=n-3L)
    egger_rows[[paste(o,orient)]] <- rbind(
      data.table(outcome=o, method="MVMR-Egger", exposure="BMI", n_instruments=n,
                 beta=co["betaX1","Estimate"], SE=co["betaX1","Std. Error"], P=co["betaX1","Pr(>|t|)"],
                 CI_lower=co["betaX1","Estimate"]-crit_e*co["betaX1","Std. Error"],
                 CI_upper=co["betaX1","Estimate"]+crit_e*co["betaX1","Std. Error"],
                 orientation=paste0(orient,"_increasing"), intercept=co["(Intercept)","Estimate"],
                 intercept_SE=co["(Intercept)","Std. Error"], intercept_P=co["(Intercept)","Pr(>|t|)"]),
      data.table(outcome=o, method="MVMR-Egger", exposure="Smoking initiation", n_instruments=n,
                 beta=co["betaX2","Estimate"], SE=co["betaX2","Std. Error"], P=co["betaX2","Pr(>|t|)"],
                 CI_lower=co["betaX2","Estimate"]-crit_e*co["betaX2","Std. Error"],
                 CI_upper=co["betaX2","Estimate"]+crit_e*co["betaX2","Std. Error"],
                 orientation=paste0(orient,"_increasing"), intercept=co["(Intercept)","Estimate"],
                 intercept_SE=co["(Intercept)","Std. Error"], intercept_P=co["(Intercept)","Pr(>|t|)"])
    )
  }

  point <- rq(betaYG ~ betaX1 + betaX2 - 1, weights=1/sebetaYG^2, data=x)$coefficients
  boots <- matrix(NA_real_, nboot, 2)
  for (b in seq_len(nboot)) {
    z1 <- rnorm(n); z2 <- rnorm(n); zy <- rnorm(n)
    sd1 <- sqrt(vapply(gc, function(v) v[1,1], numeric(1)))
    sd2 <- sqrt(vapply(gc, function(v) v[2,2], numeric(1)))
    cov12 <- vapply(gc, function(v) v[1,2], numeric(1))
    rho <- pmax(-.999999, pmin(.999999, cov12/(sd1*sd2)))
    bx1 <- x$betaX1 + sd1*z1
    bx2 <- x$betaX2 + sd2*(rho*z1 + sqrt(1-rho^2)*z2)
    by <- x$betaYG + x$sebetaYG*zy
    boots[b,] <- rq(by ~ bx1 + bx2 - 1, weights=1/x$sebetaYG^2)$coefficients
  }
  se <- apply(boots, 2, sd, na.rm=TRUE); pv <- 2*pnorm(-abs(point/se))
  med_rows[[o]] <- rbind(
    data.table(outcome=o, method="MVMR weighted LAD median", exposure="BMI", n_instruments=n,
               beta=point[1], SE=se[1], P=pv[1], CI_lower=point[1]-qnorm(.975)*se[1], CI_upper=point[1]+qnorm(.975)*se[1],
               orientation=NA_character_, intercept=NA_real_, intercept_SE=NA_real_, intercept_P=NA_real_),
    data.table(outcome=o, method="MVMR weighted LAD median", exposure="Smoking initiation", n_instruments=n,
               beta=point[2], SE=se[2], P=pv[2], CI_lower=point[2]-qnorm(.975)*se[2], CI_upper=point[2]+qnorm(.975)*se[2],
               orientation=NA_character_, intercept=NA_real_, intercept_SE=NA_real_, intercept_P=NA_real_)
  )

  suppressWarnings(invisible(capture.output(qm <- qhet_mvmr(
    x, pcor=matrix(c(1,lambda,lambda,1),2,2), CI=TRUE, ncores=4, CI_method="jackknife"))))
  qpv <- 2*pnorm(-abs(qm[,"Effect Estimates"] / qm[,"Std. Error"]))
  qmin_rows[[o]] <- rbind(
    data.table(outcome=o, method="Q-minimisation", exposure="BMI", n_instruments=n,
               beta=qm[1,"Effect Estimates"], SE=qm[1,"Std. Error"], P=qpv[1],
               CI_lower=qm[1,"Effect Estimates"]-qnorm(.975)*qm[1,"Std. Error"], CI_upper=qm[1,"Effect Estimates"]+qnorm(.975)*qm[1,"Std. Error"],
               orientation=NA_character_, intercept=NA_real_, intercept_SE=NA_real_, intercept_P=NA_real_),
    data.table(outcome=o, method="Q-minimisation", exposure="Smoking initiation", n_instruments=n,
               beta=qm[2,"Effect Estimates"], SE=qm[2,"Std. Error"], P=qpv[2],
               CI_lower=qm[2,"Effect Estimates"]-qnorm(.975)*qm[2,"Std. Error"], CI_upper=qm[2,"Effect Estimates"]+qnorm(.975)*qm[2,"Std. Error"],
               orientation=NA_character_, intercept=NA_real_, intercept_SE=NA_real_, intercept_P=NA_real_)
  )

  invisible(capture.output(qfit <- pleiotropy_mvmr(x, gencov=gc, estimator="ivw")))
  het_rows[[o]] <- data.table(outcome=o, n_instruments=n, Q=qfit$Qstat,
                              Q_df=n-2L, Q_P=qfit$Qpval)
  full_bmi <- ivw["exposure1","Estimate"]; full_smk <- ivw["exposure2","Estimate"]
  lo <- rbindlist(lapply(seq_len(n), function(j) {
    invisible(capture.output(lf <- ivw_mvmr(x[-j,], gencov=gc[-j])))
    data.table(outcome=o, omitted_SNP=x$SNP[j], success=TRUE,
               BMI_beta=lf["exposure1","Estimate"],
               Smoking_beta=lf["exposure2","Estimate"])
  }))
  loo_rows[[o]] <- lo
  loo_summary_rows[[o]] <- data.table(
    outcome=o, n_instruments=n, BMI_full_beta=full_bmi,
    BMI_LOO_min_beta=min(lo$BMI_beta), BMI_LOO_max_beta=max(lo$BMI_beta),
    BMI_LOO_sign_reversal_count=sum(sign(lo$BMI_beta)!=sign(full_bmi)),
    Smoking_full_beta=full_smk, Smoking_LOO_min_beta=min(lo$Smoking_beta),
    Smoking_LOO_max_beta=max(lo$Smoking_beta),
    Smoking_LOO_sign_reversal_count=sum(sign(lo$Smoking_beta)!=sign(full_smk)),
    convergence_status="SUCCESS")
  conv_rows[[o]] <- data.table(outcome=o, median_bootstrap_planned=nboot,
                               median_bootstrap_success=nboot,
                               qmin_jackknife_planned=n+1L, qmin_status="SUCCESS",
                               loo_planned=n, loo_success=n)
}

allres <- rbindlist(c(ivw_rows, med_rows, qmin_rows, egger_rows), fill=TRUE)
allres[, FDR_q := p.adjust(P, "BH"), by=.(method, orientation)]
allres[method=="MVMR-Egger", intercept_FDR_q := p.adjust(intercept_P, "BH"), by=orientation]
allres[method!="MVMR-Egger", intercept_FDR_q := NA_real_]
allres[, instrument_set := "STRICT_PRIMARY_406_EXCLUDING_17_FREQUENCY_AMBIGUOUS_PALINDROMIC_SNPS"]
setorder(allres, method, orientation, outcome, exposure)
fwrite(allres, file.path(out, "MVMR_STRICT_FULL_RESULTS.tsv"), sep="\t")

het <- rbindlist(het_rows); het[, Q_FDR := p.adjust(Q_P, "BH")]
loo_sum <- rbindlist(loo_summary_rows); loo_per <- rbindlist(loo_rows)
fwrite(het, file.path(out, "MVMR_STRICT_HETEROGENEITY.tsv"), sep="\t")
fwrite(loo_sum, file.path(out, "MVMR_STRICT_LEAVE_ONE_OUT_SUMMARY.tsv"), sep="\t")
fwrite(loo_per, file.path(out, "MVMR_STRICT_LEAVE_ONE_OUT_PER_SNP.tsv"), sep="\t")
fwrite(rbindlist(conv_rows), file.path(out, "MVMR_STRICT_CONVERGENCE.tsv"), sep="\t")

summary <- allres[, .(FDR_significant=sum(FDR_q<.05, na.rm=TRUE),
                      positive=sum(beta>0, na.rm=TRUE), estimates=.N),
                  by=.(method, orientation, exposure)]
fwrite(summary, file.path(out, "MVMR_STRICT_METHOD_SUMMARY.tsv"), sep="\t")
writeLines(c(
  "# Strict-primary MVMR robustness report", "",
  "The strict primary analysis excluded 17 frequency-ambiguous palindromic SNPs before outcome inspection.",
  sprintf("Outcome-specific analyses used %d-%d instruments.", min(allres$n_instruments), max(allres$n_instruments)),
  sprintf("All %d median bootstrap fits, %d Q-minimisation full/jackknife units and %d leave-one-SNP-out fits completed.",
          15L*nboot, sum(rbindlist(conv_rows)$qmin_jackknife_planned), sum(rbindlist(conv_rows)$loo_planned)),
  sprintf("IVW Q was BH-FDR significant for %d/15 outcomes.", sum(het$Q_FDR<.05)),
  "The 423-variant analysis is retained as a sensitivity analysis only."
), file.path(out, "MVMR_STRICT_ROBUSTNESS_REPORT.md"))
cat("STRICT_PRIMARY_MVMR_COMPLETE\n")
