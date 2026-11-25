suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

repo <- normalizePath(Sys.getenv("MODRISK_ROOT", unset = "."), mustWork = TRUE)
outdir <- file.path(repo, "final_major_revision", "04_univariable_mr")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ref_prefix <- file.path(repo, "data", "reference", "ld_clumping", "EUR")
bmi_file <- file.path(repo, "ModRisk-MultiNet-data", "raw", "exposures", "Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz")
smk_file <- file.path(repo, "ModRisk-MultiNet-data", "raw", "exposures", "Smoking", "SmokingInitiation.txt.gz")
outcome_dir <- file.path(repo, "ModRisk-MultiNet-data", "raw", "GEMINI")

outcomes <- c(
  "COPD", "HF", "OA", "alcohol_problems", "carpal_tunnel", "coronary_heart",
  "diabetes_type_2", "fibromyalgia", "gastritis_duodenitis", "gout",
  "hypertension", "insomnia", "periph_neuro", "sleep_apnoea", "valve_nonrheumatic"
)

comp <- function(x) chartr("ACGT", "TGCA", toupper(x))
is_pal <- function(a, b) paste0(toupper(a), toupper(b)) %in% c("AT", "TA", "CG", "GC")

message("Reading complete BMI and smoking-initiation GWAS")
bmi <- fread(bmi_file, select = c("CHR", "POS", "SNP", "Tested_Allele", "Other_Allele", "Freq_Tested_Allele_in_HRS", "BETA", "SE", "P", "N"))
setnames(bmi, c("CHR", "POS", "SNP", "Tested_Allele", "Other_Allele", "Freq_Tested_Allele_in_HRS", "BETA", "SE", "P", "N"),
         c("CHR", "POS", "SNP", "EA", "OA", "EAF", "beta", "SE", "P", "N"))
bmi[, exposure := "BMI"]
bmi <- bmi[is.finite(P) & P < 5e-8 & is.finite(beta) & is.finite(SE) & SE > 0 & !is.na(SNP) & SNP != ""]
bmi <- bmi[order(P)][!duplicated(SNP)]

smk <- fread(smk_file, select = c("CHROM", "POS", "RSID", "REF", "ALT", "AF", "BETA", "SE", "PVALUE", "N", "EFFECTIVE_N"))
setnames(smk, c("CHROM", "POS", "RSID", "REF", "ALT", "AF", "BETA", "SE", "PVALUE", "N", "EFFECTIVE_N"),
         c("CHR", "POS", "SNP", "OA", "EA", "EAF", "beta", "SE", "P", "N", "effective_N"))
smk[, exposure := "Smoking initiation"]
smk <- smk[is.finite(P) & P < 5e-8 & is.finite(beta) & is.finite(SE) & SE > 0 & !is.na(SNP) & SNP != ""]
smk <- smk[order(P)][!duplicated(SNP)]

run_clump <- function(d, tag) {
  inp <- file.path(outdir, paste0(tag, "_GWS_CLUMP_INPUT.tsv"))
  fwrite(d[, .(ID = SNP, P)], inp, sep = "\t")
  prefix <- file.path(outdir, paste0(tag, "_independent"))
  args <- c(
    "--bfile", ref_prefix, "--clump", inp,
    "--clump-id-field", "ID", "--clump-p-field", "P",
    "--clump-p1", "5e-8", "--clump-p2", "1",
    "--clump-r2", "0.001", "--clump-kb", "10000",
    "--out", prefix
  )
  status <- system2("plink2", args = args, stdout = file.path(outdir, paste0(tag, "_plink_stdout.log")), stderr = file.path(outdir, paste0(tag, "_plink_stderr.log")))
  if (status != 0) stop("PLINK clumping failed for ", tag)
  cl <- fread(paste0(prefix, ".clumps"))
  ids <- cl$ID
  list(data = d[SNP %in% ids][match(ids, SNP)], clumps = cl, prefix = prefix)
}

bmi_cl <- run_clump(bmi, "BMI")
smk_cl <- run_clump(smk, "SMK")
instruments <- rbindlist(list(bmi_cl$data, smk_cl$data), fill = TRUE)
instruments[, `:=`(EA = toupper(EA), OA = toupper(OA))]
fwrite(instruments, file.path(outdir, "INDEPENDENT_EXPOSURE_INSTRUMENTS.tsv"), sep = "\t")

union_ids <- unique(instruments$SNP)
id_file <- file.path(outdir, "independent_union_rsids.txt")
fwrite(data.table(SNP = union_ids), id_file, col.names = FALSE)

message("Extracting selected variants from the 15 GEMINI outcome files")
outcome_data <- list()
for (o in outcomes) {
  src <- file.path(outcome_dir, paste0("GEMINI.sumstats.", o, ".txt.gz"))
  # rs_id is field 9 in the frozen GEMINI v1.1 files. The header is retained.
  awk <- sprintf("awk -v ids=%s 'BEGIN{while((getline line < ids)>0){keep[line]=1}} NR==1 || ($9 in keep)'", shQuote(id_file))
  cmd <- sprintf("gzip -dc %s | %s", shQuote(src), awk)
  d <- fread(cmd = cmd, select = c("effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "rs_id"))
  setnames(d, c("effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "rs_id"),
           c("out_EA", "out_OA", "out_beta", "out_SE", "out_EAF", "SNP"))
  d[, outcome := o]
  outcome_data[[o]] <- d[!duplicated(SNP)]
}
outcome_all <- rbindlist(outcome_data, use.names = TRUE)

harmonise_one <- function(exp_d, out_d) {
  d <- merge(exp_d, out_d, by = "SNP", all.x = TRUE, sort = FALSE)
  d[, `:=`(EA = toupper(EA), OA = toupper(OA), out_EA = toupper(out_EA), out_OA = toupper(out_OA))]
  d[, `:=`(
    direct = EA == out_EA & OA == out_OA,
    swapped = EA == out_OA & OA == out_EA,
    complement_direct = EA == comp(out_EA) & OA == comp(out_OA),
    complement_swapped = EA == comp(out_OA) & OA == comp(out_EA),
    palindromic = is_pal(EA, OA)
  )]
  d[, harmonisation_action := fifelse(direct, "DIRECT",
    fifelse(swapped, "SWAPPED",
      fifelse(complement_direct, "COMPLEMENT_DIRECT",
        fifelse(complement_swapped, "COMPLEMENT_SWAPPED", "ALLELE_MISMATCH"))))]
  d[, outcome_flip := harmonisation_action %in% c("SWAPPED", "COMPLEMENT_SWAPPED")]
  d[, frequency_difference := fifelse(outcome_flip, abs(EAF - (1 - out_EAF)), abs(EAF - out_EAF))]
  d[, exclusion_reason := NA_character_]
  d[is.na(out_beta) | is.na(out_SE), exclusion_reason := "OUTCOME_VARIANT_MISSING"]
  d[harmonisation_action == "ALLELE_MISMATCH" & is.na(exclusion_reason), exclusion_reason := "ALLELE_MISMATCH"]
  d[palindromic & abs(EAF - 0.5) <= 0.08 & is.na(exclusion_reason), exclusion_reason := "PALINDROMIC_INTERMEDIATE_FREQUENCY"]
  d[palindromic & frequency_difference > 0.10 & is.na(exclusion_reason), exclusion_reason := "PALINDROMIC_FREQUENCY_DISCORDANCE"]
  d[!palindromic & frequency_difference > 0.20 & is.na(exclusion_reason), exclusion_reason := "FREQUENCY_DISCORDANCE"]
  d[!is.finite(out_SE) | out_SE <= 0, exclusion_reason := "INVALID_OUTCOME_SE"]
  d[, included := is.na(exclusion_reason)]
  d[included & outcome_flip, out_beta := -out_beta]
  d
}

methods <- list(
  IVW = mr_ivw,
  `MR-Egger` = mr_egger_regression,
  `Weighted median` = mr_weighted_median,
  `Weighted mode` = mr_weighted_mode
)
minimum_n <- c(IVW = 2L, `MR-Egger` = 3L, `Weighted median` = 3L, `Weighted mode` = 3L)

audit_rows <- list()
res_rows <- list()
k <- 1L
a <- 1L
for (ex in c("BMI", "Smoking initiation")) {
  exd <- instruments[exposure == ex]
  for (o in outcomes) {
    h <- harmonise_one(exd, outcome_all[outcome == o])
    h[, `:=`(exposure = ex, outcome = o)]
    audit_rows[[a]] <- h
    a <- a + 1L
    use <- h[included == TRUE]
    mean_f <- if (nrow(use)) mean((use$beta / use$SE)^2) else NA_real_
    for (m in names(methods)) {
      if (nrow(use) < minimum_n[[m]]) {
        res_rows[[k]] <- data.table(
          exposure = ex, outcome = o, instrument_set = "INDEPENDENT_EXPOSURE_UVMR", method = m,
          nsnp = nrow(use), mean_F = mean_f, beta = NA_real_, SE = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_,
          P = NA_real_, intercept = NA_real_, intercept_SE = NA_real_, intercept_P = NA_real_,
          exposure_unit = ifelse(ex == "BMI", "per 1-SD higher inverse-normal-transformed BMI residual", "per unit higher genetically predicted smoking-initiation liability (ever versus never regular smoking)"),
          outcome_scale = "log-odds of disease", convergence_status = "NOT_ESTIMABLE",
          warning = "Insufficient instruments for estimator", exclusion_reason = NA_character_
        )
      } else {
        set.seed(20261124 + k)
        warning_text <- character()
        fit <- withCallingHandlers(
          tryCatch(methods[[m]](use$beta, use$out_beta, use$SE, use$out_SE, default_parameters()), error = function(e) e),
          warning = function(w) { warning_text <<- c(warning_text, conditionMessage(w)); invokeRestart("muffleWarning") }
        )
        if (inherits(fit, "error")) {
          res_rows[[k]] <- data.table(
            exposure = ex, outcome = o, instrument_set = "INDEPENDENT_EXPOSURE_UVMR", method = m,
            nsnp = nrow(use), mean_F = mean_f, beta = NA_real_, SE = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_,
            P = NA_real_, intercept = NA_real_, intercept_SE = NA_real_, intercept_P = NA_real_,
            exposure_unit = ifelse(ex == "BMI", "per 1-SD higher inverse-normal-transformed BMI residual", "per unit higher genetically predicted smoking-initiation liability (ever versus never regular smoking)"),
            outcome_scale = "log-odds of disease", convergence_status = "FAILED",
            warning = conditionMessage(fit), exclusion_reason = NA_character_
          )
        } else {
          res_rows[[k]] <- data.table(
            exposure = ex, outcome = o, instrument_set = "INDEPENDENT_EXPOSURE_UVMR", method = m,
            nsnp = fit$nsnp, mean_F = mean_f, beta = fit$b, SE = fit$se,
            CI_lower = fit$b - 1.96 * fit$se, CI_upper = fit$b + 1.96 * fit$se, P = fit$pval,
            intercept = ifelse(m == "MR-Egger", fit$b_i, NA_real_),
            intercept_SE = ifelse(m == "MR-Egger", fit$se_i, NA_real_),
            intercept_P = ifelse(m == "MR-Egger", fit$pval_i, NA_real_),
            exposure_unit = ifelse(ex == "BMI", "per 1-SD higher inverse-normal-transformed BMI residual", "per unit higher genetically predicted smoking-initiation liability (ever versus never regular smoking)"),
            outcome_scale = "log-odds of disease", convergence_status = "SUCCESS",
            warning = paste(unique(warning_text), collapse = " | "), exclusion_reason = NA_character_
          )
        }
      }
      k <- k + 1L
    }
  }
}

audit <- rbindlist(audit_rows, fill = TRUE)
res <- rbindlist(res_rows, fill = TRUE)
res[, FDR_q := p.adjust(P, method = "BH"), by = .(exposure, method, instrument_set)]
setcolorder(res, c("exposure", "outcome", "instrument_set", "method", "nsnp", "mean_F", "beta", "SE", "CI_lower", "CI_upper", "P", "FDR_q", "intercept", "intercept_SE", "intercept_P", "exposure_unit", "outcome_scale", "convergence_status", "warning", "exclusion_reason"))

fwrite(audit, file.path(outdir, "UVMR_HARMONISATION_AUDIT.tsv"), sep = "\t")
fwrite(res, file.path(outdir, "UNIVARIABLE_MR_INDEPENDENT.tsv"), sep = "\t")

# Preserve the earlier union-subset sensitivity analysis as a separate instrument set.
union_src <- fread(file.path(repo, "revision_round_20261124", "analysis", "outputs", "UNIVARIABLE_MR_SENSITIVITY.tsv"))
union_strength <- fread(file.path(repo, "revision_round_20261124", "analysis", "outputs", "UNIVARIABLE_INSTRUMENT_STRENGTH.tsv"))
union_src <- merge(
  union_src,
  union_strength[, .(exposure, mean_F)],
  by = "exposure", all.x = TRUE
)
union <- union_src[, .(
  exposure, outcome, instrument_set = "UNION_SUBSET_UVMR", method,
  nsnp = n_instruments, mean_F, beta, SE, CI_lower = beta - 1.96 * SE,
  CI_upper = beta + 1.96 * SE, P, FDR_q = q, intercept, intercept_SE, intercept_P,
  exposure_unit = fifelse(exposure == "BMI", "per 1-SD higher inverse-normal-transformed BMI residual", "per unit higher genetically predicted smoking-initiation liability (ever versus never regular smoking)"),
  outcome_scale = "log-odds of disease", convergence_status = "SUCCESS",
  warning = "Historical union-subset sensitivity; q values retained within exposure-by-method family",
  exclusion_reason = NA_character_
)]
fwrite(union, file.path(outdir, "UNIVARIABLE_MR_UNION_SUBSET.tsv"), sep = "\t")
full <- rbindlist(list(union, res), fill = TRUE)
fwrite(full, file.path(outdir, "UVMR_FULL_RESULTS.tsv"), sep = "\t")

comparison <- rbindlist(list(
  instruments[, .(gws_preclump = .N, post_clump = .N), by = exposure][, gws_preclump := c(nrow(bmi), nrow(smk))],
  union_strength[, .(exposure, union_subset_n = n_instruments, union_subset_mean_F = mean_F)]
), fill = TRUE)
comparison <- merge(
  data.table(exposure = c("BMI", "Smoking initiation"),
             full_gws_preclump = c(nrow(bmi), nrow(smk)),
             independent_postclump = c(nrow(bmi_cl$data), nrow(smk_cl$data))),
  union_strength[, .(exposure, union_subset_n = n_instruments, union_subset_mean_F = mean_F)],
  by = "exposure", all = TRUE
)
comparison[, mr_ahc_n := fifelse(exposure == "Smoking initiation", 54L, NA_integer_)]
comparison[, explanation := fifelse(
  exposure == "Smoking initiation",
  "The 54-variant MR-AHC set came from the BMI-smoking union after joint clumping and required smoking-specific per-variant F>10; only 27 of those variants were smoking genome-wide significant. The 16-variant historical UVMR subset additionally required smoking genome-wide significance and complete BMI-smoking-outcome data in the corrected 423-SNP matrix. The independent exposure-specific UVMR starts from the full smoking GWAS, clumps smoking loci alone and does not require BMI completeness.",
  "Independent BMI instruments were selected directly from the full BMI GWAS; the historical subset was restricted to variants retained in the complete BMI-smoking MVMR matrix."
)]
fwrite(comparison, file.path(outdir, "UVMR_INSTRUMENT_COMPARISON.tsv"), sep = "\t")

key <- res[outcome %in% c("HF", "coronary_heart", "diabetes_type_2", "COPD", "OA") & method == "IVW"]
report <- c(
  "# Univariable MR sensitivity report", "",
  "## Design", "",
  "Independent exposure-specific instruments were selected from the complete BMI and smoking-initiation GWAS files at P < 5 × 10^-8 and clumped with the frozen 1000 Genomes European reference at r2 < 0.001 within 10 Mb. Outcome effects were independently harmonised for all 15 outcomes. Palindromic intermediate-frequency variants and unresolved allele mismatches were excluded with explicit reasons.", "",
  "The primary multiplicity family for the new analysis is 15 outcomes within each exposure-by-method combination. The historical union-subset q values are retained separately and are not mixed with the independent-analysis q values.", "",
  "## Instrument counts", "",
  paste(capture.output(print(comparison)), collapse = "\n"), "",
  "## Priority outcomes under independent IVW", "",
  paste(capture.output(print(key[, .(exposure, outcome, nsnp, mean_F, beta, SE, P, FDR_q)])), collapse = "\n"), "",
  "## Interpretation", "",
  "Univariable MR estimates total exposure associations and does not estimate the same mutually adjusted direct effects as MVMR. Agreement or disagreement is therefore informative about estimand and instrument-set dependence, not a formal validation test. The independent smoking set may equal the 54-variant MR-AHC set if the same source and clumping rules yield the same retained variants, but this was not imposed.", ""
)
writeLines(report, file.path(outdir, "UVMR_SENSITIVITY_REPORT.md"))

message("INDEPENDENT_UVMR_COMPLETE: BMI=", nrow(bmi_cl$data), ", SMK=", nrow(smk_cl$data))
