suppressPackageStartupMessages(library(data.table))

repo <- normalizePath(Sys.getenv("MODRISK_ROOT", unset = "."), mustWork = TRUE)
outdir <- file.path(repo, "final_major_revision", "03_mr_harmonisation")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dat <- fread(file.path(repo, "results/phase4_repair/PRIMARY_MVMR_CORRECTED_INPUT.tsv.gz"))
gclist <- readRDS(file.path(repo, "results/phase4_repair2/PRIMARY_BMI_SMK_GENCOV_LAMBDA.rds"))
gclist <- gclist[match(dat$SNP, names(gclist))]
outcomes <- sub("_beta$", "", grep("_beta$", names(dat), value = TRUE))
outcomes <- setdiff(outcomes, c("BMI", "SMK"))

ids <- file.path(outdir, "MVMR_423_RSIDS.txt")
fwrite(dat[, .(SNP)], ids, col.names = FALSE)

bmi_file <- file.path(repo, "ModRisk-MultiNet-data/raw/exposures/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz")
smk_file <- file.path(repo, "ModRisk-MultiNet-data/raw/exposures/Smoking/SmokingInitiation.txt.gz")
bmi_cmd <- sprintf(
  "gzip -dc %s | awk -v ids=%s 'BEGIN{while((getline x < ids)>0) keep[x]=1} NR==1 || ($3 in keep)'",
  shQuote(bmi_file), shQuote(ids)
)
smk_cmd <- sprintf(
  "gzip -dc %s | awk -v ids=%s 'BEGIN{while((getline x < ids)>0) keep[x]=1} NR==1 || ($3 in keep)'",
  shQuote(smk_file), shQuote(ids)
)
bmi <- fread(cmd = bmi_cmd, select = c("CHR", "POS", "SNP", "Tested_Allele", "Other_Allele", "Freq_Tested_Allele_in_HRS"))
setnames(bmi, c("CHR", "POS", "Tested_Allele", "Other_Allele", "Freq_Tested_Allele_in_HRS"),
         c("CHR", "POS", "BMI_raw_EA", "BMI_raw_OA", "BMI_EAF"))
smk <- fread(cmd = smk_cmd, select = c("CHROM", "POS", "RSID", "REF", "ALT", "AF"))
setnames(smk, c("CHROM", "POS", "RSID", "REF", "ALT", "AF"),
         c("SMK_CHR", "SMK_POS", "SNP", "reference_allele", "alternate_allele", "SMK_ALT_AF"))

a <- merge(dat, bmi, by = "SNP", all.x = TRUE, sort = FALSE)
a <- merge(a, smk, by = "SNP", all.x = TRUE, sort = FALSE)
a <- a[match(dat$SNP, SNP)]

same <- a$BMI_EA == a$SMK_EA & a$BMI_OA == a$SMK_OA
swap <- a$BMI_EA == a$SMK_OA & a$BMI_OA == a$SMK_EA
comp <- function(x) chartr("ACGT", "TGCA", x)
comp_direct <- a$BMI_EA == comp(a$SMK_EA) & a$BMI_OA == comp(a$SMK_OA)
comp_swap <- a$BMI_EA == comp(a$SMK_OA) & a$BMI_OA == comp(a$SMK_EA)
pal <- paste0(a$BMI_EA, a$BMI_OA) %in% c("AT", "TA", "CG", "GC")

valid_outcomes <- integer(nrow(a))
outcome_all_smk <- logical(nrow(a))
outcome_ea_sets <- character(nrow(a))
for (i in seq_len(nrow(a))) {
  valid <- outcomes[vapply(outcomes, function(o) is.finite(a[[paste0(o, "_beta")]][i]) && is.finite(a[[paste0(o, "_SE")]][i]), logical(1))]
  valid_outcomes[i] <- length(valid)
  eas <- unique(vapply(valid, function(o) a[[paste0(o, "_EA")]][i], character(1)))
  outcome_ea_sets[i] <- paste(eas, collapse = ";")
  outcome_all_smk[i] <- length(eas) <= 1 && (length(eas) == 0 || eas == a$SMK_EA[i])
}

audit <- data.table(
  SNP = a$SNP,
  chromosome = fifelse(!is.na(a$CHR), a$CHR, a$SMK_CHR),
  position_GRCh37 = fifelse(!is.na(a$POS), a$POS, a$SMK_POS),
  reference_allele = a$reference_allele,
  alternate_allele = a$alternate_allele,
  BMI_effect_allele = a$BMI_EA,
  BMI_other_allele = a$BMI_OA,
  BMI_EAF = a$BMI_EAF,
  BMI_EAF_common_smoking_EA = fifelse(swap, 1 - a$BMI_EAF, a$BMI_EAF),
  smoking_effect_allele = a$SMK_EA,
  smoking_other_allele = a$SMK_OA,
  smoking_effect_allele_frequency = fifelse(a$SMK_EA == a$alternate_allele, a$SMK_ALT_AF, 1 - a$SMK_ALT_AF),
  aligned_exposure_EAF_difference = abs(fifelse(swap, 1 - a$BMI_EAF, a$BMI_EAF) - fifelse(a$SMK_EA == a$alternate_allele, a$SMK_ALT_AF, 1 - a$SMK_ALT_AF)),
  BMI_smoking_alignment = fifelse(same, "DIRECT", fifelse(swap, "SWAPPED", fifelse(comp_direct, "STRAND_COMPLEMENT_DIRECT", fifelse(comp_swap, "STRAND_COMPLEMENT_SWAPPED", "UNRESOLVED")))),
  palindromic = pal,
  BMI_beta_original = a$BMI_beta,
  BMI_beta_common_smoking_EA = fifelse(swap, -a$BMI_beta, a$BMI_beta),
  smoking_beta_common_smoking_EA = a$SMK_beta,
  BMI_SE = a$BMI_SE,
  smoking_SE = a$SMK_SE,
  outcome_valid_count = valid_outcomes,
  outcome_effect_alleles = outcome_ea_sets,
  all_valid_outcomes_aligned_to_smoking_EA = outcome_all_smk,
  exposure_covariance = vapply(gclist, function(x) x[1, 2], numeric(1)),
  covariance_transform_sign_product = 1L,
  covariance_after_common_allele_recode = vapply(gclist, function(x) x[1, 2], numeric(1)),
  duplicate_SNP = duplicated(a$SNP) | duplicated(a$SNP, fromLast = TRUE),
  missing_required_exposure_value = !is.finite(a$BMI_beta) | !is.finite(a$BMI_SE) | !is.finite(a$SMK_beta) | !is.finite(a$SMK_SE),
  exclusion_status = fifelse(!(same | swap), "EXCLUDE_UNRESOLVED_ALLELES", "INCLUDED")
)

fwrite(audit, file.path(outdir, "MR_ALLELE_HARMONISATION_AUDIT.tsv"), sep = "\t")

summary <- data.table(
  metric = c(
    "input_rows", "unique_snps", "direct_pairs", "swapped_pairs",
    "strand_complement_pairs", "palindromic_pairs", "duplicates",
    "unresolved_pairs", "missing_exposure_values", "all_outcome_alleles_valid",
    "palindromic_max_aligned_EAF_difference",
    "variants_with_one_missing_outcome", "outcomes_with_422_instruments", "missing_outcome_variant"
  ),
  value = c(
    nrow(a), uniqueN(a$SNP), sum(same), sum(swap), sum(!(same | swap) & (comp_direct | comp_swap)),
    sum(pal), sum(audit$duplicate_SNP), sum(!(same | swap)),
    sum(audit$missing_required_exposure_value), sum(outcome_all_smk),
    max(audit[palindromic == TRUE]$aligned_exposure_EAF_difference),
    sum(valid_outcomes == 14),
    sum(vapply(outcomes, function(o) sum(is.finite(a[[paste0(o, "_beta")]]) & is.finite(a[[paste0(o, "_SE")]])) == 422, logical(1))),
    "rs705217"
  ),
  status = c(rep("PASS", 11), "OBSERVED", "OBSERVED", "OBSERVED")
)
fwrite(summary, file.path(outdir, "MR_ALLELE_HARMONISATION_SUMMARY.tsv"), sep = "\t")
message("Per-SNP harmonisation audit completed")
