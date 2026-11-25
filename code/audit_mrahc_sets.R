suppressPackageStartupMessages(library(data.table))

root <- normalizePath(Sys.getenv("MODRISK_ROOT", unset = "."), mustWork = TRUE)
out <- file.path(root, "final_revision_closeout", "01_mr_ahc_provenance")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) {
  z <- system2("shasum", c("-a", "256", path), stdout = TRUE)
  sub(" .*", "", z[[1]])
}

six_pre <- fread(file.path(root, "results/phase4/MR_INSTRUMENT_UNION_PRECLUMP.tsv"))
six <- fread(file.path(root, "results/phase4/MR_INSTRUMENT_UNION_CLUMPED.tsv"))
exp_mat <- fread(file.path(root, "results/phase4/MVMR_EXPOSURE_MATRIX.tsv"))
mvmr <- fread(file.path(root, "results/phase4_repair/PRIMARY_MVMR_CORRECTED_INPUT.tsv.gz"))
mrahc <- fread(file.path(root, "results/phase5/MRAHC_PRIMARY_INPUT.tsv"))
indep <- fread(file.path(root, "final_major_revision/04_univariable_mr/INDEPENDENT_EXPOSURE_INSTRUMENTS.tsv"))
uv_audit <- fread(file.path(root, "final_major_revision/04_univariable_mr/UVMR_HARMONISATION_AUDIT.tsv"))

# Reconstruct the BMI-smoking-only union using the same source union, LD reference,
# clumping threshold and PLINK version as the historical six-exposure clump. This
# is a new comparator, not a recovered historical file.
bms_pre <- six_pre[grepl("(^|;)(BMI|SMK)(;|$)", exposures_involved)]
clump_in <- bms_pre[, .(ID = rsID, P = min_P)]
clump_input <- file.path(out, "BMI_SMK_RECONSTRUCTED_CLUMP_INPUT.tsv")
fwrite(clump_in, clump_input, sep = "\t")
prefix <- file.path(out, "BMI_SMK_RECONSTRUCTED")
args <- c(
  "--bfile", file.path(root, "data/reference/ld_clumping/EUR"),
  "--clump", clump_input,
  "--clump-id-field", "ID", "--clump-p-field", "P",
  "--clump-p1", "5e-8", "--clump-p2", "1",
  "--clump-r2", "0.001", "--clump-kb", "10000",
  "--out", prefix
)
status <- system2("plink2", args, stdout = file.path(out, "BMI_SMK_RECONSTRUCTED.stdout.log"),
                  stderr = file.path(out, "BMI_SMK_RECONSTRUCTED.stderr.log"))
stopifnot(status == 0L)
bms_clumps <- fread(paste0(prefix, ".clumps"))
bms_ids <- unique(bms_clumps$ID)

pvals <- exp_mat[, .(SNP, BMI_beta, BMI_SE, BMI_P, SMK_beta, SMK_SE, SMK_P,
                     BMI_EA, BMI_OA, SMK_EA, SMK_OA)]
mvmr_p <- merge(mvmr[, .(SNP)], pvals[, .(SNP, BMI_P, SMK_P)], by = "SNP", all.x = TRUE)
hist_smk_ids <- mvmr_p[SMK_P < 5e-8, SNP]
mrahc_ids <- unique(mrahc$rsID)
indep_smk_ids <- indep[exposure == "Smoking initiation", unique(SNP)]
usable_counts <- uv_audit[exposure == "Smoking initiation" & included == TRUE,
                          .(n_outcomes = uniqueN(outcome)), by = SNP]
usable_smk_ids <- usable_counts[n_outcomes == uniqueN(uv_audit$outcome), SNP]

sets <- list(
  SIX_EXPOSURE_CLUMPED_UNION = unique(six$rsID),
  BMI_SMK_RECONSTRUCTED_CLUMPED_UNION = bms_ids,
  MVMR_COMPLETE_423 = unique(mvmr$SNP),
  HISTORICAL_SMK_UVMR_16 = unique(hist_smk_ids),
  MRAHC_54 = mrahc_ids,
  INDEPENDENT_SMK_CLUMPED_92 = indep_smk_ids,
  INDEPENDENT_SMK_OUTCOME_USABLE_84 = usable_smk_ids
)

all_ids <- sort(unique(unlist(sets, use.names = FALSE)))
prov <- data.table(rsID = all_ids)
for (nm in names(sets)) prov[, (nm) := rsID %chin% sets[[nm]]]
prov <- merge(prov, six[, .(rsID, six_union_source_axes = exposures_involved,
                            six_union_min_P = min_P, six_union_n_exposures = n_exposures)],
              by = "rsID", all.x = TRUE)
prov <- merge(prov, pvals, by.x = "rsID", by.y = "SNP", all.x = TRUE)
prov[, `:=`(
  BMI_F = (BMI_beta / BMI_SE)^2,
  SMK_F = (SMK_beta / SMK_SE)^2,
  BMI_GWS = is.finite(BMI_P) & BMI_P < 5e-8,
  SMK_GWS = is.finite(SMK_P) & SMK_P < 5e-8
)]
prov[, `:=`(
  MRAHC_selection_rule_verified = MRAHC_54 & SMK_F > 10,
  MRAHC_complete_four_outcomes = MRAHC_54 & rsID %chin% mrahc$rsID,
  provenance_note = fifelse(
    MRAHC_54,
    "Historical six-exposure clumped-union variant with smoking-specific F>10; four outcome profiles complete",
    "See set-membership columns"
  )
)]

source_paths <- c(
  "results/phase4/MR_INSTRUMENT_UNION_PRECLUMP.tsv",
  "results/phase4/MR_INSTRUMENT_UNION_CLUMPED.tsv",
  "results/phase4/MVMR_EXPOSURE_MATRIX.tsv",
  "results/phase4_repair/PRIMARY_MVMR_CORRECTED_INPUT.tsv.gz",
  "results/phase5/MRAHC_PRIMARY_INPUT.tsv",
  "final_major_revision/04_univariable_mr/INDEPENDENT_EXPOSURE_INSTRUMENTS.tsv",
  "final_major_revision/04_univariable_mr/UVMR_HARMONISATION_AUDIT.tsv"
)
source_manifest <- data.table(
  source_file = source_paths,
  sha256 = vapply(file.path(root, source_paths), sha256, character(1))
)
source_manifest[, source_role := c(
  "six-exposure pre-clump union", "six-exposure clumped union", "exposure effects and P values",
  "complete 423-SNP MVMR matrix", "54-SNP MR-AHC profile input",
  "independent exposure-specific clumps", "outcome harmonisation and usability"
)]
fwrite(source_manifest, file.path(out, "MR_AHC_SOURCE_FILE_MANIFEST.tsv"), sep = "\t")
fwrite(prov, file.path(out, "MR_AHC_INSTRUMENT_PROVENANCE.tsv"), sep = "\t")

overlap <- rbindlist(lapply(names(sets), function(a) {
  rbindlist(lapply(names(sets), function(b) data.table(
    set_a = a, set_b = b, n_a = length(sets[[a]]), n_b = length(sets[[b]]),
    intersection_n = length(intersect(sets[[a]], sets[[b]])),
    a_minus_b_n = length(setdiff(sets[[a]], sets[[b]])),
    b_minus_a_n = length(setdiff(sets[[b]], sets[[a]])),
    a_subset_b = all(sets[[a]] %chin% sets[[b]]),
    b_subset_a = all(sets[[b]] %chin% sets[[a]])
  )))
}))
fwrite(overlap, file.path(out, "MR_AHC_SET_OVERLAP.tsv"), sep = "\t")

stopifnot(length(sets$SIX_EXPOSURE_CLUMPED_UNION) == 633L)
stopifnot(length(sets$MVMR_COMPLETE_423) == 423L)
stopifnot(length(sets$HISTORICAL_SMK_UVMR_16) == 16L)
stopifnot(length(sets$MRAHC_54) == 54L)
stopifnot(length(sets$INDEPENDENT_SMK_CLUMPED_92) == 92L)
stopifnot(length(sets$INDEPENDENT_SMK_OUTCOME_USABLE_84) == 84L)
stopifnot(all(sets$MRAHC_54 %chin% sets$SIX_EXPOSURE_CLUMPED_UNION))
stopifnot(sum(prov$MRAHC_54 & prov$SMK_GWS, na.rm = TRUE) == 27L)
stopifnot(sum(prov$MRAHC_54 & prov$MRAHC_selection_rule_verified, na.rm = TRUE) == 54L)

summary <- c(
  "# MR-AHC historical settings audit",
  "",
  "## Verified inputs and selection",
  "",
  sprintf("- The historical six-exposure clumped union contains %d SNPs.", length(sets$SIX_EXPOSURE_CLUMPED_UNION)),
  sprintf("- The newly reconstructed BMI-smoking-only union contains %d SNPs; it is a comparator and was not substituted for the historical six-exposure union.", length(sets$BMI_SMK_RECONSTRUCTED_CLUMPED_UNION)),
  "- All 54 MR-AHC SNPs are members of the historical six-exposure clumped union and are exactly the 54 union variants with smoking-specific F > 10 and complete four-outcome profile rows.",
  sprintf("- Six of the 54 entered the union only through non-BMI/non-smoking axes; therefore the historical input cannot be described as a BMI-smoking union."),
  "- Twenty-seven of 54 were genome-wide significant for smoking initiation. The remaining 27 were selected by another exposure axis but still had smoking-specific F > 10.",
  "- Forty of 54 are retained in the complete 423-SNP BMI-smoking MVMR matrix. All 16 smoking-GWS variants in that complete matrix are among the 54 MR-AHC variants.",
  "- The 27-versus-16 difference is therefore explained by complete-matrix filtering: 11 smoking-GWS MR-AHC variants were not retained in the 423-SNP complete MVMR matrix. It is not a contradictory significance threshold.",
  "",
  "## Historical analysis settings",
  "",
  "- The primary RDS records Cluster_number = 1, Cluster_number_real = 0, no accepted AHC cluster, all 54 variants in the junk set, and F = 28.2.",
  "- The project records state that cluster acceptance followed the published Q-based stopping rule and that a four-outcome extension of the published one/two-outcome implementation was used.",
  "- The extension source code, complete per-iteration bootstrap output, bootstrap iteration count and random seed are absent from the preserved historical assets.",
  "- The recorded co-clustering summaries (median 0, mean 0.064, 2.7% of pairs >0.8) cannot be independently regenerated from the retained files.",
  "- The upstream reference example was reproduced with R 4.4.3 from commit 9b66344, but this does not validate the missing four-outcome extension or its bootstrap implementation.",
  "",
  "## Audit conclusion",
  "",
  "Instrument provenance is fully resolved. Historical analytical provenance is incomplete because the extension and bootstrap settings are not recoverable. The numerical stability claims are therefore not independently reproducible."
)
writeLines(summary, file.path(out, "MR_AHC_HISTORICAL_SETTINGS_AUDIT.md"))

disposition <- c(
  "# MR-AHC final disposition",
  "",
  "**MR_AHC_STATUS = WITHDRAWN_WITH_JUSTIFICATION**",
  "",
  "The 54-variant source is now verified as the historical six-exposure clumped union restricted to variants with smoking-specific F > 10 and complete four-outcome profiles. This resolves the instrument-set contradiction. However, the project does not preserve the source code for the four-outcome extension, the bootstrap iteration count, the random seed or the per-iteration co-clustering output. Consequently, the recorded no-stable-cluster result and stability percentages cannot be independently reproduced to the standard required for a current result.",
  "",
  "MR-AHC is therefore withdrawn from the current evidence set, manuscript results, supplementary figures and machine-readable result claims. The response letter will state the verified provenance and the reason for withdrawal. Historical files remain untouched as auditable provenance and are not presented as valid current results. No replacement clustering analysis is introduced because MR-AHC is secondary to the primary network findings."
)
writeLines(disposition, file.path(out, "MR_AHC_FINAL_DISPOSITION.md"))

cat("MR-AHC set audit complete\n")
