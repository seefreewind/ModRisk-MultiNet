suppressPackageStartupMessages(library(data.table))

root <- normalizePath(Sys.getenv("MODRISK_ROOT", unset = "."), mustWork = TRUE)
out <- file.path(root, "final_revision_closeout", "04_machine_readable_audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

pair_path <- file.path(root, "final_major_revision/07_figures_tables/Additional_file_2_Pair_level_classification_Final.csv")
uv_path <- file.path(root, "final_major_revision/04_univariable_mr/UVMR_FULL_RESULTS.tsv")
mv_path <- file.path(root, "final_major_revision/05_robust_mvmr/MVMR_FULL_RESULTS.tsv")
pair <- fread(pair_path)
uv <- fread(uv_path)
mv <- fread(mv_path)
het <- fread(file.path(root, "final_revision_closeout/02_mr_diagnostics/MVMR_HETEROGENEITY_FINAL.tsv"))
loo <- fread(file.path(root, "final_revision_closeout/02_mr_diagnostics/MVMR_LEAVE_ONE_OUT_FINAL.tsv"))
loo_per_snp <- fread(file.path(root, "final_revision_closeout/02_mr_diagnostics/MVMR_LEAVE_ONE_OUT_PER_SNP.tsv"))

pair_checks <- data.table(
  metric = c("unique_pairs","baseline_FDR_significant","joint_signed_difference_FDR_significant",
             "attenuated","enhanced","no_significant_signed_change","residual_FDR_significant",
             "near_null","other_residual_non_significant","sign_reversals","significant_reversals",
             "CI_within_0.05","CI_within_0.10","legacy_residual_Z_gate"),
  observed = c(uniqueN(pair$pair_id), sum(pair$baseline_fdr < .05), sum(pair$difference_fdr < .05),
               sum(pair$change_class == "ATTENUATED"), sum(pair$change_class == "ENHANCED"),
               sum(pair$change_class == "NO_SIGNIFICANT_CHANGE"), sum(pair$residual_significant),
               sum(pair$near_null_classified), sum(pair$residual_state == "RESIDUAL_NON_SIGNIFICANT_OTHER"),
               sum(pair$sign_reversal_flag), sum(pair$significant_reversal_flag),
               sum(pair$equivalent_005), sum(pair$equivalent_010),
               sum(pair$legacy_significant_reversal_residual_z_gate)),
  expected = c(2485,1535,1019,900,119,1466,1100,605,780,135,66,2,162,1)
)
pair_checks[, status := fifelse(observed == expected, "PASS", "FAIL")]
flag_expected <- pair$sign_reversal_flag & pair$difference_fdr < .05
flag_match <- all(pair$significant_reversal_flag == flag_expected)
mapping_ok <- all(pair[trait_1 == "alc_problems", trait_1_resource_code] == "alcohol_problems") &&
  all(pair[trait_2 == "alc_problems", trait_2_resource_code] == "alcohol_problems") &&
  all(pair[trait_1 == "deaf", trait_1_resource_code] == "hearing_loss") &&
  all(pair[trait_2 == "deaf", trait_2_resource_code] == "hearing_loss") &&
  all(pair[trait_1 == "schizo", trait_1_resource_code] == "SSD") &&
  all(pair[trait_2 == "schizo", trait_2_resource_code] == "SSD")
fwrite(pair_checks, file.path(out, "FINAL_PAIR_CLASSIFICATION_VALIDATION.tsv"), sep="\t")

uv_key <- uv[, .N, by=.(instrument_set, exposure, method, outcome)]
uv_complete <- nrow(uv) == 240L && nrow(uv_key) == 240L && all(uv_key$N == 1L)
uv_success <- all(uv$convergence_status == "SUCCESS")
mv_key <- mv[, .N, by=.(method, exposure, outcome, orientation)]
mv_complete <- nrow(mv) == 150L && nrow(mv_key) == 150L && all(mv_key$N == 1L)
egger_intercepts <- unique(mv[method == "MVMR-Egger", .(outcome, orientation, intercept, intercept_SE, intercept_P, intercept_FDR_q)])
egger_complete <- nrow(egger_intercepts) == 30L && all(complete.cases(egger_intercepts))

mr_validation <- rbindlist(list(
  data.table(component="UVMR planned effect units", observed=nrow(uv), expected=240L,
             status=if(uv_complete) "PASS" else "FAIL", note="2 instrument sets x 2 exposures x 4 methods x 15 outcomes"),
  data.table(component="UVMR successful units", observed=sum(uv$convergence_status=="SUCCESS"), expected=240L,
             status=if(uv_success) "PASS" else "FAIL", note="No failed or silently omitted rows"),
  data.table(component="MVMR planned effect units", observed=nrow(mv), expected=150L,
             status=if(mv_complete) "PASS" else "FAIL", note="IVW, median, Q-minimisation, and two Egger orientations"),
  data.table(component="Unique Egger intercepts", observed=nrow(egger_intercepts), expected=30L,
             status=if(egger_complete) "PASS" else "FAIL", note="15 outcomes x 2 orientations"),
  data.table(component="MVMR Q outcomes", observed=nrow(het), expected=15L,
             status=if(nrow(het)==15L) "PASS" else "FAIL", note="Final covariance-aware IVW Q"),
  data.table(component="MVMR LOO outcomes", observed=nrow(loo), expected=15L,
             status=if(nrow(loo)==15L) "PASS" else "FAIL", note="Full per-SNP diagnostics retained separately")
))
fwrite(mr_validation, file.path(out,"FINAL_MR_RESULT_VALIDATION.tsv"), sep="\t")

# Build a per-SNP instrument manifest from the final MVMR and independent-UVMR sources.
indep <- fread(file.path(root,"final_major_revision/04_univariable_mr/INDEPENDENT_EXPOSURE_INSTRUMENTS.tsv"))
mv_in <- fread(file.path(root,"results/phase4_repair/PRIMARY_MVMR_CORRECTED_INPUT.tsv.gz"))
mv_p <- fread(file.path(root,"results/phase4/MVMR_EXPOSURE_MATRIX.tsv"), select=c("SNP","BMI_P","SMK_P"))
mv_in <- merge(mv_in,mv_p,by="SNP",all.x=TRUE)
mrahc_prov <- fread(file.path(root,"final_revision_closeout/01_mr_ahc_provenance/MR_AHC_INSTRUMENT_PROVENANCE.tsv"))

mv_manifest <- rbindlist(list(
  mv_in[, .(SNP,analysis_set="MVMR_COMPLETE_423",exposure="BMI",EA=BMI_EA,OA=BMI_OA,
            beta=BMI_beta,SE=BMI_SE,P=BMI_P,EAF=NA_real_,selection_rule="complete BMI-smoking MVMR matrix after six-exposure union clumping")],
  mv_in[, .(SNP,analysis_set="MVMR_COMPLETE_423",exposure="Smoking initiation",EA=SMK_EA,OA=SMK_OA,
            beta=SMK_beta,SE=SMK_SE,P=SMK_P,EAF=NA_real_,selection_rule="complete BMI-smoking MVMR matrix after six-exposure union clumping")]
))
ind_manifest <- indep[, .(SNP,analysis_set=ifelse(exposure=="BMI","INDEPENDENT_BMI_521","INDEPENDENT_SMK_92"),
                          exposure,EA,OA,beta,SE,P,EAF,
                          selection_rule="exposure-specific P<5e-8 and LD clumping r2<0.001 within 10 Mb")]
manifest <- rbindlist(list(mv_manifest,ind_manifest),fill=TRUE)
manifest <- merge(manifest,mrahc_prov[,.(rsID,MRAHC_historical_54=MRAHC_54,
                                       HISTORICAL_SMK_UVMR_16,INDEPENDENT_SMK_OUTCOME_USABLE_84)],
                  by.x="SNP",by.y="rsID",all.x=TRUE)
manifest[is.na(MRAHC_historical_54), `:=`(MRAHC_historical_54=FALSE,HISTORICAL_SMK_UVMR_16=FALSE,
                                         INDEPENDENT_SMK_OUTCOME_USABLE_84=FALSE)]
manifest[, `:=`(
  exposure_unit=fifelse(exposure=="BMI","per 1-SD inverse-normal BMI residual","per smoking-initiation liability unit"),
  outcome_unit="disease log-odds",
  harmonisation_note="Effect and other alleles shown for the source instrument row; outcome-specific harmonisation is in UVMR_HARMONISATION_AUDIT or the corrected MVMR matrix",
  missing_value_meaning="Blank means unavailable or not applicable; it does not imply zero",
  mr_ahc_status=fifelse(MRAHC_historical_54,"WITHDRAWN historical exploratory set","not applicable")
)]
fwrite(manifest,file.path(out,"FINAL_MR_INSTRUMENT_MANIFEST.tsv"),sep="\t")

dictionary <- data.table(
  instrument_set=c("MVMR union","Independent BMI UVMR","Independent smoking UVMR",
                   "Historical union-subset UVMR","Historical MR-AHC set","MVMR Q/LOO diagnostics","Per-SNP manifest"),
  analysis_label=c("MVMR_COMPLETE_423","INDEPENDENT_BMI_521","INDEPENDENT_SMK_92",
                   "UNION_SUBSET_UVMR","MRAHC_HISTORICAL_WITHDRAWN","MVMR_DIAGNOSTICS_FINAL","FINAL_MR_INSTRUMENT_MANIFEST"),
  count=c("423 unique; 422 for four outcome models","521 post-clump; 504-506 used per outcome",
          "92 post-clump; 84 outcome-usable","404 BMI; 16 smoking","54 historical variants; 27 smoking GWS",
          sprintf("15 Q rows; 15 outcome summaries plus %s per-SNP LOO rows", format(nrow(loo_per_snp), big.mark=",", scientific=FALSE)),"one row per instrument/exposure-set membership"),
  selection_or_definition=c("six-exposure union clump followed by complete BMI-smoking matrix","BMI GWAS P<5e-8; BMI-only clump",
                            "smoking GWAS P<5e-8; smoking-only clump","complete MVMR matrix then exposure-specific GWS",
                            "six-exposure clumped union, smoking F>10, complete four-outcome profiles",
                            "final harmonised 423/422-SNP model with SNP-specific covariance","combined MVMR and independent-UVMR instrument identities"),
  exposure_scale=c("BMI per 1-SD inverse-normal residual; smoking liability","BMI per 1-SD inverse-normal residual",
                   "smoking-initiation liability","same as parent exposure","profile denominator; not a common MR effect set",
                   "same as final MVMR","row-specific"),
  outcome_scale="disease log-odds",
  covariance_definition=c("SNP-specific 2x2 exposure-association sampling covariance from lambda and SEs","not applicable to UVMR",
                          "not applicable to UVMR","not applicable to UVMR","not used for current effect estimation",
                          "same SNP-specific 2x2 covariance for Q; IVW LOO point estimates do not depend on gencov in MVMR 0.4.9","see analysis_set"),
  status=c("CURRENT","CURRENT","CURRENT","HISTORICAL_SENSITIVITY","WITHDRAWN_WITH_JUSTIFICATION","CURRENT","CURRENT")
)
fwrite(dictionary,file.path(out,"Additional_file_6_MR_Instrument_Data_Dictionary.tsv"),sep="\t")

report <- c(
  "# Final machine-readable audit", "",
  sprintf("- Pair file: %d unique rows; all 14 locked counts %s.",nrow(pair),if(all(pair_checks$status=="PASS")) "PASS" else "FAIL"),
  sprintf("- significant_reversal_flag equals sign_reversal_flag AND difference_fdr < 0.05 for every row: %s.",flag_match),
  sprintf("- Resource-code mappings alc_problems/alcohol_problems, deaf/hearing_loss and schizo/SSD: %s.",mapping_ok),
  sprintf("- UVMR: %d/240 planned units present; %d successful.",nrow(uv),sum(uv$convergence_status=="SUCCESS")),
  sprintf("- MVMR: %d/150 planned effect units present; %d unique Egger intercepts.",nrow(mv),nrow(egger_intercepts)),
  sprintf("- Final Q and LOO diagnostics: %d and %d outcome rows; per-SNP LOO file retained separately.",nrow(het),nrow(loo)),
  "- The previous Additional file 6 was a five-row summary and was not a complete per-SNP audit. It has been replaced by an expanded dictionary plus FINAL_MR_INSTRUMENT_MANIFEST.tsv.",
  "- MR-AHC is labelled historical and withdrawn; no historical clustering output is presented as a current valid result."
)
writeLines(report,file.path(out,"FINAL_MACHINE_READABLE_AUDIT.md"))

stopifnot(all(pair_checks$status=="PASS"),flag_match,mapping_ok,uv_complete,uv_success,mv_complete,egger_complete)
cat("Machine-readable audit complete\n")
