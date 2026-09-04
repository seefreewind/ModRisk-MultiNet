#!/usr/bin/env python3
"""Manuscript number validation — all core numbers regenerated from frozen source tables."""
import csv, sys

FAILURES = []

def check(name, got, want):
    ok = got == want
    status = "PASS" if ok else "FAIL"
    print(f"[{status}] {name}: got={got} want={want}")
    if not ok: FAILURES.append(name)

# 1. joint rows
jr = list(csv.DictReader(open("results/phase3_joint/JOINT_PARTIAL_RG.tsv"), delimiter="\t"))
check("joint_pairs", len(jr), 2485)
check("baseline_FDR_significant", sum(1 for r in jr if float(r["baseline_FDR"]) < 0.05), 1535)
check("joint_residual_FDR_significant", sum(1 for r in jr if float(r["joint_FDR"]) < 0.05), 1100)
check("significant_joint_attenuation", sum(1 for r in jr if float(r["difference_FDR"]) < 0.05), 1019)
check("persistent", sum(1 for r in jr if r["machine_class"] == "PERSISTENT"), 417)
check("attenuated", sum(1 for r in jr if r["machine_class"] == "ATTENUATED"), 1051)
check("near_null_classified", sum(1 for r in jr if r["machine_class"] == "NEAR_NULL_EQUIVALENT"), 96)
check("masked_enhanced", sum(1 for r in jr if r["machine_class"] == "MASKED_ENHANCED"), 33)
check("unstable", sum(1 for r in jr if r["machine_class"] == "UNSTABLE"), 888)

# 2. single-factor rows
sf = list(csv.DictReader(open("results/phase2_single_factor/SINGLE_EXPOSURE_RESULTS.tsv"), delimiter="\t"))
check("single_exposure_partial_tests", len(sf), 14910)

# 3. disease-exposure
de = list(csv.DictReader(open("results/phase2_disease_exposure/DISEASE_EXPOSURE_RG.tsv"), delimiter="\t"))
check("disease_exposure_rg_tests", len(de), 426)
check("disease_exposure_FDR_significant", sum(1 for r in de if float(r["p_bh"]) < 0.05), 215)

# 4. MR
fdr = list(csv.DictReader(open("results/phase4/MODEL1_MVMR_IVW_FDR.tsv"), delimiter="\t"))
check("primary_MR_tests", len(fdr), 30)
check("SMK_FDR_significant", sum(1 for r in fdr if r["exposure"] == "SMK" and float(r["q"]) < 0.05), 15)
check("BMI_FDR_significant", sum(1 for r in fdr if r["exposure"] == "BMI" and float(r["q"]) < 0.05), 0)
ivw = list(csv.DictReader(open("results/phase4/MODEL1_MVMR_IVW.tsv"), delimiter="\t"))
check("MR_unique_outcomes", len(set(r["outcome"] for r in ivw)), 15)

# 5. pair classification
pc = list(csv.DictReader(open("results/phase4/PHASE4B_PAIR_CAUSAL_CLASSIFICATION.tsv"), delimiter="\t"))
check("pair_classifications", len(pc), 24)
check("SMK_concordant_pairs", sum(1 for r in pc if r["exposure"] == "SMK" and r["classification"] == "CONCORDANT_DUAL_OUTCOME_SUPPORT"), 12)
check("BMI_no_support_pairs", sum(1 for r in pc if r["exposure"] == "BMI" and r["classification"] == "NO_DUAL_OUTCOME_SUPPORT"), 12)

# 6. MR-AHC
ahc = list(csv.DictReader(open("results/phase5/MRAHC_CLUSTER_PROFILES.tsv"), delimiter="\t"))
check("MR_AHC_stable_clusters", sum(1 for r in ahc if r["status"] == "NO ROBUST MULTI-CLUSTER (published Q rule: all merged clusters fail)" or "NO ROBUST" in r.get("status", "")), 1)

# 7. orientation audit
print("\n" + ("ALL CHECKS PASSED" if not FAILURES else f"FAILURES: {FAILURES}"))
sys.exit(0 if not FAILURES else 1)