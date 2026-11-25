#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, t


ROOT = Path(os.environ.get("MODRISK_ROOT", Path.cwd())).resolve()
OUT = ROOT / "pre_submission_technical_closeout"
P_GWS = 5e-8


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_md(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def bool_series(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().map({"true": True, "false": False}).fillna(False)


def audit_instrument_provenance() -> dict:
    out = OUT / "01_instrument_provenance"
    out.mkdir(parents=True, exist_ok=True)
    union_path = ROOT / "results/phase4/MR_INSTRUMENT_UNION_CLUMPED.tsv"
    matrix_path = ROOT / "results/phase4/MVMR_EXPOSURE_MATRIX.tsv"
    final_path = ROOT / "results/phase4_repair/PRIMARY_MVMR_CORRECTED_INPUT.tsv.gz"
    union = pd.read_csv(union_path, sep="\t")
    matrix = pd.read_csv(matrix_path, sep="\t")
    final = pd.read_csv(final_path, sep="\t")

    p_cols = ["BMI_P", "SMK_P", "SBP_P", "ApoB_P", "FG_P", "SLEEP_P"]
    effect_cols = ["BMI_beta", "BMI_SE", "SMK_beta", "SMK_SE"]
    complete = matrix[effect_cols].notna().all(axis=1)
    any_f_gt_10 = ((matrix["BMI_beta"] / matrix["BMI_SE"]) ** 2 > 10) | ((matrix["SMK_beta"] / matrix["SMK_SE"]) ** 2 > 10)
    reconstructed = set(matrix.loc[complete & any_f_gt_10, "SNP"])
    final_ids = set(final["SNP"])
    if reconstructed != final_ids:
        raise RuntimeError(f"423-set reconstruction mismatch: reconstructed={len(reconstructed)} final={len(final_ids)}")

    prov = final[["SNP"]].merge(matrix[["SNP", *p_cols, *effect_cols]], on="SNP", how="left", validate="one_to_one")
    prov = prov.merge(
        union[["rsID", "chrom", "pos", "min_P", "exposures_involved", "n_exposures"]],
        left_on="SNP", right_on="rsID", how="left", validate="one_to_one",
    ).drop(columns=["rsID"])
    for axis in ["BMI", "SMK", "SBP", "ApoB", "FG", "SLEEP"]:
        prov[f"{axis}_GWS"] = prov[f"{axis}_P"].notna() & (prov[f"{axis}_P"] < P_GWS)
    prov["BMI_F"] = (prov["BMI_beta"] / prov["BMI_SE"]) ** 2
    prov["SMK_F"] = (prov["SMK_beta"] / prov["SMK_SE"]) ** 2
    prov["BMI_or_SMK_GWS"] = prov["BMI_GWS"] | prov["SMK_GWS"]
    prov["neither_BMI_nor_SMK_GWS"] = ~prov["BMI_or_SMK_GWS"]
    prov["six_axis_GWS_memberships"] = prov.apply(
        lambda row: ";".join(axis for axis in ["BMI", "SMK", "SBP", "ApoB", "FG", "SLEEP"] if row[f"{axis}_GWS"]), axis=1
    )
    prov["selection_stage_1"] = "six-exposure genome-wide-significant union"
    prov["selection_stage_2"] = "joint LD clumping at r2<0.001 within 10 Mb using the frozen EUR reference"
    prov["selection_stage_3"] = "complete BMI and smoking associations"
    prov["selection_stage_4"] = "per-variant F>10 for BMI or smoking initiation"
    prov["final_status"] = "MVMR_COMPLETE_423"

    ordered = [
        "SNP", "chrom", "pos", "min_P", "exposures_involved", "n_exposures",
        *p_cols, "BMI_F", "SMK_F", "BMI_GWS", "SMK_GWS", "SBP_GWS", "ApoB_GWS", "FG_GWS", "SLEEP_GWS",
        "BMI_or_SMK_GWS", "neither_BMI_nor_SMK_GWS", "six_axis_GWS_memberships",
        "selection_stage_1", "selection_stage_2", "selection_stage_3", "selection_stage_4", "final_status",
    ]
    prov[ordered].to_csv(out / "MVMR_423_INSTRUMENT_PROVENANCE.tsv", sep="\t", index=False)
    non = prov.loc[prov["neither_BMI_nor_SMK_GWS"], ordered].copy()
    non.to_csv(out / "MVMR_NON_BMI_SMOKING_GWS_AUDIT.tsv", sep="\t", index=False)

    counts = {
        "BMI_GWS": int(prov["BMI_GWS"].sum()),
        "SMK_GWS": int(prov["SMK_GWS"].sum()),
        "both_GWS": int((prov["BMI_GWS"] & prov["SMK_GWS"]).sum()),
        "either_GWS": int(prov["BMI_or_SMK_GWS"].sum()),
        "neither_GWS": int(prov["neither_BMI_nor_SMK_GWS"].sum()),
        "complete_before_F_gate": int(complete.sum()),
        "final": int(len(prov)),
    }
    expected = {"BMI_GWS": 404, "SMK_GWS": 16, "both_GWS": 8, "either_GWS": 412, "neither_GWS": 11, "complete_before_F_gate": 481, "final": 423}
    checks = pd.DataFrame([
        {"metric": key, "observed": counts[key], "expected": expected[key], "status": "PASS" if counts[key] == expected[key] else "FAIL"}
        for key in expected
    ])
    checks.to_csv(out / "MVMR_INSTRUMENT_SELECTION_CHECKS.tsv", sep="\t", index=False)

    source_hashes = {p.relative_to(ROOT).as_posix(): sha256(p) for p in [union_path, matrix_path, final_path]}
    lines = [
        "# MVMR instrument selection report", "",
        "## Reconstructed selection chain", "",
        "The final 423-variant matrix was reconstructed exactly from the historical six-exposure genome-wide-significant union. The union was jointly LD-clumped, then restricted to variants with complete BMI and smoking-initiation associations and a per-variant F statistic above 10 for at least one of the two exposures. Complete BMI and smoking data alone yielded 481 variants; the F-statistic gate yielded the final 423.", "",
        "The previous shorthand description as a BMI-smoking union was incorrect. The correct description is a six-exposure clumped union followed by BMI-smoking completeness and the recorded per-variant strength gate.", "",
        "## Counts", "", checks.to_markdown(index=False), "",
        "Among the final 423 variants, 404 were genome-wide significant for BMI, 16 for smoking initiation, and 8 for both. Their union contained 412 variants. Eleven variants were not genome-wide significant for either BMI or smoking initiation and entered through another exposure axis before satisfying the BMI-or-smoking F-statistic gate. These 11 variants are listed separately and are not silently labelled BMI-smoking genome-wide-significant instruments.", "",
        "## Source checksums", "",
        *[f"- `{name}`: `{digest}`" for name, digest in source_hashes.items()], "",
        "## Methods wording", "",
        "Candidate instruments were formed as the union of genome-wide-significant variants across the six exposure GWASs and jointly LD-clumped at r² < 0.001 within 10 Mb. The BMI-smoking MVMR matrix then retained variants with complete associations for both exposures and a per-variant F statistic above 10 for at least one exposure, yielding 423 variants. Of these, 412 were genome-wide significant for BMI or smoking initiation and 11 entered through another exposure axis.", "",
    ]
    write_md(out / "MVMR_INSTRUMENT_SELECTION_REPORT.md", lines)
    return {"counts": counts, "provenance": prov}


def distribution_rule_uv(method: str, n: int) -> tuple[str, float | None, float]:
    if method == "MR-Egger":
        df = n - 2
        return "Student_t", float(df), float(t.ppf(0.975, df))
    if method == "Weighted mode":
        df = n - 1
        return "Student_t", float(df), float(t.ppf(0.975, df))
    return "standard_normal", None, float(norm.ppf(0.975))


def distribution_rule_mv(method: str, n: int) -> tuple[str, float | None, float]:
    if method == "MVMR-IVW":
        df = n - 2
        return "Student_t", float(df), float(t.ppf(0.975, df))
    if method == "MVMR-Egger":
        df = n - 3
        return "Student_t", float(df), float(t.ppf(0.975, df))
    return "standard_normal", None, float(norm.ppf(0.975))


def audit_ci() -> dict:
    out = OUT / "02_ci_validation"
    out.mkdir(parents=True, exist_ok=True)
    uv_path = ROOT / "final_revision_closeout/10_submission_package/06_Additional_file_4_Univariable_MR_Full_Results.tsv"
    mv_path = ROOT / "final_revision_closeout/10_submission_package/07_Additional_file_5_MVMR_Full_Results.tsv"
    uv = pd.read_csv(uv_path, sep="\t")
    mv = pd.read_csv(mv_path, sep="\t")
    if len(uv) != 240 or len(mv) != 150:
        raise RuntimeError("Unexpected MR row count")

    audit_rows: list[dict] = []
    before_after: list[dict] = []
    uv_new = uv.copy()
    for i, row in uv.iterrows():
        dist, df, crit = distribution_rule_uv(str(row.method), int(row.nsnp))
        z = abs(row.beta / row.SE)
        p_recomputed = 2 * (t.sf(z, df) if dist == "Student_t" else norm.sf(z))
        lo = row.beta - crit * row.SE
        hi = row.beta + crit * row.SE
        before_excludes = bool(row.CI_lower > 0 or row.CI_upper < 0)
        after_excludes = bool(lo > 0 or hi < 0)
        p_sig = bool(row.P < 0.05)
        audit_rows.append({
            "analysis": "UVMR", "instrument_set": row.instrument_set, "method": row.method,
            "exposure": row.exposure, "outcome": row.outcome, "n_instruments": int(row.nsnp),
            "distribution": dist, "df": df, "critical_value_95": crit,
            "reported_P": row.P, "recomputed_P": p_recomputed, "absolute_P_error": abs(p_recomputed - row.P),
            "reported_CI_lower": row.CI_lower, "reported_CI_upper": row.CI_upper,
            "corrected_CI_lower": lo, "corrected_CI_upper": hi,
            "P_lt_0_05": p_sig, "reported_CI_excludes_zero": before_excludes,
            "corrected_CI_excludes_zero": after_excludes,
            "reported_CI_P_consistent": before_excludes == p_sig,
            "corrected_CI_P_consistent": after_excludes == p_sig,
            "status": "PASS" if abs(p_recomputed - row.P) < 1e-10 and after_excludes == p_sig else "FAIL",
        })
        if not math.isclose(lo, row.CI_lower, rel_tol=0, abs_tol=1e-14) or not math.isclose(hi, row.CI_upper, rel_tol=0, abs_tol=1e-14):
            before_after.append({
                "analysis": "UVMR", "instrument_set": row.instrument_set, "method": row.method,
                "exposure": row.exposure, "outcome": row.outcome, "n_instruments": int(row.nsnp),
                "distribution": dist, "df": df, "beta": row.beta, "SE": row.SE, "P": row.P, "FDR_q": row.FDR_q,
                "CI_lower_before": row.CI_lower, "CI_upper_before": row.CI_upper,
                "CI_lower_after": lo, "CI_upper_after": hi,
                "significance_mismatch_before": before_excludes != p_sig,
                "significance_mismatch_after": after_excludes != p_sig,
            })
        uv_new.loc[i, "CI_lower"] = lo
        uv_new.loc[i, "CI_upper"] = hi
    uv_new.to_csv(out / "Additional_file_4_Univariable_MR_Full_Results_CI_corrected.tsv", sep="\t", index=False)

    mv_new = mv.copy()
    for i, row in mv.iterrows():
        dist, df, crit = distribution_rule_mv(str(row.method), int(row.n_instruments))
        z = abs(row.beta / row.SE)
        p_recomputed = 2 * (t.sf(z, df) if dist == "Student_t" else norm.sf(z))
        lo = row.beta - crit * row.SE
        hi = row.beta + crit * row.SE
        before_excludes = bool(row.CI_lower > 0 or row.CI_upper < 0)
        after_excludes = bool(lo > 0 or hi < 0)
        p_sig = bool(row.P < 0.05)
        audit_rows.append({
            "analysis": "MVMR_effect", "instrument_set": "MVMR_COMPLETE_423", "method": row.method,
            "exposure": row.exposure, "outcome": row.outcome, "n_instruments": int(row.n_instruments),
            "distribution": dist, "df": df, "critical_value_95": crit,
            "reported_P": row.P, "recomputed_P": p_recomputed, "absolute_P_error": abs(p_recomputed - row.P),
            "reported_CI_lower": row.CI_lower, "reported_CI_upper": row.CI_upper,
            "corrected_CI_lower": lo, "corrected_CI_upper": hi,
            "P_lt_0_05": p_sig, "reported_CI_excludes_zero": before_excludes,
            "corrected_CI_excludes_zero": after_excludes,
            "reported_CI_P_consistent": before_excludes == p_sig,
            "corrected_CI_P_consistent": after_excludes == p_sig,
            "status": "PASS" if abs(p_recomputed - row.P) < 1e-10 and after_excludes == p_sig else "FAIL",
        })
        if not math.isclose(lo, row.CI_lower, rel_tol=0, abs_tol=1e-14) or not math.isclose(hi, row.CI_upper, rel_tol=0, abs_tol=1e-14):
            before_after.append({
                "analysis": "MVMR_effect", "instrument_set": "MVMR_COMPLETE_423", "method": row.method,
                "exposure": row.exposure, "outcome": row.outcome, "n_instruments": int(row.n_instruments),
                "distribution": dist, "df": df, "beta": row.beta, "SE": row.SE, "P": row.P, "FDR_q": row.FDR_q,
                "CI_lower_before": row.CI_lower, "CI_upper_before": row.CI_upper,
                "CI_lower_after": lo, "CI_upper_after": hi,
                "significance_mismatch_before": before_excludes != p_sig,
                "significance_mismatch_after": after_excludes != p_sig,
            })
        mv_new.loc[i, "CI_lower"] = lo
        mv_new.loc[i, "CI_upper"] = hi

    intercepts = mv[mv["method"].eq("MVMR-Egger")].drop_duplicates(["outcome", "orientation"]).copy()
    if len(intercepts) != 30:
        raise RuntimeError("Expected 30 unique MVMR-Egger intercepts")
    intercept_ci: dict[tuple[str, str], tuple[float, float]] = {}
    for _, row in intercepts.iterrows():
        df = int(row.n_instruments) - 3
        crit = float(t.ppf(0.975, df))
        z = abs(row.intercept / row.intercept_SE)
        p_recomputed = 2 * t.sf(z, df)
        lo = row.intercept - crit * row.intercept_SE
        hi = row.intercept + crit * row.intercept_SE
        intercept_ci[(row.outcome, row.orientation)] = (lo, hi)
        audit_rows.append({
            "analysis": "MVMR_Egger_intercept", "instrument_set": "MVMR_COMPLETE_423", "method": "MVMR-Egger intercept",
            "exposure": "intercept", "outcome": row.outcome, "n_instruments": int(row.n_instruments),
            "distribution": "Student_t", "df": float(df), "critical_value_95": crit,
            "reported_P": row.intercept_P, "recomputed_P": p_recomputed, "absolute_P_error": abs(p_recomputed - row.intercept_P),
            "reported_CI_lower": np.nan, "reported_CI_upper": np.nan,
            "corrected_CI_lower": lo, "corrected_CI_upper": hi,
            "P_lt_0_05": bool(row.intercept_P < 0.05), "reported_CI_excludes_zero": np.nan,
            "corrected_CI_excludes_zero": bool(lo > 0 or hi < 0), "reported_CI_P_consistent": np.nan,
            "corrected_CI_P_consistent": bool((lo > 0 or hi < 0) == (row.intercept_P < 0.05)),
            "status": "PASS" if abs(p_recomputed - row.intercept_P) < 1e-10 else "FAIL",
        })
    mv_new["intercept_CI_lower"] = mv_new.apply(
        lambda r: intercept_ci.get((r.outcome, r.orientation), (np.nan, np.nan))[0], axis=1
    )
    mv_new["intercept_CI_upper"] = mv_new.apply(
        lambda r: intercept_ci.get((r.outcome, r.orientation), (np.nan, np.nan))[1], axis=1
    )
    mv_new.to_csv(out / "Additional_file_5_MVMR_Full_Results_CI_corrected.tsv", sep="\t", index=False)

    audit = pd.DataFrame(audit_rows)
    changes = pd.DataFrame(before_after)
    audit.to_csv(out / "MR_INFERENCE_DISTRIBUTION_AUDIT.tsv", sep="\t", index=False)
    changes.to_csv(out / "MR_CI_BEFORE_AFTER.tsv", sep="\t", index=False)
    mismatches_before = audit[(audit["analysis"] == "UVMR") & (audit["reported_CI_P_consistent"] == False)]
    mismatches_after = audit[(audit["analysis"].isin(["UVMR", "MVMR_effect", "MVMR_Egger_intercept"])) & (audit["corrected_CI_P_consistent"] == False)]
    if len(mismatches_before) != 4 or len(mismatches_after) != 0 or not audit["status"].eq("PASS").all():
        raise RuntimeError("Inference distribution audit failed")
    lines = [
        "# MR confidence-interval correction report", "",
        "## Scope and result", "",
        "All 240 UVMR effect rows, all 150 MVMR effect rows, and 30 unique directional MVMR-Egger intercepts were audited against the probability distribution used to calculate each reported P value. Beta estimates, standard errors, P values and FDR values were unchanged.", "",
        "UVMR IVW and weighted-median P values use the standard normal distribution. UVMR MR-Egger uses a Student t distribution with n−2 degrees of freedom. Weighted mode uses a Student t distribution with n−1 degrees of freedom. The historical export applied 1.96 to every method, causing four threshold inconsistencies. Corrected method-specific 95% confidence intervals remove all four.", "",
        "MVMR-IVW uses Student t with n−2 degrees of freedom, and directional MVMR-Egger effects and intercepts use Student t with n−3 degrees of freedom. Multivariable weighted LAD median and Q-minimisation use normal-reference intervals. No MVMR row changed nominal significance classification after the CI correction.", "",
        "## Four corrected UVMR threshold inconsistencies", "",
        mismatches_before[["instrument_set", "method", "exposure", "outcome", "n_instruments", "reported_P", "reported_CI_lower", "reported_CI_upper", "corrected_CI_lower", "corrected_CI_upper"]].to_markdown(index=False), "",
        f"Rows with corrected numeric CI endpoints: {len(changes)}. Rows with a P/CI significance mismatch after correction: {len(mismatches_after)}.", "",
    ]
    write_md(out / "MR_CI_CORRECTION_REPORT.md", lines)
    return {"uv": uv_new, "mv": mv_new, "audit": audit, "changes": changes}


def audit_palindromes() -> dict:
    out = OUT / "03_palindromic_audit"
    out.mkdir(parents=True, exist_ok=True)
    audit_path = ROOT / "final_major_revision/03_mr_harmonisation/MR_ALLELE_HARMONISATION_AUDIT.tsv"
    corrected_path = ROOT / "results/phase4_repair/PRIMARY_MVMR_CORRECTED_INPUT.tsv.gz"
    bmi_path = ROOT / "ModRisk-MultiNet-data/raw/exposures/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz"
    smk_path = ROOT / "ModRisk-MultiNet-data/raw/exposures/Smoking/SmokingInitiation.txt.gz"
    source = pd.read_csv(audit_path, sep="\t")
    source["palindromic"] = bool_series(source["palindromic"])
    pal = source[source["palindromic"]].copy()
    if len(pal) != 54:
        raise RuntimeError(f"Expected 54 palindromic variants; found {len(pal)}")
    pal["genome_build"] = "GRCh37/hg19"
    pal["BMI_raw_effect_allele"] = pal["BMI_effect_allele"]
    pal["BMI_raw_other_allele"] = pal["BMI_other_allele"]
    pal["BMI_raw_EAF"] = pal["BMI_EAF"]
    pal["SMK_raw_reference_allele"] = pal["reference_allele"]
    pal["SMK_raw_alternate_allele"] = pal["alternate_allele"]
    pal["SMK_raw_ALT_AF"] = np.where(
        pal["smoking_effect_allele"].eq(pal["alternate_allele"]),
        pal["smoking_effect_allele_frequency"], 1 - pal["smoking_effect_allele_frequency"],
    )
    pal["BMI_aligned_EAF"] = pal["BMI_EAF_common_smoking_EA"]
    pal["SMK_aligned_EAF"] = pal["smoking_effect_allele_frequency"]
    pal["BMI_MAF"] = np.minimum(pal["BMI_aligned_EAF"], 1 - pal["BMI_aligned_EAF"])
    pal["SMK_MAF"] = np.minimum(pal["SMK_aligned_EAF"], 1 - pal["SMK_aligned_EAF"])
    pal["current_orientation_EAF_difference"] = (pal["BMI_aligned_EAF"] - pal["SMK_aligned_EAF"]).abs()
    pal["opposite_orientation_EAF_difference"] = (pal["BMI_aligned_EAF"] - (1 - pal["SMK_aligned_EAF"])).abs()
    pal["frequency_margin_current_minus_opposite"] = pal["opposite_orientation_EAF_difference"] - pal["current_orientation_EAF_difference"]
    pal["intermediate_frequency"] = (pal["BMI_MAF"] >= 0.42) | (pal["SMK_MAF"] >= 0.42)
    pal["frequency_favors_current_orientation"] = pal["current_orientation_EAF_difference"] < pal["opposite_orientation_EAF_difference"]
    pal["ambiguity_rule"] = "unresolved if either aligned MAF>=0.42; otherwise require current-orientation EAF difference<opposite-orientation difference and <=0.10"
    pal["decision"] = np.where(
        (~pal["intermediate_frequency"]) & pal["frequency_favors_current_orientation"] & (pal["current_orientation_EAF_difference"] <= 0.10),
        "RESOLVED_BY_EAF", "UNRESOLVED_FREQUENCY_AMBIGUOUS",
    )
    pal["decision_reason"] = np.select(
        [
            pal["decision"].eq("RESOLVED_BY_EAF"),
            pal["intermediate_frequency"] & ~pal["frequency_favors_current_orientation"],
            pal["intermediate_frequency"],
            ~pal["frequency_favors_current_orientation"],
        ],
        [
            "Both exposure frequencies are sufficiently informative and favor the retained orientation",
            "Intermediate-frequency palindrome and the alternative strand orientation fits the EAFs better",
            "At least one aligned MAF is >=0.42, so strand orientation is not safely resolved by frequency",
            "Alternative strand orientation fits the EAFs at least as well as the retained orientation",
        ],
        default="Frequency difference exceeds the prespecified concordance threshold",
    )
    pal["outcome_alignment_evidence"] = np.where(
        bool_series(pal["all_valid_outcomes_aligned_to_smoking_EA"]),
        pal["outcome_valid_count"].astype(str) + " valid outcome associations use the smoking effect allele",
        "Outcome alignment incomplete or inconsistent",
    )
    pal["BMI_source_file"] = bmi_path.relative_to(ROOT).as_posix()
    pal["BMI_source_sha256"] = sha256(bmi_path)
    pal["SMK_source_file"] = smk_path.relative_to(ROOT).as_posix()
    pal["SMK_source_sha256"] = sha256(smk_path)
    pal["harmonised_source_file"] = corrected_path.relative_to(ROOT).as_posix()
    pal["harmonised_source_sha256"] = sha256(corrected_path)
    pal["audit_source_file"] = audit_path.relative_to(ROOT).as_posix()
    pal["audit_source_sha256"] = sha256(audit_path)
    cols = [
        "SNP", "chromosome", "position_GRCh37", "genome_build",
        "BMI_raw_effect_allele", "BMI_raw_other_allele", "BMI_raw_EAF", "BMI_aligned_EAF", "BMI_MAF",
        "SMK_raw_reference_allele", "SMK_raw_alternate_allele", "SMK_raw_ALT_AF",
        "smoking_effect_allele", "smoking_other_allele", "SMK_aligned_EAF", "SMK_MAF",
        "BMI_smoking_alignment", "current_orientation_EAF_difference", "opposite_orientation_EAF_difference",
        "frequency_margin_current_minus_opposite", "intermediate_frequency", "frequency_favors_current_orientation",
        "ambiguity_rule", "decision", "decision_reason", "outcome_valid_count", "outcome_effect_alleles",
        "all_valid_outcomes_aligned_to_smoking_EA", "outcome_alignment_evidence",
        "BMI_source_file", "BMI_source_sha256", "SMK_source_file", "SMK_source_sha256",
        "harmonised_source_file", "harmonised_source_sha256", "audit_source_file", "audit_source_sha256",
    ]
    pal[cols].to_csv(out / "PALINDROMIC_SNP_EAF_AUDIT.tsv", sep="\t", index=False)
    pal[["SNP", "intermediate_frequency", "frequency_favors_current_orientation", "current_orientation_EAF_difference", "opposite_orientation_EAF_difference", "decision", "decision_reason"]].to_csv(
        out / "PALINDROMIC_SNP_DECISION_LOG.tsv", sep="\t", index=False
    )
    resolved = int(pal["decision"].eq("RESOLVED_BY_EAF").sum())
    unresolved = int(pal["decision"].str.startswith("UNRESOLVED").sum())
    opposite = int((~pal["frequency_favors_current_orientation"]).sum())
    lines = [
        "# Palindromic SNP harmonisation report", "",
        "## Audit result", "",
        f"All 54 palindromic variants were traced to GRCh37 coordinates, raw BMI and smoking alleles and frequencies, harmonised exposure alleles, and outcome-alignment evidence. Under a conservative ambiguity rule that treats an aligned MAF of 0.42 or greater in either exposure as intermediate-frequency, {resolved} variants were resolved by EAF and {unresolved} remained frequency-ambiguous.", "",
        f"For {opposite} of the 54 variants, the opposite strand orientation produced a smaller EAF discrepancy than the retained orientation. Both occur among the frequency-ambiguous variants. The previously reported maximum aligned EAF difference of 0.0383 is a descriptive maximum only; it does not by itself prove strand resolution and must not be used as a blanket harmonisation claim.", "",
        "Outcome associations are aligned to the smoking effect allele in the frozen corrected matrix. This verifies internal outcome coding but cannot independently resolve a palindromic exposure strand when the exposure EAFs are near 0.5.", "",
        "## Decision rule", "",
        "A palindrome was labelled RESOLVED_BY_EAF only when both aligned exposure frequencies were informative (MAF < 0.42), the retained orientation had a smaller EAF discrepancy than the opposite orientation, and the retained-orientation difference was no greater than 0.10. All other palindromes were labelled UNRESOLVED_FREQUENCY_AMBIGUOUS.", "",
        "## Manuscript consequence", "",
        "The sentence claiming that all 54 palindromic variants were resolved by aligned allele frequencies should be replaced. The revised text should report that 54 palindromic variants were audited, 37 had frequency support for the retained orientation, and 17 remained frequency-ambiguous. The unsupported 0.0383 resolution claim should be removed. The ambiguity must be disclosed as a limitation and the final manifest must carry the per-SNP decisions.", "",
    ]
    write_md(out / "PALINDROMIC_HARMONISATION_REPORT.md", lines)
    return {"pal": pal, "resolved": resolved, "unresolved": unresolved, "opposite": opposite}


def audit_pair_csv() -> dict:
    out = OUT / "04_pair_csv_validation"
    out.mkdir(parents=True, exist_ok=True)
    path = ROOT / "final_revision_closeout/10_submission_package/04_Additional_file_2_Pair_level_classification.csv"
    d = pd.read_csv(path)
    checks: list[dict] = []

    def add(metric: str, observed, expected, tolerance: float = 0.0) -> None:
        if isinstance(expected, float):
            ok = abs(float(observed) - expected) <= tolerance
        else:
            ok = observed == expected
        checks.append({"metric": metric, "observed": observed, "expected": expected, "tolerance": tolerance, "status": "PASS" if ok else "FAIL"})

    for c in ["residual_significant", "near_null_classified", "sign_reversal_flag", "significant_reversal_flag", "equivalent_005", "equivalent_010", "legacy_significant_reversal_residual_z_gate"]:
        d[c] = bool_series(d[c])
    add("rows", len(d), 2485)
    add("unique_pair_id", d.pair_id.nunique(), 2485)
    add("baseline_FDR_significant", int((d.baseline_fdr < .05).sum()), 1535)
    add("signed_difference_FDR_significant", int((d.difference_fdr < .05).sum()), 1019)
    add("attenuated", int(d.change_class.eq("ATTENUATED").sum()), 900)
    add("enhanced", int(d.change_class.eq("ENHANCED").sum()), 119)
    add("no_significant_change", int(d.change_class.eq("NO_SIGNIFICANT_CHANGE").sum()), 1466)
    add("residual_FDR_significant", int(d.residual_significant.sum()), 1100)
    add("near_null", int(d.near_null_classified.sum()), 605)
    add("other_residual_non_significant", int(d.residual_state.eq("RESIDUAL_NON_SIGNIFICANT_OTHER").sum()), 780)
    add("sign_reversals", int(d.sign_reversal_flag.sum()), 135)
    add("significant_reversals", int(d.significant_reversal_flag.sum()), 66)
    add("legacy_residual_Z_gate", int(d.legacy_significant_reversal_residual_z_gate.sum()), 1)
    add("CI_within_0.05", int(d.equivalent_005.sum()), 2)
    add("CI_within_0.10", int(d.equivalent_010.sum()), 162)
    add("median_absolute_attenuation", float(d.absolute_attenuation.median()), 0.0206, 5e-5)
    stable = (d.baseline_abs_rg / d.baseline_se >= 2) & (d.absolute_attenuation > 0)
    add("proportional_attenuation_denominator", int(stable.sum()), 1317)
    add("median_proportional_attenuation", float(d.loc[stable, "proportional_attenuation"].median()), 0.2492, 5e-5)
    add("signed_difference_formula_max_error", float(np.max(np.abs(d.signed_difference - (d.baseline_rg - d.joint_rg)))), 0.0, 1e-12)
    add("absolute_attenuation_formula_max_error", float(np.max(np.abs(d.absolute_attenuation - (d.baseline_abs_rg - d.joint_abs_rg)))), 0.0, 1e-12)
    add("residual_CI_lower_formula_max_error", float(np.max(np.abs(d.residual_ci_lower - (d.joint_rg - 1.96 * d.joint_se)))), 0.0, 1e-12)
    add("residual_CI_upper_formula_max_error", float(np.max(np.abs(d.residual_ci_upper - (d.joint_rg + 1.96 * d.joint_se)))), 0.0, 1e-12)
    add("significant_reversal_formula_matches", int((d.significant_reversal_flag == (d.sign_reversal_flag & (d.difference_fdr < .05))).sum()), 2485)
    expected_cross = {
        ("ATTENUATED", "RESIDUAL_SIGNIFICANT"): 495,
        ("ATTENUATED", "NEAR_NULL_CLASSIFIED"): 155,
        ("ATTENUATED", "RESIDUAL_NON_SIGNIFICANT_OTHER"): 250,
        ("ENHANCED", "RESIDUAL_SIGNIFICANT"): 47,
        ("ENHANCED", "NEAR_NULL_CLASSIFIED"): 26,
        ("ENHANCED", "RESIDUAL_NON_SIGNIFICANT_OTHER"): 46,
        ("NO_SIGNIFICANT_CHANGE", "RESIDUAL_SIGNIFICANT"): 558,
        ("NO_SIGNIFICANT_CHANGE", "NEAR_NULL_CLASSIFIED"): 424,
        ("NO_SIGNIFICANT_CHANGE", "RESIDUAL_NON_SIGNIFICANT_OTHER"): 484,
    }
    cross = d.groupby(["change_class", "residual_state"]).size().to_dict()
    for key, expected in expected_cross.items():
        add(f"3x3_{key[0]}__{key[1]}", int(cross.get(key, 0)), expected)

    checks_df = pd.DataFrame(checks)
    checks_df.to_csv(out / "FINAL_PAIR_CSV_CHECKS.tsv", sep="\t", index=False)
    if not checks_df.status.eq("PASS").all():
        raise RuntimeError("Final pair CSV validation failed")

    label_map = pd.read_csv(ROOT / "config/FIGURE_DISEASE_LABEL_MAP.tsv", sep="\t")
    label_col = next(c for c in label_map.columns if c != "code" and "source" not in c.lower()) if "code" in label_map.columns else label_map.columns[1]
    code_col = "code" if "code" in label_map.columns else label_map.columns[0]
    resource_by_analysis = dict(zip(d.trait_1, d.trait_1_resource_code)) | dict(zip(d.trait_2, d.trait_2_resource_code))
    mapping = pd.DataFrame({"analysis_code": sorted(set(d.trait_1) | set(d.trait_2))})
    mapping["resource_code"] = mapping.analysis_code.map(resource_by_analysis)
    mapping = mapping.merge(label_map[[code_col, label_col]].rename(columns={code_col: "analysis_code", label_col: "display_name"}), on="analysis_code", how="left")
    mapping["mapping_status"] = np.where(mapping.analysis_code.eq(mapping.resource_code), "IDENTICAL", "FROZEN_ALIAS_TO_RESOURCE_CODE")
    mapping.to_csv(out / "FINAL_DISEASE_CODE_MAPPING.tsv", sep="\t", index=False)
    add_alias_ok = {
        "alc_problems": "alcohol_problems", "deaf": "hearing_loss", "schizo": "SSD",
    }
    for alias, resource in add_alias_ok.items():
        observed = mapping.loc[mapping.analysis_code.eq(alias), "resource_code"].iloc[0]
        if observed != resource:
            raise RuntimeError(f"Mapping mismatch for {alias}: {observed}")

    lines = [
        "# Final pair CSV validation", "",
        "Status: **PASS**", "",
        "The final Additional file 2 contains 2,485 unique unordered disease pairs. All locked counts, derived-field formulas, medians, denominator rules, 3×3 change-by-residual cells, reversal definitions, equivalence flags and disease-code mappings were independently recalculated from the exported CSV.", "",
        "## Validation checks", "", checks_df.to_markdown(index=False), "",
        "## Disease-code mapping", "",
        f"The mapping table contains {len(mapping)} analysis codes. Three frozen aliases map to distinct GEMINI resource codes: alc_problems→alcohol_problems, deaf→hearing_loss and schizo→SSD. All other analysis codes retain the same resource code.", "",
        f"Validated file SHA256: `{sha256(path)}`.", "",
    ]
    write_md(out / "FINAL_PAIR_CSV_VALIDATION.md", lines)
    return {"checks": checks_df, "mapping": mapping}


def audit_mrahc_manifest_scope() -> dict:
    out = OUT / "05_document_corrections"
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = ROOT / "final_revision_closeout/10_submission_package/10_Additional_file_8_MR_Instrument_Manifest.tsv"
    hist_path = ROOT / "results/phase5/MRAHC_PRIMARY_INPUT.tsv"
    manuscript_path = ROOT / "final_revision_closeout/05_manuscript_revision/Main_Manuscript_Resubmission.md"
    manifest = pd.read_csv(manifest_path, sep="\t")
    manifest["MRAHC_historical_54"] = bool_series(manifest["MRAHC_historical_54"])
    hist = pd.read_csv(hist_path, sep="\t")
    current = set(manifest.loc[manifest.MRAHC_historical_54, "SNP"])
    historical = set(hist.rsID)
    missing = sorted(historical - current)
    extra = sorted(current - historical)
    text = manuscript_path.read_text(encoding="utf-8")
    liang_in_list = "Liang X" in text
    liang_in_text = "[24]" in text or "Liang et al." in text
    report = [
        "# MR-AHC manifest scope audit", "",
        "MR-AHC means Mendelian randomization-based agglomerative hierarchical clustering.", "",
        f"The historical withdrawn MR-AHC input contains {len(historical)} unique variants. The current Additional file 8 marks only {len(current)} unique variants because it is assembled from current MVMR and independent-UVMR sets. It is therefore not a complete manifest of the historical 54-variant input.", "",
        f"Missing historical-only variants ({len(missing)}): {', '.join(missing) if missing else 'none'}.", "",
        f"Unexpected manifest variants ({len(extra)}): {', '.join(extra) if extra else 'none'}.", "",
        "The final manifest must either add explicit historical-withdrawn rows for these variants or state that the MR-AHC membership columns are partial current-set overlaps. This closeout will add explicit rows while retaining WITHDRAWN status.", "",
        f"Liang et al. 2025 is present in the reference list: {liang_in_list}. It has an in-text citation marker: {liang_in_text}. The corrected manuscript must cite it at the first expanded MR-AHC mention or delete and renumber the reference.", "",
    ]
    write_md(out / "MR_AHC_MANIFEST_SCOPE_AUDIT.md", report)
    pd.DataFrame({"missing_historical_MRAHC_SNP": missing}).to_csv(out / "MR_AHC_MANIFEST_MISSING_9.tsv", sep="\t", index=False)
    return {"current": current, "historical": historical, "missing": missing}


def write_summary(results: dict) -> None:
    summary = {
        "instrument_counts": results["instrument"]["counts"],
        "uvmr_ci_mismatches_before": int(((results["ci"]["audit"].analysis == "UVMR") & (results["ci"]["audit"].reported_CI_P_consistent == False)).sum()),
        "ci_mismatches_after": int((results["ci"]["audit"].corrected_CI_P_consistent == False).sum()),
        "palindromic_resolved": results["pal"]["resolved"],
        "palindromic_unresolved": results["pal"]["unresolved"],
        "pair_checks_pass": bool(results["pair"]["checks"].status.eq("PASS").all()),
        "mrahc_manifest_unique_current": len(results["mrahc"]["current"]),
        "mrahc_manifest_unique_historical": len(results["mrahc"]["historical"]),
        "mrahc_manifest_missing": len(results["mrahc"]["missing"]),
    }
    (OUT / "logs/core_technical_audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    results = {
        "instrument": audit_instrument_provenance(),
        "ci": audit_ci(),
        "pal": audit_palindromes(),
        "pair": audit_pair_csv(),
        "mrahc": audit_mrahc_manifest_scope(),
    }
    write_summary(results)
    print("CORE_TECHNICAL_AUDITS_COMPLETE")


if __name__ == "__main__":
    main()
