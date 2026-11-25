#!/usr/bin/env python3
"""Build strict-primary MR tables and machine-readable submission attachments."""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

from scipy.stats import t as student_t

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "final_manuscript_error_correction"
MR = WORK / "02_mr_sensitivity"
PAL = WORK / "01_palindromic_validation"
TABLES = WORK / "07_supplement/tables"
PKG = WORK / "10_submission_package"


def read_tsv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def bh(values):
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    out = [math.nan] * n
    running = 1.0
    for rank_index in range(n - 1, -1, -1):
        i = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[i] * n / rank)
        out[i] = min(1.0, running)
    return out


def num(x):
    return float(x) if x not in (None, "", "NA", "NaN") else math.nan


def fmt(x, digits=3):
    x = num(x) if isinstance(x, str) else x
    if math.isnan(x):
        return ""
    return f"{x:.{digits}f}"


def fmt_p(x):
    x = num(x) if isinstance(x, str) else x
    if math.isnan(x):
        return ""
    if x < 0.001:
        return f"{x:.2e}"
    return f"{x:.3f}"


strict_path = MR / "MVMR_STRICT_FULL_RESULTS.tsv"
strict = read_tsv(strict_path)

# Correct Egger-intercept multiplicity: 15 unique outcomes within each orientation,
# not the duplicated 30 exposure rows.
for orientation in ("BMI_increasing", "SMK_increasing"):
    unique = []
    seen = set()
    for row in strict:
        if row["method"] == "MVMR-Egger" and row["orientation"] == orientation and row["outcome"] not in seen:
            seen.add(row["outcome"])
            unique.append(row)
    qvals = bh([num(row["intercept_P"]) for row in unique])
    qmap = {row["outcome"]: q for row, q in zip(unique, qvals)}
    for row in strict:
        if row["method"] == "MVMR-Egger" and row["orientation"] == orientation:
            row["intercept_FDR_q"] = f"{qmap[row['outcome']]:.15g}"

write_tsv(strict_path, strict, list(strict[0]))

# Additional file 5: strict primary results plus complete four-model IVW sensitivity.
add5_fields = [
    "method", "exposure", "outcome", "n_instruments", "beta", "SE", "P", "FDR_q",
    "CI_lower", "CI_upper", "multiplicity_family", "uncertainty_method", "orientation",
    "intercept", "intercept_SE", "intercept_P", "intercept_FDR_q", "intercept_CI_lower",
    "intercept_CI_upper", "instrument_set", "analysis_role",
]
add5 = []
for row in strict:
    method = row["method"]
    n = int(row["n_instruments"])
    out = {field: "" for field in add5_fields}
    out.update({key: row.get(key, "") for key in ("exposure", "outcome", "n_instruments", "beta", "SE", "P", "FDR_q", "CI_lower", "CI_upper", "orientation", "intercept", "intercept_SE", "intercept_P", "intercept_FDR_q", "instrument_set")})
    out["method"] = method
    out["analysis_role"] = "STRICT_PRIMARY"
    if method == "MVMR-Egger":
        out["multiplicity_family"] = "30 exposure-by-outcome effects and 15 intercepts within orientation"
        out["uncertainty_method"] = "MVMR-Egger model-based SE; t interval with n-3 df"
        crit = student_t.ppf(0.975, n - 3)
        inter = num(row["intercept"])
        inter_se = num(row["intercept_SE"])
        out["intercept_CI_lower"] = f"{inter - crit * inter_se:.15g}"
        out["intercept_CI_upper"] = f"{inter + crit * inter_se:.15g}"
    elif method == "MVMR-IVW":
        out["multiplicity_family"] = "30 exposure-by-outcome estimates within method"
        out["uncertainty_method"] = "MVMR-IVW model-based SE; t interval with n-2 df"
    elif method == "MVMR weighted LAD median":
        out["multiplicity_family"] = "30 exposure-by-outcome estimates within method"
        out["uncertainty_method"] = "2000-replicate covariance-aware parametric bootstrap; normal interval"
    else:
        out["multiplicity_family"] = "30 exposure-by-outcome estimates within method"
        out["uncertainty_method"] = "Leave-one-SNP-out jackknife SE; normal interval"
    add5.append(out)

model_rows = read_tsv(MR / "MVMR_PALINDROME_SENSITIVITY_RESULTS.tsv")
for row in model_rows:
    out = {field: "" for field in add5_fields}
    out.update({key: row[key] for key in ("exposure", "outcome", "n_instruments", "beta", "SE", "P", "FDR_q", "CI_lower", "CI_upper")})
    out["method"] = "MVMR-IVW instrument-set sensitivity"
    out["multiplicity_family"] = "30 exposure-by-outcome estimates within instrument model"
    out["uncertainty_method"] = "MVMR-IVW model-based SE; t interval with n-2 df"
    out["orientation"] = row["model"]
    out["instrument_set"] = row["model"]
    out["analysis_role"] = "PRIMARY" if row["model"] == "C_STRICT_406" else "SENSITIVITY"
    add5.append(out)

write_tsv(PKG / "07_Additional_file_5_MVMR_Full_Results.tsv", add5, add5_fields)

# Supplementary Table S5.
ivw = [r for r in strict if r["method"] == "MVMR-IVW"]
by_outcome = defaultdict(dict)
outcome_order = []
for row in ivw:
    if row["outcome"] not in by_outcome:
        outcome_order.append(row["outcome"])
    by_outcome[row["outcome"]][row["exposure"]] = row
s5_fields = ["Outcome", "n", "BMI beta", "BMI SE", "BMI CI low", "BMI CI high", "BMI P", "BMI q", "Smoking beta", "Smoking SE", "Smoking CI low", "Smoking CI high", "Smoking P", "Smoking q"]
s5 = []
for outcome in outcome_order:
    b = by_outcome[outcome]["BMI"]
    s = by_outcome[outcome]["Smoking initiation"]
    s5.append({
        "Outcome": outcome, "n": b["n_instruments"],
        "BMI beta": fmt(b["beta"]), "BMI SE": fmt(b["SE"]), "BMI CI low": fmt(b["CI_lower"]), "BMI CI high": fmt(b["CI_upper"]), "BMI P": fmt_p(b["P"]), "BMI q": fmt_p(b["FDR_q"]),
        "Smoking beta": fmt(s["beta"]), "Smoking SE": fmt(s["SE"]), "Smoking CI low": fmt(s["CI_lower"]), "Smoking CI high": fmt(s["CI_upper"]), "Smoking P": fmt_p(s["P"]), "Smoking q": fmt_p(s["FDR_q"]),
    })
write_tsv(TABLES / "TableS5_MVMR_IVW.tsv", s5, s5_fields)

# Supplementary Table S6.
egger = [r for r in strict if r["method"] == "MVMR-Egger"]
eg = defaultdict(dict)
eg_order = []
for row in egger:
    key = (row["outcome"], row["orientation"])
    if key not in eg:
        eg_order.append(key)
    eg[key][row["exposure"]] = row
s6_fields = ["Outcome", "Orientation", "n", "BMI beta (95% CI)", "BMI P", "BMI q", "Smoking beta (95% CI)", "Smoking P", "Smoking q", "Intercept (95% CI)", "Intercept P", "Intercept q"]
s6 = []
for key in eg_order:
    b = eg[key]["BMI"]
    s = eg[key]["Smoking initiation"]
    n = int(b["n_instruments"])
    crit = student_t.ppf(0.975, n - 3)
    inter = num(b["intercept"])
    ise = num(b["intercept_SE"])
    s6.append({
        "Outcome": key[0], "Orientation": "BMI increasing" if key[1] == "BMI_increasing" else "SMK increasing", "n": n,
        "BMI beta (95% CI)": f"{fmt(b['beta'])} ({fmt(b['CI_lower'])}, {fmt(b['CI_upper'])})", "BMI P": fmt_p(b["P"]), "BMI q": fmt_p(b["FDR_q"]),
        "Smoking beta (95% CI)": f"{fmt(s['beta'])} ({fmt(s['CI_lower'])}, {fmt(s['CI_upper'])})", "Smoking P": fmt_p(s["P"]), "Smoking q": fmt_p(s["FDR_q"]),
        "Intercept (95% CI)": f"{inter:.4f} ({inter-crit*ise:.4f}, {inter+crit*ise:.4f})", "Intercept P": fmt_p(b["intercept_P"]), "Intercept q": fmt_p(b["intercept_FDR_q"]),
    })
write_tsv(TABLES / "TableS6_Directional_MVMR_Egger.tsv", s6, s6_fields)

# Supplementary Table S8.
s8_fields = ["method", "exposure", "outcomes", "positive", "FDR_significant", "median_beta"]
s8 = []
for method in ("MVMR-IVW", "MVMR weighted LAD median", "Q-minimisation"):
    for exposure in ("BMI", "Smoking initiation"):
        vals = [r for r in strict if r["method"] == method and r["exposure"] == exposure]
        betas = sorted(num(r["beta"]) for r in vals)
        median = betas[len(betas)//2]
        s8.append({"method": method, "exposure": exposure, "outcomes": len(vals), "positive": sum(num(r["beta"]) > 0 for r in vals), "FDR_significant": sum(num(r["FDR_q"]) < 0.05 for r in vals), "median_beta": f"{median:.15g}"})
write_tsv(TABLES / "TableS8_Robust_MVMR_summary.tsv", s8, s8_fields)

# Strict heterogeneity and LOO attachment/Table S9.
het = {r["outcome"]: r for r in read_tsv(MR / "MVMR_STRICT_HETEROGENEITY.tsv")}
loo_summary = {r["outcome"]: r for r in read_tsv(MR / "MVMR_STRICT_LEAVE_ONE_OUT_SUMMARY.tsv")}
loo_rows = read_tsv(MR / "MVMR_STRICT_LEAVE_ONE_OUT_PER_SNP.tsv")
loo_by = defaultdict(list)
for row in loo_rows:
    loo_by[row["outcome"]].append(row)
influential = {}
for outcome, rows in loo_by.items():
    summary = loo_summary[outcome]
    fb = num(summary["BMI_full_beta"])
    fs = num(summary["Smoking_full_beta"])
    max_b = max(abs(num(r["BMI_beta"]) - fb) for r in rows)
    max_s = max(abs(num(r["Smoking_beta"]) - fs) for r in rows)
    ids = []
    for row in rows:
        b = num(row["BMI_beta"]); s = num(row["Smoking_beta"])
        if abs(abs(b-fb)-max_b) < 1e-14 or abs(abs(s-fs)-max_s) < 1e-14 or (b * fb < 0) or (s * fs < 0):
            if row["omitted_SNP"] not in ids:
                ids.append(row["omitted_SNP"])
    influential[outcome] = ";".join(ids)

add7_fields = ["record_type", "outcome", "n_instruments", "heterogeneity_method", "Q", "Q_df", "Q_P", "Q_FDR", "nominal_significant", "FDR_significant", "BMI_full_beta", "BMI_LOO_min_beta", "BMI_LOO_max_beta", "BMI_LOO_sign_reversal_count", "Smoking_full_beta", "Smoking_LOO_min_beta", "Smoking_LOO_max_beta", "Smoking_LOO_sign_reversal_count", "influential_variant_ids", "convergence_status", "omitted_SNP", "success", "BMI_beta", "SMK_beta", "message", "instrument_set"]
add7 = []
for outcome in outcome_order:
    h = het[outcome]; l = loo_summary[outcome]
    row = {field: "" for field in add7_fields}
    row.update({"record_type": "OUTCOME_SUMMARY", "outcome": outcome, "n_instruments": h["n_instruments"], "heterogeneity_method": "MVMR::pleiotropy_mvmr IVW Q with SNP-specific exposure-association covariance", "Q": h["Q"], "Q_df": h["Q_df"], "Q_P": h["Q_P"], "Q_FDR": h["Q_FDR"], "nominal_significant": str(num(h["Q_P"]) < 0.05).upper(), "FDR_significant": str(num(h["Q_FDR"]) < 0.05).upper(), "influential_variant_ids": influential[outcome], "instrument_set": "STRICT_PRIMARY_406_EXCLUDING_17_FREQUENCY_AMBIGUOUS_PALINDROMIC_SNPS"})
    for key in ("BMI_full_beta", "BMI_LOO_min_beta", "BMI_LOO_max_beta", "BMI_LOO_sign_reversal_count", "Smoking_full_beta", "Smoking_LOO_min_beta", "Smoking_LOO_max_beta", "Smoking_LOO_sign_reversal_count", "convergence_status"):
        row[key] = l[key]
    add7.append(row)
for source in loo_rows:
    row = {field: "" for field in add7_fields}
    row.update({"record_type": "LEAVE_ONE_SNP_OUT", "outcome": source["outcome"], "n_instruments": loo_summary[source["outcome"]]["n_instruments"], "omitted_SNP": source["omitted_SNP"], "success": source["success"], "BMI_beta": source["BMI_beta"], "SMK_beta": source["Smoking_beta"], "message": source.get("message", ""), "instrument_set": "STRICT_PRIMARY_406_EXCLUDING_17_FREQUENCY_AMBIGUOUS_PALINDROMIC_SNPS"})
    add7.append(row)
write_tsv(PKG / "09_Additional_file_7_MVMR_Heterogeneity_and_LOO.tsv", add7, add7_fields)

s9_fields = ["Outcome", "SNPs", "IVW Q", "df", "Q P", "Q FDR", "BMI full beta", "BMI LOO range", "BMI sign-reversal omissions", "Smoking full beta", "Smoking LOO range", "Smoking sign-reversal omissions", "Influential variants", "Convergence"]
s9 = []
for outcome in outcome_order:
    h = het[outcome]; l = loo_summary[outcome]
    s9.append({"Outcome": outcome, "SNPs": h["n_instruments"], "IVW Q": f"{num(h['Q']):.2f}", "df": h["Q_df"], "Q P": fmt_p(h["Q_P"]), "Q FDR": fmt_p(h["Q_FDR"]), "BMI full beta": fmt(l["BMI_full_beta"]), "BMI LOO range": f"{fmt(l['BMI_LOO_min_beta'])} to {fmt(l['BMI_LOO_max_beta'])}", "BMI sign-reversal omissions": l["BMI_LOO_sign_reversal_count"], "Smoking full beta": fmt(l["Smoking_full_beta"]), "Smoking LOO range": f"{fmt(l['Smoking_LOO_min_beta'])} to {fmt(l['Smoking_LOO_max_beta'])}", "Smoking sign-reversal omissions": l["Smoking_LOO_sign_reversal_count"], "Influential variants": influential[outcome] or "None", "Convergence": l["convergence_status"]})
write_tsv(TABLES / "TableS9_MVMR_heterogeneity_and_LOO.tsv", s9, s9_fields)

# Additional file 6: instrument-set dictionary.
dict_fields = ["instrument_set", "analysis_label", "count", "selection_or_definition", "exposure_scale", "outcome_scale", "covariance_definition", "status"]
dictionary = [
    ["Strict primary MVMR", "C_STRICT_406", "406 unique; 405 for four outcomes", "Historical six-axis union after complete BMI-smoking data and F>10, excluding 17 unresolved palindromic variants", "BMI per 1-SD inverse-normal residual; smoking-initiation source-GWAS scale", "Disease log-odds", "SNP-specific 2×2 exposure-association sampling covariance", "CURRENT_PRIMARY"],
    ["Historical full MVMR", "A_FULL_423", "423 unique; 422 for four outcomes", "Six-exposure P<5×10−8 union; joint r²<0.001/10 Mb clump; complete BMI and smoking associations; per-variant BMI or smoking F>10", "Same as strict primary", "Disease log-odds", "Same covariance construction", "HISTORICAL_SENSITIVITY"],
    ["Historical exposure-GWS subset", "B_EXPOSURE_GWS_412", "412 unique; 411 for four outcomes", "A_FULL_423 restricted to BMI_GWS or SMK_GWS; not independently reclumped", "Same as strict primary", "Disease log-odds", "Same covariance construction", "SENSITIVITY"],
    ["Strict exposure-GWS subset", "D_EXPOSURE_GWS_STRICT_395", "395 unique; 394 for four outcomes", "Intersection of C_STRICT_406 with BMI_GWS or SMK_GWS", "Same as strict primary", "Disease log-odds", "Same covariance construction", "SENSITIVITY"],
    ["Independent BMI UVMR", "INDEPENDENT_BMI_521", "521 post-clump; 504–506 per outcome", "BMI GWAS P<5×10−8; BMI-only clump and outcome-specific harmonisation", "BMI per 1-SD inverse-normal residual", "Disease log-odds", "Not applicable", "CURRENT"],
    ["Independent smoking UVMR", "INDEPENDENT_SMK_92", "92 post-clump; 84 per outcome", "Smoking GWAS P<5×10−8; smoking-only clump and outcome-specific harmonisation", "Smoking-initiation source-GWAS scale", "Disease log-odds", "Not applicable", "CURRENT"],
    ["Historical union-subset UVMR", "UNION_SUBSET_UVMR", "404 BMI; 16 smoking", "Complete historical MVMR matrix then exposure-specific genome-wide significance", "Same as parent exposure", "Disease log-odds", "Not applicable", "HISTORICAL_SENSITIVITY"],
    ["Historical MR-AHC set", "MRAHC_HISTORICAL_WITHDRAWN_54", "54 explicit variants", "Six-exposure clumped union; smoking-specific F>10; complete four-outcome profiles", "Profile denominator", "Disease log-odds", "Not used for current effect estimation", "HISTORICAL_WITHDRAWN"],
]
write_tsv(PKG / "08_Additional_file_6_MR_Instrument_Data_Dictionary.tsv", [dict(zip(dict_fields, row)) for row in dictionary], dict_fields)

# Additional file 8: preserve all existing rows and add explicit A-D membership/current role.
manifest_path = PKG / "10_Additional_file_8_MR_Instrument_Manifest.tsv"
manifest = read_tsv(manifest_path)
manifest_fields = list(manifest[0])
membership_rows = read_tsv(PAL / "MVMR_STRICT_INSTRUMENT_MANIFEST.tsv")
membership = defaultdict(dict)
for row in membership_rows:
    membership[row["SNP"]][row["model"]] = row
extra = ["MODEL_A_423", "MODEL_B_412", "MODEL_C_STRICT_406", "MODEL_D_STRICT_GWS_395", "strict_primary_role", "palindromic_exclusion_reason"]
models = [("A_FULL_423", "MODEL_A_423"), ("B_EXPOSURE_GWS_412", "MODEL_B_412"), ("C_STRICT_406", "MODEL_C_STRICT_406"), ("D_EXPOSURE_GWS_STRICT_395", "MODEL_D_STRICT_GWS_395")]
for row in manifest:
    snp = row["SNP"]
    for model, field in models:
        value = membership.get(snp, {}).get(model, {}).get("included", "")
        row[field] = value
    row["strict_primary_role"] = "CURRENT_PRIMARY" if row["MODEL_C_STRICT_406"] == "TRUE" else ("EXCLUDED_UNRESOLVED_PALINDROME" if snp in membership and row["MODEL_A_423"] == "TRUE" else "NOT_APPLICABLE")
    row["palindromic_exclusion_reason"] = membership.get(snp, {}).get("C_STRICT_406", {}).get("exclusion_reason", "")
write_tsv(manifest_path, manifest, manifest_fields + extra)

print(f"strict_rows={len(strict)} add5_rows={len(add5)} add7_rows={len(add7)} manifest_rows={len(manifest)}")
