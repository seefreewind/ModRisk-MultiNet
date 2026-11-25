#!/usr/bin/env python3
"""Apply final multiplicity families and assemble auditable MR result products."""

from __future__ import annotations

from pathlib import Path
import math
import os

import pandas as pd


ROOT = Path(os.environ.get("MODRISK_ROOT", Path.cwd())).resolve()
MVMR = ROOT / "final_major_revision/05_robust_mvmr"
UVMR = ROOT / "final_major_revision/04_univariable_mr"
HARM = ROOT / "final_major_revision/03_mr_harmonisation"


def bh(values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg adjusted values, preserving missing entries and order."""
    out = pd.Series(float("nan"), index=values.index, dtype=float)
    ok = values.notna()
    p = values.loc[ok].astype(float)
    order = p.sort_values().index
    ranked = p.loc[order]
    n = len(ranked)
    adjusted = ranked * n / pd.Series(range(1, n + 1), index=ranked.index)
    adjusted = adjusted.iloc[::-1].cummin().iloc[::-1].clip(upper=1.0)
    out.loc[adjusted.index] = adjusted
    return out


def q30_wide(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path, sep="\t")
    stacked = pd.concat(
        [d[["BMI_P"]].rename(columns={"BMI_P": "P"}),
         d[["SMK_P"]].rename(columns={"SMK_P": "P"})],
        ignore_index=True,
    )
    q = bh(stacked["P"])
    n = len(d)
    d["BMI_q"] = q.iloc[:n].to_numpy()
    d["SMK_q"] = q.iloc[n:].to_numpy()
    d.to_csv(path, sep="\t", index=False)
    return d


def fmt(x: float, digits: int = 3) -> str:
    if pd.isna(x):
        return "NA"
    if x != 0 and abs(x) < 0.001:
        return f"{x:.2e}"
    return f"{x:.{digits}f}"


def postprocess_mvmr() -> None:
    ivw = q30_wide(MVMR / "MVMR_IVW_FINAL.tsv")
    median = q30_wide(MVMR / "MVMR_MEDIAN_FULL_RESULTS.tsv")
    qmin = q30_wide(MVMR / "MVMR_WEAK_IV_AWARE_RESULTS.tsv")

    egger_path = MVMR / "MVMR_EGGER_ORIENTATION_RESULTS.tsv"
    egger = pd.read_csv(egger_path, sep="\t")
    for orientation, idx in egger.groupby("orientation", sort=False).groups.items():
        g = egger.loc[idx]
        q = bh(pd.concat([g["BMI_P"], g["SMK_P"]], ignore_index=True))
        n = len(g)
        egger.loc[idx, "BMI_q"] = q.iloc[:n].to_numpy()
        egger.loc[idx, "SMK_q"] = q.iloc[n:].to_numpy()
        egger.loc[idx, "intercept_q"] = bh(g["intercept_P"]).to_numpy()
    egger.to_csv(egger_path, sep="\t", index=False)

    long = []
    for method, d in [
        ("MVMR-IVW", ivw),
        ("Multivariable weighted LAD median", median),
        ("Q-minimisation", qmin),
    ]:
        for exposure, prefix in [("BMI", "BMI"), ("Smoking initiation", "SMK")]:
            z = pd.DataFrame({
                "method": method,
                "exposure": exposure,
                "outcome": d["outcome"],
                "n_instruments": d["n_instruments"],
                "beta": d[f"{prefix}_beta"],
                "SE": d[f"{prefix}_SE"],
                "P": d[f"{prefix}_P"],
                "FDR_q": d[f"{prefix}_q"],
            })
            z["CI_lower"] = z["beta"] - 1.96 * z["SE"]
            z["CI_upper"] = z["beta"] + 1.96 * z["SE"]
            z["multiplicity_family"] = "30 exposure-by-outcome estimates within method"
            if method == "Multivariable weighted LAD median":
                z["uncertainty_method"] = "2,000-replicate covariance-aware parametric bootstrap"
            elif method == "Q-minimisation":
                z["uncertainty_method"] = "MVMR qhet_mvmr jackknife"
            else:
                z["uncertainty_method"] = "MVMR IVW model-based standard error"
            long.append(z)
    comparison = pd.concat(long, ignore_index=True)
    comparison.to_csv(MVMR / "MVMR_ROBUSTNESS_COMPARISON.tsv", sep="\t", index=False)

    full = comparison.copy()
    full["orientation"] = "common smoking effect allele"
    egger_long = []
    for exposure, prefix in [("BMI", "BMI"), ("Smoking initiation", "SMK")]:
        z = pd.DataFrame({
            "method": "MVMR-Egger",
            "exposure": exposure,
            "outcome": egger["outcome"],
            "n_instruments": egger["n_instruments"],
            "beta": egger[f"{prefix}_beta"],
            "SE": egger[f"{prefix}_SE"],
            "P": egger[f"{prefix}_P"],
            "FDR_q": egger[f"{prefix}_q"],
            "orientation": egger["orientation"],
            "intercept": egger["intercept"],
            "intercept_SE": egger["intercept_SE"],
            "intercept_P": egger["intercept_P"],
            "intercept_FDR_q": egger["intercept_q"],
            "multiplicity_family": "30 exposure-by-outcome estimates within orientation; 15 intercepts within orientation",
            "uncertainty_method": "weighted regression model-based standard error",
        })
        z["CI_lower"] = z["beta"] - 1.96 * z["SE"]
        z["CI_upper"] = z["beta"] + 1.96 * z["SE"]
        egger_long.append(z)
    full = pd.concat([full, *egger_long], ignore_index=True, sort=False)
    full.to_csv(MVMR / "MVMR_FULL_RESULTS.tsv", sep="\t", index=False)

    def counts(d: pd.DataFrame) -> tuple[int, int, int, int]:
        return (
            int((d["BMI_q"] < 0.05).sum()),
            int((d["SMK_q"] < 0.05).sum()),
            int((d["BMI_beta"] > 0).sum()),
            int((d["SMK_beta"] > 0).sum()),
        )

    rows = []
    for label, d in [("MVMR-IVW", ivw), ("Multivariable weighted LAD median", median), ("Q-minimisation", qmin)]:
        a, b, c, e = counts(d)
        rows.append((label, a, b, c, e))
    count_table = pd.DataFrame(rows, columns=["method", "BMI_FDR_significant", "SMK_FDR_significant", "BMI_positive", "SMK_positive"])
    count_table.to_csv(MVMR / "MVMR_METHOD_SUMMARY_COUNTS.tsv", sep="\t", index=False)

    priority = []
    for label, d in [("MVMR-IVW", ivw), ("Multivariable weighted LAD median", median), ("Q-minimisation", qmin)]:
        for outcome in ["HF", "coronary_heart"]:
            r = d.loc[d.outcome.eq(outcome)].iloc[0]
            priority.append(
                f"- {label}, {outcome}: BMI {fmt(r.BMI_beta)} (SE {fmt(r.BMI_SE)}, q={fmt(r.BMI_q)}); "
                f"smoking initiation {fmt(r.SMK_beta)} (SE {fmt(r.SMK_SE)}, q={fmt(r.SMK_q)})."
            )

    e_counts = []
    for orientation, g in egger.groupby("orientation", sort=False):
        e_counts.append(
            f"- {orientation}: BMI {int((g.BMI_q < .05).sum())}/15 and smoking initiation "
            f"{int((g.SMK_q < .05).sum())}/15 FDR-significant; "
            f"{int((g.intercept_P < .05).sum())}/15 intercepts had nominal P < 0.05."
        )

    report = [
        "# MVMR final sensitivity report", "",
        "## Multiplicity", "",
        "For IVW, multivariable weighted LAD median and Q-minimisation, Benjamini-Hochberg correction was applied to the complete 30-estimate exposure-by-outcome family within each method. Directional MVMR-Egger used a separate 30-estimate family for each orientation and a 15-intercept family for each orientation.", "",
        "## Method-level findings", "",
    ]
    for row in count_table.itertuples(index=False):
        report.append(
            f"- {row.method}: BMI {row.BMI_FDR_significant}/15 and smoking initiation {row.SMK_FDR_significant}/15 FDR-significant; positive directions {row.BMI_positive}/15 and {row.SMK_positive}/15, respectively."
        )
    report += ["", "## Priority cardiovascular outcomes", "", *priority, "", "## Directional Egger", "", *e_counts, "", "## Interpretation", "", "With MVMR 0.4.9 and the audited SNP-specific covariance matrices, the conditional F statistics were 25.48 for BMI and 3.85 for smoking initiation. Earlier revision-stage MVMR 0.4.8 output was 25.42 and 3.84; this small package-version difference does not change the strength classification. The low smoking conditional F can affect the entire mutually adjusted model; a stronger BMI conditional F does not guarantee that BMI estimates are free of weak-instrument distortion. Bias need not be towards the null, and unavailable exposure-outcome sampling covariance leaves additional uncertainty when source samples overlap. The robust methods use different identifying assumptions and are sensitivity analyses, not corrections that eliminate these limitations.", ""]
    (MVMR / "MVMR_FINAL_SENSITIVITY_REPORT.md").write_text("\n".join(report), encoding="utf-8")


def postprocess_uvmr() -> None:
    independent = pd.read_csv(UVMR / "UNIVARIABLE_MR_INDEPENDENT.tsv", sep="\t")
    union = pd.read_csv(UVMR / "UNIVARIABLE_MR_UNION_SUBSET.tsv", sep="\t")
    full = pd.concat([union, independent], ignore_index=True, sort=False)
    full.to_csv(UVMR / "UVMR_FULL_RESULTS.tsv", sep="\t", index=False)

    comp = pd.read_csv(UVMR / "UVMR_INSTRUMENT_COMPARISON.tsv", sep="\t")
    comp.loc[comp.exposure.eq("BMI"), "explanation"] = (
        "Independent BMI instruments were selected directly from the full BMI GWAS; the historical subset was restricted to variants retained in the complete BMI-smoking MVMR matrix."
    )
    comp.loc[comp.exposure.eq("Smoking initiation"), "explanation"] = (
        "The 54-variant MR-AHC set came from the BMI-smoking union after joint clumping and required smoking-specific per-variant F>10; only 27 were smoking genome-wide significant. The 16-variant historical UVMR subset required smoking genome-wide significance and complete joint-matrix data. Independent UVMR starts from the full smoking GWAS and clumps smoking loci alone."
    )
    comp.to_csv(UVMR / "UVMR_INSTRUMENT_COMPARISON.tsv", sep="\t", index=False)
    key = independent.loc[
        independent["outcome"].isin(["HF", "coronary_heart", "diabetes_type_2", "COPD", "OA"])
        & independent["method"].eq("IVW"),
        ["exposure", "outcome", "nsnp", "mean_F", "beta", "SE", "P", "FDR_q"],
    ]
    counts = (
        independent.assign(positive=independent.beta > 0, significant=independent.FDR_q < .05)
        .groupby(["exposure", "method"], as_index=False)
        .agg(n_outcomes=("outcome", "size"), positive=("positive", "sum"), FDR_significant=("significant", "sum"))
    )
    counts.to_csv(UVMR / "UVMR_METHOD_SUMMARY_COUNTS.tsv", sep="\t", index=False)

    lines = [
        "# Univariable MR sensitivity report", "", "## Design", "",
        "Independent exposure-specific instruments were selected from the complete BMI and smoking-initiation GWAS files at P < 5 × 10^-8 and clumped with the frozen 1000 Genomes European reference at r2 < 0.001 within 10 Mb. Outcome effects were harmonised independently for all 15 outcomes. Palindromic intermediate-frequency variants and unresolved allele mismatches were excluded with recorded reasons.", "",
        "Benjamini-Hochberg correction was applied to 15 outcomes within each exposure-by-method combination for the independent analysis. Historical union-subset q values remain in their original, separately labelled families.", "",
        "## Instrument provenance", "",
        "The independent analysis retained 521 BMI and 92 smoking-initiation variants after exposure-specific clumping; 504–506 BMI variants and 84 smoking variants were usable per outcome after harmonisation. The historical union-subset analysis used 404 BMI and 16 smoking variants because it required membership and completeness in the joint BMI-smoking matrix.", "",
        "The 54-variant MR-AHC set had a different provenance: it was selected from the BMI-smoking union after joint clumping and required smoking-specific per-variant F > 10. Only 27 of the 54 variants were genome-wide significant for smoking initiation. It is therefore neither a smoking-specific genome-wide-significant instrument set nor interchangeable with the independent UVMR set.", "",
        "## Method-level counts", "",
        counts.to_markdown(index=False), "", "## Priority outcomes under independent IVW", "", key.to_markdown(index=False), "", "## Interpretation", "",
        "Univariable MR estimates total exposure associations and does not estimate the same mutually adjusted direct effects as MVMR. Agreement or disagreement therefore reflects estimand and instrument-set dependence rather than formal validation. All 15 independent smoking IVW estimates were positive and FDR-significant, including positive estimates for heart failure and coronary heart disease. Weighted-median estimates supported those two directions, whereas weighted-mode estimates were less precise and often non-significant.", "",
    ]
    (UVMR / "UVMR_SENSITIVITY_REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    required = MVMR / "MVMR_WEAK_IV_AWARE_RESULTS.tsv"
    if not required.exists():
        raise SystemExit(f"Waiting for required Q-minimisation output: {required}")
    postprocess_mvmr()
    postprocess_uvmr()
    print("Final MR outputs postprocessed")


if __name__ == "__main__":
    main()
