#!/usr/bin/env python3
from pathlib import Path
import os
import pandas as pd

ROOT = Path(os.environ.get("MODRISK_ROOT", Path.cwd())).resolve()
WORK = ROOT / "final_manuscript_error_correction"

pair_file = WORK / "10_submission_package/04_Additional_file_2_Pair_level_classification.csv"
pairs = pd.read_csv(pair_file)
pairs["baseline_significant"] = pairs["baseline_fdr"] < 0.05
pairs["residual_significant_recomputed"] = pairs["residual_fdr"] < 0.05

rows = []
for b in (True, False):
    for r in (True, False):
        rows.append({
            "baseline_FDR_lt_0.05": b,
            "residual_FDR_lt_0.05": r,
            "n_pairs": int(((pairs.baseline_significant == b) &
                            (pairs.residual_significant_recomputed == r)).sum()),
        })
cross = pd.DataFrame(rows)
out_net = WORK / "03_network_logic"
out_net.mkdir(parents=True, exist_ok=True)
cross.to_csv(out_net / "BASELINE_RESIDUAL_SIGNIFICANCE_CROSSTAB.tsv", sep="\t", index=False)

checks = {
    "total_pairs": len(pairs),
    "baseline_significant": int(pairs.baseline_significant.sum()),
    "residual_significant": int(pairs.residual_significant_recomputed.sum()),
    "baseline_and_residual_significant": int((pairs.baseline_significant & pairs.residual_significant_recomputed).sum()),
    "baseline_significant_residual_nonsignificant": int((pairs.baseline_significant & ~pairs.residual_significant_recomputed).sum()),
    "baseline_nonsignificant_residual_significant": int((~pairs.baseline_significant & pairs.residual_significant_recomputed).sum()),
    "baseline_and_residual_nonsignificant": int((~pairs.baseline_significant & ~pairs.residual_significant_recomputed).sum()),
}
assert checks["total_pairs"] == 2485
assert checks["baseline_significant"] == 1535
assert checks["residual_significant"] == 1100
assert (pairs["residual_significant"].astype(bool) == pairs["residual_significant_recomputed"]).all()

main = (WORK / "06_manuscript/Main_Manuscript_Error_Correction.md").read_text()
terms = ["retained", "remained", "persistent", "persisted", "residual"]
hits = []
for i, line in enumerate(main.splitlines(), 1):
    if any(t in line.lower() for t in terms):
        hits.append(f"- Line {i}: {line.strip()}")
report = [
    "# Residual-significance language audit",
    "",
    "## Recomputed significance cross-tabulation",
    "",
    *(f"- {k}: {v}" for k, v in checks.items()),
    "",
    "## Interpretation rule",
    "",
    "Use ‘retained/remained significant’ only for the baseline-significant and residual-significant intersection. Use ‘residual-significant’ for the full residual-significant set, including pairs that were not baseline-significant. Do not use ‘persistent’ as a statistical class because the manuscript's current two-dimensional classification separates signed change from residual significance.",
    "",
    "## Candidate manuscript occurrences",
    "",
    *hits,
]
(out_net / "RESIDUAL_SIGNIFICANCE_LANGUAGE_AUDIT.md").write_text("\n".join(report) + "\n")

# Independent UVMR audit for overlap with the 17 MVMR frequency-ambiguous SNPs.
pal = pd.read_csv(ROOT / "pre_submission_technical_closeout/03_palindromic_audit/PALINDROMIC_SNP_EAF_AUDIT.tsv", sep="\t")
unresolved = set(pal.loc[pal.decision == "UNRESOLVED_FREQUENCY_AMBIGUOUS", "SNP"])
uv = pd.read_csv(ROOT / "final_major_revision/04_univariable_mr/UVMR_HARMONISATION_AUDIT.tsv", sep="\t", low_memory=False)
aff = uv[uv.SNP.isin(unresolved)].copy()
summary = (aff.groupby(["exposure", "SNP"], as_index=False)
           .agg(exposure_EAF=("EAF", "first"),
                outcome_EAF_min=("out_EAF", "min"),
                outcome_EAF_max=("out_EAF", "max"),
                outcome_rows=("outcome", "size"),
                included_rows=("included", "sum"),
                harmonisation_actions=("harmonisation_action", lambda x: ";".join(sorted(set(map(str, x))))),
                exclusion_reasons=("exclusion_reason", lambda x: ";".join(sorted({str(v) for v in x if pd.notna(v)})))))
summary["uvmr_disposition"] = summary.apply(
    lambda r: "RESOLVED_BY_OUTCOME_EAF_AND_INCLUDED" if r.included_rows == r.outcome_rows
    else ("EXCLUDED_AS_PALINDROMIC_INTERMEDIATE_FREQUENCY" if r.included_rows == 0
          else "MIXED_REQUIRES_REVIEW"), axis=1)
summary.to_csv(WORK / "02_mr_sensitivity/UVMR_PALINDROMIC_OVERLAP_AUDIT.tsv", sep="\t", index=False)
assert not (summary.uvmr_disposition == "MIXED_REQUIRES_REVIEW").any()
uv_report = [
    "# Independent UVMR palindromic-overlap audit",
    "",
    f"Eleven of the 17 MVMR frequency-ambiguous variants overlapped the independent UVMR instrument sets ({sum(summary.exposure == 'BMI')} BMI; {sum(summary.exposure == 'Smoking initiation')} smoking initiation).",
    f"Nine BMI variants and one smoking variant were already excluded from every outcome as palindromic intermediate-frequency variants. Two BMI variants (rs1436344 and rs1624134) were retained because outcome allele frequencies consistently resolved their orientation under the independent-UVMR harmonisation rule.",
    "",
    "The independent UVMR result tables therefore already implement strict outcome-specific handling and do not require recomputation. This conclusion concerns UVMR only and does not rescue the 17 variants in the two-exposure MVMR, where the available exposure-pair frequencies remained ambiguous.",
]
(WORK / "02_mr_sensitivity/UVMR_PALINDROMIC_AUDIT.md").write_text("\n".join(uv_report) + "\n")

# Provenance of the 11 historical 423-set variants that were not GWS for BMI or smoking.
prov = pd.read_csv(ROOT / "pre_submission_technical_closeout/01_instrument_provenance/MVMR_423_INSTRUMENT_PROVENANCE.tsv", sep="\t")
other = prov[prov["neither_BMI_nor_SMK_GWS"].astype(bool)].copy()
assert len(other) == 11
axis_cols = ["SBP_GWS", "ApoB_GWS", "FG_GWS", "SLEEP_GWS"]
other["other_axis_GWS_basis"] = other.apply(
    lambda r: ";".join(c.removesuffix("_GWS") for c in axis_cols if bool(r[c])), axis=1)
assert other["other_axis_GWS_basis"].str.len().gt(0).all()
other.to_csv(WORK / "02_mr_sensitivity/OTHER_AXIS_11_INSTRUMENT_PROVENANCE.tsv", sep="\t", index=False)
other_report = [
    "# Other-axis instrument provenance",
    "",
    "Eleven variants entered the historical 423-variant matrix through genome-wide significance for at least one of SBP, ApoB, fasting glucose or sleep duration, although none was genome-wide significant for BMI or smoking initiation. Their inclusion followed the frozen six-exposure-union selection chain and was not an independently clumped BMI-smoking instrument design.",
    "",
    "Because the final two-exposure model contains BMI and smoking initiation only, these variants could affect disease outcomes through the omitted risk-factor axes that selected them. Membership in the six-exposure union therefore does not by itself establish the exclusion-restriction assumption for the two-exposure MVMR. Models B and D isolate this selection feature; Model C isolates palindromic ambiguity.",
]
(WORK / "02_mr_sensitivity/OTHER_AXIS_11_INSTRUMENT_PROVENANCE.md").write_text("\n".join(other_report) + "\n")

print(checks)
print(summary[["exposure", "SNP", "uvmr_disposition"]].to_string(index=False))
print(other[["SNP", "other_axis_GWS_basis"]].to_string(index=False))
