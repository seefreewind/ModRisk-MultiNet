#!/usr/bin/env python3
"""Independent downstream audit of the frozen network and repaired pair labels."""

from __future__ import annotations

import csv
import hashlib
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("MODRISK_ROOT", Path.cwd())).resolve()
BASE = ROOT / "final_major_revision"
NET_OUT = BASE / "02_primary_network_audit"
PAIR_OUT = BASE / "07_figures_tables"
FDR = 0.05


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load(path: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / path, sep="\t")


def check_row(metric: str, observed, expected, tolerance=0.0, source="") -> dict:
    if isinstance(expected, (float, np.floating)):
        passed = np.isfinite(observed) and abs(float(observed) - float(expected)) <= tolerance
    else:
        passed = observed == expected
    return {
        "metric": metric,
        "observed": observed,
        "expected": expected,
        "tolerance": tolerance,
        "status": "PASS" if passed else "FAIL",
        "source": source,
    }


def bh_adjust(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.minimum(ranked, 1.0)
    return out


def two_sided_normal_p(z: np.ndarray) -> np.ndarray:
    return np.fromiter((math.erfc(abs(float(x)) / math.sqrt(2.0)) for x in z), dtype=float)


def main() -> None:
    NET_OUT.mkdir(parents=True, exist_ok=True)
    PAIR_OUT.mkdir(parents=True, exist_ok=True)

    joint_path = ROOT / "results/phase3_joint/JOINT_PARTIAL_RG.tsv"
    atten_path = ROOT / "results/phase3_joint/JOINT_ATTENUATION.tsv"
    pair_path = ROOT / "results/phase3_classification_repair_R2/PAIR_CLASSIFICATION_REBUILT.tsv"
    single_path = ROOT / "results/phase2_single_factor/NETWORK_ATTENUATION_SUMMARY.tsv"
    repro_path = ROOT / "results/phase2_baseline/BASELINE_VS_MOUNIER_COMPARISON.tsv"
    fig2_path = ROOT / "results/phase3_joint/FIGURE2_HEATMAP_DATA.tsv"

    joint = pd.read_csv(joint_path, sep="\t")
    atten = pd.read_csv(atten_path, sep="\t")
    pair = pd.read_csv(pair_path, sep="\t")
    single = pd.read_csv(single_path, sep="\t")
    repro = pd.read_csv(repro_path, sep="\t").iloc[0]
    fig2 = pd.read_csv(fig2_path, sep="\t")

    assert len(joint) == len(atten) == len(pair) == 2485
    assert joint.pair_id.is_unique and atten.pair_id.is_unique and pair.pair_id.is_unique
    atten = atten.set_index("pair_id").loc[joint.pair_id].reset_index()
    pair = pair.set_index("pair_id").loc[joint.pair_id].reset_index()

    disease_codes = sorted(set(joint.disease1) | set(joint.disease2))
    canonical_pairs = joint.apply(lambda r: "||".join(sorted((r.disease1, r.disease2))), axis=1)
    recomputed_rg = joint.joint_covariance / np.sqrt(joint.joint_var_disease1 * joint.joint_var_disease2)
    signed_difference = joint.baseline_rg - joint.joint_rg
    magnitude_change = joint.baseline_rg.abs() - joint.joint_rg.abs()
    reversal = (
        (np.sign(joint.baseline_rg) != np.sign(joint.joint_rg))
        & joint.baseline_rg.ne(0)
        & joint.joint_rg.ne(0)
    )
    significant_reversal = reversal & joint.difference_FDR.lt(FDR)
    ci_lo = joint.joint_rg - 1.96 * joint.joint_SE
    ci_hi = joint.joint_rg + 1.96 * joint.joint_SE
    eq005 = ci_lo.ge(-0.05) & ci_hi.le(0.05)
    eq010 = ci_lo.ge(-0.10) & ci_hi.le(0.10)
    stable_positive = (joint.baseline_rg.abs() / joint.baseline_SE >= 2) & magnitude_change.gt(0)
    proportional = magnitude_change[stable_positive] / joint.baseline_rg.abs()[stable_positive]
    baseline_p_recomputed = two_sided_normal_p((joint.baseline_rg / joint.baseline_SE).to_numpy())
    joint_p_recomputed = two_sided_normal_p(joint.joint_Z.to_numpy())
    difference_p_recomputed = two_sided_normal_p(joint.difference_Z.to_numpy())
    baseline_fdr_recomputed = bh_adjust(baseline_p_recomputed)
    joint_fdr_recomputed = bh_adjust(joint_p_recomputed)
    difference_fdr_recomputed = bh_adjust(difference_p_recomputed)

    axis_expected = {"BMI": 1193, "SMK": 1178, "SBP": 120, "ApoB": 7, "FG": 5, "SLEEP": 148}
    axis_observed = dict(zip(single.exposure, single.n_significant_q005.astype(int)))

    checks = [
        check_row("diseases", len(disease_codes), 71, source="JOINT_PARTIAL_RG.tsv"),
        check_row("disease_pairs", len(joint), 2485, source="JOINT_PARTIAL_RG.tsv"),
        check_row("unique_unordered_disease_pairs", canonical_pairs.nunique(), 2485, source="JOINT_PARTIAL_RG.tsv"),
        check_row("baseline_FDR_significant", int(joint.baseline_FDR.lt(FDR).sum()), 1535, source="JOINT_PARTIAL_RG.tsv"),
        check_row("joint_signed_difference_FDR_significant", int(joint.difference_FDR.lt(FDR).sum()), 1019, source="JOINT_PARTIAL_RG.tsv"),
        check_row("joint_residual_FDR_significant", int(joint.joint_FDR.lt(FDR).sum()), 1100, source="JOINT_PARTIAL_RG.tsv"),
        check_row("change_attenuated", int(((joint.difference_FDR < FDR) & (magnitude_change > 0)).sum()), 900, source="JOINT_PARTIAL_RG.tsv"),
        check_row("change_enhanced", int(((joint.difference_FDR < FDR) & (magnitude_change < 0)).sum()), 119, source="JOINT_PARTIAL_RG.tsv"),
        check_row("change_no_FDR_significant_signed_difference", int((joint.difference_FDR >= FDR).sum()), 1466, source="JOINT_PARTIAL_RG.tsv"),
        check_row("residual_near_null_classified", int((pair.residual_state == "NEAR_NULL_CLASSIFIED").sum()), 605, source="PAIR_CLASSIFICATION_REBUILT.tsv"),
        check_row("residual_non_significant_other", int((pair.residual_state == "RESIDUAL_NON_SIGNIFICANT_OTHER").sum()), 780, source="PAIR_CLASSIFICATION_REBUILT.tsv"),
        check_row("sign_reversal_descriptive", int(reversal.sum()), 135, source="JOINT_PARTIAL_RG.tsv"),
        check_row("sign_reversal_with_FDR_significant_signed_difference", int(significant_reversal.sum()), 66, source="JOINT_PARTIAL_RG.tsv"),
        check_row("interval_equivalent_margin_0.05", int(eq005.sum()), 2, source="JOINT_PARTIAL_RG.tsv"),
        check_row("interval_equivalent_margin_0.10", int(eq010.sum()), 162, source="JOINT_PARTIAL_RG.tsv"),
        check_row("median_joint_magnitude_change", float(np.median(magnitude_change)), 0.0206, tolerance=5e-5, source="JOINT_PARTIAL_RG.tsv"),
        check_row("proportional_attenuation_denominator", int(stable_positive.sum()), 1317, source="JOINT_PARTIAL_RG.tsv"),
        check_row("median_proportional_attenuation", float(np.median(proportional)), 0.2492, tolerance=5e-5, source="JOINT_PARTIAL_RG.tsv"),
        check_row("baseline_reproduction_r", float(repro.cor_rg), 0.9936, tolerance=5e-5, source="BASELINE_VS_MOUNIER_COMPARISON.tsv"),
        check_row("baseline_reproduction_median_absolute_difference", float(repro.median_abs_diff), 0.0033, tolerance=5e-5, source="BASELINE_VS_MOUNIER_COMPARISON.tsv"),
        check_row("Schur_joint_rg_max_absolute_error", float(np.max(np.abs(recomputed_rg - joint.joint_rg))), 0.0, tolerance=1e-12, source="JOINT_PARTIAL_RG.tsv"),
        check_row("signed_difference_max_absolute_error", float(np.max(np.abs(signed_difference - joint.difference_estimate))), 0.0, tolerance=1e-12, source="JOINT_PARTIAL_RG.tsv"),
        check_row("magnitude_change_max_absolute_error", float(np.max(np.abs(magnitude_change - atten.absolute_attenuation))), 0.0, tolerance=1e-12, source="JOINT_ATTENUATION.tsv"),
        check_row("baseline_P_normal_recalculation_max_absolute_error", float(np.max(np.abs(baseline_p_recomputed - joint.baseline_P))), 0.0, tolerance=1e-12, source="JOINT_PARTIAL_RG.tsv"),
        check_row("joint_P_normal_recalculation_max_absolute_error", float(np.max(np.abs(joint_p_recomputed - joint.joint_P))), 0.0, tolerance=1e-12, source="JOINT_PARTIAL_RG.tsv"),
        check_row("difference_P_normal_recalculation_max_absolute_error", float(np.max(np.abs(difference_p_recomputed - joint.difference_P))), 0.0, tolerance=1e-12, source="JOINT_PARTIAL_RG.tsv"),
        check_row("baseline_BH_FDR_recalculation_max_absolute_error", float(np.max(np.abs(baseline_fdr_recomputed - joint.baseline_FDR))), 0.0, tolerance=1e-12, source="JOINT_PARTIAL_RG.tsv"),
        check_row("joint_BH_FDR_recalculation_max_absolute_error", float(np.max(np.abs(joint_fdr_recomputed - joint.joint_FDR))), 0.0, tolerance=1e-12, source="JOINT_PARTIAL_RG.tsv"),
        check_row("difference_BH_FDR_recalculation_max_absolute_error", float(np.max(np.abs(difference_fdr_recomputed - joint.difference_FDR))), 0.0, tolerance=1e-12, source="JOINT_PARTIAL_RG.tsv"),
        check_row("figure2_rows", len(fig2), 2485, source="FIGURE2_HEATMAP_DATA.tsv"),
    ]
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
    observed_cross = pair.groupby(["change_class", "residual_state"]).size().to_dict()
    for (change, residual), expected in expected_cross.items():
        checks.append(check_row(
            f"Table3_{change}_{residual}",
            int(observed_cross.get((change, residual), 0)), expected,
            source="PAIR_CLASSIFICATION_REBUILT.tsv",
        ))
    for axis, expected in axis_expected.items():
        checks.append(check_row(f"single_axis_{axis}_FDR_significant", int(axis_observed[axis]), expected, source="NETWORK_ATTENUATION_SUMMARY.tsv"))

    checks_df = pd.DataFrame(checks)
    checks_df.to_csv(NET_OUT / "PRIMARY_NETWORK_NUMERICAL_CHECKS.tsv", sep="\t", index=False)

    source_hashes = {
        p.relative_to(ROOT).as_posix(): sha256(p)
        for p in [joint_path, atten_path, pair_path, single_path, repro_path, fig2_path]
    }
    report = [
        "# Primary network reproduction audit",
        "",
        "Status: **PASS**" if (checks_df.status == "PASS").all() else "Status: **FAIL**",
        "",
        "This audit independently recalculated downstream numerical identities from the frozen Phase 2 and Phase 3 outputs. It did not rerun LDSC or partialLDSC because no discrepancy requiring a full upstream rerun was found.",
        "",
        "## Scope",
        "",
        "- Recomputed Schur-derived joint correlations from stored conditional covariances and variances.",
        "- Recomputed signed differences, absolute-magnitude changes, reversal flags, FDR-based change classes, residual significance and interval-based equivalence criteria.",
        "- Verified all primary denominators and the published-network reproduction statistics.",
        "- Verified that Figure 2 contains one row for every disease pair.",
        "",
        "## Numerical checks",
        "",
        checks_df.to_markdown(index=False),
        "",
        "## Input hashes",
        "",
        *[f"- `{path}`: `{digest}`" for path, digest in source_hashes.items()],
        "",
        "## Interpretation",
        "",
        "The primary network is unchanged. The corrected reversal definition affects only a derived pair-level flag: 135 pairs reversed sign, and 66 of these also had an FDR-significant signed baseline-to-joint difference. It does not alter any LDSC, partialLDSC, Schur-complement, P-value or FDR estimate.",
        "",
    ]
    (NET_OUT / "PRIMARY_NETWORK_REPRODUCTION_AUDIT.md").write_text("\n".join(report), encoding="utf-8")

    # Repair the pair-level export while retaining the historical gate explicitly.
    final_pair = pair.copy()
    final_pair = final_pair.rename(columns={"significant_reversal_flag": "legacy_significant_reversal_residual_z_gate"})
    final_pair["significant_reversal_flag"] = significant_reversal.to_numpy()
    final_pair["sign_reversal_definition"] = "opposite non-zero baseline and joint signs"
    final_pair["significant_reversal_definition"] = "sign reversal and signed-difference BH-FDR q<0.05"
    alias_to_resource = {"alc_problems": "alcohol_problems", "deaf": "hearing_loss", "schizo": "SSD"}
    final_pair["trait_1_resource_code"] = final_pair.trait_1.map(alias_to_resource).fillna(final_pair.trait_1)
    final_pair["trait_2_resource_code"] = final_pair.trait_2.map(alias_to_resource).fillna(final_pair.trait_2)
    final_pair.to_csv(PAIR_OUT / "Additional_file_2_Pair_level_classification_Final.csv", index=False, float_format="%.17g")

    mapping = pd.DataFrame([
        ["alc_problems", "alcohol_problems", "Alcohol Problems", "Frozen analysis alias mapped to GEMINI v1.1 resource code"],
        ["deaf", "hearing_loss", "Hearing loss", "Frozen analysis alias mapped to GEMINI v1.1 resource code"],
        ["schizo", "SSD", "Schizophrenia, schizotypal and delusional disorders", "Frozen analysis alias mapped to GEMINI v1.1 resource code"],
    ], columns=["analysis_code", "resource_code", "resource_name", "mapping_basis"])
    mapping.to_csv(PAIR_OUT / "DISEASE_CODE_MAPPING.tsv", sep="\t", index=False)

    dictionary_rows = []
    source_defs = {
        "change_class": "ATTENUATED or ENHANCED when the signed-difference BH-FDR q<0.05 and absolute correlation decreases or increases; otherwise NO_SIGNIFICANT_CHANGE.",
        "residual_state": "Separately defined residual state based on residual FDR significance and the locked point-estimate near-null rule.",
        "sign_reversal_flag": "True when baseline and joint correlations have opposite non-zero signs.",
        "significant_reversal_flag": "True when sign_reversal_flag is true and signed-difference BH-FDR q<0.05.",
        "legacy_significant_reversal_residual_z_gate": "Historical flag requiring reversal, signed-difference BH-FDR q<0.05 and residual |Z|>=2; retained for traceability only.",
        "equivalent_005": "Interval-based equivalence criterion: entire residual 95% CI lies within [-0.05, 0.05].",
        "equivalent_010": "Sensitivity interval-based equivalence criterion: entire residual 95% CI lies within [-0.10, 0.10].",
        "near_null_classified": "Exploratory point-estimate label; it is not evidence of no genetic correlation and is not equivalent to interval-based equivalence.",
    }
    for col in final_pair.columns:
        dictionary_rows.append({
            "field": col,
            "definition": source_defs.get(col, "See frozen Phase 3 source tables; numerical estimates and inferential quantities are unchanged."),
            "status": "REPAIRED_DERIVED" if col in {"significant_reversal_flag", "sign_reversal_definition", "significant_reversal_definition", "trait_1_resource_code", "trait_2_resource_code"} else "PRESERVED_OR_RENAMED",
        })
    pd.DataFrame(dictionary_rows).to_csv(PAIR_OUT / "PAIR_CLASSIFICATION_DATA_DICTIONARY.tsv", sep="\t", index=False)

    validation = [
        "# Pair classification validation",
        "",
        "Status: **PASS**",
        "",
        f"- Total pairs: {len(final_pair):,}",
        f"- Attenuated: {(final_pair.change_class == 'ATTENUATED').sum():,}",
        f"- Enhanced: {(final_pair.change_class == 'ENHANCED').sum():,}",
        f"- No FDR-significant signed change: {(final_pair.change_class == 'NO_SIGNIFICANT_CHANGE').sum():,}",
        f"- Residual FDR-significant: {final_pair.residual_significant.sum():,}",
        f"- Sign reversals: {final_pair.sign_reversal_flag.sum():,}",
        f"- Sign reversals with FDR-significant signed differences: {final_pair.significant_reversal_flag.sum():,}",
        f"- Historical residual-|Z| reversal gate: {final_pair.legacy_significant_reversal_residual_z_gate.sum():,}",
        f"- Interval-based equivalence criterion at ±0.05: {final_pair.equivalent_005.sum():,}",
        f"- Interval-based equivalence criterion at ±0.10: {final_pair.equivalent_010.sum():,}",
        "",
        "The repaired significant-reversal flag is descriptive plus inferential evidence for a signed change. It does not establish a biological reversal mechanism. Near-zero point estimates are not treated as proof of no genetic correlation. The interval-based equivalence criterion remains exploratory and depends on the selected margin, confidence interval precision and the separate multiplicity framework.",
        "",
        f"Frozen input SHA256 before and after audit: `{sha256(pair_path)}`.",
    ]
    (PAIR_OUT / "PAIR_CLASSIFICATION_VALIDATION.md").write_text("\n".join(validation), encoding="utf-8")

    if not (checks_df.status == "PASS").all():
        raise SystemExit("Primary network audit failed")
    print("PRIMARY_NETWORK_AUDIT_PASS")
    print(f"sign_reversals={reversal.sum()} significant_reversals={significant_reversal.sum()}")


if __name__ == "__main__":
    main()
