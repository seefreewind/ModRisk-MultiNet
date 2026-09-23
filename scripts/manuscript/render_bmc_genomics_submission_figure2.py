#!/usr/bin/env python3
"""Correct Figure 2 terminology without changing its frozen values or order."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures" / "bmc_genomics_submission"
OUT.mkdir(parents=True, exist_ok=True)
AXES = ["BMI", "SMK", "SBP", "ApoB", "FG", "SLEEP", "JOINT"]
EXPECTED = {"BMI": 1193, "SMK": 1178, "SBP": 120, "ApoB": 7, "FG": 5, "SLEEP": 148, "JOINT": 1019}

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 6.5, "axes.titlesize": 8, "axes.labelsize": 6.5,
    "xtick.labelsize": 6, "ytick.labelsize": 5.1,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.bbox": "tight",
})


def main() -> None:
    pair = pd.read_csv(ROOT / "results/phase3_joint/FIGURE2_HEATMAP_DATA.tsv", sep="\t")
    rebuilt = pd.read_csv(ROOT / "results/phase3_classification_repair_R2/PAIR_CLASSIFICATION_REBUILT.tsv", sep="\t")
    old = pd.read_csv(ROOT / "results/phase3_joint/DISEASE_LEVEL_JOINT_SUMMARY.tsv", sep="\t")
    label_table = pd.read_csv(ROOT / "config/FIGURE_DISEASE_LABEL_MAP.tsv", sep="\t")
    names = dict(zip(label_table.internal_disease_id, label_table.publication_label))
    assert len(pair) == len(rebuilt) == 2485 and len(old) == 71

    summary = {}
    for axis in AXES:
        signed_magnitude_change = pair.baseline.abs() - pair[axis].abs()
        long = pd.concat([
            pd.DataFrame({"disease": pair.disease1, "change": signed_magnitude_change}),
            pd.DataFrame({"disease": pair.disease2, "change": signed_magnitude_change}),
        ], ignore_index=True)
        summary[axis] = long.groupby("disease").change.median()
    summary = pd.DataFrame(summary)

    significant = rebuilt[rebuilt.difference_fdr < .05]
    count = pd.concat([significant.trait_1, significant.trait_2]).value_counts()
    source = old.set_index("disease")
    assert (source.loc[count.index, "significantly_attenuated_edges"].to_numpy() == count.to_numpy()).all()
    # The historical rank orders by significant signed-difference count, then
    # by the median *absolute* magnitude change. Preserve this exact order.
    expected_order = old.sort_values(
        ["significantly_attenuated_edges", "median_joint_attenuation", "disease"],
        ascending=[False, False, True],
    ).disease.tolist()
    order = old.sort_values("rank_attenuation_burden").disease.tolist()
    assert order == expected_order
    summary = summary.loc[order]
    assert summary.notna().all().all()

    axis_counts = pd.read_csv(ROOT / "results/phase2_single_factor/NETWORK_ATTENUATION_SUMMARY.tsv", sep="\t")
    counts = dict(zip(axis_counts.exposure, axis_counts.n_significant_q005))
    counts["JOINT"] = int(len(significant))
    assert counts == EXPECTED

    data = pd.DataFrame({
        "disease": order,
        "publication_label": [names[x] for x in order],
        "signed_difference_significant_edge_count": source.loc[order, "significantly_attenuated_edges"].to_numpy(),
        "median_absolute_magnitude_change_tiebreak": source.loc[order, "median_joint_attenuation"].to_numpy(),
    })
    for axis in AXES:
        data[f"median_signed_magnitude_change_{axis}"] = summary[axis].to_numpy()
    data.to_csv(OUT / "Figure2_SUBMISSION_DISEASE_SUMMARY.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(5.4, 8.1))
    matrix = summary[AXES].to_numpy()
    vmax = max(float(np.nanmax(np.abs(matrix))), 0.12)
    image = ax.imshow(matrix, aspect="auto", cmap="PuOr_r",
                      norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax))
    ax.set_xticks(range(7), ["BMI", "Smoking", "SBP", "ApoB", "Fasting\nglucose", "Sleep\nduration", "Joint"])
    ax.set_yticks(range(71), data.publication_label, fontsize=5.1)
    ax.set_xlabel("Conditioning axis")
    ax.set_ylabel("Disease, ordered by significant signed-difference count")
    ax.set_title("Risk-factor-associated changes in genetic-correlation magnitude across diseases",
                 loc="left", fontweight="bold", pad=31)
    for i, axis in enumerate(AXES):
        ax.text(i, -2.3, f"{counts[axis]:,}", ha="center", va="bottom", fontsize=5.6, fontweight="bold")
    ax.text(6.48, -3.8, "FDR-significant signed baseline-to-conditioned differences",
            ha="right", va="bottom", fontsize=5.6)
    bar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    bar.set_label("Median [|baseline rg| − |conditioned rg|] across disease edges")
    fig.savefig(OUT / "Figure2_SUBMISSION.pdf")
    fig.savefig(OUT / "Figure2_SUBMISSION.svg")
    fig.savefig(OUT / "Figure2_SUBMISSION.png", dpi=600)
    fig.savefig(OUT / "Figure2_SUBMISSION.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(OUT / "Figure2_SUBMISSION.pdf")


if __name__ == "__main__":
    main()
