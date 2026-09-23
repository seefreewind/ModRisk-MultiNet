#!/usr/bin/env python3
"""Stack the frozen R2 Figure 4 panels at journal page width; data unchanged."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures/bmc_genomics_submission"
DATA = ROOT / "figures/bmc_genomics_r2/FIGURE4_R2_DISEASE_SUMMARY.tsv"

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42, "svg.fonttype": "none", "ps.fonttype": 42,
    "font.size": 6.5, "axes.titlesize": 8, "axes.labelsize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.bbox": "tight",
})


def main():
    d = pd.read_csv(DATA, sep="\t")
    assert len(d) == 71 and d.residual_significant_edges_r2.sum() == 2200
    a = d.nlargest(20, "median_absolute_magnitude_change").copy()
    b = d.nlargest(20, "residual_significant_edges_r2").copy()
    common = set(a.disease) & set(b.disease)
    fig, axes = plt.subplots(2, 1, figsize=(6.3, 8.4), gridspec_kw={"hspace": 0.35})
    for ax, frame, field, xlabel, title, panel in (
        (axes[0], a, "median_absolute_magnitude_change", "Median absolute change in |rg|", "Largest magnitude changes", "A"),
        (axes[1], b, "residual_significant_edges_r2", "Residual FDR-significant edges", "Greatest residual connectivity", "B"),
    ):
        frame = frame.sort_values(field)
        ax.barh(frame.display_name, frame[field], color=["#EE7733" if x in common else "#4477AA" for x in frame.disease], height=.72)
        ax.set_xlabel(xlabel)
        ax.set_title(f"{panel}  {title}", loc="left", fontweight="bold")
        ax.grid(axis="x", color="#D8DEE6", lw=.5)
        ax.set_axisbelow(True)
    fig.subplots_adjust(left=.47, right=.97, top=.96, bottom=.07)
    fig.text(.5, .015, "Orange: disease appears in both top-20 lists. Panel A includes increases and decreases in |rg|.", ha="center", fontsize=6)
    fig.savefig(OUT / "Figure4_SUBMISSION.pdf")
    fig.savefig(OUT / "Figure4_SUBMISSION.svg")
    fig.savefig(OUT / "Figure4_SUBMISSION.png", dpi=600)
    plt.close(fig)
    print(OUT / "Figure4_SUBMISSION.pdf")


if __name__ == "__main__":
    main()
