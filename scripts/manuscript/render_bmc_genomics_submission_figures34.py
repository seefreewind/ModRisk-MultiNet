#!/usr/bin/env python3
"""Submission wording export of R2 Figures 3/4 from frozen estimates and labels."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures" / "bmc_genomics_submission"
OUT.mkdir(parents=True, exist_ok=True)
BLUE = "#4477AA"
ORANGE = "#EE7733"
GREY = "#8A96A3"
LIGHT = "#D8DEE6"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.bbox": "tight",
})


def save(fig: plt.Figure, name: str) -> None:
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.svg")
    fig.savefig(OUT / f"{name}.png", dpi=600)
    plt.close(fig)


def figure3(d: pd.DataFrame) -> None:
    colors = {"ATTENUATED": BLUE, "ENHANCED": ORANGE, "NO_SIGNIFICANT_CHANGE": GREY}
    labels = {"ATTENUATED": "Attenuated", "ENHANCED": "Enhanced", "NO_SIGNIFICANT_CHANGE": "No FDR-significant signed change"}
    fig, ax = plt.subplots(figsize=(6.0, 5.2), constrained_layout=True)
    for cls in ("NO_SIGNIFICANT_CHANGE", "ATTENUATED", "ENHANCED"):
        for residual in (False, True):
            g = d[(d.change_class == cls) & (d.residual_significant == residual)]
            ax.scatter(
                g.baseline_abs_rg, g.joint_abs_rg, s=13 if residual else 11,
                facecolors=colors[cls] if residual else "none",
                edgecolors=colors[cls] if not residual else "white",
                linewidths=0.38 if residual else 0.6,
                alpha=0.74 if residual else 0.58,
                rasterized=True, zorder=2 if residual else 1,
            )
    limit = max(d.baseline_abs_rg.max(), d.joint_abs_rg.max()) * 1.04
    ax.plot([0, limit], [0, limit], ls="--", lw=0.8, color="#33485A", zorder=0)
    ax.set(xlim=(0, limit), ylim=(0, limit), xlabel="Baseline |rg|", ylabel="Joint residual |rg|",
           title="Baseline versus jointly conditioned genetic correlation")
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker="o", linestyle="none", markersize=5.5,
                      markerfacecolor=colors[c], markeredgecolor="white", label=f"{labels[c]} (n={sum(d.change_class == c)})")
               for c in ("ATTENUATED", "ENHANCED", "NO_SIGNIFICANT_CHANGE")]
    handles.extend([
        Line2D([0], [0], marker="o", linestyle="none", markersize=5.5,
               markerfacecolor="#444444", markeredgecolor="white", label="Residual FDR < 0.05"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=5.5,
               markerfacecolor="none", markeredgecolor="#444444", label="Residual FDR ≥ 0.05"),
    ])
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False)
    ax.text(0.98, 0.03, "2,485 pairs; 1,019 significant signed differences\n1,100 residual FDR-significant pairs",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7,
            bbox=dict(facecolor="white", edgecolor=LIGHT, boxstyle="round,pad=.3"))
    save(fig, "Figure3_SUBMISSION")


def figure4(d: pd.DataFrame) -> None:
    source = pd.read_csv(ROOT / "results/phase3_joint/DISEASE_LEVEL_JOINT_SUMMARY.tsv", sep="\t")
    registry = pd.read_csv(ROOT / "config/diseases_gemini.tsv", sep="\t")
    names = dict(zip(registry.file_prefix, registry.name))
    names.update({"COPD": "Chronic obstructive pulmonary disease", "GORD": "Gastro-oesophageal reflux disease",
                  "IDA": "Iron deficiency anaemia", "ED": "Erectile dysfunction",
                  "RA": "Rheumatoid arthritis", "HF": "Heart failure", "OA": "Osteoarthritis", "IBS": "Irritable bowel syndrome",
                  "alc_problems": "Alcohol problems", "periph_neuro": "Peripheral neuropathy",
                  "coronary_heart": "Coronary heart disease", "urine_incont": "Urinary incontinence"})
    source["display_name"] = source.disease.map(names).fillna(source.disease.str.replace("_", " "))
    incident = pd.concat([d[["trait_1", "residual_significant"]].rename(columns={"trait_1": "disease"}),
                          d[["trait_2", "residual_significant"]].rename(columns={"trait_2": "disease"})])
    counts = incident.groupby("disease").residual_significant.sum().astype(int)
    source["residual_significant_edges_r2"] = source.disease.map(counts)
    assert source.residual_significant_edges_r2.sum() == 2200
    assert (source.residual_significant_edges_r2 <= 70).all()
    source = source.rename(columns={"median_joint_attenuation": "median_absolute_magnitude_change"})
    source[["disease", "display_name", "n_edges", "median_absolute_magnitude_change", "residual_significant_edges_r2"]].to_csv(
        OUT / "FIGURE4_R2_DISEASE_SUMMARY.tsv", sep="\t", index=False)
    a = source.nlargest(20, "median_absolute_magnitude_change").copy()
    b = source.nlargest(20, "residual_significant_edges_r2").copy()
    common = set(a.disease) & set(b.disease)
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 6.4), gridspec_kw={"wspace": 1.05})
    for ax, frame, field, xlabel, title, panel in (
        (axes[0], a, "median_absolute_magnitude_change", "Median absolute change in |rg|", "Largest magnitude changes", "A"),
        (axes[1], b, "residual_significant_edges_r2", "Residual FDR-significant edges", "Greatest residual connectivity", "B"),
    ):
        frame = frame.sort_values(field)
        ax.barh(frame.display_name, frame[field],
                color=[ORANGE if x in common else BLUE for x in frame.disease], height=0.75)
        ax.set_xlabel(xlabel)
        ax.set_title(title, loc="left")
        ax.text(-0.10, 1.04, panel, transform=ax.transAxes, fontsize=10, fontweight="bold")
        ax.grid(axis="x", color=LIGHT, lw=0.5)
        ax.set_axisbelow(True)
    fig.text(0.5, 0.01, "Orange: disease appears in both top-20 lists. Panel A is unsigned and includes increases and decreases in |rg|.",
             ha="center", fontsize=7)
    save(fig, "Figure4_SUBMISSION")


if __name__ == "__main__":
    d = pd.read_csv(ROOT / "results/phase3_classification_repair_R2/PAIR_CLASSIFICATION_REBUILT.tsv", sep="\t")
    assert len(d) == 2485 and d.residual_significant.sum() == 1100
    figure3(d)
    figure4(d)
    print(OUT)
