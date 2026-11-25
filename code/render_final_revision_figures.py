#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import os
import shutil
import subprocess
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(os.environ.get("MODRISK_ROOT", Path.cwd())).resolve()
HIST = ROOT / "final_major_revision"
OUT = ROOT / "final_revision_closeout/05_manuscript_revision/figures"
OUT.mkdir(parents=True, exist_ok=True)

spec = importlib.util.spec_from_file_location(
    "previous_renderer", ROOT / "revision_round_20261124/analysis/render_revision_figures.py"
)
previous = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(previous)
previous.OUT = OUT

BLUE = "#4477AA"
ORANGE = "#EE7733"
GREY = "#7F8B99"
LIGHT = "#D8DEE6"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7.5,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6,
    "legend.fontsize": 7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.bbox": "tight",
})


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.svg")
    fig.savefig(OUT / f"{stem}.png", dpi=600)
    fig.savefig(OUT / f"{stem}.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def figure3(d):
    colours = {"ATTENUATED": BLUE, "ENHANCED": ORANGE, "NO_SIGNIFICANT_CHANGE": GREY}
    labels = {
        "ATTENUATED": "Attenuated",
        "ENHANCED": "Enhanced",
        "NO_SIGNIFICANT_CHANGE": "No FDR-significant\nsigned change",
    }
    fig, ax = plt.subplots(figsize=(8.8, 5.3))
    fig.subplots_adjust(left=.10, right=.68, top=.90, bottom=.12)
    for cls in ("NO_SIGNIFICANT_CHANGE", "ATTENUATED", "ENHANCED"):
        for residual in (False, True):
            g = d[(d.change_class == cls) & (d.residual_significant == residual)]
            ax.scatter(
                g.baseline_abs_rg, g.joint_abs_rg, s=14 if residual else 12,
                facecolors=colours[cls] if residual else "none",
                edgecolors="white" if residual else colours[cls],
                linewidths=.35 if residual else .65, alpha=.74 if residual else .62,
                rasterized=True,
            )
    limit = max(d.baseline_abs_rg.max(), d.joint_abs_rg.max()) * 1.04
    ax.plot([0, limit], [0, limit], "--", lw=.8, color="#33485A")
    ax.set(xlim=(0, limit), ylim=(0, limit), xlabel="Baseline |r$_g$|",
           ylabel="Joint residual |r$_g$|",
           title="Baseline versus jointly conditioned genetic-correlation magnitude")
    class_handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=5.8,
               markerfacecolor=colours[c], markeredgecolor="white",
               label=f"{labels[c]} (n={sum(d.change_class == c):,})")
        for c in ("ATTENUATED", "ENHANCED", "NO_SIGNIFICANT_CHANGE")
    ]
    residual_handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=5.8,
               markerfacecolor="#4B5563", markeredgecolor="white", label="Residual FDR < 0.05"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=5.8,
               markerfacecolor="none", markeredgecolor="#4B5563", label="Residual FDR ≥ 0.05"),
    ]
    first = ax.legend(handles=class_handles, title="Change class", loc="upper left",
                      bbox_to_anchor=(1.02, 1), frameon=False)
    ax.add_artist(first)
    ax.legend(handles=residual_handles, title="Residual significance (marker fill)",
              loc="upper left", bbox_to_anchor=(1.02, .66), frameon=False)
    ax.text(.98, .03, "2,485 pairs\n1,019 FDR-significant signed differences\n"
            "1,100 residual FDR-significant pairs\n135 sign reversals; 66 with significant signed change",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7,
            bbox=dict(facecolor="white", edgecolor=LIGHT, boxstyle="round,pad=.3"))
    save(fig, "Figure3_FINAL")


def figure4(d):
    registry = pd.read_csv(ROOT / "config/diseases_gemini.tsv", sep="\t")
    names = dict(zip(registry.file_prefix, registry.name))
    names.update({
        "COPD": "Chronic obstructive pulmonary disease", "GORD": "Gastro-oesophageal reflux disease",
        "IDA": "Iron deficiency anaemia", "ED": "Erectile dysfunction", "RA": "Rheumatoid arthritis",
        "HF": "Heart failure", "OA": "Osteoarthritis", "IBS": "Irritable bowel syndrome",
        "alc_problems": "Alcohol problems", "periph_neuro": "Peripheral neuropathy",
        "coronary_heart": "Coronary heart disease", "urine_incont": "Urinary incontinence",
    })
    diseases = sorted(set(d.trait_1) | set(d.trait_2))
    rows = []
    for disease in diseases:
        x = d[(d.trait_1 == disease) | (d.trait_2 == disease)]
        rows.append({
            "disease": disease,
            "display_name": names.get(disease, disease.replace("_", " ")),
            "median_absolute_magnitude_change": x.absolute_attenuation.abs().median(),
            "residual_significant_edges": int(x.residual_significant.sum()),
        })
    source = pd.DataFrame(rows)
    source.to_csv(OUT / "Figure4_FINAL_SOURCE.tsv", sep="\t", index=False)
    a = source.nlargest(20, "median_absolute_magnitude_change").copy()
    b = source.nlargest(20, "residual_significant_edges").copy()
    common = set(a.disease) & set(b.disease)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.5), gridspec_kw={"wspace": 1.05})
    panels = (
        (axes[0], a, "median_absolute_magnitude_change", "Median |Δ|r$_g$|| across 70 incident edges",
         "Top 20 by median absolute magnitude of change", "A"),
        (axes[1], b, "residual_significant_edges", "Number of residual FDR-significant edges",
         "Top 20 by residual connectivity", "B"),
    )
    for ax, frame, field, xlabel, title, panel in panels:
        frame = frame.sort_values(field)
        ax.barh(frame.display_name, frame[field],
                color=[ORANGE if x in common else BLUE for x in frame.disease], height=.75)
        ax.set_xlabel(xlabel); ax.set_title(title, loc="left")
        ax.text(-.10, 1.04, panel, transform=ax.transAxes, fontsize=10, fontweight="bold")
        ax.grid(axis="x", color=LIGHT, lw=.5); ax.set_axisbelow(True)
    fig.text(.5, .005, "Top 20 is a descriptive display, not a statistical threshold. Orange identifies diseases in both lists. "
             "Panel A uses unsigned magnitude; Figure 2 retains direction.", ha="center", fontsize=7)
    save(fig, "Figure4_FINAL")


def figure_s2():
    ivw = pd.read_csv(HIST / "05_robust_mvmr/MVMR_IVW_FINAL.tsv", sep="\t")
    med = pd.read_csv(HIST / "05_robust_mvmr/MVMR_MEDIAN_FULL_RESULTS.tsv", sep="\t")
    qh = pd.read_csv(HIST / "05_robust_mvmr/MVMR_WEAK_IV_AWARE_RESULTS.tsv", sep="\t")
    labels = {
        "COPD": "COPD", "HF": "Heart failure", "OA": "Osteoarthritis",
        "alcohol_problems": "Alcohol problems", "carpal_tunnel": "Carpal tunnel syndrome",
        "coronary_heart": "Coronary heart disease", "diabetes_type_2": "Type 2 diabetes",
        "fibromyalgia": "Fibromyalgia", "gastritis_duodenitis": "Gastritis/duodenitis",
        "gout": "Gout", "hypertension": "Hypertension", "insomnia": "Insomnia",
        "periph_neuro": "Peripheral neuropathy", "sleep_apnoea": "Sleep apnoea",
        "valve_nonrheumatic": "Non-rheumatic valve disease",
    }
    order = ivw.sort_values("SMK_beta").outcome.tolist()
    datasets = {"MVMR-IVW": ivw, "Weighted LAD median": med, "Q-minimisation": qh}
    colours = {"MVMR-IVW": BLUE, "Weighted LAD median": ORANGE, "Q-minimisation": "#228833"}
    offsets = {"MVMR-IVW": -.20, "Weighted LAD median": 0, "Q-minimisation": .20}
    y = np.arange(len(order))
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 6.3), sharey=True, constrained_layout=True)
    for ax, exposure, panel in ((axes[0], "BMI", "A"), (axes[1], "SMK", "B")):
        ax.axvline(0, color="#4B5563", lw=.8, ls="--")
        for method, df in datasets.items():
            x = df.set_index("outcome").loc[order]
            beta, se = x[f"{exposure}_beta"].to_numpy(), x[f"{exposure}_SE"].to_numpy()
            ax.errorbar(beta, y + offsets[method], xerr=1.96 * se, fmt="o", ms=3.4,
                        mfc=colours[method], mec=colours[method], ecolor=colours[method],
                        lw=.7, capsize=1.2, label=method)
        ax.set_xlabel("Beta (disease log-odds scale)")
        ax.set_title(f"{panel}  {'BMI' if exposure == 'BMI' else 'Smoking initiation'}", loc="left", fontweight="bold")
    axes[0].set_yticks(y, [labels[x] for x in order]); axes[0].set_ylabel("Outcome")
    axes[1].legend(frameon=False, loc="lower right")
    fig.suptitle("Exploratory BMI–smoking MVMR sensitivity estimates", x=.01, ha="left", fontweight="bold")
    fig.text(.01, -.01, "Points show estimates and 95% confidence intervals. Conditional F: BMI 25.48; smoking initiation 3.85. "
             "Q-minimisation uses the verified non-zero exposure-association error correlation; all estimates remain exploratory.",
             ha="left", fontsize=7)
    save(fig, "FigureS2_FINAL")


def figure_s3():
    """Rebuild the MR-AHC display without mislabelling the historical 54-variant set."""
    d = pd.read_csv(ROOT / "results/final/FIGURE6_VARIANT_PROFILES.tsv", sep="\t")
    matrix = d[["COPD", "CHD", "T2D", "OA"]].to_numpy()
    order = np.argsort(matrix[:, 0])
    matrix = matrix[order, :]
    fig = plt.figure(figsize=(11.5, 7.3), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[2.2, .85, 1.15])
    ax1, ax2, ax3 = (fig.add_subplot(gs[0, i]) for i in range(3))

    vmax = max(8.0, float(np.nanquantile(np.abs(matrix), .98)))
    im = ax1.imshow(matrix, aspect="auto", cmap="PuOr_r",
                    norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), rasterized=True)
    ax1.set_xticks(range(4), ["COPD", "Coronary\nheart disease", "Type 2\ndiabetes", "Osteoarthritis"])
    ax1.set_yticks([])
    ax1.set_ylabel("54 MR-AHC profile variants")
    ax1.set_title("A  Variant-level outcome profiles", loc="left", fontweight="bold")
    cb = fig.colorbar(im, ax=ax1, orientation="horizontal", fraction=.06, pad=.08)
    cb.set_label("Variant-specific ratio estimate")

    values = [0.0, 0.064, 0.027]
    ax2.bar([0, 1, 2], values, color=[GREY, BLUE, "#228833"], width=.8)
    ax2.set_xticks([0, 1, 2], ["Median", "Mean", ">0.8"])
    ax2.set_ylabel("Bootstrap co-clustering summary")
    ax2.set_ylim(0, .11)
    ax2.set_title("B  Stability summary", loc="left", fontweight="bold")
    ax2.text(1, .099, "No stable\nco-clustering", ha="center", va="top", color="#213547", fontsize=8)

    ax3.axis("off")
    ax3.set_title("C  Cluster-selection result", loc="left", fontweight="bold")
    ax3.text(.5, .62, "No robust\nmulti-cluster\nsolution", transform=ax3.transAxes,
             ha="center", va="center", fontsize=13, fontweight="bold", color="#213547",
             bbox=dict(boxstyle="round,pad=.8", facecolor="#FFF4E8", edgecolor=ORANGE, linewidth=1.2))
    ax3.text(.5, .25, "Published Q-based rule:\nall merged clusters failed.\n\n"
             "Descriptive step-40 profiles\nwere retained for audit only.",
             transform=ax3.transAxes, ha="center", va="center", fontsize=8)
    save(fig, "FigureS3_FINAL")


def figure1():
    """Update the workflow label while retaining the previously approved design."""
    src = ROOT / "revision_round_20261124/figures/Figure1_REVISED.svg"
    dst = OUT / "Figure1_FINAL.svg"
    text = src.read_text()
    text = text.replace(">Joint changes</text>", ">Signed differences</text>")
    text = text.replace(">1,019 pairs</text>", ">1,019 FDR-significant</text>")
    dst.write_text(text)
    try:
        import cairosvg
        cairosvg.svg2png(bytestring=text.encode(), write_to=str(OUT / "Figure1_FINAL.png"), output_width=3200)
    except Exception:
        subprocess.run(["rsvg-convert", "-w", "3200", "-o", str(OUT / "Figure1_FINAL.png"), str(dst)], check=True)
    with Image.open(OUT / "Figure1_FINAL.png") as im:
        im.save(OUT / "Figure1_FINAL.tiff", dpi=(600, 600), compression="tiff_lzw")


def main():
    pair_heat = pd.read_csv(ROOT / "results/phase3_joint/FIGURE2_HEATMAP_DATA.tsv", sep="\t")
    final_pair = pd.read_csv(HIST / "07_figures_tables/Additional_file_2_Pair_level_classification_Final.csv")
    previous.figure2(pair_heat, final_pair)
    for ext in ("png", "svg", "tiff"):
        src = OUT / f"Figure2_REVISED.{ext}"
        if src.exists():
            src.replace(OUT / f"Figure2_FINAL.{ext}")
    previous.figure_s1(pair_heat)
    for ext in ("png", "svg", "tiff"):
        src = OUT / f"FigureS1_REVISED.{ext}"
        if src.exists():
            src.replace(OUT / f"FigureS1_FINAL.{ext}")
    figure3(final_pair)
    figure4(final_pair)
    figure_s2()
    figure1()
    print(OUT)


if __name__ == "__main__":
    main()
