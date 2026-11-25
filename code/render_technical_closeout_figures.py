#!/usr/bin/env python3
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("MODRISK_ROOT", Path.cwd())).resolve()
SOURCE = ROOT / "final_revision_closeout/05_manuscript_revision/figures"
OUT = ROOT / "pre_submission_technical_closeout/05_document_corrections/figures"
OUT.mkdir(parents=True, exist_ok=True)

BLUE = "#4477AA"
ORANGE = "#EE7733"
GREEN = "#228833"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
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


def copy_unchanged() -> None:
    for stem in ["Figure2_FINAL", "Figure3_FINAL", "Figure4_FINAL", "FigureS1_FINAL"]:
        for ext in ["png", "svg", "tiff"]:
            shutil.copy2(SOURCE / f"{stem}.{ext}", OUT / f"{stem}.{ext}")


def figure1() -> None:
    src = SOURCE / "Figure1_FINAL.svg"
    text = src.read_text(encoding="utf-8")
    text = text.replace(
        '<text style="font-size: 5.7px; font-family: \'Arial\', \'Helvetica\', \'DejaVu Sans\', sans-serif; fill: #111111" transform="translate(252.000000 85.731375)">Signed differences</text>\n'
        '    <text style="font-size: 5.7px; font-family: \'Arial\', \'Helvetica\', \'DejaVu Sans\', sans-serif; fill: #111111" transform="translate(256.67137 91.199367)">1,019 FDR-significant</text>',
        '<text style="font-size: 5.2px; font-family: \'Arial\', \'Helvetica\', \'DejaVu Sans\', sans-serif; fill: #111111" transform="translate(246.900000 85.731375)">1,019 FDR-significant</text>\n'
        '    <text style="font-size: 5.2px; font-family: \'Arial\', \'Helvetica\', \'DejaVu Sans\', sans-serif; fill: #111111" transform="translate(252.100000 91.199367)">signed differences</text>',
    )
    if "1,019 FDR-significant</text>" not in text or ">signed differences</text>" not in text:
        raise RuntimeError("Figure 1 text replacement failed")
    dst = OUT / "Figure1_FINAL.svg"
    dst.write_text(text, encoding="utf-8")
    cairosvg = ROOT / "final_revision_closeout/.venv/bin/cairosvg"
    subprocess.run([str(cairosvg), str(dst), "-o", str(OUT / "Figure1_FINAL.png"), "-s", "4"], check=True)
    from PIL import Image
    with Image.open(OUT / "Figure1_FINAL.png") as image:
        image.save(OUT / "Figure1_FINAL.tiff", dpi=(600, 600), compression="tiff_lzw")


def figure_s2() -> None:
    data = pd.read_csv(
        ROOT / "pre_submission_technical_closeout/02_ci_validation/Additional_file_5_MVMR_Full_Results_CI_corrected.tsv",
        sep="\t",
    )
    data = data[data["method"].isin(["MVMR-IVW", "Multivariable weighted LAD median", "Q-minimisation"])].copy()
    labels = {
        "COPD": "COPD", "HF": "Heart failure", "OA": "Osteoarthritis",
        "alcohol_problems": "Alcohol problems", "carpal_tunnel": "Carpal tunnel syndrome",
        "coronary_heart": "Coronary heart disease", "diabetes_type_2": "Type 2 diabetes",
        "fibromyalgia": "Fibromyalgia", "gastritis_duodenitis": "Gastritis/duodenitis",
        "gout": "Gout", "hypertension": "Hypertension", "insomnia": "Insomnia",
        "periph_neuro": "Peripheral neuropathy", "sleep_apnoea": "Sleep apnoea",
        "valve_nonrheumatic": "Non-rheumatic valve disease",
    }
    ivw_smk = data[(data.method == "MVMR-IVW") & (data.exposure == "Smoking initiation")]
    order = ivw_smk.sort_values("beta").outcome.tolist()
    display = {
        "MVMR-IVW": "MVMR-IVW",
        "Multivariable weighted LAD median": "Weighted LAD median",
        "Q-minimisation": "Q-minimisation",
    }
    colours = {"MVMR-IVW": BLUE, "Multivariable weighted LAD median": ORANGE, "Q-minimisation": GREEN}
    offsets = {"MVMR-IVW": -0.20, "Multivariable weighted LAD median": 0.0, "Q-minimisation": 0.20}
    y = np.arange(len(order))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 5.6), sharey=True, constrained_layout=True)
    for ax, exposure, panel in [(axes[0], "BMI", "A"), (axes[1], "Smoking initiation", "B")]:
        ax.axvline(0, color="#4B5563", lw=0.8, ls="--")
        for method in display:
            x = data[(data.method == method) & (data.exposure == exposure)].set_index("outcome").loc[order]
            beta = x.beta.to_numpy()
            xerr = np.vstack([beta - x.CI_lower.to_numpy(), x.CI_upper.to_numpy() - beta])
            ax.errorbar(
                beta, y + offsets[method], xerr=xerr, fmt="o", ms=3.4,
                mfc=colours[method], mec=colours[method], ecolor=colours[method],
                lw=0.7, capsize=1.2, label=display[method],
            )
        ax.set_xlabel("Beta (disease log-odds scale)")
        ax.set_title(f"{panel}  {exposure}", loc="left", fontweight="bold")
    axes[0].set_yticks(y, [labels[x] for x in order])
    axes[0].set_ylabel("Outcome")
    axes[1].legend(frameon=False, loc="lower right")
    fig.suptitle("Exploratory BMI-smoking MVMR sensitivity estimates", x=0.01, ha="left", fontweight="bold")
    fig.text(
        0.01, -0.01,
        "Points show estimates and method-specific 95% confidence intervals. MVMR-IVW intervals use Student t reference values; weighted LAD median and Q-minimisation use normal-reference intervals. Conditional F: BMI 25.48; smoking initiation 3.85.",
        ha="left", fontsize=7,
    )
    fig.savefig(OUT / "FigureS2_FINAL.svg")
    fig.savefig(OUT / "FigureS2_FINAL.pdf")
    fig.savefig(OUT / "FigureS2_FINAL.png", dpi=600)
    fig.savefig(OUT / "FigureS2_FINAL.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def main() -> None:
    copy_unchanged()
    figure1()
    figure_s2()
    print(OUT)


if __name__ == "__main__":
    main()
