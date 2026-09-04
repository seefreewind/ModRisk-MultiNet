#!/usr/bin/env python3
"""Render frozen ModRisk-MultiNet figures for the Genome Medicine draft.

This script only reads frozen result tables and creates visual summaries. It
does not rerun or alter any scientific analysis.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures" / "genome_medicine_v1"
OUT.mkdir(parents=True, exist_ok=True)
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "svg.fonttype": "none", "pdf.fonttype": 42})
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 10,
                     "axes.labelsize": 8, "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
                     "savefig.bbox": "tight", "axes.spines.top": False,
                     "axes.spines.right": False})
NAVY, BLUE, TEAL, ORANGE, GREY = "#18324B", "#2F6B9A", "#2A9D8F", "#E07A3F", "#7A8793"

def save(fig, number):
    # Publication exports: .svg and .pdf vectors; .tiff and .png rasters at dpi=600.
    fig.savefig(OUT / f"Figure{number}.pdf")
    fig.savefig(OUT / f"Figure{number}.svg")
    fig.savefig(OUT / f"Figure{number}.png", dpi=600)
    fig.savefig(OUT / f"Figure{number}.tiff", dpi=600)
    plt.close(fig)

def fig1():
    fig, ax = plt.subplots(figsize=(7.2, 3.2)); ax.axis("off")
    boxes = [(0.02, .44, .15, .22, "71 diseases\n2,485 pairs", "#E8F0F5"),
             (.22, .44, .17, .22, "Baseline LDSC\nshared genetic liability", "#E8F0F5"),
             (.44, .62, .17, .22, "Six single axes\nBMI · SMK · SBP\nApoB · FG · SLEEP", "#EAF4F1"),
             (.44, .26, .17, .22, "Six-factor joint model\nresidual network", "#EAF4F1"),
             (.68, .44, .14, .22, "12 pairs\ncausal triangulation", "#FFF0E5"),
             (.87, .44, .11, .22, "54 instruments\nMR-AHC", "#F0EDF7")]
    for x,y,w,h,t,c in boxes:
        ax.add_patch(plt.Rectangle((x,y),w,h,facecolor=c,edgecolor=NAVY,lw=1.1,transform=ax.transAxes))
        ax.text(x+w/2,y+h/2,t,ha="center",va="center",color=NAVY,transform=ax.transAxes,fontsize=8)
    for a,b in [((.17,.55),(.22,.55)),((.39,.55),(.44,.70)),((.39,.55),(.44,.38)),((.61,.70),(.68,.55)),((.61,.38),(.68,.55)),((.82,.55),(.87,.55))]:
        ax.annotate("", xy=b, xytext=a, xycoords=ax.transAxes, textcoords=ax.transAxes,
                    arrowprops=dict(arrowstyle="->",color=GREY,lw=1.2))
    ax.text(.5,.93,"Analysis sequence",ha="center",va="center",transform=ax.transAxes,color=NAVY,weight="bold",fontsize=11)
    ax.text(.5,.07,"Covariance decomposition  →  selected direct-effect testing  →  variant-level heterogeneity",ha="center",transform=ax.transAxes,color=GREY,fontsize=8)
    save(fig,1)

def fig2():
    d = pd.read_csv(ROOT/"results/phase3_joint/FIGURE2_HEATMAP_DATA.tsv", sep="\t")
    axes = ["BMI","SMK","SBP","ApoB","FG","SLEEP","JOINT"]
    rows=[]
    for col in axes:
        z = d.assign(value=np.abs(d.baseline)-np.abs(d[col]))
        rows.append(pd.concat([z[["disease1","value"]].rename(columns={"disease1":"disease"}),
                               z[["disease2","value"]].rename(columns={"disease2":"disease"})]))
    s = pd.concat(rows, keys=axes, names=["axis","row"]).reset_index()
    summary = s.groupby(["disease","axis"], sort=False).value.median().unstack()
    order = pd.read_csv(ROOT/"results/phase3_joint/DISEASE_LEVEL_JOINT_SUMMARY.tsv",sep="\t").sort_values("rank_attenuation_burden").disease
    summary = summary.reindex(order)
    summary.to_csv(OUT/"Figure2_main_summary.tsv",sep="\t")
    fig, ax = plt.subplots(figsize=(6.2, 10.0))
    im=ax.imshow(summary[axes].to_numpy(),aspect="auto",cmap="RdBu_r",norm=TwoSlopeNorm(0, vmin=-.08, vmax=.16))
    ax.set_xticks(range(7),axes); ax.set_yticks(range(71),summary.index,fontsize=5)
    ax.set_xlabel("Conditioning axis"); ax.set_ylabel("Disease (ordered by joint attenuation burden)")
    ax.set_title("Disease-level summary of genetic-sharing attenuation",loc="left",weight="bold")
    fig.colorbar(im,ax=ax,label="Median |baseline rg| − |conditioned rg|",fraction=.03,pad=.02)
    save(fig,2)

def fig3():
    d=pd.read_csv(ROOT/"results/phase3_joint/FIGURE3_SCATTER_DATA.tsv",sep="\t")
    colors={"PERSISTENT":NAVY,"ATTENUATED":TEAL,"NEAR_NULL_EQUIVALENT":ORANGE,"MASKED_ENHANCED":"#9B6AB0","UNSTABLE":"#B8C0C8"}
    fig,ax=plt.subplots(figsize=(5.0,4.6))
    for cls,g in d.groupby("classification"):
        ax.scatter(g.baseline_abs_rg,g.joint_abs_rg,s=7,alpha=.65,label=cls.replace("_"," ").title(),c=colors.get(cls,GREY),rasterized=True)
    lim=max(d.baseline_abs_rg.max(),d.joint_abs_rg.max())*1.03; ax.plot([0,lim],[0,lim],ls="--",c=GREY,lw=.8)
    ax.set(xlim=(0,lim),ylim=(0,lim),xlabel="Baseline |rg|",ylabel="Joint residual |rg|",title="Baseline versus jointly conditioned genetic correlation")
    ax.legend(frameon=False,fontsize=6,loc="upper left",markerscale=1.5); save(fig,3)

def fig4():
    d=pd.read_csv(ROOT/"results/phase3_joint/DISEASE_LEVEL_JOINT_SUMMARY.tsv",sep="\t").sort_values("median_joint_attenuation",ascending=False).head(20).sort_values("median_joint_attenuation")
    fig,ax=plt.subplots(figsize=(5.6,5.0)); ax.barh(d.disease,d.median_joint_attenuation,color=BLUE)
    ax.set_xlabel("Median joint attenuation"); ax.set_ylabel(""); ax.set_title("Diseases with highest median joint attenuation",loc="left",weight="bold")
    ax.tick_params(axis="y",labelsize=7); save(fig,4)

def fig5():
    w=pd.read_csv(ROOT/"results/phase4/MODEL1_MVMR_IVW.tsv",sep="\t")
    f=pd.read_csv(ROOT/"results/phase4/MODEL1_MVMR_IVW_FDR.tsv",sep="\t")
    q=f.pivot(index="outcome",columns="exposure",values="q")
    w=w.sort_values("SMK_beta",ascending=True).reset_index(drop=True); y=np.arange(len(w))
    fig,ax=plt.subplots(figsize=(6.2,6.2))
    for i,r in w.iterrows():
        for exp, col in [("BMI",ORANGE),("SMK",TEAL)]:
            b=r[f"{exp}_beta"]; se=r[f"{exp}_SE"]; sig=q.loc[r.outcome,exp] < .05
            ax.errorbar(b,i+(-.13 if exp=="BMI" else .13),xerr=1.96*se,fmt="o",ms=4,color=col,ecolor=col,elinewidth=.9,capsize=2,alpha=1 if sig else .45)
    ax.axvline(0,color=GREY,lw=.8); ax.set_yticks(y,w.outcome); ax.set_xlabel("Direct genetically predicted effect (per SD; log-odds scale)")
    ax.set_title("Primary MVMR: BMI–SMK contrast",loc="left",weight="bold")
    ax.text(.02,.02,"BMI: 0/15 FDR-significant   SMK: 15/15 positive FDR-significant",transform=ax.transAxes,fontsize=8,color=NAVY,weight="bold")
    ax.scatter([],[],c=ORANGE,label="BMI"); ax.scatter([],[],c=TEAL,label="SMK"); ax.legend(frameon=False,loc="upper right")
    save(fig,5)

def fig6():
    d=pd.read_csv(ROOT/"results/final/FIGURE6_VARIANT_PROFILES.tsv",sep="\t").sort_values(["COPD","CHD","T2D","OA"]).reset_index(drop=True)
    cols=["COPD","CHD","T2D","OA"]
    fig,(ax,cax)=plt.subplots(1,2,figsize=(6.4,7.5),gridspec_kw={"width_ratios":[5,1],"wspace":.08})
    arr=d[cols].to_numpy(); im=ax.imshow(arr,aspect="auto",cmap="viridis")
    ax.set_xticks(range(4),cols); ax.set_yticks([]); ax.set_xlabel("Outcome profile"); ax.set_ylabel("54 independent SMK instruments (sorted)")
    ax.set_title("Variant-level smoking profiles",loc="left",weight="bold")
    cax.barh(np.arange(len(d)),d.SMK_beta,color=GREY,height=.8); cax.axvline(0,color=NAVY,lw=.7); cax.set_title("SMK\nbeta",fontsize=8); cax.set_yticks([])
    fig.colorbar(im,ax=ax,label="Variant-specific outcome estimate",fraction=.04,pad=.02)
    fig.text(.5,.02,"MR-AHC: no robust multi-cluster solution; descriptive profiles are exploratory",ha="center",fontsize=8,color=NAVY)
    save(fig,6)

if __name__ == "__main__":
    for fn in [fig1,fig2,fig3,fig4,fig5,fig6]: fn()
    print(f"Rendered six figures to {OUT}")
