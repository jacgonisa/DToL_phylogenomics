#!/usr/bin/env python3
"""Plot the de-novo best-window analysis (cenpb_denovo_bestwindow.tsv): best-window
identity to canonical in the real satellite (y) vs a dinucleotide shuffle (x). Points
above the diagonal contain a box-like window above chance. Only alpha-sat is clearly
box-like; birds and the rest sit on the null."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
SAT=Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity")
d=pd.read_csv(SAT/"figures/cenpb_denovo_bestwindow.tsv",sep="\t")
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
  "font.size":8,"axes.linewidth":0.7,"xtick.major.width":0.7,"ytick.major.width":0.7,
  "axes.labelsize":8,"legend.fontsize":6.6,"xtick.labelsize":7.5,"ytick.labelsize":7.5,"savefig.dpi":600})
ccol={"Vertebrates":"#0072B2","Invertebrate":"#E69F00","Viridiplantae":"#009E73","Fungi":"#CC79A7","Protist":"#C62828"}
dt=d[~d.clade.str.contains("Human")]; asat=d[d.name.str.contains("alpha")].iloc[0]; hsat=d[d.name.str.contains("HSat")].iloc[0]
fig,ax=plt.subplots(figsize=(5.2,4.6))
ax.plot([50,85],[50,85],ls=(0,(4,3)),color="0.5",lw=0.9,zorder=1)
for cl,g in dt.groupby("clade"):
    ax.scatter(g.best_null,g.best_obs,s=16,color=ccol.get(cl,"grey"),edgecolor="0.3",lw=0.2,alpha=0.8,zorder=3,label=cl)
ax.scatter(asat.best_null,asat.best_obs,marker="*",s=240,color="#C1272D",edgecolor="k",lw=0.6,zorder=6)
ax.annotate("α-satellite (HG002, +)\nreal box (+17)",(asat.best_null,asat.best_obs),xytext=(-6,-24),textcoords="offset points",ha="right",fontsize=6.6,fontweight="bold")
ax.scatter(hsat.best_null,hsat.best_obs,marker="X",s=70,color="0.4",edgecolor="k",lw=0.5,zorder=6)
ax.annotate("HSat (−)",(hsat.best_null,hsat.best_obs),xytext=(7,-2),textcoords="offset points",fontsize=6.4,color="0.3")
for _,r in dt[dt.name.str.contains("Accipiter|Lagopus|Porphyrio|Corvus",na=False)].iterrows():
    ax.annotate(r["name"].split()[0],(r.best_null,r.best_obs),xytext=(4,3),textcoords="offset points",fontsize=6,fontweight="bold")
ax.set_xlabel("best-window identity in shuffled satellite (%) — null")
ax.set_ylabel("best-window identity in real satellite (%)")
ax.set_title("De-novo box-like signal (no CENP-B seeding)\npoints above the diagonal contain a box-like window above chance",fontsize=8.6,fontweight="bold")
ax.text(0.03,0.97,"only α-satellite is clearly\nbox-like; birds sit on the null",transform=ax.transAxes,va="top",fontsize=6.6,color="0.3")
ax.legend(frameon=False,loc="lower right",title="DToL clade",title_fontsize=6.8)
ax.spines[["top","right"]].set_visible(False); plt.tight_layout()
for ext in ("png","pdf"): fig.savefig(SAT/f"figures/cenpb_denovo_bestwindow.{ext}",dpi=600 if ext=="png" else None,bbox_inches="tight",facecolor="white")
print("Saved figures/cenpb_denovo_bestwindow.png/pdf")
