#!/usr/bin/env python3
"""Method 2 (songbird ±5-flank) figure: box vs flank information per species, with
the human HG002 benchmarks (α-satellite positive, HSat negative) as reference
points. Left: all clades. Right: vertebrates by subgroup (size ∝ boxes/Mbp)."""
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
SAT=Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity")
df=pd.read_csv(SAT/"figures/cenpb_flank_uncapped_per_species.tsv",sep="\t"); d=df[df.n_windows>=5].copy()
b=pd.read_csv(SAT/"figures/cenpb_human_benchmark.tsv",sep="\t")
hm={r.control:(r.mean_flank_bits,r.mean_box_bits) for _,r in b.iterrows()}
ASAT=hm["alpha-satellite (positive)"]; HSAT=hm["HSat1/2/3 (negative)"]

plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
  "font.size":8,"axes.linewidth":0.7,"xtick.major.width":0.7,"ytick.major.width":0.7,
  "xtick.major.size":3,"ytick.major.size":3,"axes.labelsize":8,"legend.fontsize":6.8,
  "xtick.labelsize":7.5,"ytick.labelsize":7.5,"savefig.dpi":600})
ccol={"Vertebrates":"#0072B2","Invertebrate":"#E69F00","Viridiplantae":"#009E73","Fungi":"#CC79A7","Protist":"#D55E00"}
vcol={"Mammalia":"#8E44AD","Aves":"#0072B2","Fish":"#009E73","Reptilia":"#D55E00","Amphibia":"#56B4E9"}

def bench(ax):
    ax.plot([0,2],[0,2],ls=(0,(4,3)),color="0.5",lw=0.9,zorder=1)
    ax.scatter(*ASAT,marker="*",s=180,color="#C1272D",edgecolor="k",lw=0.5,zorder=6)
    ax.annotate("α-sat (HG002, +)",ASAT,xytext=(4,-9),textcoords="offset points",fontsize=6.6,fontweight="bold")
    ax.scatter(*HSAT,marker="s",s=42,color="0.45",edgecolor="k",lw=0.5,zorder=6)
    ax.annotate("HSat (HG002, −)",HSAT,xytext=(4,3),textcoords="offset points",fontsize=6.6,color="0.3")

fig,(a1,a2)=plt.subplots(1,2,figsize=(7.09,3.4))
# left: all clades
for cl,g in d.groupby("clade"):
    a1.scatter(g.mean_flank_bits,g.mean_box_bits,s=14,color=ccol.get(cl,"grey"),alpha=0.6,edgecolor="none",label=cl,zorder=3)
bench(a1)
a1.set_xlim(0,2); a1.set_ylim(0,2); a1.set_xlabel("flank information (bits) — control"); a1.set_ylabel("box information (bits)")
a1.set_title("All clades",fontsize=8.5,fontweight="bold",pad=3)
a1.legend(frameon=False,loc="lower right",handletextpad=0.2); a1.spines[["top","right"]].set_visible(False)
# right: vertebrates by subgroup
v=d[d.clade=="Vertebrates"].copy()
for g,sub in v.groupby("vgroup"):
    a2.scatter(sub.mean_flank_bits,sub.mean_box_bits,s=np.clip(sub.win_per_Mbp,12,320),color=vcol.get(g,"grey"),
               alpha=0.8,edgecolor="k",lw=0.4,label=g,zorder=3)
bench(a2)
for _,r in v[(v.delta>=0.4)|v.name.str.contains("Accipiter|Lagopus|Porphyrio|Cervus|Diceros")].iterrows():
    a2.annotate(r["name"].split()[0],(r.mean_flank_bits,r.mean_box_bits),fontsize=6.2,xytext=(3,2),textcoords="offset points")
a2.set_xlim(0,2); a2.set_ylim(0,2); a2.set_xlabel("flank information (bits) — control"); a2.set_ylabel("box information (bits)")
a2.set_title("Vertebrates (size ∝ boxes/Mbp)",fontsize=8.5,fontweight="bold",pad=3)
a2.legend(frameon=False,loc="lower right",handletextpad=0.2); a2.spines[["top","right"]].set_visible(False)
a2.text(0.03,0.97,"above diagonal =\nbox-enriched",transform=a2.transAxes,fontsize=6.4,color="0.35",va="top")
for ax,l in [(a1,"a"),(a2,"b")]: ax.text(-0.16,1.08,l,transform=ax.transAxes,fontsize=11,fontweight="bold",va="top")
plt.tight_layout(w_pad=2.0)
for ext in ("png","pdf"):
    fig.savefig(SAT/f"figures/cenpb_flank_uncapped_scatter.{ext}",dpi=600 if ext=="png" else None,bbox_inches="tight",facecolor="white")
print("Saved figures/cenpb_flank_uncapped_scatter.png/pdf | benchmarks alpha",ASAT,"HSat",HSAT)
