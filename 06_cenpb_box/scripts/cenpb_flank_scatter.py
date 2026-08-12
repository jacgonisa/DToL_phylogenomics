#!/usr/bin/env python3
"""Method 2 (songbird ±5-flank) figure. Points coloured by IDENTITY to the
canonical CENP-B motif (= (17 − substitutions)/17), so a 'candidate box' =
high identity AND above the box=flank diagonal (box-specific). Human HG002
α-satellite (positive) and HSat (negative) shown as reference markers.
Terminology: 'box' = functional (protein-binding); sequence hits are motifs /
candidate boxes."""
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
SAT=Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity")
df=pd.read_csv(SAT/"figures/cenpb_flank_uncapped_per_species.tsv",sep="\t"); d=df[df.n_windows>=5].copy()
# colour = identity of the species' consensus motif to the canonical CENP-B motif
#   (17 − substitutions)/17. Descriptive: how box-like the dominant motif is.
d["identity"]=100*(17-d.subs_vs_canonical)/17
b=pd.read_csv(SAT/"figures/cenpb_human_benchmark.tsv",sep="\t")
hm={r.control:(r.mean_flank_bits,r.mean_box_bits) for _,r in b.iterrows()}
ASAT=hm["alpha-satellite (positive)"]; HSAT=hm["HSat1/2/3 (negative)"]

plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
  "font.size":8,"axes.linewidth":0.7,"xtick.major.width":0.7,"ytick.major.width":0.7,
  "xtick.major.size":3,"ytick.major.size":3,"axes.labelsize":8,"legend.fontsize":6.6,
  "xtick.labelsize":7.5,"ytick.labelsize":7.5,"savefig.dpi":600})
CM="YlOrRd"; VMIN,VMAX=70,100                           # identity to canonical (%): warm = box-like
shape={"Invertebrate":"o","Vertebrates":"^","Viridiplantae":"s","Fungi":"D","Protist":"v"}

def bench(ax):
    ax.plot([0,2],[0,2],ls=(0,(4,3)),color="0.5",lw=0.9,zorder=1)
    ax.scatter(*ASAT,marker="*",s=190,color="#C1272D",edgecolor="k",lw=0.6,zorder=7)
    ax.annotate("α-sat (HG002, +)\nfunctional box",ASAT,xytext=(5,-14),textcoords="offset points",fontsize=6.2,fontweight="bold")
    ax.scatter(*HSAT,marker="X",s=55,color="0.45",edgecolor="k",lw=0.5,zorder=7)
    ax.annotate("HSat (HG002, −)",HSAT,xytext=(5,3),textcoords="offset points",fontsize=6.2,color="0.3")

fig,(a1,a2)=plt.subplots(1,2,figsize=(7.3,3.5))
# left: all clades (marker = clade), colour = identity to canonical motif
for cl,g in d.groupby("clade"):
    a1.scatter(g.mean_flank_bits,g.mean_box_bits,c=g.identity,cmap=CM,vmin=VMIN,vmax=VMAX,
               marker=shape.get(cl,"o"),s=26,edgecolor="0.3",lw=0.3,zorder=3)
bench(a1)
a1.set_xlim(0,2); a1.set_ylim(0,2); a1.set_xlabel("flank information (bits) — control"); a1.set_ylabel("motif information (bits)")
a1.set_title("All clades",fontsize=8.5,fontweight="bold",pad=3)
a1.legend(handles=[Line2D([0],[0],marker=m,ls="",mfc="0.7",mec="0.3",label=c) for c,m in shape.items() if c in d.clade.values],
          frameon=False,loc="lower right",handletextpad=0.2,title="clade",title_fontsize=6.6)
a1.spines[["top","right"]].set_visible(False)
# right: vertebrates, colour = identity, size ∝ boxes/Mbp
v=d[d.clade=="Vertebrates"].copy()
sc=a2.scatter(v.mean_flank_bits,v.mean_box_bits,c=v.identity,cmap=CM,vmin=VMIN,vmax=VMAX,
              s=np.clip(v.win_per_Mbp,14,320),edgecolor="k",lw=0.4,zorder=3)
bench(a2)
for _,r in v[(v.delta>=0.4)|v.name.str.contains("Accipiter|Lagopus|Porphyrio|Cervus|Diceros")].iterrows():
    a2.annotate(r["name"].split()[0],(r.mean_flank_bits,r.mean_box_bits),fontsize=6.2,xytext=(3,2),textcoords="offset points")
a2.set_xlim(0,2); a2.set_ylim(0,2); a2.set_xlabel("flank information (bits) — control"); a2.set_ylabel("motif information (bits)")
a2.set_title("Vertebrates (size ∝ candidate boxes/Mbp)",fontsize=8.5,fontweight="bold",pad=3)
a2.text(0.03,0.97,"warm = box-like\n(≫ shuffle chance ≈ 30%)",transform=a2.transAxes,fontsize=6.2,color="0.35",va="top")
a2.spines[["top","right"]].set_visible(False)
cb=fig.colorbar(sc,ax=a2,fraction=0.046,pad=0.03); cb.set_label("identity to canonical CENP-B motif (%)",fontsize=6.8)
cb.ax.tick_params(labelsize=6.5)
for ax,l in [(a1,"a"),(a2,"b")]: ax.text(-0.17,1.08,l,transform=ax.transAxes,fontsize=11,fontweight="bold",va="top")
plt.tight_layout(w_pad=2.2)
for ext in ("png","pdf"):
    fig.savefig(SAT/f"figures/cenpb_flank_uncapped_scatter.{ext}",dpi=600 if ext=="png" else None,bbox_inches="tight",facecolor="white")
print("Saved figures/cenpb_flank_uncapped_scatter.png/pdf")
