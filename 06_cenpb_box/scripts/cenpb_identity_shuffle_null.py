#!/usr/bin/env python3
"""Composition null for 'identity to the canonical CENP-B motif'. Short (17 bp)
motifs match partly by base composition alone, so a raw %identity is hard to read.
For each species' box consensus we compare its identity to canonical against the
MEAN identity of many shuffles of that same consensus (same composition, sequence
order destroyed). excess = observed − shuffled tells us how much of the match is
real arrangement vs composition chance."""
import random, numpy as np, pandas as pd
from pathlib import Path
random.seed(0)
SAT=Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity")
CANON=[set("CT"),{"T"},{"T"},{"C"},{"G"},{"T"},{"T"},{"G"},{"G"},{"A"},{"A"},set("AG"),{"C"},{"G"},{"G"},{"G"},{"A"}]
def idto(seq): return 100*sum(seq[i] in CANON[i] for i in range(17))/17
K=3000

f=pd.read_csv(SAT/"figures/cenpb_flank_uncapped_per_species.tsv",sep="\t")
rows=[]
for _,r in f.iterrows():
    cons=r.box_consensus
    if not isinstance(cons,str) or len(cons)!=17 or (set(cons)-set("ACGT")): continue
    obs=idto(cons); ch=list(cons); sh=np.empty(K)
    for k in range(K): random.shuffle(ch); sh[k]=idto(ch)
    rows.append(dict(species=r.species, id_obs=round(obs,1), id_shuffle=round(sh.mean(),1),
                     id_shuffle_p95=round(np.percentile(sh,95),1),
                     excess_shuffle=round(obs-sh.mean(),1)))
df=pd.DataFrame(rows); df.to_csv(SAT/"figures/cenpb_identity_shuffle_null.tsv",sep="\t",index=False)

# reference: shuffle the canonical motif itself (composition chance for the box)
cn=list("CTTCGTTGGAAACGGGA"); base=np.mean([idto(random.sample(cn,17)) for _ in range(20000)])
print(f"canonical-composition chance identity (shuffle the box): {base:.1f}%")
print(f"per-species: median id_obs={df.id_obs.median():.1f}%  shuffle={df.id_shuffle.median():.1f}%  excess={df.excess_shuffle.median():.1f}%")
fk=df.merge(f[["species","name","vgroup","clade","delta"]],on="species",how="left").reset_index(drop=True)
print("\nbirds + key:")
print(fk[fk.name.str.contains("Accipiter|Lagopus|Porphyrio|Corvus|Cervus",na=False)]
      [["name","id_obs","id_shuffle","excess_shuffle","delta"]].to_string(index=False))

# ---- plot: observed identity vs shuffle-the-motif null ----
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
  "font.size":8,"axes.linewidth":0.7,"xtick.major.width":0.7,"ytick.major.width":0.7,
  "axes.labelsize":8,"legend.fontsize":6.6,"xtick.labelsize":7.5,"ytick.labelsize":7.5,"savefig.dpi":600})
ccol={"Vertebrates":"#0072B2","Invertebrate":"#E69F00","Viridiplantae":"#009E73","Fungi":"#CC79A7","Protist":"#D55E00"}
fig,ax=plt.subplots(figsize=(4.7,4.4))
# slopegraph: link each observed motif to the mean of its own shuffled null
for _,r in fk.iterrows():
    ax.plot([0,1],[r.id_shuffle,r.id_obs],color=ccol.get(r.clade,"grey"),alpha=0.22,lw=0.5,zorder=2)
ax.scatter([0]*len(fk),fk.id_shuffle,s=11,color="0.55",edgecolor="none",alpha=0.7,zorder=3)
for cl,g in fk.groupby("clade"):
    m=(fk.clade==cl).values
    ax.scatter([1]*int(m.sum()),fk.id_obs[m],s=16,color=ccol.get(cl,"grey"),edgecolor="0.3",lw=0.2,alpha=0.9,zorder=4,label=cl)
for xx,col in [(0,"id_shuffle"),(1,"id_obs")]:                 # median tick (no text)
    ax.plot([xx-0.09,xx+0.09],[fk[col].median()]*2,color="k",lw=2,zorder=6)
# formula box
ax.text(0.02,0.99,"identity = (positions matching\ncanonical IUPAC) / 17 × 100%\nnull = same, on shuffles of the\nmotif (composition preserved)",
        transform=ax.transAxes,va="top",ha="left",fontsize=6,family="DejaVu Sans",
        bbox=dict(boxstyle="round,pad=0.4",fc="#f7f7f7",ec="0.7",lw=0.7))
ax.set_xticks([0,1]); ax.set_xticklabels(["shuffled motif\n(composition chance)","observed\nmotif"])
ax.set_xlim(-0.35,1.35); ax.set_ylim(0,100)
ax.set_ylabel("% identity to canonical CENP-B motif")
ax.set_title("Identity vs shuffle-the-motif null",fontweight="bold",fontsize=9)
ax.legend(frameon=False,loc="lower right",fontsize=6.3)
ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
for ext in ("png","pdf"): fig.savefig(SAT/f"figures/cenpb_identity_shuffle_null.{ext}",dpi=600 if ext=="png" else None,bbox_inches="tight",facecolor="white")
print("Saved figures/cenpb_identity_shuffle_null.png/pdf")
