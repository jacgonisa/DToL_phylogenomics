#!/usr/bin/env python3
"""Figure for Method 1 (Fachinetti / exact-IUPAC motif parsing). Reads the
325-restricted per-species + per-clade tables from cenpb_paper_motifs.py and
adds the human HG002 alpha-satellite as the positive benchmark (same three exact
IUPAC tiers). Panel A: per-clade enrichment over random for canonical/broad/
degenerate. Panel B: vertebrate broad/degenerate hits per Mbp."""
import re, numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
SAT=Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity")
HUMAN=SAT/"cenpb_psi/human_alpha.fasta"
MOT={"canonical":re.compile("[CT]TTCGTTGGAA[AG]CGGGA"),"broad":re.compile(".TTCG....A..CGGG."),
     "degenerated":re.compile("[CT]TTCG....A.[AG]CGGG.")}
ALLOWED={"canonical":[{"C","T"},{"T"},{"T"},{"C"},{"G"},{"T"},{"T"},{"G"},{"G"},{"A"},{"A"},{"A","G"},{"C"},{"G"},{"G"},{"G"},{"A"}],
 "broad":[set("ACGT"),{"T"},{"T"},{"C"},{"G"},set("ACGT"),set("ACGT"),set("ACGT"),set("ACGT"),{"A"},set("ACGT"),set("ACGT"),{"C"},{"G"},{"G"},{"G"},set("ACGT")],
 "degenerated":[{"C","T"},{"T"},{"T"},{"C"},{"G"},set("ACGT"),set("ACGT"),set("ACGT"),set("ACGT"),{"A"},set("ACGT"),{"A","G"},{"C"},{"G"},{"G"},{"G"},set("ACGT")]}
comp=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]

# HG002 benchmarks (alpha-sat positive, HSat negative) from cenpb_human_benchmark.py
# NULL = dinucleotide-preserving (first-order Markov / Altschul-Erikson) enrichment
bench=pd.read_csv(SAT/"figures/cenpb_human_benchmark.tsv",sep="\t")
blab={"alpha-satellite (positive)":"α-sat (HG002)","HSat1/2/3 (negative)":"HSat (HG002)"}
benchg={blab[r.control]:{t:r[f"{t}_enrich_dinuc"] for t in ["canonical","broad","degenerated"]} for _,r in bench.iterrows()}
print("benchmarks (dinuc null):",benchg)

per=pd.read_csv(SAT/"figures/cenpb_paper_motifs_per_clade.tsv",sep="\t")
tiers=["canonical","broad","degenerated"]
tlab={"canonical":"canonical","broad":"broad","degenerated":"degenerate"}
tcol={"canonical":"#D55E00","broad":"#0072B2","degenerated":"#E69F00"}   # Okabe-Ito (colour-blind safe)
benchpm={blab[r.control]:{t:r[f"{t}_perMbp"] for t in tiers} for _,r in bench.iterrows()}

# ---- Nature-style aesthetics ----
plt.rcParams.update({
  "font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
  "font.size":8,"axes.linewidth":0.7,"xtick.major.width":0.7,"ytick.major.width":0.7,
  "xtick.major.size":3,"ytick.major.size":3,"xtick.direction":"out","ytick.direction":"out",
  "axes.labelsize":8,"legend.fontsize":7,"xtick.labelsize":8,"ytick.labelsize":7.5,"savefig.dpi":600})

fig,(a1,a2)=plt.subplots(1,2,figsize=(7.09,2.9))       # ~180 mm double-column
w=0.26

# Panel a — human controls: functional-box density (per Mbp), linear
ga=["α-sat (HG002)","HSat (HG002)"]; xa=np.arange(len(ga))
for j,t in enumerate(tiers):
    a1.bar(xa+(j-1)*w,[benchpm[g][t] for g in ga],w,color=tcol[t],label=tlab[t],edgecolor="none",zorder=3)
for j,t in enumerate(tiers):                            # annotate the near-zero HSat bars
    v=benchpm["HSat (HG002)"][t]
    a1.annotate(f"{v:.2g}",(1+(j-1)*w,v),ha="center",va="bottom",fontsize=6,color=tcol[t])
a1.set_xticks(xa); a1.set_xticklabels(["α-satellite\n(positive)","HSat 1/2/3\n(negative)"])
a1.set_ylabel("CENP-B box hits per Mbp")
a1.set_title("Human HG002 controls",fontsize=8.5,fontweight="bold",pad=4)
a1.legend(frameon=False,handlelength=1.1,title="motif tier",title_fontsize=7,loc="upper right")
a1.spines[["top","right"]].set_visible(False)

# Panel b — DToL clades: enrichment over the dinucleotide-preserving null, linear
gb=list(per.clade); lab={"Vertebrates":"Vertebrates","Invertebrate":"Invertebrates","Viridiplantae":"Plants"}
xb=np.arange(len(gb))
for j,t in enumerate(tiers):
    vals=[per.loc[per.clade==g,f"{t}_enrich_dinuc"].iloc[0] for g in gb]
    a2.bar(xb+(j-1)*w,vals,w,color=tcol[t],label=tlab[t],edgecolor="none",zorder=3)
a2.axhline(1,ls=(0,(4,3)),color="0.35",lw=0.9,zorder=2)
a2.text(-0.45,1.05,"null (random)",ha="left",va="bottom",fontsize=6.5,color="0.35")
a2.annotate("canonical box = 0\nin all 325 species",xy=(0.03,0.9),xycoords="axes fraction",
            fontsize=6.8,color=tcol["canonical"],va="top")
a2.set_xticks(xb); a2.set_xticklabels([lab[g] for g in gb])
a2.set_ylabel("enrichment over dinucleotide null")
a2.set_title("DToL satellites (325 species)",fontsize=8.5,fontweight="bold",pad=4)
a2.set_ylim(0,3.4)
a2.spines[["top","right"]].set_visible(False)

for ax,l in [(a1,"a"),(a2,"b")]:
    ax.text(-0.16,1.10,l,transform=ax.transAxes,fontsize=11,fontweight="bold",va="top")
plt.tight_layout(w_pad=2.0)
for ext in ("png","pdf"):
    fig.savefig(SAT/f"figures/cenpb_paper_motifs.{ext}",dpi=600 if ext=="png" else None,bbox_inches="tight",facecolor="white")
print("Saved figures/cenpb_paper_motifs.png/pdf")
