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
n_sat=int(per.n_species.sum())                          # species that actually have satellites
# enrichment over the dinucleotide-preserving (Markov-1 / Altschul-Erikson) null
benchd={blab[r.control]:{t:r[f"{t}_enrich_dinuc"] for t in ["canonical","broad","degenerated"]} for _,r in bench.iterrows()}
tiers=["broad","degenerated"]; tlab={"broad":"broad","degenerated":"degenerate"}
tcol={"broad":"#0072B2","degenerated":"#E69F00"}        # Okabe-Ito
enr=lambda cl,t: per.loc[per.clade==cl,f"{t}_enrich_dinuc"].iloc[0]
G=[("HSat 1/2/3\n(HG002, −)", {t:benchd["HSat (HG002)"][t] for t in tiers}),
   ("Vertebrates",   {t:enr("Vertebrates",t)   for t in tiers}),
   ("Invertebrates", {t:enr("Invertebrate",t)  for t in tiers}),
   ("Plants",        {t:enr("Viridiplantae",t) for t in tiers})]
labels=[g[0] for g in G]; x=np.arange(len(G)); w=0.36

# ---- Nature-style aesthetics ----
plt.rcParams.update({
  "font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
  "font.size":8,"axes.linewidth":0.7,"xtick.major.width":0.7,"ytick.major.width":0.7,
  "xtick.major.size":3,"ytick.major.size":3,"xtick.direction":"out","ytick.direction":"out",
  "axes.labelsize":8,"legend.fontsize":7,"xtick.labelsize":7.5,"ytick.labelsize":7.5,"savefig.dpi":600})

fig,ax=plt.subplots(figsize=(5.0,3.3))
for j,t in enumerate(tiers):
    vals=[g[1][t] for g in G]
    ax.bar(x+(j-0.5)*w,vals,w,color=tcol[t],label=tlab[t],edgecolor="none",zorder=3)
    for xi,v in zip(x+(j-0.5)*w,vals):
        ax.annotate(f"{v:.2f}",(xi,v),ha="center",va="bottom",fontsize=6,color=tcol[t])
ax.axhline(1,ls=(0,(4,3)),color="0.35",lw=0.9,zorder=2)
ax.text(len(G)-0.55,1.03,"null (random)",ha="right",va="bottom",fontsize=6.5,color="0.35")
ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_ylabel("enrichment over dinucleotide null  (obs / exp)")
ax.set_ylim(0,3.35)
ax.set_title("Exact IUPAC CENP-B box enrichment",fontsize=9,fontweight="bold",pad=20)
ax.text(0.5,1.045,f"canonical (functional) box = 0 in all DToL clades ({n_sat} species w/ satellites);\n"
        "α-satellite (HG002,+) positive control off scale — broad 4,283× · degenerate 14,457×",
        transform=ax.transAxes,ha="center",va="bottom",fontsize=6.3,color="0.3")
ax.legend(frameon=False,handlelength=1.1,title="motif tier",title_fontsize=7,loc="upper center")
ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
for ext in ("png","pdf"):
    fig.savefig(SAT/f"figures/cenpb_paper_motifs.{ext}",dpi=600 if ext=="png" else None,bbox_inches="tight",facecolor="white")
print(f"Saved figures/cenpb_paper_motifs.png/pdf | species with satellites: {n_sat}")
