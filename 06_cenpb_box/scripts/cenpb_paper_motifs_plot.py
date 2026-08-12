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
n_sat=int(per.n_species.sum())                          # species that actually have satellites

# same metric (hits per Mbp) for controls AND DToL clades, so panels are comparable
G=[("α-satellite\n(HG002, +)",  {t:benchpm["α-sat (HG002)"][t] for t in tiers}),
   ("HSat 1/2/3\n(HG002, −)",   {t:benchpm["HSat (HG002)"][t]  for t in tiers}),
   ("Vertebrates",  {t:per.loc[per.clade=="Vertebrates",  f"{t}_perMbp"].iloc[0] for t in tiers}),
   ("Invertebrates",{t:per.loc[per.clade=="Invertebrate", f"{t}_perMbp"].iloc[0] for t in tiers}),
   ("Plants",       {t:per.loc[per.clade=="Viridiplantae",f"{t}_perMbp"].iloc[0] for t in tiers})]
labels=[g[0] for g in G]; x=np.arange(len(G)); w=0.26

# ---- Nature-style aesthetics ----
plt.rcParams.update({
  "font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
  "font.size":8,"axes.linewidth":0.7,"xtick.major.width":0.7,"ytick.major.width":0.7,
  "xtick.major.size":3,"ytick.major.size":3,"xtick.direction":"out","ytick.direction":"out",
  "axes.labelsize":8,"legend.fontsize":7,"xtick.labelsize":7.5,"ytick.labelsize":7.5,"savefig.dpi":600})

# broken y-axis: alpha-sat is ~500x the rest, shown on one shared per-Mbp scale
fig,(top,bot)=plt.subplots(2,1,sharex=True,figsize=(7.09,3.7),
                           gridspec_kw=dict(height_ratios=[1,1.7],hspace=0.06))
for j,t in enumerate(tiers):
    vals=[g[1][t] for g in G]
    top.bar(x+(j-1)*w,vals,w,color=tcol[t],edgecolor="none",zorder=3)
    bot.bar(x+(j-1)*w,vals,w,color=tcol[t],label=tlab[t],edgecolor="none",zorder=3)
top.set_ylim(1400,2250); bot.set_ylim(0,4.4)            # break between the two ranges
# annotate the small (near-zero) bars in the lower panel
for gi,(gname,gv) in enumerate(G):
    if gi==0: continue
    for j,t in enumerate(tiers):
        v=gv[t]
        if v<4.4: bot.annotate(f"{v:.2g}",(gi+(j-1)*w,v),ha="center",va="bottom",fontsize=5.6,color=tcol[t])
# broken-axis marks
top.spines["bottom"].set_visible(False); bot.spines["top"].set_visible(False)
top.tick_params(bottom=False)
dm=dict(marker=[(-1,-0.5),(1,0.5)],markersize=7,linestyle="none",color="k",mec="k",mew=0.8,clip_on=False)
top.plot([0,1],[0,0],transform=top.transAxes,**dm); bot.plot([0,1],[1,1],transform=bot.transAxes,**dm)
for ax in (top,bot): ax.spines[["right"]].set_visible(False)
bot.axvline(1.5,color="0.8",lw=0.8,ls=":",zorder=1)     # divider: human controls | DToL clades
top.axvline(1.5,color="0.8",lw=0.8,ls=":",zorder=1)
bot.set_xticks(x); bot.set_xticklabels(labels)
bot.set_ylabel("CENP-B box hits per Mbp"); bot.yaxis.set_label_coords(-0.075,0.72)
bot.legend(frameon=False,handlelength=1.1,title="motif tier",title_fontsize=7,loc="upper center",ncol=3)
top.set_title(f"Exact IUPAC CENP-B box density — human controls vs DToL satellites ({n_sat} species with satellites)",
              fontsize=8.5,fontweight="bold",pad=4)
top.annotate("canonical box present\nonly in α-satellite",xy=(0,1562),xytext=(0.35,1780),
             fontsize=6.6,color=tcol["canonical"],ha="left")
plt.tight_layout()
for ext in ("png","pdf"):
    fig.savefig(SAT/f"figures/cenpb_paper_motifs.{ext}",dpi=600 if ext=="png" else None,bbox_inches="tight",facecolor="white")
print(f"Saved figures/cenpb_paper_motifs.png/pdf | species with satellites: {n_sat}")
