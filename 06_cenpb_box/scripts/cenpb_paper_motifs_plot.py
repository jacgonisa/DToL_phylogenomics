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
tiers=["canonical","broad","degenerated"]; tlab={"canonical":"canonical","broad":"broad","degenerated":"degenerate"}
tcol={"canonical":"#D55E00","broad":"#0072B2","degenerated":"#E69F00"}   # Okabe-Ito
def bd(ctrl): return {t:(bench.loc[bench.control==ctrl,f"{t}_perMbp"].iloc[0], bench.loc[bench.control==ctrl,f"{t}_enrich_dinuc"].iloc[0]) for t in tiers}
def cd(cl):   return {t:(per.loc[per.clade==cl,f"{t}_perMbp"].iloc[0],        per.loc[per.clade==cl,f"{t}_enrich_dinuc"].iloc[0])       for t in tiers}
# (density hits/Mbp, enrichment over dinucleotide null) per group x tier
groups=[("α-sat (HG002, +)","*",230, bd("alpha-satellite (positive)")),
        ("HSat (HG002, −)", "X",70,  bd("HSat1/2/3 (negative)")),
        ("Vertebrates",     "^",55,  cd("Vertebrates")),
        ("Invertebrates",   "o",46,  cd("Invertebrate")),
        ("Plants",          "s",50,  cd("Viridiplantae"))]

# ---- Nature-style aesthetics ----
plt.rcParams.update({
  "font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
  "font.size":8,"axes.linewidth":0.7,"xtick.major.width":0.7,"ytick.major.width":0.7,
  "xtick.major.size":3,"ytick.major.size":3,"xtick.direction":"out","ytick.direction":"out",
  "axes.labelsize":8,"legend.fontsize":6.8,"xtick.labelsize":7.5,"ytick.labelsize":7.5,"savefig.dpi":600})

# scatter: enrichment (obs/exp) vs motif density (hits/Mbp) — a real box needs BOTH high
fig,ax=plt.subplots(figsize=(5.8,4.4))
for glab,mk,sz,dd in groups:
    for t in tiers:
        xden,yenr=dd[t]
        if xden>0 and yenr>0:
            ax.scatter(xden,yenr,marker=mk,s=sz,facecolor=tcol[t],edgecolor="k",lw=0.5,zorder=3)
ax.set_xscale("log"); ax.set_yscale("log")
ax.axhline(1,ls=(0,(4,3)),color="0.4",lw=0.9,zorder=1); ax.text(ax.get_xlim()[0]*1.3 if False else 2e-3,1.15,"enrichment = 1 (null)",fontsize=6.3,color="0.4")
# annotations
a=bd("alpha-satellite (positive)"); h=bd("HSat1/2/3 (negative)")
ax.annotate("α-satellite (functional box)\nhigh density + high enrichment",(a["broad"][0],a["broad"][1]),
            xytext=(-6,-24),textcoords="offset points",ha="right",fontsize=6.4,fontweight="bold")
ax.annotate("HSat canonical:\n5 boundary hits (0.03/Mbp)\n→ high enrichment but ~0 density",(h["canonical"][0],h["canonical"][1]),
            xytext=(8,4),textcoords="offset points",fontsize=6,color="#C62828")
for cl,lab in [("Vertebrates","Vertebrates"),("Invertebrate","Invertebrates"),("Viridiplantae","Plants")]:
    xy=cd(cl)["broad"]; ax.annotate(lab,xy,xytext=(4,-8),textcoords="offset points",fontsize=6.3)
ax.set_xlabel("CENP-B motif density (hits per Mbp)")
ax.set_ylabel("enrichment over dinucleotide null (obs / exp)")
ax.set_title(f"Exact IUPAC CENP-B motif ({n_sat} species with satellites)",fontsize=8.8,fontweight="bold",pad=4)
ax.text(0.98,0.03,"canonical = 0 in all DToL clades (not plotted)",transform=ax.transAxes,ha="right",fontsize=6.2,color="#D55E00")
# legends: tier colour + group marker
from matplotlib.lines import Line2D
tl=[Line2D([0],[0],marker="o",ls="",mfc=tcol[t],mec="k",mew=0.4,label=tlab[t]) for t in tiers]
gl=[Line2D([0],[0],marker=mk,ls="",mfc="0.7",mec="k",mew=0.4,label=glab) for glab,mk,_,_ in groups]
l1=ax.legend(handles=tl,title="motif tier",frameon=False,loc="upper left",fontsize=6.6,title_fontsize=6.8)
ax.add_artist(l1); ax.legend(handles=gl,title="group",frameon=False,loc="lower right",fontsize=6.4,title_fontsize=6.8)
ax.spines[["top","right"]].set_visible(False)
# motif-definition box: what Barra & Fachinetti define as each IUPAC tier
from matplotlib.patches import FancyBboxPatch
mdefs=[("canonical","YTTCGTTGGAARCGGGA","canonical"),
       ("broad",    "NTTCGNNNNANNCGGGN","broad"),
       ("degenerate","YTTCGNNNNANRCGGGN","degenerated")]
bx,byt=0.28,0.80
ax.add_patch(FancyBboxPatch((bx-0.012,byt-0.205),0.52,0.215,transform=ax.transAxes,
             boxstyle="round,pad=0.006",fc="white",ec="0.75",lw=0.7,zorder=8))
ax.text(bx,byt,"IUPAC tiers (Barra & Fachinetti):",transform=ax.transAxes,va="top",ha="left",
        fontsize=5.9,fontweight="bold",color="0.25",zorder=9)
for i,(nm,seq,key) in enumerate(mdefs):
    ax.text(bx,byt-0.058-i*0.05,f"{nm:<11s}{seq}",transform=ax.transAxes,va="top",ha="left",
            fontsize=5.9,family="DejaVu Sans Mono",color=tcol[key],zorder=9)
plt.tight_layout()
for ext in ("png","pdf"):
    fig.savefig(SAT/f"figures/cenpb_paper_motifs.{ext}",dpi=600 if ext=="png" else None,bbox_inches="tight",facecolor="white")
print(f"Saved figures/cenpb_paper_motifs.png/pdf | species with satellites: {n_sat}")
