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
bench=pd.read_csv(SAT/"figures/cenpb_human_benchmark.tsv",sep="\t")
blab={"alpha-satellite (positive)":"α-sat (HG002)","HSat1/2/3 (negative)":"HSat (HG002)"}
benchg={blab[r.control]:{t:r[f"{t}_enrich"] for t in ["canonical","broad","degenerated"]} for _,r in bench.iterrows()}
print("benchmarks:",benchg)

per=pd.read_csv(SAT/"figures/cenpb_paper_motifs_per_clade.tsv",sep="\t")
df=pd.read_csv(SAT/"figures/cenpb_paper_motifs_per_species.tsv",sep="\t")
tiers=["canonical","broad","degenerated"]; tcol={"canonical":"#C62828","broad":"#1565C0","degenerated":"#EF6C00"}

fig,(a1,a2)=plt.subplots(1,2,figsize=(13.5,5.6))
groups=["α-sat (HG002)","HSat (HG002)"]+list(per.clade)
enrich={g:{} for g in groups}
for g in ("α-sat (HG002)","HSat (HG002)"):
    for t in tiers: enrich[g][t]=benchg[g][t]
for _,r in per.iterrows():
    for t in tiers: enrich[r.clade][t]=r[f"{t}_enrich"]
x=np.arange(len(groups)); w=0.26
for j,t in enumerate(tiers):
    a1.bar(x+(j-1)*w,[max(enrich[g][t],0.01) for g in groups],w,color=tcol[t],label=t)
a1.axhline(1,ls="--",color="black",lw=1.2); a1.set_yscale("log")
a1.set_xticks(x); a1.set_xticklabels(groups,rotation=20,ha="right")
a1.set_ylabel("enrichment over random (obs / expected)")
a1.set_title("Method 1 — exact IUPAC motif tiers (Fachinetti)\ncanonical = functional box; broad/degenerate = looser",fontweight="bold",fontsize=10.5)
a1.legend(fontsize=9,frameon=False,title="motif"); a1.spines[["top","right"]].set_visible(False)
a1.text(0.01,0.02,"canonical = 0 in every non-human clade",transform=a1.transAxes,fontsize=8,color="grey")

v=df[df.clade=="Vertebrates"].copy(); v=v[(v.broad_perMbp>0)|(v.degenerated_perMbp>0)].sort_values("broad_perMbp")
y=np.arange(len(v))
a2.barh(y+0.2,v.broad_perMbp.clip(lower=0.001),0.4,color=tcol["broad"],label="broad")
a2.barh(y-0.2,v.degenerated_perMbp.clip(lower=0.001),0.4,color=tcol["degenerated"],label="degenerated")
a2.set_yticks(y); a2.set_yticklabels([f"{r['name']} ({r.vgroup})" for _,r in v.iterrows()],fontsize=8)
a2.set_xscale("log"); a2.set_xlabel("motif hits per Mbp (exact IUPAC, log)")
a2.set_title("Vertebrates: exact broad/degenerate hits\n(trace-level; no canonical box)",fontweight="bold",fontsize=10.5)
a2.legend(fontsize=9,frameon=False); a2.spines[["top","right"]].set_visible(False)
plt.tight_layout()
for ext in ("png","pdf"):
    fig.savefig(SAT/f"figures/cenpb_paper_motifs.{ext}",dpi=300 if ext=="png" else None,bbox_inches="tight",facecolor="white")
print("Saved figures/cenpb_paper_motifs.png/pdf")
