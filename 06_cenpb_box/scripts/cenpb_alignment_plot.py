#!/usr/bin/env python3
"""Why the goshawk scores 0 by exact IUPAC matching yet matches by the songbird
≤5-substitution method: its 2 substitutions land exactly on the 5' CpG that all
three IUPAC tiers hold fixed (TTCG anchor). Alignment panel for the report."""
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
SAT=Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity")
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"savefig.dpi":600})

rows=[("canonical (IUPAC)",  "YTTCGTTGGAARCGGGA"),
      ("broad (IUPAC)",      "YTTCGNNNNANRCGGGN"),
      ("degenerate (IUPAC)", "NTTCGNNNNANNCGGGN"),
      ("goshawk consensus",  "CTTTTTTGGAAACGGGA")]
canon=[set("CT"),{"T"},{"T"},{"C"},{"G"},{"T"},{"T"},{"G"},{"G"},{"A"},{"A"},set("AG"),{"C"},{"G"},{"G"},{"G"},{"A"}]
L=17; n=len(rows)
fig,ax=plt.subplots(figsize=(7.8,3.1))
for r,(name,seq) in enumerate(rows):
    y=n-1-r
    for i,ch in enumerate(seq):
        if name.startswith("goshawk"):
            ok = ch in canon[i]
            fc,tc = ("#C8E6C9","#1B5E20") if ok else ("#EF5350","#7f0000")
        elif ch=="N": fc,tc="#F2F2F2","0.55"
        elif ch in "YR": fc,tc="#E1BEE7","#4A148C"
        else: fc,tc="#BBDEFB","#0D47A1"
        ax.add_patch(plt.Rectangle((i,y),1,1,facecolor=fc,edgecolor="white",lw=1.2))
        ax.text(i+0.5,y+0.5,ch,ha="center",va="center",fontsize=9,fontweight="bold",color=tc,family="monospace")
    ax.text(-0.4,y+0.5,name,ha="right",va="center",fontsize=8)
# highlight the eroded CpG (positions 4–5 = index 3–4)
ax.add_patch(plt.Rectangle((3,-0.05),2,n+0.1,fill=False,edgecolor="#C62828",lw=2.2,zorder=5))
ax.annotate("eroded CpG (pos 4–5): C,G → T,T\nCENP-B contact base · methylation/deamination site",
            xy=(4,n+0.05),xytext=(4,n+0.55),ha="center",fontsize=7.2,color="#C62828",fontweight="bold")
# shared anchors bracket
for (a,b,lab) in [(1,5,"TTCG anchor"),(12,16,"CGGG anchor")]:
    ax.plot([a,b],[-0.35,-0.35],color="0.35",lw=1.2); ax.text((a+b)/2,-0.62,lab,ha="center",va="top",fontsize=6.5,color="0.35")
ax.text(-0.4,-0.48,"fixed in ALL IUPAC tiers →",ha="right",va="center",fontsize=6.6,color="0.35",style="italic")
# verdicts
ax.text(L/2,n+1.15,"Same divergence, opposite verdicts",ha="center",fontsize=10,fontweight="bold")
ax.text(L+0.4,n-1.0,"IUPAC (position-based):\nCpG anchor is fixed → 0 matches\n(canonical=broad=degenerate=0)",
        ha="left",va="center",fontsize=7.2,color="#C62828",
        bbox=dict(boxstyle="round,pad=0.4",fc="#FDECEA",ec="#C62828",lw=0.8))
ax.text(L+0.4,n-3.0,"Songbird (≤5 substitutions):\nno fixed position → 1,722 windows\n15/17 identity · flank Δ 0.85",
        ha="left",va="center",fontsize=7.2,color="#1B5E20",
        bbox=dict(boxstyle="round,pad=0.4",fc="#EAF5EA",ec="#1B5E20",lw=0.8))
ax.set_xlim(-4.2,L+7.2); ax.set_ylim(-0.9,n+1.5); ax.axis("off")
plt.tight_layout()
for ext in ("png","pdf"):
    fig.savefig(SAT/f"figures/cenpb_goshawk_alignment.{ext}",bbox_inches="tight",facecolor="white",dpi=600 if ext=="png" else None)
print("Saved figures/cenpb_goshawk_alignment.png/pdf")
