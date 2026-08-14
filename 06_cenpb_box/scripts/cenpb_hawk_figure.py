#!/usr/bin/env python3
"""Two-panel goshawk figure:
 (A) the CENP-B box definitions (canonical/broad/degenerate IUPAC) aligned to the
     goshawk satellite consensus — the 2 substitutions land on the 5' CpG all tiers fix.
 (B) the goshawk motif itself as a sequence logo (its box-like 17-bp windows + ±5 flanks):
     conserved box, random flanks, with the eroded CpG visible.
"""
import regex, re, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import logomaker
from pathlib import Path
SAT=Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity")
CORE="[CT]TTCGTTGGAA[AG]CGGGA"; FLANK=5
PRE=regex.compile("(?:TTCGTTGGAA){s<=3}"); PAT=regex.compile("(?:%s){s<=5}"%CORE)
comp=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
seqre=re.compile(r'"([ACGTNacgtn]{25,})"'); IDX={b:i for i,b in enumerate("ACGT")}

# ---- build the goshawk box±5 PFM from its arrays ----
pfm=np.zeros((27,4)); nwin=0
for line in open("/tmp/goshawk.txt"):
    m=seqre.search(line)
    if not m: continue
    s=m.group(1).upper()
    if not PRE.search(s+"#"+rc(s)): continue
    for st in (s,rc(s)):
        for mm in PAT.finditer(st,overlapped=False):
            a,b=mm.start(),mm.end()
            if b-a!=17 or a-FLANK<0 or b+FLANK>len(st): continue
            w=st[a-FLANK:b+FLANK]
            if len(w)==27 and set(w)<=set("ACGT"):
                for k,ch in enumerate(w): pfm[k,IDX[ch]]+=1
                nwin+=1
print("goshawk windows:",nwin)
import pandas as pd
prob=pd.DataFrame(pfm/pfm.sum(1,keepdims=True),columns=list("ACGT"))
info=logomaker.transform_matrix(prob,from_type="probability",to_type="information")
cons="".join("ACGT"[pfm[i+FLANK].argmax()] for i in range(17))
print("goshawk box consensus:",cons)

# ---- figure ----
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"savefig.dpi":600})
fig=plt.figure(figsize=(8.2,5.2))
gs=fig.add_gridspec(2,1,height_ratios=[1.15,1.0],hspace=0.42)
# ---- Panel A: alignment ----
axA=fig.add_subplot(gs[0]); axA.axis("off")
rows=[("canonical (IUPAC)","YTTCGTTGGAARCGGGA"),("broad (IUPAC)","YTTCGNNNNANRCGGGN"),
      ("degenerate (IUPAC)","NTTCGNNNNANNCGGGN"),("bAccGen1.1 consensus",cons)]
canon=[set("CT"),{"T"},{"T"},{"C"},{"G"},{"T"},{"T"},{"G"},{"G"},{"A"},{"A"},set("AG"),{"C"},{"G"},{"G"},{"G"},{"A"}]
nR=len(rows)
for r,(nm,seq) in enumerate(rows):
    y=nR-1-r
    for i,ch in enumerate(seq):
        if nm.endswith("consensus"):
            ok=ch in canon[i]; fc,tc=("#C8E6C9","#1B5E20") if ok else ("#EF5350","#7f0000")
        elif ch=="N": fc,tc="#F2F2F2","0.55"
        elif ch in "YR": fc,tc="#E1BEE7","#4A148C"
        else: fc,tc="#BBDEFB","#0D47A1"
        axA.add_patch(plt.Rectangle((i,y),1,1,facecolor=fc,edgecolor="white",lw=1.1))
        axA.text(i+0.5,y+0.5,ch,ha="center",va="center",fontsize=8.5,fontweight="bold",color=tc,family="monospace")
    axA.text(-0.4,y+0.5,nm,ha="right",va="center",fontsize=8)
axA.add_patch(plt.Rectangle((3,-0.05),2,nR+0.1,fill=False,edgecolor="#C62828",lw=2.2,zorder=5))
axA.annotate("eroded CpG (pos 4–5: C,G → T,T)",xy=(4,nR),xytext=(4,nR+0.5),ha="center",fontsize=7.4,color="#C62828",fontweight="bold")
axA.plot([1,5],[-0.35,-0.35],color="0.35",lw=1.1); axA.text(3,-0.62,"TTCG anchor",ha="center",va="top",fontsize=6.6,color="0.35")
axA.plot([12,16],[-0.35,-0.35],color="0.35",lw=1.1); axA.text(14,-0.62,"CGGG anchor",ha="center",va="top",fontsize=6.6,color="0.35")
axA.set_xlim(-4,17.3); axA.set_ylim(-0.95,nR+1.0)
axA.set_title("A   Box definitions vs the bAccGen1.1 (Accipiter gentilis) satellite motif",fontsize=10,fontweight="bold",loc="left")
# ---- Panel B: goshawk logo ----
axB=fig.add_subplot(gs[1])
logomaker.Logo(info,ax=axB,color_scheme="classic",show_spines=False)
axB.axvline(FLANK-0.5,color="grey",ls="--",lw=1); axB.axvline(FLANK+17-0.5,color="grey",ls="--",lw=1)
axB.axvspan(FLANK+3-0.5,FLANK+5-0.5,color="#C62828",alpha=0.10)     # eroded CpG columns
axB.set_ylim(0,2); axB.set_ylabel("bits"); axB.set_xticks([])
axB.text(FLANK+8,-0.24,"17-bp box (conserved)",ha="center",va="top",fontsize=7.5,color="0.3")
axB.text(2,-0.24,"−5 flank",ha="center",va="top",fontsize=7,color="0.45"); axB.text(24,-0.24,"+5 flank",ha="center",va="top",fontsize=7,color="0.45")
axB.set_title(f"B   bAccGen1.1 motif logo (n={nwin:,} box-like windows; consensus {cons}, 2 subs from canonical)",fontsize=10,fontweight="bold",loc="left")
for ext in ("png","pdf"):
    fig.savefig(SAT/f"figures/cenpb_hawk_box_AB.{ext}",dpi=600 if ext=="png" else None,bbox_inches="tight",facecolor="white")
print("Saved figures/cenpb_hawk_box_AB.png/pdf")
