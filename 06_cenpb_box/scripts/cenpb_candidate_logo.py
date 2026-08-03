#!/usr/bin/env python3
"""Characterise a candidate species' CENP-B-box-like motif: extract 17-mers matching the
broad box (NTTCGNNNNANNCGGGN) from its satellites, build a PWM + information-content logo.
Real conserved box => high IC at the *variable* positions too, not just the fixed core.
Usage: cenpb_candidate_logo.py <species_code>"""
import re, sys, numpy as np
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D; from matplotlib.textpath import TextPath; from matplotlib.patches import PathPatch
SP=sys.argv[1] if len(sys.argv)>1 else "drparjuda"
SAT="/home/jg2070/Desktop/dtol_review_August/2026_trees/all.satellites.txt"
OUT=Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity/figures_full")
BROAD=re.compile("[ACGT]TTCG[ACGT]{4}A[ACGT]{2}CGGG[ACGT]")   # NTTCGNNNNANNCGGGN
comp=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
sites=[]
with open(SAT) as fh:
    next(fh)
    for line in fh:
        f=line.split()
        if len(f)<15 or re.sub(r'\.\d+$','',f[14].strip('"')).lower()!=SP: continue
        s=f[11].strip('"').upper()
        for seq in (s,rc(s)):
            m=BROAD.search(seq)
            if m and set(m.group())<=set("ACGT"): sites.append(m.group())
        if len(sites)>=20000: break
print(f"{SP}: {len(sites)} broad-box sites")
idx={'A':0,'C':1,'G':2,'T':3}; W=17
M=np.zeros((4,W))
for s in sites:
    for i,ch in enumerate(s): M[idx[ch],i]+=1
f=M/M.sum(0); ic=2+(f*np.log2(f+1e-9)).sum(0)
cons="".join("ACGT"[c] for c in f.argmax(0))
print("canonical : YTTCGTTGGAARCGGGA")
print(f"{SP:9s} : {cons}   total IC {ic.sum():.1f} bits, mean {ic.mean():.2f}/pos")
print("per-pos IC:", " ".join(f"{v:.1f}" for v in ic))
COL={'A':'#2ca02c','C':'#1f77b4','G':'#ff7f0e','T':'#d62728'}
fig,ax=plt.subplots(figsize=(9,3))
for i in range(W):
    y=0
    for b in np.argsort(f[:,i]):
        h=f[b,i]*ic[i]
        if h<.01: continue
        tp=TextPath((0,0),"ACGT"[b],size=1); bb=tp.get_extents()
        t=Affine2D().translate(-bb.x0,-bb.y0).scale(1/(bb.width or 1),h/(bb.height or 1)).translate(i,y)
        ax.add_patch(PathPatch(tp.transformed(t),color=COL["ACGT"[b]],lw=0)); y+=h
ax.set_xlim(-.5,W);ax.set_ylim(0,2);ax.set_xticks(range(W));ax.set_xticklabels(range(1,W+1))
ax.set_ylabel("bits");ax.set_title(f"{SP} CENP-B-box-like motif ({len(sites)} sites, {ic.sum():.0f} bits) — consensus {cons}")
fig.tight_layout();fig.savefig(OUT/f"cenpb_candidate_{SP}_logo.png",dpi=200)
print(f"saved cenpb_candidate_{SP}_logo.png")
