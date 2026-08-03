#!/usr/bin/env python3
"""Characterise the horseshoe-bat CENP-B-box motif: extract substitution-only 17-mer
sites from its satellites, build a PWM + information-content sequence logo, and write an
HMMER-ready seed alignment (for nhmmer scans of other species)."""
import re, regex, numpy as np
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
SAT="/home/jg2070/Desktop/dtol_review_August/2026_trees/all.satellites.txt"
OUT=Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity/figures")
CORE="[CT]TTCGTTGGAA[AG]CGGGA"; P=regex.compile("(?:%s){s<=3}"%CORE)
comp=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
BAT="mrhisin"

sites=[]
with open(SAT) as fh:
    next(fh)
    for line in fh:
        f=line.split()
        if len(f)<15 or re.sub(r'\.\d+$','',f[14].strip('"')).lower()!=BAT: continue
        s=f[11].strip('"').upper()
        for seq in (s, rc(s)):
            m=P.search(seq)
            if m and len(m.group())==17 and set(m.group())<=set("ACGT"): sites.append(m.group())
        if len(sites)>=20000: break
print(f"{len(sites)} sites")

idx={'A':0,'C':1,'G':2,'T':3}
M=np.zeros((4,17))
for s in sites:
    for i,ch in enumerate(s): M[idx[ch],i]+=1
freq=M/M.sum(0); ic=2+(freq*np.log2(freq+1e-9)).sum(0)
cons="".join("ACGT"[c] for c in freq.argmax(0))
np.savetxt(OUT/"cenpb_bat_pwm.tsv",freq,delimiter="\t",header="pos1..17 rows=A,C,G,T",fmt="%.4f")
# seed alignment for hmmbuild (subsample)
with open(OUT/"cenpb_bat_seed.fa","w") as fh:
    for i,s in enumerate(sites[:500]): fh.write(f">s{i}\n{s}\n")

# information-content sequence logo
COL={'A':'#2ca02c','C':'#1f77b4','G':'#ff7f0e','T':'#d62728'}
fig,ax=plt.subplots(figsize=(9,3))
for i in range(17):
    order=np.argsort(freq[:,i]); y=0
    for b in order:
        base="ACGT"[b]; h=freq[b,i]*ic[i]
        if h<0.01: y+=h; continue
        tp=TextPath((0,0),base,size=1,prop=None)
        bb=tp.get_extents(); sx=1/(bb.width or 1); sy=h/(bb.height or 1)
        t=Affine2D().translate(-bb.x0,-bb.y0).scale(sx,sy).translate(i,y)
        ax.add_patch(PathPatch(tp.transformed(t),color=COL[base],lw=0)); y+=h
ax.set_xlim(-.5,17); ax.set_ylim(0,2); ax.set_xticks(range(17)); ax.set_xticklabels(range(1,18))
ax.set_ylabel("bits"); ax.set_title(f"Horseshoe-bat CENP-B-box motif ({len(sites)} sites)  consensus {cons}")
fig.tight_layout(); fig.savefig(OUT/"cenpb_bat_logo.png",dpi=200)
print("consensus",cons,"| total IC %.1f bits"%ic.sum())
print("saved cenpb_bat_logo.png, cenpb_bat_pwm.tsv, cenpb_bat_seed.fa")
