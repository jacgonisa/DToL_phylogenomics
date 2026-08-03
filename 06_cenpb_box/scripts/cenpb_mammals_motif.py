#!/usr/bin/env python3
"""Per-mammal-species CENP-B box motif: extract substitution-only 17-mer sites, build a
PWM per species, report consensus + total information content (real box = high IC,
canonical-matching; noise = low IC). Combined sequence-logo figure."""
import re, regex, collections, numpy as np
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
import pandas as pd
SAT="/home/jg2070/Desktop/dtol_review_August/2026_trees/all.satellites.txt"
BASE="/home/jg2070/Desktop/dtol_review_August/2026_trees/annotation_centromeres"
OUT=Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity/figures")
CORE="[CT]TTCGTTGGAA[AG]CGGGA"; P=regex.compile("(?:%s){s<=3}"%CORE)
comp=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
CANON="YTTCGTTGGAARCGGGA"

xl=pd.read_excel(BASE+"/DTOL_327_master_March.xlsx")
xl["code"]=xl["fasta"].astype(str).str.lower().str.replace(r"[0-9].*$","",regex=True).str.replace(r"[.].*$","",regex=True)
mam=set(xl["code"][xl["taxa1"]=="Mammalia"]); genus=dict(zip(xl["code"],xl.get("genus",xl["taxa1"])))
CAP=8000
sites=collections.defaultdict(list); nrec=collections.Counter()
with open(SAT) as fh:
    next(fh)
    for line in fh:
        f=line.split()
        if len(f)<15: continue
        sp=re.sub(r'\.\d+$','',f[14].strip('"')).lower()
        if sp not in mam: continue
        nrec[sp]+=1
        if len(sites[sp])>=CAP: continue
        s=f[11].strip('"').upper()
        for seq in (s,rc(s)):
            m=P.search(seq)
            if m and len(m.group())==17 and set(m.group())<=set("ACGT"): sites[sp].append(m.group())

idx={'A':0,'C':1,'G':2,'T':3}
def pwm(ss):
    M=np.zeros((4,17))
    for s in ss:
        for i,ch in enumerate(s): M[idx[ch],i]+=1
    f=M/M.sum(0); ic=2+(f*np.log2(f+1e-9)).sum(0)
    return f,ic,"".join("ACGT"[c] for c in f.argmax(0))

rows=[]
for sp in mam:
    ss=sites[sp]
    if len(ss)>=30:
        f,ic,cons=pwm(ss); rows.append((sp,genus.get(sp,''),nrec[sp],len(ss),ic.sum(),cons,f,ic))
    else:
        rows.append((sp,genus.get(sp,''),nrec[sp],len(ss),0,"-",None,None))
rows.sort(key=lambda r:-r[4])
print(f"{'sp':9s} {'genus':13s} {'n_rec':>9} {'sites':>6} {'IC_bits':>7}  consensus")
print(f"{'canonical':9s} {'':13s} {'':>9} {'':>6} {'34.0':>7}  {CANON}")
for sp,g,nr,ns,ics,cons,_,_ in rows:
    print(f"{sp:9s} {str(g):13s} {nr:>9,} {ns:>6} {ics:>7.1f}  {cons}")
pd.DataFrame([(sp,g,nr,ns,round(ics,1),cons) for sp,g,nr,ns,ics,cons,_,_ in rows],
    columns=["species","genus","n_records","n_sites","total_IC_bits","consensus"]).to_csv(OUT/"cenpb_mammals_motif.tsv",sep="\t",index=False)

# combined logo for species with a real motif (>=30 sites)
COL={'A':'#2ca02c','C':'#1f77b4','G':'#ff7f0e','T':'#d62728'}
real=[r for r in rows if r[6] is not None]
fig,axes=plt.subplots(len(real),1,figsize=(9,1.5*len(real)),squeeze=False); axes=axes[:,0]
for ax,(sp,g,nr,ns,ics,cons,f,ic) in zip(axes,real):
    for i in range(17):
        y=0
        for b in np.argsort(f[:,i]):
            base="ACGT"[b]; h=f[b,i]*ic[i]
            if h<0.01: continue
            tp=TextPath((0,0),base,size=1); bb=tp.get_extents()
            t=Affine2D().translate(-bb.x0,-bb.y0).scale(1/(bb.width or 1),h/(bb.height or 1)).translate(i,y)
            ax.add_patch(PathPatch(tp.transformed(t),color=COL[base],lw=0)); y+=h
    ax.set_xlim(-.5,17); ax.set_ylim(0,2); ax.set_xticks([])
    ax.set_ylabel(f"{g}\n{ns} sites\n{ics:.0f} bits",rotation=0,ha="right",va="center",fontsize=8)
axes[-1].set_xticks(range(17)); axes[-1].set_xticklabels(range(1,18))
fig.suptitle("Per-mammal CENP-B-box motif (high IC = real box)",fontsize=11)
fig.tight_layout(); fig.savefig(OUT/"cenpb_mammals_logos.png",dpi=200)
print("saved cenpb_mammals_motif.tsv, cenpb_mammals_logos.png")
