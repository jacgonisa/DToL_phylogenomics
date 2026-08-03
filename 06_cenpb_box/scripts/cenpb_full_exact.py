#!/usr/bin/env python3
"""Exact full-data screen (all satellite records, both strands) for the canonical
CENP-B box and its 10-bp conserved core. Fast (no fuzzy). Settles whether boxes exist
anywhere in the curated satellites."""
import re, collections, sys
import pandas as pd
SAT="/home/jg2070/Desktop/dtol_review_August/2026_trees/all.satellites.txt"
BASE="/home/jg2070/Desktop/dtol_review_August/2026_trees/annotation_centromeres"
BOX=re.compile("[CT]TTCGTTGGAA[AG]CGGGA")     # canonical 17-bp
CORE=re.compile("TTCGTTGGAA")                  # conserved core
comp=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
xl=pd.read_excel(BASE+"/DTOL_327_master_March.xlsx")
xl["code"]=xl["fasta"].astype(str).str.lower().str.replace(r"[0-9].*$","",regex=True).str.replace(r"[.].*$","",regex=True)
sp2t=dict(zip(xl["code"],xl["taxa1"]))
box=collections.Counter(); core=collections.Counter(); nrec=collections.Counter()
with open(SAT) as fh:
    next(fh)
    for i,line in enumerate(fh):
        if i%5_000_000==0: print(f"  {i:,}",flush=True)
        f=line.split()
        if len(f)<15: continue
        sp=re.sub(r'\.\d+$','',f[14].strip('"')).lower(); t=sp2t.get(sp,"NA")
        s=f[11].strip('"').upper()
        if len(s)<17: continue
        nrec[t]+=1
        if BOX.search(s) or BOX.search(rc(s)): box[t]+=1
        if CORE.search(s) or CORE.search(rc(s)): core[t]+=1
print(f"\n{'taxa1':16s} {'n_records':>10} {'canonical_box':>14} {'core10':>8}")
for t in sorted(nrec,key=lambda x:-nrec[x]):
    print(f"{t:16s} {nrec[t]:>10,} {box[t]:>14,} {core[t]:>8,}")
print(f"\nTOTAL records={sum(nrec.values()):,}  canonical boxes={sum(box.values()):,}  core10={sum(core.values()):,}")
