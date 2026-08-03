#!/usr/bin/env python3
"""Focused CENP-B box screen of the 8 mammal species (full satellite records, both
strands): exact canonical, <=2 and <=3 mismatch, and the 10-bp core, per species."""
import re, regex, collections
import pandas as pd
SAT="/home/jg2070/Desktop/dtol_review_August/2026_trees/all.satellites.txt"
BASE="/home/jg2070/Desktop/dtol_review_August/2026_trees/annotation_centromeres"
CORE="[CT]TTCGTTGGAA[AG]CGGGA"
EXACT=re.compile(CORE); P2=regex.compile("(?:%s){e<=2}"%CORE); P3=regex.compile("(?:%s){e<=3}"%CORE)
C10=re.compile("TTCGTTGGAA")
comp=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
def h(p,s): return bool(p.search(s) or p.search(rc(s)))

xl=pd.read_excel(BASE+"/DTOL_327_master_March.xlsx")
xl["code"]=xl["fasta"].astype(str).str.lower().str.replace(r"[0-9].*$","",regex=True).str.replace(r"[.].*$","",regex=True)
mam=set(xl["code"][xl["taxa1"]=="Mammalia"]); genus=dict(zip(xl["code"],xl.get("genus",xl["taxa1"])))

C=collections.defaultdict(collections.Counter)
with open(SAT) as fh:
    next(fh)
    for i,line in enumerate(fh):
        if i%5_000_000==0: print(f"  {i:,}",flush=True)
        f=line.split()
        if len(f)<15: continue
        sp=re.sub(r'\.\d+$','',f[14].strip('"')).lower()
        if sp not in mam: continue
        s=f[11].strip('"').upper()
        if len(s)<17: continue
        c=C[sp]; c["n"]+=1
        if h(EXACT,s): c["exact"]+=1
        if h(P2,s):    c["le2"]+=1
        if h(P3,s):    c["le3"]+=1
        if h(C10,s):   c["core"]+=1
print(f"\n{'species':10s} {'genus':13s} {'n':>10} {'exact':>6} {'<=2mm':>6} {'<=3mm':>6} {'core10':>7}")
tot=collections.Counter()
for sp in sorted(C,key=lambda x:-C[x]['n']):
    c=C[sp]
    for k in c: tot[k]+=c[k]
    print(f"{sp:10s} {str(genus.get(sp,'')):13s} {c['n']:>10,} {c['exact']:>6} {c['le2']:>6} {c['le3']:>6} {c['core']:>7}")
print(f"\nMAMMAL TOTAL n={tot['n']:,}  exact={tot['exact']}  <=2mm={tot['le2']}  <=3mm={tot['le3']}  core10={tot['core']}")
