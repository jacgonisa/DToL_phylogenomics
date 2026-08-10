#!/usr/bin/env python3
"""Build box +/-5 flank sequence logos for the vertebrates (uncapped arrays),
using the per-species table from cenpb_flank_uncapped.py to pick which species
to draw. Re-extracts windows for vertebrate species only (fast subset), caps at
CAP windows/species for a clean logo, and compiles one PDF ordered by Delta."""
import regex, re, collections, numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import logomaker

BASE=Path("/home/jg2070/Desktop/dtol_review_August")
SAT =BASE/"DToL_phylogenomics_publication_325genomes/05_satellite_similarity"
ALL =BASE/"2026_trees/all.satellites.txt"
CORE="[CT]TTCGTTGGAA[AG]CGGGA"; FLANK=5; CAP=20000
PRE =regex.compile("(?:TTCGTTGGAA){s<=3}"); PAT=regex.compile("(?:%s){s<=5}"%CORE)
comp=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
seqre=re.compile(r'"([ACGTNacgtn]{25,})"'); IDX={b:i for i,b in enumerate("ACGT")}

df=pd.read_csv(SAT/"figures/cenpb_flank_uncapped_per_species.tsv",sep="\t")
vert=df[(df.clade=="Vertebrates")&(df.n_windows>=5)].sort_values("delta",ascending=False)
want=set(vert.species); counts={sp:np.zeros((27,4)) for sp in want}; got={sp:0 for sp in want}
nm=dict(zip(df.species,df.name)); vg=dict(zip(df.species,df.vgroup))

for i,line in enumerate(open(ALL)):
    if i==0: continue
    if '"' not in line: continue
    sp=line.rstrip().rsplit('"',2)[-2].split(".")[0].lower()
    if sp not in want or got[sp]>=CAP: continue
    m=seqre.search(line)
    if not m: continue
    s=m.group(1).upper(); both=s+"#"+rc(s)
    if not PRE.search(both): continue
    for st in (s,rc(s)):
        for mm in PAT.finditer(st,overlapped=False):
            a,b=mm.start(),mm.end()
            if b-a!=17 or a-FLANK<0 or b+FLANK>len(st): continue
            w=st[a-FLANK:b+FLANK]
            if len(w)==27 and set(w)<=set("ACGT"):
                for k,ch in enumerate(w): counts[sp][k,IDX[ch]]+=1
                got[sp]+=1
                if got[sp]>=CAP: break

pdf=PdfPages(SAT/"figures/cenpb_box_logos_flanks_VERTEBRATES_uncapped.pdf")
for _,r in vert.iterrows():
    sp=r.species; c=counts[sp]
    if c.sum()<5*10: continue
    prob=pd.DataFrame(c/c.sum(1,keepdims=True),columns=list("ACGT"))
    info=logomaker.transform_matrix(prob,from_type="probability",to_type="information")
    fig,ax=plt.subplots(figsize=(8.5,2.6)); logomaker.Logo(info,ax=ax,color_scheme="classic",show_spines=False)
    ax.axvline(FLANK-0.5,color="grey",ls="--",lw=1); ax.axvline(FLANK+17-0.5,color="grey",ls="--",lw=1)
    ax.set_ylim(0,2); ax.set_ylabel("bits"); ax.set_xticks([])
    ax.set_title(f"{r['name']} ({sp}) | {vg.get(sp,'')} | windows={int(got[sp])} ({r.win_per_Mbp}/Mbp) | "
                 f"box {r.mean_box_bits} vs flank {r.mean_flank_bits} bits (Δ={r.delta:+.2f}) | "
                 f"consensus {r.box_consensus} ({r.subs_vs_canonical} subs)",fontsize=8,fontweight="bold")
    plt.tight_layout(); pdf.savefig(fig,bbox_inches="tight"); plt.close(fig)
pdf.close()
print("Saved:",SAT/"figures/cenpb_box_logos_flanks_VERTEBRATES_uncapped.pdf")
print(f"vertebrate logos drawn: {len(vert)}")
