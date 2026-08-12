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
# restrict to the 325 published species (tree tips)
ALLOW=set(l.strip() for l in open(SAT/"species_325.txt") if l.strip() and not l.startswith("#"))
want &= ALLOW; vert=vert[vert.species.isin(ALLOW)]

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

def count_from_fasta(path, cap=30000):
    counts=np.zeros((27,4)); nwin=0; bp=0; buf=[]; seqs=[]
    for ln in open(path):
        if ln[0]==">":
            if buf: seqs.append("".join(buf)); buf=[]
        else: buf.append(ln.strip().upper())
    if buf: seqs.append("".join(buf))
    for s in seqs:
        bp+=len(s); both=s+"#"+rc(s)
        if not PRE.search(both): continue
        for st in (s,rc(s)):
            for m in PAT.finditer(st,overlapped=False):
                a,b=m.start(),m.end()
                if b-a!=17 or a-FLANK<0 or b+FLANK>len(st): continue
                w=st[a-FLANK:b+FLANK]
                if len(w)==27 and set(w)<=set("ACGT"):
                    for k,ch in enumerate(w): counts[k,IDX[ch]]+=1
                    nwin+=1
        if nwin>=cap: break
    return counts, nwin, bp

def draw(pdf, counts, title):
    prob=pd.DataFrame(counts/counts.sum(1,keepdims=True),columns=list("ACGT"))
    info=logomaker.transform_matrix(prob,from_type="probability",to_type="information")
    box=info.iloc[FLANK:FLANK+17].sum(1).mean(); fl=pd.concat([info.iloc[:FLANK],info.iloc[FLANK+17:]]).sum(1).mean()
    fig,ax=plt.subplots(figsize=(8.5,2.6)); logomaker.Logo(info,ax=ax,color_scheme="classic",show_spines=False)
    ax.axvline(FLANK-0.5,color="grey",ls="--",lw=1); ax.axvline(FLANK+17-0.5,color="grey",ls="--",lw=1)
    ax.set_ylim(0,2); ax.set_ylabel("bits"); ax.set_xticks([])
    ax.set_title(title.format(box=round(box,2),flank=round(fl,2),delta=round(box-fl,2)),fontsize=8,fontweight="bold")
    plt.tight_layout(); pdf.savefig(fig,bbox_inches="tight"); plt.close(fig)

pdf=PdfPages(SAT/"figures/cenpb_box_logos_flanks_VERTEBRATES_uncapped.pdf")
# --- human HG002 benchmarks first: alpha-sat (positive), HSat (negative) ---
ca,na,bpa=count_from_fasta(SAT/"cenpb_psi/human_alpha.fasta")
draw(pdf,ca,f"HUMAN α-satellite (HG002, positive benchmark) | windows={{}} n={na} ({round(na/(bpa/1e6),0)}/Mbp) | box {{box}} vs flank {{flank}} bits (Δ={{delta:+.2f}})".replace("{}",str(na)))
ch,nh,bph=count_from_fasta(SAT/"cenpb_psi/human_hsat.fasta")
if nh>=10:
    draw(pdf,ch,f"HUMAN HSat1/2/3 (HG002, negative control) | n={nh} ({round(nh/(bph/1e6),1)}/Mbp) | box {{box}} vs flank {{flank}} bits (Δ={{delta:+.2f}})")
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
