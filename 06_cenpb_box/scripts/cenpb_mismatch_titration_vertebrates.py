#!/usr/bin/env python3
"""CENP-B box occurrences in vertebrate satellites as a function of allowed
mismatches (0..5), vs a mononucleotide-shuffled null. No PWM/FIMO — just the
canonical 17-bp box with increasing mismatch tolerance, both strands."""
import re, random, collections
from pathlib import Path
import regex, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

random.seed(0)
FASTA = "/mnt/ssd-8tb/satellite_dna_lm/monomers.fasta"
SAT   = Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity")
CORE  = "[CT]TTCGTTGGAA[AG]CGGGA"                      # canonical 17-bp CENP-B box
KMAX  = 5
PATS  = {k: regex.compile("(?:%s){e<=%d}" % (CORE, k)) for k in range(KMAX+1)}
comp  = str.maketrans("ACGTN","TGCAN"); rc = lambda s: s.translate(comp)[::-1]
def hit(pat, s): return bool(pat.search(s) or pat.search(rc(s)))

# vertebrate species -> subgroup (from the scanned per-species table)
tab = pd.read_csv(SAT/"figures/cenpb_box_per_species.tsv", sep="\t")
vgrp = {"Mammalia":"Mammalia","Aves":"Aves","Actinopterygii":"Fish",
        "Reptilia":"Reptilia","Amphibia":"Amphibia","Chondrichthyes":"Fish"}
vspec = {r.species: vgrp[r.taxa1] for r in tab.itertuples() if r.taxa1 in vgrp}

# load monomers per vertebrate species
seqs = collections.defaultdict(list); cur=None; buf=[]
for line in open(FASTA):
    if line[0]==">":
        if cur in vspec and buf: seqs[cur].append("".join(buf))
        cur = line[1:].split("__")[0].strip().lower(); buf=[]
    else: buf.append(line.strip().upper())
if cur in vspec and buf: seqs[cur].append("".join(buf))

# ---- human alpha-satellite (HG002) as positive-control benchmark ----
HUMAN = SAT/"cenpb_psi/human_alpha.fasta"
hs=[]; buf=[]
for line in open(HUMAN):
    if line[0]==">":
        if buf: hs.append("".join(buf)); buf=[]
    else: buf.append(line.strip().upper())
if buf: hs.append("".join(buf))
random.shuffle(hs); seqs["human"] = hs[:500]      # match the <=500/species cap
vspec["human"] = "Human"

rows=[]
for sp, mons in seqs.items():
    if not mons: continue
    shuf=["".join(random.sample(m,len(m))) for m in mons]        # mononucleotide-shuffled null
    for k in range(KMAX+1):
        obs=sum(hit(PATS[k],m) for m in mons)
        nul=sum(hit(PATS[k],m) for m in shuf)
        rows.append(dict(species=sp, group=vspec[sp], n=len(mons),
                         mm=k, obs=obs, null=nul,
                         pct_obs=100*obs/len(mons), pct_null=100*nul/len(mons)))
df=pd.DataFrame(rows)
df.to_csv(SAT/"figures/cenpb_mismatch_titration_vertebrates.tsv", sep="\t", index=False)

# ---- plot: occurrences vs mismatches, one line per species, coloured by group ----
gcol={"Mammalia":"#1565C0","Aves":"#00838F","Fish":"#2E7D32","Reptilia":"#6A1B9A","Amphibia":"#EF6C00"}
fig,(a1,a2)=plt.subplots(1,2,figsize=(13,5.2))
for sp,d in df.groupby("species"):
    d=d.sort_values("mm"); g=d.group.iloc[0]
    a1.plot(d.mm,d.obs,color=gcol[g],alpha=0.55,lw=1)
    a2.plot(d.mm,d.pct_obs,color=gcol[g],alpha=0.55,lw=1)
# mean null across vertebrates
nu=df.groupby("mm").agg(null=("null","mean"),pct_null=("pct_null","mean")).reset_index()
a1.plot(nu.mm,nu.null,color="black",ls="--",lw=1.8,label="shuffled null (mean)")
a2.plot(nu.mm,nu.pct_null,color="black",ls="--",lw=1.8,label="shuffled null (mean)")
for ax,ylab,ttl in [(a1,"# monomers with a CENP-B box","occurrence count"),
                    (a2,"% monomers with a CENP-B box","occurrence rate")]:
    ax.set_xlabel("allowed mismatches to the 17-bp box"); ax.set_ylabel(ylab)
    ax.set_title(ttl,fontweight="bold"); ax.set_xticks(range(KMAX+1))
    ax.spines[["top","right"]].set_visible(False)
from matplotlib.lines import Line2D
handles=[Line2D([0],[0],color=c,lw=2,label=g) for g,c in gcol.items() if g in df.group.values]
handles.append(Line2D([0],[0],color="black",ls="--",lw=1.8,label="shuffled null"))
a2.legend(handles=handles,fontsize=9,frameon=False)
fig.suptitle("CENP-B box in vertebrate satellites vs mismatch tolerance (canonical 17-bp box, both strands)",
             fontweight="bold")
plt.tight_layout()
for ext in ("png","pdf"):
    fig.savefig(SAT/f"figures/cenpb_mismatch_titration_vertebrates.{ext}",dpi=300 if ext=="png" else None,
                bbox_inches="tight",facecolor="white")
    print("Saved:",SAT/f"figures/cenpb_mismatch_titration_vertebrates.{ext}")
# quick summary: obs vs null at each mm, and the standout species
print(df.groupby("mm").agg(obs=("obs","sum"),null=("null","sum")).to_string())
print("\nTop species at mm=2 (obs-null):")
d2=df[df.mm==2].copy(); d2["excess"]=d2.obs-d2.null
print(d2.sort_values("excess",ascending=False)[["species","group","obs","null"]].head(8).to_string(index=False))
