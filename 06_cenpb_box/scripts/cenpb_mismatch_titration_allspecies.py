#!/usr/bin/env python3
"""CENP-B box occurrences across ALL species' satellites vs allowed mismatches
(0..5), coloured by clade, with human alpha-satellite (HG002) as the positive
benchmark and a mononucleotide-shuffled null. No PWM/FIMO - just the canonical
17-bp box with increasing mismatch tolerance, both strands."""
import random, collections
from pathlib import Path
import regex, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

random.seed(0)
FASTA = "/mnt/ssd-8tb/satellite_dna_lm/monomers.fasta"     # curated monomers, <=500/species
HUMAN = "/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity/cenpb_psi/human_alpha.fasta"
BASE  = Path("/home/jg2070/Desktop/dtol_review_August")
SAT   = BASE/"DToL_phylogenomics_publication_325genomes/05_satellite_similarity"
CORE  = "[CT]TTCGTTGGAA[AG]CGGGA"; KMAX = 5
PATS  = {k: regex.compile("(?:%s){e<=%d}" % (CORE, k)) for k in range(KMAX+1)}
comp  = str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
hit   = lambda p,s: bool(p.search(s) or p.search(rc(s)))

# species code -> broad clade
tax = pd.read_csv(BASE/"2026_trees/annotation_centromeres/centromere_code_to_species.tsv", sep="\t")
import re as _re
tax["code"]=tax["fasta"].astype(str).str.lower().str.replace(r"[0-9.].*$","",regex=True)
bmap={"Mammalia":"Vertebrates","Aves":"Vertebrates","Actinopterygii":"Vertebrates","Reptilia":"Vertebrates",
      "Amphibia":"Vertebrates","Chondrichthyes":"Vertebrates","Fungi":"Fungi","Algae":"Viridiplantae",
      "Bryophyta":"Viridiplantae","Dicots":"Viridiplantae","Monocots":"Viridiplantae","Gymnosperms":"Viridiplantae",
      "Alveolata":"Protist","Discoba":"Protist"}
clade={c:bmap.get(t,"Invertebrate") for c,t in zip(tax.code,tax.taxa1)}

# load all monomers per species
seqs=collections.defaultdict(list); cur=None; buf=[]
for line in open(FASTA):
    if line[0]==">":
        if cur and buf: seqs[cur].append("".join(buf))
        cur=line[1:].split("__")[0].strip().lower(); buf=[]
    else: buf.append(line.strip().upper())
if cur and buf: seqs[cur].append("".join(buf))
# human benchmark
hs=[]; buf=[]
for line in open(HUMAN):
    if line[0]==">":
        if buf: hs.append("".join(buf)); buf=[]
    else: buf.append(line.strip().upper())
if buf: hs.append("".join(buf))
random.shuffle(hs); seqs["human"]=hs[:500]; clade["human"]="Human"
print(f"species scanned: {len(seqs)}")

rows=[]
for sp,mons in seqs.items():
    if not mons: continue
    shuf=["".join(random.sample(m,len(m))) for m in mons]
    for k in range(KMAX+1):
        o=sum(hit(PATS[k],m) for m in mons); nu=sum(hit(PATS[k],m) for m in shuf)
        rows.append(dict(species=sp, clade=clade.get(sp,"Invertebrate"), n=len(mons), mm=k,
                         obs=o, null=nu, pct_obs=100*o/len(mons), pct_null=100*nu/len(mons)))
df=pd.DataFrame(rows); df.to_csv(SAT/"figures/cenpb_mismatch_titration_allspecies.tsv",sep="\t",index=False)

ccol={"Vertebrates":"#1565C0","Invertebrate":"#EF6C00","Viridiplantae":"#2E7D32",
      "Fungi":"#6A1B9A","Protist":"#C62828"}
fig,(a1,a2)=plt.subplots(1,2,figsize=(13,5.4))
for sp,d in df.groupby("species"):
    d=d.sort_values("mm"); cl=d.clade.iloc[0]
    if cl=="Human": continue
    a1.plot(d.mm,d.obs,color=ccol.get(cl,"grey"),alpha=0.35,lw=0.8)
    a2.plot(d.mm,d.pct_obs,color=ccol.get(cl,"grey"),alpha=0.35,lw=0.8)
# human benchmark (bold red) + mean null (black dashed)
h=df[df.species=="human"].sort_values("mm")
nu=df[df.species!="human"].groupby("mm").agg(null=("null","mean"),pct_null=("pct_null","mean")).reset_index()
for ax,yo,yn,hn in [(a1,"obs","null",h.obs),(a2,"pct_obs","pct_null",h.pct_obs)]:
    ax.plot(h.mm,hn,color="#d62728",lw=3,marker="o",ms=5,label="Human α-sat (benchmark)",zorder=5)
    ax.plot(nu.mm,nu[yn],color="black",ls="--",lw=1.8,label="shuffled null (mean)")
a1.set_ylabel("# monomers with a CENP-B box"); a2.set_ylabel("% monomers with a CENP-B box")
for ax,t in [(a1,"occurrence count (/ ≤500 monomers)"),(a2,"occurrence rate")]:
    ax.set_xlabel("allowed mismatches to the 17-bp CENP-B box"); ax.set_title(t,fontweight="bold")
    ax.set_xticks(range(KMAX+1)); ax.spines[["top","right"]].set_visible(False)
handles=[Line2D([0],[0],color=c,lw=2,label=g) for g,c in ccol.items() if g in df.clade.values]
handles=[Line2D([0],[0],color="#d62728",lw=3,label="Human α-sat (benchmark)")]+handles+\
        [Line2D([0],[0],color="black",ls="--",lw=1.8,label="shuffled null")]
a2.legend(handles=handles,fontsize=8.5,frameon=False)
fig.suptitle("CENP-B box across all DToL satellites vs mismatch tolerance (canonical 17-bp box, both strands)\n"
             "Human α-satellite = positive benchmark; each thin line = one species",fontweight="bold")
plt.tight_layout()
for ext in ("png","pdf"):
    fig.savefig(SAT/f"figures/cenpb_mismatch_titration_allspecies.{ext}",dpi=300 if ext=="png" else None,
                bbox_inches="tight",facecolor="white")
    print("Saved:",SAT/f"figures/cenpb_mismatch_titration_allspecies.{ext}")
# any non-human species with a box at <=2 mismatches above its null?
print("\nspecies with obs>null at mm<=2 (excess>=3):")
d2=df[(df.mm<=2)&(df.species!="human")].copy(); d2["ex"]=d2.obs-d2.null
print(d2[d2.ex>=3].sort_values("ex",ascending=False)[["species","clade","mm","obs","null"]].to_string(index=False) or "  none")
