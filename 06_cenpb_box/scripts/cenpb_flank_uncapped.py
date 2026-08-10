#!/usr/bin/env python3
"""Genome-wide CENP-B box +/-5 flank enrichment test for ALL species, on the
UNCAPPED satellite arrays (all.satellites.txt, ~24.5M arrays) -- not the <=500
curated monomers. Implements the songbird paper's negative control (Formenti et
al., Cell 2026, Suppl. Fig. 15): boxes are matched as fixed 17-bp windows
(<=5 substitutions, both strands), and for each species we accumulate the 27-bp
(5 + 17 + 5) position frequency matrix. A genuine CENP-B box shows high
information across the box that collapses in the flanks (Delta = box-flank > 0);
its consensus should also resemble the canonical box. We report, per species:
n_arrays, n_windows, prevalence, box vs flank bits, Delta, box consensus and its
substitutions from canonical.

CAVEAT: substitution-matched windows are constrained toward the consensus, which
inflates box bits by construction; Delta alone is not proof. The discriminators
are Delta AND consensus resemblance AND prevalence, reported together."""
import regex, re, collections, numpy as np, pandas as pd
from pathlib import Path

BASE = Path("/home/jg2070/Desktop/dtol_review_August")
SAT  = BASE/"DToL_phylogenomics_publication_325genomes/05_satellite_similarity"
ALL  = BASE/"2026_trees/all.satellites.txt"
CORE = "[CT]TTCGTTGGAA[AG]CGGGA"; FLANK = 5
PRE  = regex.compile("(?:TTCGTTGGAA){s<=3}")                 # core gate
PAT  = regex.compile("(?:%s){s<=5}" % CORE)                 # fixed 17-bp window
comp = str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
seqre= re.compile(r'"([ACGTNacgtn]{25,})"')
IDX  = {b:i for i,b in enumerate("ACGT")}

# species -> clade
tax=pd.read_csv(BASE/"2026_trees/annotation_centromeres/centromere_code_to_species.tsv",sep="\t")
tax["code"]=tax["fasta"].astype(str).str.lower().str.replace(r"[0-9.].*$","",regex=True)
name=dict(zip(tax.code,tax.genus+" "+tax.species))
vsub={"Mammalia":"Mammalia","Aves":"Aves","Actinopterygii":"Fish","Reptilia":"Reptilia",
      "Amphibia":"Amphibia","Chondrichthyes":"Fish"}
bmap={**{k:"Vertebrates" for k in vsub},"Fungi":"Fungi","Algae":"Viridiplantae","Bryophyta":"Viridiplantae",
      "Dicots":"Viridiplantae","Monocots":"Viridiplantae","Gymnosperms":"Viridiplantae",
      "Alveolata":"Protist","Discoba":"Protist"}
clade=dict(zip(tax.code,[bmap.get(t,"Invertebrate") for t in tax.taxa1]))
vgroup=dict(zip(tax.code,[vsub.get(t,"") for t in tax.taxa1]))

class Acc:
    __slots__=("counts","nwin","narr","bp")
    def __init__(self): self.counts=np.zeros((27,4)); self.nwin=0; self.narr=0; self.bp=0
acc=collections.defaultdict(Acc)

def add(sp,seq):
    a=acc[sp]; a.narr+=1; a.bp+=len(seq)
    both=seq+"#"+rc(seq)
    if not PRE.search(both): return
    for st in (seq,rc(seq)):
        for m in PAT.finditer(st,overlapped=False):
            s,e=m.start(),m.end()
            if e-s!=17 or s-FLANK<0 or e+FLANK>len(st): continue
            w=st[s-FLANK:e+FLANK]
            if len(w)!=27: continue
            ok=True
            for i,ch in enumerate(w):
                j=IDX.get(ch)
                if j is None: ok=False; break
                a.counts[i,j]+=1
            if ok: a.nwin+=1

n=0
for i,line in enumerate(open(ALL)):
    if i==0: continue
    m=seqre.search(line)
    if not m: continue
    sp=line.rstrip().rsplit('"',2)[-2].split(".")[0].lower() if '"' in line else None
    if not sp: continue
    add(sp,m.group(1).upper()); n+=1
    if n % 2000000==0: print(f"  {n:,} arrays...",flush=True)
print(f"total arrays: {n:,}",flush=True)

CANON=regex.compile("(?:%s){s<=17}"%CORE)
def ic(col_counts):                                          # bits at one position
    tot=col_counts.sum()
    if tot==0: return 0.0
    p=col_counts/tot; p=p[p>0]
    return 2.0-(-(p*np.log2(p)).sum())
rows=[]
for sp,a in acc.items():
    r=dict(species=sp, name=name.get(sp,sp), clade=clade.get(sp,"Invertebrate"),
           vgroup=vgroup.get(sp,""), n_arrays=a.narr, Mbp=round(a.bp/1e6,2), n_windows=a.nwin)
    r["win_per_Mbp"]=round(a.nwin/(a.bp/1e6),2) if a.bp else 0.0
    if a.nwin>=5:
        bits=np.array([ic(a.counts[i]) for i in range(27)])
        box=bits[FLANK:FLANK+17].mean(); fl=np.concatenate([bits[:FLANK],bits[FLANK+17:]]).mean()
        cons="".join("ACGT"[a.counts[i+FLANK].argmax()] for i in range(17))
        mm=CANON.fullmatch(cons)
        r.update(mean_box_bits=round(float(box),2), mean_flank_bits=round(float(fl),2),
                 delta=round(float(box-fl),2), box_consensus=cons,
                 subs_vs_canonical=(mm.fuzzy_counts[0] if mm else None))
    else:
        r.update(mean_box_bits=np.nan,mean_flank_bits=np.nan,delta=np.nan,box_consensus="",subs_vs_canonical=None)
    rows.append(r)
df=pd.DataFrame(rows).sort_values("delta",ascending=False,na_position="last")
df.to_csv(SAT/"figures/cenpb_flank_uncapped_per_species.tsv",sep="\t",index=False)
print("Saved:",SAT/"figures/cenpb_flank_uncapped_per_species.tsv")
print("\n=== top 25 by Delta (box - flank), all clades ===")
print(df[df.n_windows>=20][["name","species","clade","vgroup","n_windows","win_per_Mbp",
      "mean_box_bits","mean_flank_bits","delta","box_consensus","subs_vs_canonical"]].head(25).to_string(index=False))
print("\n=== VERTEBRATES (n_windows>=5), sorted by Delta ===")
v=df[(df.clade=="Vertebrates")&(df.n_windows>=5)]
print(v[["name","species","vgroup","n_arrays","n_windows","win_per_Mbp","mean_box_bits",
      "mean_flank_bits","delta","box_consensus","subs_vs_canonical"]].to_string(index=False))
