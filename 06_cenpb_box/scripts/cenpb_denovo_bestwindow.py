#!/usr/bin/env python3
"""De-novo, UNBIASED test for a box-like motif in each species' satellites.
Unlike the ±5-flank search (which seeds on the CENP-B box and so only ever collects
windows ~5 substitutions away), here we scan EVERY 17-bp window of every array
(both strands) and take the single best match to the canonical CENP-B IUPAC motif
(no seeding). The null is the same on a DINUCLEOTIDE-shuffle of the array — which
controls for the 'best of many windows' inflation. A satellite that genuinely
contains a box-like region has obs best-window identity >> shuffle; an unrelated
satellite has obs ≈ shuffle. Human α-sat (+) and HSat (−) are the controls."""
import re, random, collections, numpy as np, pandas as pd, pysam
from pathlib import Path
from ae_shuffle import dinucl_shuffle
random.seed(0)
BASE=Path("/home/jg2070/Desktop/dtol_review_August")
SAT =BASE/"DToL_phylogenomics_publication_325genomes/05_satellite_similarity"
ALL =BASE/"2026_trees/all.satellites.txt"
seqre=re.compile(r'"([ACGTNacgtn]{25,})"')
LUT=np.full(256,-1,np.int8)
for i,b in enumerate(b"ACGT"): LUT[b]=i
IDX={"A":0,"C":1,"G":2,"T":3}
CANSETS=[set("CT"),{"T"},{"T"},{"C"},{"G"},{"T"},{"T"},{"G"},{"G"},{"A"},{"A"},set("AG"),{"C"},{"G"},{"G"},{"G"},{"A"}]
ALLOW=np.zeros((17,4),bool)
for j,st in enumerate(CANSETS):
    for b in st: ALLOW[j][IDX[b]]=1
CAP=1500; MAXLEN=600                                    # arrays/species and bp/array caps

def best_id(seq):                                      # max identity (/17*100) of any window, both strands
    a=LUT[np.frombuffer(seq.encode("ascii","ignore"),np.uint8)]
    if (a<0).any(): a=a[a>=0]
    best=0
    for strand in (a,(3-a)[::-1]):
        n=len(strand)-16
        if n<1: continue
        sc=np.zeros(n,np.int16)
        for j in range(17): sc+=ALLOW[j][strand[j:j+n]]
        m=int(sc.max())
        if m>best: best=m
    return 100.0*best/17

def species_stats(seqs):
    seqs=[s[:MAXLEN] for s in seqs if set(s)<=set("ACGT") and len(s)>=25][:CAP]
    if len(seqs)<10: return None
    obs=np.mean([best_id(s) for s in seqs])
    nul=np.mean([best_id(dinucl_shuffle(s)) for s in seqs])
    return round(obs,1),round(nul,1),round(obs-nul,1),len(seqs)

def load_fa(p):
    out=[]; buf=[]
    for ln in open(p):
        if ln[0]==">":
            if buf: out.append("".join(buf)); buf=[]
        else: buf.append(ln.strip().upper())
    if buf: out.append("".join(buf)); return out
def load_hsat():
    fa=pysam.FastaFile("/mnt/ssd-8tb/HUMAN/data/assembly/genomes/hg002v1.1.fasta"); out=[]
    for ln in open("/tmp/hg002_hsat.bed"):
        c,a,b,_=ln.rstrip().split("\t")[:4]; s=fa.fetch(c,int(a),int(b)).upper()
        if s: out.append(s)
    return out

# --- collect DToL satellites per species ---
ALLOWSP=set(l.strip() for l in open(SAT/"species_325.txt") if l.strip() and not l.startswith("#"))
seqs=collections.defaultdict(list); n=0
for i,line in enumerate(open(ALL)):
    if i==0 or '"' not in line: continue
    sp=line.rstrip().rsplit('"',2)[-2].split(".")[0].lower()
    if sp not in ALLOWSP or len(seqs[sp])>=CAP: continue
    m=seqre.search(line)
    if m: seqs[sp].append(m.group(1).upper()); n+=1
    if n%4000000==0: print(f"  {n:,} arrays read...",flush=True)
print(f"collected {len(seqs)} species",flush=True)

# clade map
tax=pd.read_csv(BASE/"2026_trees/annotation_centromeres/centromere_code_to_species.tsv",sep="\t")
tax["code"]=tax["fasta"].astype(str).str.lower().str.replace(r"[0-9.].*$","",regex=True)
name=dict(zip(tax.code,tax.genus+" "+tax.species))
vsub={"Mammalia":"Mammalia","Aves":"Aves","Actinopterygii":"Fish","Reptilia":"Reptilia","Amphibia":"Amphibia","Chondrichthyes":"Fish"}
bmap={**{k:"Vertebrates" for k in vsub},"Fungi":"Fungi","Algae":"Viridiplantae","Bryophyta":"Viridiplantae",
      "Dicots":"Viridiplantae","Monocots":"Viridiplantae","Gymnosperms":"Viridiplantae","Alveolata":"Protist","Discoba":"Protist"}
clade=dict(zip(tax.code,[bmap.get(t,"Invertebrate") for t in tax.taxa1]))
vgroup=dict(zip(tax.code,[vsub.get(t,"") for t in tax.taxa1]))

rows=[]
for sp,ss in seqs.items():
    st=species_stats(ss)
    if st: rows.append(dict(species=sp,name=name.get(sp,sp),clade=clade.get(sp,"Invertebrate"),
                            vgroup=vgroup.get(sp,""),best_obs=st[0],best_null=st[1],excess=st[2],n=st[3]))
# controls
for lab,seqsc,cl in [("alpha-satellite (positive)",load_fa(SAT/"cenpb_psi/human_alpha.fasta"),"Human+"),
                     ("HSat1/2/3 (negative)",load_hsat(),"Human-")]:
    st=species_stats(seqsc)
    rows.append(dict(species=lab,name=lab,clade=cl,vgroup="",best_obs=st[0],best_null=st[1],excess=st[2],n=st[3]))
df=pd.DataFrame(rows).sort_values("excess",ascending=False)
df.to_csv(SAT/"figures/cenpb_denovo_bestwindow.tsv",sep="\t",index=False)
print("Saved figures/cenpb_denovo_bestwindow.tsv")
print("\ncontrols + top/bottom species by excess (obs best-window identity − shuffle):")
print(df[df.clade.str.contains("Human")][["name","best_obs","best_null","excess"]].to_string(index=False))
print(df[~df.clade.str.contains("Human")].head(10)[["name","clade","best_obs","best_null","excess","n"]].to_string(index=False))
print("...\n", df[df.name.str.contains("Accipiter|Lagopus|Porphyrio|Corvus",na=False)][["name","best_obs","best_null","excess"]].to_string(index=False))
d=df[~df.clade.str.contains("Human")]
print(f"\nDToL: median excess {d.excess.median():.1f} ; species with excess>=+10: {(d.excess>=10).sum()}/{len(d)}")
