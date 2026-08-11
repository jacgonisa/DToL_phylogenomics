#!/usr/bin/env python3
"""CENP-B box tiers the way the functional-box paper does it (Barra/Fachinetti,
bioRxiv 2026.05.25.727640): three EXACT IUPAC motif definitions, both strands,
counted per species on the UNCAPPED satellite arrays, compared to a random
expectation from each species' own base composition (obs/expected enrichment).
No mismatch titration -- three discrete, well-defined tiers.

  canonical   YTTCGTTGGAARCGGGA   (functional box; degenerate only at pos 1 & 12)
  broad       NTTCGNNNNANNCGGGN
  degenerated YTTCGNNNNANRCGGGN

canonical (subset of) degenerated (subset of) broad, so counts are nested."""
import re, collections, numpy as np, pandas as pd
from pathlib import Path

BASE=Path("/home/jg2070/Desktop/dtol_review_August")
SAT =BASE/"DToL_phylogenomics_publication_325genomes/05_satellite_similarity"
ALL =BASE/"2026_trees/all.satellites.txt"
seqre=re.compile(r'"([ACGTNacgtn]{25,})"')
comp=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
MOT={"canonical":re.compile("[CT]TTCGTTGGAA[AG]CGGGA"),
     "broad":re.compile(".TTCG....A..CGGG."),
     "degenerated":re.compile("[CT]TTCG....A.[AG]CGGG.")}
# allowed base-sets per position, for the random expectation
ALLOWED={
 "canonical":[{"C","T"},{"T"},{"T"},{"C"},{"G"},{"T"},{"T"},{"G"},{"G"},{"A"},{"A"},{"A","G"},{"C"},{"G"},{"G"},{"G"},{"A"}],
 "broad":[set("ACGT"),{"T"},{"T"},{"C"},{"G"},set("ACGT"),set("ACGT"),set("ACGT"),set("ACGT"),{"A"},set("ACGT"),set("ACGT"),{"C"},{"G"},{"G"},{"G"},set("ACGT")],
 "degenerated":[{"C","T"},{"T"},{"T"},{"C"},{"G"},set("ACGT"),set("ACGT"),set("ACGT"),set("ACGT"),{"A"},set("ACGT"),{"A","G"},{"C"},{"G"},{"G"},{"G"},set("ACGT")]}

tax=pd.read_csv(BASE/"2026_trees/annotation_centromeres/centromere_code_to_species.tsv",sep="\t")
tax["code"]=tax["fasta"].astype(str).str.lower().str.replace(r"[0-9.].*$","",regex=True)
name=dict(zip(tax.code,tax.genus+" "+tax.species))
vsub={"Mammalia":"Mammalia","Aves":"Aves","Actinopterygii":"Fish","Reptilia":"Reptilia","Amphibia":"Amphibia","Chondrichthyes":"Fish"}
bmap={**{k:"Vertebrates" for k in vsub},"Fungi":"Fungi","Algae":"Viridiplantae","Bryophyta":"Viridiplantae",
      "Dicots":"Viridiplantae","Monocots":"Viridiplantae","Gymnosperms":"Viridiplantae","Alveolata":"Protist","Discoba":"Protist"}
clade=dict(zip(tax.code,[bmap.get(t,"Invertebrate") for t in tax.taxa1]))
vgroup=dict(zip(tax.code,[vsub.get(t,"") for t in tax.taxa1]))

class A:
    __slots__=("c","bp","base","narr")
    def __init__(self): self.c=collections.Counter(); self.bp=0; self.narr=0; self.base=collections.Counter()
acc=collections.defaultdict(A)
BASES="ACGT"
n=0
for i,line in enumerate(open(ALL)):
    if i==0: continue
    m=seqre.search(line)
    if not m: continue
    if '"' not in line: continue
    sp=line.rstrip().rsplit('"',2)[-2].split(".")[0].lower()
    s=m.group(1).upper(); a=acc[sp]; a.narr+=1; a.bp+=len(s)
    for b in BASES: a.base[b]+=s.count(b)
    both=s+"#"+rc(s)
    for k,p in MOT.items(): a.c[k]+=len(p.findall(both))
    n+=1
    if n%4000000==0: print(f"  {n:,}...",flush=True)
print(f"total arrays: {n:,}",flush=True)

def expected(a,motif):
    tot=sum(a.base[b] for b in BASES) or 1
    f={b:a.base[b]/tot for b in BASES}
    prob=1.0
    for st in ALLOWED[motif]: prob*= sum(f[b] for b in st)
    pos=2*max(a.bp-16*a.narr,0)                 # scannable 17-mer windows, both strands
    return pos*prob

rows=[]
for sp,a in acc.items():
    r=dict(species=sp,name=name.get(sp,sp),clade=clade.get(sp,"Invertebrate"),vgroup=vgroup.get(sp,""),
           n_arrays=a.narr,Mbp=round(a.bp/1e6,2))
    for k in MOT:
        exp=expected(a,k)
        r[f"{k}"]=int(a.c[k]); r[f"{k}_perMbp"]=round(a.c[k]/(a.bp/1e6),3) if a.bp else 0.0
        r[f"{k}_enrich"]=round(a.c[k]/exp,2) if exp>0 else np.nan
    rows.append(r)
df=pd.DataFrame(rows).sort_values(["canonical","degenerated","broad"],ascending=False)
df.to_csv(SAT/"figures/cenpb_paper_motifs_per_species.tsv",sep="\t",index=False)
print("Saved:",SAT/"figures/cenpb_paper_motifs_per_species.tsv")

# per-clade summary (observed per Mbp + enrichment over random expectation)
def cl(g):
    d=df[df.clade==g]
    return dict(clade=g,n_species=len(d),
        **{f"{k}_perMbp":round(d[f"{k}_perMbp"].mean(),3) for k in MOT},
        **{f"{k}_enrich":round(d[f"{k}_enrich"].mean(),2) for k in MOT},
        n_with_canonical=int((d["canonical"]>0).sum()))
per=pd.DataFrame([cl(g) for g in ["Vertebrates","Invertebrate","Viridiplantae","Fungi","Protist"] if (df.clade==g).any()])
per.to_csv(SAT/"figures/cenpb_paper_motifs_per_clade.tsv",sep="\t",index=False)
print("\n=== per-clade (mean per Mbp | mean enrichment over random | #species with a canonical box) ===")
print(per.to_string(index=False))
print("\n=== species with >=1 EXACT canonical box (functional) ===")
print(df[df.canonical>0][["name","species","clade","vgroup","n_arrays","canonical","degenerated","broad","canonical_enrich"]].to_string(index=False) or "  none")
print("\n=== VERTEBRATES: canonical/degenerated/broad per Mbp + enrichment ===")
v=df[df.clade=="Vertebrates"].sort_values("degenerated_perMbp",ascending=False)
print(v[["name","vgroup","n_arrays","canonical","degenerated","broad","degenerated_perMbp","broad_perMbp","degenerated_enrich","broad_enrich"]].to_string(index=False))
