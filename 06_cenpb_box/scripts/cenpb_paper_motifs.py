#!/usr/bin/env python3
"""CENP-B box tiers the way the functional-box paper does it (Barra/Fachinetti,
bioRxiv 2026.05.25.727640): three EXACT IUPAC motif definitions, both strands,
counted per species on the UNCAPPED satellite arrays, compared to a random
expectation from each species' own base composition (obs/expected enrichment).
No mismatch titration -- three discrete, well-defined tiers.

  canonical   YTTCGTTGGAARCGGGA   (functional box; degenerate only at pos 1 & 12)
  broad       YTTCGNNNNANRCGGGN   (Sugimoto 1998; intermediate)
  degenerated NTTCGNNNNANNCGGGN   (most permissive; N at pos 1 & 12)

canonical (subset of) broad (subset of) degenerated, so counts are nested."""
import re, collections, numpy as np, pandas as pd
from pathlib import Path

BASE=Path("/home/jg2070/Desktop/dtol_review_August")
SAT =BASE/"DToL_phylogenomics_publication_325genomes/05_satellite_similarity"
ALL =BASE/"2026_trees/all.satellites.txt"
seqre=re.compile(r'"([ACGTNacgtn]{25,})"')
comp=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
MOT={"canonical":re.compile("[CT]TTCGTTGGAA[AG]CGGGA"),
     "broad":re.compile("[CT]TTCG....A.[AG]CGGG."),
     "degenerated":re.compile(".TTCG....A..CGGG.")}
# allowed base-sets per position, for the random expectation
ALLOWED={
 "canonical":[{"C","T"},{"T"},{"T"},{"C"},{"G"},{"T"},{"T"},{"G"},{"G"},{"A"},{"A"},{"A","G"},{"C"},{"G"},{"G"},{"G"},{"A"}],
 "broad":[{"C","T"},{"T"},{"T"},{"C"},{"G"},set("ACGT"),set("ACGT"),set("ACGT"),set("ACGT"),{"A"},set("ACGT"),{"A","G"},{"C"},{"G"},{"G"},{"G"},set("ACGT")],
 "degenerated":[set("ACGT"),{"T"},{"T"},{"C"},{"G"},set("ACGT"),set("ACGT"),set("ACGT"),set("ACGT"),{"A"},set("ACGT"),set("ACGT"),{"C"},{"G"},{"G"},{"G"},set("ACGT")]}

tax=pd.read_csv(BASE/"2026_trees/annotation_centromeres/centromere_code_to_species.tsv",sep="\t")
tax["code"]=tax["fasta"].astype(str).str.lower().str.replace(r"[0-9.].*$","",regex=True)
name=dict(zip(tax.code,tax.genus+" "+tax.species))
vsub={"Mammalia":"Mammalia","Aves":"Aves","Actinopterygii":"Fish","Reptilia":"Reptilia","Amphibia":"Amphibia","Chondrichthyes":"Fish"}
bmap={**{k:"Vertebrates" for k in vsub},"Fungi":"Fungi","Algae":"Viridiplantae","Bryophyta":"Viridiplantae",
      "Dicots":"Viridiplantae","Monocots":"Viridiplantae","Gymnosperms":"Viridiplantae","Alveolata":"Protist","Discoba":"Protist"}
clade=dict(zip(tax.code,[bmap.get(t,"Invertebrate") for t in tax.taxa1]))
vgroup=dict(zip(tax.code,[vsub.get(t,"") for t in tax.taxa1]))
# restrict to the 325 published species (tree tips); all.satellites.txt carries extras
ALLOW=set(l.strip() for l in open(SAT/"species_325.txt") if l.strip() and not l.startswith("#"))

# exact overlapping dinucleotide counting (numpy, C-level)
LUT=np.full(256,-1,dtype=np.int8)
for i,b in enumerate(b"ACGT"): LUT[b]=i
def di_counts(s):
    arr=LUT[np.frombuffer(s.encode("ascii","ignore"),dtype=np.uint8)]
    arr=arr[arr>=0]
    if arr.size<2: return np.zeros(16,dtype=np.int64)
    return np.bincount((arr[:-1].astype(np.int64)<<2)|arr[1:],minlength=16)

class A:
    __slots__=("c","bp","base","narr","di")
    def __init__(self): self.c=collections.Counter(); self.bp=0; self.narr=0; self.base=collections.Counter(); self.di=np.zeros(16,dtype=np.int64)
acc=collections.defaultdict(A)
BASES="ACGT"
n=0
for i,line in enumerate(open(ALL)):
    if i==0: continue
    m=seqre.search(line)
    if not m: continue
    if '"' not in line: continue
    sp=line.rstrip().rsplit('"',2)[-2].split(".")[0].lower()
    if sp not in ALLOW: continue
    s=m.group(1).upper(); a=acc[sp]; a.narr+=1; a.bp+=len(s)
    for b in BASES: a.base[b]+=s.count(b)
    rs=rc(s); a.di+=di_counts(s)+di_counts(rs)          # both strands
    both=s+"#"+rs
    for k,p in MOT.items(): a.c[k]+=len(p.findall(both))
    n+=1
    if n%4000000==0: print(f"  {n:,}...",flush=True)
print(f"total arrays: {n:,}",flush=True)

IB={b:i for i,b in enumerate("ACGT")}
def expected(a,motif):                                  # 0-order (mononucleotide) null
    tot=sum(a.base[b] for b in BASES) or 1
    f={b:a.base[b]/tot for b in BASES}
    prob=1.0
    for st in ALLOWED[motif]: prob*= sum(f[b] for b in st)
    return 2*max(a.bp-16*a.narr,0)*prob

def expected_dinuc(a,motif):                            # 1st-order Markov (dinucleotide-preserving) null
    D=a.di.reshape(4,4).astype(float)
    rowsum=D.sum(1)
    if (rowsum<=0).any(): return float("nan")
    T=D/rowsum[:,None]                                  # transition P(next=b | cur=a)
    tot=sum(a.base[b] for b in BASES) or 1
    pi=np.array([a.base[b]/tot for b in BASES])         # start distribution
    allowed=[ {IB[x] for x in st} for st in ALLOWED[motif] ]
    fvec=np.array([pi[b] if b in allowed[0] else 0.0 for b in range(4)])
    for pos in range(1,17):
        nf=fvec@T
        fvec=np.array([nf[b] if b in allowed[pos] else 0.0 for b in range(4)])
    return 2*max(a.bp-16*a.narr,0)*fvec.sum()

rows=[]
for sp,a in acc.items():
    r=dict(species=sp,name=name.get(sp,sp),clade=clade.get(sp,"Invertebrate"),vgroup=vgroup.get(sp,""),
           n_arrays=a.narr,Mbp=round(a.bp/1e6,2))
    for k in MOT:
        exp=expected(a,k); expd=expected_dinuc(a,k)
        r[f"{k}"]=int(a.c[k]); r[f"{k}_perMbp"]=round(a.c[k]/(a.bp/1e6),3) if a.bp else 0.0
        r[f"{k}_enrich"]=round(a.c[k]/exp,2) if exp>0 else np.nan            # mononucleotide null
        r[f"{k}_enrich_dinuc"]=round(a.c[k]/expd,2) if expd and expd>0 else np.nan  # dinucleotide (Markov-1) null
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
        **{f"{k}_enrich_dinuc":round(d[f"{k}_enrich_dinuc"].mean(),2) for k in MOT},
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
