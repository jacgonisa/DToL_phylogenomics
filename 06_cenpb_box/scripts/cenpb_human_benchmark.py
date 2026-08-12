#!/usr/bin/env python3
"""Human HG002 benchmark for BOTH CENP-B box methods, with a positive and a
negative control:
  POSITIVE  alpha-satellite  = cenpb_psi/human_alpha.fasta (50,000 monomers x171 bp)
  NEGATIVE  HSat (HSat1/2/3) = extracted from hg002v1.1.fasta via the v1.1 CenSat
            annotation (human satellites that do NOT carry CENP-B boxes)
Method 1 (Fachinetti, exact IUPAC canonical/broad/degenerate; obs/expected).
Method 2 (songbird +/-5 flank; box vs flank information Delta)."""
import re, numpy as np, pandas as pd, pysam, regex
from pathlib import Path

SAT=Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity")
ASM="/mnt/ssd-8tb/HUMAN/data/assembly/genomes/hg002v1.1.fasta"
HSATBED="/tmp/hg002_hsat.bed"
ALPHA=SAT/"cenpb_psi/human_alpha.fasta"
comp=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
# Method 1 motifs
MOT={"canonical":re.compile("[CT]TTCGTTGGAA[AG]CGGGA"),"broad":re.compile(".TTCG....A..CGGG."),
     "degenerated":re.compile("[CT]TTCG....A.[AG]CGGG.")}
ALLOWED={"canonical":[{"C","T"},{"T"},{"T"},{"C"},{"G"},{"T"},{"T"},{"G"},{"G"},{"A"},{"A"},{"A","G"},{"C"},{"G"},{"G"},{"G"},{"A"}],
 "broad":[set("ACGT"),{"T"},{"T"},{"C"},{"G"},set("ACGT"),set("ACGT"),set("ACGT"),set("ACGT"),{"A"},set("ACGT"),set("ACGT"),{"C"},{"G"},{"G"},{"G"},set("ACGT")],
 "degenerated":[{"C","T"},{"T"},{"T"},{"C"},{"G"},set("ACGT"),set("ACGT"),set("ACGT"),set("ACGT"),{"A"},set("ACGT"),{"A","G"},{"C"},{"G"},{"G"},{"G"},set("ACGT")]}
# Method 2
CORE="[CT]TTCGTTGGAA[AG]CGGGA"; FLANK=5
PRE=regex.compile("(?:TTCGTTGGAA){s<=3}"); PAT=regex.compile("(?:%s){s<=5}"%CORE); IDX={b:i for i,b in enumerate("ACGT")}

def load_alpha():
    out=[]; buf=[]
    for ln in open(ALPHA):
        if ln[0]==">":
            if buf: out.append("".join(buf)); buf=[]
        else: buf.append(ln.strip().upper())
    if buf: out.append("".join(buf)); return out

def load_hsat():
    fa=pysam.FastaFile(ASM); out=[]
    for ln in open(HSATBED):
        c,a,b,lab=ln.rstrip().split("\t")[:4]
        s=fa.fetch(c,int(a),int(b)).upper()
        if s: out.append(s)
    return out

def method1(seqs):
    base={b:0 for b in "ACGT"}; bp=0; narr=0; c={k:0 for k in MOT}
    for s in seqs:
        narr+=1; bp+=len(s)
        for b in "ACGT": base[b]+=s.count(b)
        both=s+"#"+rc(s)
        for k,p in MOT.items(): c[k]+=len(p.findall(both))
    tot=sum(base.values()) or 1; f={b:base[b]/tot for b in "ACGT"}
    r={"n_arrays":narr,"Mbp":round(bp/1e6,2)}
    for k in MOT:
        prob=1.0
        for st in ALLOWED[k]: prob*=sum(f[b] for b in st)
        exp=2*max(bp-16*narr,0)*prob
        r[f"{k}"]=c[k]; r[f"{k}_perMbp"]=round(c[k]/(bp/1e6),3) if bp else 0.0
        r[f"{k}_enrich"]=round(c[k]/exp,2) if exp>0 else np.nan
    return r

def method2(seqs):
    counts=np.zeros((27,4)); nwin=0; bp=0
    for s in seqs:
        bp+=len(s)
        both=s+"#"+rc(s)
        if not PRE.search(both): continue
        for st in (s,rc(s)):
            for m in PAT.finditer(st,overlapped=False):
                a,b=m.start(),m.end()
                if b-a!=17 or a-FLANK<0 or b+FLANK>len(st): continue
                w=st[a-FLANK:b+FLANK]
                if len(w)==27 and set(w)<=set("ACGT"):
                    for i,ch in enumerate(w): counts[i,IDX[ch]]+=1
                    nwin+=1
    def ic(col):
        t=col.sum()
        if t==0: return 0.0
        p=col/t; p=p[p>0]; return 2.0-(-(p*np.log2(p)).sum())
    r={"n_windows":nwin,"win_per_Mbp":round(nwin/(bp/1e6),2) if bp else 0.0}
    if nwin>=5:
        bits=np.array([ic(counts[i]) for i in range(27)])
        box=bits[FLANK:FLANK+17].mean(); fl=np.concatenate([bits[:FLANK],bits[FLANK+17:]]).mean()
        r["mean_box_bits"]=round(float(box),2); r["mean_flank_bits"]=round(float(fl),2); r["delta"]=round(float(box-fl),2)
    else:
        r["mean_box_bits"]=r["mean_flank_bits"]=r["delta"]=np.nan
    return r

rows=[]
for label,seqs in [("alpha-satellite (positive)",load_alpha()),("HSat1/2/3 (negative)",load_hsat())]:
    m1=method1(seqs); m2=method2(seqs)
    print(f"\n### {label}: {m1['n_arrays']:,} arrays, {m1['Mbp']} Mbp")
    print("  Method 1 (exact IUPAC):",{k:m1[k] for k in ("canonical","canonical_enrich","broad","broad_enrich","degenerated","degenerated_enrich")})
    print("  Method 2 (flank Delta):",{k:m2[k] for k in ("n_windows","win_per_Mbp","mean_box_bits","mean_flank_bits","delta")})
    rows.append({"control":label,**m1,**m2})
df=pd.DataFrame(rows)
df.to_csv(SAT/"figures/cenpb_human_benchmark.tsv",sep="\t",index=False)
print("\nSaved:",SAT/"figures/cenpb_human_benchmark.tsv")
