#!/usr/bin/env python3
"""Uncapped CENP-B box scan across ALL DToL satellite arrays (all.satellites.txt,
~24.5M arrays; species = last column {Code}.{family}). Two methods:
  (1) exact/mismatch match to the canonical 17-bp box (k=0/1/2; k=0 is null-free)
  (2) PWM log-odds scan (PWM built from the human alpha-sat box) -> catches
      diverged boxes that strict mismatch misses.
Human alpha-satellite (HG002) scanned as the positive benchmark."""
import regex, re, math, collections
from pathlib import Path
import pandas as pd

BASE = Path("/home/jg2070/Desktop/dtol_review_August")
SAT  = BASE/"DToL_phylogenomics_publication_325genomes/05_satellite_similarity"
ALL  = BASE/"2026_trees/all.satellites.txt"
HUMAN= SAT/"cenpb_psi/human_alpha.fasta"
CORE17="[CT]TTCGTTGGAA[AG]CGGGA"
PATS ={k:regex.compile("(?:%s){e<=%d}"%(CORE17,k)) for k in range(3)}
PRE  = regex.compile("(?:TTCGTTGGAA){e<=1}")           # strict pre-filter
FIND = regex.compile("(?:%s){e<=4}"%CORE17)            # window finder for PWM (loose)
comp = str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
seqre= re.compile(r'"([ACGTNacgtn]{25,})"')

# ---- build PWM from human alpha-sat box instances ----
def load_fa(p):
    out=[]; buf=[]
    for ln in open(p):
        if ln[0]==">":
            if buf: out.append("".join(buf)); buf=[]
        else: buf.append(ln.strip().upper())
    if buf: out.append("".join(buf)); return out
hum = load_fa(HUMAN)
boxes=[]
for s in hum:
    for strand in (s, rc(s)):
        m=PATS[2].search(strand)                       # <=2mm human box instance
        if m: boxes.append(strand[m.start():m.start()+17])
counts=[collections.Counter() for _ in range(17)]
for b in boxes:
    if len(b)==17:
        for i,ch in enumerate(b):
            if ch in "ACGT": counts[i][ch]+=1
bg={"A":.3,"C":.2,"G":.2,"T":.3}
PWM=[]
for i in range(17):
    tot=sum(counts[i].values()) or 1
    PWM.append({b: math.log2(((counts[i][b]+0.5)/(tot+2))/bg[b]) for b in "ACGT"})
maxscore=sum(max(PWM[i].values()) for i in range(17))
THR=0.75*maxscore                                       # PWM hit threshold (75% of max)
print(f"human box instances for PWM: {len(boxes)} | maxscore {maxscore:.1f} | THR {THR:.1f}", flush=True)
def pwm_hits(seq):
    n=0; best=-1e9
    for strand in (seq, rc(seq)):
        for m in FIND.finditer(strand):
            w=strand[m.start():m.start()+17]
            if len(w)!=17: continue
            sc=sum(PWM[i].get(w[i],-2) for i in range(17))
            best=max(best,sc)
            if sc>=THR: n+=1
    return n, best

# ---- species -> clade ----
tax=pd.read_csv(BASE/"2026_trees/annotation_centromeres/centromere_code_to_species.tsv",sep="\t")
tax["code"]=tax["fasta"].astype(str).str.lower().str.replace(r"[0-9.].*$","",regex=True)
bmap={"Mammalia":"Vertebrates","Aves":"Vertebrates","Actinopterygii":"Vertebrates","Reptilia":"Vertebrates",
      "Amphibia":"Vertebrates","Chondrichthyes":"Vertebrates","Fungi":"Fungi","Algae":"Viridiplantae",
      "Bryophyta":"Viridiplantae","Dicots":"Viridiplantae","Monocots":"Viridiplantae","Gymnosperms":"Viridiplantae",
      "Alveolata":"Protist","Discoba":"Protist"}
clade={c:bmap.get(t,"Invertebrate") for c,t in zip(tax.code,tax.taxa1)}

def blank(): return dict(n_arrays=0, bp=0, k0=0, k1=0, k2=0, pwm=0, best_pwm=-1e9)
acc=collections.defaultdict(blank)

def add(sp, seq):
    a=acc[sp]; a["n_arrays"]+=1; a["bp"]+=len(seq)
    both=seq+"#"+rc(seq)
    if PRE.search(both):
        a["k0"]+=len(PATS[0].findall(both)); a["k1"]+=len(PATS[1].findall(both)); a["k2"]+=len(PATS[2].findall(both))
        h,b=pwm_hits(seq); a["pwm"]+=h; a["best_pwm"]=max(a["best_pwm"],b)

# ---- human benchmark ----
for s in hum: add("human", s)
acc["human"]  # ensure present
# ---- stream all DToL satellites ----
n=0
for i,line in enumerate(open(ALL)):
    if i==0: continue
    m=seqre.search(line)
    if not m: continue
    sp=line.rstrip().rsplit('"',2)[-2].split(".")[0].lower() if '"' in line else None
    if not sp: continue
    add(sp, m.group(1).upper()); n+=1
    if n % 2000000==0: print(f"  {n:,} arrays...", flush=True)
print(f"total DToL arrays scanned: {n:,}", flush=True)

rows=[]
for sp,a in acc.items():
    rows.append(dict(species=sp, clade="Human" if sp=="human" else clade.get(sp,"Invertebrate"),
        n_arrays=a["n_arrays"], Mbp=round(a["bp"]/1e6,2),
        exact_k0=a["k0"], box_k1=a["k1"], box_k2=a["k2"],
        pwm_hits=a["pwm"], best_pwm=round(a["best_pwm"],1)))
df=pd.DataFrame(rows).sort_values(["exact_k0","pwm_hits"],ascending=False)
df.to_csv(SAT/"figures/cenpb_uncapped_per_species.tsv",sep="\t",index=False)
print("\n=== species with exact (k=0) canonical boxes ===")
print(df[df.exact_k0>0][["species","clade","n_arrays","exact_k0","box_k2","pwm_hits"]].to_string(index=False))
print("\n=== top 12 by PWM hits ===")
print(df.sort_values("pwm_hits",ascending=False)[["species","clade","exact_k0","box_k2","pwm_hits","best_pwm"]].head(12).to_string(index=False))
