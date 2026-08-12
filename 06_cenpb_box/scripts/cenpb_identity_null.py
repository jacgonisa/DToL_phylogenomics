#!/usr/bin/env python3
"""Null for the motif-identity metrics (per WINDOW, not consensus — the consensus
version is degenerate: selecting near-canonical windows regenerates canonical).
For each species we collect its box-like windows (17-mers within ≤5 substitutions
of canonical, both strands) and compute, per window:
  identity_canonical = mean identity of each window to the canonical CENP-B motif
  identity_motif     = mean identity of each window to the species' OWN identified
                       motif consensus (how coherent/conserved the candidate box is)
Null = the identical pipeline on a DINUCLEOTIDE-PRESERVING (Altschul–Erikson)
shuffle of the same arrays. excess = observed − null isolates real signal.
Input = candidate satellite arrays; subsampled to CAP arrays/species."""
import regex, re, collections, numpy as np, pandas as pd, random
from pathlib import Path
from ae_shuffle import dinucl_shuffle
random.seed(0)
BASE=Path("/home/jg2070/Desktop/dtol_review_August")
SAT =BASE/"DToL_phylogenomics_publication_325genomes/05_satellite_similarity"
ALL =BASE/"2026_trees/all.satellites.txt"
CANON="[CT]TTCGTTGGAA[AG]CGGGA"
PRE=regex.compile("(?:TTCGTTGGAA){s<=3}"); PAT=regex.compile("(?:%s){s<=5}"%CANON)
comp=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
seqre=re.compile(r'"([ACGTNacgtn]{25,})"'); IDX={b:i for i,b in enumerate("ACGT")}
ALLOW=set(l.strip() for l in open(SAT/"species_325.txt") if l.strip() and not l.startswith("#"))
CAP=3000

def new(): return dict(pfm=np.zeros((17,4)), N=0, subs=0)
def add(acc, seq):
    for st in (seq, rc(seq)):
        for m in PAT.finditer(st, overlapped=False):
            a,b=m.start(),m.end()
            if b-a!=17: continue
            w=st[a:b]
            if set(w)-set("ACGT"): continue
            acc["subs"]+=m.fuzzy_counts[0]; acc["N"]+=1        # subs to canonical (IUPAC)
            for i,ch in enumerate(w): acc["pfm"][i,IDX[ch]]+=1
def ids(acc):
    N=acc["N"]
    if N<10: return (np.nan,np.nan)
    id_canon=100*(17 - acc["subs"]/N)/17                       # mean per-window identity to canonical
    id_motif=100*acc["pfm"].max(1).sum()/(17*N)                # mean per-window identity to own consensus
    return (round(id_canon,1), round(id_motif,1))

obs=collections.defaultdict(new); nul=collections.defaultdict(new); narr=collections.Counter()
n=0
for i,line in enumerate(open(ALL)):
    if i==0 or '"' not in line: continue
    sp=line.rstrip().rsplit('"',2)[-2].split(".")[0].lower()
    if sp not in ALLOW or narr[sp]>=CAP: continue
    m=seqre.search(line)
    if not m: continue
    s=m.group(1).upper(); narr[sp]+=1; n+=1
    if PRE.search(s+"#"+rc(s)): add(obs[sp], s)
    sh=dinucl_shuffle(s)
    if PRE.search(sh+"#"+rc(sh)): add(nul[sp], sh)
    if n%2000000==0: print(f"  {n:,} arrays...",flush=True)
print(f"streamed {n:,} capped arrays",flush=True)

rows=[]
for sp in sorted(set(list(obs)+list(nul))):
    oc,om=ids(obs[sp]); nc,nm=ids(nul[sp])
    rows.append(dict(species=sp, n_arrays=narr[sp], n_win=obs[sp]["N"], n_win_null=nul[sp]["N"],
        id_canon_obs=oc, id_canon_null=nc, excess_canon=(round(oc-nc,1) if oc==oc and nc==nc else np.nan),
        id_motif_obs=om, id_motif_null=nm, excess_motif=(round(om-nm,1) if om==om and nm==nm else np.nan)))
df=pd.DataFrame(rows); df.to_csv(SAT/"figures/cenpb_identity_null.tsv",sep="\t",index=False)
print("Saved figures/cenpb_identity_null.tsv")
d=df.dropna(subset=["excess_canon"])
print("\n-- identity to CANONICAL (per window) --  obs %.1f | null %.1f | excess %.1f  (>0: %d/%d)"%(
    d.id_canon_obs.median(), d.id_canon_null.median(), d.excess_canon.median(), (d.excess_canon>0).sum(), len(d)))
print("-- identity to OWN MOTIF (coherence)  --  obs %.1f | null %.1f | excess %.1f  (>0: %d/%d)"%(
    d.id_motif_obs.median(), d.id_motif_null.median(), d.excess_motif.median(), (d.excess_motif>0).sum(), len(d)))
print("\ntop excess (own-motif coherence over null):")
print(d.sort_values("excess_motif",ascending=False)[["species","n_win","id_canon_obs","excess_canon","id_motif_obs","excess_motif"]].head(12).to_string(index=False))
