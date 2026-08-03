#!/usr/bin/env python3
"""Look harder for CENP-B boxes in each mammal: canonical motif allowing 0-3 SUBSTITUTIONS
(no indels), both strands, on all satellite records, with a shuffled-sequence null.
Enrichment = %real / %null. A real (diverged) box => real >> null."""
import regex, random, collections
random.seed(0)
FA="cenpb_psi/all_mammals.fasta"
CANON="[CT]TTCGTTGGAA[AG]CGGGA"
P={k:regex.compile("(?:%s){s<=%d}"%(CANON,k)) for k in (0,1,2,3)}
comp=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
def hit(p,s): return bool(p.search(s) or p.search(rc(s)))

real=collections.defaultdict(lambda:collections.Counter()); nul=collections.defaultdict(lambda:collections.Counter())
n=collections.Counter()
nm=None
for line in open(FA):
    if line[0]==">": nm=line[1:].strip().split("__")[0]; continue
    s=line.strip().upper()
    if len(s)<17: continue
    n[nm]+=1
    sh="".join(random.sample(s,len(s)))
    for k in (0,1,2,3):
        if hit(P[k],s): real[nm][k]+=1
        if hit(P[k],sh): nul[nm][k]+=1

print(f"{'species':9s} {'n':>9} | {'k=0 real/null':>14} {'k=1':>12} {'k=2':>14} {'k=3 (enrich)':>18}")
for sp in sorted(n,key=lambda x:-real[x][3]/max(1,n[x])):
    N=n[sp]
    def cell(k):
        r=100*real[sp][k]/N; u=100*nul[sp][k]/N
        e=r/u if u>0 else (float('inf') if r>0 else 0)
        return f"{r:.2f}/{u:.2f}({e:.1f}x)" if k==3 else f"{r:.2f}/{u:.2f}"
    print(f"{sp:9s} {N:>9,} | {cell(0):>14} {cell(1):>12} {cell(2):>14} {cell(3):>18}")
