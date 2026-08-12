#!/usr/bin/env python3
"""Altschul-Erikson dinucleotide-preserving shuffle (Altschul & Erikson 1985).
Produces a random sequence with EXACTLY the same dinucleotide (doublet) counts
and same first/last base as the input -- the proper null for motifs whose base
context (e.g. CpG in TTCG/CGGG) matters. Clote-style Eulerian-path implementation."""
import random

def _counts_lists(s):
    nucl={n:0 for n in "ACGT"}; di={a:{b:0 for b in "ACGT"} for a in "ACGT"}
    nucl[s[0]]=1
    for i in range(len(s)-1):
        x,y=s[i],s[i+1]; di[x][y]+=1; nucl[y]+=1
    return nucl,di

def _choose_edge(x,di):
    tot=di[x]['A']+di[x]['C']+di[x]['G']+di[x]['T']
    z=int(random.random()*tot); run=0
    for y in "ACGT":
        run+=di[x][y]
        if z<run: di[x][y]-=1; return y
    return 'T'

def _connected(last,nuclList,lastCh):
    D={x:0 for x in nuclList}
    for a,b in last:
        if b==lastCh: D[a]=1
    for _ in range(len(nuclList)):
        for a,b in last:
            if D.get(b,0)==1: D[a]=1
    return all(D[x]==1 for x in nuclList if x!=lastCh)

def dinucl_shuffle(s):
    s=s.upper()
    if len(s)<3 or set(s)-set("ACGT"): return s   # skip Ns / very short
    nucl,di=_counts_lists(s)
    nuclList=[x for x in "ACGT" if nucl[x]>0]
    firstCh,lastCh=s[0],s[-1]
    if firstCh==lastCh or len(nuclList)==1: pass
    while True:                                    # pick one "last edge" per vertex -> arborescence to lastCh
        d={a:{b:di[a][b] for b in "ACGT"} for a in "ACGT"}
        last=[(x,_choose_edge(x,d)) for x in nuclList if x!=lastCh]
        if _connected(last,nuclList,lastCh): break
    lastMap=dict(last)
    L={}
    for x in nuclList:
        lst=[]
        for y in "ACGT": lst+=[y]*d[x][y]
        random.shuffle(lst)
        if x in lastMap: lst.append(lastMap[x])
        L[x]=lst
    out=[firstCh]; cur=firstCh; ptr={x:0 for x in nuclList}
    for _ in range(len(s)-1):
        nxt=L[cur][ptr[cur]]; ptr[cur]+=1; out.append(nxt); cur=nxt
    return "".join(out)

if __name__=="__main__":                           # self-check: doublet counts + ends preserved
    random.seed(1)
    def di(s):
        from collections import Counter; return Counter(s[i:i+2] for i in range(len(s)-1))
    for t in ["ACGTACGTGGCCATATCGCGAAATTTCGGGACGT","TTTTAAAACCCCGGGG","ACGTTGCAACGT"*5]:
        sh=dinucl_shuffle(t)
        assert len(sh)==len(t) and sh[0]==t[0] and sh[-1]==t[-1], "ends/length"
        assert di(sh)==di(t), f"dinucleotide counts differ:\n{di(t)}\n{di(sh)}"
    print("ae_shuffle self-check OK: dinucleotide counts + first/last base preserved")
