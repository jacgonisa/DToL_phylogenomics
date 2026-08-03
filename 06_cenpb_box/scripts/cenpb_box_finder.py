#!/usr/bin/env python3
"""
cenpb_box_finder.py — screen satellite (or any) sequences for the CENP-B box using
several independent lines of evidence, per species, each with a shuffled-sequence null:

  1. canonical  IUPAC exact match  (Masumoto 1989)      YTTCGTTGGAARCGGGA
  2. broad      IUPAC exact match  (marmoset, PMC4843215) NTTCGNNNNANNCGGGN
  3. degenerate IUPAC exact match  (Altemose 2022 style; --degenerate to set exactly)
  4. PWM        FIMO p-value scan  (data-driven, FDR from a shuffled DB)

Input FASTA headers must be >{species}__{n}. Outputs a per-species evidence table
(counts + null counts per line of evidence) and a stacked figure. A species is called
CENP-B-box-positive when >=2 lines exceed their null.

Usage: cenpb_box_finder.py --fasta iter_db.fasta [--pwm cenpb_final.meme]
       [--degenerate NTTCGNNNNNNNNCGGGN] [--thresh 1e-7] [--out figures]
"""
import argparse, re, subprocess, collections, random, os
from pathlib import Path
import pandas as pd
random.seed(0)
FIMO="/home/jg2070/meme/bin/fimo"
IUPAC={'A':'A','C':'C','G':'G','T':'T','R':'[AG]','Y':'[CT]','S':'[GC]','W':'[AT]',
       'K':'[GT]','M':'[AC]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACGT]'}
comp=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
def iupac_re(p): return re.compile("".join(IUPAC[c] for c in p))

def read_fasta(path):
    nm=None; buf=[]
    for l in open(path):
        if l[0]==">":
            if nm is not None: yield nm,"".join(buf)
            nm=l[1:].strip(); buf=[]
        else: buf.append(l.strip())
    if nm is not None: yield nm,"".join(buf)

def shuffle_fasta(src,dst):
    with open(dst,"w") as o:
        for nm,s in read_fasta(src):
            l=list(s.upper()); random.shuffle(l); o.write(f">{nm}\n{''.join(l)}\n")

def count_regex(path,motifs):
    """per species: n records, and n records with >=1 match (both strands) per motif."""
    n=collections.Counter(); hit={m:collections.Counter() for m in motifs}
    for nm,s in read_fasta(path):
        sp=nm.split("__")[0]; s=s.upper(); r=rc(s); n[sp]+=1
        for m,rx in motifs.items():
            if rx.search(s) or rx.search(r): hit[m][sp]+=1
    return n,hit

def fimo_hits(pwm,path,thresh,bg):
    r=subprocess.run([FIMO,"--text","--thresh",thresh,"--bfile",bg,pwm,path],
                     capture_output=True,text=True)
    per=collections.Counter()
    seen=set()
    for ln in r.stdout.splitlines():
        p=ln.split("\t")
        if len(p)<10 or p[0]=="motif_id": continue
        if p[2] not in seen: seen.add(p[2]); per[p[2].split("__")[0]]+=1
    return per

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--fasta",required=True)
    ap.add_argument("--pwm",default=None,help="MEME motif for the PWM/FIMO line")
    ap.add_argument("--bg",default=None,help="FIMO background file (order-0)")
    ap.add_argument("--canonical",default="YTTCGTTGGAARCGGGA")
    ap.add_argument("--broad",default="NTTCGNNNNANNCGGGN")
    ap.add_argument("--degenerate",default="NTTCGNNNNNNNCGGGN",help="set to the exact Altemose IUPAC")
    ap.add_argument("--thresh",default="1e-7")
    ap.add_argument("--taxa",default="/home/jg2070/Desktop/dtol_review_August/2026_trees/annotation_centromeres/DTOL_327_master_March.xlsx")
    ap.add_argument("--out",default="figures")
    a=ap.parse_args()
    out=Path(a.out); out.mkdir(exist_ok=True)
    motifs={k:iupac_re(v) for k,v in
            [("canonical",a.canonical),("broad",a.broad),("degenerate",a.degenerate)]}

    shuf=str(Path(a.fasta).with_suffix(".finder_shuf.fasta"))
    if not os.path.exists(shuf): shuffle_fasta(a.fasta,shuf)

    n,hit=count_regex(a.fasta,motifs)
    _,hnull=count_regex(shuf,motifs)
    cols={"n_records":n}
    for m in motifs: cols[m]=hit[m]; cols[m+"_null"]=hnull[m]
    if a.pwm:
        bg=a.bg or str(out/"finder_bg.txt")
        if not a.bg:
            b=collections.Counter()
            for _,s in read_fasta(a.fasta):
                for ch in s.upper(): b[ch]+=1
            tot=sum(b[x] for x in "ACGT")
            open(bg,"w").write("".join(f"{x} {b[x]/tot:.5f}\n" for x in "ACGT"))
        cols["pwm"]=fimo_hits(a.pwm,a.fasta,a.thresh,bg)
        cols["pwm_null"]=fimo_hits(a.pwm,shuf,a.thresh,bg)

    # taxa map
    xl=pd.read_excel(a.taxa); xl["code"]=xl["fasta"].astype(str).str.lower().str.replace(r"[0-9].*$","",regex=True).str.replace(r"[.].*$","",regex=True)
    t=dict(zip(xl["code"],xl["taxa1"]))
    sp=sorted(n)
    ev=[m for m in list(motifs)+(["pwm"] if a.pwm else [])]
    rows=[]
    for s in sp:
        d={"species":s,"taxa1":t.get(s,"NA"),"n":n[s]}
        pos=0
        for m in ev:
            d[m]=cols[m][s]; d[m+"_null"]=cols[m+"_null"][s]
            if cols[m][s]>cols[m+"_null"][s] and cols[m][s]>=max(3,0.005*n[s]): pos+=1
        d["lines_positive"]=pos; rows.append(d)
    df=pd.DataFrame(rows).sort_values(["lines_positive","pwm" if a.pwm else "broad"],ascending=False)
    df.to_csv(out/"cenpb_box_evidence_per_species.tsv",sep="\t",index=False)
    print(f"lines of evidence: {ev}")
    print("species with >=2 lines of evidence (CENP-B-box-positive):")
    print(df[df.lines_positive>=2].to_string(index=False))
    print(f"\nfull table -> {out/'cenpb_box_evidence_per_species.tsv'}")

if __name__=="__main__": main()
