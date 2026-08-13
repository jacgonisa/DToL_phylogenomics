#!/usr/bin/env python3
"""
Per-class GroupSim (vertebrate class CENPA vs all H3) + within-class Shannon entropy.
Outputs: vertebrate_antibody/groupsim_{Class}_vs_H3.tsv for each vertebrate class.
Also outputs groupsim_Vertebrata_vs_H3.tsv (all vertebrate CENPA pooled).

Run after extract_vertebrate_cenpa.py.
"""
import math, re
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
from scipy.stats import zscore as sp_zscore

BASE    = Path(__file__).parent
ALN     = BASE / "cenpa430_H3_archaea10.aligned.clipkit.325sp.fasta"
GRP_F   = BASE / "split_entropy" / "groupsim" / "groups_gap085.txt"
TAX_F   = Path("/home/jg2070/Desktop/dtol_review_August/2026_trees/annotation_centromeres/centromere_code_to_species.tsv")
OUT_DIR = BASE / "vertebrate_antibody"
OUT_DIR.mkdir(exist_ok=True)

GAP  = "-"
AA20 = list("ACDEFGHIKLMNPQRSTVWY")
COL_GAP_CUTOFF  = 0.85  # column gap fraction for trimming (match main analysis)
GRP_GAP_CUTOFF  = 0.50  # within-group gap fraction to skip a column for that group
ARCHAEA_IDS = {
    "OLS22332.1","OLS24873.1","OLS21974.1","KKK41979.1","KXH71038.1",
    "OLS18261.1","OLS16336.1","BAD86478.1","OIO61677.1","OIO41945.1",
}

def read_fasta(path):
    seqs, cur, buf = {}, None, []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            if cur is not None: seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]; buf = []
        else: buf.append(line.upper())
    if cur is not None: seqs[cur] = "".join(buf)
    return seqs

def seq_to_code(sid):
    return re.sub(r"[0-9.].*$", "", sid.split("_")[0].lower())

def trim_columns(seqs, thr):
    ids = list(seqs); arr = [seqs[i] for i in ids]; n = len(arr); L = len(arr[0])
    keep = [c for c in range(L) if sum(1 for s in arr if s[c]==GAP)/n <= thr]
    return {sid: "".join(arr[k][c] for c in keep) for k,sid in enumerate(ids)}, keep

def aa_freqs(col):
    """Unweighted AA frequency dict (gaps excluded)."""
    aa = [c for c in col if c != GAP]
    if not aa: return None
    n = len(aa)
    return {a: cnt/n for a,cnt in Counter(aa).items()}

def groupsim_score(col, idx_a, idx_b, col_gap_cut=COL_GAP_CUTOFF, grp_gap_cut=GRP_GAP_CUTOFF):
    """Raw GroupSim score: (1-gap_frac) * (mean_within - mean_between)."""
    gap_frac = sum(1 for c in col if c==GAP) / len(col)
    if gap_frac > col_gap_cut: return None
    fa = aa_freqs([col[i] for i in idx_a])
    fb = aa_freqs([col[i] for i in idx_b])
    # skip if either group is too gappy
    if fa is None or sum(1 for i in idx_a if col[i]==GAP)/len(idx_a) > grp_gap_cut: return None
    if fb is None or sum(1 for i in idx_b if col[i]==GAP)/len(idx_b) > grp_gap_cut: return None
    within  = (sum(f**2 for f in fa.values()) + sum(f**2 for f in fb.values())) / 2
    between = sum(fa.get(a,0)*fb.get(a,0) for a in set(fa)|set(fb))
    return (1 - gap_frac) * (within - between)

def shannon_entropy(col):
    aa = [c for c in col if c != GAP]
    if not aa: return None
    cnt = Counter(aa); n = len(aa)
    return -sum((cnt[a]/n)*math.log2(cnt[a]/n) for a in cnt if cnt[a]>0)

def consensus_aa(col):
    aa = [c for c in col if c != GAP]
    if not aa: return "-"
    return Counter(aa).most_common(1)[0][0]

def js_div_conservation(col):
    """JSD-based per-column conservation vs background, for windowing correction."""
    bg = [0.074,0.052,0.045,0.054,0.025,0.034,0.054,0.074,0.026,0.068,
          0.099,0.058,0.025,0.047,0.039,0.057,0.051,0.013,0.032,0.073]
    pc = 0.001
    cnt = Counter(aa for aa in col if aa!=GAP); n = sum(cnt.values()) or 1
    fc1 = [cnt.get(aa,0)/n+pc for aa in AA20]; s1=sum(fc1); fc1=[x/s1 for x in fc1]
    fc2 = [x+pc for x in bg]; s2=sum(fc2); fc2=[x/s2 for x in fc2]
    r   = [0.5*fc1[i]+0.5*fc2[i] for i in range(20)]
    d   = sum(fc1[i]*math.log2(fc1[i]/r[i]) if fc1[i]>0 else 0 for i in range(20))
    d  += sum(fc2[i]*math.log2(fc2[i]/r[i]) if fc2[i]>0 else 0 for i in range(20))
    return (1 - sum(1 for c in col if c==GAP)/len(col)) * d/2

def norm_01(vals):
    v = [x for x in vals if x is not None]
    if not v: return vals
    lo, hi = min(v), max(v); r = hi-lo
    return [None if x is None else (1.0 if r==0 else (x-lo)/r) for x in vals]

def window_correct(raw, cols, window=3, lam=0.7):
    """Smooth GroupSim scores with local JSD conservation context."""
    cons_n = norm_01([js_div_conservation(col) for col in cols])
    out = []
    for i, s in enumerate(raw):
        if s is None: out.append(None); continue
        lo, hi = max(i-window,0), min(i+window+1,len(raw))
        terms = [cons_n[j] for j in range(lo,hi) if j!=i and cons_n[j] is not None]
        wc = sum(terms)/len(terms) if terms else 0
        out.append((1-lam)*wc + lam*s)
    return out

def run_class(name, class_ids, h3_ids, trimmed):
    ids_ord  = [s for s in class_ids+h3_ids if s in trimmed]
    cls_in   = [s for s in class_ids if s in trimmed]
    h3_in    = [s for s in h3_ids if s in trimmed]
    if not cls_in:
        print(f"  {name}: no sequences in trimmed MSA, skipping"); return None

    seqlist  = [trimmed[s] for s in ids_ord]
    idx_a    = list(range(len(cls_in)))
    idx_b    = list(range(len(cls_in), len(ids_ord)))
    L        = len(seqlist[0])
    cols     = [[seq[c] for seq in seqlist] for c in range(L)]
    cls_seqs = [trimmed[s] for s in cls_in]

    raw       = [groupsim_score(cols[c], idx_a, idx_b) for c in range(L)]
    corrected = window_correct(raw, cols)
    scores_n  = norm_01(corrected)

    valid = [(i,s) for i,s in enumerate(scores_n) if s is not None]
    if len(valid) > 1:
        zvals = sp_zscore([s for _,s in valid])
        zmap  = {i:z for (i,_),z in zip(valid, zvals)}
    else:
        zmap = {i:0.0 for i,_ in valid}

    rows = []
    for c in range(L):
        col_cls = [cls_seqs[k][c] for k in range(len(cls_seqs))]
        ent   = shannon_entropy(col_cls)
        con   = 1 - (ent/math.log2(min(20, len([x for x in col_cls if x!=GAP])+1)) if ent is not None else 0)
        gfrac = sum(1 for aa in col_cls if aa==GAP) / len(col_cls)
        rows.append({
            "col_idx":         c,
            "groupsim_z":      zmap.get(c),
            "groupsim_norm":   scores_n[c],
            "within_entropy":  ent,
            "conservation":    con if ent is not None else None,
            "consensus_aa":    consensus_aa(col_cls),
            "gap_frac_class":  round(gfrac, 4),
            "n_seqs":          len(cls_in),
        })
    df = pd.DataFrame(rows)
    sig = (df["groupsim_z"] >= 2.0).sum()
    print(f"  {name}: {len(cls_in)} seqs, {L} cols, {sig} sig positions (z≥2)")
    out = OUT_DIR / f"groupsim_{name}_vs_H3.tsv"
    df.to_csv(out, sep="\t", index=False, float_format="%.6f")
    return df

# ── main ──────────────────────────────────────────────────────────────────────
seqs_all = {k:v for k,v in read_fasta(ALN).items() if k not in ARCHAEA_IDS}

grps = {}
for line in GRP_F.read_text().splitlines():
    if ":" not in line: continue
    g, ids = line.strip().split(":",1)
    grps[g] = [x for x in ids.split(",") if x in seqs_all]
h3_ids = grps["H3"]

tax_df   = pd.read_csv(TAX_F, sep="\t")
code_map = {re.sub(r"[0-9.].*$","",r.fasta_base.lower()): r.taxa1 for _,r in tax_df.iterrows()}

cenpa_ids = grps["CENPA"]
by_class = {}
for sid in cenpa_ids:
    cls = code_map.get(seq_to_code(sid))
    if cls in {"Aves","Mammalia","Actinopterygii","Reptilia"}:
        by_class.setdefault(cls, []).append(sid)

trimmed, keep = trim_columns(seqs_all, COL_GAP_CUTOFF)
print(f"MSA: {len(seqs_all)} seqs, {len(keep)} columns after {COL_GAP_CUTOFF} gap trim")
print(f"H3: {len(h3_ids)}  Vertebrate CENPA: {sum(len(v) for v in by_class.values())}")

for cls in sorted(by_class):
    run_class(cls, by_class[cls], h3_ids, trimmed)

# pooled vertebrates
all_vert = [s for ids in by_class.values() for s in ids]
run_class("Vertebrata", all_vert, h3_ids, trimmed)

print("\nDone → vertebrate_antibody/groupsim_*_vs_H3.tsv")
