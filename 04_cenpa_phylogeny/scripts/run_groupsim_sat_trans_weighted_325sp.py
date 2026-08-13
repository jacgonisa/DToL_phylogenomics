#!/usr/bin/env python3
"""
GroupSim with two phylogenetically-aware weighting schemes — Satellite vs Transposon CENP-A.

Groups (from iTOL centromere symbols):
  Satellite  = #ff006e
  Transposon = #3a86ff

Two weighting methods (both run and compared):
  HH   — Henikoff-Henikoff position-based weights (sequence-similarity proxy)
  CLADE — clade-level weights: each sequence weighted by 1/n_clade within its
           comparison group, where n_clade = number of sequences from the same
           taxa1 clade (order-level) in that group. Every clade contributes
           equally regardless of sampling depth.

Outputs (split_entropy/groupsim_sat_trans_weighted/):
  groupsim_st_weighted_gap085.tsv / .txt        — HH weighted
  groupsim_st_clade_gap085.tsv / .txt           — clade weighted
  groupsim_st_weight_comparison.tsv             — HH vs clade vs unweighted
  groupsim_st_weighted_gap085_weights.tsv       — per-sequence weights (both methods)
"""

import math, re
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import zscore as sp_zscore

BASE         = Path(__file__).parent
ALN          = BASE / "cenpa430_H3_archaea10.aligned.clipkit.325sp.fasta"
ITOL_SYMBOLS = BASE / "iTOL_bnni" / "iTOL_centromere_symbols.txt"
TAX_F        = Path("/home/jg2070/Desktop/dtol_review_August/2026_trees/annotation_centromeres/centromere_code_to_species.tsv")
UW_DIR       = BASE / "split_entropy" / "groupsim_sat_trans"
OUT_DIR      = BASE / "split_entropy" / "groupsim_sat_trans_weighted"
OUT_DIR.mkdir(parents=True, exist_ok=True)

GAP  = '-'
AA20 = list("ACDEFGHIKLMNPQRSTVWY")

# ── I/O ───────────────────────────────────────────────────────────────────────

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

def trim_columns(seqs, thr):
    ids = list(seqs.keys())
    arr = [seqs[i] for i in ids]
    n, L = len(arr), len(arr[0])
    keep = [c for c in range(L)
            if sum(1 for s in arr if s[c] == GAP) / n <= thr]
    trimmed = {sid: "".join(arr[k][c] for c in keep)
               for k, sid in enumerate(ids)}
    return trimmed, keep

def read_uw_txt(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('#') or not line.strip(): continue
        parts = line.split('\t', 2)
        try: pos = int(parts[0])
        except (ValueError, IndexError): continue
        scor = None if parts[1].strip() in ('None', '', 'NA') else float(parts[1])
        rows.append({'pos': pos, 'groupsim_uw': scor})
    return pd.DataFrame(rows)

# ── HH weights ────────────────────────────────────────────────────────────────

def henikoff_weights(sequences):
    n, L = len(sequences), len(sequences[0])
    w = np.zeros(n)
    for c in range(L):
        col = [seq[c] for seq in sequences]
        non_gap = [(i, aa) for i, aa in enumerate(col) if aa != GAP]
        if not non_gap: continue
        cnt = Counter(aa for _, aa in non_gap)
        r   = len(cnt)
        for i, aa in non_gap:
            w[i] += 1.0 / (r * cnt[aa])
    s = w.sum()
    if s > 0:
        w *= n / s
    return w

# ── Clade weights ─────────────────────────────────────────────────────────────

def seq_to_code(seq_id):
    """Extract species code from sequence ID (prefix before first digit)."""
    base = seq_id.split("_")[0].lower()
    return re.sub(r"[0-9.].*$", "", base)

INSECTS = {"Coleoptera","Diptera","Hymenoptera","Lepidoptera","Hemiptera",
           "Ephemeroptera","Odonata","Trichoptera","Neuroptera","Dermaptera",
           "Plecoptera","Psocodea","Blattodea"}
VERTEBRATES = {"Actinopterygii","Aves","Mammalia","Reptilia","Amphibia","Chondrichthyes"}
VIRIDIPLANTAE = {"Algae","Bryophyta","Dicots","Monocots","Gymnosperms"}

def taxa1_to_broad(t):
    if t in INSECTS:      return "Insects"
    if t in VERTEBRATES:  return "Vertebrates"
    if t in VIRIDIPLANTAE: return "Viridiplantae"
    if t == "Fungi":      return "Fungi"
    return "Other_Invertebrates"

def clade_weights(seq_ids, group_ids_list, tax_df):
    """
    Weight each sequence by 1 / (number of sequences from the same broad clade
    within the same comparison group): Insects, Vertebrates, Viridiplantae,
    Fungi, Other_Invertebrates. Every major clade contributes equally.
    Weights normalised within each group so mean = 1.
    """
    code_to_taxa1 = dict(zip(
        tax_df["fasta_base"].str.lower().str.replace(r"[0-9.].*$", "", regex=True),
        tax_df["taxa1"]
    ))
    broad_arr = np.array([taxa1_to_broad(code_to_taxa1.get(seq_to_code(sid), "Unknown"))
                          for sid in seq_ids])
    w = np.ones(len(seq_ids))
    for idx in group_ids_list:
        clades = broad_arr[idx]
        counts = Counter(clades)
        for i, c in zip(idx, clades):
            w[i] = 1.0 / counts[c]
        grp_w = w[idx]
        w[idx] = grp_w / grp_w.mean()
    return w

# ── Weighted scoring ──────────────────────────────────────────────────────────

def weighted_freqs(col, indices, weights, group_gap_cutoff=0.5):
    g_col = [col[i] for i in indices]
    if sum(1 for c in g_col if c == GAP) / len(g_col) > group_gap_cutoff:
        return None
    pairs = [(aa, weights[i]) for i, aa in zip(indices, g_col) if aa != GAP]
    if not pairs: return None
    wsum = sum(wi for _, wi in pairs)
    if wsum == 0: return None
    f = {}
    for aa, wi in pairs:
        f[aa] = f.get(aa, 0.0) + wi / wsum
    return f

def groupsim_weighted_score(col, group_idx_list, weights,
                             col_gap_cutoff=0.85, group_gap_cutoff=0.5):
    n_gap = sum(1 for c in col if c == GAP)
    if n_gap / len(col) > col_gap_cutoff:
        return None
    gfreqs = [weighted_freqs(col, idx, weights, group_gap_cutoff)
              for idx in group_idx_list]
    if any(gf is None for gf in gfreqs):
        return None
    within = [sum(f**2 for f in gf.values()) for gf in gfreqs]
    mean_within = np.mean(within)
    between = []
    for g1 in range(len(gfreqs)):
        for g2 in range(g1+1, len(gfreqs)):
            all_aa = set(gfreqs[g1]) | set(gfreqs[g2])
            between.append(sum(gfreqs[g1].get(a,0)*gfreqs[g2].get(a,0)
                               for a in all_aa))
    mean_between = np.mean(between) if between else 0.0
    return (1.0 - n_gap/len(col)) * (mean_within - mean_between)

# ── JS divergence + window smoothing (same as groupsim.py) ───────────────────

def js_div(col):
    bg = [0.074,0.052,0.045,0.054,0.025,0.034,0.054,0.074,0.026,0.068,
          0.099,0.058,0.025,0.047,0.039,0.057,0.051,0.013,0.032,0.073]
    pc = 0.001
    cnt = Counter(aa for aa in col if aa != GAP)
    n   = sum(cnt.values()) or 1
    fc1 = [cnt.get(aa,0)/n + pc for aa in AA20]
    s1  = sum(fc1); fc1 = [x/s1 for x in fc1]
    fc2 = [x+pc for x in bg]; s2 = sum(fc2); fc2 = [x/s2 for x in fc2]
    r   = [0.5*fc1[i]+0.5*fc2[i] for i in range(20)]
    d   = sum(fc1[i]*math.log2(fc1[i]/r[i]) if fc1[i]>0 and r[i]>0 else 0
              for i in range(20))
    d  += sum(fc2[i]*math.log2(fc2[i]/r[i]) if fc2[i]>0 and r[i]>0 else 0
              for i in range(20))
    return (1 - sum(1 for c in col if c==GAP)/len(col)) * d/2

def norm_01(vals):
    v = [x for x in vals if x is not None]
    if not v: return vals
    lo, hi = min(v), max(v)
    r = hi - lo
    return [None if x is None else (1.0 if r==0 else (x-lo)/r) for x in vals]

def window_correct(raw, cols, window=3, lam=0.7):
    cons   = [js_div(col) if s is not None else None
              for s, col in zip(raw, cols)]
    cons_n = norm_01(cons)
    out = []
    for i, s in enumerate(raw):
        if s is None: out.append(None); continue
        lo, hi = max(i-window,0), min(i+window+1, len(raw))
        terms  = [cons_n[j] for j in range(lo,hi)
                  if j!=i and cons_n[j] is not None and cons_n[j]>=0]
        wc = sum(terms)/len(terms) if terms else 0.0
        out.append((1-lam)*wc + lam*s)
    return out

# ── Residue detail string ─────────────────────────────────────────────────────

def format_col_wt(col, group_idx_list, weights, group_names):
    parts = []
    for name, idx in zip(group_names, group_idx_list):
        pairs = [(col[i], weights[i]) for i in idx if col[i] != GAP]
        if not pairs: parts.append(f"{name}: -"); continue
        wsum = sum(w for _,w in pairs)
        f = {}
        for aa,w in pairs: f[aa] = f.get(aa,0.0) + w/wsum
        tokens = [f"{aa}({100*v:.0f}%)"
                  for aa,v in sorted(f.items(), key=lambda x:-x[1]) if v>0.005]
        parts.append(f"{name}: " + " ".join(tokens))
    return "  |  ".join(parts)

# ── Main ──────────────────────────────────────────────────────────────────────

def run(gap=0.85):
    tag = str(gap).replace('.','')
    print(f"\n── Gap threshold {gap} ──")

    # 1. Parse group membership from iTOL
    sat_tips, trans_tips = [], []
    for line in ITOL_SYMBOLS.read_text().splitlines():
        if line.startswith('#') or not line.strip() or '\t' not in line: continue
        parts = line.split('\t')
        if len(parts) < 4: continue
        tip, col = parts[0], parts[3]
        if   col == '#ff006e': sat_tips.append(tip)
        elif col == '#3a86ff': trans_tips.append(tip)
    print(f"  Satellite: {len(sat_tips)}  Transposon: {len(trans_tips)}")

    # 2. Load alignment, filter to sat+trans only
    seqs_all = read_fasta(ALN)
    keep_ids = set(sat_tips) | set(trans_tips)
    seqs_st  = {k: v for k,v in seqs_all.items() if k in keep_ids}
    print(f"  Sequences retained: {len(seqs_st)}")

    # 3. Trim on sat+trans sequences only
    trimmed, _keep = trim_columns(seqs_st, gap)
    L = len(next(iter(trimmed.values())))
    print(f"  Columns after trim: {L}")

    # 4. Build arrays
    sat_in   = [s for s in sat_tips   if s in trimmed]
    trans_in = [s for s in trans_tips if s in trimmed]
    ids_ord  = sat_in + trans_in
    seqlist  = [trimmed[sid] for sid in ids_ord]

    idx_sat   = list(range(len(sat_in)))
    idx_trans = list(range(len(sat_in), len(ids_ord)))
    group_idx = [idx_sat, idx_trans]

    # 5. Load taxonomy for clade weights
    tax_df = pd.read_csv(TAX_F, sep='\t')

    # 5a. HH weights (computed on all sat+trans together)
    print("  computing HH weights …")
    w_hh = henikoff_weights(seqlist)
    print(f"  HH weight range: {w_hh.min():.4f} – {w_hh.max():.4f}  mean: {w_hh.mean():.4f}")

    # 5b. Clade weights (broad groups: Insects / Vertebrates / Viridiplantae / Fungi / Other_Invertebrates)
    print("  computing clade weights …")
    w_cl = clade_weights(ids_ord, group_idx, tax_df)
    broad_arr = np.array([taxa1_to_broad(
        dict(zip(tax_df["fasta_base"].str.lower().str.replace(r"[0-9.].*$","",regex=True),
                 tax_df["taxa1"])).get(seq_to_code(sid), "Unknown"))
        for sid in ids_ord])
    for grp, idx in zip(['Satellite','Transposon'], group_idx):
        cnt = Counter(broad_arr[idx])
        print(f"  {grp} clades: {dict(cnt)}")
    print(f"  Clade weight range: {w_cl.min():.4f} – {w_cl.max():.4f}  mean: {w_cl.mean():.4f}")

    # Save weights
    w_df = pd.DataFrame({
        'seq_id':      ids_ord,
        'group':       ['Satellite']*len(sat_in) + ['Transposon']*len(trans_in),
        'broad_clade': broad_arr,
        'hh_weight':   w_hh,
        'clade_weight': w_cl,
    })
    w_df.to_csv(OUT_DIR / f"groupsim_st_weighted_gap{tag}_weights.tsv",
                sep='\t', index=False, float_format='%.6f')

    cols = [[seq[c] for seq in seqlist] for c in range(L)]

    def score_and_save(w, label, out_prefix):
        raw      = [groupsim_weighted_score(cols[c], group_idx, w,
                                            col_gap_cutoff=gap, group_gap_cutoff=0.5)
                    for c in range(L)]
        corrected = window_correct(raw, cols, window=3, lam=0.7)
        scores_n  = norm_01(corrected)
        valid_pairs = [(i,s) for i,s in enumerate(scores_n) if s is not None]
        vals   = np.array([s for _,s in valid_pairs])
        zs_arr = sp_zscore(vals) if len(vals) > 1 else np.zeros(len(vals))
        zs_map = {i: z for (i,_),z in zip(valid_pairs, zs_arr)}

        rows = []
        for c in range(L):
            s = scores_n[c]; z = zs_map.get(c)
            rows.append({'pos': c+1, f'groupsim_{label}': s,
                         f'z_{label}': z, f'sig_{label}': int(z is not None and z >= 2.0)})
        df = pd.DataFrame(rows)
        df.to_csv(OUT_DIR / f"{out_prefix}_gap{tag}.tsv", sep='\t', index=False, float_format='%.6f')
        n_sig = int(df[f'sig_{label}'].sum())
        print(f"  [{label}] sig positions (z≥2): {n_sig}")
        return df, scores_n, zs_map

    # 6. Score with both weight schemes
    print("  scoring with HH weights …")
    df_hh, sn_hh, zm_hh = score_and_save(w_hh, "hh",    "groupsim_st_weighted")
    print("  scoring with clade weights …")
    df_cl, sn_cl, zm_cl = score_and_save(w_cl, "clade", "groupsim_st_clade")

    # 7. Merge comparison
    comp = df_hh.merge(df_cl, on='pos')
    uw_path = UW_DIR / f"groupsim_gap{tag}.tsv"
    if uw_path.exists():
        uw = pd.read_csv(uw_path, sep='\t')[['pos','groupsim_score']].rename(
            columns={'groupsim_score': 'groupsim_unweighted'})
        comp = comp.merge(uw, on='pos', how='left')
    comp.to_csv(OUT_DIR / f"groupsim_st_weight_comparison_gap{tag}.tsv",
                sep='\t', index=False, float_format='%.6f')
    print(f"  saved comparison: groupsim_st_weight_comparison_gap{tag}.tsv")

    return df_hh, sn_hh, zm_hh, df_cl, sn_cl, zm_cl, cols, group_idx, w_hh, w_cl, ids_ord, sat_in, trans_in

# ── Comparison & summary ──────────────────────────────────────────────────────

def make_comparison(tsv_wt, uw_txt, top_n=15):
    wt = tsv_wt.dropna(subset=['groupsim_weighted']).copy()
    wt = wt.sort_values('groupsim_weighted', ascending=False).reset_index(drop=True)
    wt['rank_wt'] = wt.index + 1

    try:
        uw = read_uw_txt(uw_txt)
        uw = uw.dropna(subset=['groupsim_uw']).copy()
        uw = uw.sort_values('groupsim_uw', ascending=False).reset_index(drop=True)
        uw['rank_uw'] = uw.index + 1
        merged = wt.merge(uw[['pos','groupsim_uw','rank_uw']], on='pos', how='left')
    except Exception as e:
        print(f"  Warning: {e}")
        merged = wt.copy()
        merged['groupsim_uw'] = None
        merged['rank_uw']     = None

    return merged.sort_values('rank_wt').head(top_n)[
        ['rank_wt','pos','groupsim_weighted','z_score','sig_z2',
         'groupsim_uw','rank_uw']]


if __name__ == "__main__":
    df_hh, sn_hh, zm_hh, df_cl, sn_cl, zm_cl, cols, group_idx, w_hh, w_cl, ids_ord, sat_in, trans_in = run(gap=0.85)
    tsv_df = df_hh.rename(columns={'groupsim_hh':'groupsim_weighted','z_hh':'z_score','sig_hh':'sig_z2'})
    scores_n, zs_map, w = sn_hh, zm_hh, w_hh

    # Print top 10
    valid = tsv_df.dropna(subset=['groupsim_weighted']) \
                  .sort_values('groupsim_weighted', ascending=False)
    print(f"\n{'─'*80}")
    print(f"  Top 10 — weighted GroupSim Sat vs Trans (gap 0.85)")
    print(f"{'─'*80}")
    print(f"  {'Rk':>3}  {'Pos':>4}  {'Score':>7}  {'Z':>6}  Satellite (wt%)  |  Transposon (wt%)")
    print(f"{'─'*80}")
    for rank, (_, row) in enumerate(valid.head(10).iterrows(), 1):
        c = int(row['pos']) - 1
        detail = format_col_wt(cols[c], group_idx, w, ['Sat','Trans'])
        parts  = detail.split('  |  ')
        s_str  = parts[0].replace('Sat: ','')[:38]
        t_str  = parts[1].replace('Trans: ','') if len(parts)>1 else ''
        z_str  = f"{row['z_score']:6.2f}" if pd.notna(row['z_score']) else "    NA"
        sig    = ' *' if row['sig_z2']==1 else ''
        print(f"  {rank:>3}  {int(row['pos']):>4}  {row['groupsim_weighted']:7.3f}"
              f"  {z_str}  {s_str:<38}  {t_str}{sig}")

    # Comparison table
    uw_txt = UW_DIR / "groupsim_sat_trans_gap085.txt"
    comp   = make_comparison(tsv_df, uw_txt)
    comp_path = OUT_DIR / "groupsim_st_weighted_comparison.tsv"
    comp.to_csv(comp_path, sep='\t', index=False, float_format='%.4f')

    print(f"\n  Rank comparison (wt | pos | wt-score | uw-score | uw-rank | Δrank):")
    print(f"  {'─'*65}")
    for _, row in comp.iterrows():
        uw_r   = int(row['rank_uw'])   if pd.notna(row['rank_uw'])   else -1
        uw_s   = f"{row['groupsim_uw']:7.3f}" if pd.notna(row['groupsim_uw']) else "     NA"
        delta  = f"{uw_r - int(row['rank_wt']):+d}" if uw_r > 0 else "  NA"
        ur_str = f"{uw_r:>6}" if uw_r > 0 else "    NA"
        sig    = " *" if row['sig_z2']==1 else "  "
        print(f"  {int(row['rank_wt']):>4}  {int(row['pos']):>4}  "
              f"{row['groupsim_weighted']:7.3f}  {uw_s}  {ur_str}  {delta}{sig}")

    print(f"\nDone. Outputs: {OUT_DIR}")
