#!/usr/bin/env python3
"""
GroupSim with Henikoff-Henikoff sequence reweighting (CENP-A vs H3).

Standard GroupSim (Pazos et al. 2006) treats every sequence as an independent
observation, so over-represented clades (insects, Viridiplantae) inflate scores.
Henikoff-Henikoff (1994) position-based weights downweight sequences from dense
clades before computing within- and between-group similarity.

Scoring (identity matrix):
  within-group score  = weighted homozygosity = sum_a (w_freq_a)^2
  between-group score = sum_a w_freq_a^g1 * w_freq_a^g2
  raw score           = gap_weight * (mean_within - between)
  final score         = window-corrected (same JS-divergence window as original)

Outputs (split_entropy/groupsim_weighted/):
  groupsim_weighted_gap085.tsv / .txt   — scores + residue details
  groupsim_weighted_comparison.tsv      — side-by-side unweighted vs weighted
  groupsim_weighted_gap085_weights.tsv  — HH weights per sequence
"""

import math
import sys
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import zscore as sp_zscore

BASE     = Path(__file__).parent
H3_FASTA = BASE / "H3_all.aligned.fasta"
ALN      = BASE / "cenpa430_H3_archaea10.aligned.clipkit.325sp.fasta"
UNWEIGHTED_TXT = BASE / "split_entropy" / "groupsim" / "groupsim_gap085.txt"

OUT_DIR  = BASE / "split_entropy" / "groupsim_weighted"
OUT_DIR.mkdir(parents=True, exist_ok=True)

ARCHAEA_IDS = {
    "OLS22332.1","OLS24873.1","OLS21974.1","KKK41979.1","KXH71038.1",
    "OLS18261.1","OLS16336.1","BAD86478.1","OIO61677.1","OIO41945.1",
}
GAP = '-'

AA20 = list("ACDEFGHIKLMNPQRSTVWY")

# ── I/O helpers ───────────────────────────────────────────────────────────────

def read_fasta(path):
    seqs, cur, buf = {}, None, []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur is not None:
                seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]; buf = []
        else:
            buf.append(line.upper())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs

def trim_columns(seqs, gap_threshold=0.85):
    ids  = list(seqs.keys())
    arr  = [seqs[i] for i in ids]
    n, L = len(arr), len(arr[0])
    keep = [c for c in range(L)
            if sum(1 for s in arr if s[c] == GAP) / n <= gap_threshold]
    trimmed = {sid: "".join(arr[k][c] for c in keep)
               for k, sid in enumerate(ids)}
    return trimmed, keep

# ── Henikoff-Henikoff weights ─────────────────────────────────────────────────

def henikoff_weights(sequences):
    """
    Henikoff & Henikoff (1994) position-based sequence weights.
    Returns array of length n, normalised to sum to n.
    """
    n, L = len(sequences), len(sequences[0])
    w = np.zeros(n)
    for c in range(L):
        col = [seq[c] for seq in sequences]
        non_gap = [(i, aa) for i, aa in enumerate(col) if aa != GAP]
        if not non_gap:
            continue
        cnt = Counter(aa for _, aa in non_gap)
        r   = len(cnt)
        for i, aa in non_gap:
            w[i] += 1.0 / (r * cnt[aa])
    s = w.sum()
    if s > 0:
        w *= n / s   # normalise: sum = n (effective number of sequences)
    return w

# ── Weighted scoring ──────────────────────────────────────────────────────────

def weighted_freqs(col, indices, weights, group_gap_cutoff=0.5, with_gaps=False):
    """Weighted AA frequency dict for one group at one column, or None if too gappy.

    with_gaps=True: treat '-' as a 21st character; no group-level gap filter applied.
    """
    g_col = [col[i] for i in indices]
    n_gap = sum(1 for c in g_col if c == GAP)
    if not with_gaps and n_gap / len(g_col) > group_gap_cutoff:
        return None
    pairs = [(aa, weights[i]) for i, aa in zip(indices, g_col)
             if with_gaps or aa != GAP]
    if not pairs:
        return None
    wsum = sum(wi for _, wi in pairs)
    if wsum == 0:
        return None
    f = {}
    for aa, wi in pairs:
        f[aa] = f.get(aa, 0.0) + wi / wsum
    return f

def groupsim_weighted_score(col, group_idx_list, weights,
                            col_gap_cutoff=0.85, group_gap_cutoff=0.5,
                            with_gaps=False):
    """
    Raw weighted GroupSim score at one column.
    Returns None if column fails gap filters.

    with_gaps=True: gap treated as 21st character; no group-level gap filter;
    no gap downweighting (gap_weight fixed at 1.0).
    """
    n_gap = sum(1 for c in col if c == GAP)
    if n_gap / len(col) > col_gap_cutoff:
        return None

    gfreqs = [weighted_freqs(col, idx, weights, group_gap_cutoff, with_gaps)
              for idx in group_idx_list]
    if any(gf is None for gf in gfreqs):
        return None

    # Within-group: homozygosity (expected identity between 2 random draws)
    within = [sum(f**2 for f in gf.values()) for gf in gfreqs]
    mean_within = np.mean(within)

    # Between-group: dot product of frequency vectors
    between = []
    for g1 in range(len(gfreqs)):
        for g2 in range(g1 + 1, len(gfreqs)):
            all_aa = set(gfreqs[g1]) | set(gfreqs[g2])
            between.append(sum(gfreqs[g1].get(a, 0.0) * gfreqs[g2].get(a, 0.0)
                               for a in all_aa))
    mean_between = np.mean(between) if between else 0.0

    # When gaps are informative characters, no additional downweighting
    gap_weight = 1.0 if with_gaps else (1.0 - n_gap / len(col))
    return gap_weight * (mean_within - mean_between)

# ── JS divergence for window smoothing (matches original groupsim.py) ─────────

def js_divergence_col(col):
    bg = [0.074,0.052,0.045,0.054,0.025,0.034,0.054,0.074,0.026,0.068,
          0.099,0.058,0.025,0.047,0.039,0.057,0.051,0.013,0.032,0.073]
    pc = 0.001
    cnt = Counter(aa for aa in col if aa != GAP)
    n   = sum(cnt.values()) or 1
    fc1 = [cnt.get(aa, 0) / n + pc for aa in AA20]
    s1  = sum(fc1); fc1 = [x / s1 for x in fc1]
    fc2 = [x + pc for x in bg]; s2 = sum(fc2); fc2 = [x / s2 for x in fc2]
    r   = [0.5 * fc1[i] + 0.5 * fc2[i] for i in range(20)]
    d   = 0.0
    for i in range(20):
        if fc1[i] > 0 and r[i] > 0: d += fc1[i] * math.log2(fc1[i] / r[i])
        if fc2[i] > 0 and r[i] > 0: d += fc2[i] * math.log2(fc2[i] / r[i])
    gap_w = 1.0 - sum(1 for c in col if c == GAP) / len(col)
    return gap_w * d / 2

def norm_01(vals):
    valid = [v for v in vals if v is not None]
    if not valid: return vals
    lo, hi = min(valid), max(valid)
    r = hi - lo
    return [None if v is None else (1.0 if r == 0 else (v - lo) / r)
            for v in vals]

def window_correct(raw_scores, cols, window=3, lam=0.7):
    """Apply JS-divergence window correction (same as original groupsim.py)."""
    cons = [js_divergence_col(col) if s is not None else None
            for s, col in zip(raw_scores, cols)]
    cons_n = norm_01(cons)
    out = []
    for i, s in enumerate(raw_scores):
        if s is None:
            out.append(None)
            continue
        lo, hi = max(i - window, 0), min(i + window + 1, len(raw_scores))
        terms = [cons_n[j] for j in range(lo, hi)
                 if j != i and cons_n[j] is not None and cons_n[j] >= 0]
        wc = sum(terms) / len(terms) if terms else 0.0
        out.append((1.0 - lam) * wc + lam * s)
    return out

# ── Residue detail string (weighted percentages) ──────────────────────────────

def format_col_weighted(col, group_idx_list, weights, group_names, with_gaps=False):
    parts = []
    for name, idx in zip(group_names, group_idx_list):
        all_pairs = [(col[i], weights[i]) for i in idx]
        if with_gaps:
            pairs = all_pairs
        else:
            pairs = [(aa, w) for aa, w in all_pairs if aa != GAP]
        if not pairs:
            parts.append(f"{name}: -"); continue
        wsum = sum(w for _, w in pairs)
        f = {}
        for aa, w in pairs:
            f[aa] = f.get(aa, 0.0) + w / wsum
        tokens = [f"{aa}({100*v:.0f}%)"
                  for aa, v in sorted(f.items(), key=lambda x: (-x[1], x[0]))
                  if v > 0.005]
        parts.append(f"{name}: " + " ".join(tokens))
    return "  |  ".join(parts)

# ── Main ──────────────────────────────────────────────────────────────────────

def run(gap=0.85, with_gaps=False):
    tag     = str(gap).replace(".", "")
    gap_sfx = "_gap21" if with_gaps else ""
    print(f"\n── Gap threshold {gap}  {'(gap as 21st char)' if with_gaps else ''} ──")

    # 1. Read alignment, assign groups
    seqs_all = read_fasta(ALN)
    seqs     = {k: v for k, v in seqs_all.items() if k not in ARCHAEA_IDS}
    h3_pool  = set(read_fasta(H3_FASTA).keys())

    groups = {"CENPA": [], "H3": []}
    for sid in seqs:
        groups["H3" if sid in h3_pool else "CENPA"].append(sid)
    print(f"  CENPA: {len(groups['CENPA'])}  H3: {len(groups['H3'])}")

    # 2. Trim columns
    trimmed, _keep = trim_columns(seqs, gap)
    print(f"  columns after trim: {len(next(iter(trimmed.values())))}")

    # 3. Build arrays
    ids_ordered = groups["CENPA"] + groups["H3"]
    seqlist     = [trimmed[sid] for sid in ids_ordered if sid in trimmed]
    ids_ordered = [sid for sid in ids_ordered if sid in trimmed]
    L = len(seqlist[0])

    idx_cenpa = [i for i, sid in enumerate(ids_ordered) if sid in set(groups["CENPA"])]
    idx_h3    = [i for i, sid in enumerate(ids_ordered) if sid in set(groups["H3"])]
    group_idx = [idx_cenpa, idx_h3]

    # 4. Henikoff-Henikoff weights
    print("  computing HH weights …")
    w = henikoff_weights(seqlist)
    print(f"  weight range: {w.min():.4f} – {w.max():.4f}  mean: {w.mean():.4f}")

    # Save weights
    w_df = pd.DataFrame({"seq_id": ids_ordered,
                          "group": ["CENPA" if i in set(idx_cenpa) else "H3"
                                    for i in range(len(ids_ordered))],
                          "hh_weight": w})
    w_df.to_csv(OUT_DIR / f"groupsim_weighted_gap{tag}_weights.tsv",
                sep="\t", index=False, float_format="%.6f")

    # 5. Score each column
    cols = [[seq[c] for seq in seqlist] for c in range(L)]
    raw  = [groupsim_weighted_score(cols[c], group_idx, w,
                                    col_gap_cutoff=gap, group_gap_cutoff=0.5,
                                    with_gaps=with_gaps)
            for c in range(L)]

    # 6. Window correction + normalise
    corrected = window_correct(raw, cols, window=3, lam=0.7)
    scores_n  = norm_01(corrected)

    # 7. Z-scores and significance
    valid   = [(i, s) for i, s in enumerate(scores_n) if s is not None]
    vals    = np.array([s for _, s in valid])
    zs_vals = sp_zscore(vals) if len(vals) > 1 else np.zeros(len(vals))
    zs_map  = {i: z for (i, _), z in zip(valid, zs_vals)}

    # 8. Write TSV
    tsv_rows = []
    for c in range(L):
        s  = scores_n[c]
        z  = zs_map.get(c)
        sig = int(z is not None and z >= 2.0)
        tsv_rows.append({"pos": c, "groupsim_weighted": s,
                         "z_score": z, "sig_z2": sig})
    tsv_df = pd.DataFrame(tsv_rows)
    tsv_path = OUT_DIR / f"groupsim_weighted_gap{tag}{gap_sfx}.tsv"
    tsv_df.to_csv(tsv_path, sep="\t", index=False, float_format="%.6f")
    print(f"  saved: {tsv_path.name}")

    # 9. Write text (residue details)
    mode_str = "gap as 21st char" if with_gaps else "gap excluded"
    txt_lines = [
        f"# Weighted GroupSim (HH weights) — gap threshold {gap} — {mode_str}",
        f"# CENPA: {len(idx_cenpa)} seqs  H3: {len(idx_h3)} seqs",
        f"# col_num\tscore\tcolumn_residues",
    ]
    for c in range(L):
        s = scores_n[c]
        s_str = f"{s:.6f}" if s is not None else "None"
        details = format_col_weighted(cols[c], group_idx, w,
                                       ["CENPA", "H3"], with_gaps=with_gaps)
        txt_lines.append(f"{c}\t{s_str}\t{details}")
    txt_path = OUT_DIR / f"groupsim_weighted_gap{tag}{gap_sfx}.txt"
    txt_path.write_text("\n".join(txt_lines) + "\n")
    print(f"  saved: {txt_path.name}")

    return tsv_df, scores_n, zs_map, cols, group_idx, w, ids_ordered

# ── Comparison table ──────────────────────────────────────────────────────────

def read_unweighted_txt(path):
    """Parse groupsim .txt output: col_num, score, residue_detail."""
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t", 2)
        if len(parts) < 2:
            continue
        pos  = int(parts[0])
        scor = None if parts[1] in ("None", "", "NA") else float(parts[1])
        rows.append({"pos": pos, "groupsim_uw": scor})
    df = pd.DataFrame(rows)
    return df

def make_comparison(tsv_weighted, unweighted_txt=UNWEIGHTED_TXT, top_n=15):
    """Build side-by-side comparison of top positions."""
    wt = tsv_weighted.dropna(subset=["groupsim_weighted"]).copy()
    wt = wt.sort_values("groupsim_weighted", ascending=False).reset_index(drop=True)
    wt["rank_weighted"] = wt.index + 1

    try:
        uw = read_unweighted_txt(unweighted_txt)
        uw = uw.dropna(subset=["groupsim_uw"]).copy()
        uw = uw.sort_values("groupsim_uw", ascending=False).reset_index(drop=True)
        uw["rank_unweighted"] = uw.index + 1
        merged = wt.merge(uw[["pos", "groupsim_uw", "rank_unweighted"]],
                          on="pos", how="left")
    except Exception as e:
        print(f"  Warning: could not load unweighted scores ({e})")
        merged = wt.copy()
        merged["groupsim_uw"]     = None
        merged["rank_unweighted"] = None

    merged = merged.sort_values("rank_weighted").head(top_n)
    return merged[["rank_weighted", "pos", "groupsim_weighted", "z_score",
                   "sig_z2", "groupsim_uw", "rank_unweighted"]]

# ── Pretty-print top hits ─────────────────────────────────────────────────────

def print_top(tsv_df, scores_n, cols, group_idx, w, ids_ordered,
              top_n=10, with_gaps=False):
    label = "weighted GroupSim — gap as 21st char" if with_gaps else "weighted GroupSim (gap 0.85)"
    valid = tsv_df.dropna(subset=["groupsim_weighted"]) \
                  .sort_values("groupsim_weighted", ascending=False)
    print(f"\n{'─'*80}")
    print(f"  Top {top_n} positions — {label}")
    print(f"{'─'*80}")
    print(f"  {'Rank':>4}  {'Pos':>4}  {'Score':>7}  {'Z':>6}  {'CENPA (wt %)':40}  H3 (wt %)")
    print(f"{'─'*80}")
    for rank, (_, row) in enumerate(valid.head(top_n).iterrows(), 1):
        c = int(row["pos"])
        col = cols[c]
        detail = format_col_weighted(col, group_idx, w, ["CENPA", "H3"],
                                     with_gaps=with_gaps)
        parts = detail.split("  |  ")
        cenpa_str = parts[0].replace("CENPA: ", "")[:39]
        h3_str    = parts[1].replace("H3: ",    "") if len(parts) > 1 else ""
        z_str     = f"{row['z_score']:6.2f}" if pd.notna(row['z_score']) else "    NA"
        print(f"  {rank:>4}  {c:>4}  {row['groupsim_weighted']:7.3f}  "
              f"{z_str}  {cenpa_str:<40}  {h3_str}")


if __name__ == "__main__":
    # ── Standard runs (gap excluded from scoring) ─────────────────────────────
    tsv_df, scores_n, zs_map, cols, group_idx, w, ids_ordered = run(gap=0.85)
    print_top(tsv_df, scores_n, cols, group_idx, w, ids_ordered, top_n=10)

    run(gap=0.80)   # for Manhattan plot alignment with unweighted gap08 data

    # ── Gap-21 runs (gap treated as 21st character) ───────────────────────────
    tsv_g21, scores_g21, zs_g21, _, _, _, _ = run(gap=0.85, with_gaps=True)
    print_top(tsv_g21, scores_g21, cols, group_idx, w, ids_ordered,
              top_n=10, with_gaps=True)

    run(gap=0.80, with_gaps=True)   # for Manhattan plot

    # Comparison table (standard)
    comp = make_comparison(tsv_df, top_n=15)
    comp_path = OUT_DIR / "groupsim_weighted_comparison.tsv"
    comp.to_csv(comp_path, sep="\t", index=False, float_format="%.4f")
    print(f"\n  Comparison saved: {comp_path.name}")
    print("\n  Rank comparison (weighted rank | pos | w-score | uw-score | uw-rank):")
    print(f"  {'Wt-rank':>7}  {'Pos':>4}  {'W-score':>7}  {'UW-score':>8}  {'UW-rank':>7}  Δrank")
    print(f"  {'─'*60}")
    for _, row in comp.iterrows():
        uw_rank = int(row["rank_unweighted"]) if pd.notna(row["rank_unweighted"]) else -1
        delta   = uw_rank - int(row["rank_weighted"]) if uw_rank > 0 else 0
        uw_str  = f"{row['groupsim_uw']:8.3f}" if pd.notna(row["groupsim_uw"]) else "       NA"
        ur_str  = f"{uw_rank:>7}" if uw_rank > 0 else "      NA"
        d_str   = f"{delta:+d}" if uw_rank > 0 else "  NA"
        print(f"  {int(row['rank_weighted']):>7}  {int(row['pos']):>4}  "
              f"{row['groupsim_weighted']:7.3f}  {uw_str}  {ur_str}  {d_str}")

    print(f"\nDone. Outputs: {OUT_DIR}")
