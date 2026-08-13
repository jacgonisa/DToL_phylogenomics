#!/usr/bin/env python3
"""
GroupSim: CENH3 from Satellite species vs CENH3 from Transposon species.

Groups defined by centromere type in iTOL_bnni/iTOL_centromere_symbols.txt:
  Satellite  = #ff006e  (n=198 CenpA tips)
  Transposon = #3a86ff  (n=129 CenpA tips)

Alignment trimmed independently on just the Sat+Trans CenpA sequences.
Helix bands (STRIDE-derived) added to the output plot.

Run from:
  cd DToL_phylogenomics_publication_325genomes/04_cenpa_phylogeny/
  python3 run_groupsim_sat_trans_325sp.py
"""

import subprocess, sys, math
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker

BASE     = Path(__file__).parent
GROUPSIM = Path("/home/jg2070/Desktop/groupsim-py3/src/groupsim.py")
ENTDIR   = Path("/home/jg2070/Desktop/dtol_review_August/PhylogeneticProfiling/"
                "20_cenpa_analysis/cenpa_plus_H3_entropy")

ALN          = BASE / "cenpa430_H3_archaea10.aligned.clipkit.325sp.fasta"
H3_ALL_FASTA = BASE / "H3_all.aligned.fasta"
ITOL_SYMBOLS = BASE / "iTOL_bnni" / "iTOL_centromere_symbols.txt"
STRIDE_CENPA = ENTDIR / "AF-Q8RVQ9-F1.stride"

OUT_DIR  = BASE / "split_entropy" / "groupsim_sat_trans"
OUT_DIR.mkdir(parents=True, exist_ok=True)

ARCHAEA_IDS = {
    "OLS22332.1","OLS24873.1","OLS21974.1","KKK41979.1","KXH71038.1",
    "OLS18261.1","OLS16336.1","BAD86478.1","OIO61677.1","OIO41945.1",
}
GAP = {"-"}
MAPPED_CENPA_ID = "ddAraThal4_chr_1_000021.1"

# ── helpers ───────────────────────────────────────────────────────────────────
def read_fasta(path):
    seqs, cur, buf = {}, None, []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            if cur is not None: seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]; buf = []
        else: buf.append(line)
    if cur is not None: seqs[cur] = "".join(buf)
    return seqs

def trim_columns(seqs, thr):
    ids = list(seqs.keys())
    arr = [seqs[i] for i in ids]
    n, L = len(arr), len(arr[0])
    keep = [i for i in range(L)
            if sum(1 for s in arr if s[i] in GAP) / n <= thr]
    trimmed = {sid: "".join(arr[k][c] for c in keep) for k, sid in enumerate(ids)}
    return trimmed, keep

def write_fasta(seqs, path):
    with open(path, "w") as fh:
        for sid, seq in seqs.items():
            fh.write(f">{sid}\n{seq}\n")

def parse_stride_helices(path):
    ranges = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("LOC  AlphaHelix"):
            p = line.split()
            ranges.append((int(p[3]), int(p[6])))
    return ranges

def seqpos_to_aln(aln_seq):
    m, pos = {}, 0
    for i, aa in enumerate(aln_seq):
        if aa != "-":
            pos += 1
            m[pos] = i
    return m

def project_helices(aln_ranges, keep):
    keep_set = set(keep)
    keep_pos = {orig: i + 1 for i, orig in enumerate(keep)}
    out = []
    for s, e in aln_ranges:
        inside = [k for k in keep if s <= k <= e]
        if inside:
            out.append((keep_pos[inside[0]], keep_pos[inside[-1]]))
    return out

def read_groupsim_scores(path):
    scores = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"): continue
        parts = line.split("\t")
        if len(parts) < 2: continue
        val = parts[1]
        scores.append(np.nan if val.strip() == "None" else float(val))
    return np.array(scores, dtype=float)

# ── parse iTOL symbols → Satellite / Transposon tip IDs ──────────────────────
sat_tips, trans_tips = [], []
for line in ITOL_SYMBOLS.read_text().splitlines():
    if line.startswith("#") or not line.strip() or "\t" not in line: continue
    parts = line.split("\t")
    if len(parts) < 4: continue
    tip, colour = parts[0], parts[3]
    if colour == "#ff006e": sat_tips.append(tip)
    elif colour == "#3a86ff": trans_tips.append(tip)

print(f"Satellite CenpA tips:   {len(sat_tips)}")
print(f"Transposon CenpA tips:  {len(trans_tips)}")

# ── load alignment, keep only Sat+Trans CenpA ────────────────────────────────
seqs_all = read_fasta(ALN)
keep_ids = set(sat_tips) | set(trans_tips)
seqs_st  = {k: v for k, v in seqs_all.items() if k in keep_ids}
print(f"Sequences retained: {len(seqs_st)}")

# ── STRIDE helix mapping on CenpA reference ───────────────────────────────────
cenpa_helix_raw = parse_stride_helices(STRIDE_CENPA)
aln_ref         = seqs_all[MAPPED_CENPA_ID]
map_cenpa       = seqpos_to_aln(aln_ref)
helix_aln       = [(map_cenpa[s], map_cenpa[e])
                   for s, e in cenpa_helix_raw if s in map_cenpa and e in map_cenpa]

# ── publication-quality plot ─────────────────────────────────────────────────
import pandas as pd

# Helix labels: aN, a1, a2, a3 (N-terminal tail + 3 core helices of the
# histone fold; CENH3 has 4 predicted helices from STRIDE)
_HELIX_NAMES = ["aN", "a1", "a2", "a3", "a4", "a5"]

def _plot_groupsim(scores, pos, helix_trim, gap, sat_n, trans_n, out_dir, tag):
    """Manhattan-style GroupSim plot with Z-scores."""

    # ── Z-score transform (on valid positions only) ───────────────────────────
    valid    = ~np.isnan(scores)
    mu, sd   = np.nanmean(scores), np.nanstd(scores)
    z        = np.where(valid, (scores - mu) / (sd + 1e-10), np.nan)

    z_thresh = 2.0                               # Z ≥ 2 highlighted
    sig_mask = np.where(valid, z >= z_thresh, False)
    nan_mask = ~valid
    n_sig    = int(sig_mask.sum())
    n_nan    = int(nan_mask.sum())

    # ── Layout: two stacked panels ────────────────────────────────────────────
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(16, 6.5),
        gridspec_kw={"height_ratios": [1, 3.5], "hspace": 0.06},
        sharex=True
    )
    fig.patch.set_facecolor("white")

    # ── Top panel: helix bands ────────────────────────────────────────────────
    ax_top.set_facecolor("white")
    ax_top.set_ylim(0, 1)
    ax_top.set_xlim(pos[0] - 1, pos[-1] + 1)
    ax_top.axis("off")

    for i, (s, e) in enumerate(helix_trim):
        ax_top.add_patch(mpatches.FancyBboxPatch(
            (s - 0.4, 0.15), (e - s + 0.8), 0.70,
            boxstyle="round,pad=0.0",
            fc="#C62828", ec="#7B0000", lw=0.8, alpha=0.75
        ))
        mid = (s + e) / 2.0
        lbl = _HELIX_NAMES[i] if i < len(_HELIX_NAMES) else f"H{i+1}"
        ax_top.text(mid, 0.50, lbl, ha="center", va="center",
                    fontsize=8, fontweight="bold", color="white")

    ax_top.text(pos[0] - 1, 0.50, "CENH3\nalpha-\nhelices",
                ha="right", va="center", fontsize=7.5,
                color="#C62828", fontweight="bold")
    ax_top.set_title(
        f"GroupSim — CENH3: Satellite (n={sat_n}) vs Transposon (n={trans_n})"
        f"   |   gap threshold {gap*100:.0f}%",
        fontsize=13, fontweight="bold", pad=7, loc="left"
    )

    # ── Bottom panel: Manhattan Z-score plot ──────────────────────────────────
    ax_bot.set_facecolor("#fafafa")

    # Alternating light bands every 20 positions for readability
    for band_s in range(int(pos[0]), int(pos[-1]) + 1, 20):
        if (band_s // 20) % 2 == 0:
            ax_bot.axvspan(band_s, min(band_s + 20, pos[-1] + 1),
                           color="#f0f0f0", alpha=0.6, linewidth=0, zorder=0)

    # Threshold line
    ax_bot.axhline(z_thresh, color="#C62828", linewidth=1.1,
                   linestyle="--", alpha=0.75, zorder=2,
                   label=f"Z = {z_thresh:.0f} threshold")
    ax_bot.axhline(0, color="#888888", linewidth=0.6,
                   linestyle="-", alpha=0.5, zorder=2)

    # Masked positions: tiny grey triangles at y = 0
    if n_nan:
        ax_bot.scatter(pos[nan_mask], np.zeros(n_nan),
                       marker="v", s=18, color="#BBBBBB", alpha=0.7,
                       zorder=3, linewidths=0,
                       label=f"Masked  (n={n_nan})")

    # Below-threshold points
    below = valid & ~sig_mask
    ax_bot.scatter(pos[below], z[below],
                   s=28, color="#4A90D9", alpha=0.65,
                   linewidths=0, zorder=4,
                   label=f"Z < {z_thresh:.0f}  (n={int(below.sum())})")

    # Above-threshold points (significant)
    ax_bot.scatter(pos[sig_mask], z[sig_mask],
                   s=55, color="#d62728", alpha=0.90,
                   linewidths=0, zorder=5,
                   label=f"Z ≥ {z_thresh:.0f}  (n={n_sig})")

    # Annotate the top 5 positions
    top5_idx = np.argsort(np.where(valid, z, -np.inf))[-5:][::-1]
    for idx in top5_idx:
        ax_bot.annotate(
            str(int(pos[idx])),
            xy=(pos[idx], z[idx]),
            xytext=(0, 7), textcoords="offset points",
            ha="center", fontsize=7.5, color="#7B0000", fontweight="bold",
            arrowprops=dict(arrowstyle="-", color="#7B0000",
                            lw=0.8, shrinkB=3)
        )

    ax_bot.set_xlabel(
        f"Trimmed alignment position  (gap ≤ {gap*100:.0f}%)",
        fontsize=12
    )
    ax_bot.set_ylabel("GroupSim Z-score", fontsize=12)
    ax_bot.spines[["top", "right"]].set_visible(False)
    ax_bot.spines["left"].set_linewidth(0.8)
    ax_bot.spines["bottom"].set_linewidth(0.8)
    ax_bot.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(6))
    ax_bot.tick_params(axis="both", labelsize=10)

    handles, lbls = ax_bot.get_legend_handles_labels()
    ax_bot.legend(handles, lbls, fontsize=9, loc="upper left",
                  framealpha=0.92, edgecolor="grey", borderpad=0.8)

    for ext in [".png", ".pdf"]:
        out = out_dir / f"groupsim_sat_trans_gap{tag}{ext}"
        plt.savefig(out, dpi=200 if ext == ".png" else 100, bbox_inches="tight")
        print(f"  saved {out.name}")
    plt.close()

    # TSV for R
    df_out = pd.DataFrame({"pos": pos, "groupsim": scores, "z_score": z,
                           "sig_z2": sig_mask.astype(int)})
    df_out.to_csv(out_dir / f"groupsim_sat_trans_gap{tag}.tsv", sep="\t", index=False)


# ── run for each gap threshold ────────────────────────────────────────────────
for gap in [0.80, 0.85, 0.90]:
    tag = str(gap).replace(".", "")
    print(f"\n── Gap {gap} ──")

    trimmed, keep = trim_columns(seqs_st, gap)
    print(f"  columns kept: {len(keep)}")

    helix_trim = project_helices(helix_aln, keep)
    np.savetxt(OUT_DIR / f"helix_positions_gap{tag}.txt",
               np.array(helix_trim), fmt="%d", header="start end")

    # Write ordered FASTA: Satellite first, then Transposon
    ordered = {sid: trimmed[sid] for sid in sat_tips + trans_tips if sid in trimmed}
    aln_out = OUT_DIR / f"trimmed_sat_trans_gap{tag}.fasta"
    write_fasta(ordered, aln_out)

    # Group file
    sat_in  = [s for s in sat_tips   if s in trimmed]
    trans_in= [s for s in trans_tips if s in trimmed]
    grp_file = OUT_DIR / f"groups_sat_trans_gap{tag}.txt"
    with open(grp_file, "w") as fh:
        fh.write("Satellite:"  + ",".join(sat_in)   + "\n")
        fh.write("Transposon:" + ",".join(trans_in) + "\n")
    print(f"  Satellite: {len(sat_in)}  Transposon: {len(trans_in)}")

    # Run GroupSim (skip if output already exists)
    out_prefix = str(OUT_DIR / f"groupsim_sat_trans_gap{tag}")
    gs_out = Path(out_prefix + ".txt")
    if gs_out.exists():
        print(f"  GroupSim output exists, skipping run")
    else:
        cmd = [
            sys.executable, str(GROUPSIM),
            "-k", str(grp_file),
            "-o", out_prefix,
            "-c", str(gap),
            "-g", "0.5",
            "-w", "3", "-l", "0.7", "-n",
            str(aln_out)
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print("STDERR:", result.stderr[-1000:])
            raise RuntimeError(f"GroupSim failed (gap={gap})")
        print(f"  GroupSim done → {Path(out_prefix).name}.*")

    # ── Plot ──────────────────────────────────────────────────────────────────
    scores = read_groupsim_scores(out_prefix + ".txt")
    pos    = np.arange(1, len(scores) + 1)
    _plot_groupsim(scores, pos, helix_trim, gap, len(sat_in), len(trans_in),
                   OUT_DIR, tag)


print("\nDone. Outputs:", OUT_DIR)
