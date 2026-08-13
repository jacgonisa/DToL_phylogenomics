#!/usr/bin/env python3
"""
Clean GroupSim Manhattan plot — unweighted vs HH-weighted (gap as 21st character),
gap threshold 0.80.

Layout (top to bottom):
  ① Helix annotation track  — H3 helices (blue), CENPA-specific regions (pink)
  ② Main panel              — gray bars = unweighted; blue/red bars = gap-21 weighted
                             top-scoring positions annotated with substitution summary
  ③ Z-score track           — horizontal bar coloured by z-score magnitude

Gap-21 mode: gaps treated as 21st character in weighted GroupSim scoring, so
group-specific insertion/deletion patterns (e.g. shorter αN in CENPA, longer L1 loop)
contribute alongside amino acid composition differences.

Outputs: figures/groupsim_manhattan_v2.{pdf,png}
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy.stats import zscore as sp_zscore

BASE    = Path(__file__).parent
OUT_DIR = BASE / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Paths ─────────────────────────────────────────────────────────────────────
UW_TXT  = BASE / "split_entropy" / "groupsim"          / "groupsim_gap08.txt"
WT_TSV  = BASE / "split_entropy" / "groupsim_weighted" / "groupsim_weighted_gap08_gap21.tsv"
HELIX   = BASE / "split_entropy" / "helix_positions_gap08.tsv"

# ── Colours ───────────────────────────────────────────────────────────────────
C_UW    = "#b0bec5"   # unweighted bars  — light blue-gray
C_WT    = "#1565c0"   # weighted bars    — deep blue
C_SIG   = "#e53935"   # significant (z≥2) weighted bar
C_H3    = "#90caf9"   # H3 helix band
C_CENPA = "#f48fb1"   # CENPA-specific helix band

# Top positions to annotate (trimmed positions for gap=0.80)
# * = new significant in gap-21 mode (include gap patterns)
TOP_ANNOTATIONS = {
    91: "Phe→Trp\n(pos 91)",
    74: "Gln→diverse\n(pos 74)",
    69: "Arg→Pro\n(pos 69)",
   109: "Gly→diverse\n(L1 loop, 109)",
   126: "Ile→Leu\n(pos 126)",
   114: "Thr→Ala/Ser\n(pos 114)",
    81: "α1 + gap↑\n(pos 81)*",
    96: "α2 + gap↑\n(pos 96)*",
    85: "α1/L1 + gap↑\n(pos 85)*",
}

# ── Read unweighted scores ────────────────────────────────────────────────────
def read_uw(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t", 2)
        try:
            pos  = int(parts[0])
        except (ValueError, IndexError):
            continue   # skip summary / header lines injected by groupsim.py
        scor = None if parts[1].strip() in ("None", "", "NA") else float(parts[1])
        rows.append({"pos": pos, "score_uw": scor})
    return pd.DataFrame(rows)

# ── Load data ─────────────────────────────────────────────────────────────────
uw = read_uw(UW_TXT)
wt = pd.read_csv(WT_TSV, sep="\t")
wt = wt.rename(columns={"groupsim_weighted": "score_wt"})

df = uw.merge(wt[["pos", "score_wt", "z_score"]], on="pos", how="outer").sort_values("pos")
df["sig"] = (df["z_score"] >= 2.0).fillna(False)

helices = pd.read_csv(HELIX, sep="\t")

L = len(df)
positions = df["pos"].values

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 6))
gs  = fig.add_gridspec(3, 1, height_ratios=[0.7, 4, 0.6],
                        hspace=0.08, left=0.07, right=0.93,
                        top=0.88, bottom=0.1)

ax_helix = fig.add_subplot(gs[0])
ax_main  = fig.add_subplot(gs[1], sharex=ax_helix)
ax_z     = fig.add_subplot(gs[2], sharex=ax_helix)

# ── ① Helix track ─────────────────────────────────────────────────────────────
ax_helix.set_xlim(-1, L)
ax_helix.set_ylim(0, 1)
ax_helix.axis("off")

for _, row in helices.iterrows():
    color = C_H3 if row["histone"] == "H3" else C_CENPA
    y0    = 0.55 if row["histone"] == "H3" else 0.05
    h     = 0.35
    ax_helix.add_patch(mpatches.FancyBboxPatch(
        (row["start"] - 0.4, y0), row["end"] - row["start"] + 0.8, h,
        boxstyle="round,pad=0,rounding_size=0.5",
        facecolor=color, edgecolor="none", alpha=0.85, zorder=2
    ))

# legend
leg_handles = [
    mpatches.Rectangle((0,0), 1, 1, facecolor=C_H3,    edgecolor="none", label="H3 α-helices"),
    mpatches.Rectangle((0,0), 1, 1, facecolor=C_CENPA, edgecolor="none", label="CENP-A α-helices"),
]
ax_helix.legend(handles=leg_handles, loc="upper right", fontsize=7.5, framealpha=0.9,
                borderpad=0.4, handlelength=1.0, ncol=2, edgecolor="#cccccc")

# ── ② Main panel ─────────────────────────────────────────────────────────────
bw = 0.38   # bar half-width
ax_main.set_ylim(0, 1.18)
ax_main.set_ylabel("GroupSim score (normalised)", fontsize=9)
ax_main.axhline(0, color="#cccccc", lw=0.5, zorder=0)
ax_main.spines[["top", "right", "bottom"]].set_visible(False)
ax_main.tick_params(axis="x", bottom=False, labelbottom=False)
ax_main.tick_params(axis="y", labelsize=8)
ax_main.yaxis.grid(True, color="#eeeeee", zorder=0)

# unweighted bars (background)
ax_main.bar(df["pos"] - bw * 0.5, df["score_uw"].fillna(0),
            width=bw, color=C_UW, alpha=0.8, zorder=2, label="Unweighted")

# weighted bars (foreground) — red if significant
bar_colors = [C_SIG if s else C_WT for s in df["sig"]]
ax_main.bar(df["pos"] + bw * 0.5, df["score_wt"].fillna(0),
            width=bw, color=bar_colors, alpha=0.9, zorder=3, label="HH-weighted")

# z=2 significance line on score axis (approximate)
valid_wt = df.dropna(subset=["score_wt"])
if len(valid_wt) > 1:
    z2_threshold = valid_wt["score_wt"][valid_wt["sig"]].min() if valid_wt["sig"].any() else None
    if z2_threshold is not None:
        ax_main.axhline(z2_threshold, color=C_SIG, lw=0.8, ls="--", alpha=0.6, zorder=1,
                        label=f"z = 2 threshold")

# Annotations for top positions — fixed y-positions staggered above bar
# Sort significant TOP_ANNOTATIONS positions by x to assign y slots
sig_annot_pos = sorted(
    [p for p in TOP_ANNOTATIONS if df[df["pos"] == p]["sig"].any()],
    key=lambda p: df[df["pos"] == p]["score_wt"].values[0], reverse=True
)
y_levels = [1.03, 1.10, 1.03, 1.10, 1.03, 1.10, 1.03, 1.10]
for rank, p in enumerate(sig_annot_pos):
    row  = df[df["pos"] == p].iloc[0]
    s    = row["score_wt"]
    y_txt = y_levels[rank % len(y_levels)]
    ax_main.annotate(
        TOP_ANNOTATIONS[p],
        xy=(p + bw * 0.5, min(s, 0.99)),
        xytext=(p, y_txt),
        fontsize=6.5, color="#222222", ha="center", va="bottom",
        clip_on=False,
        arrowprops=dict(arrowstyle="-", color="#888888", lw=0.7, shrinkA=0),
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                  edgecolor="#cccccc", alpha=0.92, linewidth=0.6),
        zorder=5
    )

# Legend
handles = [
    mpatches.Patch(facecolor=C_UW,  edgecolor="none", alpha=0.8, label="Unweighted GroupSim"),
    mpatches.Patch(facecolor=C_WT,  edgecolor="none",             label="HH-weighted GroupSim (gap = 21st char)"),
    mpatches.Patch(facecolor=C_SIG, edgecolor="none",             label="Significant (z ≥ 2, weighted)"),
]
ax_main.legend(handles=handles, loc="upper left", fontsize=7.5,
               framealpha=0.9, edgecolor="#cccccc", borderpad=0.4,
               handlelength=1.0, ncol=3)

# ── ③ Z-score track ───────────────────────────────────────────────────────────
ax_z.set_ylim(0, 1)
ax_z.set_xlim(-1, L)
ax_z.spines[["top", "right", "left"]].set_visible(False)
ax_z.tick_params(axis="y", left=False, labelleft=False)
ax_z.set_xlabel("Alignment position (trimmed, gap ≤ 80%)", fontsize=9)

cmap  = plt.cm.YlOrRd
vmin, vmax = 0.0, 4.5
norm  = Normalize(vmin=vmin, vmax=vmax)

for _, row in df.iterrows():
    z = row["z_score"] if pd.notna(row["z_score"]) else 0.0
    ax_z.bar(row["pos"], 1.0, width=0.9, color=cmap(norm(max(z, 0))),
             zorder=2)

# colorbar
sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax_z, orientation="vertical",
                    fraction=0.015, pad=0.01, aspect=6)
cbar.set_label("Z-score\n(weighted)", fontsize=7.5)
cbar.ax.tick_params(labelsize=7)
cbar.set_ticks([0, 2, 4])

ax_z.tick_params(axis="x", labelsize=8)

# ── Title & save ──────────────────────────────────────────────────────────────
fig.suptitle(
    "GroupSim — CENP-A (n=418) vs H3-like (n=901)  ·  gap threshold 80%\n"
    "HH-weighted scores (gap as 21st char): 9 significant positions — "
    "6 AA substitutions + 3 involving group-specific indels (marked *)",
    fontsize=10, fontweight="bold", y=0.97
)

for ext in ("pdf", "png"):
    p = OUT_DIR / f"groupsim_manhattan_v2.{ext}"
    fig.savefig(p, dpi=300 if ext == "png" else None,
                bbox_inches="tight", facecolor="white")
    print(f"Saved: {p}")
plt.close(fig)
