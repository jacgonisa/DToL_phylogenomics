#!/usr/bin/env python3
"""
Clean GroupSim Manhattan plot — Satellite vs Transposon CENP-A,
unweighted vs HH-weighted (gap threshold 0.85).

Same layout as groupsim_manhattan_v2:
  ① CENP-A helix annotation track
  ② Main bar panel — gray = unweighted, teal = weighted, red = sig (weighted z≥2)
  ③ Z-score colour track

Key result: no positions survive z≥2 after HH reweighting, showing the
unweighted signal was driven by phylogenetic bias (plants vs animals).

Outputs: figures/groupsim_st_manhattan.{pdf,png}
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

BASE    = Path(__file__).parent
OUT_DIR = BASE / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

UW_TXT  = BASE / "split_entropy" / "groupsim_sat_trans"          / "groupsim_sat_trans_gap085.txt"
WT_TSV  = BASE / "split_entropy" / "groupsim_sat_trans_weighted" / "groupsim_st_weighted_gap085.tsv"
HELIX   = BASE / "split_entropy" / "groupsim_sat_trans"          / "helix_positions_gap085.txt"

C_UW    = "#b0bec5"   # unweighted — blue-gray
C_WT    = "#00897b"   # weighted   — teal
C_SIG   = "#e53935"   # significant weighted (z≥2) — none expected
C_HELIX = "#f48fb1"   # CENP-A helix bands — pink

# Annotate the top 6 weighted positions (none significant, so just show top scorers)
TOP_N_ANNOT = 6

# ── Read data ─────────────────────────────────────────────────────────────────

def read_uw(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip(): continue
        parts = line.split("\t", 2)
        try: pos = int(parts[0])
        except (ValueError, IndexError): continue
        scor = None if parts[1].strip() in ("None","","NA") else float(parts[1])
        rows.append({"pos": pos, "score_uw": scor})
    return pd.DataFrame(rows)

uw = read_uw(UW_TXT)
wt = pd.read_csv(WT_TSV, sep="\t").rename(columns={"groupsim_weighted": "score_wt"})

df = uw.merge(wt[["pos","score_wt","z_score"]], on="pos", how="outer").sort_values("pos")
df["sig"] = (df["z_score"] >= 2.0).fillna(False)
L = len(df)

# Read helix positions (start end, 1-based trimmed positions)
helix_ranges = []
for line in Path(HELIX).read_text().splitlines():
    if line.startswith("#") or not line.strip(): continue
    parts = line.split()
    if len(parts) >= 2:
        helix_ranges.append((int(parts[0]), int(parts[1])))

helix_labels = ["αN", "α1", "α2", "α3"]

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 6))
gs  = fig.add_gridspec(3, 1, height_ratios=[0.7, 4, 0.6],
                        hspace=0.08, left=0.07, right=0.93,
                        top=0.88, bottom=0.10)

ax_helix = fig.add_subplot(gs[0])
ax_main  = fig.add_subplot(gs[1], sharex=ax_helix)
ax_z     = fig.add_subplot(gs[2], sharex=ax_helix)

# ── ① Helix track ─────────────────────────────────────────────────────────────
ax_helix.set_xlim(0, L + 1)
ax_helix.set_ylim(0, 1)
ax_helix.axis("off")

for i, (s, e) in enumerate(helix_ranges):
    ax_helix.add_patch(mpatches.FancyBboxPatch(
        (s - 0.4, 0.08), e - s + 0.8, 0.84,
        boxstyle="round,pad=0,rounding_size=0.5",
        facecolor=C_HELIX, edgecolor="none", alpha=0.85, zorder=2
    ))
    lbl = helix_labels[i] if i < len(helix_labels) else f"α{i+1}"
    ax_helix.text((s + e) / 2, 0.50, lbl,
                  ha="center", va="center", fontsize=7.5,
                  fontweight="bold", color="#880e4f", zorder=3)

leg_h = [mpatches.Rectangle((0,0), 1, 1, facecolor=C_HELIX, edgecolor="none",
                              alpha=0.85, label="CENP-A α-helices (STRIDE)")]
ax_helix.legend(handles=leg_h, loc="upper right", fontsize=7.5,
                framealpha=0.9, borderpad=0.4, handlelength=1.0,
                edgecolor="#cccccc")

# ── ② Main panel ─────────────────────────────────────────────────────────────
bw = 0.38
ax_main.set_ylim(0, 1.22)
ax_main.set_ylabel("GroupSim score (normalised)", fontsize=9)
ax_main.axhline(0, color="#cccccc", lw=0.5, zorder=0)
ax_main.spines[["top","right","bottom"]].set_visible(False)
ax_main.tick_params(axis="x", bottom=False, labelbottom=False)
ax_main.tick_params(axis="y", labelsize=8)
ax_main.yaxis.grid(True, color="#eeeeee", zorder=0)

# unweighted bars
ax_main.bar(df["pos"] - bw * 0.5, df["score_uw"].fillna(0),
            width=bw, color=C_UW, alpha=0.8, zorder=2,
            label="Unweighted GroupSim")

# weighted bars
bar_colors = [C_SIG if s else C_WT for s in df["sig"]]
ax_main.bar(df["pos"] + bw * 0.5, df["score_wt"].fillna(0),
            width=bw, color=bar_colors, alpha=0.9, zorder=3,
            label="HH-weighted GroupSim")

# z=2 threshold line on score axis
valid_wt = df.dropna(subset=["score_wt"])
if valid_wt["sig"].any():
    z2_score = valid_wt.loc[valid_wt["sig"], "score_wt"].min()
    ax_main.axhline(z2_score, color=C_SIG, lw=0.8, ls="--", alpha=0.5,
                    zorder=1, label="z = 2 threshold (weighted)")
else:
    # draw the z=2 line anyway based on score corresponding to z=2
    mu  = valid_wt["score_wt"].mean()
    sd  = valid_wt["score_wt"].std()
    z2s = mu + 2 * sd
    ax_main.axhline(min(z2s, 1.0), color="#78909c", lw=0.8, ls="--",
                    alpha=0.55, zorder=1, label="z = 2 (no positions exceed)")

# Annotate top N weighted positions
top_rows = valid_wt.nlargest(TOP_N_ANNOT, "score_wt")
y_levels = [1.08, 1.16, 1.08, 1.16, 1.08, 1.16]
for i, (_, row) in enumerate(top_rows.sort_values("pos").iterrows()):
    p = int(row["pos"])
    s = row["score_wt"]
    z = row["z_score"] if pd.notna(row["z_score"]) else 0
    ax_main.annotate(
        f"pos {p}\n(z={z:.2f})",
        xy=(p + bw * 0.5, min(s, 0.99)),
        xytext=(p, y_levels[i % len(y_levels)]),
        fontsize=6.2, color="#004d40", ha="center", va="bottom",
        clip_on=False,
        arrowprops=dict(arrowstyle="-", color="#888888", lw=0.7, shrinkA=0),
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                  edgecolor="#80cbc4", alpha=0.92, linewidth=0.6),
        zorder=5
    )

# Legend
handles = [
    mpatches.Patch(facecolor=C_UW,  edgecolor="none", alpha=0.8,
                   label="Unweighted GroupSim"),
    mpatches.Patch(facecolor=C_WT,  edgecolor="none",
                   label="HH-weighted GroupSim"),
    mpatches.Patch(facecolor=C_SIG, edgecolor="none",
                   label="Significant (z ≥ 2, weighted) — none"),
]
ax_main.legend(handles=handles, loc="upper left", fontsize=7.5,
               framealpha=0.9, edgecolor="#cccccc", borderpad=0.4,
               handlelength=1.0, ncol=3)

# ── ③ Z-score track (weighted) ────────────────────────────────────────────────
ax_z.set_ylim(0, 1)
ax_z.set_xlim(0, L + 1)
ax_z.spines[["top","right","left"]].set_visible(False)
ax_z.tick_params(axis="y", left=False, labelleft=False)
ax_z.set_xlabel("Alignment position (trimmed, gap ≤ 85%)", fontsize=9)

cmap = plt.cm.YlOrRd
vmin, vmax = 0.0, 4.5
norm_c = Normalize(vmin=vmin, vmax=vmax)

for _, row in df.iterrows():
    z = row["z_score"] if pd.notna(row["z_score"]) else 0.0
    ax_z.bar(row["pos"], 1.0, width=0.9,
             color=cmap(norm_c(max(z, 0))), zorder=2)

sm = ScalarMappable(cmap=cmap, norm=norm_c)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax_z, orientation="vertical",
                    fraction=0.015, pad=0.01, aspect=6)
cbar.set_label("Z-score\n(weighted)", fontsize=7.5)
cbar.ax.tick_params(labelsize=7)
cbar.set_ticks([0, 2, 4])

ax_z.tick_params(axis="x", labelsize=8)

# ── Title ─────────────────────────────────────────────────────────────────────
n_uw_sig = int((df["score_uw"].fillna(0) > 0).sum())   # just for info
fig.suptitle(
    "GroupSim — Satellite CENP-A (n=198) vs Transposon CENP-A (n=129)  ·  gap threshold 85%\n"
    "After HH reweighting: 0 positions exceed z = 2  (unweighted had 4)",
    fontsize=10, fontweight="bold", y=0.97
)

for ext in ("pdf", "png"):
    p = OUT_DIR / f"groupsim_st_manhattan.{ext}"
    fig.savefig(p, dpi=300 if ext == "png" else None,
                bbox_inches="tight", facecolor="white")
    print(f"Saved: {p}")
plt.close(fig)
