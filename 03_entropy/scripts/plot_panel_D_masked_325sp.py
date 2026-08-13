#!/usr/bin/env python3
"""
Panel D — split entropy, per-group masking, filled-area aesthetic,
with dotted vertical lines delimiting helix boundaries.

Output: 03_entropy/figures/panel_D_split_entropy_325sp_masked_fills.{pdf,png}
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BASE    = Path(__file__).parent
PUB_DIR = BASE.parent
SP_DIR  = PUB_DIR / "04_cenpa_phylogeny" / "split_entropy"
FIG_DIR = BASE / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── data ──────────────────────────────────────────────────────────────────────
df      = pd.read_csv(SP_DIR / "split_entropy_bnni_pergroup_gap085.tsv",
                      sep="\t", na_values=["nan", "NaN", "NA", ""])
helices = pd.read_csv(SP_DIR / "helix_positions_gap085.tsv", sep="\t")

pos   = df["pos"].values
cenpa = df["CENPA"].values.astype(float)
h3    = df["H3"].values.astype(float)

n_cenpa = 418
n_h3    = 901

C_CENPA = "#C62828"
C_H3    = "#1565C0"

# ── figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(13, 4.5))

# Filled areas — gaps appear automatically where values are NaN (masked)
ax.fill_between(pos, 0, h3,
                color=C_H3, alpha=0.18, linewidth=0, zorder=1)
ax.fill_between(pos, 0, cenpa,
                color=C_CENPA, alpha=0.22, linewidth=0, zorder=2)

# Lines on top
ax.plot(pos, h3,    color=C_H3,    linewidth=0.75, alpha=0.90,
        label=f"H3-like  (n={n_h3})", zorder=4)
ax.plot(pos, cenpa, color=C_CENPA, linewidth=0.85, alpha=0.95,
        label=f"CENP-A  (n={n_cenpa})", zorder=5)

# Light grey background shading for CENP-A helix regions
for _, row in helices[helices["histone"] == "CENPA"].iterrows():
    ax.axvspan(row["start"], row["end"], color="#bbbbbb", alpha=0.30,
               linewidth=0, zorder=0)

# ── Helix annotation strips + dotted boundary lines ───────────────────────────
seen_boundaries = set()

for _, row in helices.iterrows():
    s, e = row["start"], row["end"]
    is_cenpa = (row["histone"] == "CENPA")
    col  = C_CENPA if is_cenpa else C_H3
    ymin = 0.935 if is_cenpa else 0.875
    ymax = 0.990 if is_cenpa else 0.930

    # Solid helix band
    ax.axvspan(s, e, ymin=ymin, ymax=ymax,
               color=col, alpha=0.45, linewidth=0, zorder=3)

    # Dotted vertical lines at helix edges (drawn once per x position)
    for x in (s, e):
        if x not in seen_boundaries:
            ax.axvline(x, color="#555555", linewidth=0.65,
                       linestyle=":", alpha=0.60, zorder=3)
            seen_boundaries.add(x)

# ── Axes & decoration ─────────────────────────────────────────────────────────
ax.set_xlim(pos.min(), pos.max())
ax.set_ylim(0, 1.01)
ax.set_xlabel("Alignment position (per-group mask, gap ≤ 85%)", fontsize=11)
ax.set_ylabel("Normalised Shannon entropy", fontsize=11)
ax.set_yticks([0, 0.25, 0.50, 0.75, 1.00])
ax.spines[["top", "right"]].set_visible(False)
ax.yaxis.grid(True, color="#eeeeee", linewidth=0.5, zorder=0)

# Legend
patch_cenpa_h  = mpatches.Patch(color=C_CENPA,  alpha=0.45, label="CENP-A α-helices")
patch_h3_h     = mpatches.Patch(color=C_H3,     alpha=0.45, label="H3 α-helices")
patch_ss       = mpatches.Patch(color="#bbbbbb", alpha=0.55, label="CENP-A secondary structure")
handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles + [patch_cenpa_h, patch_h3_h, patch_ss],
    labels  + ["CENP-A α-helices", "H3 α-helices", "CENP-A secondary structure"],
    fontsize=9, loc="upper right", framealpha=0.90,
    edgecolor="#cccccc", handlelength=1.4
)

n_masked_c = int(np.isnan(cenpa).sum())
n_masked_h = int(np.isnan(h3).sum())

ax.set_title(
    "CENP-A vs H3 split entropy — bnni tree (325-sp)",
    fontsize=12, fontweight="bold"
)
ax.text(0.5, -0.13,
        f"Per-group masking: gaps shown where group gap > 85%  "
        f"(CENP-A: {n_masked_c} masked pos · H3: {n_masked_h} masked pos).  "
        "Dotted lines = helix boundaries.",
        ha="center", va="top", fontsize=8, color="#616161",
        transform=ax.transAxes)

plt.tight_layout(rect=[0, 0.04, 1, 1])

for ext in ("pdf", "png"):
    out = FIG_DIR / f"panel_D_split_entropy_325sp_masked_fills.{ext}"
    fig.savefig(str(out),
                dpi=300 if ext == "png" else None,
                bbox_inches="tight", facecolor="white")
    print(f"Saved: {out}")

plt.close(fig)
