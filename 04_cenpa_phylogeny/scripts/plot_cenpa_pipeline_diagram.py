#!/usr/bin/env python3
"""
CENP-A identification and phylogeny pipeline diagram.

Main vertical flow (left column) + miniprot rescue branch (right column).
Saves: figures/cenpa_pipeline_diagram.{pdf, png, svg}
"""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT_DIR = Path(__file__).parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Colours ───────────────────────────────────────────────────────────────────
C_INPUT    = "#e8edf2"   # light slate  — input data
C_PROCESS  = "#d0e4f7"   # light blue   — computational step
C_CURATE   = "#fde8cc"   # light amber  — curation / filtering
C_RESCUE   = "#d6efd4"   # light green  — miniprot rescue branch
C_MERGE    = "#ede0f5"   # light purple — merge step
C_OUTPUT   = "#c8e6c9"   # mid green    — final output
C_EDGE     = "#555555"   # arrow / border colour
C_ANNOT    = "#777777"   # annotation text

# ── Canvas ────────────────────────────────────────────────────────────────────
fig_w, fig_h = 13, 18
fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.set_xlim(0, fig_w)
ax.set_ylim(0, fig_h)
ax.axis("off")

# ── Box drawing helpers ───────────────────────────────────────────────────────
BOX_W_MAIN   = 4.8
BOX_W_BRANCH = 3.6
BOX_H        = 0.78
RADIUS       = 0.12

def box(ax, cx, cy, w, h, color, label, sublabel=None, bold=False):
    """Rounded rectangle centred at (cx, cy)."""
    rect = FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle=f"round,pad=0,rounding_size={RADIUS}",
        linewidth=1.1, edgecolor=C_EDGE,
        facecolor=color, zorder=3
    )
    ax.add_patch(rect)
    fw = "bold" if bold else "normal"
    if sublabel:
        ax.text(cx, cy + 0.13, label, ha="center", va="center",
                fontsize=9.5, fontweight=fw, color="#222222", zorder=4)
        ax.text(cx, cy - 0.17, sublabel, ha="center", va="center",
                fontsize=7.8, color="#555555", zorder=4, style="italic")
    else:
        ax.text(cx, cy, label, ha="center", va="center",
                fontsize=9.5, fontweight=fw, color="#222222", zorder=4)
    return cy - h / 2   # bottom y

def arrow(ax, x, y_top, y_bot, color=C_EDGE, lw=1.5):
    """Straight downward arrow."""
    ax.annotate("", xy=(x, y_bot + 0.04), xytext=(x, y_top - 0.04),
                arrowprops=dict(arrowstyle="-|>", color=color,
                                lw=lw, mutation_scale=14),
                zorder=2)

def h_arrow(ax, x_from, x_to, y, color=C_EDGE, lw=1.5):
    """Horizontal arrow."""
    ax.annotate("", xy=(x_to, y), xytext=(x_from, y),
                arrowprops=dict(arrowstyle="-|>", color=color,
                                lw=lw, mutation_scale=14),
                zorder=2)

def annot(ax, cx, cy, text, color=C_ANNOT, size=7.5, ha="left"):
    ax.text(cx, cy, text, ha=ha, va="center",
            fontsize=size, color=color, zorder=4,
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                      edgecolor="none", alpha=0.7))

# ── Column x positions ────────────────────────────────────────────────────────
X_MAIN   = 4.2    # main pipeline column
X_BRANCH = 9.8    # rescue branch column
X_ANNOT  = 6.85   # annotation x for main column

# ── Y positions (top → bottom) ───────────────────────────────────────────────
Y = {
    "proteomes":  17.3,
    "dmnd_db":    16.1,
    "dmnd_srch":  14.7,
    "candidates": 13.4,
    "brh":        12.2,
    "merge":      10.85,
    "mafft1":      9.55,
    "curate":      8.3,
    "cenpa430":    7.05,
    "add_h3":      5.85,
    "mafft2":      4.65,
    "iqtree":      3.35,
    "final":       2.05,
}

# Rescue branch y positions (aligned between candidates and merge)
Y_R = {
    "absent":   13.4,
    "miniprot": 12.2,
    "cluster":  11.05,
    "hmm":       9.9,
}

# ── MAIN PIPELINE ─────────────────────────────────────────────────────────────

# 1 — Proteomes
box(ax, X_MAIN, Y["proteomes"], BOX_W_MAIN, BOX_H, C_INPUT,
    "325 species proteomes",
    "8,975,424 proteins", bold=True)

arrow(ax, X_MAIN, Y["proteomes"] - BOX_H/2, Y["dmnd_db"] + BOX_H/2)

# 2 — DIAMOND database
box(ax, X_MAIN, Y["dmnd_db"], BOX_W_MAIN, BOX_H, C_PROCESS,
    "DIAMOND protein database")

arrow(ax, X_MAIN, Y["dmnd_db"] - BOX_H/2, Y["dmnd_srch"] + BOX_H/2)

# 3 — DIAMOND search
box(ax, X_MAIN, Y["dmnd_srch"], BOX_W_MAIN, BOX_H * 1.25, C_PROCESS,
    "DIAMOND blastp",
    "ultra-sensitive | e ≤ 1×10⁻⁵ | qcov/scov ≥ 25%")
annot(ax, X_ANNOT, Y["dmnd_srch"],
      "query: 450 external\nCENP-A references")

arrow(ax, X_MAIN, Y["dmnd_srch"] - BOX_H*1.25/2, Y["candidates"] + BOX_H/2)

# 4 — Candidates
box(ax, X_MAIN, Y["candidates"], BOX_W_MAIN, BOX_H, C_PROCESS,
    "CENP-A candidate hits")

# Branch arrow (right) from candidates to rescue
mid_branch_y = Y["candidates"]
h_arrow(ax, X_MAIN + BOX_W_MAIN/2, X_BRANCH - BOX_W_BRANCH/2, mid_branch_y,
        color="#4a9e6b", lw=1.5)
ax.text((X_MAIN + BOX_W_MAIN/2 + X_BRANCH - BOX_W_BRANCH/2)/2,
        mid_branch_y + 0.22,
        "76 species absent", ha="center", va="bottom",
        fontsize=7.8, color="#4a9e6b", fontstyle="italic")

arrow(ax, X_MAIN, Y["candidates"] - BOX_H/2, Y["brh"] + BOX_H/2)

# 5 — BRH filter
box(ax, X_MAIN, Y["brh"], BOX_W_MAIN, BOX_H, C_CURATE,
    "Reciprocal best HMM (BRH)",
    "size filter: 50–300 aa | per-species best hit")

arrow(ax, X_MAIN, Y["brh"] - BOX_H/2, Y["merge"] + BOX_H/2)

# 6 — Merge
box(ax, X_MAIN, Y["merge"], BOX_W_MAIN, BOX_H, C_MERGE,
    "Merged CENP-A candidates", bold=True)

# Rescue rejoins here (arrow from branch bottom → merge box right side)
ax.annotate("",
    xy  =(X_MAIN + BOX_W_MAIN/2, Y["merge"]),
    xytext=(X_BRANCH - BOX_W_BRANCH/2, Y_R["hmm"] - BOX_H/2 - 0.05),
    arrowprops=dict(arrowstyle="-|>", color="#4a9e6b",
                    lw=1.5, mutation_scale=14,
                    connectionstyle="arc3,rad=-0.25"),
    zorder=2)

arrow(ax, X_MAIN, Y["merge"] - BOX_H/2, Y["mafft1"] + BOX_H/2)

# 7 — MAFFT + FastTree (round 1)
box(ax, X_MAIN, Y["mafft1"], BOX_W_MAIN, BOX_H, C_PROCESS,
    "MAFFT alignment → FastTree",
    "initial tree for curation")

arrow(ax, X_MAIN, Y["mafft1"] - BOX_H/2, Y["curate"] + BOX_H/2)

# 8 — Curation
box(ax, X_MAIN, Y["curate"], BOX_W_MAIN, BOX_H * 1.3, C_CURATE,
    "Phylogenetic curation",
    "remove non-CENP-A clusters")
annot(ax, X_ANNOT, Y["curate"],
      "removed: Dictyostelium H3v1,\nSymbiodinium, Perkinsus,\n+ 2 manual outliers")

arrow(ax, X_MAIN, Y["curate"] - BOX_H*1.3/2, Y["cenpa430"] + BOX_H/2)

# 9 — Curated CENP-A
box(ax, X_MAIN, Y["cenpa430"], BOX_W_MAIN, BOX_H, C_OUTPUT,
    "Curated CENP-A set",
    "430 sequences", bold=True)

arrow(ax, X_MAIN, Y["cenpa430"] - BOX_H/2, Y["add_h3"] + BOX_H/2)

# 10 — Add H3 + archaea
box(ax, X_MAIN, Y["add_h3"], BOX_W_MAIN, BOX_H * 1.2, C_PROCESS,
    "Add outgroups",
    "H3: ≤3/species, 99% id clustering (MMseqs2) → 897 seqs  |  Archaea: 10 seqs")

arrow(ax, X_MAIN, Y["add_h3"] - BOX_H*1.2/2, Y["mafft2"] + BOX_H/2)

# 11 — MAFFT + ClipKIT
box(ax, X_MAIN, Y["mafft2"], BOX_W_MAIN, BOX_H * 1.2, C_PROCESS,
    "MAFFT alignment → ClipKIT trimming",
    "gappy mode | 1,681 / 3,569 sites retained")

arrow(ax, X_MAIN, Y["mafft2"] - BOX_H*1.2/2, Y["iqtree"] + BOX_H/2)

# 12 — IQ-TREE2
box(ax, X_MAIN, Y["iqtree"], BOX_W_MAIN, BOX_H * 1.2, C_PROCESS,
    "IQ-TREE2 maximum-likelihood",
    "Q.pfam+R7 | UFBoot 1000 | -bnni")

arrow(ax, X_MAIN, Y["iqtree"] - BOX_H*1.2/2, Y["final"] + BOX_H/2)

# 13 — Final tree
box(ax, X_MAIN, Y["final"], BOX_W_MAIN, BOX_H, C_OUTPUT,
    "CENP-A phylogeny",
    "1,329 tips | rooted on archaeal MRCA", bold=True)

# ── RESCUE BRANCH ─────────────────────────────────────────────────────────────

box(ax, X_BRANCH, Y_R["absent"], BOX_W_BRANCH, BOX_H, C_RESCUE,
    "76 species",
    "CENP-A absent from annotation")

arrow(ax, X_BRANCH, Y_R["absent"] - BOX_H/2, Y_R["miniprot"] + BOX_H/2,
      color="#4a9e6b")

box(ax, X_BRANCH, Y_R["miniprot"], BOX_W_BRANCH, BOX_H, C_RESCUE,
    "miniprot",
    "protein-to-genome mapping | -I")

arrow(ax, X_BRANCH, Y_R["miniprot"] - BOX_H/2, Y_R["cluster"] + BOX_H/2,
      color="#4a9e6b")

box(ax, X_BRANCH, Y_R["cluster"], BOX_W_BRANCH, BOX_H, C_RESCUE,
    "1 kb window clustering",
    "remove genomic duplicates")

arrow(ax, X_BRANCH, Y_R["cluster"] - BOX_H/2, Y_R["hmm"] + BOX_H/2,
      color="#4a9e6b")

box(ax, X_BRANCH, Y_R["hmm"], BOX_W_BRANCH, BOX_H, C_RESCUE,
    "CENP-A HMM validation",
    "98.8% pass rate  |  1,299 / 1,315 seqs")

# Branch label
ax.text(X_BRANCH, Y_R["absent"] + BOX_H/2 + 0.28,
        "Miniprot rescue", ha="center", va="bottom",
        fontsize=10, fontweight="bold", color="#4a9e6b")
ax.add_patch(FancyBboxPatch(
    (X_BRANCH - BOX_W_BRANCH/2 - 0.15, Y_R["hmm"] - BOX_H/2 - 0.18),
    BOX_W_BRANCH + 0.3,
    Y_R["absent"] - Y_R["hmm"] + BOX_H + 0.55,
    boxstyle=f"round,pad=0,rounding_size={RADIUS}",
    linewidth=1.2, edgecolor="#4a9e6b",
    facecolor="none", linestyle="--", zorder=1
))

# ── TITLE & LEGEND ────────────────────────────────────────────────────────────
ax.text(fig_w / 2, fig_h - 0.28, "CENP-A identification and phylogeny pipeline",
        ha="center", va="top", fontsize=13, fontweight="bold", color="#1a1a1a")

legend_items = [
    mpatches.Patch(facecolor=C_INPUT,   edgecolor=C_EDGE, label="Input data"),
    mpatches.Patch(facecolor=C_PROCESS, edgecolor=C_EDGE, label="Computational step"),
    mpatches.Patch(facecolor=C_CURATE,  edgecolor=C_EDGE, label="Curation / filtering"),
    mpatches.Patch(facecolor=C_RESCUE,  edgecolor="#4a9e6b", label="Miniprot rescue"),
    mpatches.Patch(facecolor=C_MERGE,   edgecolor=C_EDGE, label="Merge"),
    mpatches.Patch(facecolor=C_OUTPUT,  edgecolor=C_EDGE, label="Key output"),
]
ax.legend(handles=legend_items, loc="lower right", fontsize=8,
          framealpha=0.9, edgecolor="#cccccc",
          handlelength=1.2, handleheight=1.0,
          bbox_to_anchor=(0.99, 0.01))

# ── SAVE ──────────────────────────────────────────────────────────────────────
for ext in ("pdf", "png", "svg"):
    out = OUT_DIR / f"cenpa_pipeline_diagram.{ext}"
    dpi = 300 if ext == "png" else None
    fig.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
    print(f"Saved: {out}")

plt.close(fig)
