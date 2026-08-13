#!/usr/bin/env python3
"""
Within-clade amino acid composition at the 4 significant GroupSim positions.
Shows whether Sat vs Trans differences hold independently within mixed clades
(Insects, Viridiplantae, Invertebrates) vs being a phylogenetic artefact.
"""

from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BASE    = Path(__file__).parent
OUT_DIR = BASE / "split_entropy" / "groupsim_sat_trans"

ITOL_SYM  = BASE / "iTOL_bnni" / "iTOL_centromere_symbols.txt"
ITOL_CLADE= BASE / "iTOL_bnni" / "iTOL_taxonomic_clade.txt"
ALN       = OUT_DIR / "trimmed_sat_trans_gap085.fasta"

# ── Load data ─────────────────────────────────────────────────────────────────
def read_fasta(path):
    seqs, cur, buf = {}, None, []
    for line in open(path):
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            if cur: seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]; buf = []
        else: buf.append(line)
    if cur: seqs[cur] = "".join(buf)
    return seqs

cent_type, tax_group = {}, {}
for line in open(ITOL_SYM):
    if line.startswith("#") or "\t" not in line: continue
    p = line.split("\t")
    if len(p) < 4: continue
    if p[3].strip() == "#ff006e": cent_type[p[0]] = "Satellite"
    elif p[3].strip() == "#3a86ff": cent_type[p[0]] = "Transposon"

for line in open(ITOL_CLADE):
    if line.startswith("#") or "\t" not in line: continue
    p = line.split("\t")
    if len(p) < 3: continue
    tax_group[p[0]] = p[2].strip()

seqs = read_fasta(ALN)

# ── Positions & biochemical class colours ────────────────────────────────────
SIG_POS = [
    (187, 4.27, "loop\n(α2–α3)"),
    (118, 3.49, "αN\nhelix"),
    (176, 3.18, "α2\nhelix"),
    (154, 2.45, "α2\nhelix"),
]

# Biochemical class → colour
CHEM_CLASS = {
    # Hydrophobic aliphatic
    "A": ("Hydrophobic", "#F4A261"),
    "V": ("Hydrophobic", "#F4A261"),
    "I": ("Hydrophobic", "#F4A261"),
    "L": ("Hydrophobic", "#F4A261"),
    "M": ("Hydrophobic", "#F4A261"),
    # Aromatic
    "F": ("Aromatic",    "#8338EC"),
    "Y": ("Aromatic",    "#8338EC"),
    "W": ("Aromatic",    "#8338EC"),
    # Polar uncharged
    "S": ("Polar",       "#52B788"),
    "T": ("Polar",       "#52B788"),
    "N": ("Polar",       "#52B788"),
    "Q": ("Polar",       "#52B788"),
    # Positive
    "K": ("Positive",    "#4A90D9"),
    "R": ("Positive",    "#4A90D9"),
    "H": ("Positive",    "#4A90D9"),
    # Negative
    "D": ("Negative",    "#E63946"),
    "E": ("Negative",    "#E63946"),
    # Special
    "G": ("Special",     "#BDBDBD"),
    "P": ("Special (Pro)","#FF6B6B"),
    "C": ("Special",     "#BDBDBD"),
}

def aa_colour(aa):
    return CHEM_CLASS.get(aa, ("Other", "#DDDDDD"))[1]

# ── Groups to show (order matters for display) ────────────────────────────────
GROUPS = [
    ("All",                    None),
    ("Insects",                "Insects"),
    ("Viridiplantae",          "Viridiplantae"),
    ("Invertebrates",          "Invertebrates (non-insect)"),
]

# ── Build AA count tables ─────────────────────────────────────────────────────
def aa_fracs(tips, pos_idx):
    c = Counter(seqs[t][pos_idx] for t in tips if t in seqs and seqs[t][pos_idx] != "-")
    total = sum(c.values())
    if total == 0:
        return {}, 0
    return {aa: n/total for aa, n in c.most_common()}, total

# ── Plot ──────────────────────────────────────────────────────────────────────
n_pos   = len(SIG_POS)
n_grp   = len(GROUPS)
fig, axes = plt.subplots(
    n_grp, n_pos,
    figsize=(4.2 * n_pos, 2.8 * n_grp),
    gridspec_kw={"hspace": 0.55, "wspace": 0.35}
)
fig.patch.set_facecolor("white")

SAT_EDGE  = "#C62828"
TRANS_EDGE= "#1565C0"
BAR_WIDTH = 0.38
BAR_SEP   = 0.55   # centre-to-centre distance between Sat and Trans bars

for col, (pos, z, region_lbl) in enumerate(SIG_POS):
    idx = pos - 1

    for row, (grp_label, grp_filter) in enumerate(GROUPS):
        ax = axes[row, col]
        ax.set_facecolor("white")

        # Get tips for this group
        if grp_filter is None:
            sat_tips  = [t for t, ct in cent_type.items() if ct == "Satellite"]
            trans_tips= [t for t, ct in cent_type.items() if ct == "Transposon"]
        else:
            sat_tips  = [t for t, ct in cent_type.items()
                         if ct == "Satellite" and tax_group.get(t) == grp_filter]
            trans_tips= [t for t, ct in cent_type.items()
                         if ct == "Transposon" and tax_group.get(t) == grp_filter]

        sat_frac,  n_sat  = aa_fracs(sat_tips,  idx)
        trans_frac,n_trans= aa_fracs(trans_tips, idx)

        # All AAs present in either group, ordered by Sat fraction desc
        all_aas = sorted(set(sat_frac) | set(trans_frac),
                         key=lambda a: -sat_frac.get(a, 0))

        # Draw stacked bars
        for bar_x, fracs, edge, n_label in [
            (0,        sat_frac,   SAT_EDGE,   f"Sat\n(n={n_sat})"),
            (BAR_SEP,  trans_frac, TRANS_EDGE, f"Trans\n(n={n_trans})"),
        ]:
            bottom = 0.0
            for aa in all_aas:
                h = fracs.get(aa, 0.0)
                if h < 0.003: continue
                col_c = aa_colour(aa)
                ax.bar(bar_x, h, BAR_WIDTH, bottom=bottom,
                       color=col_c, edgecolor="white", linewidth=0.4, zorder=3)
                if h >= 0.10:
                    ax.text(bar_x, bottom + h/2, aa,
                            ha="center", va="center",
                            fontsize=7.5, fontweight="bold", color="white",
                            zorder=4)
                bottom += h
            # border
            ax.bar(bar_x, 1.0, BAR_WIDTH, bottom=0,
                   color="none", edgecolor=edge, linewidth=1.4, zorder=5)
            # n label below
            ax.text(bar_x, -0.12, n_label,
                    ha="center", va="top", fontsize=6.5, color=edge,
                    fontweight="bold")

        ax.set_xlim(-BAR_WIDTH, BAR_SEP + BAR_WIDTH)
        ax.set_ylim(0, 1.05)
        ax.set_xticks([])
        ax.set_yticks([0, 0.5, 1.0])
        ax.tick_params(axis="y", labelsize=7)
        ax.spines[["top", "right", "bottom"]].set_visible(False)
        ax.spines["left"].set_linewidth(0.7)

        # Row label (group name) on left edge
        if col == 0:
            ax.set_ylabel(grp_label, fontsize=9, fontweight="bold", labelpad=6)

        # Column header (pos + Z + region) on top row only
        if row == 0:
            ax.set_title(
                f"pos {pos}  (Z = {z:.2f})\n{region_lbl}",
                fontsize=9, fontweight="bold", pad=6
            )

        # Verdict annotation bottom-right of each panel
        if row > 0:  # skip "All" row
            if   (col == 0 and grp_filter in ("Insects", "Viridiplantae")):
                verdict, vc = "≈ artefact", "#888888"
            elif (col == 1 and grp_filter in ("Insects",)):
                verdict, vc = "≈ artefact", "#888888"
            elif (col == 1 and grp_filter == "Viridiplantae"):
                verdict, vc = "subtle", "#E65100"
            elif (col == 2 and grp_filter == "Insects"):
                verdict, vc = "REVERSED", "#C62828"
            elif (col == 2 and grp_filter == "Viridiplantae"):
                verdict, vc = "holds", "#2E7D32"
            elif (col == 3 and grp_filter in ("Insects", "Viridiplantae")):
                verdict, vc = "holds", "#2E7D32"
            else:
                verdict, vc = "", "black"
            if verdict:
                ax.text(0.97, 0.97, verdict,
                        transform=ax.transAxes,
                        ha="right", va="top", fontsize=7, color=vc,
                        fontweight="bold")

# ── Legend (biochemical classes) ─────────────────────────────────────────────
class_colours = {
    "Hydrophobic": "#F4A261",
    "Aromatic":    "#8338EC",
    "Polar":       "#52B788",
    "Positive":    "#4A90D9",
    "Negative":    "#E63946",
    "Pro":         "#FF6B6B",
    "Other":       "#BDBDBD",
}
legend_patches = [mpatches.Patch(facecolor=c, edgecolor="grey", linewidth=0.5, label=lbl)
                  for lbl, c in class_colours.items()]
fig.legend(
    handles=legend_patches, title="Biochemical class",
    loc="lower center", ncol=7, fontsize=8, title_fontsize=8.5,
    framealpha=0.9, edgecolor="grey",
    bbox_to_anchor=(0.5, -0.01)
)

fig.suptitle(
    "Within-clade AA composition at GroupSim significant positions\n"
    "Satellite (red border) vs Transposon (blue border) CENH3  —  gap ≤ 85% trimmed alignment",
    fontsize=11, fontweight="bold", y=1.01
)

for ext in [".png", ".pdf"]:
    out = OUT_DIR / f"groupsim_within_clade{ext}"
    plt.savefig(out, dpi=200 if ext == ".png" else 100,
                bbox_inches="tight", facecolor="white")
    print(f"Saved {out.name}")
plt.close()
