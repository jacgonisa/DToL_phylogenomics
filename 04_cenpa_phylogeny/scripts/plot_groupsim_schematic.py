#!/usr/bin/env python3
"""
GroupSim algorithm schematic.

Layout:
  A: MSA input (two groups, one column highlighted)
  B: AA frequency profiles for that column
  C: Jensen-Shannon divergence formula + score
  D: Full per-position score profile

Saves: groupsim_schematic.{svg, pdf, png}
SVG is fully editable in Inkscape / Illustrator.
"""

from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.ticker
from matplotlib.patches import FancyArrowPatch

OUT_DIR = Path(__file__).parent / "split_entropy" / "groupsim_sat_trans"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Chemical-class colour map ─────────────────────────────────────────────────
_CLASS = {}
for a in "ACFGILMPVWY": _CLASS[a] = "hydrophobic"
for a in "NQST":         _CLASS[a] = "polar"
for a in "HKR":          _CLASS[a] = "positive"
for a in "DE":            _CLASS[a] = "negative"
_CLASS["G"] = "special"
_CLASS["-"] = "gap"

_CCOL = {
    "hydrophobic": "#F4A261",
    "polar":       "#52B788",
    "positive":    "#4A90D9",
    "negative":    "#E63946",
    "special":     "#BDBDBD",
}

def aa_col(aa):
    return _CCOL.get(_CLASS.get(aa.upper(), "gap"), "#EEEEEE")

# ── Mock alignment (controlled to illustrate the concept clearly) ─────────────
rng = np.random.default_rng(42)

N_SAT, N_TRANS, N_COLS = 10, 8, 22
HCOL = 13          # highlighted "divergent" column (0-indexed)
MOD_COLS = [4, 18] # moderately divergent columns

HYD  = list("LIVMFA")    # Satellite: hydrophobic at HCOL
POS  = list("KRH")        # Transposon: positively charged at HCOL
MIX  = list("ACDEFGILMNPQRSTVWY")

def gen_seqs(n, hi_pool):
    seqs = []
    for _ in range(n):
        row = [rng.choice(MIX) for _ in range(N_COLS)]
        row[HCOL] = rng.choice(hi_pool)
        for mc in MOD_COLS:
            row[mc] = rng.choice(hi_pool + MIX[:5])
        seqs.append(row)
    return seqs

sat_seqs   = gen_seqs(N_SAT,   HYD)
trans_seqs = gen_seqs(N_TRANS, POS)

# ── JS divergence helpers ─────────────────────────────────────────────────────
AAS = list("ACDEFGHIKLMNPQRSTVWY")

def freq_vec(col_aas):
    c = {a: 0 for a in AAS}
    n = 0
    for a in col_aas:
        if a in c:
            c[a] += 1; n += 1
    v = np.array([c[a] for a in AAS], dtype=float)
    return v / max(n, 1)

def js_div(p, q, eps=1e-9):
    p = p + eps; p /= p.sum()
    q = q + eps; q /= q.sum()
    m = (p + q) / 2.0
    kl = lambda a, b: np.sum(a * np.log(a / b))
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)

sat_cols   = [[sat_seqs[r][c]   for r in range(N_SAT)]   for c in range(N_COLS)]
trans_cols = [[trans_seqs[r][c] for r in range(N_TRANS)] for c in range(N_COLS)]

raw = np.array([js_div(freq_vec(sat_cols[c]), freq_vec(trans_cols[c]))
                for c in range(N_COLS)])
scores = (raw - raw.min()) / (raw.max() - raw.min() + 1e-10)

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(17, 9.5), facecolor="white")

gs_main = gridspec.GridSpec(
    2, 1, figure=fig,
    height_ratios=[2.8, 1.9],
    hspace=0.40
)
gs_top = gridspec.GridSpec(
    1, 4, figure=fig,
    width_ratios=[3.6, 0.25, 3.2, 1.8],
    wspace=0.08,
    left=0.07, right=0.97, top=0.91, bottom=0.42
)
gs_prof = gridspec.GridSpecFromSubplotSpec(
    2, 1, subplot_spec=gs_top[2], hspace=0.10
)
gs_bot = gridspec.GridSpec(
    1, 1, figure=fig,
    left=0.07, right=0.97, top=0.37, bottom=0.08
)

ax_msa     = fig.add_subplot(gs_top[0])
ax_arr1    = fig.add_subplot(gs_top[1])   # arrow placeholder
ax_sat_p   = fig.add_subplot(gs_prof[0])
ax_trans_p = fig.add_subplot(gs_prof[1])
ax_score   = fig.add_subplot(gs_top[3])
ax_full    = fig.add_subplot(gs_bot[0])

ax_arr1.axis("off")

# ─────────────────────────────────────────────────────────────────────────────
# Panel A — Mini MSA
# ─────────────────────────────────────────────────────────────────────────────
SEP = 1           # blank separator row between groups
n_draw = N_SAT + SEP + N_TRANS

ax_msa.set_xlim(-0.6, N_COLS + 0.1)
ax_msa.set_ylim(-0.8, n_draw + 0.2)
ax_msa.axis("off")

# Group background bands
ax_msa.add_patch(mpatches.FancyBboxPatch(
    (-0.45, N_TRANS + SEP - 0.05), N_COLS + 0.4, N_SAT + 0.15,
    boxstyle="round,pad=0.08", fc="#FFF0F3", ec="#ff006e",
    lw=1.6, zorder=0, clip_on=False
))
ax_msa.add_patch(mpatches.FancyBboxPatch(
    (-0.45, -0.05), N_COLS + 0.4, N_TRANS + 0.1,
    boxstyle="round,pad=0.08", fc="#EEF4FF", ec="#3a86ff",
    lw=1.6, zorder=0, clip_on=False
))

# Highlighted column
ax_msa.add_patch(mpatches.FancyBboxPatch(
    (HCOL - 0.08, -0.15), 1.18, n_draw + 0.30,
    boxstyle="round,pad=0.04", fc="#FFE566", ec="#C62828",
    lw=2.2, alpha=0.70, zorder=1, clip_on=False
))

# Draw cells
r_abs = 0
for group_seqs, group_offset in [(sat_seqs, N_TRANS + SEP), (trans_seqs, 0)]:
    for r_idx, row_data in enumerate(group_seqs):
        y = group_offset + len(group_seqs) - r_idx - 1
        for c_idx, aa in enumerate(row_data):
            rect = mpatches.Rectangle(
                (c_idx + 0.05, y + 0.05), 0.90, 0.90,
                fc=aa_col(aa), ec="white", lw=0.35, zorder=2
            )
            ax_msa.add_patch(rect)
            ax_msa.text(c_idx + 0.50, y + 0.50, aa,
                        ha="center", va="center", fontsize=5.2,
                        fontweight="bold", color="white",
                        fontfamily="monospace", zorder=3)

# Group labels (left margin)
ax_msa.text(-0.55, N_TRANS + SEP + N_SAT / 2,
            "Satellite\n(n=198)", ha="right", va="center",
            fontsize=9, fontweight="bold", color="#ff006e")
ax_msa.text(-0.55, N_TRANS / 2,
            "Transposon\n(n=129)", ha="right", va="center",
            fontsize=9, fontweight="bold", color="#3a86ff")

# Column label
ax_msa.text(HCOL + 0.5, -0.55, f"pos {HCOL + 1}",
            ha="center", va="top", fontsize=8, color="#C62828",
            fontweight="bold")

# Column index tick marks
for c in range(0, N_COLS, 5):
    ax_msa.text(c + 0.5, -0.3, str(c + 1),
                ha="center", va="top", fontsize=5.5, color="#888888")

ax_msa.set_title(
    "A   Input: MSA of CENH3 sequences — two functional groups",
    fontsize=10, fontweight="bold", loc="left", pad=5
)

# ─────────────────────────────────────────────────────────────────────────────
# Arrow A → B
# ─────────────────────────────────────────────────────────────────────────────
ax_arr1.text(0.5, 0.50, "→", ha="center", va="center",
             fontsize=26, color="#555555", transform=ax_arr1.transAxes)

# ─────────────────────────────────────────────────────────────────────────────
# Panel B — AA frequency profiles
# ─────────────────────────────────────────────────────────────────────────────
p_vec = freq_vec(sat_cols[HCOL])
q_vec = freq_vec(trans_cols[HCOL])

nz      = np.where((p_vec > 0.02) | (q_vec > 0.02))[0]
aa_lbls = [AAS[i] for i in nz]
p_v     = p_vec[nz]
q_v     = q_vec[nz]
bar_col = [aa_col(a) for a in aa_lbls]
xi      = np.arange(len(aa_lbls))

for ax, vals, grp_label, grp_col, bg_col, is_top in [
    (ax_sat_p,   p_v, "P   (Satellite)",  "#ff006e", "#FFF0F3", True),
    (ax_trans_p, q_v, "Q   (Transposon)", "#3a86ff", "#EEF4FF", False),
]:
    ax.bar(xi, vals, color=bar_col, edgecolor="white", linewidth=0.5, width=0.82)
    ax.set_xticks(xi)
    ax.set_xticklabels(aa_lbls, fontsize=6.5, fontfamily="monospace")
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("Freq.", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", length=0, bottom=False)
    ax.text(0.97, 0.88, grp_label, transform=ax.transAxes,
            ha="right", va="top", fontsize=8.5, fontweight="bold", color=grp_col,
            bbox=dict(fc=bg_col, ec=grp_col, lw=0.8, pad=2.5))
    if is_top:
        ax.set_title(f"B   AA profiles at highlighted position {HCOL + 1}",
                     fontsize=10, fontweight="bold", loc="left", pad=4)
        ax.tick_params(labelbottom=False)

# Chemical class legend on bottom profile panel
patches = [mpatches.Patch(fc=c, label=k.capitalize())
           for k, c in _CCOL.items()]
ax_trans_p.legend(handles=patches, fontsize=6, ncol=3,
                   loc="upper left", framealpha=0.85,
                   title="Chemical class", title_fontsize=6.5,
                   borderpad=0.4, handlelength=1.0)

# ─────────────────────────────────────────────────────────────────────────────
# Panel C — JS divergence formula + score
# ─────────────────────────────────────────────────────────────────────────────
ax_score.axis("off")
sv = scores[HCOL]

# Step label
ax_score.text(0.05, 0.99, "C   Divergence score",
              transform=ax_score.transAxes, ha="left", va="top",
              fontsize=10, fontweight="bold")

# Formula box
formula_txt = (
    "For each column $i$:\n\n"
    r"$M_i = \frac{P_i + Q_i}{2}$" + "\n\n"
    r"$\mathrm{JS}(P_i \| Q_i) =$" + "\n"
    r"$\;\frac{1}{2}D_{KL}(P_i\|M_i)$" + "\n"
    r"$+\frac{1}{2}D_{KL}(Q_i\|M_i)$"
)
ax_score.text(0.50, 0.64, formula_txt,
              transform=ax_score.transAxes,
              ha="center", va="center",
              fontsize=9.0, linespacing=1.75,
              bbox=dict(boxstyle="round,pad=0.55",
                        fc="#FFFDF0", ec="#AAAAAA", lw=1.2))

# Score badge
badge_col = "#d62728" if sv > 0.65 else ("#F4A261" if sv > 0.35 else "#4A90D9")
ax_score.add_patch(mpatches.FancyBboxPatch(
    (0.08, 0.04), 0.84, 0.20,
    transform=ax_score.transAxes,
    boxstyle="round,pad=0.06",
    fc=badge_col, ec="#444444", lw=1.0, alpha=0.92
))
ax_score.text(0.50, 0.14, f"Score pos {HCOL+1} = {sv:.2f}",
              transform=ax_score.transAxes,
              ha="center", va="center", fontsize=9.5,
              fontweight="bold", color="white")

# ─────────────────────────────────────────────────────────────────────────────
# Panel D — Full score profile
# ─────────────────────────────────────────────────────────────────────────────
pos    = np.arange(1, N_COLS + 1)
thresh = np.percentile(scores, 70)
top    = scores >= thresh

ax_full.fill_between(pos, scores, where=~top,
                     color="#9ecae1", alpha=0.45, linewidth=0)
ax_full.fill_between(pos, scores, where=top,
                     color="#d62728", alpha=0.40, linewidth=0)
ax_full.plot(pos, scores, color="#1a1a2e", linewidth=1.9, zorder=4)
ax_full.scatter(pos[top], scores[top], s=60, color="#7B0000",
                zorder=5, linewidths=0, label="High-scoring positions")

# Mark highlighted column
ax_full.axvline(HCOL + 1, color="#F4A261", linewidth=2.2,
                linestyle="--", alpha=0.90, zorder=3)
ax_full.annotate(
    f"pos {HCOL+1}  (score = {sv:.2f})\n← computed in panels B+C",
    xy=(HCOL + 1, sv),
    xytext=(HCOL + 1 + 2.5, sv + 0.12),
    fontsize=8.5, color="#C62828", fontweight="bold",
    arrowprops=dict(arrowstyle="->", color="#C62828", lw=1.3)
)

# Mock helix band (show concept of structural annotation)
helix_regions = [(6, 9), (11, 15), (17, 21)]  # illustrative
for hs, he in helix_regions:
    ax_full.axvspan(hs + 0.5, he + 0.5, ymin=0.88, ymax=0.97,
                    color="#C62828", alpha=0.45, linewidth=0)
ax_full.text(helix_regions[0][0] + 0.5, 1.04, "α-helices (STRIDE)",
             fontsize=7.5, color="#C62828", fontweight="bold", va="bottom")

ax_full.set_xlim(0.3, N_COLS + 0.7)
ax_full.set_ylim(-0.06, 1.16)
ax_full.set_xlabel("Alignment position", fontsize=11)
ax_full.set_ylabel("GroupSim score\n(normalised JS divergence)", fontsize=10)
ax_full.set_title(
    "D   Output: per-position score — high score = position differentiates the two groups",
    fontsize=10, fontweight="bold", loc="left"
)
ax_full.spines[["top", "right"]].set_visible(False)
ax_full.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.25))
ax_full.grid(axis="y", color="grey", alpha=0.18, linewidth=0.6)
ax_full.legend(fontsize=8.5, loc="upper right", framealpha=0.9)

# ─────────────────────────────────────────────────────────────────────────────
# Connecting arrow from score panel (C) down to full profile (D)
# ─────────────────────────────────────────────────────────────────────────────
pos_c = ax_score.get_position()
pos_d = ax_full.get_position()
x_mid = (pos_c.x0 + pos_c.x1) / 2

arr = FancyArrowPatch(
    (x_mid, pos_c.y0 - 0.005),
    (x_mid, pos_d.y1 + 0.005),
    transform=fig.transFigure, figure=fig,
    arrowstyle="-|>", mutation_scale=16,
    color="#555555", linewidth=1.6
)
fig.add_artist(arr)

# ─────────────────────────────────────────────────────────────────────────────
# Main title
# ─────────────────────────────────────────────────────────────────────────────
fig.text(0.50, 0.975,
         "GroupSim — identifying positions that differentiate two functional groups"
         " by amino acid composition",
         ha="center", va="top", fontsize=12, fontweight="bold")

fig.text(0.50, 0.955,
         "Pazos et al. 2006  |  per-column Jensen-Shannon divergence between "
         "group amino acid frequency profiles",
         ha="center", va="top", fontsize=8.5, color="#555555")

# ─────────────────────────────────────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────────────────────────────────────
for ext in [".svg", ".pdf", ".png"]:
    kwargs = {"bbox_inches": "tight"}
    if ext == ".png":
        kwargs["dpi"] = 200
    out = OUT_DIR / f"groupsim_schematic{ext}"
    plt.savefig(out, **kwargs)
    print(f"Saved {out.name}")
plt.close()
