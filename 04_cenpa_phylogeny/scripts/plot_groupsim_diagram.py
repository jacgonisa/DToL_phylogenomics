#!/usr/bin/env python3
"""
Schematic diagram of the GroupSim SDP algorithm including HH weighting.
Outputs: figures/groupsim_diagram.{pdf,png}
"""

from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OUT_DIR = Path(__file__).parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Palette ───────────────────────────────────────────────────────────────────
C_CENPA   = "#f8bbd0"
C_H3      = "#bbdefb"
C_CENPA_D = "#c62828"
C_H3_D    = "#0d47a1"
C_HL      = "#fff9c4"
C_WITHIN  = "#2e7d32"
C_BETWEEN = "#e65100"
C_SDP     = "#e53935"
C_NOSIG   = "#90a4ae"
C_HH      = "#6a1b9a"   # purple for HH weighting

fig = plt.figure(figsize=(15, 12))
fig.patch.set_facecolor("white")

ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.axis("off")

# ── Helpers ───────────────────────────────────────────────────────────────────
def box(x, y, w, h, fc="#eceff1", ec="#bbbbbb", lw=0.8, zorder=1):
    ax.add_patch(mpatches.FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.01",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=zorder,
        transform=ax.transAxes, clip_on=False))

def txt(x, y, s, **kw):
    ax.text(x, y, s, transform=ax.transAxes, clip_on=False, **kw)

def arr(x0, y0, x1, y1, color="#9e9e9e", lw=1.3):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                xycoords="axes fraction", textcoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=color,
                               lw=lw, mutation_scale=10))

# ── Row y-coordinates (bottom of each row) ────────────────────────────────────
# Row 0 — Alignment        0.80 – 0.97
# Row 1 — HH weighting     0.55 – 0.76
# Row 2 — Score calc       0.28 – 0.50
# Row 3 — Pipeline         0.01 – 0.24

R = {
    0: (0.80, 0.17),   # (y_bottom, height)
    1: (0.55, 0.21),
    2: (0.27, 0.23),
    3: (0.01, 0.22),
}

# ── Title ─────────────────────────────────────────────────────────────────────
txt(0.5, 0.985, "GroupSim — Specificity Determining Position (SDP) Algorithm",
    ha="center", va="center", fontsize=13, fontweight="bold", color="#212121")

# ═══════════════════════════════════════════════════════════════════════════════
# ROW 0 — ALIGNMENT
# ═══════════════════════════════════════════════════════════════════════════════
r0y, r0h = R[0]
box(0.02, r0y - 0.01, 0.96, r0h + 0.02, fc="white", ec="#cccccc")
txt(0.5, r0y + r0h + 0.005,
    "① Multiple Sequence Alignment  (two groups)",
    ha="center", va="bottom", fontsize=9.5, fontweight="bold", color="#424242")

seqs_cenpa = [list("IWAFQSLL"), list("WVAFRSLA"),
              list("IWAYQSLV"), list("WWAFESLL")]
seqs_h3    = [list("FFAGQSAV"), list("FFAGPSAV"),
              list("FFAGQSAI"), list("LFAGQSAV")]
n_pos      = 8
col_hl     = 0
label_x0   = 0.06
seq_x0     = 0.18
col_w      = 0.072
row_h      = 0.022
seq_y0     = r0y + r0h - 0.01

txt(label_x0, seq_y0 - 1.5 * row_h, "CENP-A\ngroup",
    ha="center", va="center", fontsize=8.5, fontweight="bold", color=C_CENPA_D)
txt(label_x0, seq_y0 - (len(seqs_cenpa) + 1.5) * row_h - 0.008, "H3\ngroup",
    ha="center", va="center", fontsize=8.5, fontweight="bold", color=C_H3_D)

for si, seq in enumerate(seqs_cenpa + seqs_h3):
    is_cenpa = si < len(seqs_cenpa)
    bg  = C_CENPA if is_cenpa else C_H3
    sep = 0.008 if si == len(seqs_cenpa) else 0
    y_s = seq_y0 - si * row_h - sep
    for pi, aa in enumerate(seq):
        xc = seq_x0 + pi * col_w
        hl = (pi == col_hl)
        fc = C_HL if hl else bg
        ec = "#f57f17" if hl else "white"
        lw = 1.5 if hl else 0.3
        ax.add_patch(mpatches.FancyBboxPatch(
            (xc, y_s - row_h + 0.002), col_w - 0.004, row_h - 0.003,
            boxstyle="round,pad=0.003",
            facecolor=fc, edgecolor=ec, linewidth=lw,
            transform=ax.transAxes, clip_on=False, zorder=2))
        ax.text(xc + col_w/2, y_s - row_h/2 + 0.001, aa,
                ha="center", va="center", fontsize=8,
                fontweight="bold" if hl else "normal",
                color="#b71c1c" if hl else "#212121",
                transform=ax.transAxes, clip_on=False, zorder=3)

for pi in range(n_pos):
    txt(seq_x0 + pi * col_w + col_w/2, seq_y0 + 0.009, str(pi + 1),
        ha="center", va="bottom", fontsize=7, color="#757575")

x_hl = seq_x0 + col_hl * col_w
txt(x_hl + col_w/2, r0y - 0.022, "scored\ncolumn",
    ha="center", va="top", fontsize=7, color="#f57f17", fontstyle="italic")
arr(x_hl + col_w/2, r0y - 0.012, x_hl + col_w/2, r0y + 0.003, "#f57f17")

txt(seq_x0 + n_pos * col_w + 0.02, seq_y0 - 1.5 * row_h,
    "n=418", ha="left", va="center", fontsize=7, color=C_CENPA_D)
txt(seq_x0 + n_pos * col_w + 0.02,
    seq_y0 - (len(seqs_cenpa) + 1.5) * row_h - 0.008,
    "n=901", ha="left", va="center", fontsize=7, color=C_H3_D)

# ═══════════════════════════════════════════════════════════════════════════════
# ROW 1 — HH WEIGHTING
# ═══════════════════════════════════════════════════════════════════════════════
r1y, r1h = R[1]
box(0.02, r1y - 0.01, 0.96, r1h + 0.02, fc="white", ec="#cccccc")
txt(0.5, r1y + r1h + 0.005,
    "② Henikoff-Henikoff sequence weights  (correct for phylogenetic over-representation)",
    ha="center", va="bottom", fontsize=9.5, fontweight="bold", color="#424242")

arr(0.5, r0y - 0.025, 0.5, r1y + r1h + 0.012)

# ── Left: problem — dense clade ───────────────────────────────────────────────
prob_ax = fig.add_axes([0.05, r1y + 0.03, 0.25, r1h - 0.05])
prob_ax.set_xlim(0, 1); prob_ax.set_ylim(0, 1); prob_ax.axis("off")
prob_ax.set_facecolor("#fce4ec")
prob_ax.add_patch(mpatches.FancyBboxPatch(
    (0.02, 0.02), 0.96, 0.96, boxstyle="round,pad=0.02",
    facecolor="#fce4ec", edgecolor="#ef9a9a", linewidth=1.2))
prob_ax.text(0.5, 0.91, "Without weighting", ha="center", va="center",
             fontsize=8, fontweight="bold", color=C_CENPA_D)

# show 5 identical "insect" seqs and 1 unique "plant" seq
rows_prob = [("W", "Insect 1", True), ("W", "Insect 2", True),
             ("W", "Insect 3", True), ("W", "Insect 4", True),
             ("W", "Insect 5", True), ("A", "Plant",    False)]
yy = 0.79
for aa, label, is_ins in rows_prob:
    c = "#ffcdd2" if is_ins else "#c8e6c9"
    prob_ax.add_patch(mpatches.FancyBboxPatch(
        (0.05, yy - 0.07), 0.12, 0.07,
        boxstyle="round,pad=0.01", facecolor=c, edgecolor="white", lw=0.5))
    prob_ax.text(0.11, yy - 0.035, aa, ha="center", va="center",
                 fontsize=9, fontweight="bold", color="#212121")
    prob_ax.text(0.25, yy - 0.035, label, ha="left", va="center",
                 fontsize=7, color="#424242")
    prob_ax.text(0.88, yy - 0.035, "w = 1", ha="right", va="center",
                 fontsize=7, color="#757575")
    yy -= 0.10

prob_ax.text(0.5, 0.04,
             "f(W) = 5/6 = 0.83  ← biased!",
             ha="center", va="center", fontsize=7.5,
             fontweight="bold", color=C_CENPA_D)

# ── Right: solution — HH weights ─────────────────────────────────────────────
sol_ax = fig.add_axes([0.38, r1y + 0.03, 0.56, r1h - 0.05])
sol_ax.set_xlim(0, 1); sol_ax.set_ylim(0, 1); sol_ax.axis("off")
sol_ax.set_facecolor("#f3e5f5")
sol_ax.add_patch(mpatches.FancyBboxPatch(
    (0.02, 0.02), 0.96, 0.96, boxstyle="round,pad=0.02",
    facecolor="#f3e5f5", edgecolor="#ce93d8", linewidth=1.2))
sol_ax.text(0.5, 0.91, "With HH weighting", ha="center", va="center",
            fontsize=8, fontweight="bold", color=C_HH)

# formula
sol_ax.text(0.5, 0.80,
            r"$w_i = \sum_c \frac{1}{r_c \times n_{a_c}}$",
            ha="center", va="center", fontsize=9, color="#212121")
sol_ax.text(0.5, 0.68,
            r"$r_c$ = distinct residues in col  |  $n_{a_c}$ = seqs sharing that residue",
            ha="center", va="center", fontsize=7, color="#555555")

# example: col has W×5, A×1  → r=2
# W: 1/(2×5)=0.1 each;  A: 1/(2×1)=0.5
rows_sol = [("W", "Insect 1", "1/(2×5) = 0.10"),
            ("W", "Insect 2", "1/(2×5) = 0.10"),
            ("W", "Insect 3", "1/(2×5) = 0.10"),
            ("W", "Insect 4", "1/(2×5) = 0.10"),
            ("W", "Insect 5", "1/(2×5) = 0.10"),
            ("A", "Plant",    "1/(2×1) = 0.50")]
yy = 0.57
for aa, label, w_str in rows_sol:
    c = "#ffcdd2" if aa == "W" else "#c8e6c9"
    sol_ax.add_patch(mpatches.FancyBboxPatch(
        (0.03, yy - 0.065), 0.07, 0.065,
        boxstyle="round,pad=0.01", facecolor=c, edgecolor="white", lw=0.5))
    sol_ax.text(0.065, yy - 0.032, aa, ha="center", va="center",
                fontsize=9, fontweight="bold", color="#212121")
    sol_ax.text(0.14, yy - 0.032, label, ha="left", va="center",
                fontsize=7, color="#424242")
    sol_ax.text(0.70, yy - 0.032, w_str, ha="left", va="center",
                fontsize=7.5, color=C_HH, family="monospace")
    yy -= 0.087

# normalise note
sol_ax.text(0.5, 0.04,
            "Normalised: Σw = n  →  f(W) = (5×0.10)/(5×0.10 + 0.50) = 0.50  ✓",
            ha="center", va="center", fontsize=7.5,
            fontweight="bold", color=C_HH)

arr(0.305, r1y + r1h/2 + 0.01, 0.37, r1y + r1h/2 + 0.01, C_HH)

# ═══════════════════════════════════════════════════════════════════════════════
# ROW 2 — SCORE CALCULATION
# ═══════════════════════════════════════════════════════════════════════════════
r2y, r2h = R[2]
box(0.02, r2y - 0.01, 0.96, r2h + 0.02, fc="white", ec="#cccccc")
txt(0.5, r2y + r2h + 0.005,
    "③ GroupSim score at each column  (identity matrix, using weighted frequencies)",
    ha="center", va="bottom", fontsize=9.5, fontweight="bold", color="#424242")

arr(0.5, r1y - 0.022, 0.5, r2y + r2h + 0.012)

# CENPA frequency bar
bar_c = fig.add_axes([0.06, r2y + 0.05, 0.17, 0.14])
cenpa_aa = ["W", "I/V", "other"]
cenpa_f  = [0.55, 0.24, 0.21]
colors_c = [C_CENPA_D, "#ef9a9a", "#ffcdd2"]
bars = bar_c.bar(cenpa_aa, cenpa_f, color=colors_c, edgecolor="white", width=0.55)
bar_c.set_ylim(0, 1.0)
bar_c.set_ylabel("Weighted\nfrequency", fontsize=7)
bar_c.set_title("CENP-A", fontsize=8.5, fontweight="bold", color=C_CENPA_D, pad=3)
bar_c.tick_params(labelsize=7)
bar_c.spines[["top", "right"]].set_visible(False)
bar_c.set_facecolor(C_CENPA + "55")
for b, f in zip(bars, cenpa_f):
    bar_c.text(b.get_x() + b.get_width()/2, f + 0.02, f"{f:.0%}",
               ha="center", va="bottom", fontsize=6.5, color="#424242")

within_c = sum(f**2 for f in cenpa_f)
txt(0.155, r2y + 0.015,
    r"Within$_{CENPA}$ = $\sum_a (f_a)^2$" + f" = {within_c:.3f}",
    ha="center", va="bottom", fontsize=8, color=C_WITHIN)

# H3 frequency bar
bar_h = fig.add_axes([0.30, r2y + 0.05, 0.17, 0.14])
h3_aa = ["F", "L/M", "other"]
h3_f  = [0.89, 0.06, 0.05]
colors_h = [C_H3_D, "#64b5f6", "#bbdefb"]
bars_h = bar_h.bar(h3_aa, h3_f, color=colors_h, edgecolor="white", width=0.55)
bar_h.set_ylim(0, 1.0)
bar_h.set_title("H3", fontsize=8.5, fontweight="bold", color=C_H3_D, pad=3)
bar_h.tick_params(labelsize=7)
bar_h.spines[["top", "right"]].set_visible(False)
bar_h.set_facecolor(C_H3 + "55")
for b, f in zip(bars_h, h3_f):
    bar_h.text(b.get_x() + b.get_width()/2, f + 0.02, f"{f:.0%}",
               ha="center", va="bottom", fontsize=6.5, color="#424242")

within_h = sum(f**2 for f in h3_f)
txt(0.385, r2y + 0.015,
    r"Within$_{H3}$ = $\sum_a (f_a)^2$" + f" = {within_h:.3f}",
    ha="center", va="bottom", fontsize=8, color=C_WITHIN)

# Between-group box
between = sum(fc * fh for fc, fh in zip(cenpa_f, h3_f))
bet_ax = fig.add_axes([0.55, r2y + 0.05, 0.16, 0.14])
bet_ax.set_xlim(0, 1); bet_ax.set_ylim(0, 1); bet_ax.axis("off")
bet_ax.add_patch(mpatches.FancyBboxPatch(
    (0.02, 0.02), 0.96, 0.96, boxstyle="round,pad=0.02",
    facecolor="#fff8e1", edgecolor=C_BETWEEN, linewidth=1.2))
bet_ax.text(0.5, 0.88, "Between groups", ha="center", va="center",
            fontsize=8, fontweight="bold", color=C_BETWEEN)
bet_ax.text(0.5, 0.70, r"$\sum_a f_a^{CENPA} \times f_a^{H3}$",
            ha="center", va="center", fontsize=8.5, color="#212121")
bet_ax.text(0.5, 0.50, "W×F: 0.55×0.00 = 0.00", ha="center", va="center",
            fontsize=6.5, color="#666666", family="monospace")
bet_ax.text(0.5, 0.37, "I×L: 0.24×0.06 = 0.01", ha="center", va="center",
            fontsize=6.5, color="#666666", family="monospace")
bet_ax.text(0.5, 0.18, f"Total = {between:.3f}", ha="center", va="center",
            fontsize=8, fontweight="bold", color=C_BETWEEN)

# Score box
mean_within = (within_c + within_h) / 2
raw = mean_within - between
sc_ax = fig.add_axes([0.79, r2y + 0.05, 0.16, 0.14])
sc_ax.set_xlim(0, 1); sc_ax.set_ylim(0, 1); sc_ax.axis("off")
sc_ax.add_patch(mpatches.FancyBboxPatch(
    (0.02, 0.02), 0.96, 0.96, boxstyle="round,pad=0.02",
    facecolor="#e8f5e9", edgecolor=C_WITHIN, linewidth=1.5))
sc_ax.text(0.5, 0.88, "Raw GroupSim", ha="center", va="center",
           fontsize=8, fontweight="bold", color=C_WITHIN)
sc_ax.text(0.5, 0.65, r"$\overline{\mathrm{within}} - \mathrm{between}$",
           ha="center", va="center", fontsize=11, color="#212121")
sc_ax.text(0.5, 0.44, f"= {mean_within:.3f} − {between:.3f}",
           ha="center", va="center", fontsize=8.5, color="#424242")
sc_ax.text(0.5, 0.22, f"= {raw:.3f}", ha="center", va="center",
           fontsize=11, fontweight="bold", color=C_WITHIN)
sc_ax.text(0.5, 0.08, "(high = discriminating)", ha="center", va="center",
           fontsize=6.5, color="#666666", fontstyle="italic")

# Arrows row2
arr(0.24,  r2y + 0.125, 0.545, r2y + 0.125, C_BETWEEN)
arr(0.475, r2y + 0.125, 0.545, r2y + 0.125, C_BETWEEN)
arr(0.715, r2y + 0.125, 0.785, r2y + 0.125, C_WITHIN)
txt(0.748, r2y + 0.165,
    r"$\overline{\mathrm{within}}$" + f" = {mean_within:.3f}",
    ha="center", va="bottom", fontsize=8, color=C_WITHIN)

# ═══════════════════════════════════════════════════════════════════════════════
# ROW 3 — PIPELINE
# ═══════════════════════════════════════════════════════════════════════════════
r3y, r3h = R[3]
box(0.02, r3y, 0.96, r3h, fc="white", ec="#cccccc")
txt(0.5, r3y + r3h + 0.005,
    "④ From raw score to SDP",
    ha="center", va="bottom", fontsize=9.5, fontweight="bold", color="#424242")

arr(0.5, r2y - 0.022, 0.5, r3y + r3h + 0.012)

steps = [
    ("Raw score\nper column",        "#fff9c4", "#f9a825"),
    ("Window correction\n(JS div, w=3)", "#fce4ec", "#c62828"),
    ("Normalise\nto [0, 1]",        "#e8eaf6", "#283593"),
    ("Z-score across\nall columns",  "#e0f2f1", "#00695c"),
    ("SDP\n(z ≥ 2)",                 "#ffebee", C_SDP),
]
n_s   = len(steps)
bw    = 0.12
bh    = 0.12
gx    = (0.96 - 0.04 - n_s * bw) / (n_s - 1)
sx    = 0.04
cy    = r3y + r3h/2 - bh/2 + 0.01

for i, (label, fc, ec) in enumerate(steps):
    bx = sx + i * (bw + gx)
    box(bx, cy, bw, bh, fc=fc, ec=ec, lw=1.5, zorder=2)
    txt(bx + bw/2, cy + bh/2, label,
        ha="center", va="center", fontsize=8.5, fontweight="bold",
        color=ec, zorder=3)
    if i < n_s - 1:
        arr(bx + bw + 0.004, cy + bh/2,
            bx + bw + gx - 0.004, cy + bh/2)

# Mini Manhattan
mini = fig.add_axes([0.06, r3y + 0.01, 0.88, 0.055])
np.random.seed(42)
nc = 60
sc = np.abs(np.random.normal(0.4, 0.18, nc))
sc = np.clip(sc, 0, 1)
for p in [8, 21, 35, 47]:
    sc[p] = np.random.uniform(0.80, 0.99)
zz = (sc - sc.mean()) / sc.std()
colors_m = [C_SDP if z >= 2 else C_NOSIG for z in zz]
mini.bar(np.arange(nc), sc, color=colors_m, width=0.85)
z2s = sc.mean() + 2 * sc.std()
mini.axhline(z2s, color=C_SDP, lw=1.0, ls="--", alpha=0.7)
mini.text(nc - 0.5, z2s + 0.01, "z = 2", ha="right", va="bottom",
          fontsize=6.5, color=C_SDP)
mini.spines[["top", "right", "left"]].set_visible(False)
mini.set_xticks([]); mini.set_yticks([])
mini.set_xlabel("Alignment position", fontsize=7.5, labelpad=2)
mini.set_ylabel("Score", fontsize=7)
mini.set_facecolor("white")
mini.legend(handles=[
    mpatches.Patch(facecolor=C_SDP,   edgecolor="none", label="SDP (z ≥ 2)"),
    mpatches.Patch(facecolor=C_NOSIG, edgecolor="none", label="Non-significant"),
], loc="upper left", fontsize=6.5, framealpha=0.9,
   edgecolor="#cccccc", handlelength=1.0, borderpad=0.3, ncol=2)

# ── Save ──────────────────────────────────────────────────────────────────────
for ext in ("pdf", "png"):
    p = OUT_DIR / f"groupsim_diagram.{ext}"
    fig.savefig(p, dpi=300 if ext == "png" else None,
                bbox_inches="tight", facecolor="white")
    print(f"Saved: {p}")
plt.close(fig)
