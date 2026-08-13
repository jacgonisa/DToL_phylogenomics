#!/usr/bin/env python3
"""
Publication-ready GroupSim Manhattan plots — unified style.

Panel A: CENP-A (n=418) vs H3-like (n=901), gap 80%
Panel B: Satellite CENP-A (n=198) vs Transposon CENP-A (n=129), gap 85%

Both use:
  • Top track   — helix annotation blocks
  • Main track  — gray = unweighted, coloured = HH-weighted, red = z ≥ 2
  • Bottom track — z-score heatmap strip

Outputs: figures/groupsim_cenpa_vs_h3_pub.{pdf,png}
         figures/groupsim_sat_vs_trans_pub.{pdf,png}
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
import matplotlib.ticker as ticker

BASE    = Path(__file__).parent
OUT_DIR = BASE / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Shared style ──────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.size":       10,
    "axes.labelsize":  10,
    "axes.titlesize":  11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
})

C_UW    = "#cfd8dc"   # unweighted — very light blue-grey
C_WT_A  = "#1565c0"   # weighted — deep blue (both plots)
C_WT_B  = "#1565c0"
C_SIG   = "#c62828"   # significant (z ≥ 2)
C_H3    = "#90caf9"   # H3 helix
C_CENPA = "#f48fb1"   # CENP-A helix


# ── Helpers ───────────────────────────────────────────────────────────────────
def read_uw(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t", 2)
        try:
            pos = int(parts[0])
        except (ValueError, IndexError):
            continue
        scor = None if parts[1].strip() in ("None","","NA") else float(parts[1])
        rows.append({"pos": pos + 1, "score_uw": scor})  # convert 0-based to 1-based
    return pd.DataFrame(rows)


C_CLADE = "#e65100"   # clade-weighted — deep orange

def make_groupsim_plot(
    uw_path, wt_path, helix_path,
    title, xlabel, out_stem,
    c_wt, annotations,
    wt_label="Weighted",     # label for weighted track in legend
    helix_mode="two_rows",   # "two_rows" (H3+CENPA) | "one_row" (CENPA only)
    top_n_auto=0,            # if >0, auto-annotate top N weighted positions
    wt2_path=None,           # optional second weighted method (clade)
    c_wt2=None,
    fig_w=13, fig_h=6.5,
):
    # load
    uw  = read_uw(uw_path)
    wt  = pd.read_csv(wt_path, sep="\t").rename(columns={
              "groupsim_weighted": "score_wt", "groupsim_hh": "score_wt",
              "groupsim_clade": "score_wt",
              "z_score": "z_score", "z_hh": "z_score", "z_clade": "z_score"})
    df  = uw.merge(wt[["pos","score_wt","z_score"]], on="pos", how="outer").sort_values("pos")
    df["sig"] = (df["z_score"] >= 2.0).fillna(False)
    if wt2_path is not None:
        wt2 = pd.read_csv(wt2_path, sep="\t").rename(columns={"groupsim_clade":"score_wt2",
                                                                "z_clade":"z_score2"})
        df = df.merge(wt2[["pos","score_wt2","z_score2"]], on="pos", how="left")
        df["sig2"] = (df["z_score2"] >= 2.0).fillna(False)
    L    = len(df)
    xmin = df["pos"].min() - 1
    xmax = df["pos"].max() + 1

    # helix file — two formats: TSV with histone col, or plain start/end
    raw = [l for l in Path(helix_path).read_text().splitlines()
           if l.strip() and not l.startswith("#")]
    if "\t" in raw[0] and "histone" in raw[0]:
        helix_df = pd.read_csv(helix_path, sep="\t")
        hx_h3    = helix_df[helix_df["histone"] == "H3"]
        hx_cenpa = helix_df[helix_df["histone"] == "CENPA"]
    else:
        rows = [r.split() for r in raw if not r.startswith("start")]
        hx_cenpa = pd.DataFrame(rows, columns=["start","end"]).astype(int)
        hx_h3    = None

    # layout
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs  = fig.add_gridspec(3, 1, height_ratios=[0.55, 4.5, 0.5],
                            hspace=0.06, left=0.07, right=0.94,
                            top=0.91, bottom=0.09)
    ax_hx = fig.add_subplot(gs[0])
    ax_mn = fig.add_subplot(gs[1], sharex=ax_hx)
    ax_z  = fig.add_subplot(gs[2], sharex=ax_hx)

    # ── helix track ──────────────────────────────────────────────────────────
    ax_hx.set_xlim(xmin, xmax)
    ax_hx.set_ylim(0, 1)
    ax_hx.axis("off")

    helix_names = ["αN","α1","α2","α3"]

    def draw_helix_row(hx, color, y0, h, label_color, show_labels):
        for i, (_, row) in enumerate(hx.iterrows()):
            s, e = row["start"], row["end"]
            ax_hx.add_patch(mpatches.FancyBboxPatch(
                (s - 0.4, y0), e - s + 0.8, h,
                boxstyle="round,pad=0,rounding_size=0.6",
                facecolor=color, edgecolor="none", alpha=0.88, zorder=2
            ))
            if show_labels:
                lbl = helix_names[i] if i < len(helix_names) else f"α{i+1}"
                ax_hx.text((s + e) / 2, y0 + h / 2, lbl,
                           ha="center", va="center", fontsize=8,
                           fontweight="bold", color=label_color, zorder=3)

    if helix_mode == "two_rows" and hx_h3 is not None:
        draw_helix_row(hx_h3,    C_H3,    0.54, 0.38, "#0d47a1", False)
        draw_helix_row(hx_cenpa, C_CENPA, 0.08, 0.38, "#880e4f", False)
        leg_handles = [
            mpatches.Patch(facecolor=C_H3,    edgecolor="none", label="H3 α-helices"),
            mpatches.Patch(facecolor=C_CENPA, edgecolor="none", label="CENP-A α-helices"),
        ]
    else:
        draw_helix_row(hx_cenpa, C_CENPA, 0.08, 0.84, "#880e4f", True)
        leg_handles = [
            mpatches.Patch(facecolor=C_CENPA, edgecolor="none", label="CENP-A α-helices"),
        ]

    ax_hx.legend(handles=leg_handles, loc="upper right", fontsize=8,
                 framealpha=0.9, edgecolor="#dddddd", borderpad=0.4,
                 handlelength=1.0, ncol=len(leg_handles))

    # ── main panel ───────────────────────────────────────────────────────────
    ax_mn.set_ylim(0, 1.30)
    ax_mn.set_ylabel("GroupSim score (normalised)", fontsize=10)
    ax_mn.axhline(0, color="#e0e0e0", lw=0.5, zorder=0)
    ax_mn.spines[["top","right","bottom"]].set_visible(False)
    ax_mn.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_mn.yaxis.grid(True, color="#f0f0f0", lw=0.6, zorder=0)
    ax_mn.yaxis.set_major_locator(ticker.MultipleLocator(0.25))

    # bar widths: 2 or 3 tracks
    n_tracks = 3 if wt2_path is not None else 2
    bw_each  = 0.80 / n_tracks

    if n_tracks == 2:
        offsets = [-bw_each/2, bw_each/2]
    else:
        offsets = [-bw_each, 0, bw_each]

    # unweighted (background) — only where score exists
    uw_v = df.dropna(subset=["score_uw"])
    ax_mn.bar(uw_v["pos"] + offsets[0], uw_v["score_uw"],
              width=bw_each, color=C_UW, alpha=1.0, zorder=2, label="Unweighted")

    # weighted (foreground) — only where score exists
    wt_v = df.dropna(subset=["score_wt"])
    bar_colors = [C_SIG if s else c_wt for s in wt_v["sig"]]
    ax_mn.bar(wt_v["pos"] + offsets[1], wt_v["score_wt"],
              width=bw_each, color=bar_colors, alpha=0.92, zorder=3, label=wt_label)

    # clade-weighted (third track, optional) — only where score exists
    if wt2_path is not None:
        c2   = c_wt2 if c_wt2 else C_CLADE
        wt2_v = df.dropna(subset=["score_wt2"])
        bar_colors2 = [C_SIG if s else c2 for s in wt2_v["sig2"]]
        ax_mn.bar(wt2_v["pos"] + offsets[2], wt2_v["score_wt2"],
                  width=bw_each, color=bar_colors2, alpha=0.92, zorder=3, label="Clade-weighted")

    # z=2 dashed line
    valid = df.dropna(subset=["score_wt"])
    if valid["sig"].any():
        z2s = valid.loc[valid["sig"], "score_wt"].min()
        ax_mn.axhline(z2s, color=C_SIG, lw=0.9, ls="--", alpha=0.55, zorder=1)
    else:
        mu, sd = valid["score_wt"].mean(), valid["score_wt"].std()
        ax_mn.axhline(min(mu + 2*sd, 1.0), color="#9e9e9e", lw=0.8,
                      ls="--", alpha=0.5, zorder=1)

    # annotations — explicit dict takes priority; else auto top-N
    annot_positions = annotations if annotations else {}
    if top_n_auto > 0 and not annotations:
        top_rows = valid.nlargest(top_n_auto, "score_wt")
        for _, r in top_rows.iterrows():
            p = int(r["pos"])
            z = r["z_score"] if pd.notna(r["z_score"]) else 0
            annot_positions[p] = f"pos {p}"

    # annotations: stagger y to avoid overlap based on x proximity
    sorted_annots = sorted(annot_positions.items())
    n_annots = len(sorted_annots)
    # assign alternating y levels, spacing by position proximity
    y_base, y_step = 1.10, 0.07
    y_levels = []
    for i in range(n_annots):
        # alternate high/low within groups of nearby positions
        y_levels.append(y_base + y_step * (i % 2))

    for i, (p, label) in enumerate(sorted_annots):
        row = df[df["pos"] == p]
        if row.empty:
            continue
        s       = float(row["score_wt"].fillna(row["score_uw"]).values[0])
        is_sig  = bool(row["sig"].values[0])
        txt_col = C_SIG if is_sig else "#333333"
        bar_x   = p + offsets[1]   # tip of weighted bar
        ax_mn.annotate(
            label,
            xy=(bar_x, min(s, 1.0)),
            xytext=(bar_x, y_levels[i]),
            fontsize=7, fontstyle="italic", color=txt_col,
            ha="center", va="bottom", clip_on=False,
            arrowprops=dict(arrowstyle="-", color="#bbbbbb", lw=0.7,
                            shrinkA=0, shrinkB=2),
            bbox=dict(boxstyle="round,pad=0.20", facecolor="white",
                      edgecolor=txt_col, alpha=0.97, linewidth=0.6),
            zorder=6
        )

    # legend
    leg_patches = [
        mpatches.Patch(facecolor=C_UW,  edgecolor="none", alpha=0.9, label="Unweighted"),
        mpatches.Patch(facecolor=c_wt,  edgecolor="none",             label=wt_label),
        *([] if wt2_path is None else [
          mpatches.Patch(facecolor=c_wt2 or C_CLADE, edgecolor="none", label="Clade-weighted")]),
        mpatches.Patch(facecolor=C_SIG, edgecolor="none",             label="Significant  z ≥ 2"),
    ]
    ax_mn.legend(handles=leg_patches, loc="upper left", fontsize=8,
                 framealpha=0.95, edgecolor="#dddddd", borderpad=0.5,
                 handlelength=1.1, ncol=3)

    # ── z-score strip ─────────────────────────────────────────────────────────
    ax_z.set_ylim(0, 1)
    ax_z.set_xlim(xmin, xmax)
    ax_z.spines[["top","right","left"]].set_visible(False)
    ax_z.tick_params(axis="y", left=False, labelleft=False)
    ax_z.set_xlabel(xlabel, fontsize=10)

    # diverging colormap centred at 0: blue = conserved between groups, red = specific
    z_vals = df["z_score"].dropna()
    z_abs  = max(abs(z_vals.min()), abs(z_vals.max()))
    cmap   = plt.cm.RdBu_r
    norm   = Normalize(vmin=-z_abs, vmax=z_abs)
    # pcolormesh for pixel-perfect alignment with main panel bars
    pos_min = int(df["pos"].min())
    pos_max = int(df["pos"].max())
    z_arr   = np.full(pos_max - pos_min + 1, np.nan)
    for _, row in df.iterrows():
        if pd.notna(row["z_score"]):
            z_arr[int(row["pos"]) - pos_min] = float(row["z_score"])
    x_edges  = np.arange(pos_min - 0.5, pos_max + 1.5)
    z_masked = np.ma.masked_invalid(z_arr).reshape(1, -1)
    ax_z.pcolormesh(x_edges, [0, 1], z_masked, cmap=cmap, norm=norm,
                    shading="flat", zorder=2)
    ax_z.set_facecolor("#cccccc")  # grey for NaN/masked cells

    # attach colorbar to main panel so it steals from ax_mn, not ax_z
    # (keeping ax_z the same width as ax_mn for alignment)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax_mn, orientation="vertical",
                        fraction=0.014, pad=0.01, aspect=20)
    cbar.set_label("z-score", fontsize=8)
    cbar.ax.tick_params(labelsize=7.5)
    cbar.set_ticks([-round(z_abs,1), 0, round(z_abs,1)])

    ax_z.tick_params(axis="x", labelsize=9)

    ax_hx.set_title(title, fontsize=11, fontweight="bold", pad=4, loc="left")

    for ext in ("pdf", "png"):
        out = OUT_DIR / f"{out_stem}.{ext}"
        fig.savefig(str(out), dpi=300 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
        print(f"Saved: {out.name}")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# Plot A — CENP-A vs H3
# ══════════════════════════════════════════════════════════════════════════════
ANNOT_A = {
    75:  "H3:Gln→CENPA:diverse",
    92:  "H3:Phe→CENPA:Trp",
    97:  "H3:Val→CENPA:Leu",
   110:  "H3:Gly→CENPA:diverse",
   115:  "H3:Thr→CENPA:Ala/Ser",
   127:  "H3:Ile→CENPA:Leu",
}

make_groupsim_plot(
    uw_path    = BASE / "split_entropy/groupsim/groupsim_gap08.txt",
    wt_path    = BASE / "split_entropy/groupsim_weighted/groupsim_cenpa_h3_clade_gap08.tsv",
    helix_path = BASE / "split_entropy/helix_positions_gap08.tsv",
    title      = "CENP-A (n=418) vs H3-like (n=901)",
    xlabel     = "Alignment position (trimmed, gap ≤ 80%)",
    out_stem   = "groupsim_cenpa_vs_h3_pub",
    c_wt       = C_WT_A,
    wt_label   = "Clade-weighted",
    annotations= ANNOT_A,
    helix_mode = "two_rows",
)

# ══════════════════════════════════════════════════════════════════════════════
# Plot B — Satellite vs Transposon CENP-A
# ══════════════════════════════════════════════════════════════════════════════
make_groupsim_plot(
    uw_path    = BASE / "split_entropy/groupsim_sat_trans/groupsim_sat_trans_gap085.txt",
    wt_path    = BASE / "split_entropy/groupsim_sat_trans_weighted/groupsim_st_clade_gap085.tsv",
    helix_path = BASE / "split_entropy/groupsim_sat_trans/helix_positions_gap085.txt",
    title      = "Satellite CENP-A (n=198) vs Transposon CENP-A (n=129)",
    xlabel     = "Alignment position (trimmed, gap ≤ 85%)",
    out_stem   = "groupsim_sat_vs_trans_pub",
    c_wt       = C_CLADE,
    wt_label   = "Clade-weighted",
    annotations= {},
    helix_mode = "one_row",
    top_n_auto = 6,
)
