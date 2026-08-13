#!/usr/bin/env python3
"""
1. Bar chart of CENPA copy number per species (n_loci_pass40pct), coloured by centromere type
2. Boxplot + Mann-Whitney test: monotypic vs ditypic satellite species
3. Gene duplication calls from the CENPA gene tree (same-species paralog pairs)
"""
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from Bio import Phylo

BASE  = "/home/jg2070/Desktop/dtol_review_August"
PUB   = f"{BASE}/DToL_phylogenomics_publication_325genomes"
TREE_DIR = f"{PUB}/04_cenpa_phylogeny"

cn  = pd.read_csv(f"{TREE_DIR}/cenpa_copy_number_table.tsv", sep="\t")
sat = pd.read_csv(f"{BASE}/2026_trees/annotation_centromeres/organized/plots/"
                  "centromere_length/genome_vs_centromere_length_satellite_species_table.tsv",
                  sep="\t")

# ── normalise keys for join ───────────────────────────────────────────────────
def norm(s):
    s = str(s).lower()
    s = re.sub(r'\.(grcg7b|hap\d+(\.\d+)?)$', '', s)
    s = re.sub(r'\.[\d]+$', '', s)
    return s

cn["_key"]  = cn["species_code_full"].map(norm)
sat["_key"] = sat["species_id"].map(norm)
df = cn.merge(sat[["_key","n_sat_families"]], on="_key", how="left")

# working copy number column
df["copies"] = pd.to_numeric(df["n_loci_pass40pct"], errors="coerce").fillna(0).astype(int)

# centromere type palette
pal = {"Satellite":"#ff2d87","Transposon":"#4f84e8","Holocentric":"#2d7d32",
       "Satellite/transposon":"#8338ec","Unknown":"#aaaaaa","Monocentric sequence unknown":"#cccccc"}
df["ct"] = df["centromere_type"].fillna("Unknown")

# ── Figure 1: copy-number bar chart ──────────────────────────────────────────
plot_df = (df[df["copies"] > 0]
           .sort_values(["copies","ct"], ascending=[False,True])
           .reset_index(drop=True))

fig, ax = plt.subplots(figsize=(18, 5))
for i, row in plot_df.iterrows():
    ax.bar(i, row["copies"], color=pal.get(row["ct"],"#aaaaaa"),
           width=0.85, linewidth=0)

ax.set_xticks(range(len(plot_df)))
ax.set_xticklabels(plot_df["species_code_full"], rotation=90, fontsize=4.5)
ax.set_ylabel("CENPA copy number\n(unique loci, ≥40% ref identity)", fontsize=10)
ax.set_title("CENPA copy number per species — 325sp dataset (curated gene tree)", fontsize=11, fontweight="bold")
ax.spines[["top","right"]].set_visible(False)
ax.yaxis.grid(True, color="#eeeeee", zorder=0)

patches = [plt.matplotlib.patches.Patch(color=c, label=l) for l,c in pal.items()
           if l in df["ct"].unique()]
ax.legend(handles=patches, fontsize=8, loc="upper right",
          framealpha=0.9, edgecolor="#cccccc")

plt.tight_layout()
fig.savefig(f"{TREE_DIR}/figures/cenpa_copy_number_bar.pdf", bbox_inches="tight")
fig.savefig(f"{TREE_DIR}/figures/cenpa_copy_number_bar.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print("Saved: cenpa_copy_number_bar")

# ── Figure 2: ditypic vs monotypic test ──────────────────────────────────────
sat_only = df[(df["ct"] == "Satellite") & df["n_sat_families"].notna()].copy()
sat_only["group"] = np.where(sat_only["n_sat_families"] >= 2, "Ditypic\n(≥2 sat families)", "Monotypic\n(1 sat family)")
mono = sat_only.loc[sat_only["group"].str.startswith("Mono"), "copies"]
dity = sat_only.loc[sat_only["group"].str.startswith("Dity"), "copies"]

stat, pval = mannwhitneyu(dity, mono, alternative="greater")
print(f"\nMann-Whitney (ditypic > monotypic): U={stat:.0f}  p={pval:.4f}")
print(f"  Monotypic n={len(mono)}  median={mono.median():.1f}  mean={mono.mean():.2f}")
print(f"  Ditypic   n={len(dity)}  median={dity.median():.1f}  mean={dity.mean():.2f}")

# also Wilcoxon rank-sum on all satellite species with n_sat_families known
print(f"\n  n_sat_families distribution among satellite species:")
print(sat_only["n_sat_families"].value_counts().sort_index().to_string())

fig2, axes = plt.subplots(1, 2, figsize=(11, 5),
                           gridspec_kw={"width_ratios":[1.6,1]})

# left: boxplot + jitter
ax = axes[0]
groups  = ["Monotypic\n(1 sat family)", "Ditypic\n(≥2 sat families)"]
data    = [mono.values, dity.values]
bp = ax.boxplot(data, positions=[0,1], widths=0.45, patch_artist=True,
                medianprops=dict(color="black", lw=1.8),
                boxprops=dict(facecolor="#f0f0f0"), showfliers=False)
for i, (vals, col) in enumerate(zip(data, ["#ff2d87","#c2185b"])):
    bp["boxes"][i].set_facecolor(col); bp["boxes"][i].set_alpha(0.35)
    jx = np.random.uniform(-0.18, 0.18, len(vals))
    ax.scatter(i + jx, vals, s=22, color=col, alpha=0.7, zorder=3)

sig = "***" if pval<0.001 else "**" if pval<0.01 else "*" if pval<0.05 else f"p={pval:.3f}"
ymax = max(max(mono), max(dity)) + 0.5
ax.plot([0,0,1,1],[ymax-0.2,ymax,ymax,ymax-0.2],"k-",lw=0.9)
ax.text(0.5, ymax+0.05, sig, ha="center", va="bottom", fontsize=11)
ax.set_xticks([0,1]); ax.set_xticklabels(groups, fontsize=10)
ax.set_ylabel("CENPA copy number (unique loci, ≥40%)", fontsize=10)
ax.set_title("CENPA copies: monotypic vs ditypic satellite species", fontsize=10, fontweight="bold")
ax.spines[["top","right"]].set_visible(False)
ax.yaxis.grid(True, color="#eeeeee", zorder=0)
ax.text(0.98, 0.02, f"Mann-Whitney U={stat:.0f}, p={pval:.4f}\n(one-sided: ditypic > monotypic)",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="#555555")

# right: copy number distribution by n_sat_families (scatter + mean)
ax2 = axes[1]
for fam, grp in sat_only.groupby("n_sat_families"):
    ax2.scatter([fam]*len(grp), grp["copies"],
                s=25, alpha=0.6, color="#ff2d87", zorder=2)
    ax2.plot(fam, grp["copies"].mean(), "D", color="#7b0000",
             ms=8, zorder=3)
ax2.set_xlabel("Number of satellite families", fontsize=10)
ax2.set_ylabel("CENPA copy number", fontsize=10)
ax2.set_title("Copies vs satellite family richness", fontsize=10, fontweight="bold")
ax2.set_xticks(sorted(sat_only["n_sat_families"].unique()))
ax2.spines[["top","right"]].set_visible(False)
ax2.yaxis.grid(True, color="#eeeeee", zorder=0)
ax2.plot([], [], "D", color="#7b0000", ms=7, label="mean")
ax2.legend(fontsize=8)

plt.tight_layout()
fig2.savefig(f"{TREE_DIR}/figures/cenpa_ditypic_test.pdf", bbox_inches="tight")
fig2.savefig(f"{TREE_DIR}/figures/cenpa_ditypic_test.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print("Saved: cenpa_ditypic_test")

# ── Gene duplications from tree ───────────────────────────────────────────────
print("\n── Gene duplications from CENPA tree ──")
cenpa_tips_set = set(open(f"{TREE_DIR}/cenpa_tips_from_contree.txt").read().split())

def sp_of(sid):
    if "_rescue" in sid: return sid.split("_rescue")[0]
    for sep in ("_GRCg7b_chr_","_chr_","_SUPER_","_ENA|","_scaffold","_tig","_ptg","_HiC","_LG","_NC_"):
        if sep in sid: return sid.split(sep)[0]
    m = re.match(r'^([a-zA-Z0-9._]+?)_[A-Z]', sid)
    return m.group(1) if m else sid.split("_")[0]

from collections import defaultdict
sp_to_seqs = defaultdict(list)
for t in cenpa_tips_set:
    sp_to_seqs[sp_of(t)].append(t)

dup_rows = [(sp, len(seqs), seqs) for sp, seqs in sp_to_seqs.items() if len(seqs) > 1]
dup_rows.sort(key=lambda x: -x[1])
print(f"Species with ≥2 CENPA copies in tree: {len(dup_rows)}")
for sp, n, seqs in dup_rows[:10]:
    print(f"  {sp}: {n} paralogs")

# write duplication table
dup_df = pd.DataFrame([
    {"species": sp, "n_paralogs": n,
     "seq_ids": "|".join(seqs)}
    for sp, n, seqs in dup_rows
])
dup_df.to_csv(f"{TREE_DIR}/cenpa_gene_duplications.tsv", sep="\t", index=False)
print(f"Saved: cenpa_gene_duplications.tsv ({len(dup_df)} species with duplications)")
