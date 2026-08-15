#!/usr/bin/env python3
"""trait_pairplot_seaborn_325sp.py
Seaborn pairplot of per-species centromere/genome traits (dominant-array table from
trait_pgls_correlations_325sp.R), with BOTH naive Pearson r and phylogenetically
corrected r (PIC/PGLS) annotated in the upper triangle, plus a two-triangle
correlation heatmap (upper = naive Pearson, lower = phylo PIC)."""

import os, numpy as np, pandas as pd, seaborn as sns, matplotlib.pyplot as plt
from itertools import combinations

AGG  = os.environ.get("AGG", "ian").lower()
DATA = "/home/jg2070/Desktop/DToL_phylogenomics/01_species_tree/data"
FIG  = "/home/jg2070/Desktop/DToL_phylogenomics/01_species_tree/figures"
tab  = pd.read_csv(f"{DATA}/centromere_trait_table_325sp_{AGG}.tsv", sep="\t")
cor  = pd.read_csv(f"{DATA}/pgls_trait_correlations_325sp_{AGG}.tsv", sep="\t")

# log-transform the two skewed size traits for display + correlation consistency
tab["mono_len_bp"] = np.log10(tab["mono_len_bp"])
tab["genome_mb"]   = np.log10(tab["genome_mb"])
traits = ["regi","hor","mono_len_bp","sat_gc","genome_gc","genome_mb","chr_n"]
labels = {"regi":"HOR regim.","hor":"HOR score","mono_len_bp":"monomer length (bp)",
          "sat_gc":"sat GC%","genome_gc":"genome GC%","genome_mb":"genome size (Mb)",
          "chr_n":"chrom. number"}
# these two are plotted as log10; show real units on the axes
logticks = {"mono_len_bp":[30,100,300,1000], "genome_mb":[200,1000,5000]}
pal = {"Vertebrata":"#1f78b4","Invertebrata":"#e6550d","Viridiplantae":"#31a354",
       "Fungi":"#7a0177","Protist":"#525252","Other":"#999999"}

# lookup dicts for annotation (both orderings)
def mk(col_r, col_p):
    d = {}
    for _,r in cor.iterrows():
        d[(r.trait_x,r.trait_y)] = (r[col_r], r[col_p]); d[(r.trait_y,r.trait_x)] = (r[col_r], r[col_p])
    return d
naive = mk("pearson_r","pearson_p"); phylo = mk("phylo_r_pic","phylo_p_pic")
def star(p): return "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<0.05 else "ns"

d = tab[["clade"]+traits].dropna(subset=traits, how="all")
g = sns.PairGrid(d, vars=traits, hue="clade", palette=pal, diag_sharey=False, corner=False)
g.map_lower(sns.scatterplot, s=14, alpha=0.7, edgecolor="none")
g.map_diag(sns.kdeplot, fill=True, warn_singular=False, color="grey", lw=1)

def annot(x, y, **kw):
    ax = plt.gca(); ax.set_axis_off()
    xi, yi = traits[kw["_c"]], traits[kw["_r"]]
    rn,pn = naive.get((xi,yi),(np.nan,np.nan)); rp,pp = phylo.get((xi,yi),(np.nan,np.nan))
    ax.text(0.5,0.62,f"Pearson r = {rn:+.2f} {star(pn)}",ha="center",va="center",
            fontsize=9,transform=ax.transAxes,color="#444444")
    ax.text(0.5,0.34,f"phylo r = {rp:+.2f} {star(pp)}",ha="center",va="center",
            fontsize=9,fontweight="bold",transform=ax.transAxes,
            color="#b2182b" if abs(rp)>=0.2 and pp<0.05 else "#888888")
n=len(traits)
for r in range(n):
    for c in range(n):
        if c>r:
            plt.sca(g.axes[r][c]); annot(None,None,_r=r,_c=c)
for i,t in enumerate(traits):
    g.axes[-1][i].set_xlabel(labels[t]); g.axes[i][0].set_ylabel(labels[t])
# relabel the log10-plotted axes (monomer length, genome size) with real units
for t,vals in logticks.items():
    i = traits.index(t); pos = [np.log10(v) for v in vals]
    g.axes[-1][i].set_xticks(pos); g.axes[-1][i].set_xticklabels([str(v) for v in vals])
    g.axes[i][0].set_yticks(pos); g.axes[i][0].set_yticklabels([str(v) for v in vals])
g.add_legend(title="Clade", bbox_to_anchor=(1.01,0.5))
agg_lab = "dominant array" if AGG=="ian" else "freq-weighted mean"
g.figure.suptitle(f"Centromere / genome trait pairplot ({agg_lab}, {len(d)} species)\n"
                  "upper: naive Pearson vs phylogenetically corrected (PIC) r",
                  y=1.02, fontsize=12, fontweight="bold")
g.savefig(f"{FIG}/trait_pairplot_seaborn_325sp_{AGG}.png", dpi=200, bbox_inches="tight")
g.savefig(f"{FIG}/trait_pairplot_seaborn_325sp_{AGG}.pdf", bbox_inches="tight")
print("Wrote pairplot")

# ── two-triangle correlation heatmap: upper naive, lower phylo ───────────────────
M = np.full((n,n), np.nan); S = np.empty((n,n), dtype=object)
for i in range(n): S[i,i]=""; M[i,i]=np.nan
for (i,a),(j,b) in combinations(enumerate(traits),2):
    rn,pn = naive[(a,b)]; rp,pp = phylo[(a,b)]
    M[i,j]=rn; S[i,j]=f"{rn:+.2f}\n{star(pn)}"           # upper = naive
    M[j,i]=rp; S[j,i]=f"{rp:+.2f}\n{star(pp)}"           # lower = phylo
fig,ax=plt.subplots(figsize=(8.5,7.5))
sns.heatmap(M, annot=S, fmt="", cmap="RdBu_r", center=0, vmin=-0.6, vmax=0.6,
            xticklabels=[labels[t] for t in traits], yticklabels=[labels[t] for t in traits],
            linewidths=1, linecolor="white", cbar_kws={"label":"correlation r"},
            annot_kws={"fontsize":8}, ax=ax)
ax.set_title("Trait correlations: upper = naive Pearson, lower = phylogenetic (PIC)\n"
             "*** p<0.001  ** p<0.01  * p<0.05", fontsize=11, fontweight="bold")
plt.xticks(rotation=40,ha="right"); plt.yticks(rotation=0); plt.tight_layout()
fig.savefig(f"{FIG}/trait_corr_heatmap_naive_vs_phylo_325sp_{AGG}.png", dpi=200, bbox_inches="tight")
fig.savefig(f"{FIG}/trait_corr_heatmap_naive_vs_phylo_325sp_{AGG}.pdf", bbox_inches="tight")
print("Wrote heatmap")

# ── partial vs pairwise (phylogenetic) for the regi/HOR/monomer trio ─────────────
trio = pd.read_csv(f"{DATA}/phylo_pairwise_vs_partial_trio_{AGG}.tsv", sep="\t")
xlab = [f"{p}\n(control: {c})" for p,c in zip(trio.pair, trio.control_for)]
xi = np.arange(len(trio)); w = 0.38
fig,ax = plt.subplots(figsize=(8,4.8))
ax.bar(xi-w/2, trio.pic_r_pairwise, w, label="pairwise (phylo, PIC)", color="#9ecae1", edgecolor="k", lw=.5)
ax.bar(xi+w/2, trio.pic_r_partial,  w, label="partial (control 3rd trait)", color="#fc9272", edgecolor="k", lw=.5)
for k,(rp,pp) in enumerate(zip(trio.pic_r_partial, trio.partial_p)):
    ax.text(xi[k]+w/2, rp+0.02*np.sign(rp+1e-9), star(pp), ha="center",
            va="bottom" if rp>=0 else "top", fontsize=11, fontweight="bold")
ax.axhline(0, color="grey", lw=.8); ax.set_xticks(xi); ax.set_xticklabels(xlab, fontsize=9)
ax.set_ylabel("phylogenetic correlation r"); ax.set_ylim(-0.15, max(0.75, trio.pic_r_pairwise.max()+0.1))
ax.set_title(f"Pairwise vs partial phylogenetic correlation — regi/HOR/monomer trio ({agg_lab})\n"
             "partial holds the 3rd trait constant; stars = PGLS partial-slope p", fontsize=10, fontweight="bold")
ax.legend(frameon=False, fontsize=9); plt.tight_layout()
fig.savefig(f"{FIG}/trait_partial_vs_pairwise_trio_325sp_{AGG}.png", dpi=200, bbox_inches="tight")
fig.savefig(f"{FIG}/trait_partial_vs_pairwise_trio_325sp_{AGG}.pdf", bbox_inches="tight")
print("Wrote partial-vs-pairwise trio panel")

# ── phylogenetic signal per trait (Pagel lambda) ────────────────────────────────
sig = pd.read_csv(f"{DATA}/trait_phylo_signal_{AGG}.tsv", sep="\t").sort_values("lambda")
col = ["#fc9272" if s=="weak" else "#3182bd" for s in sig.structured]
fig,ax = plt.subplots(figsize=(7,4.3))
ax.barh([labels[t] for t in sig.trait], sig["lambda"], color=col, edgecolor="k", lw=.5)
ax.axvline(0.5, color="grey", ls="--", lw=1)
for y,(l,k) in enumerate(zip(sig["lambda"], sig["K"])):
    ax.text(l+0.02, y, f"λ={l:.2f} (K={k:.2f})", va="center", fontsize=8)
ax.set_xlim(0,1.2); ax.set_xlabel("Pagel's λ  (0 = no phylogenetic signal, 1 = Brownian)")
ax.set_title(f"Phylogenetic signal per trait ({agg_lab})\n"
             "blue = tree-structured (λ≥0.5) → correlations can be shared-ancestry artifacts;\n"
             "red = weak signal → correlations are already phylogeny-independent",
             fontsize=9.5, fontweight="bold", loc="left")
plt.tight_layout()
fig.savefig(f"{FIG}/trait_phylo_signal_325sp_{AGG}.png", dpi=200, bbox_inches="tight")
fig.savefig(f"{FIG}/trait_phylo_signal_325sp_{AGG}.pdf", bbox_inches="tight")
print("Wrote phylo-signal barplot")
