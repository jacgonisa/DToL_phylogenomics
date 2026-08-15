#!/usr/bin/env python3
"""trait_pairplot_seaborn_325sp.py
Seaborn pairplot of per-species centromere/genome traits (dominant-array table from
trait_pgls_correlations_325sp.R), with BOTH naive Pearson r and phylogenetically
corrected r (PIC/PGLS) annotated in the upper triangle, plus a two-triangle
correlation heatmap (upper = naive Pearson, lower = phylo PIC)."""

import numpy as np, pandas as pd, seaborn as sns, matplotlib.pyplot as plt
from itertools import combinations

DATA = "/home/jg2070/Desktop/DToL_phylogenomics/01_species_tree/data"
FIG  = "/home/jg2070/Desktop/DToL_phylogenomics/01_species_tree/figures"
tab  = pd.read_csv(f"{DATA}/centromere_trait_table_325sp.tsv", sep="\t")
cor  = pd.read_csv(f"{DATA}/pgls_trait_correlations_325sp.tsv", sep="\t")

# log-transform the two skewed size traits for display + correlation consistency
tab["mono_len_bp"] = np.log10(tab["mono_len_bp"])
tab["genome_mb"]   = np.log10(tab["genome_mb"])
traits = ["regi","hor","mono_len_bp","sat_gc","genome_gc","genome_mb","chr_n"]
labels = {"regi":"HOR regim.","hor":"HOR score","mono_len_bp":"log10 monomer bp",
          "sat_gc":"sat GC%","genome_gc":"genome GC%","genome_mb":"log10 genome Mb",
          "chr_n":"chrom. number"}
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
g.add_legend(title="Clade", bbox_to_anchor=(1.01,0.5))
g.figure.suptitle("Centromere / genome trait pairplot (dominant array, n=173 species)\n"
                  "upper: naive Pearson vs phylogenetically corrected (PIC) r",
                  y=1.02, fontsize=12, fontweight="bold")
g.savefig(f"{FIG}/trait_pairplot_seaborn_325sp.png", dpi=200, bbox_inches="tight")
g.savefig(f"{FIG}/trait_pairplot_seaborn_325sp.pdf", bbox_inches="tight")
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
fig.savefig(f"{FIG}/trait_corr_heatmap_naive_vs_phylo_325sp.png", dpi=200, bbox_inches="tight")
fig.savefig(f"{FIG}/trait_corr_heatmap_naive_vs_phylo_325sp.pdf", bbox_inches="tight")
print("Wrote heatmap")
