#!/usr/bin/env python3
"""Bar views of the exact-IUPAC CENP-B motif result, to accompany the
enrichment-vs-density scatter (cenpb_paper_motifs_plot.py):
  (1) enrichment over the dinucleotide null (obs/exp)   -> cenpb_bars_enrichment
  (2) motif density (hits per Mbp)                       -> cenpb_bars_density
Both use a broken y-axis (α-satellite is ~1000× the DToL clades). HSat canonical is
5 boundary hits, hatched + flagged."""
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
SAT=Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity")
bench=pd.read_csv(SAT/"figures/cenpb_human_benchmark.tsv",sep="\t")
per=pd.read_csv(SAT/"figures/cenpb_paper_motifs_per_clade.tsv",sep="\t"); n_sat=int(per.n_species.sum())
tiers=["canonical","broad","degenerated"]; tlab={"canonical":"canonical","broad":"broad","degenerated":"degenerate"}
tcol={"canonical":"#D55E00","broad":"#0072B2","degenerated":"#E69F00"}
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
  "font.size":8,"axes.linewidth":0.7,"xtick.major.width":0.7,"ytick.major.width":0.7,
  "xtick.major.size":3,"ytick.major.size":3,"axes.labelsize":8,"legend.fontsize":7,
  "xtick.labelsize":7.3,"ytick.labelsize":7.5,"savefig.dpi":600})
def val(col):
    def g(row_ctrl=None,clade=None):
        if row_ctrl is not None: return bench.loc[bench.control==row_ctrl,col].iloc[0]
        return per.loc[per.clade==clade,col].iloc[0]
    return g
labels=["α-satellite\n(HG002, +)","HSat 1/2/3\n(HG002, −)","Vertebrates","Invertebrates","Plants"]

def make(col_tpl, ylabel, fname, top_ylim, bot_ylim, nullline, unit):
    G=[]
    for lab,ctrl,cl in [(labels[0],"alpha-satellite (positive)",None),(labels[1],"HSat1/2/3 (negative)",None),
                        (labels[2],None,"Vertebrates"),(labels[3],None,"Invertebrate"),(labels[4],None,"Viridiplantae")]:
        d={}
        for t in tiers:
            c=col_tpl.format(t=t)
            d[t]= bench.loc[bench.control==ctrl,c].iloc[0] if ctrl else per.loc[per.clade==cl,c].iloc[0]
        G.append((lab,d))
    x=np.arange(len(G)); w=0.26
    fig,(top,bot)=plt.subplots(2,1,sharex=True,figsize=(5.9,4.1),gridspec_kw=dict(height_ratios=[1,1.5],hspace=0.07))
    for j,t in enumerate(tiers):
        vals=[g[1][t] for g in G]
        top.bar(x+(j-1)*w,vals,w,color=tcol[t],edgecolor="none",zorder=3)
        bot.bar(x+(j-1)*w,vals,w,color=tcol[t],label=tlab[t],edgecolor="none",zorder=3)
    top.set_ylim(*top_ylim); bot.set_ylim(*bot_ylim)
    # HSat canonical = 5 boundary hits -> hatch + flag
    top.bar(1-w,G[1][1]["canonical"],w,facecolor="white",edgecolor="#C62828",hatch="////",lw=0.9,zorder=4)
    top.annotate("HSat canonical:\n5 boundary hits",(1-w,G[1][1]["canonical"]),ha="center",va="bottom",fontsize=5.2,color="#C62828")
    for gi,(gn,gv) in enumerate(G):
        for j,t in enumerate(tiers):
            v=gv[t]; xx=gi+(j-1)*w
            if gi==1 and t=="canonical": continue
            if top_ylim[0]<v<top_ylim[1]: top.annotate(f"{v:,.0f}{unit}",(xx,v),ha="center",va="bottom",fontsize=5.2,color=tcol[t])
            elif 0<v<bot_ylim[1]: bot.annotate(f"{v:.2f}",(xx,v),ha="center",va="bottom",fontsize=5.2,color=tcol[t])
    if nullline: bot.axhline(1,ls=(0,(4,3)),color="0.4",lw=0.9); bot.text(1.55,1.05,"null",ha="left",fontsize=6.2,color="0.4")
    top.spines["bottom"].set_visible(False); bot.spines["top"].set_visible(False); top.tick_params(bottom=False)
    dm=dict(marker=[(-1,-0.5),(1,0.5)],markersize=7,linestyle="none",color="k",mec="k",mew=0.8,clip_on=False)
    top.plot([0,1],[0,0],transform=top.transAxes,**dm); bot.plot([0,1],[1,1],transform=bot.transAxes,**dm)
    for ax in (top,bot): ax.spines[["right"]].set_visible(False)
    bot.set_xticks(x); bot.set_xticklabels(labels)
    bot.set_ylabel(ylabel); bot.yaxis.set_label_coords(-0.10,0.72)
    bot.legend(frameon=False,handlelength=1.1,title="motif tier",title_fontsize=7,loc="upper center",ncol=3)
    top.set_title(f"Exact IUPAC CENP-B motif — {fname.split('_')[-1]} ({n_sat} species with satellites)",fontsize=8.4,fontweight="bold",pad=4)
    top.annotate(f"α-sat canonical {'off top ↑' if 'enrich' in fname else ''}",xy=(0-w,top_ylim[1]*0.95),ha="left",va="top",fontsize=5.8,color="#D55E00")
    if nullline: bot.annotate("canonical = 0 in all DToL clades",xy=(0.5,0.62),xycoords="axes fraction",fontsize=6,color="#D55E00",va="top",ha="center")
    plt.tight_layout()
    for ext in ("png","pdf"): fig.savefig(SAT/f"figures/{fname}.{ext}",dpi=600 if ext=="png" else None,bbox_inches="tight",facecolor="white")
    plt.close(fig); print("Saved",fname)

make("{t}_enrich_dinuc","enrichment over dinucleotide null (obs/exp)","cenpb_bars_enrichment",(3.6,16000),(0,3.35),True,"×")
make("{t}_perMbp","motif density (hits per Mbp)","cenpb_bars_density",(5,2300),(0,4.6),False,"")
