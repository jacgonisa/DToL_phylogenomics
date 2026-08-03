#!/usr/bin/env python3
"""Validation figure for cenpb_box_finder: human alpha-sat (positive control) vs
horseshoe bat vs the rest of DToL, across the 4 lines of evidence (as % of records,
real vs shuffled null)."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D="/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity"
df=pd.read_csv(f"{D}/figures_sweep/cenpb_box_evidence_per_species.tsv",sep="\t")
METH=["canonical","broad","degenerate","pwm"]
for m in METH:
    df[m+"_pct"]=100*df[m]/df["n"]; df[m+"_npct"]=100*df[m+"_null"]/df["n"]

def row(sp): return df[df.species==sp].iloc[0]
human=row("human"); bat=row("mrhisin")
rest=df[~df.species.isin(["human","mrhisin"])]
restmax=rest[[m+"_pct" for m in METH]].max()

groups=[("Human α-sat\n(positive control)",human,"#1f77b4"),
        ("Horseshoe bat\n(Rhinolophus)",bat,"#F72485"),
        ("Other 174 DToL\n(max across spp)",None,"#999999")]
fig,ax=plt.subplots(figsize=(10,5.5))
x=np.arange(len(METH)); w=0.26
for gi,(lab,r,c) in enumerate(groups):
    vals=[ (r[m+"_pct"] if r is not None else restmax[m+"_pct"]) for m in METH]
    nvals=[ (r[m+"_npct"] if r is not None else 0) for m in METH]
    ax.bar(x+(gi-1)*w, vals, w, color=c, label=lab, zorder=3)
    ax.bar(x+(gi-1)*w, nvals, w, color="none", edgecolor="k", hatch="////", lw=0.5, zorder=4)
ax.set_xticks(x); ax.set_xticklabels(["canonical\nIUPAC","broad\nIUPAC","degenerate\nIUPAC","PWM\n(FIMO)"])
ax.set_ylabel("% of satellite records with a CENP-B box")
ax.set_title("CENP-B box finder — validation across lines of evidence\n(hatched = shuffled-sequence null)")
ax.legend(frameon=False); ax.spines[["top","right"]].set_visible(False)
fig.tight_layout(); fig.savefig(f"{D}/figures_sweep/cenpb_validation.png",dpi=200)
print("human:",{m:round(human[m+'_pct'],1) for m in METH})
print("bat  :",{m:round(bat[m+'_pct'],1) for m in METH})
print("restmax:",{m:round(restmax[m+'_pct'],2) for m in METH})
print("saved figures_sweep/cenpb_validation.png")
