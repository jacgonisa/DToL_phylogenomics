#!/usr/bin/env python3
"""Identity of a species' box-like windows to the canonical CENP-B motif, computed
PER WINDOW (not on the consensus). The consensus version is inflated: for a
high-copy satellite whose windows deviate at random positions, the deviations
cancel on averaging and the consensus returns to canonical (~100%). Per-window
identity avoids this — a random-scatter satellite sits near the ≤5-substitution
floor (~70%), a genuinely box-like satellite sits higher.

Null: a 17-bp motif matches partly by composition, so for each window we also
shuffle its bases (composition preserved, order destroyed) and score identity.
excess = mean(observed per-window identity) − mean(shuffled per-window identity).
Subsampled to CAP arrays/species."""
import regex, re, collections, numpy as np, pandas as pd, random
from pathlib import Path
random.seed(0)
BASE=Path("/home/jg2070/Desktop/dtol_review_August")
SAT =BASE/"DToL_phylogenomics_publication_325genomes/05_satellite_similarity"
ALL =BASE/"2026_trees/all.satellites.txt"
PAT =regex.compile("(?:[CT]TTCGTTGGAA[AG]CGGGA){s<=5}")
PRE =regex.compile("(?:TTCGTTGGAA){s<=3}")
comp=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
seqre=re.compile(r'"([ACGTNacgtn]{25,})"')
CAN=[set("CT"),{"T"},{"T"},{"C"},{"G"},{"T"},{"T"},{"G"},{"G"},{"A"},{"A"},set("AG"),{"C"},{"G"},{"G"},{"G"},{"A"}]
def subs(w): return sum(w[i] not in CAN[i] for i in range(17))     # substitutions vs canonical IUPAC
ALLOW=set(l.strip() for l in open(SAT/"species_325.txt") if l.strip() and not l.startswith("#"))
CAP=3000

acc=collections.defaultdict(lambda:[0,0,0])                        # [N, sum_obs_subs, sum_shuf_subs]
narr=collections.Counter(); n=0
for i,line in enumerate(open(ALL)):
    if i==0 or '"' not in line: continue
    sp=line.rstrip().rsplit('"',2)[-2].split(".")[0].lower()
    if sp not in ALLOW or narr[sp]>=CAP: continue
    m=seqre.search(line)
    if not m: continue
    s=m.group(1).upper(); narr[sp]+=1; n+=1
    if not PRE.search(s+"#"+rc(s)): continue
    a=acc[sp]
    for st in (s,rc(s)):
        for mm in PAT.finditer(st,overlapped=False):
            x,y=mm.start(),mm.end()
            if y-x!=17: continue
            w=st[x:y]
            if set(w)-set("ACGT"): continue
            a[0]+=1; a[1]+=subs(w)
            ch=list(w); random.shuffle(ch); a[2]+=subs(ch)        # per-window composition shuffle
    if n%2000000==0: print(f"  {n:,} arrays...",flush=True)
print(f"streamed {n:,} capped arrays",flush=True)

rows=[]
for sp,(N,so,ss) in acc.items():
    if N<10: continue
    rows.append(dict(species=sp, n_windows=N,
        id_obs=round(100*(17-so/N)/17,1),       # mean per-window identity to canonical
        id_shuffle=round(100*(17-ss/N)/17,1),   # mean per-window identity of shuffled windows
        excess=round(100*(ss-so)/(17*N),1)))
df=pd.DataFrame(rows)
f=pd.read_csv(SAT/"figures/cenpb_flank_uncapped_per_species.tsv",sep="\t")[["species","name","vgroup","clade"]]
fk=df.merge(f,on="species",how="left").reset_index(drop=True)
fk.to_csv(SAT/"figures/cenpb_identity_shuffle_null.tsv",sep="\t",index=False)
print("median per-window id_obs=%.1f  shuffle=%.1f  excess=%.1f | max id_obs=%.1f"%(
    fk.id_obs.median(),fk.id_shuffle.median(),fk.excess.median(),fk.id_obs.max()))
print(fk[fk.name.str.contains("Accipiter|Lagopus|Porphyrio|Triglochin|Berberis",na=False)]
      [["name","n_windows","id_obs","id_shuffle"]].to_string(index=False))

# ---- slopegraph (jittered): each species linked to its own per-window shuffle null ----
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
  "font.size":8,"axes.linewidth":0.7,"xtick.major.width":0.7,"ytick.major.width":0.7,
  "axes.labelsize":8,"legend.fontsize":6.6,"xtick.labelsize":7.5,"ytick.labelsize":7.5,"savefig.dpi":600})
ccol={"Vertebrates":"#0072B2","Invertebrate":"#E69F00","Viridiplantae":"#009E73","Fungi":"#CC79A7","Protist":"#D55E00"}
rng=np.random.default_rng(0); jx=rng.uniform(-0.075,0.075,len(fk)); xs=0+jx; xo=1+jx
fig,ax=plt.subplots(figsize=(5.0,4.4))
for i in range(len(fk)):
    r=fk.iloc[i]; ax.plot([xs[i],xo[i]],[r.id_shuffle,r.id_obs],color=ccol.get(r.clade,"grey"),alpha=0.16,lw=0.4,zorder=2)
ax.scatter(xs,fk.id_shuffle,s=10,color="0.55",edgecolor="none",alpha=0.6,zorder=3)
for cl in [c for c in ccol if c in set(fk.clade)]:
    m=(fk.clade==cl).values
    ax.scatter(xo[m],fk.id_obs[m],s=15,color=ccol[cl],edgecolor="0.3",lw=0.2,alpha=0.85,zorder=4,label=cl)
for xx,col in [(0,"id_shuffle"),(1,"id_obs")]:
    ax.plot([xx-0.13,xx+0.13],[fk[col].median()]*2,color="k",lw=2,zorder=6)
ax.set_xticks([0,1]); ax.set_xticklabels(["shuffled windows\n(composition chance)","observed windows\n(mean per-window)"])
ax.set_xlim(-0.4,1.4); ax.set_ylim(0,103)
ax.set_ylabel("mean per-window identity to canonical (%)")
ax.set_title("Per-window identity vs shuffle-the-window null",fontweight="bold",fontsize=9)
ax.legend(frameon=False,loc="lower right",fontsize=6.3); ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
for ext in ("png","pdf"): fig.savefig(SAT/f"figures/cenpb_identity_shuffle_null.{ext}",dpi=600 if ext=="png" else None,bbox_inches="tight",facecolor="white")
print("Saved figures/cenpb_identity_shuffle_null.png/pdf")
