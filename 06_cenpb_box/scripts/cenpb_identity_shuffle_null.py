#!/usr/bin/env python3
"""Composition null for 'identity to the canonical CENP-B motif'. Short (17 bp)
motifs match partly by base composition alone, so a raw %identity is hard to read.
For each species' box consensus we compare its identity to canonical against the
MEAN identity of many shuffles of that same consensus (same composition, sequence
order destroyed). excess = observed − shuffled tells us how much of the match is
real arrangement vs composition chance."""
import random, numpy as np, pandas as pd
from pathlib import Path
random.seed(0)
SAT=Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity")
CANON=[set("CT"),{"T"},{"T"},{"C"},{"G"},{"T"},{"T"},{"G"},{"G"},{"A"},{"A"},set("AG"),{"C"},{"G"},{"G"},{"G"},{"A"}]
def idto(seq): return 100*sum(seq[i] in CANON[i] for i in range(17))/17
K=3000

f=pd.read_csv(SAT/"figures/cenpb_flank_uncapped_per_species.tsv",sep="\t")
rows=[]
for _,r in f.iterrows():
    cons=r.box_consensus
    if not isinstance(cons,str) or len(cons)!=17 or (set(cons)-set("ACGT")): continue
    obs=idto(cons); ch=list(cons); sh=np.empty(K)
    for k in range(K): random.shuffle(ch); sh[k]=idto(ch)
    rows.append(dict(species=r.species, id_obs=round(obs,1), id_shuffle=round(sh.mean(),1),
                     id_shuffle_p95=round(np.percentile(sh,95),1),
                     excess_shuffle=round(obs-sh.mean(),1)))
df=pd.DataFrame(rows); df.to_csv(SAT/"figures/cenpb_identity_shuffle_null.tsv",sep="\t",index=False)

# reference: shuffle the canonical motif itself (composition chance for the box)
cn=list("CTTCGTTGGAAACGGGA"); base=np.mean([idto(random.sample(cn,17)) for _ in range(20000)])
print(f"canonical-composition chance identity (shuffle the box): {base:.1f}%")
print(f"per-species: median id_obs={df.id_obs.median():.1f}%  shuffle={df.id_shuffle.median():.1f}%  excess={df.excess_shuffle.median():.1f}%")
fk=df.merge(f[["species","name","vgroup","delta"]],on="species",how="left")
print("\nbirds + key:")
print(fk[fk.name.str.contains("Accipiter|Lagopus|Porphyrio|Corvus|Cervus",na=False)]
      [["name","id_obs","id_shuffle","excess_shuffle","delta"]].to_string(index=False))
