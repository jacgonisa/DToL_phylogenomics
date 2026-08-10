#!/usr/bin/env python3
"""Songbird 'labile' CENP-B-like box method (Formenti et al., Cell 2026, zebra
finch T2T; Suppl. Fig. 15). Unlike the strict Fachinetti substitution scan and
the PWM log-odds scan, the labile method searches the canonical 17-bp human
CENP-B box by EDIT DISTANCE (Levenshtein: indels allowed, not just
substitutions), with the paper's thresholds:
    edit distance <=5  for cross-species  (all non-human DToL)
    edit distance <=2  for the same species (human HG002 benchmark)
Candidates are called against a mononucleotide-shuffled null (the paper's
+/-5 flanking negative control), enrichment = obs/null.

The per-species edit-distance sweep (e<=0..5, both strands, shuffled null) was
already computed by cenpb_mismatch_titration_allspecies.py (regex {e<=k} == edit
distance), so this just reads that table and applies the songbird thresholds."""
from pathlib import Path
import pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SAT = Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity")
tit = pd.read_csv(SAT/"figures/cenpb_mismatch_titration_allspecies.tsv", sep="\t")

CROSS_E, SAME_E = 5, 2                       # songbird thresholds
def at(sp_df, e):                            # obs/null/pct at a given edit distance
    r = sp_df[sp_df.mm == e]
    return (r.obs.iloc[0], r.null.iloc[0], r.pct_obs.iloc[0]) if len(r) else (0, 0, 0.0)

rows = []
for sp, d in tit.groupby("species"):
    clade = d.clade.iloc[0]; n = int(d.n.iloc[0])
    e = SAME_E if sp == "human" else CROSS_E
    obs, null, pct = at(d, e)
    rows.append(dict(species=sp, clade=clade, n_monomers=n, edit_threshold=e,
                     labile_obs=int(obs), labile_null=int(null),
                     pct_labile=round(pct, 2),
                     enrichment=round(obs/null, 2) if null else (float("inf") if obs else 0.0)))
res = pd.DataFrame(rows).sort_values(["enrichment", "labile_obs"], ascending=False)
res.to_csv(SAT/"figures/cenpb_songbird_labile_per_species.tsv", sep="\t", index=False)

# candidate = labile box clearly above its shuffled null (obs>=null+3 and >=2x)
cand = res[(res.species != "human") & (res.labile_obs >= res.labile_null + 3)
           & (res.labile_obs >= 2*res.labile_null.clip(lower=1))]
print("=== songbird labile method (edit distance <=5 cross-species, <=2 human) ===")
print(res[["species","clade","edit_threshold","labile_obs","labile_null","enrichment"]]
      .head(15).to_string(index=False))
print(f"\ncross-species candidates above null (obs>=null+3 & >=2x null): {len(cand)}")
print(cand[["species","clade","labile_obs","labile_null","enrichment"]].to_string(index=False) or "  none")

# ---- plot: edit-distance sweep, both strands, human benchmark, thresholds marked ----
ccol = {"Vertebrates":"#1565C0","Invertebrate":"#EF6C00","Viridiplantae":"#2E7D32",
        "Fungi":"#6A1B9A","Protist":"#C62828"}
fig, ax = plt.subplots(figsize=(7.6, 5.6))
for sp, d in tit.groupby("species"):
    if d.clade.iloc[0] == "Human": continue
    d = d.sort_values("mm")
    ax.plot(d.mm, d.pct_obs, color=ccol.get(d.clade.iloc[0], "grey"), alpha=0.35, lw=0.8)
nu = tit[tit.species != "human"].groupby("mm").pct_null.mean().reset_index()
ax.plot(nu.mm, nu.pct_null, color="black", ls="--", lw=1.8, label="shuffled null (mean)")
h = tit[tit.species == "human"].sort_values("mm")
if len(h):
    ax.plot(h.mm, h.pct_obs, color="#d62728", lw=3, marker="o", ms=5,
            label="Human α-sat (HG002 benchmark)", zorder=5)
for x, lab in [(SAME_E, "same-species (e≤2)"), (CROSS_E, "cross-species (e≤5)")]:
    ax.axvline(x, color="grey", ls=":", lw=1)
    ax.text(x, ax.get_ylim()[1]*0.97, lab, rotation=90, va="top", ha="right", fontsize=8, color="grey")
ax.set_xlabel("edit distance to the 17-bp CENP-B box (Levenshtein, indels allowed)")
ax.set_ylabel("% monomers with a CENP-B-like box")
ax.set_title("Songbird 'labile' CENP-B-like box method\n(Formenti et al., Cell 2026 — edit-distance search, both strands)",
             fontweight="bold", fontsize=11)
ax.set_xticks(range(6)); ax.spines[["top","right"]].set_visible(False)
handles = [Line2D([0],[0],color="#d62728",lw=3,label="Human α-sat (HG002 benchmark)")] + \
          [Line2D([0],[0],color=c,lw=2,label=g) for g,c in ccol.items() if g in tit.clade.values] + \
          [Line2D([0],[0],color="black",ls="--",lw=1.8,label="shuffled null")]
ax.legend(handles=handles, fontsize=8.5, frameon=False)
plt.tight_layout()
for ext in ("png","pdf"):
    fig.savefig(SAT/f"figures/cenpb_songbird_labile.{ext}", dpi=300 if ext=="png" else None,
                bbox_inches="tight", facecolor="white")
    print("Saved:", SAT/f"figures/cenpb_songbird_labile.{ext}")
