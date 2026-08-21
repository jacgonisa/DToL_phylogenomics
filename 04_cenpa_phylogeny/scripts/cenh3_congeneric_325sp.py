#!/usr/bin/env python3
"""Congeneric CENH3 test: for every pair of same-genus species, how identical
is their CENP-A? Congeneric pairs are shallow divergences, so this isolates
lineage-specific rate variation (the Geum question: Geum barely changed, yet
Viridiplantae is the fastest clade on average — who else is Geum-like vs fast?).

Genus is taken from the real binomial (cenpb name table); the tolID 3-letter
genus abbreviation collides (e.g. Macrophya vs Macropis both 'iyMac'), so
tolID prefix is only a *candidate* filter — confirmed by name, or by a sane
divergence time when a name is missing.
"""
import re
from pathlib import Path
from itertools import combinations
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from Bio import Phylo

HERE = Path(__file__).resolve().parent; MOD = HERE.parent; REPO = MOD.parent
PAIRS = MOD/"analysis/cenh3_halflife/cenh3_seqsim_pairs_325sp.tsv"
TREE  = REPO/"01_species_tree/data/full_325sp_calibrated_correlatedlambda01.nwk"
NAMES = REPO/"06_cenpb_box/data/cenpb_paper_motifs_per_species.tsv"
FIG   = MOD/"figures"; OUTD = MOD/"analysis/cenh3_halflife"
MAX_CONGENERIC_MYA = 45   # plausibility cap when a real name is missing

norm = lambda s: re.sub(r"[0-9.].*$", "", str(s).lower())

def main():
    t = Phylo.read(str(TREE), "newick")
    tolid = {norm(x.name): x.name.replace(".fa","") for x in t.get_terminals()}
    # candidate genus key from tolID: prefix(1-2 lc) + 3-letter genus block
    def gkey(code):
        m = re.match(r"^([a-z]{1,2}[A-Z][a-z]{2})[A-Z]", tolid.get(code,""))
        return m.group(1) if m else None
    nm = pd.read_csv(NAMES, sep="\t")
    name = {norm(c): n for c, n in zip(nm["species"], nm["name"])}
    genus = lambda code: (name[code].split()[0] if isinstance(name.get(code), str) else None)
    # propagate a genus name to every tolID-key that has one known member,
    # + a small manual fill for genera absent from the name table entirely
    key_genus = {}
    for code, tid in tolid.items():
        g = genus(code); k = gkey(code)
        if g and k: key_genus.setdefault(k, g)
    key_genus.setdefault("drMal", "Malus")   # drMalDome/Sylv: no binomial in name table

    p = pd.read_csv(PAIRS, sep="\t")
    p["gk"] = p.spA.map(gkey); p["gkB"] = p.spB.map(gkey)
    cand = p[(p.gk == p.gkB) & p.gk.notna()].copy()

    def classify(r):
        gA, gB = genus(r.spA), genus(r.spB)
        lbl = key_genus.get(r.gk, r.gk)
        if gA and gB:
            return ("same" if gA == gB else "diff"), (gA if gA==gB else f"{gA}/{gB}")
        # name missing on ≥1 side: trust tolID prefix only if divergence is plausible
        return ("same" if r.mya <= MAX_CONGENERIC_MYA else "diff"), lbl
    cand[["verdict","genus_lbl"]] = cand.apply(
        lambda r: pd.Series(classify(r)), axis=1)
    cong = cand[cand.verdict == "same"].copy().sort_values("pct_id")

    keep = ["genus_lbl","spA","spB","clade","mya","pct_id"]
    out = cong[keep].rename(columns={"genus_lbl":"genus"})
    out.to_csv(OUTD/"cenh3_congeneric_pairs.tsv", sep="\t", index=False)
    print(f"{len(out)} confirmed congeneric pairs "
          f"({cand[cand.verdict=='diff'].shape[0]} tolID-collisions dropped)\n")
    show = out.copy(); show["mya"] = show.mya.round(1); show["pct_id"] = show.pct_id.round(1)
    print(show.to_string(index=False))

    # ── figure: mya vs identity, labelled, coloured by clade ──────────────────
    pal = {"Vertebrates":"#1565C0","Invertebrates":"#E65100","Viridiplantae":"#2E7D32"}
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    for _, r in out.iterrows():
        ax.scatter(r.mya, r.pct_id, s=70, color=pal.get(r.clade,"#777"),
                   alpha=0.85, zorder=3, edgecolor="white", lw=0.8)
        g = r.genus.split()[0]
        hi = g == "Geum"
        ax.annotate(g, (r.mya, r.pct_id), xytext=(5, 3), textcoords="offset points",
                    fontsize=8, fontweight="bold" if hi else "normal",
                    color="black" if hi else "#444")
    ax.set_xlabel("Divergence time of the two species (My)")
    ax.set_ylabel("Best-copy CENP-A % identity")
    ax.set_title("Congeneric CENH3 conservation is strongly lineage-specific\n"
                 "same divergence age → very different identity (Geum unchanged; Juncus/Luzula fast)",
                 fontsize=10)
    ax.set_ylim(40, 102); ax.spines[["top","right"]].set_visible(False)
    ax.yaxis.grid(True, color="#f0f0f0")
    handles = [plt.Line2D([0],[0], marker='o', ls='', color=c, label=g)
               for g, c in pal.items() if g in set(out.clade)]
    ax.legend(handles=handles, fontsize=8, frameon=False, title="clade")
    plt.tight_layout()
    for ext in ("png","pdf"):
        fig.savefig(FIG/f"cenh3_congeneric_325sp.{ext}",
                    dpi=300 if ext=="png" else None, bbox_inches="tight", facecolor="white")
        print("Saved:", FIG/f"cenh3_congeneric_325sp.{ext}")

if __name__ == "__main__":
    main()
