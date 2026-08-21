#!/usr/bin/env python3
"""CENH3/CENP-A protein-similarity decay & half-life — the satellite
`seqsim` / Melters analysis (05_satellite_similarity) replicated for the
CENP-A protein instead of centromeric satellite DNA.

Difference from the satellite pipeline: CENP-A copies are already in one
protein MSA, so pairwise % identity is read straight off the alignment
(no all-vs-all BLAST). Everything downstream is identical in spirit:
  - one point per species pair -> node-averaged by MRCA age
  - divergence time = MRCA age in the calibrated chronogram (chronos
    correlated lambda=0.1), mya = tree.distance(a,b)/2
  - fit H = A*exp(-lam*t) + C  (empirical floor C free); t_half = ln2/lam
  - per clade (Vertebrates / Invertebrates / Viridiplantae)

The biological question: how fast does the CENP-A protein diverge, what is
its floor (conserved histone fold), and is the decay lineage-specific?

Outputs (figures/ + analysis/cenh3_halflife/):
  cenh3_halflife_325sp.{png,pdf}   3-panel decay + overlay
  cenh3_seqsim_pairs_325sp.tsv     per species-pair identities
  cenh3_halflife_325sp.tsv         per-clade half-life summary
"""
import re
from pathlib import Path
from itertools import combinations
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from Bio import Phylo

HERE  = Path(__file__).resolve().parent
MOD   = HERE.parent                      # 04_cenpa_phylogeny
REPO  = MOD.parent
MSA   = MOD/"data/cenpa430_H3_archaea10.aligned.clipkit.325sp.fasta"
GROUP = MOD/"data/groupsim/groups_gap085.txt"
TREE  = REPO/"01_species_tree/data/full_325sp_calibrated_correlatedlambda01.nwk"
TRAIT = REPO/"01_species_tree/data/centromere_trait_table_325sp_piotr.tsv"
SAT   = REPO/"05_satellite_similarity/data/seqsim_blastn_melters_325sp.tsv"  # satellite DNA pairs, for contrast
FIG   = MOD/"figures"; FIG.mkdir(exist_ok=True)
OUTD  = MOD/"analysis/cenh3_halflife"; OUTD.mkdir(parents=True, exist_ok=True)

MIN_OVERLAP = 30          # aligned non-gap columns required to score a copy pair
GAP = set("-.Xx")

# ── helpers ───────────────────────────────────────────────────────────────────
def read_fasta(p):
    seqs, cur, buf = {}, None, []
    for line in Path(p).read_text().splitlines():
        if line.startswith(">"):
            if cur: seqs[cur] = "".join(buf)
            cur, buf = line[1:].split()[0], []
        else: buf.append(line.strip())
    if cur: seqs[cur] = "".join(buf)
    return seqs

def parse_sp(sid):
    """species code from a CENP-A sequence id (mirrors copy-number script)."""
    if "_rescue" in sid: return sid.split("_rescue")[0]
    for sep in ("_GRCg7b_chr_","_chr_","_SUPER_"):
        if sep in sid: return sid.split(sep)[0]
    if "_ENA|" in sid: return sid.split("_ENA|")[0]
    return sid.split("_")[0]

# CENP-A ids that normalise to a code the tree/trait tables don't use
REMAP = {"roscan": "roscan_s"}      # Rosa canina: rosCan_subseq1_* -> tip rosCan_S27_v1

def norm(s):                        # unify code / tip / label: cut at first digit or dot
    n = re.sub(r"[0-9.].*$", "", str(s).lower())
    return REMAP.get(n, n)

def pct_id(a, b):
    """% identity over columns where BOTH aligned seqs have a residue."""
    m = n = 0
    for x, y in zip(a, b):
        if x in GAP or y in GAP: continue
        n += 1; m += (x == y)
    return (100.0*m/n, n) if n else (np.nan, 0)

def pair_identity(seqsA, seqsB):
    """species-pair identity = best copy-pair identity (the orthologous CENH3
    pair). Max over all cross-species copy pairs, so extra divergent paralogs
    don't drag the score down — we want how fast *the* CENP-A protein diverges.
    # ponytail: best-pair; switch to mean-of-best-per-copy if you want the
    #           paralog-inclusive Melters-style score instead."""
    best = [p for a in seqsA for b in seqsB
            for p, n in [pct_id(a, b)] if n >= MIN_OVERLAP]
    return max(best) if best else np.nan

def exp_free(t, A, lam, C): return A*np.exp(-lam*t) + C

def fit_decay(mya, sim):
    """fit H=A·e^(−λt)+C. Returns (A, lam, C, half_life, rate0) where
    rate0 = A·λ = %-identity lost per My at divergence 0 (the tangent slope —
    a cleaner 'how fast does it diverge' number than half-life, which hides the
    amplitude). Half-life is only meaningful once the floor C is sampled."""
    C0 = float(np.quantile(sim, 0.10)); A0 = float(np.max(sim)) - C0
    popt, _ = curve_fit(exp_free, mya, sim, p0=[max(A0,1), 0.02, C0],
                        bounds=([0,1e-5,0],[100,5,100]), maxfev=20000)
    A, lam, C = popt
    return A, lam, C, np.log(2)/lam, A*lam

# ── load CENP-A alignment, grouped by species ─────────────────────────────────
def main():
    grp = GROUP.read_text()
    cenpa_ids = [x for line in grp.splitlines() if line.startswith("CENPA:")
                 for x in line[6:].split(",") if x]
    aln = read_fasta(MSA)
    sp_seqs = {}
    for sid in cenpa_ids:
        if sid not in aln: continue
        sp_seqs.setdefault(norm(parse_sp(sid)), []).append(aln[sid])
    print(f"{len(cenpa_ids)} CENP-A ids -> {len(sp_seqs)} species")

    tr = Phylo.read(str(TREE), "newick")
    for c in tr.get_terminals(): c.name = c.name.replace(".fa", "")
    tip = {norm(t.name): t.name for t in tr.get_terminals()}
    trait = pd.read_csv(TRAIT, sep="\t")
    clade_raw = {norm(l): c for l, c in zip(trait["label"], trait["clade"])}
    CMAP = {"Vertebrata":"Vertebrates", "Invertebrata":"Invertebrates",
            "Viridiplantae":"Viridiplantae", "Fungi":"Fungi", "Protist":"Protist"}
    clade = {k: CMAP.get(v, v) for k, v in clade_raw.items()}

    species = [s for s in sp_seqs if s in tip and s in clade]
    print(f"{len(species)} species with tree tip + clade")

    # ── all-vs-all species pairs -> identity + MRCA age ───────────────────────
    rows = []
    for a, b in combinations(species, 2):
        if clade[a] != clade[b]: continue          # within-clade only (as satellites)
        try: mya = tr.distance(tip[a], tip[b]) / 2
        except Exception: continue
        if not mya or mya <= 0: continue
        pid = pair_identity(sp_seqs[a], sp_seqs[b])
        if np.isnan(pid): continue
        rows.append(dict(spA=a, spB=b, mya=mya, pct_id=pid, clade=clade[a]))
    df = pd.DataFrame(rows)
    df.to_csv(OUTD/"cenh3_seqsim_pairs_325sp.tsv", sep="\t", index=False)
    print(f"{len(df)} within-clade species pairs scored")

    # ── node-average (one point per MRCA) ─────────────────────────────────────
    df["mya"] = df["mya"].round(3)
    navg = df.groupby(["mya","clade"]).agg(
        sim=("pct_id","mean"), n=("pct_id","count")).reset_index()

    # ── fit + plot (3 panels + overlay), mirroring the satellite figure ───────
    pal = {"Vertebrates":"#1565C0","Invertebrates":"#E65100","Viridiplantae":"#2E7D32"}
    groups = ["Vertebrates","Invertebrates","Viridiplantae"]
    fig, axes = plt.subplots(1, 4, figsize=(20, 4.6))
    summary = []
    for col, grp_ in enumerate(groups):
        cp = pal[grp_]; ax = axes[col]
        d = navg[navg["clade"]==grp_].dropna(subset=["sim"])
        ax.scatter(d["mya"], d["sim"], s=d["n"]*0.8, color=cp, alpha=0.55, zorder=3)
        hl = C = rate0 = np.nan
        try:
            A, lam, C, hl, rate0 = fit_decay(d["mya"].values, d["sim"].values)
            tmax = d["mya"].max()
            # floor is only identified if the plateau is actually sampled; if the
            # half-life exceeds the deepest node, C rails and t½ is an extrapolation.
            constrained = hl <= tmax and C > 1
            t = np.linspace(0, tmax, 400)
            ax.plot(t, exp_free(t,A,lam,C), color=cp, lw=1.8, ls="-" if constrained else "--")
            ax.axhline(C, color=cp, lw=0.7, ls=":", alpha=0.6)
            lab_txt = ((f"t½ = {hl:.0f} My\nfloor = {C:.0f}%\n" if constrained
                        else f"t½ > {tmax:.0f} My (floor not reached)\n")
                       + f"initial rate = {rate0:.2f} %/My")
            ax.text(0.97,0.05, lab_txt, transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=10, color=cp, fontweight="bold")
            leg = (f"{grp_} (t½={hl:.0f} My)" if constrained
                   else f"{grp_} (still declining)")
            axes[3].plot(t, exp_free(t,A,lam,C), color=cp, lw=2,
                         ls="-" if constrained else "--", label=leg)
        except Exception as e:
            print(grp_, "fit failed", e)
        axes[3].scatter(d["mya"], d["sim"], s=d["n"]*0.5, color=cp, alpha=0.25, zorder=2)
        ax.set_title(grp_, fontsize=11, fontweight="bold", color=cp)
        ax.set_ylabel("Mean CENP-A % identity" if col==0 else "", fontsize=10)
        ax.set_xlabel("Divergence time (My)", fontsize=10)
        ax.set_ylim(40, 102)
        ax.spines[["top","right"]].set_visible(False); ax.yaxis.grid(True, color="#f0f0f0")
        summary.append(dict(clade=grp_, half_life_My=round(hl,1), floor_pct=round(C,1),
                            init_rate_pct_per_My=round(rate0,3),
                            floor_reached=bool(hl <= d["mya"].max() and C > 1),
                            n_nodes=len(d), n_pairs=int(df[df.clade==grp_].shape[0])))
    axes[3].set_title("All lineages", fontsize=11, fontweight="bold")
    axes[3].set_xlabel("Divergence time (My)", fontsize=10); axes[3].set_ylim(40,102)
    axes[3].legend(fontsize=8, frameon=False, loc="upper right")
    axes[3].spines[["top","right"]].set_visible(False); axes[3].yaxis.grid(True, color="#f0f0f0")
    fig.suptitle("CENP-A / CENH3 protein-similarity decay — all-vs-all, node-averaged (325-sp)\n"
                 "MRCA ages: chronos correlated λ=0.1.  Fit H = A·e^(−λt)+C (floor free).",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    for ext in ("png","pdf"):
        fig.savefig(FIG/f"cenh3_halflife_325sp.{ext}",
                    dpi=300 if ext=="png" else None, bbox_inches="tight", facecolor="white")
        print("Saved:", FIG/f"cenh3_halflife_325sp.{ext}")
    sdf = pd.DataFrame(summary)
    sdf.to_csv(OUTD/"cenh3_halflife_325sp.tsv", sep="\t", index=False)
    print(sdf.to_string(index=False))

    build_comparison(navg)

def _node_avg(df, valcol):
    df = df.copy(); df["mya"] = df["mya"].round(3)
    return df.groupby(["mya","clade"]).agg(sim=(valcol,"mean"),
                                           n=(valcol,"count")).reset_index()

def build_comparison(cenh3_navg):
    """CENP-A protein vs centromeric-satellite DNA: overlaid decay curves +
    a %/My initial-rate bar chart. Answers 'how fast, and is it lineage-specific,
    compared with the satellite it binds'."""
    pal = {"Vertebrates":"#1565C0","Invertebrates":"#E65100","Viridiplantae":"#2E7D32"}
    groups = ["Vertebrates","Invertebrates","Viridiplantae"]
    if not SAT.exists():
        print("satellite tsv absent — skipping comparison"); return
    sat = pd.read_csv(SAT, sep="\t").rename(columns={"group":"clade"})
    sat = sat[sat["clade"].isin(groups)]
    sat_navg = _node_avg(sat, "mean_pct_id")

    fig, (axC, axB) = plt.subplots(1, 2, figsize=(13, 4.8),
                                   gridspec_kw={"width_ratios":[1.6,1]})
    rate_rows = []
    for src, navg, ls, lab in [("CENP-A protein", cenh3_navg, "-", "CENP-A"),
                               ("satellite DNA", sat_navg, "--", "satellite")]:
        for grp_ in groups:
            d = navg[navg["clade"]==grp_].dropna(subset=["sim"])
            if len(d) < 4: continue
            try: A, lam, C, hl, r0 = fit_decay(d["mya"].values, d["sim"].values)
            except Exception as e: print(src, grp_, "fit failed", e); continue
            t = np.linspace(0, d["mya"].max(), 400)
            axC.plot(t, exp_free(t,A,lam,C), color=pal[grp_], lw=2, ls=ls,
                     label=f"{grp_} — {lab}")
            rate_rows.append(dict(source=src, clade=grp_, init_rate=r0, floor=C))
    axC.set_xlabel("Divergence time (My)"); axC.set_ylabel("Mean % identity")
    axC.set_title("Decay: CENP-A protein (solid) vs satellite DNA (dashed)", fontsize=10)
    axC.set_xlim(0, 120)   # zoom the fast early regime where the contrast lives
    axC.set_ylim(20, 100); axC.spines[["top","right"]].set_visible(False)
    axC.legend(fontsize=7, frameon=False, ncol=2)

    rr = pd.DataFrame(rate_rows)
    x = np.arange(len(groups)); w = 0.38
    for i, src in enumerate(["CENP-A protein","satellite DNA"]):
        vals = [rr[(rr.source==src)&(rr.clade==g)]["init_rate"].squeeze() if
                not rr[(rr.source==src)&(rr.clade==g)].empty else 0 for g in groups]
        axB.bar(x + (i-0.5)*w, vals, w, label=src,
                color=[pal[g] for g in groups], alpha=0.55 if i else 1.0,
                hatch="" if i==0 else "//", edgecolor="white")
    axB.set_xticks(x); axB.set_xticklabels(groups, rotation=20, ha="right", fontsize=8)
    axB.set_ylabel("Initial rate  (% identity lost / My)")
    axB.set_title("How fast identity erodes (tangent at t=0)\nsolid = CENP-A, hatched = satellite", fontsize=9)
    axB.spines[["top","right"]].set_visible(False)
    fig.suptitle("CENP-A protein diverges far slower than its satellite DNA — and the rate is lineage-specific",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    for ext in ("png","pdf"):
        fig.savefig(FIG/f"cenh3_vs_satellite_halflife_325sp.{ext}",
                    dpi=300 if ext=="png" else None, bbox_inches="tight", facecolor="white")
        print("Saved:", FIG/f"cenh3_vs_satellite_halflife_325sp.{ext}")
    rr.to_csv(OUTD/"cenh3_vs_satellite_rates_325sp.tsv", sep="\t", index=False)
    print(rr.round(3).to_string(index=False))

def _selftest():
    core = "MARTKQTARKSTGGKAPRKQLATKAARKSAP"    # 31 aa, > MIN_OVERLAP
    a = core[:5] + "--" + core[5:]
    b = core[:5] + "GG" + core[5:]              # 31 shared non-gap cols, identical -> 100%
    p, n = pct_id(a, b); assert n == 31 and abs(p-100) < 1e-9, (p, n)
    c = "W" + a[1:]                             # one mismatch over 31 -> 30/31
    p, n = pct_id(a, c); assert n == 31 and abs(p-3000/31) < 1e-9, (p, n)
    # pair_identity best-per-copy: identical copy present -> 100 despite a junk copy
    assert abs(pair_identity([a],[b, "G"*len(a)]) - 100) < 1e-9
    print("selftest ok")

if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv: _selftest()
    else: main()
