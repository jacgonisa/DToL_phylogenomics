#!/usr/bin/env python3
"""Screen curated satellite families for the canonical CENP-B box (both strands),
allowing mismatches (boxes are diverged in non-model species), with a composition-
shuffled null. Uses the representative satellite set (<=500 seqs/species). Prediction:
strongest in mammals."""
import regex, re, collections, random
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
random.seed(0)

FASTA = "/mnt/ssd-8tb/satellite_dna_lm/monomers.fasta"   # curated satellites, <=500/sp, >=20bp, upper
BASE  = Path("/home/jg2070/Desktop/dtol_review_August/2026_trees/annotation_centromeres")
OUT   = Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity/figures")

CORE = "[CT]TTCGTTGGAA[AG]CGGGA"                 # canonical 17-bp CENP-B box (Masumoto 1989)
P2, P3 = (regex.compile("(?:%s){e<=%d}" % (CORE, k)) for k in (2, 3))
comp = str.maketrans("ACGTN", "TGCAN"); rc = lambda s: s.translate(comp)[::-1]
def hit(p, s): return bool(p.search(s) or p.search(rc(s)))

xl = pd.read_excel(BASE/"DTOL_327_master_March.xlsx")
xl["code"] = xl["fasta"].astype(str).str.lower().str.replace(r"[0-9].*$","",regex=True).str.replace(r"[.].*$","",regex=True)
sp2taxa = dict(zip(xl["code"], xl["taxa1"]))

C = collections.defaultdict(collections.Counter)
name = None; seq = []
def flush():
    if not name: return
    s = "".join(seq)
    if len(s) < 17: return
    sp = name.split("__")[0].lower(); c = C[sp]; c["n"] += 1
    if hit(P2, s): c["k2"] += 1
    if hit(P3, s): c["k3"] += 1
    if hit(P3, "".join(random.sample(s, len(s)))): c["null"] += 1   # composition null
for line in open(FASTA):
    if line[0] == ">": flush(); name = line[1:].strip(); seq = []
    else: seq.append(line.strip())
flush()

rows = [dict(species=sp, taxa1=sp2taxa.get(sp,"NA"), n=c["n"],
             pct_le2mm=100*c["k2"]/c["n"], pct_le3mm=100*c["k3"]/c["n"],
             pct_null_le3mm=100*c["null"]/c["n"]) for sp,c in C.items() if c["n"]]
df = pd.DataFrame(rows)
df.sort_values("pct_le3mm", ascending=False).to_csv(OUT/"cenpb_box_per_species.tsv", sep="\t", index=False)

def grp(t):
    if t == "Mammalia": return "Mammalia"
    if t in ("Actinopterygii","Aves","Reptilia","Tunicata"): return "other Chordata"
    if t in ("Algae","Bryophyta","Dicots","Monocots"): return "Viridiplantae"
    if t == "Fungi": return "Fungi"
    return "Invertebrates"
df["group"] = df["taxa1"].map(grp)
g = df.groupby("group").apply(lambda x: pd.Series({
    "n_species": len(x), "n_seqs": x["n"].sum(),
    "pct_le2mm": np.average(x["pct_le2mm"], weights=x["n"]),
    "pct_le3mm": np.average(x["pct_le3mm"], weights=x["n"]),
    "pct_null_le3mm": np.average(x["pct_null_le3mm"], weights=x["n"])}), include_groups=False).reset_index()
g["enrichment"] = g["pct_le3mm"] / g["pct_null_le3mm"].replace(0, np.nan)
g = g.sort_values("pct_le3mm", ascending=False)
g.to_csv(OUT/"cenpb_box_per_clade.tsv", sep="\t", index=False)
print(g.to_string(index=False))

fig, ax = plt.subplots(figsize=(8,5))
x = np.arange(len(g)); cols = ["#F72485" if t=="Mammalia" else "#999" for t in g["group"]]
ax.bar(x-0.2, g["pct_le3mm"], 0.4, color=cols, label="satellites (≤3 mismatch)")
ax.bar(x+0.2, g["pct_null_le3mm"], 0.4, color="#ccc", label="composition-shuffled null")
ax.set_xticks(x); ax.set_xticklabels(g["group"], rotation=20, ha="right")
ax.set_ylabel("% satellite monomers with a CENP-B box"); ax.legend()
ax.set_title("Canonical CENP-B box in curated satellites (Masumoto motif, ≤3 mm)")
fig.tight_layout(); fig.savefig(OUT/"cenpb_box_by_clade.png", dpi=200)
print("saved cenpb_box_{per_species,per_clade}.tsv, cenpb_box_by_clade.png")
