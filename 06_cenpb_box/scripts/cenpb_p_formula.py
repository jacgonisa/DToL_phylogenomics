#!/usr/bin/env python3
"""Render the combinatorial formula for p (probability a random 17-mer is <=5 subs
from the CENP-B box) as a small image, so the self-contained report can show it."""
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
SAT=Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity")
plt.rcParams.update({"mathtext.fontset":"cm"})
lines=[
 r"$N_k=\sum_{i+j=k}\ \binom{15}{i}\,3^{\,i}\ \binom{2}{j}\,2^{\,j}\,2^{\,2-j}\qquad(0\leq i\leq 15,\ 0\leq j\leq 2)$",
 r"$p=\frac{1}{4^{17}}\sum_{k=0}^{5} N_k=\frac{4\,458\,112}{4^{17}}\approx 2.6\times10^{-4}\approx\frac{1}{3854}$",
 r"$E=D\cdot p$",
]
fig=plt.figure(figsize=(7.6,2.2)); fig.patch.set_facecolor("white")
y=0.90
for ln in lines:
    fig.text(0.02,y,ln,fontsize=15,va="top",ha="left"); y-=0.36
for ext in ("png","pdf"):
    fig.savefig(SAT/f"figures/cenpb_p_formula.{ext}",dpi=300,bbox_inches="tight",facecolor="white",pad_inches=0.12)
print("Saved figures/cenpb_p_formula.png/pdf")
