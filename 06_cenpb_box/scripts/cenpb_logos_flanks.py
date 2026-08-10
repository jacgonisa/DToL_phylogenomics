#!/usr/bin/env python3
"""CENP-B box sequence logos with +/-5 bp flanks, as the songbird paper did
(Formenti et al., Cell 2026, Suppl. Fig. 15) to check whether a candidate is a
genuinely ENRICHED motif: the conserved 17-bp box should carry high information
content while the +/-5 flanking positions collapse to ~random (a built-in
negative control -> 'sharp transition from box to surrounding random sequence').

One logo per species: the human HG002 alpha-sat benchmark, the bat Rhinolophus
sinicus (closest non-human by the strict scan), and every labile candidate
(edit-distance candidate above its shuffled null). Compiled into one PDF.

Boxes are located with fixed 17-bp (substitution) windows so the flanks are
well-defined: <=2 mismatches for the same species (human), <=5 for cross-species
(the songbird thresholds). Both strands. Vertical lines mark the box boundaries."""
import regex, math, collections
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import logomaker

SAT   = Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity")
FASTA = "/mnt/ssd-8tb/satellite_dna_lm/monomers.fasta"
HUMAN = SAT/"cenpb_psi/human_alpha.fasta"
CORE  = "[CT]TTCGTTGGAA[AG]CGGGA"                     # canonical 17-bp box
FLANK = 5
PAT   = {k: regex.compile("(?:%s){s<=%d}" % (CORE, k)) for k in (2, 5)}  # substitutions -> fixed 17bp
comp  = str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(comp)[::-1]
MAXWIN= 3000                                          # cap windows per species for a clean logo

def windows(seqs, kmax):
    out=[]
    for seq in seqs:
        for strand in (seq, rc(seq)):
            for m in PAT[kmax].finditer(strand, overlapped=False):
                s,e=m.start(),m.end()
                if e-s!=17 or s-FLANK<0 or e+FLANK>len(strand): continue
                w=strand[s-FLANK:e+FLANK]
                if len(w)==27 and set(w)<=set("ACGT"): out.append(w)
                if len(out)>=MAXWIN: return out
    return out

# ---- species -> readable name & clade ----
tax=pd.read_csv("/home/jg2070/Desktop/dtol_review_August/2026_trees/annotation_centromeres/centromere_code_to_species.tsv",sep="\t")
tax["code"]=tax["fasta"].astype(str).str.lower().str.replace(r"[0-9.].*$","",regex=True)
name={c:f"{g} {s}" for c,g,s in zip(tax.code,tax.genus,tax.species)}

# ---- candidate list from the labile table ----
lab=pd.read_csv(SAT/"figures/cenpb_songbird_labile_per_species.tsv",sep="\t")
nh=lab[lab.species!="human"].copy()
nh["cand"]=(nh.labile_obs>=nh.labile_null+3)&(nh.labile_obs>=2*nh.labile_null.clip(lower=1))
cands=nh[nh.cand].sort_values("enrichment",ascending=False).species.tolist()
order=["human","mrhisin"]+[c for c in cands if c!="mrhisin"]   # benchmark, bat, then candidates
clade_of=dict(zip(lab.species,lab.clade))

# ---- load monomers only for species we need ----
need=set(order)
seqs=collections.defaultdict(list); cur=None; buf=[]
for line in open(FASTA):
    if line[0]==">":
        if cur in need and buf: seqs[cur].append("".join(buf))
        cur=line[1:].split("__")[0].strip().lower(); buf=[]
    else: buf.append(line.strip().upper())
if cur in need and buf: seqs[cur].append("".join(buf))
# human benchmark
hs=[]; buf=[]
for line in open(HUMAN):
    if line[0]==">":
        if buf: hs.append("".join(buf)); buf=[]
    else: buf.append(line.strip().upper())
if buf: hs.append("".join(buf))
seqs["human"]=hs

def logo_matrix(wins):
    counts=np.zeros((27,4)); idx={b:i for i,b in enumerate("ACGT")}
    for w in wins:
        for i,ch in enumerate(w): counts[i,idx[ch]]+=1
    prob=counts/counts.sum(1,keepdims=True)
    prob=pd.DataFrame(prob,columns=list("ACGT"))
    return logomaker.transform_matrix(prob,from_type="probability",to_type="information")

# ---- first pass: compute box vs flank information per species ----
data={}
for sp in order:
    kmax=2 if sp=="human" else 5
    wins=windows(seqs.get(sp,[]),kmax)
    rec=dict(species=sp,name=name.get(sp,sp),
             clade=("Human" if sp=="human" else clade_of.get(sp,"?")),n_windows=len(wins))
    if len(wins)>=5:
        info=logo_matrix(wins)
        box=info.iloc[FLANK:FLANK+17].sum(1); flank=pd.concat([info.iloc[:FLANK],info.iloc[FLANK+17:]]).sum(1)
        rec["mean_box_bits"]=round(float(box.mean()),2); rec["mean_flank_bits"]=round(float(flank.mean()),2)
        rec["box_minus_flank"]=round(float(box.mean()-flank.mean()),2); rec["_info"]=info
    else:
        rec["mean_box_bits"]=rec["mean_flank_bits"]=rec["box_minus_flank"]=float("nan")
    data[sp]=rec

# order pages: human & bat first, then by box-minus-flank enrichment (most box-specific first)
drawn=[sp for sp in order if data[sp]["n_windows"]>=5]
lead=[sp for sp in ("human","mrhisin") if sp in drawn]
rest=sorted([sp for sp in drawn if sp not in lead], key=lambda s:-data[s]["box_minus_flank"])
pages=lead+rest

# ---- summary scatter: box vs flank information (the negative-control test) ----
pdf=PdfPages(SAT/"figures/cenpb_box_logos_flanks_all_candidates.pdf")
ccol={"Human":"#d62728","Vertebrates":"#1565C0","Invertebrate":"#EF6C00","Viridiplantae":"#2E7D32"}
figs,axs=plt.subplots(figsize=(7,6.4))
for sp in drawn:
    r=data[sp]; axs.scatter(r["mean_flank_bits"],r["mean_box_bits"],s=45,
        color=ccol.get(r["clade"],"grey"),edgecolor="k",lw=0.4,zorder=3,alpha=0.9)
axs.plot([0,2],[0,2],ls="--",color="grey",lw=1.2,label="box = flank (no enrichment)")
for sp in ("human","mrhisin"):
    if sp in drawn:
        r=data[sp]; axs.annotate(r["name"].split()[0],(r["mean_flank_bits"],r["mean_box_bits"]),
            fontsize=8,fontweight="bold",xytext=(4,4),textcoords="offset points")
axs.set_xlabel("mean information in ±5 flank (bits)  — negative control")
axs.set_ylabel("mean information in 17-bp box (bits)")
axs.set_title("CENP-B box vs flank information (songbird ±5-flank test)\nonly points ABOVE the diagonal are box-enriched motifs",
              fontsize=10.5,fontweight="bold")
axs.set_xlim(0,2); axs.set_ylim(0,2)
from matplotlib.lines import Line2D
axs.legend(handles=[Line2D([0],[0],marker='o',ls='',mfc=c,mec='k',label=g) for g,c in ccol.items() if g in [data[s]['clade'] for s in drawn]]
           +[Line2D([0],[0],ls='--',color='grey',label='box = flank')],fontsize=8.5,frameon=False)
axs.spines[["top","right"]].set_visible(False); plt.tight_layout()
pdf.savefig(figs,bbox_inches="tight"); plt.close(figs)

# ---- one logo per species ----
for sp in pages:
    r=data[sp]; info=r["_info"]
    fig,ax=plt.subplots(figsize=(8.5,2.6))
    logomaker.Logo(info,ax=ax,color_scheme="classic",show_spines=False)
    ax.axvline(FLANK-0.5,color="grey",ls="--",lw=1); ax.axvline(FLANK+17-0.5,color="grey",ls="--",lw=1)
    ax.set_ylim(0,2); ax.set_ylabel("bits"); ax.set_xticks([])
    lbl="Human α-sat (HG002 benchmark, ≤2 mm)" if sp=="human" else f"{r['name']} ({sp}, ≤5 mm)"
    ax.set_title(f"{lbl}  |  {r['clade']}  |  n={r['n_windows']}  |  "
                 f"box {r['mean_box_bits']} vs flank {r['mean_flank_bits']} bits  (Δ={r['box_minus_flank']:+.2f})",
                 fontsize=9,fontweight="bold")
    ax.text(FLANK+8,-0.28,"17-bp CENP-B box",ha="center",va="top",fontsize=8,color="grey")
    ax.text(2,-0.28,"−5 flank",ha="center",va="top",fontsize=7.5,color="grey")
    ax.text(24,-0.28,"+5 flank",ha="center",va="top",fontsize=7.5,color="grey")
    plt.tight_layout(); pdf.savefig(fig,bbox_inches="tight"); plt.close(fig)
pdf.close()
summ=pd.DataFrame([{k:v for k,v in data[sp].items() if k!="_info"} for sp in order])
summ=summ.sort_values("box_minus_flank",ascending=False)
summ.to_csv(SAT/"figures/cenpb_box_logos_flanks_summary.tsv",sep="\t",index=False)
print("Saved PDF:",SAT/"figures/cenpb_box_logos_flanks_all_candidates.pdf")
print(f"species logos drawn: {len(drawn)} / {len(order)}")
print(summ[["name","species","clade","n_windows","mean_box_bits","mean_flank_bits","box_minus_flank"]].to_string(index=False))
