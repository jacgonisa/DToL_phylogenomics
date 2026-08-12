#!/usr/bin/env python3
"""Build a self-contained HTML report for the CENP-B box analysis (two methods +
HG002 controls). Figures are base64-embedded so the single .html works anywhere."""
import base64, datetime
from pathlib import Path
import pandas as pd
SAT=Path("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity")
FIG=SAT/"figures"

def img(p, w="100%"):
    b=base64.b64encode(Path(p).read_bytes()).decode()
    return f'<img style="max-width:{w};height:auto;border:1px solid #e5e5e5;border-radius:6px" src="data:image/png;base64,{b}">'

bench=pd.read_csv(FIG/"cenpb_human_benchmark.tsv",sep="\t").set_index("control")
per=pd.read_csv(FIG/"cenpb_paper_motifs_per_clade.tsv",sep="\t")
flank=pd.read_csv(FIG/"cenpb_flank_uncapped_per_species.tsv",sep="\t")
n_sat=int(per.n_species.sum()); n_arr="20.3 million"

def tbl(df, cls="tbl"):
    return df.to_html(index=False, classes=cls, border=0, justify="center")

# benchmark table
a=bench.loc["alpha-satellite (positive)"]; h=bench.loc["HSat1/2/3 (negative)"]
bench_df=pd.DataFrame([
 ["α-satellite (positive)","50,000 monomers · 8.55 Mbp",f"{int(a.canonical):,} ({a.canonical_perMbp:.0f}/Mbp)",
  f"{a.broad_enrich_dinuc:,.0f}×",f"{a.delta:.2f} ({a.win_per_Mbp:.0f}/Mbp)"],
 ["HSat 1/2/3 (negative)","7,982 arrays · 161.8 Mbp",f"{int(h.canonical)} ({h.canonical_perMbp:.2g}/Mbp)",
  f"{h.broad_enrich_dinuc:.2f}× (depleted)",f"{h.delta:.2f} ({h.win_per_Mbp:.0f}/Mbp)"]],
 columns=["HG002 control","size","canonical box (M1)","broad enrichment, dinuc null","flank Δ (M2)"])

# per-clade enrichment table (dinucleotide null)
clade_df=per.rename(columns={"clade":"clade","broad_enrich_dinuc":"broad","degenerated_enrich_dinuc":"degenerate",
    "canonical_enrich_dinuc":"canonical"})[["clade","n_species","canonical","broad","degenerate"]].copy()
clade_df["clade"]=clade_df["clade"].replace({"Viridiplantae":"Plants","Invertebrate":"Invertebrates"})

# top vertebrates by flank delta
v=flank[flank.clade=="Vertebrates"].sort_values("delta",ascending=False).head(8)
vert_df=v[["name","vgroup","n_windows","win_per_Mbp","mean_box_bits","mean_flank_bits","delta","box_consensus","subs_vs_canonical"]]\
    .rename(columns={"name":"species","vgroup":"group","win_per_Mbp":"boxes/Mbp","mean_box_bits":"box bits",
                     "mean_flank_bits":"flank bits","delta":"Δ","box_consensus":"consensus","subs_vs_canonical":"subs"})

CSS="""
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;color:#1a1a1a;line-height:1.55;
 max-width:960px;margin:32px auto;padding:0 20px}
h1{font-size:26px;margin:0 0 4px} h2{font-size:19px;margin:30px 0 8px;border-bottom:2px solid #0072B2;padding-bottom:4px}
h3{font-size:15px;margin:18px 0 6px;color:#333} .sub{color:#666;font-size:13px;margin-top:0}
.key{background:#eef6fc;border-left:4px solid #0072B2;padding:12px 16px;border-radius:4px;margin:16px 0}
.tbl{border-collapse:collapse;width:100%;font-size:13px;margin:10px 0}
.tbl th{background:#37474f;color:#fff;padding:7px 9px;text-align:center;font-weight:600}
.tbl td{padding:6px 9px;border-bottom:1px solid #eee;text-align:center}
.tbl tr:nth-child(even){background:#fafafa}
.fig{margin:14px 0;text-align:center} .cap{color:#555;font-size:12.5px;margin-top:6px}
code{background:#f3f3f3;padding:1px 5px;border-radius:3px;font-size:12.5px}
.tag{display:inline-block;background:#0072B2;color:#fff;font-size:11px;padding:2px 8px;border-radius:10px;margin-right:6px}
footer{color:#888;font-size:12px;margin-top:40px;border-top:1px solid #eee;padding-top:12px}
"""

html=f"""<!doctype html><html><head><meta charset="utf-8"><title>CENP-B box screen — DToL</title>
<style>{CSS}</style></head><body>
<h1>CENP-B box screen across the DToL satellite repertoire</h1>
<p class="sub">Generated {datetime.date.today().isoformat()} · module <code>06_cenpb_box</code></p>

<div class="key"><b>Terminology.</b> A CENP-B <b>box</b> is a <i>functional</i> motif that
binds the CENP-B protein. Here we detect <i>sequence</i> matches only, so we call them
<b>motifs</b> / <b>candidate boxes</b>; none is a confirmed functional box (protein binding
not tested).<br><br>
<b>Bottom line.</b> The canonical CENP-B motif (<code>[CT]TTCGTTGGAA[AG]CGGGA</code>) as an
exact match is confined to human α-satellite — <b>0 exact canonical motifs in any of the
{n_sat} DToL species with satellites</b>. A permissive ±5-flank test flags box-specific,
near-canonical <b>candidate boxes</b> in 20 species across clades, strongest in <b>birds</b>
(goshawk, ptarmigan, takahē); these have lost the CENP-B-essential 5′ CG — a likely
TIGD4-type avian candidate box, not a mammalian CENP-B box.</div>

<h2>What was searched</h2>
<p><span class="tag">input</span><code>all.satellites.txt</code> — every annotated satellite
<i>array</i> across the DToL assemblies. Restricted to the <b>325 published species</b>
(tree tips); 12 out-of-set species excluded. <b>"Uncapped"</b> = every array of each
species, both strands (~{n_arr} arrays across the <b>{n_sat} of 325 species</b> that carry
satellite annotations) — not a ≤500-monomer subsample.</p>

<h2>Human HG002 benchmarks</h2>
<p>Both methods were validated on human satellites: α-satellite (positive; carries the box)
and HSat 1/2/3 (negative; lacks it).</p>
{tbl(bench_df)}
<p class="sub">Canonical boxes are ~50,000× denser per Mbp in α-satellite than HSat; both
methods separate the controls cleanly.</p>

<h2>Method 1 — exact IUPAC motif parsing (Fachinetti)</h2>
<p>Three exact IUPAC tiers (canonical / broad / degenerate; Barra &amp; Fachinetti,
bioRxiv 2026.05.25.727640), both strands, scored as enrichment over a
<b>dinucleotide-preserving null</b> (first-order Markov ≈ Altschul–Erikson doublet shuffle;
validated against an explicit shuffle on the α-sat benchmark, 4,283× vs 4,524×).</p>
<div class="fig">{img(FIG/'cenpb_paper_motifs.png')}<div class="cap">
Enrichment over the dinucleotide null. Canonical box = 0 in all DToL clades; α-satellite
positive control is off scale. Vertebrates lead the broad tier (1.39×); plants only lead the
degenerate/stochastic tier (3.01×).</div></div>
{tbl(clade_df.round(2))}
<p class="sub"><b>Note.</b> This uses composition-corrected <i>enrichment</i>, not raw
hits/Mbp — raw density inflates plants (AT-rich, more satellite bp); under the dinucleotide
null vertebrates actually edge plants on the broad tier.</p>

<h2>Method 2 — songbird ±5-flank test (Formenti et al., Cell 2026)</h2>
<p>The box matched as 17-bp windows (≤5 substitutions, both strands); a real motif has high
information across the box that <b>collapses in the ±5 flanks</b> (Δ = box − flank). Reported
with the box consensus, substitutions from canonical, and prevalence.</p>
<div class="fig">{img(FIG/'cenpb_flank_uncapped_scatter.png')}<div class="cap">
Box-motif vs flank information, coloured by <b>identity to the canonical CENP-B motif</b>. A
<b>candidate box</b> = high identity <i>and</i> above the box=flank diagonal (box-specific).
α-satellite (red star) = the functional box; HSat (grey ×) = negative. Birds (goshawk, takahē,
ptarmigan) are the most box-specific vertebrates; falcons sit on the diagonal (conserved
satellite, not a candidate box).</div></div>
<h3>Top vertebrates by flank Δ (box − flank)</h3>
{tbl(vert_df.round(2))}
<p class="sub">Candidate boxes (Δ≥0.5, ≤2 substitutions from canonical) occur in <b>20 species
across clades</b>, not only birds — but birds carry the strongest signal.</p>

<h2>CENP-B box across the chronogram</h2>
<div class="fig">{img(FIG/'cenpb_box_tree_325sp.png')}<div class="cap">
★ = candidate box (box-specific + near-canonical; 20 species, strongest in birds). Ring 1
(green) = candidate-box signal (flank Δ, Method 2); ring 2 (purple) = CENP-B-like motif density
(broad matches/Mbp, Method 1). Tips coloured by clade. No functional box confirmed.</div></div>
<p class="sub">Per-species logos (box ±5 flank), including the human α-sat and HSat controls:
<code>figures/cenpb_box_logos_flanks_VERTEBRATES_uncapped.pdf</code>.</p>

<h2>Conclusion</h2>
<p>The canonical functional CENP-B box is a <b>human/mammalian α-satellite feature</b>; by the
strict motif definition it is absent from all {n_sat} DToL species with satellites. The
±5-flank test detects a <b>diverged, box-shaped motif strongest in birds</b>, but it lacks the
CENP-B-essential 5′ CG — consistent with a <b>TIGD4-type avian box</b> rather than mammalian
CENP-B (cf. the songbird paper). Broad low-signal hits in plants/invertebrates fail the
composition-corrected controls and reflect homogenised satellite, not a functional box.</p>

<footer>Reproducible scripts in <code>06_cenpb_box/scripts/</code>:
<code>cenpb_paper_motifs.py</code>, <code>cenpb_flank_uncapped.py</code>,
<code>cenpb_human_benchmark.py</code>, <code>ae_shuffle.py</code>,
<code>cenpb_paper_motifs_plot.py</code>, <code>cenpb_flank_scatter.py</code>,
<code>cenpb_tree_325sp.R</code>. Allowlist: <code>species_325.txt</code>.</footer>
</body></html>"""
out=SAT/"cenpb_box_report.html"; out.write_text(html)
print("Saved:",out,"|",round(len(html)/1e6,2),"MB")
