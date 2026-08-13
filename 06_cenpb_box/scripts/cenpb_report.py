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
{n_sat} DToL species with satellites</b>. A ±5-flank test flags <b>box-like</b> motifs in 20
species, strongest in <b>birds</b>: the goshawk consensus <code>CTTTTTTGGAAACGGGA</code> is
<b>15/17 identical</b> to the canonical box, differing only at the CpG (positions 4–5,
CpG→TT). These are strong <b>candidate divergent boxes</b> — near-canonical but lacking the
CENP-B-essential CpG (the exact broad motif, which fixes that CpG, therefore scores 0). Their
functionality is untested (no protein-binding assay).</div>

<h2>What was searched</h2>
<p><span class="tag">input</span><code>all.satellites.txt</code> — every <b>candidate</b>
satellite-repeat <i>array</i> called across the DToL assemblies (candidates, not curated
satellites). Restricted to the <b>325 published species</b> (tree tips); 12 out-of-set
species excluded. <b>"Uncapped"</b> = every candidate array of each species, both strands
(~{n_arr} arrays across the <b>{n_sat} of 325 species</b> that carry satellite annotations)
— not a ≤500-monomer subsample.</p>

<h2>Human HG002 benchmarks</h2>
<p>Both methods were validated on human satellites: α-satellite (positive; carries the box)
and HSat 1/2/3 (negative; lacks it).</p>
{tbl(bench_df)}
<p class="sub">Canonical boxes are ~50,000× denser per Mbp in α-satellite than HSat; both
methods separate the controls cleanly.</p>

<h2>Method 1 — exact IUPAC motif parsing (Fachinetti)</h2>
<p>Three <b>exact IUPAC</b> motif tiers (Barra &amp; Fachinetti, bioRxiv 2026.05.25.727640),
matched on both strands (Y = [C/T], R = [A/G], N = any base):</p>
<table class="tbl" style="max-width:560px">
<tr><th>tier</th><th>IUPAC motif</th><th>note</th></tr>
<tr><td>canonical</td><td><code>YTTCGTTGGAARCGGGA</code></td><td>the functional box</td></tr>
<tr><td>broad</td><td><code>YTTCGNNNNANRCGGGN</code></td><td>Sugimoto 1998; intermediate</td></tr>
<tr><td>degenerated</td><td><code>NTTCGNNNNANNCGGGN</code></td><td>most permissive (N at pos 1 &amp; 12)</td></tr>
</table>
<p>Each is counted per species and scored as <b>enrichment = observed / expected</b>, where the
expected count comes from a <b>dinucleotide-preserving null</b> (first-order Markov ≈
Altschul–Erikson doublet shuffle; validated against an explicit shuffle on the α-sat benchmark,
4,283× vs 4,524×). Nesting: canonical ⊂ broad ⊂ degenerated.</p>
<p>The two quantities are shown first as bar charts, then combined in a scatter.</p>
<div class="fig">{img(FIG/'cenpb_bars_enrichment.png',w="70%")}<div class="cap">
<b>(i) Enrichment</b> over the dinucleotide null (obs/exp), per clade × tier (broken y-axis;
α-satellite ≫ DToL). Canonical = 0 in all DToL clades; HSat canonical = 5 boundary hits (hatched).</div></div>
<div class="fig">{img(FIG/'cenpb_bars_density.png',w="70%")}<div class="cap">
<b>(ii) Motif density</b> — exact IUPAC hits per Mbp, per clade × tier (broken y-axis).</div></div>
<div class="fig">{img(FIG/'cenpb_paper_motifs.png',w="75%")}<div class="cap">
<b>(iii) Scatter</b> of enrichment (i) vs density (ii). A real box needs <b>both</b>: only
<b>α-satellite</b> sits top-right. HSat canonical is a lone point at ~0.03/Mbp (5 boundary hits →
high enrichment but ~0 density — not a box). DToL clades sit near the null (enrichment ≈ 1).</div></div>
{tbl(clade_df.round(2))}
<p class="sub"><b>Note.</b> This uses composition-corrected <i>enrichment</i>, not raw
hits/Mbp — raw density inflates plants (AT-rich, more satellite bp). Under the dinucleotide
null the looser-tier enrichments are modest (plants lead the broad tier, 3.0×; vertebrates the
degenerate tier, 1.4×) and the apparent plant broad signal deflates (11.7× → 3.0×).</p>

<h2>Method 2 — songbird ±5-flank test (Formenti et al., Cell 2026)</h2>
<p>This is a permissive, distance-based search (following the zebra-finch T2T paper, Formenti et
al. 2026, Suppl. Fig. 15). It runs per species on the uncapped satellite arrays, in six steps:</p>
<ol style="font-size:13.5px;line-height:1.5">
<li><b>Find windows.</b> Slide along every satellite array, on <b>both strands</b>, and record every
17-bp stretch that matches the canonical box <code>[CT]TTCGTTGGAA[AG]CGGGA</code> within
<b>≤5 substitutions</b> (fuzzy matching). Each such 17-bp hit is a <b>window</b>. We record
<code>n_windows</code> (total per species) and <b>prevalence</b> = windows per Mbp of satellite.</li>
<li><b>Add flanks.</b> Extend each window by 5 bp on each side → a <b>27-bp window</b> (5 + 17 + 5).
The flanks are a built-in negative control.</li>
<li><b>Stack into a position matrix.</b> Pile all of a species' 27-bp windows and count, at each
column, how often each base (A/C/G/T) occurs — a position frequency matrix.</li>
<li><b>Information content (conservation).</b> Per column, IC = 2 + Σ<sub>b</sub> p<sub>b</sub>
log₂ p<sub>b</sub> bits (0 = random, 2 = one fixed base). Then <b>motif information</b> = mean IC
over the 17 motif columns, <b>flank information</b> = mean IC over the 10 flank columns, and
<b>Δ = motif − flank</b>. A real motif is <b>conserved within itself but random in the flanks</b>
→ high motif information, low flank information → Δ &gt; 0 (the scatter's y and x axes).</li>
<li><b>Consensus &amp; identity.</b> The <b>consensus</b> = the most common base at each of the 17
motif columns. <b>Identity to canonical</b> = the fraction of the 17 positions whose consensus
base is allowed by the canonical IUPAC, × 100% (the scatter's colour; equivalently 17 − substitutions).</li>
<li><b>Shuffle-the-motif null.</b> Because a 17-bp motif matches partly by base composition, we
shuffle the consensus (same bases, random order) many times and recompute identity → the chance
level (≈30%). A motif is "real" when observed identity ≫ its own shuffled null (see below).</li>
</ol>
<p class="sub"><b>Definitions.</b> <i>window</i> = one located 17-bp box-like hit; <i>consensus</i>
= the per-position majority base over a species' windows; <i>Δ</i> = motif − flank information;
<i>prevalence</i> = windows/Mbp. <b>Caveat:</b> windows are selected to be ≤5 subs from canonical,
so identity has a floor; and for very high-copy satellites the windows deviate at <i>random</i>
positions that cancel on averaging, so the consensus returns to canonical (≈100%) — the
degenerate-consensus effect, not a true perfect box. A genuine motif deviates <i>consistently</i>
(e.g. the goshawk's eroded CpG), keeping its consensus off-canonical.</p>
<div class="fig">{img(FIG/'cenpb_flank_uncapped_scatter.png')}<div class="cap">
Box vs flank information, coloured by <b>identity to the canonical CENP-B motif</b>. A
<b>candidate box</b> = high identity <i>and</i> above the motif=flank diagonal (motif-specific).
α-satellite (red star) = the functional box; HSat (grey ×) = negative. Birds (goshawk, takahē,
ptarmigan) are the most box-specific vertebrates; falcons sit on the diagonal (conserved
satellite, not a candidate box).</div></div>
<h3>Top vertebrates by flank Δ (box − flank)</h3>
{tbl(vert_df.round(2))}
<p class="sub">Box-like motifs (Δ≥0.5, ≤2 substitutions from canonical) occur in <b>20 species
across clades</b>, strongest in birds; the scatter is coloured by identity to the canonical
motif (warm = box-like).</p>
<h3>Is the box-like identity real? (shuffle-the-motif null)</h3>
<p>A short 17-bp motif matches partly by base composition alone, so we calibrated identity by
<b>shuffling the motif</b> (same composition, order destroyed). <b>identity</b> = (positions
matching the canonical IUPAC) / 17 × 100%; <b>null</b> = the same on shuffles of that motif.</p>
<div class="fig">{img(FIG/'cenpb_identity_shuffle_null.png',w="52%")}<div class="cap">
Each observed motif (right) is linked to the mean of its own shuffled null (left). A motif of this
composition matches canonical only <b>≈30% by chance</b>, whereas observed consensus identities are
<b>65–100%</b> (median <b>71%</b>, +44 above chance). So the matches are <b>real arrangement
similarity, not a short-motif composition artifact</b>. (Points near 100% are high-copy satellites
whose windows scatter around canonical, so their consensus averages back to canonical — the
degenerate-consensus effect, not a perfect box.) <code>cenpb_identity_shuffle_null.py</code></div></div>

<h3>Why the two methods disagree on birds</h3>
<div class="fig">{img(FIG/'cenpb_goshawk_alignment.png')}<div class="cap">
The goshawk motif is 15/17 identical to canonical; its only two substitutions land on the 5′
<b>CpG</b> (positions 4–5) that <i>all three</i> IUPAC tiers hold fixed (the TTCG anchor). So
position-based IUPAC matching scores <b>0</b>, while the distance-based songbird ≤5-substitution
search (no fixed position) finds 1,722 windows. The eroded position is a CpG — a CENP-B contact
base and a methylation/deamination hotspot — which is exactly why it is a <b>candidate divergent
box</b>, not a canonical one.</div></div>

<h2>CENP-B box across the chronogram</h2>
<div class="fig">{img(FIG/'cenpb_box_tree_325sp.png')}<div class="cap">
★ = candidate box (box-specific + near-canonical; 20 species, strongest in birds). Ring 1
(green) = candidate-box signal (flank Δ, Method 2); ring 2 (purple) = CENP-B-like motif density
(broad matches/Mbp, Method 1). Tips coloured by clade. No functional box confirmed.</div></div>
<p class="sub">Per-species logos (box ±5 flank), including the human α-sat and HSat controls:
<code>figures/cenpb_box_logos_flanks_VERTEBRATES_uncapped.pdf</code>.</p>

<footer>Reproducible scripts in <code>06_cenpb_box/scripts/</code>:
<code>cenpb_paper_motifs.py</code>, <code>cenpb_flank_uncapped.py</code>,
<code>cenpb_human_benchmark.py</code>, <code>ae_shuffle.py</code>,
<code>cenpb_paper_motifs_plot.py</code>, <code>cenpb_flank_scatter.py</code>,
<code>cenpb_tree_325sp.R</code>. Allowlist: <code>species_325.txt</code>.</footer>
</body></html>"""
out=SAT/"cenpb_box_report.html"; out.write_text(html)
print("Saved:",out,"|",round(len(html)/1e6,2),"MB")
