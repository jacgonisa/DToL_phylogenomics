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
<li><b>Consensus &amp; identity.</b> The <b>consensus</b> = the most common base per motif column
(used for the scatter colour and the reported substitutions-from-canonical). For the null test
we score identity <b>per window</b> — the mean, over all of a species' 17-bp windows, of
(positions matching the canonical IUPAC)/17 × 100%. Doing this per window (rather than on the
consensus) avoids the <b>degenerate-consensus artifact</b>: a high-copy satellite whose windows
deviate at <i>random</i> positions has a canonical-looking consensus (≈100%), yet its individual
windows sit at the ≤5-substitution floor (~70%).</li>
<li><b>Shuffle-the-window null.</b> Each window's bases are shuffled (composition preserved, order
destroyed) and identity re-scored → the chance level (≈29%). Observed per-window identity
(~70–81%) sits far above chance (windows are compositionally box-like), though near the floor.</li>
</ol>
<p class="sub"><b>Definitions.</b> <i>window</i> = one located 17-bp box-like hit; <i>consensus</i>
= the per-position majority base over a species' windows; <i>Δ</i> = motif − flank information;
<i>prevalence</i> = windows/Mbp. <b>Note:</b> per-window identity is ~70% for every species (the
≤5-substitution floor), so it does not single out particular species. The species-specific
box-likeness (e.g. the goshawk's 15/17 consensus) is a <b>consensus/coherence</b> property — its
windows deviate <i>consistently</i> (the eroded CpG) so the average keeps that off-canonical
signature — and is best seen in the logos and the alignment below, together with the flank Δ.</p>
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
<h3>How many species "recover" a box-like motif?</h3>
<p>Of the {n_sat} species with satellites, <b>152 (94%)</b> contain at least one ≤5-substitution
box-like window; 6 have only 1–4 weak hits; and <b>4 recover nothing</b> (<i>Bibio marci</i>,
<i>Limnoperna fortunei</i>, <i>Pholis gunnellus</i>, <i>Lemur catta</i> — AT-rich satellites lacking
the GC-rich <code>TTCG…CGGG</code> core). But <b>≤5 substitutions on 17 bp is very permissive</b> —
nearly any GC-containing satellite has such a window by chance — so 94% recovery is <i>expected</i>,
not evidence of real boxes. This is precisely why the unbiased test below is needed.</p>
<p class="sub"><b>Why ≤5 subs matches so easily — the numbers game.</b> The threshold is loose: a
window needs only ≥ 12/17 = 70.6% identity. From a Monte-Carlo over random 17-mers, the probability
that one matches the box is:</p>
<table class="tbl" style="max-width:440px">
<tr><th>threshold</th><th>P(random 17-mer matches)</th><th>≈ 1 in</th></tr>
<tr><td>≤2 substitutions</td><td>&lt; 5×10⁻⁷</td><td><b>&gt; 2,000,000</b></td></tr>
<tr><td>≤3 substitutions</td><td>4.5×10⁻⁶</td><td>222,000</td></tr>
<tr><td>≤4 substitutions</td><td>4.3×10⁻⁵</td><td>23,000</td></tr>
<tr><td><b>≤5 substitutions</b></td><td><b>2.9×10⁻⁴</b></td><td><b>~3,500</b></td></tr>
</table>
<p class="sub">So ~1 in 3,500 random 17-mers already matches at ≤5 subs. Now multiply by the number
of windows: a species with ~100,000 satellite arrays × ~300 bp × 2 strands ≈ <b>60 million windows</b>
→ <b>~17,000 chance matches expected</b> — recovery is essentially guaranteed for any GC-containing
satellite. The meaningful flip side: <b>≤2 substitutions is &lt; 1 in 2 million</b> (essentially
never by chance), which is exactly why <b>Method 1 finding ~0 exact canonical/broad/degenerate hits
in DToL is real signal</b>, whereas the ≤5-substitution abundance is just the loose threshold.
(Matching code: <code>cenpb_flank_uncapped.py</code> — a 10-bp core gate <code>TTCGTTGGAA{{s&lt;=3}}</code>
then the full box <code>[CT]TTCGTTGGAA[AG]CGGGA{{s&lt;=5}}</code>, both strands.)</p>
<p class="sub"><b>The CENP-B box is functional, not abundant.</b> Tellingly, the box is <i>not</i> the
most common sequence even in human α-satellite: the most frequent 17-mer there is
<code>CAAAAAGAGTGTTTCA</code> (~18% identity to the box) — a different conserved part of the 171-bp
monomer. The box is a specific, functionally conserved ~17 bp region, so it must be searched for
<i>as the box</i>; a de-novo "dominant motif" scan would miss it entirely. That is why the test
below scans for the best <i>box-matching</i> window at any position, rather than the satellite's
single most abundant motif.</p>

<h3>Unbiased de-novo test (no CENP-B seeding) — the decisive check</h3>
<p>The ±5-flank search is <b>seeded</b> on the CENP-B box, so it only ever collects windows that are
already ≤5 substitutions from it — it cannot tell whether a satellite's <i>own</i> motif is box-like
(indeed ~99% of the windows it finds sit at exactly 5 substitutions, the matching boundary). To
remove this bias we scan <b>every</b> 17-bp window of every array (both strands, no seeding), take
the single <b>best</b> match to the canonical box, and compare to the same on a <b>dinucleotide
shuffle</b> of the array (which controls for the "best of many windows" inflation). A satellite that
genuinely carries a box in its repeat unit — like α-satellite — scores well above its shuffle; an
unrelated satellite scores on the null.</p>
<div class="fig">{img(FIG/'cenpb_denovo_bestwindow.png',w="60%")}<div class="cap">
Best-window identity in the real satellite (y) vs a dinucleotide shuffle (x). Note the <b>null is
already high (~61%)</b> — scanning ~300 windows of <i>any</i> sequence finds a decent match to a
17-bp motif — so absolute values are close and the <b>excess</b> (distance above the diagonal) is
the signal. <b>α-satellite</b> (red star) is clearly above (78% vs 61%, <b>+17</b>) — a real box in
every monomer. <b>HSat</b> (−1) is on the null. DToL clades — including the <b>birds</b> (goshawk
+2.4, ptarmigan +0.2, takahē −3.8) — sit near the diagonal (median excess −0.8; 5/162 clear +10, all
well below α-sat). <code>cenpb_denovo_bestwindow.py</code></div></div>
<p class="sub"><b>Reading the two methods together.</b> Both are reported because they answer
different questions. The <b>seeded ±5-flank</b> search highlights a genuine, interesting observation:
several birds carry a satellite consensus that is <b>15/17 identical to the CENP-B box</b> (goshawk
`CTTTTTTGGAAACGGGA`, missing only the CpG). The <b>unbiased de-novo</b> test then asks whether that
resemblance exceeds chance <i>per array</i> — and it does not (the birds sit on the null; the 15/17
consensus arises because the box-matched windows are pre-selected and their random deviations
average out). So the honest reading: birds have a <b>candidate CENP-B-like motif worth following up</b>
(near-canonical consensus, box-shaped by flank Δ, TIGD-family binder plausible), but it is
<b>not confirmed</b> as a box-carrying satellite the way human α-satellite is.</p>

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
