# CENP-B box search — response to reviewers

We screened the centromeric satellite repertoire of the 325 DToL species for the 17-bp
**CENP-B box** (canonical consensus `[CT]TTCGTTGGAA[AG]CGGGA`; Masumoto et al. 1989), on
both strands, across all putative centromeric satellite arrays (162 of 325 species carry
satellites), using two complementary approaches. Human HG002 α-satellite was the positive
control and human HSat the negative control.

## Method 1 — exact IUPAC motif parsing (Fachinetti)
We counted three **exact IUPAC** definitions of the box — canonical `YTTCGTTGGAARCGGGA`, **broad**,
and **degenerate** (Barra & Fachinetti, bioRxiv 2026.05.25.727640) — per species across the full
uncapped set of ~24.5 million satellite arrays (both strands), scoring each as **enrichment =
observed / expected** against a **dinucleotide-preserving null**. A genuine box requires *both* a
high motif density (hits/Mbp) *and* high enrichment.

**Result.** The exact canonical box occurs **only in human α-satellite** (positive control; 13,359
hits in the 50,000-monomer HG002 α-sat set). **No DToL species carries a single exact canonical box**
(0 in all 162 species with satellites), and the **broad and degenerate** tiers are likewise
essentially absent — at or below the null in every clade. Only α-satellite is high on *both* density
and enrichment; every DToL clade sits at the null.

## Method 2 — songbird ±5-flank approach (Formenti et al., Cell 2026)
Following the zebra-finch T2T paper (Suppl. Fig. 15), the box was matched as fixed 17-bp windows
(≤5 substitutions, both strands) and, for every species on the full ~24.5 M arrays, we built the
position-frequency matrix of the box plus ±5 flanking bp. A genuine motif shows high information
across the 17-bp box that collapses in the flanks (Δ = box − flank information > 0); we also report
the box consensus, its substitutions from canonical, and the prevalence (boxes/Mbp).

**Result (seeded search).** Among the box-seeded hits, birds give the strongest signal: the goshawk
carries a satellite consensus **15/17 identical** to the canonical box — differing only at the central
CpG — with the highest box-vs-flank Δ (0.85) in the 325-species dataset; ptarmigan and takahē follow.

**Unbiased de-novo check.** Because the ±5-flank search is *seeded* on the box (it only ever collects
windows already ≤5 substitutions away), we added an unbiased test: scan every window of every array
(no seeding), take the best match to the box, and compare to a dinucleotide-shuffle of the *same*
satellite. Human α-satellite scores well above its shuffle (a real box in every monomer); by
contrast **every DToL clade — birds included — sits on the null.** So the near-canonical bird
consensus, while striking, is a *suggestive candidate* that does not exceed chance the way human
α-satellite does.

## Summary
The canonical, functional CENP-B box is confined to **human/mammalian α-satellite**. In DToL we find
**no confirmed CENP-B box**: the exact box (and its broad/degenerate variants) is absent, and the
relaxed ±5-flank hits do not exceed a composition-matched null. The one signal worth following up is
in **birds** — a diverged, box-shaped satellite motif 15/17 identical to the box but lacking the
central CpG (goshawk, ptarmigan, takahē) — a plausible **candidate** CENP-B-like / TIGD-family box
that we **flag but do not confirm** (protein binding not tested; it does not beat the de-novo null).
This is consistent with, but does not yet establish, the recently reported avian CENP-B-like system.

## Suggested figures (vector PDFs, Inkscape-editable)
- Method 1: **`cenpb_paper_motifs.pdf`** — enrichment vs motif density; only α-satellite is a real
  box (top-right), all DToL clades sit near the null.
- Method 2: **`cenpb_denovo_bestwindow.pdf`** — best box-like window vs a shuffled null; α-satellite
  is above the diagonal, all DToL (birds included) sit on it.
- Overview: **`cenpb_box_tree_325sp.pdf`** — the signal mapped onto the chronogram.
- Optional (the bird candidate): **`cenpb_flank_uncapped_scatter.pdf`** — box vs flank information,
  showing the goshawk's near-canonical motif.

*(Full methods, all figures and per-species data: `cenpb_box_report.html`.)*
