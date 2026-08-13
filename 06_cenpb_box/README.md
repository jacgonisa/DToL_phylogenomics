# 06 — CENP-B box screen across the 325-species satellites

Screens the centromeric satellite repertoire of the DToL species for the **CENP-B
box** (17 bp, `[CT]TTCGTTGGAA[AG]CGGGA`; Masumoto et al. 1989), using **two
methods**, with human HG002 satellites as positive/negative controls.

## Which satellites were searched
- **Input:** `all.satellites.txt` — every **putative centromeric satellite** array (curated
  candidate centromeric satellites) across
  the DToL assemblies (one sequence per array, tagged by species).
- **Species set:** restricted to the **325 published species** (`species_325.txt`,
  derived from the calibrated species-tree tips). 12 extra species present in
  `all.satellites.txt` but absent from the final tree were excluded.
- **"Uncapped":** every satellite array of each species is searched (both strands),
  not a ≤500-monomer subsample — **20.3 M arrays / 3.58 Gbp across the 162 of 325
  species with satellite annotations**.

## Human HG002 benchmarks
- **Positive — α-satellite:** 50,000 monomers (171 bp; 8.55 Mbp) — carries the box.
- **Negative — HSat1/2/3:** 7,982 arrays (161.8 Mbp), extracted from `hg002v1.1`
  via the v1.1 CenSat annotation — human satellite that lacks the box.

Canonical boxes are ~50,000× denser per Mbp in α-sat than HSat, and both methods
separate the two controls cleanly (`cenpb_human_benchmark.py`,
`data/cenpb_human_benchmark.tsv`).

## Method 1 — exact IUPAC motif parsing (Fachinetti)
`cenpb_paper_motifs.py` (+ `cenpb_paper_motifs_plot.py`). Three **exact IUPAC**
motif tiers (both strands): canonical `YTTCGTTGGAARCGGGA`, broad `YTTCGNNNNANRCGGGN`,
degenerated `NTTCGNNNNANNCGGGN` (Barra & Fachinetti, bioRxiv 2026.05.25.727640).
Enrichment is obs/expected vs a **dinucleotide-preserving null** (first-order
Markov ≈ Altschul–Erikson doublet shuffle, `ae_shuffle.py`; a 0-order null is also
reported). **Result:** no DToL species has an exact *canonical* box; looser-tier
hits are trace and the plant broad signal deflates under the dinucleotide null
(11.7×→3.0×). *Figure:* `figures/cenpb_paper_motifs.png`.

## Method 2 — songbird ±5-flank test
`cenpb_flank_uncapped.py` (+ `cenpb_flank_uncapped_logos.py`). The box matched as
17-bp windows (≤5 substitutions, both strands), then the position frequency matrix
of the box **plus ±5 flanking bp**; a real motif has high box information that
collapses in the flanks (Δ = box − flank). Reports Δ, box consensus, subs from
canonical, prevalence (Formenti et al., *Cell* 2026, zebra-finch T2T).
**Result:** the diverged box signal is strongest in **birds** (goshawk Δ 0.85,
2 subs; ptarmigan 1 sub; takahē), but these lose the essential 5′ CG — a likely
TIGD4-type avian box, not a mammalian CENP-B box. *Figures:*
`figures/cenpb_flank_uncapped_scatter.png`, `figures/cenpb_box_logos_flanks_VERTEBRATES_uncapped.pdf`.

## Report & reviewer response
- **`cenpb_box_report.html`** — self-contained HTML report (figures embedded), built by
  `scripts/cenpb_report.py`.
- `RESPONSE_TO_REVIEWERS.md` — the two-method write-up with benchmarks and conclusions.

> **Terminology:** a CENP-B *box* is functional (binds protein). We detect sequence
> matches only → **motifs / candidate boxes**; none is a confirmed functional box.

## Data (`data/`)
`cenpb_human_benchmark.tsv`, `cenpb_paper_motifs_per_{species,clade}.tsv`,
`cenpb_flank_uncapped_per_species.tsv`, `species_325.txt`.
