# 06 — CENP-B box screen across the 325-species satellites

Screens the centromeric satellite repertoire of each species for the **CENP-B
box** — the 17-bp motif bound by the CENP-B protein
(`[CT]TTCGTTGGAA[AG]CGGGA`, canonical; Masumoto et al. 1989) — and quantifies its
abundance per species and per clade, mapped onto the calibrated chronogram.

> **Paths:** scripts carry absolute paths from the original analysis dir; edit the
> path variables to run elsewhere. Satellite input = curated monomer set
> (≤500 seqs/species). 175 species have a scanned satellite set; 163 are also
> tips in the 325-sp tree.

## What was done — two complementary detection methods

**1. Mismatch scan** (`scripts/cenpb_box_scan.py`)
Scans every satellite monomer (both strands) for the canonical 17-bp box allowing
**≤2 and ≤3 mismatches** (boxes are diverged outside mammals), against a
composition-**shuffled null**. Reports, per species: `pct_box_le2mm`,
`pct_box_le3mm`, the null rate, and enrichment (box/null).

**2. Lines-of-evidence** (`scripts/cenpb_box_finder.py`)
Four independent detectors on 2000 sampled monomers, each with its own shuffled
null: **canonical** IUPAC exact (`YTTCGTTGGAARCGGGA`), **broad** IUPAC
(`NTTCGNNNNANNCGGGN`), **degenerate** IUPAC, and a data-driven **PWM** scanned
with FIMO. A species is called box-positive when ≥2 lines exceed their null.

**Focused mammal / bat analyses** (`cenpb_mammals*.py`, `cenpb_bat_motif.py`):
exact + mismatch + core-motif screen of the mammals (where CENP-B boxes are
canonical) and a bat-specific PWM/logo.

## Abundance table (`data/` + `figures/`)
- `cenpb_abundance_per_species_325sp.tsv` — merged per-species table, all methods
  (n_monomers, %box ≤2mm/≤3mm, null, enrichment, box count, canonical/broad/
  degenerate/PWM-FIMO counts, lines_positive).
- `cenpb_abundance_per_clade_325sp.tsv` + `figures/cenpb_abundance_by_clade_325sp.{png,pdf}`
  — per-clade summary rendered as a table.

| clade | n sp | mean %box ≤2mm | mean %box ≤3mm | enrichment | canonical | broad | degen. | PWM/FIMO |
|---|---|---|---|---|---|---|---|---|
| Vertebrates | 24 | 0.041 | 1.38 | 1.38 | 0 | 21.0 | 27.0 | **35.8** |
| Invertebrate | 102 | 0.031 | 0.83 | 2.35 | 0 | 0.6 | 3.6 | 0 |
| Viridiplantae | 49 | 0.002 | 0.28 | 0.84 | 0 | 2.9 | 5.7 | 0 |

**Key result:** the stringent PWM/FIMO detector finds the CENP-B box almost
exclusively in **Vertebrates** (mainly mammals; PWM 35.8 per 2000 vs 0 elsewhere)
— the classic functional CENP-B box. The looser ≤3-mismatch scan yields a higher
*apparent* enrichment in invertebrates, but these are short degenerate matches
with a low shuffled-null baseline, not the canonical functional box.

## Tree figure
`figures/cenpb_box_occurrences_tree_325sp.{png,pdf}` — fan chronogram, tips
coloured by clade, outer bars = number of satellite monomers carrying a CENP-B box
(≤2 mismatches) per species (`scripts/plot_cenpb_tree_325sp.R`).

## Other figures
`cenpb_box_by_clade.png` (scan overview), `cenpb_mammals_logos.png` / `cenpb_mammals.png`
(mammal motif + logos), `cenpb_bat_logo.png`, `cenpb_fimo_logo.png`.
