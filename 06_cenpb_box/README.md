# 06 — CENP-B box screen across the 325-species satellites

Screens the centromeric satellite repertoire of each species for the **CENP-B
box** — the 17-bp motif bound by the CENP-B protein
(`[CT]TTCGTTGGAA[AG]CGGGA`, canonical; Masumoto et al. 1989) — and asks which
species **genuinely** carry it, mapped onto the calibrated chronogram. The
**human HG002 α-satellite** annotation is scanned throughout as the positive
benchmark.

> **Paths:** scripts carry absolute paths from the original analysis dir; edit
> the path variables to run elsewhere. Two satellite inputs are used: the curated
> monomer set (≤500 seqs/species) for the mismatch/labile titrations, and the
> full uncapped array table (`all.satellites.txt`, ~24.5M arrays) for the
> genome-wide exact scan.

## Three complementary detection methods

**1. Fachinetti / canonical mismatch** (`cenpb_box_scan.py`,
`cenpb_mismatch_titration_*.py`, `cenpb_uncapped_scan.py`)
The canonical 17-bp box, both strands, allowing an increasing number of
**mismatches (0→5)**, against a composition-**shuffled null**. This is the strict
functional-box test (Fachinetti-style). Titrations are provided for all species
and for vertebrates specifically, each with the human α-sat benchmark and null.

**2. PWM / log-odds** (`cenpb_box_finder.py` + FIMO; `cenpb_uncapped_scan.py`)
A position-weight matrix built from the human α-sat box instances, scanned as a
log-odds score with a 75%-of-max threshold. Catches diverged-but-real boxes that
a hard mismatch cap misses, while staying far stricter than the labile method.

**3. Songbird "labile" method** (`cenpb_songbird_labile.py`)
From Formenti et al., *Cell* 2026 (zebra finch T2T; Suppl. Fig. 15): the human
CENP-B box searched by **edit distance (Levenshtein — indels allowed, not just
substitutions)** with the paper's generous thresholds — **edit distance ≤5 for
cross-species**, **≤2 for the same species (human)** — against a shuffled null
(their ±5 flanking negative control). This is the deliberately *labile* axis for
CENP-B-**like** box candidates in non-mammals, where the box is diverged. Their
bird box sits in the Tgut716A satellite, with the pogo/Tigger transposase
**TIGD4** as the putative CENP-B replacement.

## Key results

**Exact functional box = human only.** The uncapped genome-wide scan (22.5M
arrays, `data/cenpb_uncapped_per_species.tsv`) finds the exact canonical box
(0 mismatches) **only in human** (13,359 arrays; benchmark). No non-human DToL
species carries a single exact canonical box.

| species | clade | exact (0 mm) | box ≤2 mm | PWM hits | best PWM |
|---|---|---|---|---|---|
| **human** (benchmark) | Human | 13,359 | 21,026 | 437 | 26.6 |
| *Rhinolophus sinicus* (mrhisin, bat) | Vertebrate | 0 | 390 | 8 | 24.0 |
| ibectpall | Invertebrate | 0 | 248 | 7 | 22.4 |
| all others | — | 0 | ≤4 | ≤2 | — |

The single **closest non-human** by every strict metric is the bat *Rhinolophus
sinicus* — the lone vertebrate approaching a canonical box, though still with
zero exact matches.

**Labile method: broad, low-baseline candidates.** With edit distance ≤5
(`data/cenpb_songbird_labile_per_species.tsv`), the human α-sat is a clean
benchmark (219 boxes at e≤2, null 0). Across non-human species, **38 show
labile CENP-B-like boxes above their shuffled null** (enrichment 2–7×), spread
across invertebrates and plants. These are diverged, indel-tolerant candidates
with a low null baseline — consistent with the paper's view that *labile*
CENP-B-like boxes are widespread, but they are **not** the strict functional box
(which remains human-only here).

## Figures
- `figures/cenpb_uncapped_per_species.tsv` scan → exact/PWM story (human-only box).
- `figures/cenpb_mismatch_titration_allspecies.{png,pdf}` and
  `..._vertebrates.{png,pdf}` — occurrences vs. mismatch tolerance, one line per
  species, human α-sat benchmark (bold red) + shuffled null (Fachinetti axis).
- `figures/cenpb_songbird_labile.{png,pdf}` — edit-distance sweep with the e≤2
  (same-species) and e≤5 (cross-species) thresholds marked, human benchmark, null.
- `figures/cenpb_box_occurrences_tree_325sp.{png,pdf}` — fan chronogram, tips by
  clade, outer bars = # satellite monomers carrying a CENP-B box (≤2 mm) per
  species (`scripts/plot_cenpb_tree_325sp.R`).
- `figures/cenpb_abundance_by_clade_325sp.{png,pdf}` — per-clade summary table.
- `cenpb_mammals*.png`, `cenpb_bat_logo.png`, `cenpb_fimo_logo.png` — mammal/bat
  motif logos.

## Data tables (`data/`)
- `cenpb_uncapped_per_species.tsv` — genome-wide exact + PWM per species.
- `cenpb_mismatch_titration_allspecies.tsv`, `..._vertebrates.tsv` — obs/null per
  species × mismatch level.
- `cenpb_songbird_labile_per_species.tsv` — labile obs/null/enrichment per species.
- `cenpb_abundance_per_species_325sp.tsv`, `cenpb_abundance_per_clade_325sp.tsv`
  — merged abundance tables.
