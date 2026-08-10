# CENP-B box search — text for response to reviewers

## How we searched
We screened the centromeric satellite repertoire of all 325 DToL species for the
17-bp **CENP-B box** (canonical consensus `[CT]TTCGTTGGAA[AG]CGGGA`; Masumoto et
al. 1989), on **both strands**, using three complementary methods of increasing
permissiveness, each against a composition-matched **shuffled null**:

1. **Strict / Fachinetti mismatch scan** — the canonical box allowing 0–5
   *substitutions*. Run both on a curated monomer set (≤500 monomers/species) and
   on the full uncapped array set (**22.5 million** satellite arrays).
2. **PWM log-odds** — a position-weight matrix built from the human α-satellite
   box instances, scanned at ≥75 % of the maximum score.
3. **Labile edit-distance method** (following Formenti et al., *Cell* 2026, zebra
   finch T2T) — the box searched by **edit distance (indels allowed)** ≤5 for
   cross-species and ≤2 within species. For every candidate we then built
   **sequence logos of the box ± 5 flanking bp**, exactly as in that paper: a
   genuine motif shows high information content across the 17-bp box that
   *collapses sharply to near-random in the flanks*. If the flanks are as
   conserved as the box, the "hit" is merely conserved satellite sequence, not an
   enriched CENP-B box (built-in negative control).

## What we found
- The **exact canonical box occurs only in human** α-satellite (positive control;
  13,359 of 22.5 M arrays). **No non-human DToL species carries a single exact
  box.**
- The closest non-human is the bat ***Rhinolophus sinicus*** (0 exact matches, but
  the highest ≤2-mismatch and PWM scores of any non-human species).
- The labile method flags 38 species above their shuffled null, but (i) these are
  **not overrepresented in vertebrates** — vertebrate mean enrichment is 0.89,
  *below* the null, versus ~1.4× in plants/invertebrates — and (ii) under the
  ±5-flank logo test only **human** (box − flank information Δ = **+0.64 bits**)
  and, secondarily, the **bat** (Δ = +0.46) show the expected box>flank
  transition. Most labile "candidates" have flanks as informative as the box
  (Δ ≈ 0 or negative) and collapse to a handful of gap-free windows — i.e. the
  labile counts were largely **indel-driven** and reflect conserved satellite
  sequence rather than a functional CENP-B box.

## Conclusion
A canonical, functional CENP-B box is confined to **human (mammalian)
α-satellite** in this dataset; the bat *R. sinicus* is the only non-human that
even approaches it. Outside mammals we find **no evidence for an enriched CENP-B
box motif** — consistent with CENP-B being a mammal-specific feature.

## Figures / data referenced
- `figures/cenpb_mismatch_titration_allspecies.png` — occurrence vs. mismatch
  tolerance, all species, human α-sat benchmark + shuffled null (Fachinetti axis).
- `figures/cenpb_songbird_labile.png`, `figures/cenpb_songbird_labile_byclade.png`
  — labile edit-distance method and its per-clade breakdown (no vertebrate excess).
- `figures/cenpb_box_logos_flanks_all_candidates.pdf` — box ±5-flank sequence
  logos for **all** candidates (page 1 = box-vs-flank information scatter; only
  human and the bat sit clearly above the diagonal).
- `data/cenpb_uncapped_per_species.tsv` — genome-wide exact + PWM counts.
- `data/cenpb_box_logos_flanks_summary.tsv` — per-species box vs flank information.
