# 05 — Satellite similarity decay & half-life (325 species)

How fast centromeric satellite sequences diverge, as a function of species
divergence time (from the calibrated chronogram in `01_species_tree/`).

> **Paths:** scripts carry absolute paths from the original analysis directory;
> edit the `ROOT`/`PUB` variables to run elsewhere.

## 1. Pairwise satellite similarity (BLASTN all-vs-all)
```
python3 scripts/seqsim_blastn_melters_325sp.py
```
BLASTN all-vs-all of the sampled satellite monomers between species (Melters
et al. 2013 style). Output: `data/seqsim_blastn_melters_325sp.tsv`
(`spA, spB, mya, mean_pct_id, frac_sig, n_hits, group`).

## 2. Decay curve & half-life
```
python3 scripts/halflife_chronos_correlated_325sp.py
```
- Divergence time `mya` recomputed per pair as the MRCA age in the
  **chronos-correlated** chronogram (`01_species_tree/.../full_325sp_chronos_over_correlated.nwk`).
- Node-averaged (one point per MRCA), then fit `H = A·exp(−λt) + C`
  (empirical floor `C` free); half-life `t½ = ln(2)/λ`.

Output: `figures/seqsim_halflife_chronos_correlated_325sp.{png,pdf}` +
`data/halflife_chronos_correlated_325sp.tsv`.

## Result (converged tree)
| Clade | half-life | floor |
|---|---|---|
| Vertebrates | 10.2 My | ~29% |
| Invertebrates | 8.3 My | ~31% |
| Viridiplantae | 6.9 My | ~31% |

Satellite identity decays with a half-life of **~7–10 My** across all clades.

**Note:** an earlier "oldest node still ≥60% identity" marker was removed — it was
driven by over-dated shallow nodes (e.g. *Schoenoplectus lacustris × tabernaemontani*
is dated 17.6 My in-tree but ~3.4 My in TimeTree), so it overstated satellite
persistence.
