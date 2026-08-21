# CENP-A / CENH3 protein-similarity decay & half-life

The satellite-DNA `seqsim` / half-life analysis (`05_satellite_similarity/`)
replicated for the **CENP-A protein**: how fast does CENH3 diverge with species
divergence time, and is it lineage-specific?

Run: `python3 ../../scripts/cenh3_halflife_325sp.py` (`--selftest` for the unit check).

## Method (mirrors the satellite pipeline)
- Pairwise **protein % identity** read straight from the curated 422-seq CENP-A
  MSA (`data/cenpa430_H3_archaea10.aligned.clipkit.325sp.fasta`); no BLAST needed.
  Species-pair score = **best copy-pair identity** (the orthologous CENH3 pair),
  over columns where both copies have a residue (≥30 shared columns).
- Divergence time = MRCA age in the calibrated chronogram
  (`full_325sp_calibrated_correlatedlambda01.nwk`), `mya = distance/2`.
- One point per MRCA (node-averaged), within-clade pairs only.
- Fit `H = A·e^(−λt) + C` (floor `C` free); `t½ = ln2/λ`, per clade.

## Result
| Clade | half-life | floor | note |
|---|---|---|---|
| Viridiplantae | **17.5 My** | ~63% | fastest; clean plateau |
| Invertebrates | ~101 My | ~45% | shallow, scattered |
| Vertebrates | **> ~440 My** | — | still declining at deepest node; floor not reached, t½ not identifiable (do not quote 844 My) |

**Takeaways**
1. **Lineage-specific, strongly.** Plant CENP-A turns over ~6× faster than
   invertebrate and ≫ vertebrate — the opposite of centromeric *satellite* DNA,
   whose half-life was ~7–10 My and roughly uniform across all three clades.
   So the CENP-A protein diverges far slower than the satellite it binds, and
   its rate is clade-dependent.
2. **High floor = conserved histone fold (CATD).** Even between distant species
   CENP-A stays ~45–63% identical; the "half-life" measures decay of the
   variable regions (N-tail, loop1) above that floor.
3. **Vertebrate CENP-A is the most conserved** here — no plateau within ~450 My,
   so its half-life is only a lower bound. Interpret the vertebrate curve as
   "slow, still declining", not a fitted t½.

## Files
- `cenh3_seqsim_pairs_325sp.tsv` — per species-pair identity + MRCA age + clade
- `cenh3_halflife_325sp.tsv` — per-clade fit summary (incl. `floor_reached` flag)
- figure: `../../figures/cenh3_halflife_325sp.{png,pdf}`
