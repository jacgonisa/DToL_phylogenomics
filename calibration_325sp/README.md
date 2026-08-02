# 325-species DToL chronogram — calibration pipeline

Reproducible pipeline for dating the 325-species FastSpeciesTree (supermatrix ML
tree) and validating it against TimeTree.

> **Note on paths:** the scripts carry absolute paths from the original analysis
> directory (`/home/jg2070/Desktop/dtol_review_August/...`). To run them
> elsewhere, edit the `BASE`/`SP`/`SRC`/`ROOT` path variables at the top of each
> script to point at this `calibration_325sp/` folder. The `data/`, `trees/` and
> `figures/` here are the committed outputs of that pipeline.

## 0. Input tree
- `fast_species_tree_325sp_renamed.nwk` — ML supermatrix tree (325 tips, branch
  lengths in substitutions/site). Midpoint-rooted + `multi2di` inside the
  calibration script.

## 1. TimeTree reference data
```
python3 fetch_timetree_ci.py
```
Resolves each calibration node's representative species pair → NCBI taxid →
`timetree.org/api/pairwise/{a}/{b}/json`. Captures, per node:
- `tt_median`   = TimeTree "Median Time" (`precomputed_age`)
- `tt_adjusted` = TimeTree "Adjusted Time" (`adjusted_age`; 0 → none)
- `tt_ci_low/high` = TimeTree "Range" (`precomputed_ci_low/high`)
- `source` = TimeTree pairwise API | literature (2 nodes) | unresolved

Outputs: `outputs/calibration_qc/calibration_nodes_timetree.tsv`
(+ cache `outputs/calibration_qc/timetree_api_cache.json`).
NCBI synonym aliases handled in the `ALIAS` map (e.g. Potentilla anserina →
Argentina anserina).

## 2. Calibration constraints
`over_calib.tsv` — 64 nodes (`label, tip_a, tip_b, age_min, age_max`).
**Uniform rule:** constraint = TimeTree **CI range** where a real range exists
(ci_low>0 & ci_high>ci_low), else **median ± 20%**; 2 nodes are
literature-derived (no TimeTree): `Limnephilidae` (2026 Limnephilidae
phylogenomics study) and `Luzula_TT` (crown ~10 Ma, literature).
Rebuild the uniform constraints from the TimeTree table + a source of tip pairs;
current file already encodes the rule. Backups: `over_calib.tsv.bak_*`.

## 3. Date the tree (chronos, penalised likelihood)
```
Rscript calibrate_chronos_correlated_325sp.R over_calib.tsv outputs/full_325sp_chronos_over_correlated.nwk
```
chronos `lambda=1, model="correlated"`, high eval budget, **restarts until a
run converges** (breaks on first `converged=TRUE`; reports PHIIC). Node ages
resolved by `getMRCA(tip_a, tip_b)`; root pinned specially. `.fa` suffixes
stripped from output tips.

Other dating options (for the benchmark):
- `gen_ratesmoothed.R`     → `full_325sp_chronos_rootonly_ratesmoothed.nwk` (root-only PL)
- `gen_uncalib_parrett.R`  → `full_325sp_chronos_nocalib_relaxed.nwk` (no calibration, relaxed)
- chronos relaxed, treePL, RelTime trees also in `outputs/`.

## 4. QC vs TimeTree (PAReTT-style node-age concordance)
```
Rscript gen_parrett_all_325sp.R
```
Node-averaged MRCA age (our tree vs TimeTree, 210 shared species) for all four
dating options → `outputs/calibration_qc/parrett_{method}.tsv`.

## 5. Supplementary calibration table (all 64 nodes)
```
Rscript make_calibration_table_64_325sp.R
```
→ `outputs/calibration_table_325sp_supp_64.tsv` + `figures/…supp_64.{pdf,png}`
(Clade, Node, Taxa, constraint bounds, calibrated age, TimeTree median/adjusted/range, Source).

## 6. Benchmark validation figure
```
Rscript plot_calibration_combined_benchmark_publication_325sp.R
```
→ `figures/calibration_combined_qc_benchmark_325sp_publication.{pdf,png}`
- Panel A: 64 calibration nodes — red diamond = TimeTree median, blue diamond =
  adjusted, red line = TimeTree range, black circle = calibrated age; grey line
  + `*` = literature (not TimeTree); y-labels coloured by clade.
- Panels B–E: node-age concordance for chronos correlated / relaxed /
  rate-smoothed (root only) / uncalibrated.

## 7. Tree figures
`plot_tree_v5_chronos.R` / `plot_tree_v5_treepl.R` →
`figures/centromere_annotation_tree_FASTSPECIES_325sp_v1_{chronos,treePL}.{pdf,png}`.

## 8. Ancestral state reconstruction
Separate pipeline: `2026_trees/annotation_centromeres/ASR_March2026_327species/`.
Input trees `inputs/{full,metazoa,viridiplantae}_chronos_correlated/tree_renamed.nw`
are pruned from `outputs/full_325sp_chronos_over_correlated.nwk`. Run scripts
40 (model table) → 44 (formatted table), 43 (ASR tree plots).
