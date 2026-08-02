# 01 — Species tree: building, calibration & plotting (325 species)

Builds the 325-species **FastSpeciesTree**, time-calibrates it with `chronos`,
validates against TimeTree, and renders the annotated tree figures.

> **Note on paths:** the R/python scripts carry absolute paths from the original
> analysis directory (`/home/jg2070/Desktop/dtol_review_August/...`). To run them
> elsewhere, edit the `BASE`/`SP`/`SRC`/`ROOT` path variables at the top of each
> script to point at this `01_species_tree/` folder. The `data/`, `trees/` and
> `figures/` here are the committed outputs.

## 0. Tree building — FastSpeciesTree (supermatrix, IQ-TREE)

The species tree is a **supermatrix ML tree**, NOT the BUSCO strategy in the
repo root README. It was built with the `FastSpeciesTree` pipeline (mode
SENSITIVE): BLAST-based single-copy gene selection across the annotated
proteomes → per-gene MAFFT alignment → concatenated pseudo-alignment (+ IQ-TREE
partition file) → IQ-TREE with per-partition model selection. Full run log:
[`tree_building/fast_species_tree_iqtree.log`](tree_building/fast_species_tree_iqtree.log).

The tree-inference command (both untrimmed and trimmed pseudo-alignments):
```
iqtree -T 32 \
  -s Results/psuedo_alignment.fasta \
  -p Results/IQTree_Partition_file.partitions \
  -B 1000 --alrt 1000 -st AA \
  -mset LG,JTT,Q.BIRD,Q.MAMMAL,Q.INSECT,Q.PLANT,Q.YEAST \
  -mrate I,G,I+G -m MFP \
  --prefix Results/fast_species_tree_iqtree
```
Result: `trees/fast_species_tree_325sp_renamed.nwk` — 325 tips, branch lengths in
substitutions/site (midpoint-rooted + `multi2di` inside the calibration script).

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
- Panel A: 64 calibration nodes — red diamond = TimeTree median, red line =
  TimeTree range (= the imposed constraint for nodes with a real range, else
  median ±20%), black circle = calibrated age; grey line + `*` = literature
  (not TimeTree); y-labels coloured by clade.
- Panels B–E: node-age concordance vs TimeTree (Spearman + R², 1:1 line) for
  chronos correlated / relaxed / rate-smoothed (root only) / uncalibrated.

## 7. Annotated tree figures
```
Rscript plot_tree_v5_chronos.R     # chronos-correlated chronogram (Figure 1 style)
Rscript plot_tree_v5_treepl.R      # treePL chronogram
```
Fan chronogram coloured by clade + centromere architecture (Satellite /
Transposon / Mixed / Holocentric). Inputs: the calibrated `*_fa.nwk` tree
(the calibration script writes a `.fa`-tipped sibling automatically),
`data/DTOL_327_master_March.xlsx` (species metadata) and
`data/branch_symbol_anno.tsv` (iTOL-style architecture symbols).
Outputs: `figures/centromere_annotation_tree_FASTSPECIES_325sp_v1_{chronos,treePL}.{pdf,png}`.

## 8. Ancestral state reconstruction
Separate pipeline: `2026_trees/annotation_centromeres/ASR_March2026_327species/`.
Input trees `inputs/{full,metazoa,viridiplantae}_chronos_correlated/tree_renamed.nw`
are pruned from `outputs/full_325sp_chronos_over_correlated.nwk`. Run scripts
40 (model table) → 44 (formatted table), 43 (ASR tree plots).
