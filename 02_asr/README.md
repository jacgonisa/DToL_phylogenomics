# 02 — Ancestral state reconstruction of centromere architecture (325 species)

Fits Markov (Mk) models of centromere-type evolution on the calibrated
chronogram and tests the α–ω cyclical hypothesis (can centromere type cycle
between Satellite- and Transposon-based states?). Run on the **chronos-correlated**
tree from `01_species_tree/`.

States: **H** Holocentric · **Mixed** (Sat/Transposon) · **Sat** Satellite · **Trans** Transposon.

> **Paths:** scripts use `getwd()` = parent of the ASR root and carry absolute
> paths from the original analysis dir; edit accordingly to run elsewhere.

## Inputs (`inputs/{full,metazoa,viridiplantae}_chronos_correlated/`)
- `tree_renamed.nw` — the chronos-correlated chronogram, pruned to each clade.
- `branch_symbol_anno.tsv` — per-species architecture annotation.
Built by `scripts/00_prepare_inputs.R`.

## Pipeline (`scripts/`, run in order)
- `37_test_cyclical_ARD_irrevH.R` — fit ER/ARD/ARD_irrevH/TerminalTrans, LRT, stochastic maps
- `38_find_reversals.R`, `39_independent_cycles.R` — detect X→Y→X reversals, count independent cycles
- `40_all_models_table.R` — AICc table for ER/SYM/ARD/ARD_irrevH/TerminalTrans (+ `44_model_table.R` formats it)
- `41_Qmatrix_plots.R` — rate-matrix visualisations
- `43_cycles_tree_mk_parsimony.R` — Mk (phytools `ancr`) + Fitch parsimony ASR on the tree
- `02_custom_models.R` — design matrices (incl. `ARD_irrevH`: H is an irreversible sink)

## Results on the converged chronos-correlated tree (`outputs/`)
Model comparison (Akaike weights):

| Dataset | Best model | weight | note |
|---|---|---|---|
| Full tree | **ARD_irrevH** | 0.76 | H irreversible |
| Metazoa | **ARD_irrevH** | 0.97 | H irreversible |
| Viridiplantae | **SYM** | 0.59 | ARD_irrevH within ΔAICc 1.7 (comparable) |

Cycles (Mk ASR):
- **Metazoa**: 1 broad Sat→Trans→Sat (Hymenoptera/Diptera/Coleoptera/…, 67 tips) + 1 Trans→Sat→Trans (12 tips).
- **Viridiplantae**: 2 Trans→Sat→Trans (Monocots+Dicots; Bryophyta), 0 Sat→Trans→Sat.

Figures: `outputs/all_models/model_comparison_table.{png,pdf}`,
`outputs/cycles_mk_parsimony/{dataset}_chronos_correlated_{mk,parsimony}_{rectangular,fan}.{png,pdf}`.
