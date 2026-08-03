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

## 3-state model set (H / Sat / Trans) — full 8-model comparison

The 4-state table above (`40_all_models_table.R`) fits only 5 models. For the
paper we use the **3-state** analysis: Mixed (Satellite/transposon, n=11) and
Unknown/Monocentric-unknown taxa are pruned to NA, leaving H/Sat/Trans, and the
**full 8-model set** is fit:
ER, SYM, ARD, ARD_irrevH, ARD_irrevH_noDirectST, ARD_irrevH_symST,
ARD_irrevH_noSatToTrans, ARD_irrevH_noTransToSat.

Script: `scripts/40b_all_models_3state.R` (reads this repo's own
`inputs/{ds}_chronos_correlated/`; AICc + AIC Akaike weights).

### Results on the chronos-correlated tree (`outputs/all_models_3state/`)
| Dataset | Best (AICc) | w | note |
|---|---|---|---|
| Full tree | **ARD_irrevH_symST** | 0.39 | ARD_irrevH within ΔAICc≈0.4 |
| Metazoa | **ARD_irrevH** | 0.63 | H irreversible |
| Viridiplantae | **SYM** | 0.55 | ARD_irrevH_symST within ΔAICc 1.9 |

The two one-directional models (noSatToTrans / noTransToSat) always carry
≈0 weight — evidence for **bidirectional Satellite<->Transposon cycling**.
Figures: `outputs/all_models_3state/all_models_3state_weights.{png,pdf}`.
