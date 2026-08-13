# 03 — Positional entropy of the CENP-A / CENH3 alignment (325 species)

Quantifies **per-position sequence conservation** of the centromeric histone
CENP-A/CENH3 across the DToL species, as the positional **Shannon entropy** of the
protein multiple-sequence alignment. Produces the entropy panels of Figure 1.

## Inputs
- CENP-A / CENH3 protein sequences from the DToL species, aligned with **MAFFT**
  and curated (the same alignment underlying `04_cenpa_phylogeny`).
- The calibrated species tree (`01_species_tree`) for the on-tree plots.

## What was done
**Panel C — positional Shannon entropy** (`scripts/plot_entropy_325sp.R`,
`scripts/entropy_analysis.Rmd`)
Per alignment column, H = −Σ pᵢ log pᵢ over the residue frequencies. Low entropy =
conserved position (e.g. the histone-fold core), high entropy = variable position
(e.g. the N-terminal tail). Provided **raw**, **smoothed**, and **combined**
versions (`panel_C_entropy_{raw,smooth,combined}_325sp`).

**Panel D — split entropy** (`scripts/plot_panel_D_nomasking_325sp.py`,
`scripts/plot_panel_D_masked_325sp.py`)
Positional entropy computed **separately for CENP-A/CENH3 vs H3-like** sequences,
highlighting positions that are conserved within each group but differ between them
(the specificity-determining regions expanded in `04_cenpa_phylogeny` GroupSim).
> **Use the no-masking, filled version** (`panel_D_split_entropy_325sp_nomasking`)
> as the publication panel; `_masked_fills` and the base version are alternatives.

**On-tree views** (`scripts/plot_kmer_on_tree.R`, `scripts/plot_seqsim_on_tree.R`)
Map per-species k-mer / sequence-similarity summaries onto the chronogram.

## Figures (`figures/`)
- `panel_C_entropy_{raw,smooth,combined}_325sp.{png,pdf}` — positional entropy.
- `panel_D_split_entropy_325sp_nomasking.{png,pdf}` — split entropy (publication).
- `panel_D_split_entropy_325sp{,_masked_fills}.{png,pdf}` — alternatives.
- `panel_CD_entropy_combined_325sp.{png,pdf}` — combined C+D panel.
