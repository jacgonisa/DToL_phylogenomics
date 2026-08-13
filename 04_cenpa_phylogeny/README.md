# 04 — CENP-A/CENH3 phylogeny & GroupSim (CENP-A vs H3, 325 species)

Phylogeny of the centromeric histone **CENP-A/CENH3** together with canonical
**H3** across the DToL species, and a **GroupSim** analysis of the
**specificity-determining positions** that distinguish CENP-A/CENH3 from H3.

## Alignment & tree
- **Alignment:** `data/cenpa430_H3_archaea10.aligned.clipkit.325sp.fasta` — CENP-A/CENH3
  + H3 + 10 archaeal histone outgroups (**1,329 tips** = 418 CENP-A + 901 H3 + 10
  archaea), MAFFT-aligned and ClipKit-trimmed.
- **ML tree:** IQ-TREE, model **VT+G4** (ModelFinder/BIC), `-B 1000 -bnni`:
  `data/iqtree2_cenpa430_H3_archaea10_325sp_VTG4_bnni.{treefile,contree,iqtree}`.
- **Curated groups:** **422 CENP-A/CENH3 vs 897 H3-like** — starting from
  418/901 the final VT+G4 tree was inspected (iTOL) and 6 H3-labelled sequences that
  fall inside the CENP-A clade were moved to CENP-A, 2 non-clustering to H3
  (`data/groupsim/groups_gap08.txt`, `groups_gap085.txt`; see
  `data/HANDOFF_CENH3_H3_groupsim_entropy.md` for the exact edits).
- Plot: `scripts/plot_cenpa_phylogeny_325sp.R` → `figures/panel_D_cenpa_phylogeny_325sp.{png,pdf}`;
  phylogenetic signal / root-to-tip: `scripts/plot_cenpa_phylosignal.R`.

## GroupSim — specificity-determining positions
`scripts/run_groupsim_325sp.py` (+ `_weighted`, `_cenpa_h3_clade`) computes, per
alignment column, the **GroupSim** score separating CENP-A/CENH3 from H3 (gap
thresholds 0.80 / 0.85 / 0.90; unweighted and clade-weighted = 1/n per broad taxon,
not Henikoff–Henikoff). High-scoring columns are the residues that specify CENP-A
identity (concentrated in the histone-fold, mapped onto the α-helices).
- Results: `data/groupsim/groupsim_gap0{8,85}.txt`.
- Figures: `figures/groupsim_cenpa_vs_h3_pub.{png,pdf}`,
  `figures/groupsim_cenpa_vs_h3_gap085_pub.{png,pdf}`,
  `figures/groupsim_manhattan_v2.{png,pdf}`,
  `figures/groupsim_gap085_with_helices.{png,pdf}`, `figures/groupsim_diagram.{png,pdf}`.
- **Variants** (same method, different groupings): satellite vs transposon
  centromeres (`run_groupsim_sat_trans_325sp.py`, `figures/groupsim_sat_vs_trans_*`),
  per-phylum (`run_groupsim_phylum_325sp.py`), holocentric vs monocentric, vertebrate,
  and within-clade (`plot_groupsim_within_clade.py`).

## Split entropy
`scripts/split_entropy_325sp.py` — positional Shannon entropy computed **separately**
for CENP-A/CENH3 vs H3 (the split-entropy view; complements the GroupSim scores and
`03_entropy`).

## CENP-A copy number
`scripts/cenpa_copy_number_analysis.py` + `plot_cenpa_copy_number.py` — per-species
CENP-A/CENH3 gene copy number and duplications, mapped on the species tree.
- Data: `data/cenpa_copy_number_table.tsv`, `cenpa_copies_by_order.tsv`,
  `cenpa_gene_duplications.tsv`, `cenpa_blast_vs_refs_contree.tsv`.
- Figures: `figures/cenpa_copies_species_tree.{png,pdf}`, `cenpa_copy_number_bar.{png,pdf}`.

> **Not included here:** IQ-TREE intermediates (`.ckp.gz`, `.ufboot`, `.splits.nex`,
> ModelFinder files) and the separate CENH3-antibody design collaboration — kept out
> of the repository to avoid bloat.
