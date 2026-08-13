# Handoff: CENH3/H3 entropy-MSA / GroupSim / terminal-branch analyses

For a fresh Claude session picking up the CENP-A/CENH3 vs H3 comparative work
in the 325-genome DToL centromere paper. Everything below is current state.

## THE ONE NUMBER THAT MATTERS
Final curated groups: **422 CENP-A/CENH3  vs  897 H3-like**.
Source of truth (both files identical split, different gap trims):
`04_cenpa_phylogeny/split_entropy/groupsim/groups_gap08.txt`  (gap 0.80)
`04_cenpa_phylogeny/split_entropy/groupsim/groups_gap085.txt` (gap 0.85)
Format: two lines, `CENPA:id1,id2,...` and `H3:id1,id2,...`

### How 422/897 was reached (do NOT recompute from scratch)
- Started 418 CENP-A / 901 H3 (from the IQ-TREE input fasta).
- Inspected the final VT+G4 tree in iTOL; moved **6 H3-labelled seqs that fall
  inside the CENP-A clade → CENP-A**: tzMemMemb1_chr_1_002492.1,
  gfChlBrun_ENA|OY756145|OY756145.1_000837.1, qeDicMinu_ENA|OX461808|OX461808.1_020535.1,
  ieEcdTorr_ENA|OX439129|OX439129.1_002001.1, icPyrSerr1_chr_X_000403.1,
  icHalSede1_chr_9_000164.1.
- Removed 2 that did NOT cluster with CENP-A → H3: idSicFerr1_chr_X_002574.1,
  idSicFerr1_chr_X_002608.1.
- Net: 418−2+6 = **422 CENP-A**, 901+2−6 = **897 H3**.

## THE TREE
Final ML tree: `iqtree2_cenpa430_H3_archaea10_325sp_VTG4_bnni.{treefile,contree}`
in `04_cenpa_phylogeny/`.
- Model **VT+G4** (selected by ModelFinder/BIC on the actual 325sp alignment;
  replaced the previously hard-coded Q.pfam+R7 which was never justified here).
- Command: `iqtree2 -m VT+G4 -B 1000 -bnni -nm 2000 -T 6`. Converged iter 404/505.
- Alignment: `cenpa430_H3_archaea10.aligned.clipkit.325sp.fasta`
  (418 CENP-A + 901 H3 + 10 archaea = 1329 tips). Archaea IDs list is in every
  script (OLS22332.1, OLS24873.1, OLS21974.1, KKK41979.1, KXH71038.1, OLS18261.1,
  OLS16336.1, BAD86478.1, OIO61677.1, OIO41945.1).

## GROUPSIM (specificity-determining positions)
Scripts in `04_cenpa_phylogeny/`:
- `run_groupsim_325sp.py` — unweighted GroupSim, gaps 0.80/0.85/0.90.
- `run_groupsim_cenpa_h3_clade_325sp.py` — clade-weighted (weight = 1/n_clade per
  broad taxon group; NOT Henikoff-Henikoff — HH was tried then dropped).
- `run_groupsim_sat_trans_325sp.py` + `run_groupsim_sat_trans_weighted_325sp.py`
  — Satellite vs Transposon CENP-A comparison.
- Plot scripts: `plot_groupsim_pub.R` (main, gap085), `plot_groupsim_cenpa_h3_gap085.R`,
  `plot_groupsim_st_pub.R` (+ `plot_groupsim_st_gap08_pub.R`).

### CRITICAL gotcha — group source
All the run_* scripts were patched to read the curated `groups_gap{tag}.txt`
FIRST, and only fall back to inferring H3 vs CENPA from `H3_all.aligned.fasta`
if that file is missing. If you re-add sequences, edit the groups_gap files and
rerun; do not trust the fasta-based inference (it gives the old 418/901).
Same fix was applied to `split_entropy_325sp.py` and `plot_panel_D_nomasking_325sp.py`.

### Current results
- CENP-A vs H3: **6 significant positions (z≥2)** at both gap 0.80 and 0.85. Robust.
- Satellite vs Transposon: **0 significant positions** at both thresholds.
- Holo vs Sat vs Trans (3-group, supervisor question): 7 sig positions but mostly
  monocot-conservation bias (all 20 holo CENP-A are Cyperaceae/Juncaceae). Output:
  `split_entropy/groupsim_holo_vs_mono/`.
- Plots (n=422/897 in titles): `figures/cenpa422/groupsim_*_pub.{pdf,png}`.

## ENTROPY / MSA
- `split_entropy_325sp.py` → split-entropy + per-group + MSA-matrix plots for
  gaps 0.80/0.85/0.90. Reads curated groups file (patched).
- Publication panels in `03_entropy/`: `plot_entropy_325sp.R` (panel C),
  `plot_panel_D_nomasking_325sp.py` (panel D). **User strongly prefers the
  no-masking + filled panel D** (`panel_D_split_entropy_325sp_nomasking.{pdf,png}`),
  NOT the masked version. Panel D reads counts dynamically from groups_gap085.txt.
- Diagnostic plots for the 422 split live in `03_entropy/figures/cenpa422/`.

## TERMINAL BRANCH LENGTH (CENP-A vs H3 evolutionary rate)
`plot_terminal_branch_325sp.R` — preferred over root-to-tip (RTD conflates a
tip's own rate with clade depth; terminal branch = lineage-specific rate).
- Uses the VT+G4 **.contree**, groups_gap08.txt.
- Result: **249/262 paired species have CENP-A evolving faster than H3**,
  paired Wilcoxon p≈1.2e-40. Only 12/262 the other way.
- y-axis clamped `limits=c(0,NA)` — the earlier "whiskers below 0" was a ggplot
  boxplot artefact, not negative branch lengths (min TBL is exactly 0).
- Style: violin + boxplot + median crossbar + grey paired-species connectors,
  matching `plot_root_to_tip_325sp.R`. Output `figures/cenpa422/terminal_branch_*`.
- Note: 606 tips have zero terminal branch (near-identical H3 paralogs) — expected.

## OTHER CONTEXT
- *Eimeria praecox* (pxEimPrae) is the one CENP-A species with no detectable H3
  (Apicomplexa histones too divergent for the H3 HMM even at E=10) — excluded
  from paired comparisons. 257–262 paired species depending on group set.
- iTOL annotation files for validating the tree: `iTOL_VTG4_416/` and
  `iTOL_VTG4_bnni/` (colour strips for architecture, CENP-A presence, taxonomy).

## CONVENTIONS
- Publication figures for the 422/897 split go in `.../figures/cenpa422/`
  subfolders so the earlier 418/901 versions aren't overwritten.
- Colours: CENP-A `#C62828` (red), H3 `#1565C0` (blue); architecture Sat `#E53935`,
  Trans `#FB8C00`, Holo `#2d7d32`.
