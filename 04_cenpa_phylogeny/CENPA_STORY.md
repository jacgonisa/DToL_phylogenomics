# CENP-A / CENH3 story — onboarding for a fresh Claude session

Read this first if you're picking up the **CENP-A/CENH3 phylogenetic** work in the
325-genome DToL centromere paper. It's the narrative + the decisions that aren't obvious
from the code. For the module file-map see `README.md`; for the detailed current state
see `data/HANDOFF_CENH3_H3_groupsim_entropy.md`. This doc ties them together.

Everything lives in `04_cenpa_phylogeny/` (repo: `github.com/jacgonisa/DToL_phylogenomics`,
branch `main`). Commit as **jacgonisa**, no Claude attribution.

---

## TL;DR — the numbers that matter
- **422 CENP-A/CENH3  vs  897 H3-like** sequences (final curated split). This is the
  source of truth for every downstream analysis. Files:
  `data/groupsim/groups_gap085.txt` and `groups_gap08.txt` (two lines:
  `CENPA:id1,id2,...` / `H3:id1,id2,...`; identical split, different gap-trim of the MSA).
- **Tree:** `data/iqtree2_cenpa430_H3_archaea10_325sp_VTG4_bnni.{treefile,contree,iqtree}`
  — 1,329 tips (422 CENP-A + 897 H3 + 10 archaeal-histone outgroups), model **VT+G4**.
- **Headline results:** CENP-A vs H3 has **6 specificity-determining positions** (GroupSim,
  robust); **CENP-A evolves faster than H3 in 249/262 species** (terminal branch length,
  Wilcoxon p≈1e-40); the specificity residues sit in the **CATD / histone-fold α-helices**.

---

## The story, in order

### 1. Sequence retrieval (two-stage homology search + rescue)
CENP-A/CENH3 proteins were pulled from the 325 Helixer-predicted proteomes by a two-stage
strategy:
- **Stage 1 — DIAMOND (v2.1.11, `--sensitive`)** of all proteomes against a curated DB of
  **450 external CENP-A/CENH3 references** (published eukaryotic CENH3 DB + KEGG K11495
  CENP-A / K11253 CENH3). Hits kept if they passed a **CENP-A HMM profile** (HMMER v3.4)
  → 7,846 candidates across 253 species. **Reciprocal-best-hit** confirmed orthology.
- **Stage 2 — genome rescue (miniprot v0.15)** for species with no DIAMOND hit: references
  aligned directly to assemblies; after collapsing overlaps in 1-kb windows, 1,299/1,315
  candidates passed the HMM filter → extended detection to 21 more species.
- **H3 set:** HMMER `hmmsearch` against an H3-specific HMM, **top-3 per species** → 901.
- **69 species** have no confirmed CENP-A/CENH3 after both stages.

### 2. Curation to bona-fide orthologs (the 418 → 422 step)
The candidate pool + 5 well-characterised reference CENP-A (human P49450, mouse O35216,
zebrafish Q803H4, *Drosophila* cid Q9V6Q2, *C. elegans* hcp-3 P34470) were aligned
(MAFFT v7.526), trimmed (ClipKit v2.3.0), and a **FastTree** guide tree built. Sequences
forming a monophyletic clade with the 5 references (clearly separated from H3) were kept →
**418 CENP-A** across 256 species after removing contaminants.

The **final ML tree** (IQ-TREE2 v2.3.4, `-m VT+G4 -B 1000 -bnni -nm 2000`; VT+G4 chosen by
ModelFinder/BIC — it replaced a previously hard-coded Q.pfam+R7 that was never justified;
converged iter 404) was then re-inspected in iTOL and **6 H3-labelled sequences inside the
CENP-A clade were moved to CENP-A**, and **2 non-clustering ones moved to H3**:
- → CENP-A: `tzMemMemb1_chr_1_002492.1`, `gfChlBrun_…OY756145…_000837.1`,
  `qeDicMinu_…OX461808…_020535.1`, `ieEcdTorr_…OX439129…_002001.1`,
  `icPyrSerr1_chr_X_000403.1`, `icHalSede1_chr_9_000164.1`.
- → H3: `idSicFerr1_chr_X_002574.1`, `idSicFerr1_chr_X_002608.1`.
- Net **418 − 2 + 6 = 422 CENP-A**, **901 + 2 − 6 = 897 H3**. Do NOT recompute from scratch.

### 3. Downstream analyses (all keyed off the 422/897 split)
- **GroupSim** (specificity-determining positions, CENP-A vs H3): `scripts/run_groupsim_325sp.py`
  (+ clade-weighted `run_groupsim_cenpa_h3_clade_325sp.py`). **6 significant columns (z≥2)**
  at gap 0.80 and 0.85 — robust. Variants: Satellite-vs-Transposon CENP-A (**0 significant**),
  per-phylum, holo-vs-mono, vertebrate.
- **Split entropy:** `scripts/split_entropy_325sp.py` (positional Shannon entropy per group).
  Publication panels live in `03_entropy/` — **panel D uses the no-masking + filled version**
  (`plot_panel_D_nomasking_325sp.py`), the user's strong preference.
- **Terminal branch length** (evolutionary rate, CENP-A vs H3): the preferred rate metric
  (not root-to-tip). **249/262 species CENP-A faster**, paired Wilcoxon p≈1.2e-40.
- **CENP-A copy number / duplications** across the species tree:
  `scripts/cenpa_copy_number_analysis.py`, `plot_cenpa_copy_number.py`.
- **VGL vertebrate antibody spin-off** (kept OUT of the repo): 33 vertebrate CENP-A from 25
  species → per-class antibody peptides; top windows fall in the CATD (~MSA col 86–103).
  See memory `project_vgl_antibody_design`; outputs were in `vertebrate_antibody/`.

---

## Gotchas a fresh session MUST know
1. **Always read the curated `groups_gap*.txt` for the CENP-A/H3 labels.** Every `run_*`
   / entropy / terminal-branch script was patched to read those FIRST and only fall back to
   fasta-based inference if missing — the fasta inference gives the stale **418/901**. If
   you add/move sequences, edit `groups_gap08.txt` **and** `groups_gap085.txt` and rerun.
2. **VT+G4**, not Q.pfam+R7 (the latter was an unjustified leftover).
3. **Colours:** CENP-A `#C62828` (red), H3 `#1565C0` (blue); architecture Sat `#E53935`,
   Trans `#FB8C00`, Holo `#2d7d32`.
4. **Zero terminal branches** (606 tips) are near-identical H3 paralogs — expected, not a bug.
   The old "whiskers below 0" was a ggplot boxplot artefact; min TBL is exactly 0.
5. ***Eimeria praecox*** (`pxEimPrae`) is the one CENP-A species with **no detectable H3**
   (Apicomplexa histones too divergent for the H3 HMM even at E=10) → excluded from paired
   comparisons (so paired n is 257–262, not 325).
6. Archaeal outgroup IDs (10) are hard-coded in the scripts: OLS22332.1, OLS24873.1,
   OLS21974.1, KKK41979.1, KXH71038.1, OLS18261.1, OLS16336.1, BAD86478.1, OIO61677.1,
   OIO41945.1.

## Where things are (current repo paths)
- Alignment: `data/cenpa430_H3_archaea10.aligned.clipkit.325sp.fasta`
- Tree: `data/iqtree2_cenpa430_H3_archaea10_325sp_VTG4_bnni.{treefile,contree,iqtree}`
- Curated groups: `data/groupsim/groups_gap0{8,85}.txt`; GroupSim scores `data/groupsim/groupsim_gap0{8,85}.txt`
- Copy-number tables: `data/cenpa_copy_number_table.tsv`, `cenpa_copies_by_order.tsv`, `cenpa_gene_duplications.tsv`
- Scripts: `scripts/` (GroupSim runners + plotters, split entropy, phylogeny/phylosignal, copy number)
- Figures: `figures/`
- Detailed state: `data/HANDOFF_CENH3_H3_groupsim_entropy.md`; module map: `README.md`

## Open / likely next questions
- Re-root or re-examine any tips near the CENP-A↔H3 boundary if new genomes are added
  (re-inspect in iTOL, update the groups files, rerun).
- Map the 6 specificity positions precisely onto the histone-fold structure (CATD).
- Whether CENP-A copy-number expansions correlate with centromere architecture
  (satellite vs transposon vs holocentric) — ties into the α–ω story in `02_asr/`.
- The species-tree used elsewhere in the paper is now the **chronos correlated λ=0.1,
  62-calibration-point** tree (`01_species_tree/outputs/full_325sp_calibrated_correlatedlambda01.nwk`);
  if you map CENP-A traits onto the species tree, use that one for consistency.
