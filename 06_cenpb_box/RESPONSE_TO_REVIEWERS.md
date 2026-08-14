# CENP-B box search — response to reviewers

**How we searched.** We screened the centromeric satellites of the 325 DToL species
for the 17-bp **CENP-B box** (`[CT]TTCGTTGGAA[AG]CGGGA`; Masumoto et al. 1989), on
both strands, across all annotated **putative centromeric satellite** arrays (162 of 325 species
have satellites). We used two approaches, with human HG002 α-satellite as a positive
control and human HSat as a negative control:

1. **Exact motif match** — the canonical box and two looser IUPAC variants (Barra &
   Fachinetti 2026), counted against a composition-matched (dinucleotide-shuffle) null.
2. **Relaxed / de-novo match** — box-like windows within ≤5 substitutions, plus an
   unbiased scan of every window against a shuffled-sequence null.

**What we found.**
- The **exact CENP-B box occurs only in human α-satellite** (positive control). It is
  **absent from all 325 DToL species** — no species carries a canonical box.
- Relaxed (≤5-substitution) "box-like" hits appear in ~94% of species, but this is
  expected by chance (≤5 subs on 17 bp is very permissive) and does **not** exceed a
  shuffled-sequence null in any species — whereas human α-satellite clearly does.
- The strongest relaxed signal is in **birds** (e.g. the goshawk has a satellite motif
  15/17 identical to the box, lacking only the central CpG). We report this as a
  **suggestive candidate** CENP-B-like motif, not a confirmed box; protein binding was
  not tested.

**In short:** a functional CENP-B box is confined to human/mammalian α-satellite; we
find no confirmed CENP-B box elsewhere in the DToL set, with birds as a candidate worth
follow-up.

**Suggested figures.**
- Method 1 (exact match): **`cenpb_paper_motifs.png`** — enrichment vs motif density; only
  α-satellite is a real box (top-right), all DToL clades sit near the null.
- Method 2 (de-novo): **`cenpb_denovo_bestwindow.png`** — best box-like window vs a shuffled
  null; α-satellite is above the diagonal, all DToL (birds included) sit on it.
- Overview: **`cenpb_box_tree_325sp.png`** — the signal mapped onto the chronogram.

*(Full methods, all figures and per-species data: `cenpb_box_report.html`.)*
