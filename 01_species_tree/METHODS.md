# 01 — Species tree & time calibration (methods)

Verified against the FastSpeciesTree run
`…/2026_trees/fast_species_tree_iqtree_only325/Results/` and the calibration tables
in `data/` (`calibration_table_325sp_supp_64.tsv`, `over_calib.tsv`). Superscript
numbers are manuscript reference markers.

## Species-tree inference
A maximum-likelihood species tree was inferred from the 325 DToL proteomes using
FastSpeciesTree⁶⁰ (v1.0) in a single "sensitive" run, which performs ortholog
identification, alignment and tree inference end-to-end. Each proteome was searched with
DIAMOND¹¹¹ blastp against the eukaryota_odb10 BUSCO reference orthologs; for each
single-copy ortholog the best-hit local alignment was anchored to the reference-query
coordinates and gap-padded to the query length, and these per-ortholog alignments were
concatenated into a query-anchored pseudo-alignment (orthologs absent from a proteome
encoded as gaps). Orthologs present in at least 80% of species (no more than 20% missing
taxa) were retained, yielding **194 single-copy loci**, and species with excessive gaps
were excluded. After FastSpeciesTree's internal column trimming, the supermatrix
comprised **325 taxa × 52,754 amino-acid sites (4.22% missing data)**. Within the same
run, a partitioned maximum-likelihood tree was inferred with IQ-TREE¹¹³ (**v2.3.4**)
under an **edge-linked proportional partition model** (`-p`), with ModelFinder selecting
the best-fit substitution model independently per locus from LG, JTT, Q.INSECT, Q.YEAST,
Q.BIRD, Q.MAMMAL and Q.PLANT, and among +I, +G and +I+G rate models. Across the 194
partitions, Q.INSECT (171 loci) and LG (15) were selected most often, with Q.YEAST (6),
Q.PLANT (1) and JTT (1) elsewhere (Q.BIRD and Q.MAMMAL were offered but never selected).
Branch support was assessed with 1,000 ultrafast bootstrap replicates (UFBoot) and 1,000
SH-aLRT replicates.

FastSpeciesTree was invoked as (sensitive mode → IQ-TREE):
```bash
python FastSpeciesTree.py -f all_proteomes_only325/ -o fast_species_tree_iqtree_only325 \
  -s sensitive -t 32
```
which internally ran the partitioned IQ-TREE command:
```bash
iqtree -T 32 -s trim_psuedo_alignment.fasta -p trim_IQTree_Partition_file.partitions \
  -B 1000 --alrt 1000 -st AA \
  -mset LG,JTT,Q.BIRD,Q.MAMMAL,Q.INSECT,Q.PLANT,Q.YEAST -mrate I,G,I+G -m MFP
```

> **Note (no MAFFT/trimAl).** FastSpeciesTree builds a DIAMOND-anchored pseudo-alignment
> and trims columns itself; MAFFT and trimAl are *not* part of this pipeline (its env is
> diamond + iqtree + veryfasttree only). The env pins iqtree 3.0.1, but the actual run
> used the `iqtree` on PATH — **v2.3.4** per the run log.

## Time calibration
The inferred tree topology was time-calibrated using the `chronos` function of the
**ape** R package⁶². Calibration nodes were placed at well-supported clades whose
divergence could be defined by the most recent common ancestor of two sampled taxa, and
each node was assigned a minimum–maximum age constraint. **Sixty-four constraints** were
used in total: **one root constraint** (Eukaryota, 1,085–1,671 Mya) and **63 internal
nodes** distributed across **Opisthokonta (n=1), Metazoa (n=32), Viridiplantae (n=22)
and Fungi (n=8)**. For **62** of these, the bounds were taken from the corresponding
TimeTree 5⁶³ divergence-time confidence interval, retrieved automatically with PAReTT
(https://github.com/LSLeClercq/PAReTT). The remaining **two** nodes had no TimeTree
record and were constrained from the literature: the caddisfly family **Limnephilidae**
(10–18 Mya) and the rush genus **Luzula** (crown ~6–10 Mya). Constraint depths spanned
the full tree: **39 of the 64 constraints were older than 100 Mya** (23 older than 300
Mya), while at the shallow end **six were congeneric (genus-level) splits** and **8 were
younger than 20 Mya** — the shallowest being the *Falco* (~2 Mya) and *Thunnus*
(~3–4 Mya) divergences.

To choose a dating model, four approaches were compared and benchmarked against TimeTree
by node-age concordance (`parrett_*` tables; Supplementary Fig. X): a **correlated-rates**
penalised-likelihood model (`chronos`, `model="correlated"`, λ=1), a **relaxed-rates**
model (`model="relaxed"`, λ=10), a **rate-smoothed** tree constrained only at the root,
and an **uncalibrated** ultrametric tree (treePL and RelTime were additionally examined).
The correlated-rates model gave the closest agreement with TimeTree node ages and was
retained; the final chronogram was estimated under this model with all 64 constraints
applied simultaneously to the full tree.

```r
chronos(tree_full, lambda = 1, model = "correlated", calibration = calib_df)
```
