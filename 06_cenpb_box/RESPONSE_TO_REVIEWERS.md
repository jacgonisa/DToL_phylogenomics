# CENP-B box search — text for response to reviewers

We screened the centromeric satellite repertoire of the DToL species for the
17-bp **CENP-B box** (canonical consensus `[CT]TTCGTTGGAA[AG]CGGGA`; Masumoto et
al. 1989), on **both strands**, using two complementary methods. Human HG002
α-satellite was scanned throughout as the positive benchmark.

## Which satellites were searched
- **Input:** `all.satellites.txt` — every satellite-repeat *array* annotated across
  the DToL genome assemblies (one sequence per array, tagged by species).
- **Species set:** we restricted to the **325 published species** (allowlist
  `species_325.txt`, derived directly from the tips of the calibrated species
  tree). `all.satellites.txt` contained 12 extra species that are not in the
  final tree; these were excluded.
- **"Uncapped":** for each species we searched **every satellite array** (both
  strands) — **not** a ≤500-monomer subsample. This covers **20.3 million arrays
  (3.58 Gbp) across the 162 of 325 species that carry satellite annotations**.

## Human HG002 benchmarks (positive + negative control)
Both methods were validated on human HG002 satellites:
- **Positive — α-satellite:** 50,000 α-satellite monomers (171 bp each; 8.55 Mbp),
  which carry the functional CENP-B box.
- **Negative — HSat1/2/3:** 7,982 HSat arrays (161.8 Mbp) extracted from the
  hg002v1.1 assembly via the v1.1 CenSat annotation — human satellites that do
  **not** carry CENP-B boxes.

| control | Method 1: canonical hits (per Mbp) | broad enrichment | Method 2: flank Δ (windows/Mbp) |
|---|---|---|---|
| **α-sat (positive)** | **13,359 (1,562 /Mbp)** | 676× | **0.56 (3,247 /Mbp)** |
| **HSat (negative)** | 5 (0.03 /Mbp) | 0.16× (below random) | 0.23 (93 /Mbp) |

Both methods cleanly separate the controls: canonical boxes are **~50,000× denser
per Mbp in α-satellite than in HSat**, HSat's looser-motif hits fall *below* random
expectation, and the flank Δ and box prevalence are far higher in α-satellite. This
confirms the pipeline detects the functional box where it exists and stays quiet on
human satellite that lacks it.

## Method 1 — exact IUPAC motif parsing (Fachinetti)
Following Barra & Fachinetti et al. (bioRxiv 2026.05.25.727640), we counted three
**exact IUPAC** CENP-B box definitions and compared each species' observed count to
a random expectation from its own base composition (obs/expected enrichment):

| tier | motif | note |
|---|---|---|
| canonical | `YTTCGTTGGAARCGGGA` | the functional box (degenerate only at pos 1 & 12) |
| broad | `NTTCGNNNNANNCGGGN` | keeps the essential `TTCG…CGGG` core |
| degenerated | `YTTCGNNNNANRCGGGN` | looser interior |

**Result.** In the human benchmark the method works overwhelmingly (canonical
enrichment ~5.7 × 10⁶; broad 676×; degenerate 2655× over random). Across the 325
DToL species, **no species carries a single exact *canonical* box** (canonical = 0
in every clade). The looser tiers yield only **trace** exact hits, at or below
random expectation in vertebrates as a group (mean broad enrichment 0.45×): the
largest per-species counts are *Gallus gallus* (12 degenerate; 19× over random),
*Canis lupus* (16), *Trachurus trachurus* (12 broad) — a handful of arrays, not an
enriched motif. **By the strict motif definition, essentially nothing is found
outside the human benchmark.** *Figures:* `cenpb_paper_motifs.png`;
*data:* `cenpb_paper_motifs_per_species.tsv`, `cenpb_paper_motifs_per_clade.tsv`.

## Method 2 — songbird ±5-flank test (Formenti et al., Cell 2026)
Following the zebra-finch T2T paper (Suppl. Fig. 15), the box was matched as fixed
17-bp windows (≤5 substitutions, both strands) and, for **every species on the full
uncapped set**, we built the position frequency matrix of the box **plus ±5
flanking bp**. A genuine motif shows high information across the 17-bp box that
**collapses in the flanks** (Δ = box − flank > 0); we also report the box consensus,
its substitutions from canonical, and prevalence (boxes/Mbp). *Caveat:*
substitution-matched windows are constrained toward the consensus, giving a modest
baseline Δ in every clade, so the signals of interest are those well above baseline
**and** with a near-canonical consensus at realistic prevalence.

**Result — the signal is in birds.** The strongest box-vs-flank enrichment in the
whole dataset is in the goshawk, with several other birds close behind:

| species | group | Δ (box−flank) | box consensus | subs vs canonical | boxes/Mbp |
|---|---|---|---|---|---|
| *Accipiter gentilis* (goshawk) | Aves | **0.85** | `CTTTTTTGGAAACGGGA` | 2 | 88 |
| *Porphyrio hochstetteri* (takahē) | Aves | 0.68 | `TTTCCTTGGAAACGGAA` | 2 | 26 |
| *Lagopus muta* (ptarmigan) | Aves | 0.66 | `CTTTGTTGGAAACGGGA` | **1** | 95 |
| *Diceros bicornis* (rhino) | Mammalia | 0.52 | `CTTCCTTAGAAGCAGGA` | 3 | 7 |
| *Cervus elaphus* (red deer) | Mammalia | 0.45 | `TTTCGTGGGAAGGGGGA` | 2 | 217 |

The goshawk box (Δ 0.85, 2 substitutions, intact `CGGGA` tail) is box-specific with
random flanks; several birds carry a near-canonical box (ptarmigan 1 substitution).
The three *Falco* species instead have flanks as conserved as the box (Δ ≈ 0.1) —
conserved satellite, not a box (internal negative control). The zebra finch itself
was not represented in our satellite set. *Figures:*
`cenpb_flank_uncapped_scatter.png`, `cenpb_box_logos_flanks_VERTEBRATES_uncapped.pdf`;
*data:* `cenpb_flank_uncapped_per_species.tsv`.

## Reconciling the two methods
The methods are complementary and disagree informatively on birds:
- **Method 1 (exact IUPAC)** is conservative and *position-aware* — it requires the
  essential `TTCG…CGGG` core. No DToL species passes it.
- **Method 2 (±5-flank)** is permissive. Birds pass it, but their divergence falls
  **on the essential 5′ `CG`** (goshawk `TTCG`→`TTTT` at positions 4–5), which is
  exactly why they score 0 under the exact broad motif.

## Conclusion
The **canonical, functional CENP-B box** is confined to **human (mammalian)
α-satellite**; by the strict exact-motif definition it is absent from all 325 DToL
species. The **±5-flank test does reveal a diverged, box-shaped motif that is
strongest in birds** (goshawk, ptarmigan, takahē) — but this motif has lost the
CENP-B-essential 5′ CG, so it is a plausible **TIGD4-type avian box** rather than a
mammalian CENP-B box, consistent with the songbird paper's own model. Whether it is
functional cannot be settled from the mammalian motif alone.
