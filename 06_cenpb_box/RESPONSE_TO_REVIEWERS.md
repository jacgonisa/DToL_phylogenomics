# CENP-B box search — text for response to reviewers

We screened the centromeric satellite repertoire of the DToL species for the
17-bp **CENP-B box** (canonical consensus `[CT]TTCGTTGGAA[AG]CGGGA`; Masumoto et
al. 1989), on **both strands**, using two complementary methods. Human HG002
α-satellite was scanned throughout as the positive benchmark.

**Terminology.** A CENP-B *box* is a **functional** motif that binds the CENP-B
protein. We detect *sequence* matches only, so below we call them **motifs** or
**candidate boxes**; none is a confirmed functional box (protein binding was not
tested). A "candidate box" = a box-specific match (flank information > flank
control) that is near-canonical (≤2 substitutions).

A self-contained HTML version with all figures is in `cenpb_box_report.html`.

## Which satellites were searched
- **Input:** `all.satellites.txt` — every **candidate** satellite-repeat *array*
  called across the DToL assemblies (candidates, not a curated satellite set; one
  sequence per array, tagged by species).
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

| control | canonical hits (/Mbp) | degenerated enrichment (mono / **dinuc** / AE-shuffle) | flank Δ (windows/Mbp) |
|---|---|---|---|
| **α-sat (positive)** | **13,359 (1,562 /Mbp)** | 676× / **4,283×** / 4,524× | **0.56 (3,247 /Mbp)** |
| **HSat (negative)** | 5 (0.03 /Mbp) | 0.16× / **0.06×** / 0.44× (below random) | 0.23 (93 /Mbp) |

Both methods cleanly separate the controls: canonical boxes are **~50,000× denser
per Mbp in α-satellite than in HSat**, HSat's looser-motif hits fall *below* random,
and the flank Δ and box prevalence are far higher in α-satellite.

**Null model.** Enrichment is reported against a **dinucleotide-preserving null**
(first-order Markov expectation, the analytic form of an **Altschul–Erikson doublet
shuffle**), which controls for base *context* (e.g. the CpG in `TTCG`/`CGGG`) rather
than base composition alone. On the α-sat benchmark (degenerated tier) the analytic dinucleotide null
(4,283×) matches an explicit Altschul–Erikson shuffle (4,524×), validating it. A
0-order (mononucleotide) null and per-Mbp counts are also reported
(`cenpb_human_benchmark.tsv`, `ae_shuffle.py`).

## Method 1 — exact IUPAC motif parsing (Fachinetti)
Following Barra & Fachinetti et al. (bioRxiv 2026.05.25.727640), we counted three
**exact IUPAC** CENP-B box definitions and compared each species' observed count to
a random expectation from its own base composition (obs/expected enrichment):

| tier | motif | note |
|---|---|---|
| canonical | `YTTCGTTGGAARCGGGA` | the functional box (degenerate only at pos 1 & 12) |
| broad | `YTTCGNNNNANRCGGGN` | Sugimoto 1998; intermediate |
| degenerated | `NTTCGNNNNANNCGGGN` | most permissive (N at pos 1 & 12) |

**Result.** In the human benchmark the method works overwhelmingly (canonical
dinucleotide-null enrichment ~1.3 × 10⁷; broad 14,457×; degenerated 4,283×). Across the 325 DToL
species, **no species carries a single exact *canonical* box** (canonical = 0 in
every clade). The looser tiers yield only **trace** exact hits — a handful of
arrays (*Gallus gallus* 12, *Canis lupus* 17, *Trachurus trachurus* 12 — all in the
degenerate tier), not an enriched motif. Importantly, the apparent broad-motif
enrichment in plants **deflates under the dinucleotide null** (Viridiplantae
broad 11.7× → **3.0×**), i.e. most of it was base-context structure, not the
box. **By the strict motif definition, essentially nothing is found outside the
human benchmark.** *Figures:* `cenpb_paper_motifs.png`;
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

**Is the box-like identity real?** We score identity **per window** (not on the
consensus — the consensus is inflated for high-copy satellites, where random per-window
deviations average back to canonical, ≈100%). Against a per-window composition shuffle,
observed mean per-window identity is **~71%** vs **≈29% chance** — the windows are
compositionally box-like. Per-window identity sits near the ≤5-substitution floor for
*all* species, so it does not single out any clade; the species-specific box-likeness
(the goshawk's 15/17 consensus) is a **consensus/coherence** property — its windows deviate
*consistently* at the eroded CpG — best seen in the logos and alignment
(`cenpb_identity_shuffle_null.py`).

## Reconciling the two methods (the goshawk case)
The methods disagree on birds in a diagnostic way. The **goshawk consensus is 15/17
identical to canonical**; its only two substitutions land **exactly on the 5′ CpG**
(positions 4–5, `C,G`→`T,T`) that *all three* IUPAC tiers hold fixed (the `TTCG`
anchor):

```
 canonical  Y T T C G T T G G A A R C G G G A
 goshawk    C T T T T T T G G A A A C G G G A
                  ^ ^  (eroded CpG)
```

- **Method 1 (exact IUPAC)** is *position-aware*: the CpG anchor is non-negotiable →
  **0 matches** at every tier (canonical = broad = degenerate = 0).
- **Method 2 (≤5-substitution)** is *distance-based*: no fixed position → **1,722
  windows**, 15/17 identity, Δ 0.85.

The eroded position is a **CpG** — both a CENP-B protein-contact base and a
methylation/deamination hotspot — so this is a **candidate divergent box** (plausibly
read by a relaxed-specificity TIGD-family transposase, per the songbird model), not a
canonical CENP-B box. See the alignment figure `cenpb_goshawk_alignment.png`.

## Conclusion
The **canonical, functional CENP-B box** is confined to **human (mammalian)
α-satellite**; by the strict exact-motif definition it is absent from all 325 DToL
species. The **±5-flank test does reveal a diverged, box-shaped motif that is
strongest in birds** (goshawk, ptarmigan, takahē) — but this motif has lost the
CENP-B-essential 5′ CG, so it is a plausible **TIGD4-type avian box** rather than a
mammalian CENP-B box, consistent with the songbird paper's own model. Whether it is
functional cannot be settled from the mammalian motif alone.
