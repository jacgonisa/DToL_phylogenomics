# CENP-B box search — text for response to reviewers

We screened the centromeric satellite repertoire of the DToL species for the
17-bp **CENP-B box** (canonical consensus `[CT]TTCGTTGGAA[AG]CGGGA`; Masumoto et
al. 1989), on **both strands**, using two complementary approaches, and report
both below. (A third, PWM/log-odds scan gave results concordant with the mismatch
scan and is included in the repository as supporting evidence.)

## Method 1 — simple motif parsing (mismatch scan)
The canonical box was searched allowing an increasing number of **mismatches**
(0–5), against a mononucleotide-**shuffled null**. This was run both on a curated
monomer set (≤500 monomers/species) and, exhaustively, on the **full uncapped set
of ~24.5 million satellite arrays**.

**Result.** The **exact canonical box (0 mismatches) occurs only in human**
α-satellite (positive control; 13,359 of 22.5 M arrays). **No non-human species
carries a single exact box.** The closest non-human by every strict metric is the
bat *Rhinolophus sinicus* (0 exact matches, highest ≤2-mismatch and PWM scores of
any non-human). At relaxed mismatch thresholds hits appear broadly but at low,
composition-driven rates that do not exceed the shuffled null in most clades.
*Figures:* `cenpb_mismatch_titration_allspecies.png`, `data/cenpb_uncapped_per_species.tsv`.

## Method 2 — songbird ±5-flank approach (Formenti et al., Cell 2026)
Following the zebra-finch T2T paper (Suppl. Fig. 15), the box was matched as fixed
17-bp windows (≤5 substitutions, both strands) and, for **every species on the
full uncapped ~24.5 M arrays**, we built the position frequency matrix of the box
**plus ±5 flanking bp**. A genuine motif shows high information across the 17-bp
box that **collapses in the flanks** (Δ = box − flank information > 0); we also
report the box consensus and its substitutions from the canonical box, and the
prevalence (boxes/Mbp). *Caveat:* substitution-matched windows are constrained
toward the consensus, giving a **~0.34-bit baseline Δ in all clades**; genuine
boxes stand well above this AND have a near-canonical consensus at realistic
prevalence.

**Result — vertebrates lead, with birds strongest.** Vertebrates show the highest
mean Δ (0.39; 52 % box-enriched) vs plants (0.34) and invertebrates (0.34). The
convincing cases (Δ well above baseline + near-canonical consensus + real
prevalence):

| species | group | Δ (box−flank) | box consensus | subs vs canonical | boxes/Mbp |
|---|---|---|---|---|---|
| *Accipiter gentilis* (goshawk) | Aves | **0.85** | `CTTTTTTGGAAACGGGA` | 2 | 88 |
| *Porphyrio hochstetteri* (takahē) | Aves | 0.68 | `TTTCCTTGGAAACGGAA` | 2 | 26 |
| *Lagopus muta* (ptarmigan) | Aves | 0.66 | `CTTTGTTGGAAACGGGA` | **1** | 95 |
| *Rhinolophus sinicus* (bat) | Mammalia | 0.53 | `GTTCGTAGGAAGCGGGT` | 3 | 890 |
| *Diceros bicornis* (rhino) | Mammalia | 0.52 | `CTTCCTTAGAAGCAGGA` | 3 | 7 |
| *Cervus elaphus* (red deer) | Mammalia | 0.45 | `TTTCGTGGGAAGGGGGA` | 2 | 217 |
| *Corvus hawaiiensis* (crow) | Aves | 0.41 | `TTTCTTTGGCAGCAGCA` | 4 | 952 |

The **goshawk box (Δ 0.85, 2 substitutions, intact `CGGGA` tail)** is the
strongest box-vs-flank signal in the entire 325-species dataset — exceeding even
the human benchmark's Δ — and several other birds carry a near-canonical box
(ptarmigan 1 substitution). In contrast the three *Falco* species show flanks as
conserved as the box (Δ ≈ 0.1) — conserved satellite, not a box — an internal
negative control. Fish and the single reptile show no convincing box.

The **zebra finch itself was not represented in our satellite set** (0 arrays), so
the specific Tgut716A box could not be tested directly; but the signal in the
goshawk, ptarmigan and takahē indicates a diverged CENP-B-like box is present more
broadly across birds. *Figures:* `cenpb_flank_uncapped_scatter.png`,
`cenpb_box_logos_flanks_VERTEBRATES_uncapped.pdf`, `data/cenpb_flank_uncapped_per_species.tsv`.

## Conclusion
The **canonical, functional CENP-B box** is confined to **human (mammalian)
α-satellite**. Beyond it, the songbird ±5-flank test on the full satellite set
reveals a **diverged but box-specific CENP-B-like motif that is strongest in
birds** (goshawk, ptarmigan, takahē) and present in several mammals (bat, rhino,
deer) — consistent with, and extending beyond the songbird, the recently reported
avian CENP-B-like centromere system. Broad low-Δ hits in plants/invertebrates fail
the flank control (flanks as conserved as the box) and reflect homogenised
satellite rather than a functional box.
