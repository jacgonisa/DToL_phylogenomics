# CENP-B box search — text for response to reviewers

We screened the centromeric satellite repertoire of the DToL species for the
17-bp **CENP-B box** (canonical consensus `[CT]TTCGTTGGAA[AG]CGGGA`; Masumoto et
al. 1989), on **both strands**, using two complementary approaches, and report
both below. (A third, PWM/log-odds scan gave results concordant with the mismatch
scan and is included in the repository as supporting evidence.)

## Method 1 — exact IUPAC motif tiers (as in the functional-box paper)
Following Barra/Fachinetti et al. (bioRxiv 2026.05.25.727640), we searched three
**exact IUPAC** CENP-B box definitions (both strands) across the **full uncapped
~24.5 M satellite arrays**, and compared each species' observed counts to a random
expectation from its own base composition (obs/expected enrichment):

| tier | motif | note |
|---|---|---|
| canonical | `YTTCGTTGGAARCGGGA` | the functional box (degenerate only at pos 1 & 12) |
| broad | `NTTCGNNNNANNCGGGN` | keeps the essential `TTCG…CGGG` core |
| degenerated | `YTTCGNNNNANRCGGGN` | looser interior |

**Result.** **No non-human species carries a single exact *canonical* box**
(canonical = 0 in every clade; human α-satellite is the positive control, from the
exact scan = human-only). For the looser tiers, one vertebrate stands out sharply:
the bat **_Rhinolophus sinicus_** — **22,677 exact _broad_-motif hits (37× over
random)** plus 289 degenerate — an order of magnitude above any other non-human.
All other vertebrates have only trace exact hits (*Gallus* 12, *Trachurus* 12,
*Lagopus* 10). Plants show degenerate-motif enrichment but zero canonical,
matching the paper's observation that degenerate motifs occur stochastically and
reflect neutral variation. *Figures:* `cenpb_paper_motifs.png`,
`data/cenpb_paper_motifs_per_species.tsv`, `data/cenpb_paper_motifs_per_clade.tsv`.
(An earlier mismatch-titration, `cenpb_mismatch_titration_allspecies.png`, gives a
concordant but less clearly-tiered view and is retained as supporting.)

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

## Reconciling the two methods
The two methods agree on the bat and disagree on birds, informatively:
- The **exact-IUPAC method (Method 1)** is conservative and *position-aware*: it
  requires the essential `TTCG…CGGG` core. Only the bat *R. sinicus* passes.
- The **±5-flank method (Method 2)** is permissive: it tolerates substitutions
  anywhere. Birds pass here, but their divergence falls **on the essential 5′ `CG`**
  (goshawk consensus `CTTTTTTGGAAACGGGA` — the `CG` at positions 4–5 has become
  `TT`), which is exactly why they score 0 under the exact broad motif.

So the **bat _Rhinolophus sinicus_ is the only vertebrate robust to both methods**.
The avian signal is real as a *diverged, box-shaped* motif but does not preserve
the mammalian CENP-B essential positions — consistent with the songbird paper's
own premise that birds may use a **TIGD4** homolog whose sequence preference need
not match mammalian CENP-B. Whether the avian box is functional therefore cannot be
settled from the mammalian motif alone.

## Conclusion
The **canonical, functional CENP-B box** is confined to **human (mammalian)
α-satellite**. Among all other species, the strongest and most robust signal — by
both the exact-IUPAC and the flank method — is the bat **_Rhinolophus sinicus_**
(37× broad-motif enrichment). Birds (goshawk, ptarmigan, takahē) carry a
**diverged, box-shaped motif** that passes the permissive flank test but not the
position-aware exact motif, because it has lost the CENP-B-essential 5′ CG — a
plausible **TIGD4-type avian box** rather than a mammalian CENP-B box. Broad
low-signal hits in plants/invertebrates fail both controls and reflect homogenised
satellite rather than a functional box.
