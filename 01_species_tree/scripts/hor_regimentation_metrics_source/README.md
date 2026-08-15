# HOR score & regimentation — source scripts (Piotr Włodzimierz)

Provenance for the two per-family centromere metrics used in the ringed tree
(`plot_tree_regimentation_chronos.R`) and the continuous-trait ASR
(`asr_continuous_regimentation_325sp.R`), i.e. the `HOR_score` and
`regimentation_score` columns of `cen_families.csv`.

These are the original scripts as sent by Piotr; the hard-coded `/home/pwlodzimierz/...`
paths are his cluster and are kept verbatim for provenance (not runnable here).

Both start from **TRASH self-HOR output**: a table of higher-order-repeat (HOR) blocks
per chromosome/satellite, each row a pair of matching blocks `A` and `B` (each block is a
run of consecutive repeat monomers), with `start_A/end_A/start_B/end_B` (monomer indices)
and `block.size.in.units`.

## `HOR_score.R` — how *much* HOR structure a monomer participates in
Per **repeat monomer**, the fraction of all other monomers on the chromosome that it
shares any HOR pattern with.

1. Expand every HOR block-pair into monomer-level `element → partner` pairs (both
   A→B and B→A), dedupe, drop self-pairs.
2. For each monomer (`element`), count its **unique HOR partners** (`num_interactors`).
3. `HOR_score = 100 * num_interactors / (total monomers on chr)`  → 0–100.

So HOR score = **connectedness / abundance of HOR relationships**. High = the monomer is
part of an HOR motif repeated widely across the array; 0 = not in any HOR. It is optimised
/ terse but the logic is just "per repeat, fraction of other repeats it HOR-pairs with."

## `67_regimentation.R` — how *regularly spaced* the HORs are along the array
A spectral-periodicity pipeline asking whether HOR blocks recur at a **regular period**
(regimented, e.g. human α-satellite) vs irregularly (heterogeneous).

1. `dist = start_B - start_A` — spacing (in monomers) between the two blocks of each HOR.
2. `build_count_signal` — histogram of distances, **edge-corrected** by `/(N-d)` possible
   pairs so array-length combinatorics don't masquerade as structure.
3. `run_periodogram` — `spec.pgram` power spectrum (frequency = cycles/repeat).
4. `peak_prominence` + **Harmonic Product Spectrum** — pick the fundamental period,
   resisting octave (harmonic) errors; `estimated_period = 1/peak_freq`.
5. `spectral_entropy` (0 = periodic, 1 = spread) and **Fisher's g-test** p-value
   → `classify_periodicity`: `periodic` vs `heterogeneous` at `alpha`.
6. `window_track` slides along the array classifying each window; `summarize_regions`
   merges contiguous same-period runs (smoothing short "flickers"); `summarise_chromosome`
   tallies the **fraction of the array in periodic regions**.

So **regimentation = the % of the centromeric array that sits in regularly-periodic HOR
regions** (`100 − fraction with period 0`; see line ~519). High = ordered/regimented HORs;
low = disordered.

## The two are complementary
- **HOR score** = amount/degree of HOR sharing (can be high even if irregular).
- **Regimentation** = regularity/order of the HOR spacing (periodicity).

Aggregation to the tree: per species we take the **dominant array** (family with the most
copies) — see `plot_tree_regimentation_chronos.R`. Holocentrics are not scored (grey/NA).
