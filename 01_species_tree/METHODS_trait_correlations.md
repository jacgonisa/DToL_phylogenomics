# Methods — trait correlations & phylogenetic comparative analysis

Paragraph for the manuscript Methods (matches the CENP-A/CENH3 retrieval style).

---

**Trait correlations and phylogenetic comparative analysis.** For each species, seven
centromere/genome traits were assembled — HOR regimentation, HOR score, satellite monomer
length, satellite GC content, host genome GC content, genome size and chromosome number —
from the per-family TRASH annotations and the HOR metrics (regimentation and HOR score).
For poly- and di-typic species carrying multiple centromeric repeat families, each trait
was taken from the dominant array (the family with the highest monomer copy number);
genome-level traits (genome size, GC content, chromosome number) are family-invariant.
Monomer length and genome size were log₁₀-transformed. Pairwise associations were computed
both naively (Pearson correlation) and with correction for phylogenetic non-independence
on the chronos-calibrated 325-species time tree, using phylogenetic generalised least
squares (PGLS) with Pagel's λ estimated by maximum likelihood (caper119 v1.0.4 `pgls`,
`lambda = "ML"`; Brownian-motion fit used as a fallback where the λ optimisation failed to
converge). The phylogenetically corrected correlation coefficient for each pair was
obtained from the PGLS slope t-statistic (r = t/√(t² + df); effectsize122 v1.0.1
`t_to_r`), which is λ-adjusted — collapsing to the ordinary correlation when λ = 0 (no
phylogenetic signal) and to the Brownian-motion correlation when λ = 1. As a cross-check,
the Brownian-motion correlation was also computed from the evolutionary covariance matrix
(phytools121 v2.5.2 `phyl.vcv`; identical to a Felsenstein independent-contrasts
correlation, ape120 v5.8.1 `pic`) and was concordant. Phylogenetic signal for each trait
was quantified with Blomberg's *K* and Pagel's λ (phytools `phylosig`, 999 permutations).
To distinguish direct from mediated associations among HOR regimentation, HOR score and
monomer length, partial correlations were computed from multiple-predictor PGLS (each
predictor's partial slope tested, and converted to a partial correlation via `t_to_r`,
with the third trait held constant). All analyses were performed in R v4.3.3 and were
repeated under an alternative aggregation (copy-number-weighted genomic mean across
families), yielding equivalent results. Visualisations used seaborn122 v0.13.2.

---

Refs to add: caper (Orme et al.), ape (Paradis & Schliep), phytools (Revell), seaborn (Waskom).
Scripts: `trait_pgls_correlations_325sp.R` (stats), `trait_pairplot_seaborn_325sp.py` (figures).
