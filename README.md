# DToL phylogenomics — centromere evolution across 325 species

Reproducible analyses for the Darwin Tree of Life centromere-evolution study: a
time-calibrated 325-species phylogeny and, on it, the evolution of centromere
architecture, the centromeric histone CENP-A/CENH3, and centromeric satellite DNA
(including the CENP-B box). Each module is self-contained with its own `README.md`,
`scripts/`, `data/` and `figures/`.

## Modules
| | Analysis | What it does | Key output |
|---|---|---|---|
| **01** | [Species tree](01_species_tree/) | 325-species tree via FastSpeciesTree (DIAMOND pseudo-alignment → partitioned IQ-TREE) + `chronos` time-calibration (62 TimeTree points, correlated λ=0.1) | calibrated chronogram |
| **02** | [ASR](02_asr/) | Mk models of centromere-architecture evolution; alpha–omega cyclical hypothesis | ASR trees, model table |
| **03** | [Entropy](03_entropy/) | per-position conservation of the CENP-A/CENH3 alignment | positional-entropy profile |
| **04** | [CENP-A phylogeny](04_cenpa_phylogeny/) | CENP-A/CENH3 vs H3 phylogeny (IQ-TREE VT+G4) + GroupSim | 422 CENP-A / 897 H3 tree |
| **05** | [Satellite similarity](05_satellite_similarity/) | satellite divergence / similarity decay & half-life | decay curves |
| **06** | [CENP-B box](06_cenpb_box/) | CENP-B box screen of the satellites (2 methods, HG002 controls) | per-species motif scatter, HTML report |

Each module's README lists its run order and inputs. Analyses target the **325
published species** (the calibrated species-tree tips).

## Data
- Repo: `github.com/jacgonisa/DToL_phylogenomics`.
- Large inputs (assemblies, `all.satellites.txt`, HG002) live outside the repo; each
  module README states the paths and how its `data/` tables were produced.

The species tree is built with **FastSpeciesTree** (a DIAMOND-anchored pseudo-alignment
of BUSCO orthologs → partitioned IQ-TREE); see `01_species_tree/METHODS.md`. This
replaces an earlier BUSCO-supermatrix workflow, which has been removed.
