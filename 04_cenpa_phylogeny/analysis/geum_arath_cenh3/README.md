# Geum vs Arabidopsis CENH3 — small side analysis

Compares the CENH3/CENP-A copies of the two *Geum* species in the 325-genome set
against *Arabidopsis thaliana* CENH3 (HTR12).

## Sequences (`*.input.fasta`)
- `Riva_680`, `Riva_678` — *G. rivale* (`drGeuRiva`), the 2 curated CENP-A copies
- `Urba_S4`, `Urba_S1`, `Urba_S3` — *G. urbanum* (`drGeuUrba1`), the 3 curated CENP-A copies
- `AraTh_CENH3` — dataset's *A. thaliana* CENH3 (`ddAraThal4_chr_1_000021.1`, curated CENP-A)
- `AraTh_HTR12_ref` — canonical UniProt HTR12 (Q8RVQ9), sanity-check reference
- `Maize_CENH3` — *Zea mays* CENH3 (monocot), used as outgroup to root the tree
  (`geum_arath_maize_*` files; user-supplied sequence)

IDs taken from the curated CENP-A group (`data/groupsim/groups_gap085.txt`); Geum/AraTh
seqs degapped from the full 325sp MSA, HTR12 fetched from UniProt.

## Files
- `*.mafft_linsi.aln` — untrimmed alignment (`mafft --localpair --maxiterate 1000`)
- `*.VTG4.treefile` / `.contree` / `.iqtree` — ML tree, `iqtree2 -m VT+G4 -B 1000`
  (VT+G4 to match the main CENP-A tree)
- `geum_arath_maize_*.maize_rooted.*` — same, with maize as an explicit `-o` outgroup
  (rooted version; eudicot Geum + Arabidopsis on one side, monocot maize as outgroup)

## Result
- `AraTh_CENH3` ≡ `AraTh_HTR12_ref` (~1 residue) → dataset's Arabidopsis call is the real CENH3.
- `Riva_680` = `Urba_S4` — 100% AA identical: the bona-fide CENH3 is conserved between the two Geum.
- `Urba_S1` ~94% id (sister, 94% UFBoot); `Riva_678` + `Urba_S3` are divergent paralogs.
- Divergence sits in the N-terminal tail; the histone fold is conserved across all, incl. Arabidopsis.
