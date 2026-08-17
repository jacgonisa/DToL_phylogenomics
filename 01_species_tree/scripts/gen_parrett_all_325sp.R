#!/usr/bin/env Rscript
# gen_parrett_all_325sp.R
# Regenerate PAReTT node-age-vs-TimeTree tables for all four dating options
# with a single, consistent taxon matching (full 325-species metadata).
suppressPackageStartupMessages({library(ape); library(dplyr)})

BASE <- "/home/jg2070/Desktop/dtol_review_August"
PUB  <- file.path(BASE, "DToL_phylogenomics_publication_325genomes/01_species_tree")
QC   <- file.path(PUB, "outputs/calibration_qc")
tt_file   <- file.path(BASE, "2026_trees/TimeTree_tree/species_list_timetree.nwk")
tax_file  <- file.path(BASE, "2026_trees/annotation_centromeres/centromere_code_to_species.tsv")

# full 325-sp code<->species map (matches the most TimeTree species)
tax <- read.delim(tax_file, check.names = FALSE)
tax$sp_name  <- paste(tax$genus, tax$species, sep = "_")
tax$tip_code <- sub("[.].*$","", sub("[0-9].*$","", tolower(tax$fasta)))
sp_to_code   <- setNames(tax$tip_code, tax$sp_name)
broad_map <- c(Actinopterygii="Vertebrates",Aves="Vertebrates",Mammalia="Vertebrates",Reptilia="Vertebrates",
  Annelida="Invertebrates",Arthropoda="Invertebrates",Blattodea="Invertebrates",Bryozoa="Invertebrates",
  Chelicerata="Invertebrates",Cnidaria="Invertebrates",Coleoptera="Invertebrates",Dermaptera="Invertebrates",
  Diptera="Invertebrates",Echinodermata="Invertebrates",Ephemeroptera="Invertebrates",Hemiptera="Invertebrates",
  Hymenoptera="Invertebrates",Lepidoptera="Invertebrates",Mollusca="Invertebrates",Nematoda="Invertebrates",
  Nemertea="Invertebrates",Neuroptera="Invertebrates",Odonata="Invertebrates",Plecoptera="Invertebrates",
  Porifera="Invertebrates",Psocodea="Invertebrates",Thecostraca="Invertebrates",Trichoptera="Invertebrates",
  Tunicata="Invertebrates",Fungi="Fungi",Algae="Viridiplantae",Bryophyta="Viridiplantae",
  Dicots="Viridiplantae",Monocots="Viridiplantae",Alveolata="Protists",Discoba="Protists")
code_to_broad <- setNames(broad_map[tax$taxa1], tax$tip_code)
tt_corrections <- c("Ailanthus_altissimus"="Ailanthus_altissima")

tr_tt <- read.tree(tt_file)
tr_tt$tip.label <- ifelse(tr_tt$tip.label %in% names(tt_corrections),
                          tt_corrections[tr_tt$tip.label], tr_tt$tip.label)

make_parrett <- function(tree_file, out_name) {
  tr_our <- read.tree(file.path(PUB, "outputs", tree_file))
  otc <- sub("[0-9].*$","", tolower(tr_our$tip.label))
  matched <- tr_tt$tip.label[tr_tt$tip.label %in% names(sp_to_code)]
  mc <- sp_to_code[matched]; mo <- setNames(tr_our$tip.label, otc)[mc]
  ok <- !is.na(mo) & mo %in% tr_our$tip.label
  matched <- matched[ok]; mo <- mo[ok]
  ttp <- keep.tip(tr_tt, matched)
  op  <- keep.tip(tr_our, mo); op$tip.label <- matched[match(op$tip.label, mo)]
  ctt <- cophenetic.phylo(ttp)/2; cou <- cophenetic.phylo(op)/2
  ord <- rownames(ctt); cou <- cou[ord, ord]
  n <- length(ord); idx <- which(upper.tri(matrix(0,n,n)), arr.ind=TRUE)
  dp <- data.frame(sp1=ord[idx[,1]], sp2=ord[idx[,2]], timetree=ctt[idx], ours=cou[idx])
  dp$tt_node_age <- round(dp$timetree, 2)
  na <- dp %>% group_by(tt_node_age) %>%
    summarise(our_age=mean(ours), n_pairs=n(), sp1=first(sp1), sp2=first(sp2), .groups="drop")
  c1 <- sp_to_code[na$sp1]; c2 <- sp_to_code[na$sp2]
  na$broad <- ifelse(!is.na(code_to_broad[c1]) & !is.na(code_to_broad[c2]) &
                     code_to_broad[c1]==code_to_broad[c2], code_to_broad[c1], "Cross-group")
  out <- na[, c("tt_node_age","broad","our_age","n_pairs")]
  write.table(out, file.path(QC, out_name), sep="\t", quote=FALSE, row.names=FALSE)
  cat(sprintf("%-32s %d matched sp, %d nodes\n", out_name, length(matched), nrow(out)))
}

make_parrett("full_325sp_chronos_over_correlated.nwk",    "parrett_chronos_correlated.tsv")
make_parrett("full_325sp_chronos_over_relaxed.nwk",       "parrett_chronos_relaxed.tsv")
make_parrett("full_325sp_chronos_over_correlated_l01.nwk","parrett_correlated_l01.tsv")   # 62-cal correlated, lambda=0.1
make_parrett("full_325sp_chronos_over_relaxed_l01.nwk",   "parrett_relaxed_l01.tsv")      # 62-cal relaxed, lambda=0.1
make_parrett("full_325sp_chronos_over_correlated_l0.nwk", "parrett_correlated_l0.tsv")    # 62-cal correlated, lambda=0
make_parrett("full_325sp_chronos_over_relaxed_l0.nwk",    "parrett_relaxed_l0.tsv")       # 62-cal relaxed, lambda=0
make_parrett("full_325sp_chronos_over_correlated_l10.nwk","parrett_correlated_l10.tsv")   # 62-cal correlated, lambda=10
make_parrett("full_325sp_chronos_over_relaxed_l10.nwk",   "parrett_relaxed_l10.tsv")      # 62-cal relaxed, lambda=10
make_parrett("full_325sp_chronos_over_clock.nwk",         "parrett_clock.tsv")            # strict clock (no lambda)
make_parrett("full_325sp_chronos_over_discrete_l0.nwk",   "parrett_discrete_l0.tsv")      # discrete, lambda=0
make_parrett("full_325sp_chronos_over_discrete_l01.nwk",  "parrett_discrete_l01.tsv")     # discrete, lambda=0.1
make_parrett("full_325sp_chronos_over_discrete_l1.nwk",   "parrett_discrete_l1.tsv")      # discrete, lambda=1
make_parrett("full_325sp_chronos_over_discrete_l10.nwk",  "parrett_discrete_l10.tsv")     # discrete, lambda=10
make_parrett("full_325sp_chronos_rootonly_ratesmoothed.nwk","parrett_ratesmoothed.tsv")   # root-only, correlated
make_parrett("full_325sp_chronos_rootonly_relaxed.nwk",   "parrett_rootonly_relaxed.tsv") # root-only, relaxed
make_parrett("full_325sp_chronos_nocalib_relaxed.nwk",    "parrett_uncalibrated.tsv")        # no calibration, relaxed
make_parrett("full_325sp_chronos_nocalib_correlated.nwk", "parrett_nocalib_correlated.tsv")  # no calibration, correlated
cat("Done.\n")
