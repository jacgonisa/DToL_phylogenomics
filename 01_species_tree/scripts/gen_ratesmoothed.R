#!/usr/bin/env Rscript
# Rate-smoothed chronogram: chronos (penalized likelihood) with ONLY the root
# age constrained (Eukaryota_root 1085-1671 My). No other calibrations.
suppressPackageStartupMessages({library(ape); library(dplyr)})

BASE <- "/home/jg2070/Desktop/dtol_review_August"
PUB  <- file.path(BASE, "DToL_phylogenomics_publication_325genomes/01_species_tree")
QC   <- file.path(PUB, "outputs/calibration_qc")
tt_file <- file.path(BASE, "2026_trees/TimeTree_tree/species_list_timetree.nwk")
meta_file <- file.path(BASE, "our_DToL.tsv")

tr <- read.tree(file.path(PUB, "outputs/full_325sp_prepped_nofa.nwk"))
root_node <- Ntip(tr) + 1L
calib <- data.frame(node = root_node, age.min = 1085, age.max = 1671,
                    soft.bounds = FALSE)
ch <- chronos(tr, lambda = 1, model = "correlated", calibration = calib,
              control = chronos.control(iter.max = 1e4, eval.max = 1e4))
class(ch) <- "phylo"
write.tree(ch, file.path(PUB, "outputs/full_325sp_chronos_rootonly_ratesmoothed.nwk"))
cat("root age:", max(node.depth.edgelength(ch)), "\n")

# ── PAReTT tsv vs TimeTree (same logic as the other methods) ──────────────────
meta <- read.delim(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
meta$sp_name  <- paste(meta$Genus, meta$Species, sep = "_")
meta$tip_code <- sub("[.].*$","", sub("[0-9].*$","", tolower(meta$fasta)))
broad_map <- c(Actinopterygii="Vertebrates",Aves="Vertebrates",Mammalia="Vertebrates",Reptilia="Vertebrates",
  Annelida="Invertebrates",Arthropoda="Invertebrates",Blattodea="Invertebrates",Bryozoa="Invertebrates",
  Chelicerata="Invertebrates",Cnidaria="Invertebrates",Coleoptera="Invertebrates",Dermaptera="Invertebrates",
  Diptera="Invertebrates",Echinodermata="Invertebrates",Ephemeroptera="Invertebrates",Hemiptera="Invertebrates",
  Hymenoptera="Invertebrates",Lepidoptera="Invertebrates",Mollusca="Invertebrates",Nematoda="Invertebrates",
  Nemertea="Invertebrates",Neuroptera="Invertebrates",Odonata="Invertebrates",Plecoptera="Invertebrates",
  Porifera="Invertebrates",Psocodea="Invertebrates",Thecostraca="Invertebrates",Trichoptera="Invertebrates",
  Tunicata="Invertebrates",Fungi="Fungi",Algae="Viridiplantae",Bryophyta="Viridiplantae",
  Dicots="Viridiplantae",Monocots="Viridiplantae",Alveolata="Protists",Discoba="Protists")
meta$broad <- broad_map[meta$Group]
sp_to_code <- setNames(meta$tip_code, meta$sp_name)
code_to_broad <- setNames(meta$broad, meta$tip_code)
tt_corrections <- c("Ailanthus_altissimus"="Ailanthus_altissima")

tr_tt <- read.tree(tt_file); tr_our <- ch
tr_tt$tip.label <- ifelse(tr_tt$tip.label %in% names(tt_corrections),
                          tt_corrections[tr_tt$tip.label], tr_tt$tip.label)
our_tip_code <- sub("[0-9].*$","", tolower(tr_our$tip.label))
matched <- tr_tt$tip.label[tr_tt$tip.label %in% names(sp_to_code)]
match_code <- sp_to_code[matched]
match_ourtip <- setNames(tr_our$tip.label, our_tip_code)[match_code]
ok <- !is.na(match_ourtip) & match_ourtip %in% tr_our$tip.label
matched <- matched[ok]; match_ourtip <- match_ourtip[ok]
tr_tt_p <- keep.tip(tr_tt, matched); tr_our_p <- keep.tip(tr_our, match_ourtip)
tr_our_p$tip.label <- matched[match(tr_our_p$tip.label, match_ourtip)]
coph_tt <- cophenetic.phylo(tr_tt_p)/2; coph_our <- cophenetic.phylo(tr_our_p)/2
ord <- rownames(coph_tt); coph_our <- coph_our[ord, ord]
n <- length(ord); idx <- which(upper.tri(matrix(0,n,n)), arr.ind=TRUE)
dp <- data.frame(sp1=ord[idx[,1]], sp2=ord[idx[,2]], timetree=coph_tt[idx], ours=coph_our[idx])
dp$tt_node_age <- round(dp$timetree, 2)
na <- dp %>% group_by(tt_node_age) %>% summarise(our_age=mean(ours), n_pairs=n(),
                                                 sp1=first(sp1), sp2=first(sp2), .groups="drop")
c1 <- sp_to_code[na$sp1]; c2 <- sp_to_code[na$sp2]
na$broad <- ifelse(!is.na(code_to_broad[c1]) & !is.na(code_to_broad[c2]) &
                   code_to_broad[c1]==code_to_broad[c2], code_to_broad[c1], "Cross-group")
out <- na[, c("tt_node_age","broad","our_age","n_pairs")]
write.table(out, file.path(QC,"parrett_ratesmoothed.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("Wrote parrett_ratesmoothed.tsv:", nrow(out), "nodes\n")
cat("shallow R2 (<100My):", round(summary(lm(our_age~tt_node_age, subset(out, tt_node_age<100)))$r.squared,3), "\n")
