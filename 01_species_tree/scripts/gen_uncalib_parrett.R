#!/usr/bin/env Rscript
# Generate parrett_uncalibrated.tsv: raw ML tree node "ages" vs TimeTree.
# Uncalibrated ML branch lengths (subst/site) are root-scaled to the chronos
# root age (strict-clock rescaling) so they plot on the same My axis.
suppressPackageStartupMessages({library(ape); library(dplyr)})

BASE    <- "/home/jg2070/Desktop/dtol_review_August"
PUB_DIR <- file.path(BASE, "DToL_phylogenomics_publication_325genomes/01_species_tree")
QC_DIR  <- file.path(PUB_DIR, "outputs/calibration_qc")
tt_file <- file.path(BASE, "2026_trees/TimeTree_tree/species_list_timetree.nwk")
meta_file <- file.path(BASE, "our_DToL.tsv")
uncal_file <- file.path(PUB_DIR, "outputs/full_325sp_prepped_nofa.nwk")
cal_file   <- file.path(PUB_DIR, "outputs/full_325sp_chronos_over_correlated.nwk")

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

tr_tt   <- read.tree(tt_file)
ml.tree <- read.tree(uncal_file)
# penalized-likelihood rate-smoothing with NO calibration points at all.
# model="relaxed" (cf. Revell 2024, phytools blog "Obtaining a time-calibrated
# ultrametric tree"): here relaxed yields node ages that correlate better with the
# calibrated tree than correlated does.
tr_our <- chronos(ml.tree, lambda = 1, model = "relaxed")
class(tr_our) <- "phylo"
# Keep the node ages in RELATIVE units (root = 1), as in the phytools blog: we do NOT
# rescale to Mya. Only the correlation with TimeTree is reported, which is scale-invariant.
write.tree(tr_our, file.path(PUB_DIR, "outputs/full_325sp_chronos_nocalib_relaxed.nwk"))

tr_tt$tip.label <- ifelse(tr_tt$tip.label %in% names(tt_corrections),
                          tt_corrections[tr_tt$tip.label], tr_tt$tip.label)
our_tip_code <- sub("[0-9].*$","", tolower(tr_our$tip.label))
matched <- tr_tt$tip.label[tr_tt$tip.label %in% names(sp_to_code)]
match_code <- sp_to_code[matched]
match_ourtip <- setNames(tr_our$tip.label, our_tip_code)[match_code]
ok <- !is.na(match_ourtip) & match_ourtip %in% tr_our$tip.label
matched <- matched[ok]; match_ourtip <- match_ourtip[ok]

tr_tt_p  <- keep.tip(tr_tt, matched)
tr_our_p <- keep.tip(tr_our, match_ourtip)
tr_our_p$tip.label <- matched[match(tr_our_p$tip.label, match_ourtip)]

coph_tt  <- cophenetic.phylo(tr_tt_p)/2
coph_our <- cophenetic.phylo(tr_our_p)/2
ord <- rownames(coph_tt); coph_our <- coph_our[ord, ord]
n <- length(ord); idx <- which(upper.tri(matrix(0,n,n)), arr.ind=TRUE)
df_pairs <- data.frame(sp1=ord[idx[,1]], sp2=ord[idx[,2]],
                       timetree=coph_tt[idx], ours=coph_our[idx], stringsAsFactors=FALSE)
df_pairs$tt_node_age <- round(df_pairs$timetree, 2)
node_avg <- df_pairs %>% group_by(tt_node_age) %>%
  summarise(our_age=mean(ours), n_pairs=n(), sp1=first(sp1), sp2=first(sp2), .groups="drop")
c1 <- sp_to_code[node_avg$sp1]; c2 <- sp_to_code[node_avg$sp2]
node_avg$broad <- ifelse(!is.na(code_to_broad[c1]) & !is.na(code_to_broad[c2]) &
                         code_to_broad[c1]==code_to_broad[c2], code_to_broad[c1], "Cross-group")
out <- node_avg[, c("tt_node_age","broad","our_age","n_pairs")]
write.table(out, file.path(QC_DIR,"parrett_uncalibrated.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("Wrote parrett_uncalibrated.tsv:", nrow(out), "nodes; root_age", round(root_age), "\n")
cat("shallow R2 (<100My):", round(summary(lm(our_age~tt_node_age, subset(out, tt_node_age<100)))$r.squared,3), "\n")
