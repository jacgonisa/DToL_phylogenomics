#!/usr/bin/env Rscript
# plot_cenpb_tree_325sp.R
# CENP-B box occurrences across the 325-sp chronogram: fan tree coloured by
# clade with an outer ring = number of satellite monomers carrying a CENP-B box.
suppressPackageStartupMessages({
  library(ape); library(ggtree); library(ggtreeExtra); library(ggplot2)
  library(dplyr); library(ggnewscale)
})

BASE <- "/home/jg2070/Desktop/dtol_review_August"
PUB  <- file.path(BASE, "DToL_phylogenomics_publication_325genomes")
ST   <- file.path(PUB, "01_species_tree")
SAT  <- file.path(PUB, "05_satellite_similarity")

tr  <- read.tree(file.path(ST, "outputs/full_325sp_chronos_over_correlated.nwk"))
box <- read.delim(file.path(SAT, "figures/cenpb_box_per_species.tsv"), stringsAsFactors = FALSE)
tax <- read.delim(file.path(BASE, "2026_trees/annotation_centromeres/centromere_code_to_species.tsv"),
                  stringsAsFactors = FALSE)

# tip code -> species/clade
code   <- function(x) sub("[0-9.].*$", "", tolower(x))
tax$code <- code(tax$fasta)
broad_map <- c(Actinopterygii="Vertebrates",Aves="Vertebrates",Mammalia="Vertebrates",
  Reptilia="Vertebrates",Amphibia="Vertebrates",Chondrichthyes="Vertebrates",
  Fungi="Fungi",Algae="Viridiplantae",Bryophyta="Viridiplantae",Dicots="Viridiplantae",
  Monocots="Viridiplantae",Gymnosperms="Viridiplantae",Alveolata="Protist",Discoba="Protist")
tax$clade <- ifelse(tax$taxa1 %in% names(broad_map), broad_map[tax$taxa1], "Invertebrate")

tip_df <- data.frame(label = tr$tip.label, code = code(tr$tip.label), stringsAsFactors = FALSE)
tip_df$clade <- tax$clade[match(tip_df$code, tax$code)]
tip_df$clade[is.na(tip_df$clade)] <- "Invertebrate"
# CENP-B box occurrence count = monomers with a <=2-mismatch box
box$n_box <- round(box$n_monomers * box$pct_box_le2mm / 100)
tip_df$n_box <- box$n_box[match(tip_df$code, box$species)]
tip_df$pct  <- box$pct_box_le2mm[match(tip_df$code, box$species)]

clade_pal <- c(Vertebrates="#1565C0", Invertebrate="#EF6C00",
               Viridiplantae="#2E7D32", Fungi="#6A1B9A", Protist="#C62828")

p <- ggtree(tr, layout = "fan", open.angle = 12, size = 0.25) %<+% tip_df +
  geom_tippoint(aes(colour = clade), size = 0.7) +
  scale_colour_manual(values = clade_pal, name = "Clade", na.translate = FALSE) +
  new_scale_fill() +
  geom_fruit(
    geom = geom_col,
    mapping = aes(y = label, x = n_box, fill = clade),
    orientation = "y", pwidth = 0.38, offset = 0.06,
    axis.params = list(axis = "x", text.size = 2, nbreak = 3),
    grid.params = list()
  ) +
  scale_fill_manual(values = clade_pal, guide = "none") +
  labs(title = "CENP-B box occurrences across the 325-species DToL chronogram",
       subtitle = "Outer bars = no. of satellite monomers carrying a CENP-B box (≤2 mismatches); tips coloured by clade") +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, colour = "grey35"),
        legend.position = c(0.5, 0.52), legend.title = element_text(face="bold"))

for (ext in c("pdf", "png")) {
  out <- file.path(SAT, "figures", paste0("cenpb_box_occurrences_tree_325sp.", ext))
  ggsave(out, p, width = 11, height = 11, dpi = 300, bg = "white")
  cat("Saved:", out, "\n")
}
cat(sprintf("species with box data on tree: %d ; top count: %s (%d boxes)\n",
    sum(!is.na(tip_df$n_box)), tip_df$label[which.max(tip_df$n_box)], max(tip_df$n_box, na.rm=TRUE)))
