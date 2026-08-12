#!/usr/bin/env Rscript
# Map the two CENP-B box methods onto the 325-sp chronogram.
#   ring 1 = Method 2 (songbird +/-5 flank): Delta = box - flank information (bits)
#   ring 2 = Method 1 (Fachinetti exact IUPAC): broad-motif hits per Mbp (log1p)
# Tips coloured by clade; the top birds are labelled. Species without satellite
# data (163/325) are left blank.
suppressPackageStartupMessages({library(ape); library(ggtree); library(ggtreeExtra)
  library(ggplot2); library(dplyr); library(treeio)})
SAT <- "/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity"
TREE <- "/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/01_species_tree/outputs/full_325sp_chronos_over_correlated.nwk"
tr <- read.tree(TREE)
norm <- function(s) sub("[0-9.].*$","",tolower(s))

flank <- read.delim(file.path(SAT,"figures/cenpb_flank_uncapped_per_species.tsv"))
paper <- read.delim(file.path(SAT,"figures/cenpb_paper_motifs_per_species.tsv"))
dat <- data.frame(label=tr$tip.label, code=norm(tr$tip.label)) %>%
  left_join(flank %>% transmute(code=species, clade, vgroup, delta, box_consensus, subs=subs_vs_canonical, winMbp=win_per_Mbp), by="code") %>%
  left_join(paper %>% transmute(code=species, broadMbp=broad_perMbp), by="code")
dat$clade[is.na(dat$clade)] <- "no data"
dat$delta_pos <- pmax(dat$delta, 0)                       # only positive box>flank shown

ccol <- c(Vertebrates="#0072B2", Invertebrate="#E69F00", Viridiplantae="#009E73",
          Fungi="#CC79A7", Protist="#D55E00", `no data`="grey88")   # Okabe-Ito
# top birds to label
lab <- flank %>% filter(vgroup=="Aves", delta>=0.4) %>% arrange(desc(delta))
labtips <- dat$label[dat$code %in% lab$species]

p <- ggtree(tr, layout="fan", open.angle=12, size=0.2) %<+% dat
p <- p + geom_tippoint(aes(color=clade), size=0.9) +
  scale_color_manual(values=ccol, name="clade") +
  # ring 1: Method 2 flank Delta
  geom_fruit(geom=geom_col, mapping=aes(y=label, x=delta_pos, fill="Δ box−flank (Method 2)"),
             pwidth=0.30, offset=0.06, axis.params=list(axis="x", text.size=1.6, nbreak=3),
             grid.params=list()) +
  # ring 2: Method 1 broad hits per Mbp (log1p)
  geom_fruit(geom=geom_col, mapping=aes(y=label, x=log1p(broadMbp), fill="broad hits/Mbp (Method 1, log)"),
             pwidth=0.30, offset=0.10, axis.params=list(axis="x", text.size=1.6, nbreak=3)) +
  scale_fill_manual(values=c("Δ box−flank (Method 2)"="#117733","broad hits/Mbp (Method 1, log)"="#882255"),
                    name="CENP-B box signal") +
  geom_tiplab2(aes(subset=(label %in% labtips), label=code), size=2.4, offset=0.35, fontface="bold") +
  ggtitle("CENP-B box across the DToL chronogram (162 of 325 species have satellites)",
          subtitle="ring 1 = songbird flank Δ (box>flank, Method 2); ring 2 = Fachinetti broad-motif hits/Mbp (Method 1); birds lead") +
  theme(legend.position="right", plot.title=element_text(face="bold"))

for (ext in c("pdf","png")) {
  ggsave(file.path(SAT, paste0("figures/cenpb_box_tree_325sp.", ext)), p,
         width=11, height=10, dpi=300, bg="white", limitsize=FALSE)
  cat("Saved figures/cenpb_box_tree_325sp.", ext, "\n", sep="")
}
cat("tips with data:", sum(dat$clade!="no data"), "/", nrow(dat), "\n")
cat("labelled birds:", paste(lab$species, round(lab$delta,2), collapse=" | "), "\n")
