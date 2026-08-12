#!/usr/bin/env Rscript
# Map the two CENP-B box methods onto the 325-sp chronogram.
#   ring 1 = Method 2 (songbird +/-5 flank): Delta = box - flank information (bits)
#   ring 2 = Method 1 (Fachinetti exact IUPAC): broad-motif hits per Mbp (log1p)
# Tips coloured by clade; the top birds are labelled. Species without satellite
# data (163/325) are left blank.
suppressPackageStartupMessages({library(ape); library(ggtree); library(ggtreeExtra)
  library(ggplot2); library(dplyr); library(treeio); library(ggnewscale)})
SAT <- "/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity"
TREE <- "/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/01_species_tree/outputs/full_325sp_chronos_over_correlated.nwk"
tr <- read.tree(TREE)
norm <- function(s) sub("[0-9.].*$","",tolower(s))

flank <- read.delim(file.path(SAT,"figures/cenpb_flank_uncapped_per_species.tsv"))
paper <- read.delim(file.path(SAT,"figures/cenpb_paper_motifs_per_species.tsv"))
idn   <- read.delim(file.path(SAT,"figures/cenpb_identity_null.tsv"))   # excess identity over shuffle null
dat <- data.frame(label=tr$tip.label, code=norm(tr$tip.label)) %>%
  left_join(flank %>% transmute(code=species, clade, vgroup, delta, subs=subs_vs_canonical), by="code") %>%
  left_join(paper %>% transmute(code=species, canonMbp=canonical_perMbp, broadMbp=broad_perMbp,
                                degenMbp=degenerated_perMbp), by="code") %>%
  left_join(idn %>% transmute(code=species, excess=excess_canon), by="code")
dat$clade[is.na(dat$clade)] <- "no data"
for (cc in c("canonMbp","broadMbp","degenMbp")) dat[[cc]][is.na(dat[[cc]])] <- 0
dat$delta_pos <- ifelse(is.na(dat$delta), 0, pmax(dat$delta, 0))   # motif signal (flank Δ), Method 2
dat$identity  <- ifelse(is.na(dat$subs), NA, 100*(17-dat$subs)/17) # % identity of consensus motif to canonical
# 'box' = functional (protein-binding); sequence hits are motifs / candidate boxes.
#   candidate box = box-specific (flank Δ≥0.5) AND near-canonical motif (≤2 subs)
dat$candidate <- !is.na(dat$delta) & dat$delta>=0.5 & !is.na(dat$subs) & dat$subs<=2

ccol <- c(Vertebrates="#0072B2", Invertebrate="#E69F00", Viridiplantae="#009E73",
          Fungi="#CC79A7", Protist="#D55E00", `no data`="grey88")   # Okabe-Ito
candtips <- dat$label[dat$candidate]
toptips  <- dat$label[dat$code %in% (flank %>% filter(delta>=0.5 & subs_vs_canonical<=2) %>%
              arrange(desc(delta)) %>% head(6) %>% pull(species))]

p <- ggtree(tr, layout="fan", open.angle=14, size=0.2) %<+% dat +
  geom_tippoint(aes(color=clade), size=0.9) + scale_color_manual(values=ccol, name="clade")
# ring 1 — Method 2 MOTIF signal: bar height = flank Δ, colour = identity to canonical motif
p <- p + new_scale_fill() +
  geom_fruit(geom=geom_col, mapping=aes(y=label, x=delta_pos, fill=identity),
             pwidth=0.26, offset=0.06, axis.params=list(axis="x", text.size=1.5, nbreak=3)) +
  scale_fill_gradientn(colours=c("#FFFFB2","#FD8D3C","#BD0026"), limits=c(70,100), na.value="grey90",
                       name="motif identity to\ncanonical (%)")
# ★ candidate box (box-specific + near-canonical motif)
p <- p + geom_tippoint(aes(subset=(label %in% candtips)), shape=8, size=1.6, stroke=0.5, color="black")
# rings 2-4 — Method 1 exact-motif density (hits/Mbp, log1p): canonical / broad / degenerate
p <- p + new_scale_fill() +
  geom_fruit(geom=geom_col, mapping=aes(y=label, x=log1p(canonMbp), fill="canonical"),
             pwidth=0.11, offset=0.10, axis.params=list(axis="none")) +
  geom_fruit(geom=geom_col, mapping=aes(y=label, x=log1p(broadMbp), fill="broad"),
             pwidth=0.11, offset=0.015, axis.params=list(axis="none")) +
  geom_fruit(geom=geom_col, mapping=aes(y=label, x=log1p(degenMbp), fill="degenerate"),
             pwidth=0.11, offset=0.015, axis.params=list(axis="x", text.size=1.5, nbreak=2)) +
  scale_fill_manual(values=c(canonical="#D55E00", broad="#0072B2", degenerate="#E69F00"),
                    breaks=c("canonical","broad","degenerate"), name="M1 exact motif\n(hits/Mbp, log)") +
  geom_tiplab2(aes(subset=(label %in% toptips), label=code), size=2.3, offset=0.42, fontface="bold") +
  ggtitle("CENP-B motif signal across the DToL chronogram",
          subtitle=paste0("Ring 1 = motif signal (flank Δ, M2), coloured by identity of the consensus motif to canonical. ",
                          "Rings 2–4 = exact-motif density (canonical/broad/degenerate, M1). ",
                          "★ = box-like candidate (",length(candtips)," sp; strongest in birds — near-canonical but missing the CpG).")) +
  theme(legend.position="right", plot.title=element_text(face="bold", size=12),
        plot.subtitle=element_text(size=7.6, colour="grey30"))

for (ext in c("pdf","png")) {
  ggsave(file.path(SAT, paste0("figures/cenpb_box_tree_325sp.", ext)), p,
         width=11, height=10, dpi=300, bg="white", limitsize=FALSE)
  cat("Saved figures/cenpb_box_tree_325sp.", ext, "\n", sep="")
}
cat("tips with data:", sum(dat$clade!="no data"), "/", nrow(dat), "\n")
cat("candidate-box species:", length(candtips), "\n")
print(table(dat$boxclass))
