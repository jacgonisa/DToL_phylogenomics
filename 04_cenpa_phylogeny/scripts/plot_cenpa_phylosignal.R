#!/usr/bin/env Rscript
# Species tree fan plot with CENPA copy number + Blomberg's K / Pagel's lambda
suppressPackageStartupMessages({
  library(ape); library(ggtree); library(ggtreeExtra); library(ggplot2)
  library(dplyr); library(phytools); library(ggnewscale)
})

BASE    <- "/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes"
TREE_F  <- file.path(BASE, "01_species_tree/outputs/full_325sp_calibrated.nwk")
CN_F    <- file.path(BASE, "04_cenpa_phylogeny/cenpa_copy_number_table.tsv")
FIG_DIR <- file.path(BASE, "04_cenpa_phylogeny/figures")

tr <- read.tree(TREE_F)
cn <- read.delim(CN_F, stringsAsFactors = FALSE)

# ── match tip labels (strip .fa, normalise) ───────────────────────────────────
strip <- function(x) sub("\\.fa$", "", x)
tr$tip.label <- strip(tr$tip.label)

norm_key <- function(x) {
  x <- tolower(sub("\\.fa$", "", x))
  x <- sub("\\.(grcg7b|hap[0-9]+(\\.[0-9]+)?)$", "", x)
  sub("\\.[0-9]+$", "", x)
}
cn$key <- norm_key(cn$species_code_full)
tip_df <- data.frame(label = tr$tip.label, key = norm_key(tr$tip.label),
                     stringsAsFactors = FALSE)
tip_df <- left_join(tip_df, cn[, c("key","n_loci_pass40pct","centromere_type")],
                    by = "key")
tip_df$copies <- as.integer(tip_df$n_loci_pass40pct)
tip_df$copies[is.na(tip_df$copies)] <- 0L
tip_df$ct <- ifelse(is.na(tip_df$centromere_type), "Unknown", tip_df$centromere_type)

cat(sprintf("Tips matched: %d / %d  (copies > 1: %d)\n",
            sum(!is.na(tip_df$n_loci_pass40pct)), nrow(tip_df),
            sum(tip_df$copies > 1, na.rm = TRUE)))

ct_pal <- c(Satellite="#ff2d87", Transposon="#4f84e8", Holocentric="#2d7d32",
            "Satellite/transposon"="#8338ec", Unknown="#aaaaaa",
            "Monocentric sequence unknown"="#cccccc")

# ── Blomberg's K + Pagel's lambda ────────────────────────────────────────────
copies_vec <- setNames(tip_df$copies, tip_df$label)

# prune to species with data
has_data <- !is.na(copies_vec)
tr_pruned <- drop.tip(tr, tr$tip.label[!has_data])
x <- copies_vec[tr_pruned$tip.label]
cat(sprintf("\nPhylogenetic signal on %d tips\n", length(x)))

K_res <- phylosig(tr_pruned, x, method = "K", test = TRUE, nsim = 999)
L_res <- phylosig(tr_pruned, x, method = "lambda", test = TRUE)
cat(sprintf("Blomberg's K = %.4f  (p = %.4f, %d permutations)\n",
            K_res$K, K_res$P, K_res$nsim))
cat(sprintf("Pagel's λ   = %.4f  (LRT p = %.4g)\n",
            L_res$lambda, L_res$P))

signal_note <- sprintf(
  "Blomberg's K = %.3f (p = %.3f)   Pagel's λ = %.3f (p = %.4f)",
  K_res$K, K_res$P, L_res$lambda, L_res$P
)

# ── tree plot ─────────────────────────────────────────────────────────────────
# copy-number category for ring
tip_df$cn_cat <- cut(
  tip_df$copies,
  breaks = c(-1, 0, 1, 2, 3, 4, Inf),
  labels = c("0", "1", "2", "3", "4", ">4")
)
ring_pal <- c("0" = "#dddddd", "1" = "#ffffb2",
              "2" = "#fecc5c", "3" = "#fd8d3c",
              "4" = "#e31a1c", ">4" = "#7b0000")

p <- ggtree(tr, layout = "fan", size = 0.15, color = "grey55") %<+% tip_df +
  geom_tippoint(aes(color = ct), size = 1.2, alpha = 0.75) +
  scale_color_manual(values = ct_pal, name = "Centromere type",
                     na.value = "#aaaaaa", guide = guide_legend(order = 1)) +
  new_scale_color() +
  geom_fruit(
    geom     = geom_tile,
    mapping  = aes(y = label, x = 1, fill = cn_cat),
    offset   = 0.07, pwidth = 0.06
  ) +
  scale_fill_manual(values = ring_pal, name = "CENPA copies",
                    guide  = guide_legend(order = 2)) +
  labs(
    title    = "CENPA copy number on the 325-sp calibrated species tree",
    subtitle = signal_note
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40"),
    legend.position = "right", legend.text = element_text(size = 8)
  )

ggsave(file.path(FIG_DIR, "cenpa_copies_species_tree.pdf"), p, width = 14, height = 14)
ggsave(file.path(FIG_DIR, "cenpa_copies_species_tree.png"), p, width = 14, height = 14, dpi = 300)
cat("Saved: cenpa_copies_species_tree.{pdf,png}\n")
