#!/usr/bin/env Rscript
# plot_cenpa_phylogeny_325sp.R
# Figure 1D — CenpA + H3 + Archaea phylogeny for the 325-sp dataset.
#
# ── Tree ──────────────────────────────────────────────────────────────────────
#   iqtree2_cenpa430_H3_archaea10_325sp_bnni.treefile
#   1329 tips: 418 CenpA | 901 H3-like | 10 archaea
#   Model: Q.pfam+R7   |   UFBoot 1000 with -bnni (converged, ~38h)
#
# ── Tip classification source ─────────────────────────────────────────────────
#   iTOL_bnni/iTOL_centromere_symbols.txt
#     → per CenpA tip: colour → centromere architecture
#       #ff006e = Satellite        #3a86ff = Transposon
#       #2d7d32 = Holocentric      #8338ec = Satellite/transposon
#       #808080 = Cryptic (Unknown + Monocentric seq. unknown)
#
#   iTOL_bnni/iTOL_taxonomic_clade.txt
#     → per tip: colour → broad taxonomic group
#       #FF6F00 = Insects          #1565C0 = Vertebrates
#       #EF6C00 = Invertebrates    #2E7D32 = Viridiplantae
#       #6A1B9A = Fungi            #C62828 = Protists
#       #BDBDBD = H3               #212121 = Archaea
#
#   H3 and Archaea tips have no centromere symbol in the symbols file.
#
# ── Archaea (outgroup) ────────────────────────────────────────────────────────
#   10 hardcoded accessions used to root the tree at their MRCA.
#
# Outputs → figures/panel_D_cenpa_phylogeny_325sp.{pdf,png}
#           figures/ (top-level)

suppressPackageStartupMessages({
  library(ape); library(ggtree); library(ggplot2)
  library(dplyr); library(ggnewscale)
})

BASE     <- "/home/jg2070/Desktop/dtol_review_August"
PUB_DIR  <- file.path(BASE, "DToL_phylogenomics_publication_325genomes")
TREE_DIR <- file.path(PUB_DIR, "04_cenpa_phylogeny")
ITOL_DIR <- file.path(TREE_DIR, "iTOL_bnni")
FIG_DIR  <- file.path(TREE_DIR, "figures")
TOP_FIG  <- file.path(PUB_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TOP_FIG, showWarnings = FALSE, recursive = TRUE)

# ── Input files ───────────────────────────────────────────────────────────────
TREEFILE <- file.path(TREE_DIR, "iqtree2_cenpa430_H3_archaea10_325sp_bnni.treefile")
SYMS_TXT <- file.path(ITOL_DIR, "iTOL_centromere_symbols.txt")
CLADE_TXT<- file.path(ITOL_DIR, "iTOL_taxonomic_clade.txt")

# ── Archaea accessions (outgroup) ─────────────────────────────────────────────
archaea_ids <- c(
  "BAD86478.1", "KKK41979.1", "KXH71038.1", "OIO41945.1", "OIO61677.1",
  "OLS16336.1", "OLS18261.1", "OLS21974.1", "OLS22332.1", "OLS24873.1"
)

# ── Load and root tree ────────────────────────────────────────────────────────
cat("Loading tree:", basename(TREEFILE), "\n")
tr_raw <- read.tree(TREEFILE)
cat(sprintf("  %d tips, %d internal nodes\n", Ntip(tr_raw), tr_raw$Nnode))
tr <- root(tr_raw,
           node         = getMRCA(tr_raw, archaea_ids),
           resolve.root = TRUE)
cat(sprintf("  Rooted at archaea MRCA\n"))

# ── Parse iTOL centromere symbols (CenpA tips only) ────────────────────────────
parse_itol_data <- function(path, col_names) {
  lines <- readLines(path)
  data_start <- which(lines == "DATA") + 1L
  data_lines <- lines[data_start:length(lines)]
  data_lines <- data_lines[nzchar(data_lines) & !startsWith(data_lines, "#")]
  as.data.frame(do.call(rbind, strsplit(data_lines, "\t")),
                stringsAsFactors = FALSE) %>%
    setNames(col_names[seq_len(ncol(.))])
}

sym_df <- parse_itol_data(SYMS_TXT, c("label", "size", "shape", "colour", "fill", "alpha"))
# colour → centromere type
colour_to_type <- c(
  "#ff006e" = "Satellite",
  "#3a86ff" = "Transposon",
  "#2d7d32" = "Holocentric",
  "#8338ec" = "Satellite/transposon",
  "#808080" = "Cryptic"
)
sym_df$cenpa_type <- colour_to_type[sym_df$colour]
cat(sprintf("  CenpA symbols loaded: %d tips\n", nrow(sym_df)))

# ── Parse iTOL taxonomic clade ────────────────────────────────────────────────
clade_df <- parse_itol_data(CLADE_TXT, c("label", "colour", "taxon_name"))
colour_to_broad <- c(
  "#FF6F00" = "Insects",
  "#1565C0" = "Vertebrates",
  "#EF6C00" = "Invertebrates",
  "#2E7D32" = "Viridiplantae",
  "#6A1B9A" = "Fungi",
  "#C62828" = "Protists",
  "#BDBDBD" = "H3",
  "#212121" = "Archaea",
  "#9E9E9E" = "Other"
)
clade_df$broad <- colour_to_broad[clade_df$colour]
cat(sprintf("  Taxonomic clade annotations: %d tips\n", nrow(clade_df)))

# ── Build per-tip metadata ─────────────────────────────────────────────────────
tips <- tr$tip.label

tip_df <- data.frame(label = tips, stringsAsFactors = FALSE) %>%
  left_join(sym_df[, c("label", "cenpa_type")], by = "label") %>%
  left_join(clade_df[, c("label", "broad")],    by = "label") %>%
  mutate(
    # Tips not in symbols file = H3 or Archaea
    tip_class = case_when(
      !is.na(cenpa_type)       ~ cenpa_type,
      label %in% archaea_ids   ~ "Archaea",
      TRUE                     ~ "H3"
    )
  )

cat("Tip classification summary:\n")
print(sort(table(tip_df$tip_class), decreasing = TRUE))
cat("\nBroad group summary (CenpA only):\n")
cenpa_broad <- tip_df %>% filter(!is.na(cenpa_type))
print(sort(table(cenpa_broad$broad), decreasing = TRUE))

cat("Note: gene tree — CenpA sequences are NOT monophyletic by taxon.\n")
cat("Using outer ring for taxonomic group instead of clade highlights.\n")

# ── Colour/shape maps ─────────────────────────────────────────────────────────
pal_cenpa <- c(
  Satellite             = "#ff006e",
  Transposon            = "#3a86ff",
  Holocentric           = "#2d7d32",
  "Satellite/transposon"= "#8338ec",
  Cryptic               = "#808080"
)
cenpa_levels <- names(pal_cenpa)

count_type <- function(tp) sum(tip_df$tip_class == tp)

# ── Build ggtree ──────────────────────────────────────────────────────────────
cat("Building circular ggtree...\n")

p <- ggtree(tr, layout = "circular", linewidth = 0.08, colour = "grey70") %<+% tip_df

# Layer 1: H3-like (faint grey background)
p <- p +
  geom_tippoint(
    data  = function(d) filter(d, tip_class == "H3"),
    color = "#BDBDBD", shape = 16, size = 0.28, alpha = 0.30
  )

# Layer 2: Archaea outgroup
p <- p +
  geom_tippoint(
    data  = function(d) filter(d, tip_class == "Archaea"),
    color = "#212121", shape = 15, size = 1.2, alpha = 0.9
  )

# Layer 3: CenpA — coloured by centromere type, triangles for Holocentric
p <- p + new_scale_colour() +
  geom_tippoint(
    data = function(d) filter(d, tip_class %in% cenpa_levels) %>%
                       mutate(tip_class = factor(tip_class, levels = cenpa_levels),
                              pt_shape  = ifelse(tip_class == "Holocentric", 17L, 16L)),
    aes(colour = tip_class, shape = I(pt_shape)),
    size = 0.85, alpha = 0.95
  ) +
  scale_colour_manual(
    values = pal_cenpa,
    name   = "Centromere type",
    labels = c(
      Satellite             = sprintf("Satellite (n=%d)",             count_type("Satellite")),
      Transposon            = sprintf("Transposon (n=%d)",            count_type("Transposon")),
      Holocentric           = sprintf("Holocentric (n=%d)",           count_type("Holocentric")),
      "Satellite/transposon"= sprintf("Sat/Trans (n=%d)",             count_type("Satellite/transposon")),
      Cryptic               = sprintf("Cryptic/Unknown (n=%d)",       count_type("Cryptic"))
    ),
    guide = guide_legend(
      order        = 1,
      override.aes = list(size = 2.5, shape = c(16, 16, 17, 16, 16))
    )
  )

# ── Outer ring: broad taxonomic group (CenpA tips only) ───────────────────────
# Gene tree — CenpA sequences are NOT monophyletic by taxon, so clade highlights
# spanning the whole tree are misleading. Instead, add a second outer ring layer.
pal_broad <- c(
  Insects       = "#FF6F00",
  Vertebrates   = "#1565C0",
  Invertebrates = "#EF6C00",
  Viridiplantae = "#2E7D32",
  Fungi         = "#6A1B9A",
  Protists      = "#C62828",
  Other         = "#9E9E9E"
)

p <- p + new_scale_colour() +
  geom_tippoint(
    data = function(d) filter(d, !is.na(cenpa_type), !is.na(broad), broad %in% names(pal_broad)) %>%
                       mutate(broad = factor(broad, levels = names(pal_broad))),
    aes(colour = broad),
    size = 2.2, alpha = 0.55, shape = 16
  ) +
  scale_colour_manual(
    values = pal_broad,
    name   = "Taxonomic group\n(CenpA tips)",
    guide  = guide_legend(order = 2, override.aes = list(size = 2.5, alpha = 0.9))
  )

# ── Legend entries for H3 and Archaea ────────────────────────────────────────
# Manual legend text (appended via annotate on a blank panel — use caption instead)
n_cenpa   <- sum(tip_df$tip_class %in% cenpa_levels)
n_h3      <- sum(tip_df$tip_class == "H3")
n_archaea <- sum(tip_df$tip_class == "Archaea")

p <- p +
  labs(
    title   = sprintf(
      "CenpA + H3 + Archaea phylogeny  (IQ-TREE2 bnni, Q.pfam+R7, UFBoot 1000)"),
    caption = sprintf(
      "325-sp dataset  |  %d tips total: %d CenpA | %d H3-like (grey) | %d Archaea (black squares, outgroup)\n"
      , Ntip(tr), n_cenpa, n_h3, n_archaea)
  ) +
  theme(
    plot.title        = element_text(size = 9.5, face = "bold", hjust = 0.5,
                                     margin = margin(b = 2)),
    plot.caption      = element_text(size = 7, hjust = 0.5, color = "grey40"),
    legend.position   = c(0.86, 0.43),
    legend.key.size   = unit(0.30, "cm"),
    legend.text       = element_text(size = 6.5),
    legend.title      = element_text(size = 7.5, face = "bold"),
    legend.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.3),
    plot.background   = element_rect(fill = "white", colour = NA)
  )

# ── Save ──────────────────────────────────────────────────────────────────────
cat("Saving figures...\n")
ggsave(file.path(FIG_DIR, "panel_D_cenpa_phylogeny_325sp.pdf"), p,
       width = 17, height = 17)
ggsave(file.path(FIG_DIR, "panel_D_cenpa_phylogeny_325sp.png"), p,
       width = 17, height = 17, dpi = 200)
file.copy(file.path(FIG_DIR, "panel_D_cenpa_phylogeny_325sp.pdf"),
          file.path(TOP_FIG, "panel_D_cenpa_phylogeny_325sp.pdf"), overwrite = TRUE)
file.copy(file.path(FIG_DIR, "panel_D_cenpa_phylogeny_325sp.png"),
          file.path(TOP_FIG, "panel_D_cenpa_phylogeny_325sp.png"), overwrite = TRUE)
cat("Saved: panel_D_cenpa_phylogeny_325sp.{pdf,png}\n")
cat("Output:", FIG_DIR, "\n")
