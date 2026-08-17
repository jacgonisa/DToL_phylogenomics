#!/usr/bin/env Rscript
# plot_calibration_combined_benchmark_publication_325sp.R
# Publication version of the combined calibration validation figure:
#   A: Calibration intervals as thin lines (no dot markers on names)
#   B-E: PAReTT node-age concordance for chronos correlated, chronos relaxed,
#        rate-smoothed (treePL) and uncalibrated (strict clock) vs TimeTree.

suppressPackageStartupMessages({
  library(ape); library(dplyr); library(ggplot2); library(patchwork); library(forcats)
  library(ggnewscale)
})

BASE    <- "/home/jg2070/Desktop/dtol_review_August"
PUB_DIR <- file.path(BASE, "DToL_phylogenomics_publication_325genomes/01_species_tree")
QC_DIR  <- file.path(PUB_DIR, "outputs/calibration_qc")
FIG_DIR <- file.path(PUB_DIR, "figures")
TOP_FIG <- file.path(BASE, "DToL_phylogenomics_publication_325genomes/figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TOP_FIG, showWarnings = FALSE, recursive = TRUE)

clade_pal <- c(Root="#212121", Opisthokonta="#607D8B", Metazoa="#1565C0",
               Viridiplantae="#2E7D32", Fungi="#6A1B9A")
broad_pal <- c(Vertebrates="#1565C0", Invertebrates="#EF6C00", Viridiplantae="#2E7D32",
               Fungi="#6A1B9A", Protists="#C62828", `Cross-group`="#BDBDBD")

# ── Load data ─────────────────────────────────────────────────────────────────
clade_order <- c("Root", "Opisthokonta", "Metazoa", "Viridiplantae", "Fungi")

# Real TimeTree data fetched from the TimeTree pairwise API (fetch_timetree_ci.py):
#   tt_median   = "Median Time"  (precomputed_age)
#   tt_ci_low/high = "Range"     (precomputed_ci_low / _high)
df <- read.delim(file.path(QC_DIR, "calibration_nodes_timetree.tsv"),
                 stringsAsFactors = FALSE)
c64 <- read.delim(file.path(QC_DIR, "calib64_constraints_for_panelA.tsv"),
                  stringsAsFactors = FALSE)
df$disp_label <- c64$disp_label[match(df$label, c64$label)]
df$age_min <- as.numeric(c64$age_min[match(df$label, c64$label)])
df$age_max <- as.numeric(c64$age_max[match(df$label, c64$label)])
# all calibration nodes are now TimeTree-sourced (no literature-only constraints)
df$clade <- factor(df$clade, levels = clade_order)
num <- function(x) { x <- as.character(x)
  suppressWarnings(as.numeric(ifelse(x %in% c("None","NA",""), NA, x))) }
df$tt_median   <- num(df$tt_median)
df$tt_adjusted <- num(df$tt_adjusted)
df$tt_ci_low  <- num(df$tt_ci_low)
df$tt_ci_high <- num(df$tt_ci_high)
df$our_age    <- num(df$our_age)
# 0/0 or degenerate ranges are "no range reported" -> drop the line
bad <- is.na(df$tt_ci_low) | is.na(df$tt_ci_high) |
       (df$tt_ci_low == 0 & df$tt_ci_high == 0) | (df$tt_ci_low >= df$tt_ci_high)
df$tt_ci_low[bad] <- NA; df$tt_ci_high[bad] <- NA

clean_label <- function(x) {
  x <- gsub("Metazoa_node16","Avemetatarsalia node",x)
  x <- gsub("Metazoa_node17","Aculeata/Tipuloidea node",x)
  x <- gsub("Vir_node10","Alismatales/Poales",x)
  x <- gsub("Vir_node11","Poaceae/Asteraceae",x)
  x <- gsub("Vir_node12","Fabales/Caryophyllales",x)
  x <- gsub("Vir_node16","Asteraceae node",x)
  x <- gsub("Vir_node17","Moss node",x)
  x <- gsub("Vir_node1","Potamogetonaceae clade",x)
  x <- gsub("Marchantiophyta","Liverwort node",x)
  x <- gsub("Vir_node5","Cyperaceae clade",x)
  x <- gsub("Vir_node8","Caryophyllales/Poales",x)
  gsub("_", " ", x)
}
df$clean_label <- clean_label(df$disp_label)   # no dot markers
df$ord_age <- ifelse(is.na(df$tt_median), df$our_age, df$tt_median)
df <- df %>% arrange(clade, ord_age)
df$node_f <- fct_rev(factor(df$clean_label, levels = df$clean_label))
# y-axis label colours by clade (aligned to factor levels, bottom->top)
ylev  <- levels(df$node_f)
ycols <- unname(clade_pal[as.character(df$clade[match(ylev, df$clean_label)])])

# ── Panel A: TimeTree range (line) + median (diamond) + calibrated age (circle)
p_a <- ggplot(df) +
  # invisible layer purely to build the clade colour legend (y-axis is coloured too)
  geom_point(aes(x = ord_age, y = node_f, colour = clade), alpha = 0, na.rm = TRUE) +
  scale_colour_manual(values = clade_pal, name = "Clade",
                      guide = guide_legend(order = 1,
                        override.aes = list(alpha = 1, shape = 15, size = 4.2))) +
  new_scale_colour() +
  # TimeTree range (CI): one line per node, where TimeTree reports a range
  geom_segment(aes(x = tt_ci_low, xend = tt_ci_high, y = node_f, yend = node_f,
                   colour = "TimeTree range (CI)"),
               linewidth = 0.95, lineend = "butt", na.rm = TRUE) +
  scale_colour_manual(name = NULL,
                      values = c("TimeTree range (CI)" = "#d62728"),
                      guide = guide_legend(order = 2,
                        override.aes = list(linewidth = 1.5))) +
  # TimeTree median (red diamond) + calibrated age (black circle)
  geom_point(aes(x = tt_median, y = node_f, fill = "TimeTree median"),
             shape = 23, size = 2.9, colour = "white", stroke = 0.35, na.rm = TRUE) +
  geom_point(aes(x = our_age, y = node_f, fill = "Calibrated age (chronos correlated)"),
             shape = 21, size = 2.9, colour = "white", stroke = 0.35, na.rm = TRUE) +
  scale_fill_manual(name = NULL,
                    breaks = c("Calibrated age (chronos correlated)", "TimeTree median"),
                    values = c("Calibrated age (chronos correlated)" = "black",
                               "TimeTree median" = "#d62728"),
                    guide = guide_legend(order = 3, override.aes = list(
                      shape = c(21, 23), size = 3.6))) +
  scale_x_continuous(name = "Age (Mya)", expand = expansion(mult = c(0.01, 0.10))) +
  scale_y_discrete(name = NULL) +
  labs(tag = "A",
       title    = sprintf("Calibration nodes: TimeTree range vs calibrated age (%d nodes)", nrow(df)),
       subtitle = "Red diamond = TimeTree median  |  red line = TimeTree range (constraint for most nodes)  |  black circle = calibrated age  |  y labels by clade") +
  theme_bw(base_size = 13) +
  theme(axis.text.y = element_text(size = 9, colour = ycols),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        legend.position = "right", legend.text = element_text(size = 11.5),
        legend.title = element_text(size = 13, face = "bold"),
        plot.tag = element_text(face = "bold", size = 19),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 10, colour = "grey35"),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())

# ── Panels B-E: PAReTT node-age comparison, one panel per dating method ───────
build_parrett <- function(fname, tag, ttl) {
  d <- read.delim(file.path(QC_DIR, fname), stringsAsFactors = FALSE)
  d$broad <- factor(d$broad, levels = names(broad_pal))
  r2  <- round(summary(lm(our_age ~ tt_node_age, d))$r.squared, 3)
  pear <- round(cor(d$tt_node_age, d$our_age, method = "pearson"), 3)
  mx  <- max(c(d$tt_node_age, d$our_age), na.rm = TRUE)
  ggplot(d, aes(tt_node_age, our_age, colour = broad)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey20", linewidth = 0.95) +
    geom_smooth(aes(tt_node_age, our_age), inherit.aes = FALSE, method = "lm",
                se = TRUE, colour = "#1565C0", linewidth = 0.85, alpha = 0.12) +
    geom_point(aes(size = log10(n_pairs + 1)), alpha = 0.72, shape = 16) +
    annotate("text", x = mx * 0.04, y = mx * 0.90, hjust = 0, size = 4.3,
             colour = "grey20", lineheight = 1.4,
             label = sprintf("Pearson = %.3f\nR² = %.3f", pear, r2)) +
    scale_colour_manual(values = broad_pal, name = "Group", na.value = "grey60") +
    scale_size_continuous(name = "log10(pairs)", range = c(1.8, 5.5)) +
    scale_x_continuous(name = "TimeTree node age (Mya)", limits = c(0, mx * 1.03), expand = expansion(0)) +
    scale_y_continuous(name = "Estimated node age (Mya)", limits = c(0, mx * 1.03), expand = expansion(0)) +
    coord_fixed(ratio = 1) +
    labs(tag = tag, title = ttl,
         subtitle = sprintf("%d node-age comparisons; dashed = 1:1 line", nrow(d))) +
    theme_bw(base_size = 13) +
    theme(plot.tag = element_text(face = "bold", size = 19),
          plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(size = 10, colour = "grey35"),
          axis.title = element_text(size = 14, face = "bold"),
          legend.position = "right", legend.text = element_text(size = 11.5),
          legend.title = element_text(size = 13, face = "bold"),
          axis.text = element_text(size = 12), panel.grid.minor = element_blank())
}
p_b <- build_parrett("parrett_chronos_correlated.tsv", "B", "chronos correlated")
p_c <- build_parrett("parrett_chronos_relaxed.tsv",    "C", "chronos relaxed")
p_d <- build_parrett("parrett_ratesmoothed.tsv",       "D", "rate-smoothed (root only)")
p_e <- build_parrett("parrett_uncalibrated.tsv",       "E", "uncalibrated (no calibration)")

# ── Combine ───────────────────────────────────────────────────────────────────
# collect the shared Group / log10(pairs) legends of the B-E row into one
p_row <- (p_b | p_c | p_d | p_e) + plot_layout(guides = "collect")
p_combined <- (p_a / p_row) +
  plot_layout(heights = c(2.2, 2.0)) +
  plot_annotation(
    title   = "Supplementary Figure - Calibration validation: 325-species DToL chronogram",
    caption = paste0(
      "A: The 62 calibration nodes (y-axis labels coloured by clade). Red diamond = TimeTree ",
      "median (precomputed_age); red line = TimeTree range (precomputed_ci_low-high) - the ",
      "imposed constraint for nodes with a real range (else median +/-20%); black circle = ",
      "calibrated age (chronos correlated). All 62 constraints are TimeTree-derived.  ",
      "B-E: PAReTT node-age concordance with TimeTree (all shared species pairs) for ",
      "chronos correlated (B), chronos relaxed (C), rate-smoothed chronos with only the ",
      "root age fixed (D) and uncalibrated penalized-likelihood rate-smoothing with no ",
      "calibration points (E, chronos relaxed, lambda=1); dashed = 1:1."
    ),
    theme = theme(plot.title = element_text(face = "bold", size = 17),
                  plot.caption = element_text(size = 10, colour = "grey35", hjust = 0),
                  plot.background = element_rect(fill = "white", colour = NA))
  )

for (ext in c(".pdf", ".png")) {
  out <- file.path(FIG_DIR, paste0("calibration_combined_qc_benchmark_325sp_publication", ext))
  ggsave(out, p_combined, width = 22, height = 23,
         dpi = if (ext == ".png") 320 else 100, bg = "white")
  cat("Saved:", basename(out), "\n")
}
file.copy(file.path(FIG_DIR, "calibration_combined_qc_benchmark_325sp_publication.pdf"),
          file.path(TOP_FIG, "calibration_combined_qc_benchmark_325sp_publication.pdf"), overwrite = TRUE)
file.copy(file.path(FIG_DIR, "calibration_combined_qc_benchmark_325sp_publication.png"),
          file.path(TOP_FIG, "calibration_combined_qc_benchmark_325sp_publication.png"), overwrite = TRUE)
cat("\nDone.\n")
