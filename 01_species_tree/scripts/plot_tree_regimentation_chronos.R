#!/usr/bin/env Rscript
# plot_tree_v5_325sp.R
# Produces the publication-style v5 circular annotated tree for the 325-sp dataset.
# Replicates centromere_annotation_tree_FASTSPECIES_327sp_withTetJugo_v5 style:
#   - Background clade highlights
#   - Branch subclade colour (pure-descendant edges)
#   - Tip symbols: circle=Sat/Trans/Unk, triangle=Holocentric, rectangle=Mixed
#   - Subclade arch labels
#   - Concentric time rings (200 Mya intervals)
#   - Calibration point dots with index numbers
#
# Output -> 01_species_tree/figures/centromere_regimentation_HOR_tree_325sp.{pdf,png}

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(ggtree)
  library(ggtreeExtra)
  library(ggplot2)
  library(ggnewscale)
  library(readxl)
  library(dplyr)
  library(tibble)
  library(stringr)
})

BASE      <- "/home/jg2070/Desktop/dtol_review_August"
PUB_DIR   <- file.path(BASE, "DToL_phylogenomics_publication_325genomes")
TR_DIR    <- file.path(PUB_DIR, "01_species_tree/outputs")
FIG_DIR   <- file.path(PUB_DIR, "01_species_tree/figures")
TOP_FIG   <- file.path(PUB_DIR, "figures")
excel_path <- file.path(BASE, "2026_trees/annotation_centromeres/DTOL_327_master_March.xlsx")
iTOL_path  <- file.path(PUB_DIR, "02_ASR/inputs/full/branch_symbol_anno.tsv")
tree_path  <- file.path(TR_DIR, "full_325sp_calibrated_correlatedlambda01.nwk")  # chronos correlated, lambda=0.1, 62 calibration points
# Aggregation for poly-/di-typic species (multiple centromeric families):
#   ian   = dominant family (highest copy number)          [Ian's preference]
#   piotr = genomic average weighted by family frequency   [Piotr's preference]
AGG <- tolower(Sys.getenv("AGG", "ian"))
stopifnot(AGG %in% c("ian", "piotr"))
agg_lab <- if (AGG == "ian") "dominant array" else "freq-weighted mean"
out_pdf    <- file.path(FIG_DIR, sprintf("centromere_regimentation_HOR_tree_325sp_%s.pdf", AGG))
out_png    <- file.path(FIG_DIR, sprintf("centromere_regimentation_HOR_tree_325sp_%s.png", AGG))

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TOP_FIG, showWarnings = FALSE, recursive = TRUE)

# ── Helpers ───────────────────────────────────────────────────────────────────
strip_ext <- function(x) str_replace(x, "\\.(fa|fasta)$", "")
label_fix_map <- c("daTanVulg1.hap2.1.fasta" = "daTanVulg1.hap1.1.fa")
normalize_id <- function(x) {
  x |>
    str_replace("\\.fasta\\.F2B\\.ChrOnly\\.fa$", ".fa") |>
    strip_ext() |>
    tolower()
}

# ── Load tree ─────────────────────────────────────────────────────────────────
tree_calibrated <- read.tree(tree_path)
# .nwk is already in Mya; no branch-length rescaling needed
tree_scale <- 1

fixed <- ifelse(
  tree_calibrated$tip.label %in% names(label_fix_map),
  unname(label_fix_map[tree_calibrated$tip.label]),
  tree_calibrated$tip.label
)
tree_calibrated$tip.label <- normalize_id(fixed)
norm_fix <- setNames(normalize_id(unname(label_fix_map)), normalize_id(names(label_fix_map)))
tree_calibrated$tip.label <- ifelse(
  tree_calibrated$tip.label %in% names(norm_fix),
  unname(norm_fix[tree_calibrated$tip.label]),
  tree_calibrated$tip.label
)
tree_plot <- tree_calibrated
# tree_scale = 1, so tree_plot == tree_calibrated (no scaling)

cat("Tips:", Ntip(tree_plot), "\n")

# ── Load metadata ─────────────────────────────────────────────────────────────
meta <- suppressMessages(read_excel(excel_path, sheet = "Sheet1")) %>%
  select(fasta, genus, species, classification, taxa1, taxa2) %>%
  mutate(
    fasta_fixed = ifelse(fasta %in% names(label_fix_map), unname(label_fix_map[fasta]), fasta),
    fasta_base  = normalize_id(fasta_fixed)
  )

itol <- read.delim(iTOL_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
arch_map <- itol %>% select(architecture, symbol, colour) %>% distinct()

class_to_arch <- tibble(
  classification = c("Unknown", "Transposon", "Satellite", "Holocentric",
                     "Satellite/transposon", "Monocentric sequence unknown"),
  architecture   = c("Unknown", "Transposon", "Satellite", "Holocentric no satellite",
                     "Mixed", "Unknown")
)

class_style <- class_to_arch %>% left_join(arch_map, by = "architecture")
meta2 <- meta %>% left_join(class_style, by = "classification")

plot_df <- tibble(label = tree_calibrated$tip.label) %>%
  left_join(meta2, by = c("label" = "fasta_base")) %>%
  mutate(
    classification = ifelse(is.na(classification), "Unknown", classification),
    architecture   = ifelse(is.na(architecture),   "Unknown", architecture),
    symbol         = ifelse(is.na(symbol),          "rectangle", symbol),
    colour         = ifelse(is.na(colour),           "#808080", colour),
    display_label  = ifelse(
      !is.na(genus) & !is.na(species) & genus != "" & species != "",
      paste(genus, species), label
    )
  ) %>%
  mutate(
    colour = ifelse(classification == "Holocentric", "#2d7d32", colour),
    symbol = ifelse(classification == "Holocentric", "triangle left", symbol)
  )

shape_map <- c("circle" = 21, "rectangle" = 22, "triangle left" = 24)

# ── Broad group + subclade assignment ─────────────────────────────────────────
plot_df <- plot_df %>%
  mutate(
    broad_group = case_when(
      taxa1 %in% c("Actinopterygii", "Aves", "Reptilia", "Mammalia", "Amphibia") ~ "Vertebrates",
      taxa1 %in% c("Monocots", "Dicots", "Bryophyta", "Algae") ~ "Viridiplantae",
      taxa1 == "Fungi" ~ "Fungi",
      taxa1 %in% c("Alveolata", "Discoba") ~ "Protists",
      is.na(taxa1) | taxa1 == "" ~ "Other",
      TRUE ~ "Invertebrates"
    ),
    fungi_lineage = ifelse(
      broad_group == "Fungi" & !is.na(taxa2) & taxa2 != "",
      str_trim(sub(",.*$", "", taxa2)),
      NA_character_
    ),
    inv_detail_raw = case_when(
      broad_group == "Invertebrates" & taxa1 == "Arthropoda" & !is.na(taxa2) & taxa2 != "" ~ taxa2,
      broad_group == "Invertebrates" & !is.na(taxa1) & taxa1 != "" ~ taxa1,
      TRUE ~ NA_character_
    )
  )

inv_counts  <- sort(table(plot_df$inv_detail_raw[plot_df$broad_group == "Invertebrates" & !is.na(plot_df$inv_detail_raw)]), decreasing = TRUE)
inv_top     <- names(inv_counts)[seq_len(min(8, length(inv_counts)))]
fungi_counts <- sort(table(plot_df$fungi_lineage[plot_df$broad_group == "Fungi" & !is.na(plot_df$fungi_lineage)]), decreasing = TRUE)
fungi_top   <- names(fungi_counts)[seq_len(min(4, length(fungi_counts)))]
inv_force   <- c(
  "Chelicerata", "Thecostraca", "Entomobryomorpha", "Symphypleona",
  "Odonata", "Ephemeroptera", "Plecoptera", "Dermaptera", "Blattodea",
  "Psocodea", "Neuroptera", "Hemiptera", "Trichoptera",
  "Coleoptera", "Diptera", "Lepidoptera", "Hymenoptera"
)

plot_df <- plot_df %>%
  mutate(
    subclade = case_when(
      broad_group %in% c("Vertebrates", "Viridiplantae", "Protists") ~ taxa1,
      broad_group == "Invertebrates" & !is.na(inv_detail_raw) ~ inv_detail_raw,
      broad_group == "Invertebrates" ~ "Other invertebrates",
      broad_group == "Fungi" & !is.na(fungi_lineage) ~ fungi_lineage,
      broad_group == "Fungi" ~ "Other fungi",
      TRUE ~ "Other"
    )
  )

# ── Time rings — geological age (0 Mya at tips, root_age Mya at root) ────────
tree_df  <- ggtree::fortify(tree_plot)
max_x    <- max(tree_df$x[tree_df$isTip], na.rm = TRUE)   # root-to-tip = root age
ring_step <- 200
# Geological ages to mark (0 = present at tips, counting back toward root)
geo_ages  <- seq(0, floor(max_x / ring_step) * ring_step, by = ring_step)
# ggtree x = distance from root; age_from_present = max_x - x
ring_x      <- max_x - geo_ages
ring_labels <- paste0(geo_ages, " Mya")

# ── Background and branch colour palettes ─────────────────────────────────────
bg_palette <- c(
  "Invertebrates" = "#f6ddce",
  "Vertebrates"   = "#d8edf5",
  "Viridiplantae" = "#d9efdd",
  "Fungi"         = "#e3dcf4",
  "Protists"      = "#e5e5e5",
  "Other"         = "#f0f0f0"
)
pal_bank <- list(
  "Invertebrates" = c("#e66101", "#fdb863", "#b2abd2", "#5e3c99", "#1b9e77",
                      "#d95f02", "#7570b3", "#66a61e", "#e7298a"),
  "Vertebrates"   = c("#1f78b4", "#33a02c", "#6a3d9a", "#ff7f00", "#b15928"),
  "Viridiplantae" = c("#006d2c", "#31a354", "#74c476", "#a1d99b"),
  "Fungi"         = c("#7a0177", "#c51b8a", "#f768a1", "#fa9fb5"),
  "Protists"      = c("#525252", "#969696"),
  "Other"         = c("#8c8c8c")
)

# ── Edge group assignments (pure-descendant branch colouring) ─────────────────
tip_major_map <- setNames(plot_df$broad_group, plot_df$label)
tip_sub_map   <- setNames(plot_df$subclade,    plot_df$label)

edge_group_for_node <- function(node_id) {
  if (node_id <= Ntip(tree_plot)) {
    return(unname(tip_sub_map[tree_plot$tip.label[node_id]]))
  }
  d <- phytools::getDescendants(tree_plot, node_id)
  d <- d[d <= Ntip(tree_plot)]
  if (!length(d)) return("Other")
  tip_labels <- tree_plot$tip.label[d]
  groups <- unique(na.omit(tip_sub_map[tip_labels]))
  majors <- unique(na.omit(tip_major_map[tip_labels]))
  if (!length(groups)) return("Other")
  if (length(groups) == 1) return(groups[[1]])
  if (length(majors) == 1) return(paste0("Mixed_", majors[[1]]))
  "Mixed_cross"
}

edge_df <- tibble(
  node       = tree_plot$edge[, 2],
  edge_group = vapply(tree_plot$edge[, 2], edge_group_for_node, character(1))
)

sub_to_major <- plot_df %>% distinct(subclade, broad_group)
legend_df <- sub_to_major %>%
  group_by(broad_group) %>%
  arrange(subclade, .by_group = TRUE) %>%
  mutate(
    idx          = row_number(),
    branch_label = paste0(broad_group, ": ", subclade),
    branch_col   = {
      bank <- pal_bank[[first(broad_group)]]
      bank[((idx - 1) %% length(bank)) + 1]
    },
    branch_col = ifelse(subclade == "Other invertebrates", "#9e9e9e", branch_col)
  ) %>%
  ungroup()

detail_to_label <- setNames(legend_df$branch_label, legend_df$subclade)
edge_df$edge_label <- ifelse(
  grepl("^Mixed_", edge_df$edge_group) | edge_df$edge_group == "Mixed_cross",
  edge_df$edge_group,
  unname(detail_to_label[edge_df$edge_group])
)

branch_palette <- setNames(legend_df$branch_col, legend_df$branch_label)
branch_breaks  <- legend_df$branch_label
branch_levels  <- legend_df$branch_label
for (mx in c("Mixed_Invertebrates", "Mixed_Vertebrates", "Mixed_Viridiplantae",
             "Mixed_Fungi", "Mixed_Protists", "Mixed_Other", "Mixed_cross")) {
  if (mx %in% edge_df$edge_label) branch_levels <- c(branch_levels, mx)
}
branch_palette <- c(
  branch_palette,
  setNames(rep("#000000", 7),
           c("Mixed_Invertebrates", "Mixed_Vertebrates", "Mixed_Viridiplantae",
             "Mixed_Fungi", "Mixed_Protists", "Mixed_Other", "Mixed_cross"))
)
edge_df$edge_label <- factor(edge_df$edge_label, levels = branch_levels)

# ── Tip symbols dataframe ─────────────────────────────────────────────────────
tip_df <- ggtree::fortify(tree_plot) %>% filter(isTip) %>% select(label, x, y, angle)
tip_df <- left_join(tip_df, plot_df %>% select(label, symbol, colour, broad_group), by = "label")
tip_df$shape_num <- unname(shape_map[tip_df$symbol])
tip_df$tip_col   <- tip_df$colour
tip_df <- tip_df %>% mutate(x = x + max_x * 0.035)

# ── Singleton subclade tip labels ─────────────────────────────────────────────
singleton_sub  <- plot_df %>% filter(subclade %in% inv_force) %>%
  count(subclade) %>% filter(n == 1) %>% pull(subclade)
singleton_tips <- plot_df %>% filter(subclade %in% singleton_sub) %>%
  select(label, subclade) %>%
  left_join(tip_df %>% select(label, x, y, angle), by = "label") %>%
  mutate(hjust = ifelse(angle > 90 & angle < 270, 1.05, -0.05))

bg_legend_df <- tibble(
  broad_group = factor(names(bg_palette), levels = names(bg_palette)),
  x = Inf, y = Inf
)

# ── Build tree ────────────────────────────────────────────────────────────────
p <- ggtree(tree_plot, layout = "circular", linewidth = 0.2, alpha = 0.8)
p$data <- left_join(p$data, edge_df, by = "node")
p <- p +
  geom_tree(aes(color = edge_label), linewidth = 0.38, alpha = 0.9) +
  geom_point(
    data = tip_df,
    aes(x = x, y = y, shape = symbol),
    inherit.aes = FALSE,
    color = tip_df$tip_col, fill = tip_df$tip_col,
    size = 1.6, stroke = 0.4
  ) +
  # per-tip (species) name labels removed for now
  scale_shape_manual(values = shape_map, drop = FALSE, guide = "none") +
  geom_point(
    data = bg_legend_df,
    aes(x = x, y = y, fill = broad_group),
    inherit.aes = FALSE, shape = 22, size = 3, alpha = 0
  ) +
  scale_fill_manual(
    values = bg_palette, name = "Background clades", drop = TRUE,
    guide = guide_legend(order = 1, override.aes = list(alpha = 0.9))
  ) +
  new_scale_fill() +
  geom_point(
    data = data.frame(
      x = rep(0, 4), y = rep(1, 4),
      ctype = factor(c("Satellite", "Transposon", "Mixed (Sat/Trans)", "Holocentric"),
                     levels = c("Satellite", "Transposon", "Mixed (Sat/Trans)", "Holocentric"))
    ),
    aes(x = x, y = y, fill = ctype),
    inherit.aes = FALSE, shape = 21, size = 0, alpha = 0
  ) +
  scale_fill_manual(
    values = c("Satellite" = "#ff006e", "Transposon" = "#3a86ff",
               "Mixed (Sat/Trans)" = "#8338ec", "Holocentric" = "#2d7d32"),
    name = "Centromere type",
    guide = guide_legend(
      order = 3,
      override.aes = list(
        shape  = c(21L, 21L, 21L, 24L), size = 3, alpha = 1,
        color  = c("#ff006e", "#3a86ff", "#8338ec", "#2d7d32"),
        fill   = c("#ff006e", "#3a86ff", "#8338ec", "#2d7d32"),
        stroke = 0.4
      )
    )
  ) +
  scale_color_manual(
    values = branch_palette, drop = TRUE,
    name = "Branch subclades (within background clades)",
    na.translate = FALSE, breaks = branch_breaks,
    guide = "none"                 # branch colours kept on tree, legend hidden (phylogeny only)
  ) +
  theme_tree2() +
  theme(legend.position = "right", plot.margin = margin(10, 10, 10, 10))

# ── Background clade highlights ───────────────────────────────────────────────
major_clades <- plot_df %>%
  filter(!is.na(broad_group) & broad_group != "") %>%
  group_by(broad_group) %>%
  summarize(tips = list(label), .groups = "drop")

for (i in seq_len(nrow(major_clades))) {
  tips <- major_clades$tips[[i]]
  mg   <- major_clades$broad_group[[i]]
  if (length(tips) >= 2 && mg %in% names(bg_palette)) {
    mrca_node <- getMRCA(tree_plot, tips)
    if (!is.null(mrca_node) && !is.na(mrca_node)) {
      p <- p + geom_hilight(node = mrca_node, fill = bg_palette[[mg]],
                            alpha = 0.2, extendto = max_x)
    }
  }
}

# ── Subclade arch labels ──────────────────────────────────────────────────────
subclade_anno <- legend_df %>%
  transmute(subclade, branch_label, branch_col) %>%
  distinct() %>%
  filter(!subclade %in% c("Other"))

# clade-arc labels removed (branch colours + legend already convey subclades)

# ── Concentric time rings ─────────────────────────────────────────────────────
if (length(ring_x) > 0) {
  p <- p +
    geom_vline(xintercept = ring_x, color = "#999999",
               linewidth = 0.2, alpha = 0.5, linetype = "dotted") +
    geom_text(
      data = data.frame(x = ring_x, y = 0, lab = ring_labels),
      aes(x = x, y = y, label = lab),
      inherit.aes = FALSE, size = 2.2, hjust = -0.05, color = "#666666"
    )
}
p <- p + xlim_tree(max_x * 1.30)

# ── Calibration points ────────────────────────────────────────────────────────
# Build calibration table inline (same constraints as calibrate_325sp.R)
tip_fn_norm <- function(x) {
  if (is.na(x)) return(NA_character_)
  hits <- tree_plot$tip.label[
    tree_plot$tip.label == normalize_id(x) |
    tree_plot$tip.label == normalize_id(paste0(x, ".fa")) |
    tree_plot$tip.label == normalize_id(paste0(x, ".fasta"))
  ]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

cal_raw <- read.delim("/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/01_species_tree/over_calib.tsv", stringsAsFactors = FALSE) %>%
  mutate(
    source       = "calibration",
    index        = row_number(),
    age_interval = paste0(age_min, "-", age_max),
    node         = NA_integer_
  )

# Exclude the cross-kingdom Alveolata->Viridiplantae calibration (used only for root age estimation)
cal <- cal_raw %>% filter(label != "Eimeria_tenella_to_Viridiplantae")

# Normalize tip IDs to match tree labels
cal$tip_a_n <- sapply(cal$tip_a, function(x) {
  if (is.na(x)) return(NA_character_)
  normalize_id(x)
})
cal$tip_b_n <- sapply(cal$tip_b, function(x) {
  if (is.na(x)) return(NA_character_)
  normalize_id(x)
})

# Also try with .fa suffix
resolve_tip <- function(x_raw) {
  if (is.na(x_raw)) return(NA_character_)
  xn <- normalize_id(x_raw)
  if (xn %in% tree_plot$tip.label) return(xn)
  xfa  <- normalize_id(paste0(x_raw, ".fa"))
  if (xfa %in% tree_plot$tip.label) return(xfa)
  xfas <- normalize_id(paste0(x_raw, ".fasta"))
  if (xfas %in% tree_plot$tip.label) return(xfas)
  NA_character_
}
cal$tip_a_n <- sapply(cal$tip_a, resolve_tip)
cal$tip_b_n <- sapply(cal$tip_b, resolve_tip)

cal$node_recomputed <- mapply(function(a, b) {
  if (is.na(a) || is.na(b)) return(NA_integer_)
  if (!(a %in% tree_plot$tip.label) || !(b %in% tree_plot$tip.label)) return(NA_integer_)
  getMRCA(tree_plot, c(a, b))
}, cal$tip_a_n, cal$tip_b_n)

cal$node_plot <- cal$node_recomputed
# Root calibration fixed to tree root
cal$node_plot[cal$label == "Eukaryota_root"] <- as.integer(Ntip(tree_plot) + 1L)
# Viridiplantae: MRCA of land plants + algae representative
algae_tips <- plot_df$label[!is.na(plot_df$taxa1) & plot_df$taxa1 == "Algae" & plot_df$label %in% tree_plot$tip.label]
rose_tip   <- resolve_tip("rosCan_S27_v1")
if (!is.null(rose_tip) && !is.na(rose_tip) && length(algae_tips) >= 1) {
  cal$node_plot[cal$label == "Viridiplantae"] <- getMRCA(tree_plot, c(rose_tip, algae_tips[[1]]))
}
# Metazoa: crown from Eumetazoa-vs-sponge split
met_a <- resolve_tip("idNepFlae1.1")
met_b <- resolve_tip("ocSycCili1.1")
if (!is.na(met_a) && !is.na(met_b) && met_a %in% tree_plot$tip.label && met_b %in% tree_plot$tip.label) {
  cal$node_plot[cal$label == "Metazoa"] <- getMRCA(tree_plot, c(met_a, met_b))
}

resolved_n <- sum(!is.na(cal$node_plot))
cat("Calibration points resolved:", resolved_n, "/", nrow(cal), "\n")

tree_xy <- ggtree::fortify(tree_plot) %>% select(node, x, y)
cal_xy  <- dplyr::left_join(cal, tree_xy, by = c("node_plot" = "node"))
cal_xy  <- cal_xy[!is.na(cal_xy$x) & !is.na(cal_xy$y), ]

if (FALSE) {                      # calibration points removed from the figure (per request)
  cal_xy <- cal_xy %>%
    group_by(node_plot) %>%
    mutate(dup_n = n(), dup_i = row_number()) %>%
    ungroup() %>%
    mutate(
      r      = sqrt(x^2 + y^2),
      ux     = ifelse(r > 0, x / r, cos(pi / 4)),
      uy     = ifelse(r > 0, y / r, sin(pi / 4)),
      delta  = (dup_i - (dup_n + 1) / 2) * max_x * 0.04,
      x_plot = x + ux * delta,
      y_plot = y + uy * delta
    )

  p <- p +
    geom_point(
      data = cal_xy,
      aes(x = x_plot, y = y_plot),
      inherit.aes = FALSE, shape = 21, size = 2.25, stroke = 0.35,
      fill = "#ffd54f", color = "black", alpha = 0.95
    ) +
    geom_text(
      data = cal_xy,
      aes(x = x_plot, y = y_plot, label = index),
      inherit.aes = FALSE, size = 1.7, color = "black",
      fontface = "bold", nudge_y = 0.002
    )
}

# ── Outer rings: HOR regimentation + HOR score (Piotr, cen_families.csv) ──────
cf <- read.csv("/home/jg2070/Desktop/DToL_phylogenomics/01_species_tree/data/cen_families.csv", stringsAsFactors = FALSE)
cf$regimentation_score <- suppressWarnings(as.numeric(cf$regimentation_score))
cf$HOR_score           <- suppressWarnings(as.numeric(cf$HOR_score))
cf$count_total         <- suppressWarnings(as.numeric(cf$count_total))
cf$pfx  <- sub("\\d.*$", "", tolower(cf$fasta))
tip_pfx <- sub("\\d.*$", "", tolower(tree_plot$tip.label))
# Per-species aggregation (see AGG at top): dominant family by copy number, or a
# copy-number-weighted genomic mean across all centromeric families.
if (AGG == "ian") {
  agg <- cf %>% group_by(pfx) %>%
    slice_max(count_total, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(pfx, is_holocentric, regi = regimentation_score, hor = HOR_score)
} else {
  agg <- cf %>% group_by(pfx) %>%
    summarise(is_holocentric = first(is_holocentric),
              regi = weighted.mean(regimentation_score, count_total, na.rm = TRUE),
              hor  = weighted.mean(HOR_score,           count_total, na.rm = TRUE),
              .groups = "drop")
}
lk <- function(pp) { h <- tree_plot$tip.label[tip_pfx == pp]; if (length(h) == 1) return(h)
  h <- tree_plot$tip.label[startsWith(tree_plot$tip.label, pp)]; if (length(h) >= 1) return(h[1]); NA_character_ }
agg$tip <- vapply(agg$pfx, lk, character(1))
cat("cen_families species matched to tree:", sum(!is.na(agg$tip)), "/", nrow(agg), "\n")
regi_df <- tibble(label = tree_plot$tip.label) %>%
  left_join(agg %>% filter(!is.na(tip)) %>%
              transmute(label = tip, regi, hor, holo = is_holocentric == "TRUE"),
            by = "label") %>%
  mutate(regi = ifelse(holo %in% TRUE | !is.finite(regi), NA, regi),
         hor  = ifelse(holo %in% TRUE | !is.finite(hor),  NA, hor))

# (all rings are added together below, in Ian's outside->inside order)

# ── Ring data: satellite monomer length + total amount + genome GC + size + chr# ──
# read.csv renames gc% -> gc.; genome.bp / chrs / genome.gc are genome-level (constant
# per species) so aggregation only affects the satellite-family metrics.
sd <- read.csv("/home/jg2070/Desktop/DToL_phylogenomics/01_species_tree/data/dtol.sat.dat.csv",
               stringsAsFactors = FALSE)
sd$pfx <- sub("\\d.*$", "", tolower(sd$fasta))
if (AGG == "ian") {
  sd_ag <- sd %>% group_by(pfx) %>% slice_max(n, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(pfx, sat_len = width, sat_totbp = tot.bp,
              genome_gc = genome.gc, genome_bp = genome.bp, chrs = chrs)
} else {
  sd_ag <- sd %>% group_by(pfx) %>%
    summarise(sat_len   = weighted.mean(width, n, na.rm = TRUE),
              sat_totbp = sum(tot.bp, na.rm = TRUE),          # total centromeric satellite bp
              genome_gc = first(genome.gc), genome_bp = first(genome.bp),
              chrs = first(chrs), .groups = "drop")
}
sd_ag$tip <- vapply(sd_ag$pfx, lk, character(1))
sat_df <- tibble(label = tree_plot$tip.label) %>%
  left_join(sd_ag %>% filter(!is.na(tip)) %>% select(-pfx) %>% rename(label = tip),
            by = "label")

# ── Genome-level rings for ALL species (from the master table, not just satellite spp) ──
# genome size / chromosome number / host GC are genome properties available for every
# species, so these three rings are populated for all 325 tips (satellite metrics stay
# satellite-only above).
gm <- suppressMessages(read_excel(excel_path, sheet = "Sheet1"))
gm$pfx <- sub("\\d.*$", "", tolower(gm$fasta))
for (cc in c("chrs", "genome.bp", "genome.gc")) gm[[cc]] <- suppressWarnings(as.numeric(gm[[cc]]))
gm <- gm %>% group_by(pfx) %>% slice(1) %>% ungroup()
gm$tip <- vapply(gm$pfx, lk, character(1))
gen_df <- tibble(label = tree_plot$tip.label) %>%
  left_join(gm %>% filter(!is.na(tip)) %>%
              transmute(label = tip, genome_gc = genome.gc, genome_bp = genome.bp, chrs = chrs),
            by = "label")
cat("genome-level metrics (master) matched to tree:", sum(!is.na(gen_df$genome_bp)), "/", Ntip(tree_plot), "\n")

# ── Rings, added INSIDE -> OUTSIDE = Ian's outside->inside order reversed ──────
# (centromere architecture is already shown by the tip symbols, so no ring for it)
ro <- 0.085                        # spacing between rings
p <- p +
  # (1 innermost) HOR regimentation
  ggnewscale::new_scale_fill() +
  geom_fruit(data = regi_df, geom = geom_tile, mapping = aes(y = label, fill = regi),
             width = 0.045 * max_x, offset = 0.10, axis.params = list(axis = "none")) +
  # 7 distinct, colour-blind-safe ColorBrewer sequential ramps (one hue family each)
  scale_fill_gradientn(colours = c("#fff5f0","#fcbba1","#fb6a4a","#cb181d","#67000d"),  # Reds
    na.value = "#e0e0e0", limits = c(0, 100), name = sprintf("HOR regimentation\n(%s)", agg_lab)) +
  # (3) HOR score
  ggnewscale::new_scale_fill() +
  geom_fruit(data = regi_df, geom = geom_tile, mapping = aes(y = label, fill = hor),
             width = 0.045 * max_x, offset = ro, axis.params = list(axis = "none")) +
  scale_fill_gradientn(colours = c("#f7fcf5","#c7e9c0","#74c476","#238b45","#00441b"),  # Greens
    na.value = "#e0e0e0", name = sprintf("HOR score\n(%s)", agg_lab)) +
  # (4) Satellite total amount (bp)
  ggnewscale::new_scale_fill() +
  geom_fruit(data = sat_df, geom = geom_tile, mapping = aes(y = label, fill = sat_totbp),
             width = 0.045 * max_x, offset = ro, axis.params = list(axis = "none")) +
  scale_fill_gradientn(colours = c("#f7fbff","#c6dbef","#6baed6","#2171b5","#08306b"),  # Blues
    na.value = "#e0e0e0", trans = "log10", name = "Satellite total\namount (bp)") +
  # (5) Satellite monomer length (bp)
  ggnewscale::new_scale_fill() +
  geom_fruit(data = sat_df, geom = geom_tile, mapping = aes(y = label, fill = sat_len),
             width = 0.045 * max_x, offset = ro, axis.params = list(axis = "none")) +
  scale_fill_gradientn(colours = c("#fcfbfd","#dadaeb","#9e9ac8","#6a51a3","#3f007d"),  # Purples
    na.value = "#e0e0e0", trans = "log10", name = "Satellite\nmonomer length (bp)") +
  # (6) Host genome GC%  (all species, from master)
  ggnewscale::new_scale_fill() +
  geom_fruit(data = gen_df, geom = geom_tile, mapping = aes(y = label, fill = genome_gc),
             width = 0.045 * max_x, offset = ro, axis.params = list(axis = "none")) +
  scale_fill_gradientn(colours = c("#ffffe5","#fee391","#fe9929","#d95f0e","#993404"),  # YlOrBr
    na.value = "#e0e0e0", limits = c(20, 60), name = "Host genome\nGC %") +
  # (7) Chromosome number  (all species, from master)
  ggnewscale::new_scale_fill() +
  geom_fruit(data = gen_df, geom = geom_tile, mapping = aes(y = label, fill = chrs),
             width = 0.045 * max_x, offset = ro, axis.params = list(axis = "none")) +
  scale_fill_gradientn(colours = c("#edf8fb","#b2e2e2","#66c2a4","#2ca25f","#006d2c"),  # BuGn teal (grey reserved for NA)
    na.value = "#e0e0e0", name = "Chromosome\nnumber (n)") +
  # (8 outermost) Genome size (Mb)  (all species, from master)
  ggnewscale::new_scale_fill() +
  geom_fruit(data = gen_df, geom = geom_tile, mapping = aes(y = label, fill = genome_bp / 1e6),
             width = 0.045 * max_x, offset = ro, axis.params = list(axis = "none")) +
  scale_fill_gradientn(colours = c("#f7f4f9","#d4b9da","#df65b0","#ce1256","#67001f"),  # PuRd
    na.value = "#e0e0e0", trans = "log10", name = "Genome size\n(Mb)")

# ── Save outputs ──────────────────────────────────────────────────────────────
pdf(out_pdf, width = 15, height = 18)
print(p)
dev.off()

png(out_png, width = 2000, height = 2400, res = 250)
print(p)
dev.off()

# Copy to top-level figures/
file.copy(out_pdf, file.path(TOP_FIG, basename(out_pdf)), overwrite = TRUE)
file.copy(out_png, file.path(TOP_FIG, basename(out_png)), overwrite = TRUE)

match_n <- sum(!is.na(plot_df$genus) & !is.na(plot_df$species))
cat("Tree tips:", Ntip(tree_calibrated), "\n")
cat("Matched genus/species labels:", match_n, "\n")
cat("Classification counts:\n")
print(sort(table(plot_df$classification), decreasing = TRUE))
cat("Wrote:\n", out_pdf, "\n", out_png, "\n", sep = "")
