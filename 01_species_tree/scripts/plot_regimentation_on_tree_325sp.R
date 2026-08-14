#!/usr/bin/env Rscript
# plot_regimentation_on_tree_325sp.R
# HOR regimentation score (Piotr, data/cen_families.csv) mapped onto the 325-sp
# chronos-correlated tree. Per repeat family (chromosomes with >=10 HORs),
# aggregated per species = MAX regimentation (most-regimented array).
# Holocentrics carry no regimentation (method doesn't translate) -> grey.
# Self-contained: paths are relative to this repo's 01_species_tree/.

suppressPackageStartupMessages({
  library(ape); library(ggtree); library(ggtreeExtra); library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = FALSE)
sf   <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(sf)) normalizePath(file.path(dirname(sf), "..")) else getwd()
TREE_F  <- file.path(root, "trees", "full_325sp_chronos_over_correlated.nwk")
CSV_F   <- file.path(root, "data",  "cen_families.csv")
FIG_DIR <- file.path(root, "figures"); dir.create(FIG_DIR, showWarnings = FALSE)
DATA_DIR<- file.path(root, "data")

# ── tree ──────────────────────────────────────────────────────────────────────
tr <- read.tree(TREE_F)
norm <- function(x) tolower(sub("\\.(fa|fasta)$", "", x))
tr$tip.label <- norm(tr$tip.label)
pfx <- function(x) sub("\\d.*$", "", tolower(x))
tip_pfx <- pfx(tr$tip.label)

# ── data: regimentation per family -> per species (max) ───────────────────────
d <- read.csv(CSV_F, stringsAsFactors = FALSE)
d$regimentation_score <- suppressWarnings(as.numeric(d$regimentation_score))
d$HOR_score           <- suppressWarnings(as.numeric(d$HOR_score))
d$pfx <- pfx(d$fasta)

sp <- d %>%
  group_by(pfx, is_holocentric) %>%
  summarise(regi_max = suppressWarnings(max(regimentation_score, na.rm = TRUE)),
            hor_max  = suppressWarnings(max(HOR_score, na.rm = TRUE)),
            n_fam    = n(), .groups = "drop")

tip_lookup <- function(p) {
  hit <- tr$tip.label[tip_pfx == p]
  if (length(hit) == 1) return(hit)
  hit <- tr$tip.label[startsWith(tr$tip.label, p)]
  if (length(hit) >= 1) return(hit[1])
  NA_character_
}
sp$tip <- vapply(sp$pfx, tip_lookup, character(1))
cat(sprintf("Species matched to tree: %d / %d\n", sum(!is.na(sp$tip)), nrow(sp)))

ring <- data.frame(label = tr$tip.label) %>%
  left_join(sp %>% filter(!is.na(tip)) %>%
              transmute(label = tip, regi = regi_max,
                        holo = is_holocentric == "TRUE"),
            by = "label") %>%
  mutate(regi = ifelse(holo %in% TRUE | !is.finite(regi), NA, regi))

# ── plot: fan tree + regimentation heatmap ring ───────────────────────────────
p <- ggtree(tr, layout = "fan", open.angle = 10, linewidth = 0.22) +
  geom_fruit(data = ring, geom = geom_tile,
             mapping = aes(y = label, fill = regi),
             width = 35, offset = 0.04,
             axis.params = list(axis = "none")) +
  scale_fill_gradientn(
    colours  = c("#fff5eb", "#fdd0a2", "#fd8d3c", "#e6550d", "#a63603", "#4d0000"),
    na.value = "#e0e0e0",
    name     = "HOR regimentation\nscore (max/species)",
    limits   = c(0, 100)) +
  labs(title = "Centromere HOR regimentation across the 325-species DToL tree",
       subtitle = "Per repeat family (>=10 HORs/chr), aggregated per species as max regimentation. Grey = holocentric (not scored) or no data.") +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 9, colour = "grey35"))

for (ext in c("pdf", "png")) {
  ggsave(file.path(FIG_DIR, paste0("regimentation_on_tree_325sp.", ext)),
         p, width = 12, height = 12, dpi = if (ext == "png") 300 else 100, bg = "white")
}

# per-species summary table (max family, non-holocentric)
best <- d %>% filter(is_holocentric != "TRUE", is.finite(regimentation_score)) %>%
  group_by(fasta) %>% slice_max(regimentation_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(fasta, regimentation_max = round(regimentation_score, 3),
            top_family = TRASH_new_class, HOR_score = round(HOR_score, 3)) %>%
  arrange(desc(regimentation_max))
write.table(best, file.path(DATA_DIR, "regimentation_per_species.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved: figures/regimentation_on_tree_325sp.{pdf,png} + data/regimentation_per_species.tsv\n")
