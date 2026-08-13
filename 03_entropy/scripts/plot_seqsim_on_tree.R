#!/usr/bin/env Rscript
# Plot satellite sequence similarity (BLASTN global % identity) mapped onto
# the 325-sp calibrated species tree.
#
# For each internal node: mean similarity between all satellite species pairs
# that span that node (one species from each daughter subtree).
# Only satellite-bearing species contribute; tips without satellites are grey.
#
# Outputs: figures/satellite_halflife/seqsim_on_tree.{pdf,png}

suppressPackageStartupMessages({
  library(ape); library(phangorn); library(ggtree); library(ggplot2)
  library(dplyr); library(tidyr); library(ggnewscale)
})

BASE    <- "/home/jg2070/Desktop/dtol_review_August"
PUB     <- file.path(BASE, "DToL_phylogenomics_publication_325genomes")
FIG_DIR <- file.path(PUB, "03_entropy/figures/satellite_halflife")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

TREE_F  <- file.path(PUB, "01_species_tree/outputs/full_325sp_calibrated.nwk")
# mode: "edlib" (default) or "blastn" — set via Rscript --args blastn
args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args) && args[1] == "blastn") "blastn" else "edlib"

SIM_F   <- file.path(BASE, "2026_trees/annotation_centromeres/repeat_similarity",
                     if (MODE == "blastn") "seqsim_blastn_325sp.tsv"
                     else                 "seqsim_allvsall_325sp.tsv")
SIM_FALL<- file.path(BASE, "2026_trees/annotation_centromeres/repeat_similarity",
                     if (MODE == "blastn") "seqsim_allvsall_325sp.tsv"
                     else                 "seqsim_blastn_325sp.tsv")
SAT_F   <- file.path(BASE, "2026_trees/annotation_centromeres/organized/plots",
                     "centromere_length/genome_vs_centromere_length_satellite_species_table.tsv")
TAX_F   <- file.path(BASE, "2026_trees/annotation_centromeres/centromere_code_to_species.tsv")

# ── tree ──────────────────────────────────────────────────────────────────────
tr <- read.tree(TREE_F)
tr$tip.label <- sub("\\.fa$", "", tr$tip.label)
# species code = lowercase prefix before digits
tip_code <- sub("[0-9.].*$", "", tolower(tr$tip.label))
names(tip_code) <- tr$tip.label   # tip.label → code

# ── similarity data — all-vs-all edlib primary, BLASTN fallback ──────────────
build_lookup <- function(df, col) {
  d <- data.frame(a = pmin(df$spA, df$spB),
                  b = pmax(df$spA, df$spB),
                  v = df[[col]], stringsAsFactors = FALSE)
  setNames(d$v, paste(d$a, d$b, sep = "|"))
}

sim_allvsall <- read.delim(SIM_F,    stringsAsFactors = FALSE)
sim_blastn   <- read.delim(SIM_FALL, stringsAsFactors = FALSE)

lk_avsa   <- build_lookup(sim_allvsall, "mean_pct_id")
lk_blastn <- build_lookup(sim_blastn,   "mean_pct_id")

cat(sprintf("all-vs-all pairs: %d  |  BLASTN fallback pairs: %d\n",
            length(lk_avsa), length(lk_blastn)))

get_sim <- function(cA, cB) {
  key <- paste(sort(c(cA, cB)), collapse = "|")
  v   <- lk_avsa[key];   if (!is.na(v)) return(as.numeric(v))
  v   <- lk_blastn[key]; if (!is.na(v)) return(as.numeric(v))
  NA_real_
}

# ── satellite species set ──────────────────────────────────────────────────────
sat <- read.delim(SAT_F, stringsAsFactors = FALSE)
sat_codes <- sub("[0-9.].*$", "", tolower(sat$species_id[sat$is_satellite_species]))

# ── taxonomy for tip colours ──────────────────────────────────────────────────
tax <- read.delim(TAX_F, stringsAsFactors = FALSE)
tax$code <- sub("[0-9.].*$", "", tolower(tax$fasta_base))
vert  <- c("Actinopterygii","Aves","Mammalia","Reptilia","Amphibia","Chondrichthyes")
virid <- c("Algae","Bryophyta","Dicots","Monocots","Gymnosperms")
tax$broad <- ifelse(tax$taxa1 %in% vert, "Vertebrates",
             ifelse(tax$taxa1 %in% virid, "Viridiplantae",
             ifelse(tax$taxa1 == "Fungi", "Fungi", "Invertebrates")))
code2broad <- setNames(tax$broad, tax$code)

# ── compute per-node mean similarity ─────────────────────────────────────────
# For each internal node: get descendant tips from left vs right subtrees,
# keep only satellite species, compute mean pairwise similarity across subtrees.
cat("Computing per-node satellite similarity …\n")

n_tips  <- Ntip(tr)
n_nodes <- Nnode(tr)

# descendants of each node (tip indices)
desc_tips <- lapply(seq_len(n_tips + n_nodes), function(nd) {
  tips <- unlist(Descendants(tr, nd, type = "tips"))
  tips
})

# for each internal node, compute mean cross-subtree satellite similarity
node_sim   <- rep(NA_real_, n_nodes)
node_npair <- rep(NA_integer_, n_nodes)

for (i in seq_len(n_nodes)) {
  nd       <- n_tips + i
  children <- tr$edge[tr$edge[,1] == nd, 2]
  if (length(children) < 2) next

  left_idx  <- intersect(unlist(desc_tips[[children[1]]]), seq_len(n_tips))
  right_idx <- intersect(unlist(desc_tips[[children[2]]]), seq_len(n_tips))

  left_codes  <- unique(tip_code[left_idx])
  right_codes <- unique(tip_code[right_idx])
  left_codes  <- left_codes[left_codes %in% sat_codes]
  right_codes <- right_codes[right_codes %in% sat_codes]

  if (length(left_codes) == 0 || length(right_codes) == 0) next

  vals <- na.omit(as.vector(outer(left_codes, right_codes,
                                   Vectorize(get_sim))))
  if (length(vals) == 0) next
  node_sim[i]   <- mean(vals)
  node_npair[i] <- length(vals)
}

cat(sprintf("  Nodes with similarity data: %d / %d\n",
            sum(!is.na(node_sim)), n_nodes))
cat(sprintf("  Range: %.1f – %.1f%%  |  pairs per node: %d – %d\n",
            min(node_sim, na.rm=TRUE), max(node_sim, na.rm=TRUE),
            min(node_npair, na.rm=TRUE), max(node_npair, na.rm=TRUE)))
cat(sprintf("  Uncoloured: nodes where one/both subtrees lack satellite species\n"))

# ── tip metadata ───────────────────────────────────────────────────────────────
tip_df <- data.frame(
  label  = tr$tip.label,
  code   = tip_code,
  is_sat = tip_code %in% sat_codes,
  broad  = code2broad[tip_code],
  stringsAsFactors = FALSE
)

broad_pal <- c(Vertebrates  = "#1565C0", Invertebrates = "#E65100",
               Viridiplantae= "#2E7D32", Fungi         = "#6A1B9A")

# node data frame
node_df <- data.frame(
  node    = seq(n_tips + 1, n_tips + n_nodes),
  sat_sim = node_sim,
  n_pairs = node_npair
)

# ── plot ───────────────────────────────────────────────────────────────────────
# ── clade MRCA nodes ──────────────────────────────────────────────────────────
get_mrca <- function(taxa_vec) {
  codes <- tax$code[tax$taxa1 %in% taxa_vec]
  tips  <- tr$tip.label[tip_code %in% codes]
  if (length(tips) < 2) return(NA_integer_)
  getMRCA(tr, tips)
}

clade_labels <- list(
  "Vertebrates"   = get_mrca(c("Actinopterygii","Aves","Mammalia","Reptilia","Amphibia","Chondrichthyes")),
  "Aves"          = get_mrca("Aves"),
  "Mammalia"      = get_mrca("Mammalia"),
  "Actinopterygii"= get_mrca("Actinopterygii"),
  "Insects"       = get_mrca(c("Coleoptera","Diptera","Hymenoptera","Lepidoptera",
                                "Hemiptera","Ephemeroptera","Odonata","Trichoptera",
                                "Neuroptera","Dermaptera","Plecoptera","Psocodea","Blattodea")),
  "Viridiplantae" = get_mrca(c("Bryophyta","Dicots","Monocots","Algae","Gymnosperms")),
  "Monocots"      = get_mrca("Monocots"),
  "Dicots"        = get_mrca("Dicots"),
  "Bryophyta"     = get_mrca("Bryophyta"),
  "Fungi"         = get_mrca("Fungi")
)
clade_df <- data.frame(
  node  = unlist(clade_labels),
  label = names(clade_labels),
  stringsAsFactors = FALSE
) %>% filter(!is.na(node))
cat("Clade labels:\n"); print(clade_df)

# ── merge tip + node data ──────────────────────────────────────────────────────
all_df <- bind_rows(
  tip_df  %>% mutate(node = match(label, tr$tip.label),
                     sat_sim = NA_real_, n_pairs = NA_integer_),
  node_df %>% mutate(label = NA_character_, code = NA_character_,
                     is_sat = NA, broad = NA_character_)
)

# ── plot ───────────────────────────────────────────────────────────────────────
p <- ggtree(tr, layout = "fan", linewidth = 0.12, color = "grey70") %<+% all_df +
  # satellite tip points
  geom_tippoint(aes(color = broad, subset = !is.na(is_sat) & is_sat),
                size = 1.0, alpha = 0.85) +
  geom_tiplab(aes(subset = !is.na(is_sat) & is_sat,
                  label = code),
              size = 1.4, offset = 1.5, align = FALSE) +
  scale_color_manual(values = broad_pal, name = "Group (satellite spp.)",
                     na.value = "#cccccc") +
  ggnewscale::new_scale_color() +
  # internal nodes: colour = similarity, size = number of species pairs
  geom_nodepoint(aes(color = sat_sim,
                     size  = n_pairs,
                     subset = !is.na(sat_sim)),
                 alpha = 0.90) +
  scale_size_continuous(name = "Species pairs\nat node",
                        range = c(0.8, 7),
                        breaks = c(1, 5, 20, 50),
                        guide  = guide_legend(override.aes = list(alpha=1,
                                                                   color="#fc8d59"))) +
  scale_color_gradientn(
    colors   = c("#ffffcc","#fed976","#fd8d3c","#e31a1c","#800026"),
    limits   = c(floor(min(node_sim, na.rm=TRUE) / 5) * 5,
                 ceiling(max(node_sim, na.rm=TRUE) / 5) * 5),
    name     = paste0("Mean satellite\nsimilarity (%)\n[", MODE, "]"),
    na.value = NA,
    guide    = guide_colorbar(title.position = "top", barwidth = 0.8, barheight = 5)
  ) +
  # clade labels
  geom_cladelab(
    data     = clade_df,
    mapping  = aes(node = node, label = label),
    fontsize = 2.8, offset = 3, offset.text = 1.5,
    barsize  = 0.4, barcolor = "grey50", textcolor = "grey20",
    align    = TRUE, hjust = 0.5
  ) +
  labs(
    title    = "Satellite sequence similarity mapped onto the 325-sp tree",
    subtitle = paste0(
      "Node colour: mean BLASTN global % identity between satellite species\n",
      "in daughter subtrees (all-vs-all edlib NW; fitting: A·exp(-lt) + C).\n",
      "Node size: number of cross-subtree species pairs compared.\n",
      "Uncoloured nodes: one or both subtrees contain no satellite species.")
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 7, hjust = 0.5, color = "grey40"),
    legend.position = "right"
  )

out_stem <- file.path(FIG_DIR, paste0("seqsim_on_tree_", MODE))
ggsave(paste0(out_stem, ".pdf"), p, width = 14, height = 14)
ggsave(paste0(out_stem, ".png"), p, width = 14, height = 14, dpi = 300)
cat(sprintf("Saved: seqsim_on_tree_%s.{pdf,png}\n", MODE))
