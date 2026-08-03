#!/usr/bin/env Rscript
# 42_plot_cycles_on_tree.R
#
# Plot the independent evolutionary cycles on the tree:
#   - Tips coloured by observed state
#   - Internal nodes: MAP ancestral state (pie charts from simmap)
#   - Cycle entry nodes highlighted (numbered stars)
#   - Rectangular + fan layouts for metazoa and viridiplantae

suppressPackageStartupMessages({
  library(ape); library(phytools); library(dplyr); library(tidyr)
  library(ggplot2); library(ggtree); library(readxl); library(stringr)
  library(patchwork)
})
set.seed(42)

root      <- normalizePath(file.path(getwd(), "ASR_March2026_327species"))
cache_dir <- file.path(root, "outputs", "cyclical_ARD_irrevH", "cache")
out_dir   <- file.path(root, "outputs", "cycles_on_tree")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

state_levels <- c("H","Mixed","Sat","Trans")
state_colors <- c(H="#f46d43", Mixed="#8e44ad", Sat="#ff2d87", Trans="#4f84e8")
state_labels <- c(H="Holocentric", Mixed="Sat/Transposon", Sat="Satellite", Trans="Transposon")
datasets       <- c("metazoa","viridiplantae")
dataset_labels <- c(metazoa="Metazoa", viridiplantae="Viridiplantae")

# ── Annotations ──────────────────────────────────────────────────────────────
anno <- read.csv(file.path(root,"inputs","full","branch_symbol_anno.tsv"),
                 sep="\t", header=TRUE, stringsAsFactors=FALSE)
anno$tree_label <- gsub(" ","",anno$Species)
map_arch <- c("Holocentric"="H","Satellite/transposon"="Mixed",
              "Satellite"="Sat","Transposon"="Trans")
label_fix_map <- c("daTanVulg1.hap2.1.fasta"="daTanVulg1.hap1.1.fa",
                   "rosCan_S27_v1.fasta.F2B.ChrOnly.fa"="rosCan_S27_v1")
meta <- suppressMessages(
  read_excel(file.path(dirname(root),"DTOL_327_master_March.xlsx"), sheet="Sheet1")
) %>%
  select(fasta, genus, species, taxa1) %>%
  mutate(fasta_fixed = ifelse(fasta %in% names(label_fix_map),
                              unname(label_fix_map[fasta]), fasta),
         fasta_base  = str_replace(fasta_fixed,"\\.(fa|fasta)$",""),
         sp_name     = paste(genus, species)) %>%
  distinct(fasta_base, .keep_all=TRUE)
sp_name_map <- setNames(meta$sp_name,  meta$fasta_base)
taxa1_map   <- setNames(meta$taxa1,    meta$fasta_base)

# ── Helpers ───────────────────────────────────────────────────────────────────
load_dataset <- function(ds) {
  tr <- read.tree(file.path(root,"inputs",ds,"tree_renamed.nw"))
  tr <- midpoint.root(tr); tr <- multi2di(tr, random=FALSE)
  ch <- map_arch[anno$architecture[match(tr$tip.label,anno$tree_label)]]
  names(ch) <- tr$tip.label
  ch <- factor(ch, levels=state_levels); kp <- !is.na(ch)
  tr <- drop.tip(tr, tr$tip.label[!kp])
  ch <- factor(as.character(ch[kp]), levels=state_levels); names(ch) <- tr$tip.label
  list(tree=tr, char=ch, n=length(ch))
}

build_parent <- function(tree) {
  n <- length(tree$tip.label) + tree$Nnode
  p <- setNames(rep(NA_integer_,n), seq_len(n))
  for (i in seq_len(nrow(tree$edge))) p[tree$edge[i,2]] <- tree$edge[i,1]
  p
}

find_entry_nodes <- function(tree, all_states, tip_state, mid_state) {
  n_tips <- length(tree$tip.label)
  parent <- build_parent(tree)
  entry_node_map <- list()
  for (ti in seq_len(n_tips)) {
    if (as.character(all_states[ti]) != tip_state) next
    path <- integer(0); curr <- ti
    while (!is.na(curr)) { path <- c(curr, path); curr <- parent[curr] }
    st  <- as.character(all_states[path])
    rle_res <- rle(st); rv <- rle_res$values
    if (length(rv) < 3) next
    found <- any(rv[seq_len(length(rv)-2)]==tip_state &
                 rv[seq_len(length(rv)-2)+1]==mid_state &
                 rv[seq_len(length(rv)-2)+2]==tip_state)
    if (!found) next
    first_mid <- which(st==mid_state)[1]
    if (!is.na(first_mid)) entry_node_map[[tree$tip.label[ti]]] <- path[first_mid]
  }
  entry_node_map
}

# ── Main loop ─────────────────────────────────────────────────────────────────
for (ds in datasets) {
  cat("\n====", dataset_labels[ds], "====\n")
  dat    <- load_dataset(ds)
  cf     <- file.path(cache_dir, paste0("simmap_ARD_irrevH_",ds,".rds"))
  if (!file.exists(cf)) { cat("No simmap cache\n"); next }
  sm_list <- readRDS(cf)
  sm_sum  <- summary(sm_list)
  ace_mat <- sm_sum$ace
  n_tips  <- length(dat$tree$tip.label)

  # Handle ace_mat dimensions
  if (nrow(ace_mat) == n_tips + dat$tree$Nnode) {
    ace_internal <- ace_mat[(n_tips+1):nrow(ace_mat),,drop=FALSE]
  } else {
    ace_internal <- ace_mat
  }

  tip_states  <- as.character(dat$char)
  node_states <- apply(ace_internal, 1, function(r) colnames(ace_internal)[which.max(r)])
  node_probs  <- apply(ace_internal, 1, max)
  all_states  <- factor(c(tip_states, node_states), levels=state_levels)

  # Find independent cycles
  sts_map <- find_entry_nodes(dat$tree, all_states, "Sat",   "Trans")
  tst_map <- find_entry_nodes(dat$tree, all_states, "Trans", "Sat")

  sts_entry_nodes <- unique(unlist(sts_map))
  tst_entry_nodes <- unique(unlist(tst_map))

  cat(sprintf("  STS: %d tips, %d independent cycles\n",
              length(sts_map), length(sts_entry_nodes)))
  cat(sprintf("  TST: %d tips, %d independent cycles\n",
              length(tst_map), length(tst_entry_nodes)))

  # ── Build tip annotation ──────────────────────────────────────────────────
  tip_ann <- data.frame(
    label     = dat$tree$tip.label,
    obs_state = as.character(dat$char),
    sp_name   = coalesce(sp_name_map[dat$tree$tip.label], dat$tree$tip.label),
    taxa1     = coalesce(taxa1_map[dat$tree$tip.label], ""),
    is_sts    = dat$tree$tip.label %in% names(sts_map),
    is_tst    = dat$tree$tip.label %in% names(tst_map),
    stringsAsFactors=FALSE
  )

  # ── Build node annotation ─────────────────────────────────────────────────
  # All internal nodes: MAP state + whether they are cycle entry nodes
  node_ids <- (n_tips+1):(n_tips+dat$tree$Nnode)

  # Number cycles: STS cycles 1..N, TST cycles N+1..M
  sts_cycle_num <- setNames(seq_along(sts_entry_nodes), as.character(sts_entry_nodes))
  tst_cycle_num <- setNames(seq_along(tst_entry_nodes) + length(sts_entry_nodes),
                             as.character(tst_entry_nodes))

  node_ann <- data.frame(
    node      = node_ids,
    map_state = node_states,
    map_prob  = node_probs,
    stringsAsFactors=FALSE
  ) %>%
  mutate(
    is_sts_entry = node %in% sts_entry_nodes,
    is_tst_entry = node %in% tst_entry_nodes,
    cycle_num    = case_when(
      is_sts_entry ~ sts_cycle_num[as.character(node)],
      is_tst_entry ~ tst_cycle_num[as.character(node)],
      TRUE ~ NA_integer_
    ),
    cycle_label = ifelse(!is.na(cycle_num), as.character(cycle_num), NA_character_),
    cycle_type  = case_when(
      is_sts_entry ~ "Sat→Trans→Sat",
      is_tst_entry ~ "Trans→Sat→Trans",
      TRUE ~ NA_character_
    )
  )

  # ── Get tree layout coordinates ───────────────────────────────────────────
  tr_rect <- ggtree(dat$tree, size=0.25, color="grey50")
  layout_rect <- tr_rect$data   # has: node, x, y, isTip, label, ...

  # Join node_ann with layout coordinates
  node_ann_xy <- node_ann %>%
    left_join(layout_rect %>% select(node, x, y), by="node")

  # ── Rectangular tree ──────────────────────────────────────────────────────
  p_rect <- tr_rect %<+% tip_ann +
    geom_tippoint(aes(color=obs_state), size=1.8, alpha=0.9) +
    # MAP state dots at internal nodes (confident ones)
    geom_point(data = node_ann_xy %>% filter(map_prob > 0.7),
               aes(x=x, y=y, color=map_state),
               size=0.8, alpha=0.4, inherit.aes=FALSE) +
    # Cycle entry nodes: numbered filled circles
    geom_point(data = node_ann_xy %>% filter(!is.na(cycle_num)),
               aes(x=x, y=y, fill=cycle_type),
               shape=21, size=6, color="white", stroke=1.5,
               inherit.aes=FALSE) +
    geom_text(data = node_ann_xy %>% filter(!is.na(cycle_num)),
              aes(x=x, y=y, label=cycle_label),
              size=2.8, fontface="bold", color="white", inherit.aes=FALSE) +
    scale_color_manual(values=state_colors, name="State",
                       labels=state_labels, breaks=names(state_labels)) +
    scale_fill_manual(values=c("Sat→Trans→Sat"="#d73027",
                                "Trans→Sat→Trans"="#313695"),
                      name="Cycle type", na.value=NA) +
    labs(title=paste(dataset_labels[ds],"— independent evolutionary cycles on tree"),
         subtitle=paste0("Numbered circles = cycle entry nodes  ",
                         "(red=Sat→Trans→Sat, blue=Trans→Sat→Trans)\n",
                         "Small dots = MAP ancestral state (posterior > 0.7)")) +
    theme_tree2() +
    theme(plot.title=element_text(face="bold", size=12),
          legend.position="right")

  h_rect <- max(8, n_tips * 0.06 + 2)
  ggsave(file.path(out_dir, paste0(ds,"_cycles_rectangular.pdf")),
         p_rect, width=14, height=h_rect)
  ggsave(file.path(out_dir, paste0(ds,"_cycles_rectangular.png")),
         p_rect, width=14, height=h_rect, dpi=200)
  cat("  Saved rectangular tree\n")

  # ── Fan/circular tree ─────────────────────────────────────────────────────
  tr_fan <- ggtree(dat$tree, layout="fan", open.angle=15, size=0.25, color="grey50")
  layout_fan <- tr_fan$data
  node_ann_fan <- node_ann %>%
    left_join(layout_fan %>% select(node, x, y), by="node")

  p_fan <- tr_fan %<+% tip_ann +
    geom_tippoint(aes(color=obs_state), size=1.5, alpha=0.9) +
    geom_tiplab(aes(label=sp_name), size=1.6, offset=0.5, hjust=0) +
    geom_point(data = node_ann_fan %>% filter(!is.na(cycle_num)),
               aes(x=x, y=y, fill=cycle_type),
               shape=21, size=6, color="white", stroke=1.5,
               inherit.aes=FALSE) +
    geom_text(data = node_ann_fan %>% filter(!is.na(cycle_num)),
              aes(x=x, y=y, label=cycle_label),
              size=2.8, fontface="bold", color="white", inherit.aes=FALSE) +
    scale_color_manual(values=state_colors, name="State",
                       labels=state_labels, breaks=names(state_labels)) +
    scale_fill_manual(values=c("Sat→Trans→Sat"="#d73027",
                                "Trans→Sat→Trans"="#313695"),
                      name="Cycle type", na.value=NA) +
    labs(title=paste(dataset_labels[ds],"— cycles (fan layout)")) +
    theme_tree2() +
    theme(plot.title=element_text(face="bold", size=11),
          legend.position="bottom")

  ggsave(file.path(out_dir, paste0(ds,"_cycles_fan.pdf")),
         p_fan, width=16, height=16)
  ggsave(file.path(out_dir, paste0(ds,"_cycles_fan.png")),
         p_fan, width=16, height=16, dpi=200)
  cat("  Saved fan tree\n")

  # ── Cycle legend table ────────────────────────────────────────────────────
  cat("\n  Cycle key:\n")
  for (i in seq_along(sts_entry_nodes)) {
    en   <- sts_entry_nodes[i]
    tips <- names(sts_map)[unlist(sts_map) == en]
    taxa <- unique(coalesce(taxa1_map[tips], tips))
    cat(sprintf("    Cycle %d (STS): %s [%d tips]\n", i,
                paste(taxa, collapse=", "), length(tips)))
  }
  for (i in seq_along(tst_entry_nodes)) {
    en   <- tst_entry_nodes[i]
    tips <- names(tst_map)[unlist(tst_map) == en]
    taxa <- unique(coalesce(taxa1_map[tips], tips))
    cat(sprintf("    Cycle %d (TST): %s [%d tips]\n",
                i + length(sts_entry_nodes),
                paste(taxa, collapse=", "), length(tips)))
  }
}

cat("\nOutputs in:", out_dir, "\n")
