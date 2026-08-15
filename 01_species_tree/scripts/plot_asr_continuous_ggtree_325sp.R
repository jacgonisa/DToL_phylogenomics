#!/usr/bin/env Rscript
# plot_asr_continuous_ggtree_325sp.R
# Publication-style ancestral-state reconstruction of HOR regimentation and HOR score
# on the calibrated 325-sp chronos tree, rendered with ggtree (smooth gradient along
# branches via continuous="colour"; ancestral node states = phytools::fastAnc, ML/BM).
# Annotated with phylogenetic-signal (Blomberg K) and OU-vs-BM AICc support.

suppressPackageStartupMessages({
  library(ape); library(phytools); library(ggtree); library(ggplot2)
  library(dplyr); library(stringr); library(patchwork)
})

BASE    <- "/home/jg2070/Desktop/dtol_review_August"
PUB_DIR <- file.path(BASE, "DToL_phylogenomics_publication_325genomes")
TREE_F  <- file.path(PUB_DIR, "01_species_tree/outputs/full_325sp_chronos_over_correlated_fa.nwk")
CF_F    <- "/home/jg2070/Downloads/cen_families.csv"
FIG_DIR <- file.path(PUB_DIR, "01_species_tree/figures")
SIG_F   <- file.path(FIG_DIR, "asr_continuous_phylosig.tsv")          # from asr_continuous_*.R
MOD_F   <- file.path(FIG_DIR, "asr_continuous_model_selection.tsv")

strip_ext    <- function(x) str_replace(x, "\\.(fa|fasta)$", "")
normalize_id <- function(x) x |> str_replace("\\.fasta\\.F2B\\.ChrOnly\\.fa$",".fa") |> strip_ext() |> tolower()
label_fix_map <- c("daTanVulg1.hap2.1.fasta" = "daTanVulg1.hap1.1.fa")

tr <- read.tree(TREE_F)
fixed <- ifelse(tr$tip.label %in% names(label_fix_map), unname(label_fix_map[tr$tip.label]), tr$tip.label)
tr$tip.label <- normalize_id(fixed)
tip_pfx <- sub("\\d.*$", "", tolower(tr$tip.label))

cf <- read.csv(CF_F, stringsAsFactors = FALSE)
cf$regimentation_score <- suppressWarnings(as.numeric(cf$regimentation_score))
cf$HOR_score           <- suppressWarnings(as.numeric(cf$HOR_score))
cf$count_total         <- suppressWarnings(as.numeric(cf$count_total))
cf$pfx <- sub("\\d.*$", "", tolower(cf$fasta))
agg <- cf %>% group_by(pfx) %>% slice_max(count_total, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(pfx, is_holo = is_holocentric == "TRUE", regi = regimentation_score, hor = HOR_score)
lk <- function(pp){h<-tr$tip.label[tip_pfx==pp];if(length(h)==1)return(h)
  h<-tr$tip.label[startsWith(tr$tip.label,pp)];if(length(h)>=1)return(h[1]);NA_character_}
agg$tip <- vapply(agg$pfx, lk, character(1))
agg <- agg %>% filter(!is.na(tip), !is_holo)

# optional stat annotations (may be absent if the analysis script wasn't run first)
sig <- tryCatch(read.delim(SIG_F), error = function(e) NULL)
mod <- tryCatch(read.delim(MOD_F), error = function(e) NULL)
stat_sub <- function(nm) {
  s <- if (!is.null(sig)) sig[sig$trait == nm, ] else NULL
  m <- if (!is.null(mod)) mod[mod$trait == nm, ] else NULL
  bits <- c()
  if (!is.null(s) && nrow(s)) bits <- c(bits, sprintf("Blomberg K = %.2f", s$K[1]),
                                              sprintf("Pagel λ = %.2f", s$lambda[1]))
  if (!is.null(m) && nrow(m)) { best <- m[order(m$AICc), ][1, ]
    bits <- c(bits, sprintf("best model: %s (w=%.2f, BM ΔAICc=+%.0f)",
                            best$model, best$weight, m$dAICc[m$model=="BM"][1])) }
  paste(bits, collapse = "   ·   ")
}

vec_of <- function(col){v<-agg[[col]][match(tr$tip.label,agg$tip)];names(v)<-tr$tip.label;v[is.finite(v)]}

one_panel <- function(x, cols, title, nm) {
  drop <- setdiff(tr$tip.label, names(x))
  trx  <- if (length(drop)) drop.tip(tr, drop) else tr
  x    <- x[trx$tip.label]
  anc  <- fastAnc(trx, x)                                  # ML/BM ancestral node states
  # node-indexed state vector: tips 1..Ntip, internal nodes (Ntip+1)..
  st <- numeric(Ntip(trx) + trx$Nnode)
  st[1:Ntip(trx)] <- x[trx$tip.label]
  st[as.integer(names(anc))] <- anc
  dd <- data.frame(node = seq_along(st), state = st)

  p <- ggtree(trx, layout = "circular", aes(color = state),
              continuous = "colour", size = 0.55) %<+% dd +
    scale_color_gradientn(colours = cols, name = title,
                          guide = guide_colourbar(barheight = 8)) +
    labs(title = sprintf("Ancestral %s", title), subtitle = stat_sub(nm)) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 8, colour = "grey35"),
          plot.margin = margin(6, 6, 6, 6))
  list(p = p, trx = trx)
}

pr <- one_panel(vec_of("regi"), c("#fff5eb","#fdd0a2","#fd8d3c","#e6550d","#a63603","#4d0000"),
                "HOR regimentation", "regimentation")
ph <- one_panel(vec_of("hor"),  c("#f7fcf5","#c7e9c0","#74c476","#238b45","#00441b"),
                "HOR score", "HORscore")

# individual + combined outputs
save_fig <- function(plot, stem, w, h) {
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), plot, width = w, height = h)
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), plot, width = w, height = h, dpi = 300, bg = "white")
  cat("Wrote:", file.path(FIG_DIR, paste0(stem, ".{pdf,png}")), "\n")
}
save_fig(pr$p, "asr_ggtree_regimentation_325sp", 9, 9)
save_fig(ph$p, "asr_ggtree_HORscore_325sp",      9, 9)
save_fig(pr$p + ph$p + patchwork::plot_layout(ncol = 2),
         "asr_ggtree_combined_325sp", 17, 9)
