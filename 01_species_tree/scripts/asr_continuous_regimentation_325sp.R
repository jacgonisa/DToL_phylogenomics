#!/usr/bin/env Rscript
# asr_continuous_regimentation_325sp.R
# Continuous-trait ancestral-state reconstruction (phytools contMap, ML/fastAnc under BM)
# of HOR regimentation and HOR score on the calibrated 325-sp chronos tree, plus
# phylogenetic-signal statistics (Blomberg's K, Pagel's lambda).
#
# Motivation: the PI observes regimentation/HOR "flickering" up and down across the
# phylogeny with no clear phylogenetic signal (consistent with the alpha-omega cyclical
# centromere-evolution hypothesis). contMap visualises the ancestral trajectory and
# phylosig quantifies the (expected weak) signal.
#
# ponytail: contMap's ML anc estimation assumes Brownian motion; low K/lambda is the
# result we want to report, not a modelling failure -> report plainly, no model-shopping.

suppressPackageStartupMessages({
  library(ape); library(phytools); library(dplyr); library(stringr)
})

BASE     <- "/home/jg2070/Desktop/dtol_review_August"
PUB_DIR  <- file.path(BASE, "DToL_phylogenomics_publication_325genomes")
TREE_F   <- file.path(PUB_DIR, "01_species_tree/outputs/full_325sp_chronos_over_correlated_fa.nwk")
CF_F     <- "/home/jg2070/Downloads/cen_families.csv"
FIG_DIR  <- file.path(PUB_DIR, "01_species_tree/figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── tree + tip normalisation (same logic as plot_tree_regimentation_chronos.R) ──
strip_ext    <- function(x) str_replace(x, "\\.(fa|fasta)$", "")
normalize_id <- function(x) {
  x |> str_replace("\\.fasta\\.F2B\\.ChrOnly\\.fa$", ".fa") |> strip_ext() |> tolower()
}
label_fix_map <- c("daTanVulg1.hap2.1.fasta" = "daTanVulg1.hap1.1.fa")

tr <- read.tree(TREE_F)
fixed <- ifelse(tr$tip.label %in% names(label_fix_map),
                unname(label_fix_map[tr$tip.label]), tr$tip.label)
tr$tip.label <- normalize_id(fixed)
tip_pfx <- sub("\\d.*$", "", tolower(tr$tip.label))

# ── dominant-array regimentation / HOR per species (matches the ring aggregation) ──
cf <- read.csv(CF_F, stringsAsFactors = FALSE)
cf$regimentation_score <- suppressWarnings(as.numeric(cf$regimentation_score))
cf$HOR_score           <- suppressWarnings(as.numeric(cf$HOR_score))
cf$count_total         <- suppressWarnings(as.numeric(cf$count_total))
cf$pfx <- sub("\\d.*$", "", tolower(cf$fasta))
agg <- cf %>% group_by(pfx) %>%
  slice_max(count_total, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(pfx, is_holo = is_holocentric == "TRUE",
            regi = regimentation_score, hor = HOR_score)

lk <- function(pp) { h <- tr$tip.label[tip_pfx == pp]; if (length(h) == 1) return(h)
  h <- tr$tip.label[startsWith(tr$tip.label, pp)]; if (length(h) >= 1) return(h[1]); NA_character_ }
agg$tip <- vapply(agg$pfx, lk, character(1))
agg <- agg %>% filter(!is.na(tip), !is_holo)   # holocentrics not scored -> excluded

# named numeric vectors keyed by tip label
vec_of <- function(col) { v <- agg[[col]][match(tr$tip.label, agg$tip)]
  names(v) <- tr$tip.label; v[is.finite(v)] }

traits <- list(
  regimentation = list(x = vec_of("regi"),
                       cols = c("#fff5eb","#fdd0a2","#fd8d3c","#e6550d","#a63603","#4d0000"),
                       title = "HOR regimentation"),
  HORscore      = list(x = vec_of("hor"),
                       cols = c("#f7fcf5","#c7e9c0","#74c476","#238b45","#00441b"),
                       title = "HOR score")
)

# ── contMap + phylosig per trait ────────────────────────────────────────────────
sig <- list()
for (nm in names(traits)) {
  tt   <- traits[[nm]]
  x    <- tt$x
  drop <- setdiff(tr$tip.label, names(x))
  trx  <- if (length(drop)) drop.tip(tr, drop) else tr
  x    <- x[trx$tip.label]
  cat(sprintf("[%s] tips with data: %d\n", nm, length(x)))

  cm <- contMap(trx, x, plot = FALSE, res = 400)
  cm <- setMap(cm, colors = tt$cols)

  for (ext in c("pdf", "png")) {
    fn <- file.path(FIG_DIR, sprintf("asr_contMap_%s_325sp.%s", nm, ext))
    if (ext == "pdf") pdf(fn, width = 11, height = 11)
    else png(fn, width = 2200, height = 2200, res = 250)
    plot(cm, type = "fan", legend = 0.7 * max(nodeHeights(trx)),
         fsize = c(0.3, 0.8), lwd = 1.6, outline = FALSE,
         ftype = "off")   # ftype off: 325 tips -> per-tip names unreadable
    title(main = sprintf("Ancestral %s (contMap, ML/BM)", tt$title))
    dev.off()
    cat("Wrote:", fn, "\n")
  }

  kK <- phylosig(trx, x, method = "K",      test = TRUE, nsim = 999)
  kL <- phylosig(trx, x, method = "lambda", test = TRUE)
  sig[[nm]] <- data.frame(
    trait      = nm,
    n_tips     = length(x),
    K          = round(unname(kK$K), 4),
    K.p        = signif(kK$P, 3),
    lambda     = round(unname(kL$lambda), 4),
    lambda.p   = signif(kL$P, 3),
    stringsAsFactors = FALSE)
  cat(sprintf("[%s] Blomberg K=%.3f (p=%.3g)  Pagel lambda=%.3f (p=%.3g)\n",
              nm, kK$K, kK$P, kL$lambda, kL$P))
}

sig_tbl <- do.call(rbind, sig)
out_tsv <- file.path(FIG_DIR, "asr_continuous_phylosig.tsv")
write.table(sig_tbl, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nPhylogenetic-signal summary:\n"); print(sig_tbl)
cat("Wrote:", out_tsv, "\n")
