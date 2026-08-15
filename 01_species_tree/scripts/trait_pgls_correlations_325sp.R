#!/usr/bin/env Rscript
# trait_pgls_correlations_325sp.R
# Build a per-species centromere/genome trait table (dominant-array = Ian aggregation)
# and compute pairwise correlations BOTH naively (Pearson) and phylogenetically
# corrected (PGLS via nlme::gls + Pagel's lambda; phylo r from Felsenstein PICs).
# Outputs feed the seaborn pairplot (trait_pairplot_seaborn_325sp.py).

suppressPackageStartupMessages({
  library(ape); library(nlme); library(readxl); library(dplyr); library(stringr)
})

BASE    <- "/home/jg2070/Desktop/dtol_review_August"
PUB     <- file.path(BASE, "DToL_phylogenomics_publication_325genomes")
TREE_F  <- file.path(PUB, "01_species_tree/outputs/full_325sp_chronos_over_correlated_fa.nwk")
CF_F    <- "/home/jg2070/Downloads/cen_families.csv"
SD_F    <- "/home/jg2070/Desktop/DToL_phylogenomics/01_species_tree/data/dtol.sat.dat.csv"
XLSX    <- file.path(BASE, "2026_trees/annotation_centromeres/DTOL_327_master_March.xlsx")
OUT_DIR <- "/home/jg2070/Desktop/DToL_phylogenomics/01_species_tree/data"
# Aggregation for poly-/di-typic species: ian = dominant family (highest copy number)
#                                         piotr = copy-number-weighted genomic mean
AGG     <- tolower(Sys.getenv("AGG", "ian")); stopifnot(AGG %in% c("ian","piotr"))
TAB_F   <- file.path(OUT_DIR, sprintf("centromere_trait_table_325sp_%s.tsv", AGG))
COR_F   <- file.path(OUT_DIR, sprintf("pgls_trait_correlations_325sp_%s.tsv", AGG))
PART_F  <- file.path(OUT_DIR, sprintf("pgls_partial_regi_hor_monomer_%s.tsv", AGG))

strip_ext    <- function(x) str_replace(x, "\\.(fa|fasta)$", "")
normalize_id <- function(x) x |> str_replace("\\.fasta\\.F2B\\.ChrOnly\\.fa$",".fa") |> strip_ext() |> tolower()
fixm <- c("daTanVulg1.hap2.1.fasta" = "daTanVulg1.hap1.1.fa")

tr <- read.tree(TREE_F)
tr$tip.label <- normalize_id(ifelse(tr$tip.label %in% names(fixm), unname(fixm[tr$tip.label]), tr$tip.label))
tip_pfx <- sub("\\d.*$", "", tolower(tr$tip.label))
lk <- function(pp){h<-tr$tip.label[tip_pfx==pp];if(length(h)==1)return(h)
  h<-tr$tip.label[startsWith(tr$tip.label,pp)];if(length(h)>=1)return(h[1]);NA_character_}

# ── trait assembly (AGG: ian = dominant family; piotr = copy-number-weighted mean) ──
cf <- read.csv(CF_F, stringsAsFactors = FALSE)
for (c in c("regimentation_score","HOR_score","count_total")) cf[[c]] <- suppressWarnings(as.numeric(cf[[c]]))
cf$pfx <- sub("\\d.*$","",tolower(cf$fasta))
sd <- read.csv(SD_F, stringsAsFactors = FALSE)
sd$pfx <- sub("\\d.*$","",tolower(sd$fasta))
if (AGG == "ian") {
  cfa <- cf %>% group_by(pfx) %>% slice_max(count_total, n=1, with_ties=FALSE) %>% ungroup() %>%
    transmute(pfx, is_holo = is_holocentric=="TRUE", regi=regimentation_score, hor=HOR_score)
  sda <- sd %>% group_by(pfx) %>% slice_max(n, n=1, with_ties=FALSE) %>% ungroup() %>%
    transmute(pfx, mono_len_bp=width, sat_gc=100*gc., genome_gc=genome.gc,
              genome_mb=genome.bp/1e6, chr_n=chrs)
} else {
  cfa <- cf %>% group_by(pfx) %>%
    summarise(is_holo = first(is_holocentric=="TRUE"),
              regi = weighted.mean(regimentation_score, count_total, na.rm=TRUE),
              hor  = weighted.mean(HOR_score,           count_total, na.rm=TRUE), .groups="drop")
  sda <- sd %>% group_by(pfx) %>%
    summarise(mono_len_bp = weighted.mean(width,     n, na.rm=TRUE),
              sat_gc      = weighted.mean(100*gc.,   n, na.rm=TRUE),
              genome_gc = first(genome.gc), genome_mb = first(genome.bp)/1e6,
              chr_n = first(chrs), .groups="drop")
}

meta <- suppressMessages(read_excel(XLSX, sheet="Sheet1"))
meta$base <- normalize_id(ifelse(meta$fasta %in% names(fixm), unname(fixm[meta$fasta]), meta$fasta))
t1 <- setNames(meta$taxa1, meta$base)
clade_of <- function(lbl){t<-t1[lbl]
  ifelse(t %in% c("Actinopterygii","Aves","Reptilia","Mammalia","Amphibia"),"Vertebrata",
  ifelse(t %in% c("Monocots","Dicots","Bryophyta","Algae"),"Viridiplantae",
  ifelse(t=="Fungi","Fungi",
  ifelse(t %in% c("Alveolata","Discoba"),"Protist",
  ifelse(is.na(t),"Other","Invertebrata")))))}

tab <- tibble(label = tr$tip.label) %>%
  left_join(cfa %>% mutate(tip=vapply(pfx,lk,character(1))) %>% filter(!is.na(tip)) %>%
              transmute(label=tip, is_holo, regi, hor), by="label") %>%
  left_join(sda %>% mutate(tip=vapply(pfx,lk,character(1))) %>% filter(!is.na(tip)) %>%
              transmute(label=tip, mono_len_bp, sat_gc, genome_gc, genome_mb, chr_n), by="label") %>%
  mutate(clade = clade_of(label),
         regi = ifelse(is_holo %in% TRUE, NA, regi),
         hor  = ifelse(is_holo %in% TRUE, NA, hor))

write.table(tab %>% select(-is_holo), TAB_F, sep="\t", quote=FALSE, row.names=FALSE)
cat("Wrote trait table:", TAB_F, " (", sum(rowSums(!is.na(tab[,c('regi','hor','mono_len_bp')]))>0), "species with data )\n")

# ── correlations: naive Pearson + phylo (PGLS lambda) + PIC r ────────────────────
traits <- c("regi","hor","mono_len_bp","sat_gc","genome_gc","genome_mb","chr_n")
logv   <- c(mono_len_bp=TRUE, genome_mb=TRUE)       # log10 these before correlating
trans  <- function(v, nm) if (isTRUE(logv[nm])) log10(v) else v

pic_rp <- function(phy, x, y) {          # Felsenstein PIC correlation + test (through origin)
  cx <- pic(x, phy); cy <- pic(y, phy)
  r  <- sum(cx*cy)/sqrt(sum(cx^2)*sum(cy^2))
  p  <- summary(lm(cy ~ cx + 0))$coefficients["cx","Pr(>|t|)"]   # contrasts regression through origin
  c(r=r, p=p)
}
pgls_fit <- function(phy, d) {           # returns lambda, slope p; corPagel, fallback BM
  out <- tryCatch({
    m <- gls(y ~ x, data=d, correlation=corPagel(1, phy, form=~sp, fixed=FALSE), method="ML")
    c(lambda=as.numeric(coef(m$modelStruct$corStruct, unconstrained=FALSE)),
      p=summary(m)$tTable["x","p-value"])
  }, error=function(e) tryCatch({
    m <- gls(y ~ x, data=d, correlation=corBrownian(1, phy, form=~sp), method="ML")
    c(lambda=1, p=summary(m)$tTable["x","p-value"])
  }, error=function(e) c(lambda=NA, p=NA)))
  out
}

res <- list()
for (i in 1:(length(traits)-1)) for (j in (i+1):length(traits)) {
  a <- traits[i]; b <- traits[j]
  dd <- tab %>% select(label, xa=all_of(a), xb=all_of(b)) %>%
    filter(is.finite(xa), is.finite(xb), label %in% tr$tip.label)
  if (nrow(dd) < 15) next
  x <- trans(dd$xa, a); y <- trans(dd$xb, b); names(x) <- names(y) <- dd$label
  phy <- keep.tip(tr, dd$label)
  x <- x[phy$tip.label]; y <- y[phy$tip.label]
  ct <- cor.test(x, y)
  rp <- pic_rp(phy, x, y)
  d  <- data.frame(sp=phy$tip.label, x=x, y=y, row.names=phy$tip.label)
  pg <- pgls_fit(phy, d)
  lam <- unname(pg["lambda"]); pgp <- unname(pg["p"])
  if (!is.finite(pgp) || pgp <= 0 || (is.finite(lam) && lam > 1.001)) pgp <- NA  # guard corPagel artifacts
  res[[paste(a,b)]] <- data.frame(
    trait_x=a, trait_y=b, n=nrow(dd),
    pearson_r=round(unname(ct$estimate),3), pearson_p=signif(ct$p.value,3),
    phylo_r_pic=round(unname(rp["r"]),3), phylo_p_pic=signif(unname(rp["p"]),3),
    pgls_lambda=round(lam,3), pgls_p=signif(pgp,3), stringsAsFactors=FALSE)
}
cor_tbl <- do.call(rbind, res); rownames(cor_tbl) <- NULL
write.table(cor_tbl, COR_F, sep="\t", quote=FALSE, row.names=FALSE)
cat("\nPairwise correlations (naive vs phylogenetic):\n"); print(cor_tbl)
cat("\nWrote:", COR_F, "\n")

# ── Partial (multiple) PGLS: regi / hor / monomer length trio ───────────────────
# Does regimentation<->HOR survive controlling for monomer length, and vice versa?
# Multiple-predictor PGLS (corPagel); read each predictor's PARTIAL slope + p.
mono_log <- function(df) { df$mono_len_bp <- log10(df$mono_len_bp); df }
trio <- tab %>% select(label, regi, hor, mono_len_bp) %>%
  filter(is.finite(regi), is.finite(hor), is.finite(mono_len_bp), label %in% tr$tip.label) %>%
  mono_log()
phy3 <- keep.tip(tr, trio$label); trio <- trio[match(phy3$tip.label, trio$label), ]
trio$sp <- phy3$tip.label

partial_pgls <- function(form) {
  m <- tryCatch(gls(form, data=trio, correlation=corPagel(1, phy3, form=~sp, fixed=FALSE), method="ML"),
                error=function(e) gls(form, data=trio, correlation=corBrownian(1, phy3, form=~sp), method="ML"))
  tt <- summary(m)$tTable
  lam <- tryCatch({v <- as.numeric(coef(m$modelStruct$corStruct, unconstrained=FALSE))
                   if (length(v) == 0) 1 else v[1]}, error=function(e) NA)  # corBrownian fallback -> lambda=1
  preds <- rownames(tt)[rownames(tt) != "(Intercept)"]
  do.call(rbind, lapply(preds, function(p) data.frame(
    response = as.character(form)[2], predictor = p, n = nrow(trio),
    partial_slope = round(tt[p,"Value"],4), partial_p = signif(tt[p,"p-value"],3),
    pgls_lambda = round(lam,3), stringsAsFactors=FALSE)))
}
part <- rbind(
  partial_pgls(regi ~ hor + mono_len_bp),          # is HOR still assoc. w/ regi given monomer length?
  partial_pgls(hor  ~ regi + mono_len_bp),
  partial_pgls(mono_len_bp ~ regi + hor)           # is monomer length assoc. w/ HOR-ness independently?
)
rownames(part) <- NULL
write.table(part, PART_F, sep="\t", quote=FALSE, row.names=FALSE)
cat("\nPartial PGLS (regi/hor/log10 monomer trio, n=", nrow(trio), "):\n", sep=""); print(part)
cat("\nWrote:", PART_F, "\n")
