#!/usr/bin/env Rscript
# trait_pgls_correlations_325sp.R
# Build a per-species centromere/genome trait table (dominant-array = Ian aggregation)
# and compute pairwise correlations BOTH naively (Pearson) and phylogenetically
# corrected (PGLS via nlme::gls + Pagel's lambda; phylo r from Felsenstein PICs).
# Outputs feed the seaborn pairplot (trait_pairplot_seaborn_325sp.py).

suppressPackageStartupMessages({
  library(ape); library(caper); library(phytools); library(readxl); library(dplyr); library(stringr)
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

pic_rp <- function(phy, x, y) {          # phylo correlation from BM covariance (phytools::phyl.vcv)
  X <- cbind(x[phy$tip.label], y[phy$tip.label]); rownames(X) <- phy$tip.label
  r <- cov2cor(phyl.vcv(X, vcv(phy), lambda = 1)$R)[1, 2]         # symmetric evolutionary correlation
  cx <- pic(x, phy); cy <- pic(y, phy)                           # contrasts only for the significance test
  p  <- summary(lm(cy ~ cx + 0))$coefficients["cx", "Pr(>|t|)"]   # BM-correlation p (through origin)
  c(r = r, p = p)
}
pgls_fit <- function(phy, d) {           # caper::pgls, Pagel's lambda by ML (Brownian fallback)
  cd <- comparative.data(phy, d, names.col="sp", vcv=TRUE, warn.dropped=FALSE)
  m  <- tryCatch(pgls(y ~ x, data=cd, lambda="ML"), error=function(e) NULL)
  if (is.null(m)) m <- tryCatch(pgls(y ~ x, data=cd, lambda=1), error=function(e) NULL)  # BM if ML fails
  if (is.null(m)) return(c(lambda=NA, p=NA))
  p <- summary(m)$coefficients["x","Pr(>|t|)"]
  if (is.finite(p) && p == 0) p <- 2.2e-16                 # floor underflow (p ~ machine eps)
  c(lambda=unname(m$param["lambda"]), p=p)
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
  if (!is.finite(pgp)) pgp <- NA        # caper bounds lambda to [0,1]; only guard failed fits
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

cd3 <- comparative.data(phy3, as.data.frame(trio), names.col="sp", vcv=TRUE, warn.dropped=FALSE)
partial_pgls <- function(form) {         # caper::pgls multiple regression; partial slope + p per predictor
  m  <- pgls(form, data=cd3, lambda="ML")
  co <- summary(m)$coefficients
  lam <- unname(m$param["lambda"])
  preds <- rownames(co)[rownames(co) != "(Intercept)"]
  do.call(rbind, lapply(preds, function(p) data.frame(
    response = as.character(form)[2], predictor = p, n = nrow(trio),
    partial_slope = round(co[p,"Estimate"],4), partial_p = signif(co[p,"Pr(>|t|)"],3),
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

# ── Phylogenetic pairwise vs partial correlation for the trio (phytools::phyl.vcv) ──
# Joint BM correlation matrix R3; pairwise = R3[i,j], partial = precision-matrix formula
# (-P[i,j]/sqrt(P[i,i]P[j,j])). Significance stars from the PGLS partial-slope p-values.
X3 <- cbind(regi = trio$regi, hor = trio$hor, mono = trio$mono_len_bp); rownames(X3) <- trio$sp
R3 <- cov2cor(phyl.vcv(X3, vcv(phy3), lambda = 1)$R)
P3 <- solve(R3)
pcor3 <- function(i, j) -P3[i, j] / sqrt(P3[i, i] * P3[j, j])
pget <- function(resp,pred) part$partial_p[part$response==resp & part$predictor==pred][1]
trio_tbl <- data.frame(
  pair       = c("regimentation~HOR","regimentation~monomer","HOR~monomer"),
  control_for= c("monomer length","HOR score","regimentation"),
  n          = nrow(trio),
  pic_r_pairwise = round(c(R3["regi","hor"], R3["regi","mono"], R3["hor","mono"]), 3),
  pic_r_partial  = round(c(pcor3("regi","hor"), pcor3("regi","mono"), pcor3("hor","mono")), 3),
  partial_p  = signif(c(pget("regi","hor"), pget("regi","mono_len_bp"), pget("hor","mono_len_bp")), 3),
  stringsAsFactors = FALSE)
PTRIO_F <- file.path(OUT_DIR, sprintf("phylo_pairwise_vs_partial_trio_%s.tsv", AGG))
write.table(trio_tbl, PTRIO_F, sep="\t", quote=FALSE, row.names=FALSE)
cat("\nPhylo pairwise vs partial r (trio):\n"); print(trio_tbl)
cat("\nWrote:", PTRIO_F, "\n")

# ── Phylogenetic signal per trait (K, Pagel lambda) ─────────────────────────────
# Explains WHICH correlations phylogenetic correction affects: only tree-structured
# (high-lambda) traits can carry shared-ancestry artifacts. Low-lambda HOR traits
# barely change under correction.
sig_tab <- do.call(rbind, lapply(traits, function(nm){
  v <- tab[[nm]]; names(v) <- tab$label; v <- v[is.finite(v)]; v <- trans(v, nm)
  phy <- keep.tip(tr, intersect(names(v), tr$tip.label)); v <- v[phy$tip.label]
  data.frame(trait=nm, n=length(v),
             K=round(as.numeric(phylosig(phy, v, method="K")),3),
             lambda=round(phylosig(phy, v, method="lambda")$lambda,3),
             structured=ifelse(phylosig(phy, v, method="lambda")$lambda>=0.5,"strong","weak"),
             stringsAsFactors=FALSE)
}))
SIG_F <- file.path(OUT_DIR, sprintf("trait_phylo_signal_%s.tsv", AGG))
write.table(sig_tab, SIG_F, sep="\t", quote=FALSE, row.names=FALSE)
cat("\nPhylogenetic signal per trait:\n"); print(sig_tab)
cat("\nWrote:", SIG_F, "\n")
