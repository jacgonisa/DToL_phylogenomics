#!/usr/bin/env Rscript
# 40c_all_models_3state_mcmctree.R
# Same 3-state ASR model comparison as 40b, but on the MCMCtree 325-sp chronogram
# (Bayesian relaxed-clock dating) instead of chronos-correlated lambda=0.1.
# Checks whether the model-selection / bidirectional-cycling result is robust to the
# dating method. Tips matched to the state annotation by alphabetic ToLID prefix
# (325/325), since MCMCtree uses different version suffixes than tree_renamed.

suppressPackageStartupMessages({
  library(ape); library(phytools); library(dplyr); library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
sf   <- sub("^--file=", "", args[grep("^--file=", args)])
asr_root <- if (length(sf)) normalizePath(file.path(dirname(sf), "..")) else getwd()
in_dir   <- file.path(asr_root, "inputs")
out_dir  <- file.path(asr_root, "outputs", "all_models_3state_mcmctree")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
MC_TREE  <- file.path(in_dir, "mcmctree_325sp_chronogram.nwk")   # MCMCtree Bayesian chronogram (copied into repo)

states3   <- c("H","Sat","Trans")
map_arch3 <- c(Holocentric = "H", Satellite = "Sat", Transposon = "Trans")  # rest -> NA
pfx <- function(x) sub("[0-9].*$", "", tolower(x))   # alphabetic ToLID stub (unique per species)

aicc_calc <- function(logL, k, n) 2*k - 2*logL + 2*k*(k+1)/(n - k - 1)

# ── MCMCtree chronogram (rooted, dated); fix rounding non-ultrametricity ───────
mc <- read.tree(MC_TREE)
mc <- multi2di(mc, random = FALSE)
if (!is.ultrametric(mc)) mc <- force.ultrametric(mc, method = "nnls")  # rounding only
mc_pfx <- pfx(mc$tip.label)

# state per MCMCtree tip, via prefix -> full annotation
ann <- read.delim(file.path(in_dir, "full_chronos_correlated", "branch_symbol_anno.tsv"),
                  stringsAsFactors = FALSE)
colnames(ann)[c(1,8)] <- c("tip","architecture")
ann_pfx <- pfx(ann$tip)
arch_by_pfx <- setNames(ann$architecture, ann_pfx)

# dataset membership (which prefixes belong to metazoa / viridiplantae) from 40b inputs
sub_pfx <- function(ds) pfx(read.tree(file.path(in_dir, ds, "tree_renamed.nw"))$tip.label)
members <- list(
  `Full tree`   = mc$tip.label,
  Metazoa       = mc$tip.label[mc_pfx %in% sub_pfx("metazoa_chronos_correlated")],
  Viridiplantae = mc$tip.label[mc_pfx %in% sub_pfx("viridiplantae_chronos_correlated")]
)

load3 <- function(ds) {
  keep_tips <- members[[ds]]
  tr <- keep.tip(mc, keep_tips)
  ch <- factor(map_arch3[arch_by_pfx[pfx(tr$tip.label)]], levels = states3)
  names(ch) <- tr$tip.label
  kp <- !is.na(ch); tr <- drop.tip(tr, tr$tip.label[!kp])
  ch <- factor(as.character(ch[kp]), levels = states3); names(ch) <- tr$tip.label
  list(tree = tr, char = ch, n = length(ch))
}

# ── design matrices (identical to 40b) ────────────────────────────────────────
mk_ard <- function() { dm <- matrix(0L,3,3,dimnames=list(states3,states3)); k<-1L
  for(i in states3) for(j in states3) if(i!=j){dm[i,j]<-k; k<-k+1L}; dm }
mk_irrevH <- function(){ dm<-mk_ard(); dm["H","Sat"]<-0L; dm["H","Trans"]<-0L; dm }
mk_irrevH_noDirectST <- function(){ dm<-mk_irrevH(); dm["Sat","Trans"]<-0L; dm["Trans","Sat"]<-0L; dm }
mk_irrevH_noSatToTrans <- function(){ dm<-mk_irrevH(); dm["Sat","Trans"]<-0L; dm }
mk_irrevH_noTransToSat <- function(){ dm<-mk_irrevH(); dm["Trans","Sat"]<-0L; dm }
models <- list(ER="ER", SYM="SYM", ARD="ARD", ARD_irrevH=mk_irrevH(),
  ARD_irrevH_noDirectST=mk_irrevH_noDirectST(),
  ARD_irrevH_noSatToTrans=mk_irrevH_noSatToTrans(), ARD_irrevH_noTransToSat=mk_irrevH_noTransToSat())

datasets <- c("Full tree","Metazoa","Viridiplantae")
res <- list()
for (ds in datasets) {
  cat("\n===", ds, "(MCMCtree) ===\n")
  dat <- load3(ds)
  cat(sprintf("  n=%d  (H=%d Sat=%d Trans=%d)  root age=%.0f Mya\n", dat$n,
              sum(dat$char=="H"), sum(dat$char=="Sat"), sum(dat$char=="Trans"),
              max(node.depth.edgelength(dat$tree))))
  for (mn in names(models)) {
    fit <- tryCatch(fitMk(dat$tree, dat$char, model=models[[mn]], states=states3,
                          control=list(maxit=3000)), error=function(e) NULL)
    if (is.null(fit)) { cat(sprintf("  %-24s FAILED\n", mn)); next }
    logL <- as.numeric(fit$logLik)
    idx  <- fit$index.matrix; k <- length(unique(idx[!is.na(idx) & idx > 0]))
    ac   <- aicc_calc(logL, k, dat$n)
    cat(sprintf("  %-24s k=%d  logL=%.2f  AICc=%.2f\n", mn, k, logL, ac))
    res[[paste(ds,mn)]] <- data.frame(dataset=ds, model=mn, k=k, logL=round(logL,3),
      AIC=round(2*k-2*logL,2), AICc=round(ac,2), stringsAsFactors=FALSE)
  }
}

tbl <- bind_rows(res) %>% group_by(dataset) %>%
  mutate(dAICc = round(AICc-min(AICc),2),
         w_AICc = round(exp(-0.5*(AICc-min(AICc)))/sum(exp(-0.5*(AICc-min(AICc)))),3),
         best = dAICc==0) %>% ungroup() %>% arrange(dataset, AICc)
write.table(tbl, file.path(out_dir,"all_models_3state_aicc_mcmctree.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("\n========= 3-STATE MODEL TABLE (MCMCtree) =========\n")
print(as.data.frame(tbl %>% select(dataset,model,k,logL,AICc,dAICc,w_AICc)))
cat("\nWrote:", file.path(out_dir,"all_models_3state_aicc_mcmctree.tsv"), "\n")
