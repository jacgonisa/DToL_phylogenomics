#!/usr/bin/env Rscript
# Select the chronos smoothing parameter lambda by PHIIC (Paradis 2013): fit the
# 62-constraint chronogram at a range of lambda for both rate models, extract the
# penalized log-likelihood and PHIIC (best of several restarts), and report.
suppressPackageStartupMessages({library(ape); library(phytools)})
SP <- "/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/01_species_tree"
CTRL <- chronos.control(iter.max=1e6, eval.max=1e6, dual.iter.max=1e4, tol=1e-8)
tree <- multi2di(midpoint.root(read.tree(file.path(SP,"fast_species_tree_325sp_renamed.nwk"))), random=FALSE)
tip <- function(b){h<-grep(paste0("^",gsub("\\.","\\\\.",b),"(\\.fa|\\.fasta)?$"),tree$tip.label,value=TRUE); if(length(h)==0) NA_character_ else h[1]}
node_of <- function(a,b) if(is.na(tip(a))||is.na(tip(b))) NA_integer_ else as.integer(getMRCA(tree,c(tip(a),tip(b))))
cal <- read.delim(file.path(SP,"over_calib.tsv"), stringsAsFactors=FALSE)
cal$node <- mapply(function(l,a,b){ if(l=="Eukaryota_root") return(as.integer(Ntip(tree)+1L)); node_of(a,b) }, cal$label, cal$tip_a, cal$tip_b)
cal <- cal[!is.na(cal$node),]
mn <- tapply(cal$age_min, cal$node, max); mx <- tapply(cal$age_max, cal$node, min)
nodes <- as.integer(names(mn)); keep <- mn<=mx
calib <- makeChronosCalib(tree, node=nodes[keep], age.min=mn[keep], age.max=mx[keep])

get_ph <- function(f){p<-attr(f,"PHIIC"); if(is.list(p)) p$PHIIC else p}
get_ll <- function(f){p<-attr(f,"PHIIC"); if(is.list(p)&&!is.null(p$logLik)) p$logLik else attr(f,"ploglik")}

res <- data.frame()
for (model in c("correlated","relaxed")) {
  for (lam in c(0, 0.1, 1, 10)) {
    best_ph<-Inf; best_ll<-NA
    for (i in 1:5) {
      f <- tryCatch(suppressWarnings(chronos(tree, lambda=lam, model=model, calibration=calib, control=CTRL)), error=function(e) NULL)
      if (is.null(f)) next
      ph<-get_ph(f); if (!is.null(ph) && ph<best_ph){best_ph<-ph; best_ll<-get_ll(f)}
    }
    res <- rbind(res, data.frame(model=model, lambda=lam, logLik=round(best_ll,3), PHIIC=round(best_ph,3)))
    cat(sprintf("%-11s lambda=%-4s  logLik=%.3f  PHIIC=%.3f\n", model, lam, best_ll, best_ph))
  }
}
write.table(res, file.path(SP,"outputs/calibration_qc/lambda_selection_phiic.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("\nLowest PHIIC per model:\n")
for (m in unique(res$model)) { s<-res[res$model==m,]; b<-s[which.min(s$PHIIC),]; cat(sprintf("  %-11s -> lambda=%s (PHIIC=%.3f)\n", m, b$lambda, b$PHIIC)) }
