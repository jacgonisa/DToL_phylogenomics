#!/usr/bin/env Rscript
# 40b_all_models_3state.R
# 3-state ASR model comparison (H / Sat / Trans; Mixed + Unknown pruned to NA)
# on the chronos-correlated chronogram, full 8-model set:
#   ER, SYM, ARD, ARD_irrevH, ARD_irrevH_noDirectST, ARD_irrevH_symST,
#   ARD_irrevH_noSatToTrans, ARD_irrevH_noTransToSat
# Reads this repo's own 02_asr/inputs/; writes outputs/all_models_3state/.

suppressPackageStartupMessages({
  library(ape); library(phytools); library(dplyr); library(ggplot2)
})

# repo 02_asr root = dir containing this script's parent
args <- commandArgs(trailingOnly = FALSE)
sf   <- sub("^--file=", "", args[grep("^--file=", args)])
asr_root <- if (length(sf)) normalizePath(file.path(dirname(sf), "..")) else getwd()
in_dir   <- file.path(asr_root, "inputs")
out_dir  <- file.path(asr_root, "outputs", "all_models_3state")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

states3   <- c("H","Sat","Trans")
map_arch3 <- c(Holocentric = "H", Satellite = "Sat", Transposon = "Trans")  # rest -> NA
datasets  <- c(full_chronos_correlated = "Full tree",
               metazoa_chronos_correlated = "Metazoa",
               viridiplantae_chronos_correlated = "Viridiplantae")

aicc_calc <- function(logL, k, n) 2*k - 2*logL + 2*k*(k+1)/(n - k - 1)

load3 <- function(ds) {
  ann <- read.delim(file.path(in_dir, ds, "branch_symbol_anno.tsv"), stringsAsFactors = FALSE)
  colnames(ann)[c(1,8)] <- c("tip","architecture")
  tr  <- read.tree(file.path(in_dir, ds, "tree_renamed.nw"))
  tr  <- midpoint.root(tr); tr <- multi2di(tr, random = FALSE)
  ch  <- factor(map_arch3[ann$architecture[match(tr$tip.label, ann$tip)]], levels = states3)
  names(ch) <- tr$tip.label
  kp  <- !is.na(ch); tr <- drop.tip(tr, tr$tip.label[!kp])
  ch  <- factor(as.character(ch[kp]), levels = states3); names(ch) <- tr$tip.label
  list(tree = tr, char = ch, n = length(ch))
}

# ── design matrices (unique index per free rate) ──────────────────────────────
mk_ard <- function() { dm <- matrix(0L,3,3,dimnames=list(states3,states3)); k<-1L
  for(i in states3) for(j in states3) if(i!=j){dm[i,j]<-k; k<-k+1L}; dm }
mk_irrevH <- function(){ dm<-mk_ard(); dm["H","Sat"]<-0L; dm["H","Trans"]<-0L; dm }
mk_irrevH_noDirectST <- function(){ dm<-mk_irrevH(); dm["Sat","Trans"]<-0L; dm["Trans","Sat"]<-0L; dm }
mk_irrevH_symST <- function(){ dm<-mk_irrevH(); dm["Trans","Sat"]<-dm["Sat","Trans"]; dm }
mk_irrevH_noSatToTrans <- function(){ dm<-mk_irrevH(); dm["Sat","Trans"]<-0L; dm }
mk_irrevH_noTransToSat <- function(){ dm<-mk_irrevH(); dm["Trans","Sat"]<-0L; dm }

models <- list(
  ER                      = "ER",
  SYM                     = "SYM",
  ARD                     = "ARD",
  ARD_irrevH              = mk_irrevH(),
  ARD_irrevH_noDirectST   = mk_irrevH_noDirectST(),
  ARD_irrevH_symST        = mk_irrevH_symST(),
  ARD_irrevH_noSatToTrans = mk_irrevH_noSatToTrans(),
  ARD_irrevH_noTransToSat = mk_irrevH_noTransToSat()
)

# ── fit ───────────────────────────────────────────────────────────────────────
res <- list()
for (ds in names(datasets)) {
  cat("\n===", datasets[ds], "===\n")
  dat <- load3(ds)
  cat(sprintf("  n=%d  (H=%d Sat=%d Trans=%d)\n", dat$n,
              sum(dat$char=="H"), sum(dat$char=="Sat"), sum(dat$char=="Trans")))
  for (mn in names(models)) {
    fit <- tryCatch(fitMk(dat$tree, dat$char, model = models[[mn]],
                          states = states3, control = list(maxit = 3000)),
                    error = function(e) NULL)
    if (is.null(fit)) { cat(sprintf("  %-24s FAILED\n", mn)); next }
    logL <- as.numeric(fit$logLik)
    idx  <- fit$index.matrix; k <- length(unique(idx[!is.na(idx) & idx > 0]))
    ac   <- aicc_calc(logL, k, dat$n)
    cat(sprintf("  %-24s k=%d  logL=%.2f  AICc=%.2f\n", mn, k, logL, ac))
    res[[paste(ds,mn)]] <- data.frame(dataset = datasets[ds], model = mn,
      k = k, logL = round(logL,3), AIC = round(2*k-2*logL,2), AICc = round(ac,2),
      stringsAsFactors = FALSE)
  }
}

tbl <- bind_rows(res) %>%
  group_by(dataset) %>%
  mutate(dAICc  = round(AICc - min(AICc), 2),
         w_AICc = round(exp(-0.5*(AICc-min(AICc))) / sum(exp(-0.5*(AICc-min(AICc)))), 3),
         dAIC   = round(AIC - min(AIC), 2),
         w_AIC  = round(exp(-0.5*(AIC-min(AIC)))  / sum(exp(-0.5*(AIC-min(AIC)))),  3),
         best   = dAICc == 0) %>%
  ungroup() %>% arrange(dataset, AICc)

write.table(tbl, file.path(out_dir, "all_models_3state_aicc.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n========= 3-STATE MODEL TABLE =========\n")
print(as.data.frame(tbl %>% select(dataset, model, k, logL, AICc, dAICc, w_AICc)))

# ── weight dot plot (AIC + AICc facets, red one-directional labels) ───────────
model_order <- c("ER","SYM","ARD","ARD_irrevH_noDirectST","ARD_irrevH_symST",
                 "ARD_irrevH_noSatToTrans","ARD_irrevH_noTransToSat","ARD_irrevH")
red_models  <- c("ARD_irrevH_noSatToTrans","ARD_irrevH_noTransToSat")
ds_pal   <- c(`Full tree`="#E41A1C", Metazoa="#377EB8", Viridiplantae="#2E7D32")
ds_shape <- c(`Full tree`=16, Metazoa=17, Viridiplantae=15)
lab_cols <- ifelse(model_order %in% red_models, "#c62828", "black")

plt <- tbl %>%
  tidyr::pivot_longer(c(w_AIC, w_AICc), names_to="criterion", values_to="w") %>%
  mutate(criterion = ifelse(criterion=="w_AIC","AIC","AICc"),
         model     = factor(model, levels=model_order),
         dataset   = factor(dataset, levels=names(ds_pal))) %>%
  ggplot(aes(w, model, colour=dataset, shape=dataset)) +
  geom_vline(xintercept=0, colour="grey85") +
  geom_point(size=3.4, alpha=0.92) +
  scale_colour_manual(values=ds_pal, name=NULL) +
  scale_shape_manual(values=ds_shape, name=NULL) +
  scale_x_continuous(limits=c(0,1.05), breaks=c(0,.25,.5,.75,1)) +
  scale_y_discrete(limits=rev(model_order)) +
  facet_wrap(~criterion, ncol=2) +
  labs(title="3-state ASR model weights — chronos-correlated (325-sp)",
       subtitle="Red labels = one-directional models (near-zero weight supports bidirectional Sat<->Trans cycling)",
       x="Akaike weight", y=NULL) +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
        legend.position="bottom", strip.text=element_text(face="bold"),
        axis.text.y=element_text(colour=rev(lab_cols)),
        plot.subtitle=element_text(size=9, colour="grey40"))

ggsave(file.path(out_dir,"all_models_3state_weights.pdf"), plt, width=10, height=5)
ggsave(file.path(out_dir,"all_models_3state_weights.png"), plt, width=10, height=5, dpi=300)
cat("\nOutputs in:", out_dir, "\n")
