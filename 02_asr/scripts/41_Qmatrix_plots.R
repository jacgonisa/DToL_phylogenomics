#!/usr/bin/env Rscript
# 41_Qmatrix_plots.R
# Publication-quality Q matrix visualisations for ARD_irrevH (best model):
#   A. Heatmap of off-diagonal rates (log scale, per dataset)
#   B. ggraph rate-flow diagrams (arrow width ∝ sqrt(rate), labelled)
#   C. Combined panel figure

suppressPackageStartupMessages({
  library(ape); library(phytools); library(dplyr); library(tidyr)
  library(ggplot2); library(readxl); library(stringr)
  library(igraph); library(ggraph); library(patchwork); library(scales)
})
set.seed(42)

root    <- normalizePath(file.path(getwd(), "ASR_March2026_327species"))
out_dir <- file.path(root, "outputs", "Qmatrix_plots")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

state_levels <- c("H","Mixed","Sat","Trans")
state_colors <- c(H="#f46d43", Mixed="#8e44ad", Sat="#ff2d87", Trans="#4f84e8")
state_labels <- c(H="Holocentric", Mixed="Sat/Trans\n(Mixed)", Sat="Satellite", Trans="Transposon")
datasets       <- c("full","metazoa","viridiplantae")
dataset_labels <- c(full="Full tree", metazoa="Metazoa", viridiplantae="Viridiplantae")

# ── Annotations ──────────────────────────────────────────────────────────────
anno <- read.csv(file.path(root,"inputs","full","branch_symbol_anno.tsv"),
                 sep="\t", header=TRUE, stringsAsFactors=FALSE)
anno$tree_label <- gsub(" ","",anno$Species)
map_arch <- c("Holocentric"="H","Satellite/transposon"="Mixed",
              "Satellite"="Sat","Transposon"="Trans")
label_fix_map <- c("daTanVulg1.hap2.1.fasta"="daTanVulg1.hap1.1.fa",
                   "rosCan_S27_v1.fasta.F2B.ChrOnly.fa"="rosCan_S27_v1")

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

# ── Design matrix: ARD_irrevH ─────────────────────────────────────────────────
dm_cyc <- matrix(0L, 4, 4, dimnames=list(state_levels, state_levels))
k <- 1L
for (i in state_levels) for (j in state_levels)
  if (i != j) { dm_cyc[i,j] <- k; k <- k+1L }
dm_cyc["H","Mixed"] <- 0L; dm_cyc["H","Sat"] <- 0L; dm_cyc["H","Trans"] <- 0L

# ── Fit models and extract Q matrices ────────────────────────────────────────
cat("Fitting ARD_irrevH for all datasets...\n")
Q_list <- list()

for (ds in datasets) {
  cat(" ", ds, "...")
  dat <- load_dataset(ds)
  fit <- fitMk(dat$tree, dat$char, model=dm_cyc,
               states=state_levels, control=list(maxit=3000))
  Q <- matrix(0, 4, 4, dimnames=list(state_levels, state_levels))
  for (i in 1:4) for (j in 1:4)
    if (dm_cyc[i,j] > 0) Q[i,j] <- fit$rates[dm_cyc[i,j]]
  diag(Q) <- -rowSums(Q)
  Q_list[[ds]] <- Q
  cat(sprintf(" logL=%.3f\n", fit$logLik))
}

# ── Long-format rate table ────────────────────────────────────────────────────
rate_df <- bind_rows(lapply(datasets, function(ds) {
  Q <- Q_list[[ds]]
  expand.grid(from=state_levels, to=state_levels, stringsAsFactors=FALSE) %>%
    mutate(dataset = dataset_labels[ds],
           rate    = mapply(function(f,t) Q[f,t], from, to),
           is_diag = from == to)
})) %>%
  mutate(
    from    = factor(from,    levels=state_levels),
    to      = factor(to,      levels=state_levels),
    dataset = factor(dataset, levels=dataset_labels)
  )

# Save full Q matrices
write.table(rate_df %>% filter(!is_diag) %>% select(dataset,from,to,rate),
            file.path(out_dir,"Q_matrix_ARD_irrevH.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ── A. Heatmap ────────────────────────────────────────────────────────────────
# Off-diagonal only; diagonal shown as grey
hm_df <- rate_df %>%
  mutate(
    rate_plot  = ifelse(is_diag, NA_real_, pmax(rate, 0)),
    label_text = case_when(
      is_diag  ~ "",
      rate <= 0 ~ "0",
      TRUE      ~ sprintf("%.4f", rate)
    )
  )

p_heatmap <- ggplot(hm_df, aes(x=to, y=from, fill=rate_plot)) +
  facet_wrap(~dataset, nrow=1) +
  geom_tile(color="white", linewidth=0.9) +
  geom_text(aes(label=label_text), size=3, color="white", fontface="bold") +
  scale_fill_gradientn(
    colors  = c("#f7fbff","#6baed6","#2171b5","#08306b"),
    na.value= "grey88",
    trans   = "log1p",
    name    = "Rate\n(per MY)",
    labels  = label_number(accuracy=0.0001)
  ) +
  scale_x_discrete(labels=state_labels) +
  scale_y_discrete(limits=rev(state_levels), labels=rev(state_labels)) +
  labs(title="Q matrix — ARD_irrevH (best-fit model)",
       subtitle="Off-diagonal rates (per million years). Grey = not estimated (H irreversible, diagonal excluded).",
       x="To state", y="From state") +
  theme_bw(base_size=12) +
  theme(
    strip.text      = element_text(face="bold", size=12),
    axis.text.x     = element_text(angle=30, hjust=1, size=9),
    axis.text.y     = element_text(size=9),
    legend.position = "right",
    panel.grid      = element_blank()
  )

ggsave(file.path(out_dir,"A_Qmatrix_heatmap.pdf"), p_heatmap, width=13, height=4.5)
ggsave(file.path(out_dir,"A_Qmatrix_heatmap.png"), p_heatmap, width=13, height=4.5, dpi=300)
cat("Saved heatmap\n")

# ── B. ggraph rate-flow diagrams ──────────────────────────────────────────────
# Node layout: Sat (left), Trans (right), Mixed (bottom-centre), H (top-centre)
node_pos <- data.frame(
  state = state_levels,
  x     = c( 0,  0, -1,  1),   # H, Mixed, Sat, Trans
  y     = c( 1, -1,  0,  0),
  label = state_labels,
  color = state_colors,
  stringsAsFactors = FALSE
)

flow_plots <- lapply(datasets, function(ds) {
  Q <- Q_list[[ds]]

  edges <- rate_df %>%
    filter(dataset == dataset_labels[ds], !is_diag, rate > 1e-10) %>%
    select(from, to, rate) %>%
    mutate(
      log_rate   = log10(rate),
      rate_label = sprintf("%.4f", rate)
    )

  g <- graph_from_data_frame(edges, directed=TRUE,
                              vertices=node_pos %>% rename(name=state))

  ggraph(g, layout="manual", x=node_pos$x, y=node_pos$y) +
    geom_edge_arc(
      aes(width=rate, color=log_rate,
          label=sprintf("%.4f", rate)),
      strength=0.25,
      arrow=arrow(length=unit(3,"mm"), type="closed"),
      end_cap=circle(7,"mm"), start_cap=circle(7,"mm"),
      label_size=2.8, label_colour="grey20",
      angle_calc="along", label_dodge=unit(2.5,"mm"),
      show.legend=c(width=FALSE, color=TRUE)
    ) +
    geom_node_point(aes(color=I(color)), size=14) +
    geom_node_text(aes(label=label), size=2.8, fontface="bold",
                   color="white", lineheight=0.85) +
    scale_edge_width_continuous(range=c(0.5, 4)) +
    scale_edge_color_gradient2(
      low="#74add1", mid="#ffffbf", high="#d73027",
      midpoint=median(log10(edges$rate)),
      name="log₁₀(rate)",
      guide=guide_edge_colorbar(title.position="top")
    ) +
    labs(title=dataset_labels[ds]) +
    coord_cartesian(xlim=c(-1.7,1.7), ylim=c(-1.6,1.6)) +
    theme_void(base_size=11) +
    theme(
      plot.title       = element_text(face="bold", size=13, hjust=0.5),
      legend.position  = "bottom",
      legend.key.width = unit(1.2,"cm")
    )
})
names(flow_plots) <- datasets

p_flow <- wrap_plots(flow_plots, nrow=1) +
  plot_annotation(
    title    = "Rate-flow diagrams — ARD_irrevH",
    subtitle = "Arrow width ∝ rate (per MY). Colour = log₁₀(rate): blue=slow → red=fast.",
    theme    = theme(plot.title    = element_text(face="bold", size=14),
                     plot.subtitle = element_text(size=10, color="grey40"))
  )

ggsave(file.path(out_dir,"B_flow_diagrams.pdf"), p_flow, width=15, height=6)
ggsave(file.path(out_dir,"B_flow_diagrams.png"), p_flow, width=15, height=6, dpi=300)
cat("Saved flow diagrams\n")

# ── C. Bar chart: rates side-by-side per transition ──────────────────────────
bar_df <- rate_df %>%
  filter(!is_diag, rate > 0) %>%
  mutate(
    transition = paste0(from," → ",to),
    is_direct  = (from=="Sat" & to=="Trans") | (from=="Trans" & to=="Sat"),
    direction  = case_when(
      from %in% c("Sat","Mixed") & to %in% c("Trans","H") ~ "Forward / exit",
      TRUE ~ "Return"
    )
  ) %>%
  filter(!(from=="H"))   # H→X = 0 always, skip

p_bar <- ggplot(bar_df, aes(x=reorder(transition, rate), y=rate,
                             fill=as.character(from))) +
  facet_wrap(~dataset, nrow=1, scales="free_y") +
  geom_col(width=0.7) +
  geom_text(aes(label=sprintf("%.4f", rate)), hjust=-0.1, size=2.8) +
  coord_flip() +
  scale_fill_manual(values=state_colors, name="From state",
                    labels=state_labels[names(state_colors)]) +
  scale_y_continuous(expand=expansion(mult=c(0,0.3))) +
  labs(title="C: Individual rates from ARD_irrevH Q matrix",
       x=NULL, y="Rate (per MY)") +
  theme_bw(base_size=11) +
  theme(strip.text=element_text(face="bold"),
        panel.grid.minor=element_blank(),
        panel.grid.major.y=element_blank(),
        legend.position="bottom")

ggsave(file.path(out_dir,"C_rate_barplot.pdf"), p_bar, width=14, height=7)
ggsave(file.path(out_dir,"C_rate_barplot.png"), p_bar, width=14, height=7, dpi=300)
cat("Saved rate barplot\n")

# ── Combined ─────────────────────────────────────────────────────────────────
p_combined <- p_heatmap / p_flow / p_bar +
  plot_annotation(
    title = "Centromere transition rate matrices — ARD_irrevH best-fit model",
    theme = theme(plot.title=element_text(face="bold", size=15))
  ) +
  plot_layout(heights=c(1, 1.3, 1.4))

ggsave(file.path(out_dir,"combined_Qmatrix.pdf"), p_combined, width=15, height=18)
ggsave(file.path(out_dir,"combined_Qmatrix.png"), p_combined, width=15, height=18, dpi=300)
cat("Saved combined figure\n")

cat("\nOutputs in:", out_dir, "\n")
