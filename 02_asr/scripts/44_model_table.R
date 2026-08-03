#!/usr/bin/env Rscript
# 44_model_table.R
# Clean formatted table of all model fits across datasets.
# Reads existing all_models_aicc.tsv and produces a publication-quality PDF table.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(gridExtra); library(grid)
})

root    <- normalizePath(file.path(getwd(), "ASR_March2026_327species"))
tsv     <- file.path(root, "outputs", "all_models", "all_models_aicc.tsv")
out_dir <- file.path(root, "outputs", "all_models")

df <- read.table(tsv, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Build display table
tbl <- df %>%
  arrange(dataset, dAICc) %>%
  mutate(
    logL    = sprintf("%.2f",  logL),
    AICc    = sprintf("%.2f",  AICc),
    dAICc   = sprintf("%.2f",  dAICc),
    weight  = sprintf("%.3f",  weight),
    best    = ifelse(best, "✓", "")
  ) %>%
  select(Dataset=dataset, Model=model, Description=desc,
         k, logL, AICc, `ΔAICc`=dAICc, Weight=weight, Best=best)

# Split into one sub-table per dataset for readability
datasets <- unique(tbl$Dataset)

make_gtable <- function(sub_df, ds_name) {
  # Shade best-fit row
  best_row <- which(sub_df$Best == "✓")
  n_rows <- nrow(sub_df)
  n_cols <- ncol(sub_df)

  fill_mat <- matrix("white", nrow=n_rows+1, ncol=n_cols)
  fill_mat[1, ]          <- "#2c3e50"          # header
  if (length(best_row)) fill_mat[best_row+1, ] <- "#d5e8d4"  # best row light green

  txt_mat <- matrix("black", nrow=n_rows+1, ncol=n_cols)
  txt_mat[1, ] <- "white"

  tt <- ttheme_minimal(
    core    = list(bg_params = list(fill=fill_mat[-1,], col=NA),
                   fg_params = list(col=txt_mat[-1,], fontsize=9)),
    colhead = list(bg_params = list(fill=fill_mat[1,], col=NA),
                   fg_params = list(col=txt_mat[1,], fontsize=9, fontface="bold"))
  )

  gt <- tableGrob(sub_df %>% select(-Dataset), rows=NULL, theme=tt)

  title_grob <- textGrob(ds_name,
                          gp=gpar(fontsize=11, fontface="bold", col="#2c3e50"),
                          just="left", x=0.01)
  arrangeGrob(title_grob, gt, heights=unit(c(0.6,1)*cm(1), c("cm","null")),
              nrow=2)
}

grob_list <- lapply(datasets, function(ds) {
  make_gtable(tbl %>% filter(Dataset==ds), ds)
})

# Add a note at the bottom
note <- textGrob(
  paste0("Models fitted with fitMk() (phytools). AICc = AIC + 2k(k+1)/(n-k-1).\n",
         "ARD_irrevH: all-rates-different + H->X = 0 (H irreversible).\n",
         "TerminalTrans: ARD_irrevH + Trans->Sat = Trans->Mixed = 0 (Trans is a sink)."),
  gp=gpar(fontsize=7.5, col="grey40"), just="left", x=0.01
)

full_grob <- arrangeGrob(
  grobs = c(
    list(textGrob("Centromere evolution — model comparison (AICc)",
                  gp=gpar(fontsize=14, fontface="bold"), just="left", x=0.01)),
    grob_list,
    list(note)
  ),
  nrow = length(grob_list) + 2,
  heights = unit(c(1, rep(5, length(grob_list)), 1.8), "cm")
)

ggsave(file.path(out_dir, "model_comparison_table.pdf"),
       full_grob, width=13, height=4 + length(grob_list)*4.5)
ggsave(file.path(out_dir, "model_comparison_table.png"),
       full_grob, width=13, height=4 + length(grob_list)*4.5, dpi=300)

cat("Saved model_comparison_table.pdf/.png\n")
cat("Output:", out_dir, "\n")
