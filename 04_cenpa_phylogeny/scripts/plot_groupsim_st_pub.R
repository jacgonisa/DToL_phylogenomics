#!/usr/bin/env Rscript
# GroupSim — Satellite vs Transposon CENP-A (clade-weighted)
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

BASE    <- "/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/04_cenpa_phylogeny"
SP_DIR  <- file.path(BASE, "split_entropy")
FIG_DIR <- file.path(BASE, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

C_UW    <- "#cfd8dc"
C_WT    <- "#1565c0"   # same blue as CENP-A vs H3
C_SIG   <- "#c62828"
C_CENPA <- "#f48fb1"

# ── Load data ─────────────────────────────────────────────────────────────────
# Unweighted (already 1-based for sat_trans)
uw_raw <- readLines(file.path(SP_DIR, "groupsim_sat_trans/groupsim_sat_trans_gap085.txt"))
uw_raw <- uw_raw[!grepl("^#", uw_raw) & nchar(trimws(uw_raw)) > 0]
uw <- do.call(rbind, lapply(uw_raw, function(l) {
  p <- strsplit(l, "\t")[[1]]
  if (length(p) < 2) return(NULL)
  data.frame(pos      = as.integer(p[1]) + 1L,  # 0-based → 1-based
             score_uw = suppressWarnings(as.numeric(p[2])))
})) %>% filter(!is.na(pos))

wt <- read.delim(file.path(SP_DIR,
                            "groupsim_sat_trans_weighted/groupsim_st_clade_gap085.tsv"),
                 stringsAsFactors = FALSE)
names(wt)[names(wt)=="groupsim_clade"] <- "score_wt"
names(wt)[names(wt)=="z_clade"]        <- "z_score"

df <- merge(uw, wt[,c("pos","score_wt","z_score","sig_clade")],
            by="pos", all=TRUE)
df <- df[order(df$pos),]
cat("Positions:", nrow(df), "| sig (z>=2):", sum(df$sig_clade==1, na.rm=TRUE), "\n")

# Helix — plain "start end" format, CENPA only
helix_raw <- readLines(file.path(SP_DIR,
                                  "groupsim_sat_trans/helix_positions_gap085.txt"))
helix_raw <- helix_raw[!grepl("^#", helix_raw) & nchar(trimws(helix_raw)) > 0]
helix <- do.call(rbind, lapply(helix_raw, function(l) {
  p <- as.integer(strsplit(trimws(l), "\\s+")[[1]])
  if (length(p) < 2) return(NULL)
  data.frame(start=p[1], end=p[2])
}))

# no significant positions — no annotations

x_min <- min(df$pos, na.rm=TRUE) - 0.5
x_max <- max(df$pos, na.rm=TRUE) + 0.5

# ── Panel 1: helix track (CENPA only) ────────────────────────────────────────
p_helix <- ggplot(helix) +
  geom_rect(aes(xmin=start-0.5, xmax=end+0.5, ymin=0.1, ymax=0.9),
            fill=C_CENPA, alpha=0.85, colour=NA) +
  scale_x_continuous(limits=c(x_min,x_max), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  annotate("text", x=x_min+2, y=0.5, label="CENP-A helices",
           hjust=0, size=3.2, colour="#880e4f") +
  theme_void(base_size=9) +
  theme(plot.margin=margin(2,2,0,2))

# ── Panel 2: GroupSim bars ────────────────────────────────────────────────────
df_uw <- df[!is.na(df$score_uw),]; df_uw$fill_col <- "uw"
df_wt <- df[!is.na(df$score_wt),]; df_wt$fill_col <- "wt"

p_main <- ggplot() +
  geom_col(data=df_uw, aes(x=pos-0.22, y=score_uw, fill=fill_col),
           width=0.40, alpha=1.0) +
  geom_col(data=df_wt, aes(x=pos+0.22, y=score_wt, fill=fill_col),
           width=0.40, alpha=0.92) +
  scale_fill_manual(
    values = c(uw=C_UW, wt=C_WT),
    labels = c(uw="Unweighted", wt="Clade-weighted"),
    name   = NULL,
    breaks = c("uw","wt")
  ) +
  scale_x_continuous(limits=c(x_min,x_max), expand=c(0,0)) +
  scale_y_continuous(name="GroupSim score", limits=c(0,1.10),
                     breaks=c(0,0.25,0.5,0.75,1.0)) +
  labs(x=NULL) +
  theme_bw(base_size=11) +
  theme(axis.text.x        = element_blank(),
        axis.ticks.x       = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position    = "right",
        legend.text        = element_text(size=10),
        legend.key.size    = unit(0.5,"cm"),
        plot.margin        = margin(0,2,0,2))

# ── Panel 3: z-score strip ────────────────────────────────────────────────────
df_z <- df[!is.na(df$z_score),]

p_zscore <- ggplot(df_z, aes(x=pos, y=0.5, fill=z_score)) +
  geom_tile(width=1, height=1) +
  scale_fill_gradient2(low="#2166ac", mid="white", high="#c62828",
                       midpoint=0, limits=c(-4,4), name="z-score",
                       guide=guide_colorbar(barwidth=1.0, barheight=5,
                                            title.position="top",
                                            title.theme=element_text(size=10),
                                            label.theme=element_text(size=9))) +
  scale_x_continuous(limits=c(x_min,x_max), expand=c(0,0),
                     name="Alignment position (trimmed, gap <= 85%)") +
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  theme_void(base_size=9) +
  theme(axis.title.x      = element_text(size=10, margin=margin(t=4)),
        axis.text.x       = element_text(size=9),
        axis.ticks.x      = element_line(linewidth=0.3),
        legend.position   = "right",
        legend.text       = element_text(size=9),
        legend.title      = element_text(size=10),
        panel.background  = element_rect(fill="#cccccc", colour=NA),
        plot.margin       = margin(0,2,2,2))

# ── Combine ───────────────────────────────────────────────────────────────────
p_out <- p_helix / p_main / p_zscore +
  plot_layout(guides="collect", heights=c(0.5, 6, 0.5)) +
  plot_annotation(
    title    = "Satellite CENP-A (n=198) vs Transposon CENP-A (n=129)  |  clade-weighted GroupSim",
    subtitle = "No positions reach z >= 2. Top 6 candidates labelled. Grey = unweighted | Orange = clade-weighted.",
    theme    = theme(plot.title    = element_text(face="bold", size=11),
                     plot.subtitle = element_text(size=8, colour="grey40"),
                     plot.background = element_rect(fill="white", colour=NA))
  )

ggsave(file.path(FIG_DIR,"groupsim_sat_vs_trans_pub.pdf"),
       p_out, width=13, height=6.5, bg="white")
ggsave(file.path(FIG_DIR,"groupsim_sat_vs_trans_pub.png"),
       p_out, width=13, height=6.5, dpi=300, bg="white")
cat("Saved: groupsim_sat_vs_trans_pub.{pdf,png}\n")
