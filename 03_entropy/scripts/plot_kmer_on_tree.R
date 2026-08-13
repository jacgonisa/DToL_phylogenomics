#!/usr/bin/env Rscript
# Plot k-mer Jaccard similarity (k=15) mapped onto the 325-sp tree.
# Same layout as plot_seqsim_on_tree.R — internal nodes coloured by
# mean cross-subtree Jaccard; size = number of species pairs compared.
# Jaccard = shared_hashes / total_hashes (col5 of mash matrix, NOT col3).

suppressPackageStartupMessages({
  library(ape); library(phangorn); library(ggtree); library(ggplot2)
  library(dplyr); library(ggnewscale)
})

BASE    <- "/home/jg2070/Desktop/dtol_review_August"
PUB     <- file.path(BASE, "DToL_phylogenomics_publication_325genomes")
FIG_DIR <- file.path(PUB, "03_entropy/figures/satellite_halflife")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

TREE_F  <- file.path(PUB, "01_species_tree/outputs/full_325sp_calibrated.nwk")
KMER_F  <- file.path(BASE, "matrix_katie_k15.txt")
SAT_F   <- file.path(BASE, "2026_trees/annotation_centromeres/organized/plots",
                     "centromere_length/genome_vs_centromere_length_satellite_species_table.tsv")
TAX_F   <- file.path(BASE, "2026_trees/annotation_centromeres/centromere_code_to_species.tsv")

# ── tree ──────────────────────────────────────────────────────────────────────
tr <- read.tree(TREE_F)
tr$tip.label <- sub("\\.fa$", "", tr$tip.label)
tip_code <- sub("[0-9.].*$", "", tolower(tr$tip.label))

# ── parse k-mer matrix ────────────────────────────────────────────────────────
cat("Parsing k-mer matrix …\n")
raw <- readLines(KMER_F)
fields <- strsplit(raw, "\t")
ok     <- sapply(fields, length) == 5
mat    <- as.data.frame(do.call(rbind, fields[ok]), stringsAsFactors = FALSE)
colnames(mat) <- c("fileA","fileB","mash_dist","pval","frac")

# species code from filename: FASTAS/bcorhaw_1398.fasta → bcorhaw
parse_sp <- function(x) tolower(sub("_[0-9]+$", "", sub(".*/([^/]+)\\.fasta?$","\\1",x)))
mat$spA <- parse_sp(mat$fileA)
mat$spB <- parse_sp(mat$fileB)
mat     <- mat[mat$spA != mat$spB, ]

# true Jaccard = num/denom from col5
frac_parts   <- strsplit(mat$frac, "/")
mat$jaccard  <- as.numeric(sapply(frac_parts, `[`, 1)) /
                as.numeric(sapply(frac_parts, `[`, 2))
mat          <- mat[!is.na(mat$jaccard), ]
cat(sprintf("  %d valid comparisons, %d unique species\n",
            nrow(mat), length(unique(c(mat$spA, mat$spB)))))

# per-species-pair mean Jaccard
pair_df <- mat %>%
  mutate(a = pmin(spA,spB), b = pmax(spA,spB)) %>%
  group_by(a, b) %>%
  summarise(mean_jaccard = mean(jaccard, na.rm=TRUE),
            n_comp = n(), .groups="drop")

# lookup function
jac_lookup <- setNames(pair_df$mean_jaccard,
                       paste(pair_df$a, pair_df$b, sep="|"))
get_jac <- function(cA, cB) {
  key <- paste(sort(c(cA, cB)), collapse="|")
  v   <- jac_lookup[key]
  if (is.na(v)) NA_real_ else as.numeric(v)
}

# ── taxonomy ──────────────────────────────────────────────────────────────────
tax <- read.delim(TAX_F, stringsAsFactors=FALSE)
tax$code <- sub("[0-9.].*$","",tolower(tax$fasta_base))
vert  <- c("Actinopterygii","Aves","Mammalia","Reptilia","Amphibia","Chondrichthyes")
virid <- c("Algae","Bryophyta","Dicots","Monocots","Gymnosperms")
tax$broad <- ifelse(tax$taxa1 %in% vert, "Vertebrates",
             ifelse(tax$taxa1 %in% virid, "Viridiplantae",
             ifelse(tax$taxa1=="Fungi","Fungi","Invertebrates")))
code2broad <- setNames(tax$broad, tax$code)

sat     <- read.delim(SAT_F, stringsAsFactors=FALSE)
sat_codes <- sub("[0-9.].*$","",tolower(sat$species_id[sat$is_satellite_species]))

# ── per-node k-mer similarity ─────────────────────────────────────────────────
cat("Computing per-node Jaccard similarity …\n")
n_tips  <- Ntip(tr); n_nodes <- Nnode(tr)
desc_tips <- lapply(seq_len(n_tips+n_nodes), function(nd)
  unlist(Descendants(tr, nd, type="tips")))

node_jac   <- rep(NA_real_, n_nodes)
node_npair <- rep(NA_integer_, n_nodes)

for (i in seq_len(n_nodes)) {
  nd  <- n_tips + i
  ch  <- tr$edge[tr$edge[,1]==nd, 2]
  if (length(ch) < 2) next

  lc <- intersect(unlist(desc_tips[[ch[1]]]), seq_len(n_tips))
  rc <- intersect(unlist(desc_tips[[ch[2]]]), seq_len(n_tips))

  lc <- unique(tip_code[lc]); lc <- lc[lc %in% sat_codes]
  rc <- unique(tip_code[rc]); rc <- rc[rc %in% sat_codes]
  if (!length(lc) || !length(rc)) next

  vals <- na.omit(as.vector(outer(lc, rc, Vectorize(get_jac))))
  if (!length(vals)) next
  node_jac[i]   <- mean(vals)
  node_npair[i] <- length(vals)
}

n_col <- sum(!is.na(node_jac))
cat(sprintf("  Nodes with data: %d / %d  (uncoloured = no satellite spp on one side)\n",
            n_col, n_nodes))
cat(sprintf("  Jaccard range: %.4f – %.4f\n",
            min(node_jac, na.rm=TRUE), max(node_jac, na.rm=TRUE)))

# ── clade MRCA ────────────────────────────────────────────────────────────────
get_mrca <- function(taxa_vec) {
  codes <- tax$code[tax$taxa1 %in% taxa_vec]
  tips  <- tr$tip.label[tip_code %in% codes]
  if (length(tips) < 2) return(NA_integer_)
  getMRCA(tr, tips)
}
clade_df <- data.frame(
  node = c(
    get_mrca(c("Actinopterygii","Aves","Mammalia","Reptilia","Amphibia","Chondrichthyes")),
    get_mrca("Aves"), get_mrca("Mammalia"), get_mrca("Actinopterygii"),
    get_mrca(c("Coleoptera","Diptera","Hymenoptera","Lepidoptera","Hemiptera",
               "Ephemeroptera","Odonata","Trichoptera","Neuroptera","Blattodea")),
    get_mrca(c("Bryophyta","Dicots","Monocots","Algae","Gymnosperms")),
    get_mrca("Monocots"), get_mrca("Dicots"), get_mrca("Fungi")
  ),
  label = c("Vertebrates","Aves","Mammalia","Actinopterygii",
            "Insects","Viridiplantae","Monocots","Dicots","Fungi"),
  stringsAsFactors = FALSE
) %>% filter(!is.na(node))

# ── build combined data for ggtree ────────────────────────────────────────────
tip_df <- data.frame(
  label   = tr$tip.label,
  code    = tip_code,
  is_sat  = tip_code %in% sat_codes,
  broad   = code2broad[tip_code],
  stringsAsFactors = FALSE
)
node_df <- data.frame(
  node    = seq(n_tips+1, n_tips+n_nodes),
  jac_sim = node_jac,
  n_pairs = node_npair
)
all_df <- bind_rows(
  tip_df  %>% mutate(node = match(label, tr$tip.label),
                     jac_sim = NA_real_, n_pairs = NA_integer_),
  node_df %>% mutate(label=NA_character_, code=NA_character_,
                     is_sat=NA, broad=NA_character_)
)

broad_pal <- c(Vertebrates="#1565C0", Invertebrates="#E65100",
               Viridiplantae="#2E7D32", Fungi="#6A1B9A")

# ── plot ──────────────────────────────────────────────────────────────────────
p <- ggtree(tr, layout="fan", linewidth=0.12, color="grey70") %<+% all_df +
  geom_tippoint(aes(color=broad, subset=!is.na(is_sat) & is_sat),
                size=1.0, alpha=0.85) +
  scale_color_manual(values=broad_pal, name="Group (satellite spp.)",
                     na.value="#cccccc") +
  new_scale_color() +
  geom_nodepoint(aes(color=jac_sim, size=n_pairs,
                     subset=!is.na(jac_sim)),
                 alpha=0.90) +
  scale_size_continuous(name="Species pairs\nat node", range=c(0.8,7),
                        breaks=c(1,5,20,100),
                        guide=guide_legend(override.aes=list(alpha=1,color="#fc8d59"))) +
  scale_color_gradientn(
    colors  = c("#f7f7f7","#fee08b","#fc8d59","#d73027","#a50026"),
    limits  = c(0, 0.30),
    name    = "Mean k-mer\nJaccard (k=15)",
    na.value= NA,
    guide   = guide_colorbar(title.position="top", barwidth=0.8, barheight=5)
  ) +
  geom_cladelab(data=clade_df, mapping=aes(node=node, label=label),
                fontsize=2.8, offset=3, offset.text=1.5,
                barsize=0.4, barcolor="grey50", textcolor="grey20",
                align=TRUE, hjust=0.5) +
  labs(
    title    = "Satellite k-mer similarity (Jaccard, k=15) mapped onto the 325-sp tree",
    subtitle = paste0(
      "Node colour: mean Jaccard (shared k-mer hashes / total) between satellite species\n",
      "in daughter subtrees. Node size: number of species pairs. ",
      "Uncoloured = no satellite spp. on one/both sides.")
  ) +
  theme(plot.title=element_text(face="bold",size=11,hjust=0.5),
        plot.subtitle=element_text(size=7,hjust=0.5,color="grey40"),
        legend.position="right")

ggsave(file.path(FIG_DIR,"kmer_on_tree.pdf"), p, width=14, height=14)
ggsave(file.path(FIG_DIR,"kmer_on_tree.png"), p, width=14, height=14, dpi=300)
cat("Saved: kmer_on_tree.{pdf,png}\n")
