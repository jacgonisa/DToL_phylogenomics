#!/usr/bin/env Rscript
# cenpb_abundance_table.R
# Merge the two CENP-B box detection methods into one per-species abundance
# table + a per-clade summary figure.
#   scan    (cenpb_box_scan.py):   canonical box with <=2 / <=3 mismatches, % of ALL monomers, vs shuffled null
#   evidence(cenpb_box_finder.py): 4 independent lines (canonical/broad/degenerate IUPAC + PWM/FIMO), counts per 2000 sampled
suppressPackageStartupMessages({library(dplyr); library(gridExtra); library(grid); library(ggplot2)})

SAT <- "/home/jg2070/Desktop/dtol_review_August/DToL_phylogenomics_publication_325genomes/05_satellite_similarity"
scan <- read.delim(file.path(SAT,"figures/cenpb_box_per_species.tsv"), stringsAsFactors=FALSE)
evid <- read.delim(file.path(SAT,"figures/cenpb_box_evidence_per_species.tsv"), stringsAsFactors=FALSE)

broad_map <- c(Mammalia="Vertebrates",Aves="Vertebrates",Actinopterygii="Vertebrates",
  Reptilia="Vertebrates",Amphibia="Vertebrates",Chondrichthyes="Vertebrates",
  Fungi="Fungi",Algae="Viridiplantae",Bryophyta="Viridiplantae",Dicots="Viridiplantae",
  Monocots="Viridiplantae",Gymnosperms="Viridiplantae",Alveolata="Protist",Discoba="Protist")
clade_of <- function(t) ifelse(t %in% names(broad_map), broad_map[t], "Invertebrate")

m <- full_join(scan, evid %>% select(-taxa1), by="species")
m$clade <- clade_of(m$taxa1)
m$box_le2mm_count <- round(m$n_monomers * m$pct_box_le2mm/100)

# ---- per-species merged table (all methods) ----
per_sp <- m %>% transmute(
  species, clade, taxa1,
  n_monomers,
  pct_box_le2mm = round(pct_box_le2mm,2),
  pct_box_le3mm = round(pct_box_le3mm,2),
  pct_null_le3mm= round(pct_null_le3mm,2),
  enrichment_le3mm = round(enrichment_le3mm,2),
  box_le2mm_count,
  ev_canonical = canonical, ev_broad = broad, ev_degenerate = degenerate, ev_pwm_fimo = pwm,
  lines_positive) %>%
  arrange(desc(pct_box_le2mm))
write.table(per_sp, file.path(SAT,"figures/cenpb_abundance_per_species_325sp.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE, na="")
cat("per-species rows:", nrow(per_sp), "\n")

# ---- per-clade summary ----
per_clade <- m %>% group_by(clade) %>% summarise(
  n_species          = n(),
  `mean %box <=2mm`  = round(mean(pct_box_le2mm, na.rm=TRUE),3),
  `mean %box <=3mm`  = round(mean(pct_box_le3mm, na.rm=TRUE),3),
  `mean %null <=3mm` = round(mean(pct_null_le3mm, na.rm=TRUE),3),
  `enrichment`       = round(mean(pct_box_le3mm,na.rm=TRUE)/mean(pct_null_le3mm,na.rm=TRUE),2),
  `canonical (/2000)`= round(mean(canonical, na.rm=TRUE),1),
  `broad`            = round(mean(broad, na.rm=TRUE),1),
  `degenerate`       = round(mean(degenerate, na.rm=TRUE),1),
  `PWM/FIMO`         = round(mean(pwm, na.rm=TRUE),1),
  `n >=2 lines pos`  = sum(lines_positive>=2, na.rm=TRUE),
  .groups="drop") %>%
  arrange(desc(`enrichment`))
write.table(per_clade, file.path(SAT,"figures/cenpb_abundance_per_clade_325sp.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
print(as.data.frame(per_clade))

# ---- render per-clade table figure ----
clade_bg <- c(Vertebrates="#E3F2FD",Invertebrate="#FFF3E0",Viridiplantae="#E8F5E9",
              Fungi="#F3E5F5",Protist="#FFEBEE")
th <- ttheme_minimal(base_size=9,
  core=list(fg_params=list(hjust=0.5),
            bg_params=list(fill=clade_bg[per_clade$clade], col="grey85", lwd=0.4)),
  colhead=list(fg_params=list(fontface=2, col="white", parse=FALSE),
               bg_params=list(fill="#37474F", col=NA)))
g <- tableGrob(per_clade, rows=NULL, theme=th)
title <- textGrob("CENP-B box abundance across the 325-sp DToL satellites, by clade",
                  gp=gpar(fontface="bold", fontsize=12), x=0.01, hjust=0)
sub <- textGrob(paste("Scan method: canonical 17-bp box [CT]TTCGTTGGAA[AG]CGGGA with <=2/<=3 mismatches, both strands, % of all monomers vs shuffled null.",
                "Evidence method: mean counts per 2000 sampled monomers for 4 independent lines (canonical/broad/degenerate IUPAC exact + PWM FIMO)."),
                gp=gpar(fontsize=7, col="grey35"), x=0.01, hjust=0)
full <- arrangeGrob(title, sub, g, heights=unit.c(unit(1.4,"lines"), unit(2.2,"lines"), unit(1,"null")))
for (ext in c("pdf","png")) {
  out <- file.path(SAT,"figures", paste0("cenpb_abundance_by_clade_325sp.", ext))
  ggsave(out, full, width=12, height=3.2, dpi=300, bg="white")
  cat("Saved:", out, "\n")
}
