hor_df_file <- ""

hors <- read.csv(hor_df_file)

hors_dt <- as.data.table(hors)[, .(
  start_A = as.integer(start_A),
  end_A   = as.integer(end_A),
  start_B = as.integer(start_B),
  end_B   = as.integer(end_B),
  bs      = as.integer(block.size.in.units)
)]

# Build the pair table by expanding only the *ranges*, not the block size
# (we only care about unique element-partner pairs, so multiplicity is irrelevant)
make_pairs <- function(sA, eA, sB, eB) {
  # one row per HOR
  data.table(
    element = Map(seq.int, sA, eA),
    partner = Map(seq.int, sB, eB)
  )[, .(
    element = unlist(element),
    partner = unlist(partner)
  )]
}

# Direction A→B
pairs1 <- make_pairs(hors_dt$start_A, hors_dt$end_A,
                     hors_dt$start_B, hors_dt$end_B)
# Direction B→A
pairs2 <- make_pairs(hors_dt$start_B, hors_dt$end_B,
                     hors_dt$start_A, hors_dt$end_A)

long_data <- unique(rbind(pairs1, pairs2))
rm(pairs1, pairs2, hors_dt); gc()

# Drop self-pairs
long_data <- long_data[element != partner]

# Count unique partners per element
unique_interactions <- long_data[, .(num_interactors = uniqueN(partner)), by = element]
rm(long_data); gc()

# Write into the result vector
result_vector <- integer(nrow(repeats_chr))
result_vector[unique_interactions$element] <- unique_interactions$num_interactors
rm(unique_interactions)

repeats_chr$HOR_score <- 100 * result_vector / nrow(repeats_chr)