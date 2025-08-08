# =============================================================================
# Build Proteotyping FASTA Database
# =============================================================================
# This script builds a FASTA database for proteotyping by filtering peptides
# that have a Lowest Common Ancestor (LCA) at the species level.

# =============================================================================
# Setup and Logging
# =============================================================================

# Open log file for both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

# Record start time
start_time <- Sys.time()
conduitR::log_with_timestamp("Starting build_proteotyping_fasta_db.R script")

# =============================================================================
# Input/Output File Definitions
# =============================================================================

# Input files
sequences_file <- snakemake@input[["sequences.tsv.lz4"]]
taxons_file <- snakemake@input[["taxons.tsv.lz4"]]

# Output file
output_fasta_file <- snakemake@output[["first_pass_fasta"]]

conduitR::log_with_timestamp("Input files: %s, %s", sequences_file, taxons_file)
conduitR::log_with_timestamp("Output file: %s", output_fasta_file)

# =============================================================================
# Data Loading
# =============================================================================

conduitR::log_with_timestamp("Reading sequence file")

# Read the LZ4 compressed sequence file
sequences <- readr::read_tsv(
  pipe(paste("lz4 -d -c", shQuote(sequences_file))),
  col_names = c("id", "sequence", "lca", "lca_il", "fa", "fa_il"),
  col_types = readr::cols(
    id = "c", sequence = "c", lca = "c", 
    lca_il = "c", fa = "c", fa_il = "c"
  )
)

conduitR::log_with_timestamp("Completed reading sequence file (%d rows)", nrow(sequences))

conduitR::log_with_timestamp("Reading taxon file")

# Read the LZ4 compressed taxon file
taxons <- readr::read_tsv(
  pipe(paste("lz4 -d -c", shQuote(taxons_file))),
  col_names = c("id", "name", "rank", "parent_id"),
  col_types = readr::cols(
    id = "c", name = "c", rank = "c", parent_id = "c"
  )
)

conduitR::log_with_timestamp("Completed reading taxon file (%d rows)", nrow(taxons))

# =============================================================================
# Data Processing
# =============================================================================

conduitR::log_with_timestamp("Filtering peptides to species-level LCA")

# Filter to only contain peptides that have an LCA at the species level
species_peptides <- dplyr::left_join(
  sequences, 
  taxons, 
  by = c("lca_il" = "id")
) |>
  dplyr::filter(rank == "species") |>
  dplyr::mutate(fasta_header = paste(id, lca_il, sep = "|"))

conduitR::log_with_timestamp(
  "Filtered to %d peptides with species-level LCA (from %d total)", 
  nrow(species_peptides), 
  nrow(sequences)
)

# =============================================================================
# FASTA File Generation
# =============================================================================

conduitR::log_with_timestamp("Creating FASTA file")

# Create Biostrings object and write to file
fasta_sequences <- Biostrings::AAStringSet(species_peptides$sequence)
names(fasta_sequences) <- species_peptides$fasta_header

Biostrings::writeXStringSet(fasta_sequences, filepath = output_fasta_file)

conduitR::log_with_timestamp("FASTA file saved to: %s", output_fasta_file)

# =============================================================================
# Cleanup and Logging
# =============================================================================

end_time <- Sys.time()
elapsed_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))

conduitR::log_with_timestamp(
  "Completed build_proteotyping_fasta_db.R script. Time taken: %.2f minutes", 
  elapsed_minutes
)

# Close log file connections
sink(type = "message")
sink()
close(log_connection)