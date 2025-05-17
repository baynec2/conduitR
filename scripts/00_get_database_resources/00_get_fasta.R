################################################################################
# Getting Fasta Files from organism IDS
################################################################################
## Defining inputs and outputs from snakemake workflow
 input_file <- snakemake@input[[1]]
 output_file <- snakemake@output[[1]]

# Make the database_resources_directory if it doesn't exist.
if (!dir.exists(dirname(output_file))) {
  dir.create(dirname(output_file))
}
# Read organism IDs from the input file
organism_ids <- readr::read_delim(input_file) |>
  dplyr::pull(organism_id)

# Create a fasta file by hitting uniprot API via R.
# Warnings are suppressed because the df containing reference proteomes has
# a number of rows that could change with future updates
suppressWarnings(conduitR::download_fasta_from_organism_ids(organism_ids,
                                        destination_fp = output_file))

