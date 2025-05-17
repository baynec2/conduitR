################################################################################
# Extract Protein Info From Fasta
################################################################################
# Get input and output files from Snakemake workflow
database_fasta <- snakemake@input[["database_fasta"]]
taxonomy_txt <- snakemake@input[["taxonomy_txt"]]
output_file <- snakemake@output[[1]]

# #Testing
# database_fasta <- "user_input/00_database_resources/00_database.fasta"
# taxonomy_txt <- "user_input/00_database_resources/01_taxonomy.txt"
# output_file <- "user_input/00_database_resources/02_protein_info.txt"

# Loading organism annotation
taxonomy_txt <- readr::read_delim(taxonomy_txt)

# Extracting protein information from fasta database
protein_info <- conduitR::extract_fasta_info(database_fasta)

# Adding organism annotation to protein info
protein_info <- dplyr::inner_join(protein_info, taxonomy_txt,
  by = "organism_id"
)

readr::write_delim(protein_info, output_file, delim = "\t")
